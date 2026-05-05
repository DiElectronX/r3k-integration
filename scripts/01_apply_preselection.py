import ROOT
import os
import yaml
import numpy as np
import pandas as pd
import argparse

# --------------------------------------------------------------------------------
# DEBUG HELPERS
# --------------------------------------------------------------------------------
DEBUG_SAMPLES = {"jpsi_resonant", "kaon_psi2s_resonant", "kaon_lowq_resonant"}


def should_debug(sample_name):
    return sample_name in DEBUG_SAMPLES


def n_unique_events(rdf):
    arr = rdf.AsNumpy(["event_key"])["event_key"]
    return len(np.unique(arr))


def debug_stage(sample_name, label, rdf):
    cols = [str(c) for c in rdf.GetColumnNames()]
    n_cand = rdf.Count().GetValue()
    n_evt = n_unique_events(rdf)

    msg = f"[DEBUG][{sample_name}][{label}] candidates={n_cand}, unique_events={n_evt}"
    if "FONLLweight" in cols:
        msg += f", sum(FONLLweight)={rdf.Sum('FONLLweight').GetValue()}"
    if "total_weight" in cols:
        msg += f", sum(total_weight)={rdf.Sum('total_weight').GetValue()}"
    print(msg)


def debug_trigger_paths(sample_name, rdf, trigger_data):
    frac_sum = sum(v[0] for v in trigger_data.values())
    print(f"[DEBUG][{sample_name}] trigger fraction sum = {frac_sum}")

    for key in trigger_data:
        n_assigned = rdf.Filter(f'assigned_path == "{key}"').Count().GetValue()
        n_pass = rdf.Filter(f'assigned_path == "{key}" && {key}_pass').Count().GetValue()
        w_pass = rdf.Filter(f'assigned_path == "{key}" && {key}_pass').Sum('total_weight').GetValue()

        print(
            f"[DEBUG][{sample_name}][path={key}] "
            f"cfg_frac={trigger_data[key][0]:.6f}, "
            f"assigned_cands={n_assigned}, "
            f"pass_cands={n_pass}, "
            f"pass_sumw={w_pass}"
        )


# --------------------------------------------------------------------------------
# C++ HELPER
# --------------------------------------------------------------------------------
def declare_cpp_helper(trigger_data):
    if hasattr(ROOT, 'assign_path') and hasattr(ROOT, 'det_rand'):
        return

    paths = list(trigger_data.keys())
    fractions = [trigger_data[k][0] for k in paths]
    norm_fracs = np.array(fractions) / sum(fractions)
    cdf = np.cumsum(norm_fracs)

    cpp_code = f'''
    #include <vector>
    #include <string>
    #include <cstdint>

    double det_rand(unsigned long long key, unsigned int seed) {{
        uint64_t x = (uint64_t)key + 0x9e3779b97f4a7c15ULL + (uint64_t)seed;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        x = (x ^ (x >> 31));
        return (x >> 11) * (1.0 / 9007199254740992.0);
    }}

    std::string assign_path(double rand) {{
        const std::vector<std::string> paths = {{ {", ".join(f'"{p}"' for p in paths)} }};
        const std::vector<double> cdf = {{ {", ".join(str(x) for x in cdf)} }};
        for (size_t i = 0; i < cdf.size(); ++i) {{
            if (rand < cdf[i]) return paths[i];
        }}
        return paths.back();
    }}
    '''
    ROOT.gInterpreter.Declare(cpp_code)


def build_sf_expression(sf_bins):
    val_parts, unc_parts = [], []
    for ptmin, ptmax, val, unc in sf_bins:
        ptmax = float(ptmax)
        if ptmax == float('inf'):
            cond = f'(BToKEE_fit_l2_pt >= {ptmin})'
        else:
            cond = f'(BToKEE_fit_l2_pt >= {ptmin} && BToKEE_fit_l2_pt < {ptmax})'
        val_parts.append(f'{cond} ? {val}')
        unc_parts.append(f'{cond} ? {unc}')
    return ' : '.join(val_parts) + ' : 1.0', ' : '.join(unc_parts) + ' : 0.0'


# --------------------------------------------------------------------------------
# MAIN PROCESSING LOGIC
# --------------------------------------------------------------------------------
def process_file(file_info, config, base_output_dir, args):
    input_path = file_info['path']
    sample_name = file_info['name']
    is_data = file_info.get('is_data', False)
    target_regions = file_info.get('regions', ['none'])

    if not os.path.exists(input_path):
        print(f"Skipping missing: {input_path}")
        return []

    print(f"\nProcessing {sample_name}...")

    # 0. Load Dataframe + define event_key immediately for all samples
    df = (
        ROOT.RDataFrame("Events", input_path)
        .Define(
            'event_key',
            '(static_cast<unsigned long long>(run) << 42) ^ '
            '(static_cast<unsigned long long>(luminosityBlock) << 20) ^ '
            'static_cast<unsigned long long>(event)'
        )
    )

    # -----------------------------------------------------------------------
    # STEP 1: GLOBAL DEFINITIONS
    # -----------------------------------------------------------------------
    if is_data:
        df = (
            df.Define('FONLLweight', '1.0')
              .Define('trigger_sf_value', '1.0')
              .Define('trigger_sf_error', '0.0')
              .Define('total_weight', '1.0')
        )
    else:
        trig_sf_val, trig_sf_unc = build_sf_expression(config['trigger_sf_params']['bins'])

        if 'FONLLweight' not in [str(c) for c in df.GetColumnNames()]:
            df = df.Define('FONLLweight', '0.7878582119')

        df = (
            df.Define('trigger_sf_value', trig_sf_val)
              .Define('trigger_sf_error', trig_sf_unc)
              .Define('total_weight', 'FONLLweight * trigger_sf_value')
        )

    # -----------------------------------------------------------------------
    # STEP 2: DEFINE TRIGGER LOGIC
    # -----------------------------------------------------------------------
    trigger_data = config['trigger']['fractions']

    if not is_data and args.mode == 'mix':
        declare_cpp_helper(trigger_data)
 
        seed = int(config['random_seed'])

        # IMPORTANT FIX:
        # use event_key, not rdfentry_, so all candidates from the same event
        # get the same assigned trigger path
        df = (
            df.Define('rand', f'det_rand(event_key, {seed})')
              .Define('assigned_path', 'assign_path(rand)')
        )

        pass_cols = []
        for key, val in trigger_data.items():
            l1, hlt = val[1], val[2]
            col_name = f'{key}_pass'
            df = df.Define(col_name, f'(assigned_path == "{key}") && ({l1} && {hlt})')
            pass_cols.append(col_name)

        trig_filter_string = ' || '.join(pass_cols)

        if should_debug(sample_name):
            debug_trigger_paths(sample_name, df, trigger_data)

    elif not is_data and args.mode == 'single':
        target = args.target_trigger
        l1_t, hlt_t = trigger_data[target][1], trigger_data[target][2]
        trig_filter_string = f'({l1_t} && {hlt_t})'

    else:
        if args.mode == 'single':
            target = args.target_trigger
            trig_filter_string = f'{trigger_data[target][1]} && {trigger_data[target][2]}'
        else:
            others = [f'({v[1]} && {v[2]})' for _, v in trigger_data.items()]
            trig_filter_string = ' || '.join(others)

    if should_debug(sample_name):
        print(f"\n[DEBUG][{sample_name}] input_path = {input_path}")
        print(f"[DEBUG][{sample_name}] triplet cut    = {config['preselection']['triplet']}")
        print(f"[DEBUG][{sample_name}] trigger filter = {trig_filter_string}")
        print(f"[DEBUG][{sample_name}] presel BDT cut = {config['preselection']['bdt_score']}")
        print(f"[DEBUG][{sample_name}] anti-D0 cut    = {config['preselection']['anti_d0']}")

    # -----------------------------------------------------------------------
    # STEP 3: SEQUENTIAL FILTERING
    # -----------------------------------------------------------------------
    df_triplet = df.Filter(config['preselection']['triplet'])
    df_trig = df_triplet.Filter(trig_filter_string)
    df_presel = df_trig.Filter(config['preselection']['bdt_score'])
    df_antid0 = df_presel.Filter(config['preselection']['anti_d0'])

    if should_debug(sample_name):
        debug_stage(sample_name, "RAW", df)
        debug_stage(sample_name, "AFTER_TRIPLET", df_triplet)
        debug_stage(sample_name, "AFTER_TRIGGER", df_trig)
        debug_stage(sample_name, "AFTER_PRESEL_BDT", df_presel)
        debug_stage(sample_name, "AFTER_ANTI_D0", df_antid0)

    # -----------------------------------------------------------------------
    # STEP 4: COMMON COUNTS
    # -----------------------------------------------------------------------
    file_results = []

    count_before = df.Sum('FONLLweight').GetValue()
    count_triplet = df_triplet.Sum('FONLLweight').GetValue()
    count_trig = df_trig.Sum('FONLLweight').GetValue()
    sow_trig = df_trig.Sum('total_weight').GetValue()
    sow_presel = df_presel.Sum('total_weight').GetValue()
    sow_antid0 = df_antid0.Sum('total_weight').GetValue()

    # -----------------------------------------------------------------------
    # STEP 5: REGION LOOP
    # -----------------------------------------------------------------------
    for region in target_regions:
        region_dir = os.path.join(base_output_dir, region)
        print(region_dir)
        os.makedirs(region_dir, exist_ok=True)

        mode_suffix = f"_{args.target_trigger}" if args.mode == 'single' else "_mix"
        out_name = os.path.basename(input_path).replace('.root', f'_skimmed{mode_suffix}.root')
        out_path = os.path.join(region_dir, out_name)

        q2_cut_str = config['q2_cuts'].get(region, "1")
        df_final = df_antid0.Filter(q2_cut_str)

        if should_debug(sample_name):
            print(f"[DEBUG][{sample_name}][{region}] q2 cut = {q2_cut_str}")
            debug_stage(sample_name, f"AFTER_Q2_{region}", df_final)

        sow_final = df_final.Sum('total_weight').GetValue()

        file_results.append({
            'Sample Name': sample_name,
            'Region': region,
            'Trigger_Mode': args.target_trigger if args.mode == 'single' else 'Mixture',
            'Before': count_before,
            'After_Triplet': count_triplet,
            'After_Trigger': count_trig,
            'SumWeights_Trigger': sow_trig,
            'After_PreselBDT': sow_presel,
            'After_AntiD0': sow_antid0,
            'After_Q2': sow_final
        })

        print(f"  -> Saving {region} region to {out_path}")
        opts = ROOT.RDF.RSnapshotOptions()
        opts.fMode = "RECREATE"
        df_final.Snapshot("Events", out_path, "", opts)

    return file_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='../config/cuts.yml')
    parser.add_argument('--samples', default='../config/samples.yml')
    parser.add_argument('--mode', choices=['mix', 'single'], default='mix')
    parser.add_argument('--target_trigger', default=None)
    parser.add_argument('--output', default='../data/pre_bdt_ntuples')
    args = parser.parse_args()

    #ROOT.ROOT.EnableImplicitMT()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    with open(args.samples) as f:
        sample_list = yaml.safe_load(f)

    if args.mode == 'single' and not args.target_trigger:
        raise ValueError("Single mode requires --target_trigger")

    all_results = []
    for entry in sample_list['samples']:
        res_list = process_file(entry, config, args.output, args)
        all_results.extend(res_list)

    os.makedirs('../data/logs', exist_ok=True)

    if args.mode == 'single':
        log_name = f'cutflow_step1_{args.target_trigger}.csv'
    else:
        log_name = 'cutflow_step1_mix.csv'

    pd.DataFrame(all_results).to_csv(os.path.join('../data/logs', log_name), index=False)