import ROOT
import os
import yaml
import numpy as np
import pandas as pd
import argparse

# --------------------------------------------------------------------------------
# C++ HELPER (Unchanged)
# --------------------------------------------------------------------------------
def declare_cpp_helper(trigger_data):
    if hasattr(ROOT, 'assign_path'): return
    paths = list(trigger_data.keys())
    fractions = [trigger_data[k][0] for k in paths]
    norm_fracs = np.array(fractions) / sum(fractions)
    cdf = np.cumsum(norm_fracs)

    cpp_code = f'''
    #include <vector>
    #include <string>
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
        if ptmax == float('inf'): cond = f'(BToKEE_fit_l2_pt >= {ptmin})'
        else: cond = f'(BToKEE_fit_l2_pt >= {ptmin} && BToKEE_fit_l2_pt < {ptmax})'
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
    target_regions = file_info.get('regions', ['none']) # Default to list

    if not os.path.exists(input_path):
        print(f"Skipping missing: {input_path}")
        return []

    print(f"\nProcessing {sample_name}...")
    
    # 0. Load Dataframe
    df = ROOT.RDataFrame("Events", input_path)

    # -----------------------------------------------------------------------
    # CRITICAL STEP 1: GLOBAL DEFINITIONS (Must happen BEFORE filtering)
    # This ensures random numbers are assigned to the same Event IDs regardless of cuts
    # -----------------------------------------------------------------------
    
    # Define Weights
    if is_data:
        df = df.Define('FONLLweight', '1.0') \
               .Define('trigger_sf_value', '1.0') \
               .Define('trigger_sf_error', '0.0') \
               .Define('total_weight', '1.0')
    else:
        trig_sf_val, trig_sf_unc = build_sf_expression(config['trigger_sf_params']['bins'])
        
        if 'FONLLweight' not in [str(c) for c in df.GetColumnNames()]:
             df = df.Define('FONLLweight', '1.0')

        df = df.Define('trigger_sf_value', trig_sf_val) \
               .Define('trigger_sf_error', trig_sf_unc) \
               .Define('total_weight', 'FONLLweight * trigger_sf_value')

    # Define Trigger Columns (The "Fragile" Part)
    trigger_data = config['trigger']['fractions']
    
    # Case A: MC Mixture (Requires Random Generation on raw DF)
    if not is_data and args.mode == 'mix':
        ROOT.gRandom.SetSeed(config['random_seed'])
        declare_cpp_helper(trigger_data)
        
        # Define random assignment on the FULL dataset
        df = df.Define('rand', 'gRandom->Rndm()') \
               .Define('assigned_path', 'assign_path(rand)')

        # Define individual pass columns
        pass_cols = []
        for key, val in trigger_data.items():
            l1, hlt = val[1], val[2]
            col_name = f'{key}_pass'
            df = df.Define(col_name, f'(assigned_path == "{key}") && ({l1} && {hlt})')
            pass_cols.append(col_name)
        
        # The filter string is just OR of these columns
        trig_filter_string = ' || '.join(pass_cols)

    # Case B: MC Single / Data (Logic-based, no new columns needed immediately)
    elif not is_data and args.mode == 'single':
        target = args.target_trigger
        l1_t, hlt_t = trigger_data[target][1], trigger_data[target][2]
        
        # OLD LOGIC (N-1):
        # reject = f'!({l1_t} && {hlt_t})'
        # others = [f'({v[1]} && {v[2]})' for k,v in trigger_data.items()]
        # pass_any = '(' + ' || '.join(others) + ')'
        # trig_filter_string = f'{reject} && {pass_any}'
        
        # NEW LOGIC (Pure Selection):
        # "Keep events where THIS trigger fired."
        # (We do not care if other triggers also fired)
        trig_filter_string = f'({l1_t} && {hlt_t})'

    else: # Data
        if args.mode == 'single':
            target = args.target_trigger
            trig_filter_string = f'{trigger_data[target][1]} && {trigger_data[target][2]}'
        else:
            others = [f'({v[1]} && {v[2]})' for k,v in trigger_data.items()]
            trig_filter_string = ' || '.join(others)

    # -----------------------------------------------------------------------
    # STEP 2: SEQUENTIAL FILTERING
    # -----------------------------------------------------------------------
    
    # A. Triplet Cuts
    df_triplet = df.Filter(config['preselection']['triplet'])
    
    # B. Trigger Filter (Using the string prepared above)
    df_trig = df_triplet.Filter(trig_filter_string)

    # C. Preselection & Anti-D0
    df_presel = df_trig.Filter(config['preselection']['bdt_score'])
    df_antid0 = df_presel.Filter(config['preselection']['anti_d0'])

    # -----------------------------------------------------------------------
    # STEP 3: REGION LOOP (Branching Paths)
    # -----------------------------------------------------------------------
    file_results = []

    # Calculate common counts once
    count_before = df.Sum('FONLLweight').GetValue()
    count_triplet = df_triplet.Sum('FONLLweight').GetValue()
    count_trig = df_trig.Sum('FONLLweight').GetValue()
    sow_trig = df_trig.Sum('total_weight').GetValue()
    sow_presel = df_presel.Sum('total_weight').GetValue()
    sow_antid0 = df_antid0.Sum('total_weight').GetValue()

    for region in target_regions:
        # Determine output path
        region_dir = os.path.join(base_output_dir, region)
        os.makedirs(region_dir, exist_ok=True)
        
        mode_suffix = f"_{args.target_trigger}" if args.mode == 'single' else "_mix"
        out_name = os.path.basename(input_path).replace('.root', f'_skimmed{mode_suffix}.root')
        out_path = os.path.join(region_dir, out_name)

        # Apply Region-Specific Q2 Cut
        q2_cut_str = config['q2_cuts'].get(region, "1")
        df_final = df_antid0.Filter(q2_cut_str)
        
        # Count Final
        sow_final = df_final.Sum('total_weight').GetValue()

        # Log
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

    ROOT.ROOT.EnableImplicitMT()

    with open(args.config) as f: config = yaml.safe_load(f)
    with open(args.samples) as f: sample_list = yaml.safe_load(f)

    if args.mode == 'single' and not args.target_trigger:
        raise ValueError("Single mode requires --target_trigger")

    all_results = []
    for entry in sample_list['samples']:
        res_list = process_file(entry, config, args.output, args)
        all_results.extend(res_list)

    # Log Saving
    os.makedirs('../data/logs', exist_ok=True)
    log_name = f'cutflow_step1_{args.mode}.csv'
    if args.mode == 'single':
        log_name = f'cutflow_step1_{args.target_trigger}.csv' 
    else:
        log_name = f'cutflow_step1_mix.csv'
    
    pd.DataFrame(all_results).to_csv(os.path.join('../data/logs', log_name), index=False)
