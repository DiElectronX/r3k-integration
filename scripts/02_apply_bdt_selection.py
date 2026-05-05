import yaml
import ROOT
import os
import pandas as pd
import argparse
import numpy as np

# --------------------------------------------------------------------------------
# DEBUG HELPERS
# --------------------------------------------------------------------------------
DEBUG_SAMPLES = {"jpsi_resonant", "kaon_psi2s_resonant", "kaon_lowq_resonant"}


def should_debug(sample_name):
    return sample_name in DEBUG_SAMPLES


def n_unique_events(rdf, event_key_col):
    arr = rdf.AsNumpy([event_key_col])[event_key_col]
    return len(np.unique(arr))


def maybe_define_event_key(df, sample_name, region, input_path):
    cols = [str(c) for c in df.GetColumnNames()]

    if all(c in cols for c in ["run", "luminosityBlock", "event"]):
        df = df.Define(
            "event_key",
            '(static_cast<unsigned long long>(run) << 42) ^ '
            '(static_cast<unsigned long long>(luminosityBlock) << 20) ^ '
            'static_cast<unsigned long long>(event)'
        )
        return df, "event_key"

    if all(c in cols for c in ["run_", "luminosityBlock_", "event_"]):
        df = df.Define(
            "event_key",
            '(static_cast<unsigned long long>(run_) << 42) ^ '
            '(static_cast<unsigned long long>(luminosityBlock_) << 20) ^ '
            'static_cast<unsigned long long>(event_)'
        )
        return df, "event_key"

    if should_debug(sample_name):
        print(f"[DEBUG][{sample_name}][{region}] no run/lumi/event branches in {input_path}")
        print(f"[DEBUG][{sample_name}][{region}] unique-event debug disabled")
        print(f"[DEBUG][{sample_name}][{region}] available columns (first 40):")
        print([str(c) for c in cols[:40]])

    return df, None


def debug_bdt_stage(sample_name, region, label, rdf, bdt_branch, event_key_col=None):
    cols = [str(c) for c in rdf.GetColumnNames()]
    n_cand = rdf.Count().GetValue()

    msg = f"[DEBUG][{sample_name}][{region}][{label}] candidates={n_cand}"

    if event_key_col is not None and event_key_col in cols:
        n_evt = n_unique_events(rdf, event_key_col)
        msg += f", unique_events={n_evt}"
    else:
        msg += ", unique_events=NA"

    if "total_weight" in cols:
        msg += f", sum(total_weight)={rdf.Sum('total_weight').GetValue()}"

    print(msg)

    if bdt_branch in cols and n_cand > 0:
        print(
            f"[DEBUG][{sample_name}][{region}] "
            f"{bdt_branch} min={rdf.Min(bdt_branch).GetValue()}, "
            f"max={rdf.Max(bdt_branch).GetValue()}"
        )

        for thr in [0.998, 0.999, 0.9995, 0.9998, 0.9999, 0.99995, 0.99998, 0.99999]:
            df_thr = rdf.Filter(f"{bdt_branch} > {thr}")
            n_thr = df_thr.Count().GetValue()
            w_thr = df_thr.Sum("total_weight").GetValue() if "total_weight" in cols else float("nan")
            print(
                f"[DEBUG][{sample_name}][{region}] "
                f"{bdt_branch} > {thr}: candidates={n_thr}, sumw={w_thr}"
            )


# --------------------------------------------------------------------------------
# FILE FINDER
# --------------------------------------------------------------------------------
def find_scored_file(input_base_dir, region, sample_name, core_filename, suffix, trig_tag):
    """
    Intelligently hunts for the scored file in various subdirectory structures
    and matching variable intermediate suffixes.
    """
    candidate_dirs = [
        os.path.join(input_base_dir, region),
        os.path.join(input_base_dir, f"{region}_wScores"),
        os.path.join(input_base_dir, sample_name, region),
        os.path.join(input_base_dir, sample_name, f"{region}_wScores"),
    ]

    found_file = None

    for d in candidate_dirs:
        if not os.path.exists(d):
            continue

        files = os.listdir(d)

        matches = [
            f for f in files
            if f.startswith(core_filename)
            and f.endswith(suffix + ".root")
            and trig_tag in f
        ]

        if len(matches) == 1:
            found_file = os.path.join(d, matches[0])
            break
        elif len(matches) > 1:
            print(f"  [Warning] Ambiguous match in {d}: {matches}. Using {matches[0]}")
            found_file = os.path.join(d, matches[0])
            break

    return found_file


# --------------------------------------------------------------------------------
# MAIN STEP-2 PROCESSING
# --------------------------------------------------------------------------------
def process_step2(file_info, config, input_base_dir, args):
    sample_name = file_info["name"]
    target_regions = file_info.get("regions", ["none"])

    bdt_branch = config["final_selection"]["bdt_branch"]
    global_cut = config["final_selection"]["cut_value"]

    results = []

    print(f"\nProcessing {sample_name}...")

    raw_path = file_info["path"]
    core_filename = os.path.splitext(os.path.basename(raw_path))[0]

    for region in target_regions:
        input_path = find_scored_file(
            input_base_dir,
            region,
            sample_name,
            core_filename,
            args.bdt_suffix,
            args.trigger_tag,
        )

        print(f"    [DEBUG] region={region} using input_path={input_path}")

        if not input_path:
            print(f"  [Skipping] Could not find scored file for region '{region}' matching core '{core_filename}'")
            continue

        # Load the scored file
        df = ROOT.RDataFrame("Events", input_path)

        if should_debug(sample_name):
            tmp_cols = [str(c) for c in df.GetColumnNames()]
            print(f"[DEBUG][{sample_name}][{region}] n_columns = {len(tmp_cols)}")

        # Define event_key only if the scored file kept run/lumi/event
        df, event_key_col = maybe_define_event_key(df, sample_name, region, input_path)

        cols = [str(c) for c in df.GetColumnNames()]
        if bdt_branch not in cols:
            print(f"  [Error] Branch '{bdt_branch}' not found in {input_path}!")
            print(f"  [Error] Available columns (first 50): {cols[:50]}")
            continue

        cut_expression = f"{bdt_branch} > {global_cut}"

        if should_debug(sample_name):
            print(f"\n[DEBUG][{sample_name}][{region}] input_path = {input_path}")
            print(f"[DEBUG][{sample_name}][{region}] final BDT branch = {bdt_branch}")
            print(f"[DEBUG][{sample_name}][{region}] final BDT cut    = {cut_expression}")
            debug_bdt_stage(sample_name, region, "PRE_BDT", df, bdt_branch, event_key_col)

        df_final = df.Filter(cut_expression)

        if should_debug(sample_name):
            debug_bdt_stage(sample_name, region, "POST_BDT", df_final, bdt_branch, event_key_col)

        if "total_weight" not in cols:
            print(f"  [Error] 'total_weight' missing in {input_path}. Counts will be raw!")
            sow_before = df.Count().GetValue()
            sow_final = df_final.Count().GetValue()
        else:
            sow_before = df.Sum("total_weight").GetValue()
            sow_final = df_final.Sum("total_weight").GetValue()

        eff = 0.0
        if sow_before > 0:
            eff = sow_final / sow_before

        results.append({
            "Sample Name": sample_name,
            "Region": region,
            "Trigger_Mode": args.trigger_mode,
            "BDT_Cut_Value": global_cut,
            "SoW_Pre_BDT": sow_before,
            "SoW_Post_BDT": sow_final,
        })

        print(f"  -> Region: {region} | Found: {os.path.basename(input_path)} | Eff: {eff:.2%}")

    return results


# --------------------------------------------------------------------------------
# ENTRY POINT
# --------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 2: Apply Final BDT Selection")
    parser.add_argument("--config", default="../config/cuts.yml")
    parser.add_argument("--samples", default="../config/samples.yml")
    parser.add_argument("--input_dir", required=True, help="Base directory containing scored files")
    parser.add_argument("--trigger_mode", default="mix", help="mix or single")
    parser.add_argument("--target_trigger", default=None, help="Required if mode=single")
    parser.add_argument("--bdt_suffix", default="_wScores", help="Suffix to match (without .root)")

    args = parser.parse_args()

    ROOT.ROOT.EnableImplicitMT()

    if args.trigger_mode == "single":
        if not args.target_trigger:
            raise ValueError("Must specify --target_trigger for single mode")
        args.trigger_tag = args.target_trigger
    else:
        args.trigger_tag = "mix"

    with open(args.config) as f:
        config = yaml.safe_load(f)
    with open(args.samples) as f:
        sample_list = yaml.safe_load(f)

    all_results = []
    for entry in sample_list["samples"]:
        res = process_step2(entry, config, args.input_dir, args)
        all_results.extend(res)

    os.makedirs("../data/logs", exist_ok=True)
    log_name = f"cutflow_step2_{args.trigger_tag}.csv"
    pd.DataFrame(all_results).to_csv(os.path.join("../data/logs", log_name), index=False)
    print(f"\nStep 2 Log saved to ../data/logs/{log_name}")