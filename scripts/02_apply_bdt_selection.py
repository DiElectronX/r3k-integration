import ROOT
import os
import yaml
import pandas as pd
import argparse

def find_scored_file(input_base_dir, region, sample_name, core_filename, suffix, trig_tag):
    """
    Intelligently hunts for the scored file in various subdirectory structures
    and matching variable intermediate suffixes.
    """
    # 1. Define Candidate Directories
    # We check these locations in order of likelihood
    candidate_dirs = [
        os.path.join(input_base_dir, region),                     # Standard: /jpsi/
        os.path.join(input_base_dir, f"{region}_wScores"),        # Dir Suffix: /jpsi_wScores/
        os.path.join(input_base_dir, sample_name, region),        # Sample folder: /Bd_JpsiK/jpsi/
        os.path.join(input_base_dir, sample_name, f"{region}_wScores") # Combo: /Bd_JpsiK/jpsi_wScores/
    ]

    found_file = None
    
    for d in candidate_dirs:
        if not os.path.exists(d):
            continue
            
        # List all files in this candidate directory
        files = os.listdir(d)
        
        # 2. Filter files
        # Criteria:
        # a) Starts with the Core Name (e.g. "Bu_kaon_jpsi_resonant")
        # b) Ends with the BDT suffix (e.g. "_wScores.root")
        # c) Contains the trigger tag (e.g. "mix" or "L1_8p0") to avoid mixing modes
        
        matches = [
            f for f in files 
            if f.startswith(core_filename) 
            and f.endswith(suffix + '.root')
            and trig_tag in f
        ]
        
        if len(matches) == 1:
            found_file = os.path.join(d, matches[0])
            break # Found it!
        elif len(matches) > 1:
            print(f"  [Warning] Ambiguous match in {d}: {matches}. Using {matches[0]}")
            found_file = os.path.join(d, matches[0])
            break

    return found_file

def process_step2(file_info, config, input_base_dir, args):
    sample_name = file_info['name']
    target_regions = file_info.get('regions', ['none'])
    
    # Read Global Config
    bdt_branch = config['final_selection']['bdt_branch']
    global_cut = config['final_selection']['cut_value']
    
    results = []

    print(f"\nProcessing {sample_name}...")
    
    # Extract "Core" filename from the raw input path
    # e.g. "/eos/.../Bu_kaon_jpsi_resonant.root" -> "Bu_kaon_jpsi_resonant"
    raw_path = file_info['path']
    core_filename = os.path.splitext(os.path.basename(raw_path))[0]

    for region in target_regions:
        
        # USE SMART FINDER
        input_path = find_scored_file(
            input_base_dir, 
            region, 
            sample_name, 
            core_filename, 
            args.bdt_suffix,
            args.trigger_tag
        )
        
        if not input_path:
            print(f"  [Skipping] Could not find scored file for region '{region}' matching core '{core_filename}'")
            continue
            
        # 2. Apply Selection
        df = ROOT.RDataFrame("mytree", input_path)
        
        cols = [str(c) for c in df.GetColumnNames()]
        if bdt_branch not in cols:
            print(f"  [Error] Branch '{bdt_branch}' not found in {input_path}!")
            continue

        # Apply Global Cut
        cut_expression = f"{bdt_branch} > {global_cut}"
        df_final = df.Filter(cut_expression)
        
        # 3. Count
        if 'total_weight' not in cols:
             # Fallback if BDT tool dropped the weight column but kept relative events
             # This is risky, but allows code to run if weights are missing.
             # Better to warn.
             print(f"  [Error] 'total_weight' missing in {input_path}. Counts will be raw!")
             sow_before = df.Count().GetValue()
             sow_final = df_final.Count().GetValue()
        else:
             sow_before = df.Sum('total_weight').GetValue()
             sow_final = df_final.Sum('total_weight').GetValue()
        
        # Efficiency Check
        eff = 0.0
        if sow_before > 0:
            eff = sow_final / sow_before

        results.append({
            'Sample Name': sample_name,
            'Region': region,
            'Trigger_Mode': args.trigger_mode,
            'BDT_Cut_Value': global_cut,
            'SoW_Pre_BDT': sow_before,     
            'SoW_Post_BDT': sow_final
        })
        
        print(f"  -> Region: {region} | Found: {os.path.basename(input_path)} | Eff: {eff:.2%}")

    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 2: Apply Final BDT Selection")
    parser.add_argument('--config', default='../config/cuts.yml')
    parser.add_argument('--samples', default='../config/samples.yml')
    
    parser.add_argument('--input_dir', required=True, help="Base directory containing scored files")
    
    parser.add_argument('--trigger_mode', default='mix', help="mix or single")
    parser.add_argument('--target_trigger', default=None, help="Required if mode=single")
    
    # Default matches your request
    parser.add_argument('--bdt_suffix', default='_wScores', help="Suffix to match (without .root)")
    
    args = parser.parse_args()

    ROOT.ROOT.EnableImplicitMT()

    if args.trigger_mode == 'single':
        if not args.target_trigger: raise ValueError("Must specify --target_trigger for single mode")
        # Trigger tag used for fuzzy matching
        args.trigger_tag = args.target_trigger 
    else:
        args.trigger_tag = "mix"

    with open(args.config) as f: config = yaml.safe_load(f)
    with open(args.samples) as f: sample_list = yaml.safe_load(f)

    all_results = []
    for entry in sample_list['samples']:
        res = process_step2(entry, config, args.input_dir, args)
        all_results.extend(res)

    os.makedirs('../data/logs', exist_ok=True)
    log_name = f'cutflow_step2_{args.trigger_tag}.csv'
    pd.DataFrame(all_results).to_csv(os.path.join('../data/logs', log_name), index=False)
    print(f"\nStep 2 Log saved to ../data/logs/{log_name}")
