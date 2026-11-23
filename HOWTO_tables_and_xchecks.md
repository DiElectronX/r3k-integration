# Efficiency tables and cross-checks

All commands from the r3k-integration directory. Environment for the scripts:

    source /cvmfs/sft.cern.ch/lcg/views/LCG_106_swan/x86_64-el9-gcc13-opt/setup.sh

(the notebook runner sources it itself; the r3k-int conda kernel also works interactively)

## 1. Preselection (step 1)

    python3 scripts/01_apply_preselection.py --config config/cuts.yml --samples config/01_samples.yml --output <output_loc>
    python3 combine_data.py data/logs/cutflow_step1_mix.csv     # hadd data + same-sign data counts

Trigger mixture fractions and trigger scale factors live in the cuts file. `config/cuts_passfail_v2.yml` holds the
current pass/fail-fit trigger SFs with the exclusive fractions including run 362439;
`config/01_samples_data_withG4.yml` is the 2022 data list including era G ParkingDoubleElectronLowMass4.

## 2. BDT scores

Run `bdt_inference.py` in r3k-bdttools on the step-1 skims (see that repo), stage the wScores files next to the skims.

## 3. BDT selection (step 2)

    python3 scripts/02_apply_bdt_selection.py --config config/cuts.yml --samples config/02_samples.yml --input_dir <wScores_loc>

## 4. Efficiency table

    notebooks/run_tables_mix.sh            # executes notebooks/tables_maker.ipynb in place (TRIGGER_TAG='mix')

Inputs: `data/logs/cutflow_step{1,2}_mix.csv`. The final acceptance x efficiency for J/psi and psi(2S) goes into the
cross-check config.

## 5. Cross-checks

Enter the do_fit signal yields (r3k-fitter) and the table efficiencies in `config/xcheck_params.yml`, then

    python3 scripts/rk_xchecks_v2.py -c config/xcheck_params.yml

Note: any custom fitter input file must span the full B-mass observable range 4.5-5.7 GeV (do_fit hard-codes the
domain; the config fit_range is only the window).
