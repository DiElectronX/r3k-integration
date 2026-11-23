#!/bin/bash
# One-command run of the MIXTURE efficiency table (tables_maker.ipynb).
# Executes the notebook in place (open it afterwards to see the rendered tables) and
# also writes an HTML snapshot next to it.
# The notebook reads:
#   data/logs/cutflow_step1_L1_9p0_HLT_6p0.csv                (pipeline, unchanged)
#   data/logs/cutflow_step2_L1_9p0_HLT_6p0.csv                (CORRECTED: analysis-BDT
#        scores at WP 0.998; original 0.999/wrong-model file kept as *_run5model_0p999.bak.csv)
#   meeting_plots/npz/single_path_result.json                 (fit results, last cell)
set -e
cd "$(dirname "$0")"
source /cvmfs/sft.cern.ch/lcg/views/LCG_106_swan/x86_64-el9-gcc13-opt/setup.sh
python3 -m nbconvert --to notebook --execute --inplace tables_maker.ipynb --ExecutePreprocessor.timeout=900
python3 -m nbconvert --to html tables_maker.ipynb
echo "DONE: open tables_maker.ipynb (or .html) for the mixture table"
