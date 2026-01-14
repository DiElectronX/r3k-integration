Here is the complete `README.md` file content. You can copy and paste the code block below directly into your file.

```markdown
# R(K) Analysis Integration Framework

This repository contains the analysis integration tools for the B meson precision measurement. It includes modules for generating cutflow/efficiency tables and performing rapid physics cross-checks (R(K), Double Ratios) against Standard Model expectations.

## 📦 Setup & Dependencies

Ensure you have a Python environment (e.g., LCG View or Conda) with the following libraries:

```bash
pip install pandas numpy matplotlib seaborn uncertainties pyyaml

```

## 📂 Repository Structure

```text
.
├── config/
│   └── xcheck_params.yaml      # Configuration for yields, effs, and PDG constants
├── data/
│   └── logs/                   # Input CSV logs from Step 1/Step 2 processing
├── notebooks/
│   └── tables_maker.ipynb      # Main script for Acc x Eff tables and plots
├── scripts/
│   └── quick_xchecks.py        # CLI script for R(K) ratios and tension calculations
└── README.md

```

---

## 🚀 Module 1: Cutflow & Efficiency Tables

**File:** `notebooks/tables_maker.ipynb`

This Jupyter notebook processes the raw log files from the analysis pipeline (Step 1 & Step 2) to produce publication-ready tables.

### **Key Features**

* **Physics Normalization:** Injects `N_miniaod` and `FONLL` weights to calculate absolute acceptance.
* **Detector Acceptance Injection:** Manually injects pre-calculated  values (Step 0) for signal and control channels.
* **Uncertainty Propagation:** Uses the `uncertainties` library to rigorously propagate errors through efficiencies and ratios.
* **Duplicate Handling:** Robustly filters duplicate sample entries across analysis regions (e.g., J/ leakage into Low-).

### **Usage**

1. Open the notebook:
```bash
jupyter notebook notebooks/tables_maker.ipynb

```


2. Configure the `TRIGGER_TAG` (e.g., `'mix'`, `'L1_8p0_HLT_5p0'`) in **Cell 1**.
3. Run all cells.
4. **Outputs:**
* Cutflow Plots (Log scale).
* Detailed Efficiency Table (Step-by-step relative efficiency).
* Final **Acceptance  Efficiency** Table (Formatted for Note/Paper).



---

## 🔬 Module 2: Physics Cross-Checks

**File:** `scripts/quick_xchecks.py`

A command-line tool to calculate  ratios, double ratios, and statistical pulls (tension with SM) using yields and efficiencies defined in a YAML config.

### **Key Features**

* **Era Comparisons:** Checks consistency between Run 2 (2018) and Run 3 (2022).
* **Double Ratios:** Calculates  to cancel systematic uncertainties.
* **Config Driven:** All inputs (Yields, Lumi, Efficiencies) are separated from code in `config/xcheck_params.yaml`.

### **Configuration**

Edit `config/xcheck_params.yaml` to update your latest numbers:

```yaml
# --- Physics Constants (PDG) ---
pdg:
  br_b_to_jpsik:    {value: 1.020e-3, error: 0.019e-3}
  # ... (other constants)

# --- Analysis Data ---
channels:
  mumu:
    2022:
      yields:
        n_jpsi: {value: 3921547, error: 1641}
      effs:
        eff_jpsi: {value: 0.00453, error: 0.00001}

```

### **Usage**

Run the script from the root directory:

```bash
python scripts/quick_xchecks.py -c config/xcheck_params.yaml

```

**Example Output:**

```text
================================================================================
CROSS-CHECK RESULTS
================================================================================

--- 2022 Era (Muon + Electron) ---
R(J/ψ) [Control]                     : 1.0312+/-0.0450      (Pull from 1.0: 0.7σ)
Double Ratio                         : 0.9850+/-0.0620      (Pull from 1.0: 0.2σ)

--- Mixed Era (2018 Muon + 2022 Electron) ---
R(J/ψ) [Control]                     : 1.0120+/-0.0410      (Pull from 1.0: 0.3σ)

--- Internal Ratios (ψ(2S) / J/ψ) ---
Muons (2022)                         : 0.0836+/-0.0003      (PDG: 0.0821+/-0.0069, Pull: 0.2σ)
================================================================================

```

```

```
