import numpy as np
from uncertainties import ufloat
import yaml
import argparse
import sys

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def get_ufloat(d):
    """Parses a dict with 'value' and 'error' keys into a ufloat."""
    return ufloat(float(d['value']), float(d['error']))

def get_sigma_unc(val1, val2, rounded=1):
    """Calculates the compatibility (in sigmas) between two values."""
    diff = abs(val1.n - val2.n)
    unc = np.sqrt(val1.s**2 + val2.s**2)
    return round(diff / unc, rounded)

def main():
    parser = argparse.ArgumentParser(description="Quick R(K) cross-checks using yields and efficiencies.")
    parser.add_argument("-c", "--config", type=str, required=True, help="Path to YAML config file")
    args = parser.parse_args()

    # 1. Load Configuration
    try:
        with open(args.config, 'r') as f:
            cfg = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: Config file '{args.config}' not found.")
        sys.exit(1)

    # 2. Parse PDG Constants (Standard Model Baselines)
    try:
        pdg_vals = dotdict({k: get_ufloat(v) for k, v in cfg['pdg'].items()})
    except KeyError:
        print("Error: 'pdg' section missing from configuration file.")
        sys.exit(1)

    # 3. Parse Analysis Data into Structured Data
    try:
        channels = dotdict({
            ch: dotdict({
                # FIX: Force str(year) here to handle YAML integer parsing
                str(year): dotdict({
                    'lumi': get_ufloat(cfg['channels'][ch][year]['lumi']),
                    'yields': dotdict({k: get_ufloat(v) for k, v in cfg['channels'][ch][year].get('yields', {}).items()}),
                    'effs': dotdict({k: get_ufloat(v) for k, v in cfg['channels'][ch][year].get('effs', {}).items()}),
                }) for year in cfg['channels'][ch]
            }) for ch in cfg['channels']
        })
    except KeyError as e:
        print(f"Error parsing channels in config file. Missing key: {e}")
        sys.exit(1)

    # 4. Assign Variables for Calculation
    # Note: We use 'mu18' and 'mu22' to compare run periods, and 'ee22' for the ratio.
    mu18 = channels.mumu['2018']
    mu22 = channels.mumu['2022']
    ee22 = channels.ee['2022']

    # --- Calculations ---

    # R(K) J/psi (2022)
    # Ratio of (Yield/Lumi)_mu / (Yield/Lumi)_e * (Eff_e / Eff_mu)
    _rk_jpsi_2022 = (
        (mu22.yields.n_jpsi / mu22.lumi) / (ee22.yields.n_jpsi / ee22.lumi) *
        (ee22.effs.eff_jpsi / mu22.effs.eff_jpsi)
    )

    # R(K) Psi(2S) (2022)
    _rk_psi2s_2022 = (
        (mu22.yields.n_psi2s / mu22.lumi) / (ee22.yields.n_psi2s / ee22.lumi) *
        (ee22.effs.eff_psi2s / mu22.effs.eff_psi2s)
    )

    # Double Ratio 2022: R(K) Psi(2S) / R(K) J/psi
    # This cancels out luminosity and many systematic uncertainties
    _rk_psi2s_jpsi_2022 = _rk_psi2s_2022 / _rk_jpsi_2022

    # Mixed Era Checks (2018 Muons vs 2022 Electrons)
    # Useful for stability checks if 2022 muon data is suspect
    _rk_jpsi_2018 = (
        (mu18.yields.n_jpsi / mu18.lumi) / (ee22.yields.n_jpsi / ee22.lumi) *
        (ee22.effs.eff_jpsi / mu18.effs.eff_jpsi)
    )
    _rk_psi2s_2018 = (
        (mu18.yields.n_psi2s / mu18.lumi) / (ee22.yields.n_psi2s / ee22.lumi) *
        (ee22.effs.eff_psi2s / mu18.effs.eff_psi2s)
    )
    _rk_psi2s_jpsi_2018 = _rk_psi2s_2018 / _rk_jpsi_2018

    # Branching Fraction Ratios (Data vs PDG)
    _ratio_psi2s_jpsi_mu_2022 = (mu22.yields.n_psi2s / mu22.yields.n_jpsi) * (mu22.effs.eff_jpsi / mu22.effs.eff_psi2s)
    _ratio_psi2s_jpsi_el_2022 = (ee22.yields.n_psi2s / ee22.yields.n_jpsi) * (ee22.effs.eff_jpsi / ee22.effs.eff_psi2s)
    
    # PDG Expectations
    _pdg_ratio_mu = (pdg_vals.br_b_to_psi2sk * pdg_vals.br_psi2s_to_mumu) / (pdg_vals.br_b_to_jpsik * pdg_vals.br_jpsi_to_mumu)
    _pdg_ratio_ee = (pdg_vals.br_b_to_psi2sk * pdg_vals.br_psi2s_to_ee)   / (pdg_vals.br_b_to_jpsik * pdg_vals.br_jpsi_to_ee)

    # --- Output ---
    # Width configurations
    L_WIDTH = 38 
    V_WIDTH = 25

    l_rk_jpsi   = "R(J/\u03C8) [Control]"
    l_rk_psi2s  = "R(\u03C8(2S)) [Signal]"
    l_double    = "Double Ratio"
    
    l_mu_ratio  = "Muons (2022)"
    l_el_ratio  = "Electrons (2022)"

    print(f"{'='*80}")
    print(f"CROSS-CHECK RESULTS")
    print(f"{'='*80}")
    
    print('\n--- 2022 Era (Muon + Electron) ---')
    print(f"{l_rk_jpsi:<{L_WIDTH}} : {str(_rk_jpsi_2022):<{V_WIDTH}} (Pull from 1.0: {get_sigma_unc(ufloat(1,0), _rk_jpsi_2022)}σ)")
    print(f"{l_rk_psi2s:<{L_WIDTH}} : {str(_rk_psi2s_2022):<{V_WIDTH}} (Pull from 1.0: {get_sigma_unc(ufloat(1,0), _rk_psi2s_2022)}σ)")
    print(f"{l_double:<{L_WIDTH}} : {str(_rk_psi2s_jpsi_2022):<{V_WIDTH}} (Pull from 1.0: {get_sigma_unc(ufloat(1,0), _rk_psi2s_jpsi_2022)}σ)")

    print('\n--- Mixed Era (2018 Muon + 2022 Electron) ---')
    print(f"{l_rk_jpsi:<{L_WIDTH}} : {str(_rk_jpsi_2018):<{V_WIDTH}} (Pull from 1.0: {get_sigma_unc(ufloat(1,0), _rk_jpsi_2018)}σ)")
    print(f"{l_double:<{L_WIDTH}} : {str(_rk_psi2s_jpsi_2018):<{V_WIDTH}} (Pull from 1.0: {get_sigma_unc(ufloat(1,0), _rk_psi2s_jpsi_2018)}σ)")

    print('\n--- Internal Ratios (\u03C8(2S) / J/\u03C8) ---')
    print(f"{l_mu_ratio:<{L_WIDTH}} : {str(_ratio_psi2s_jpsi_mu_2022):<{V_WIDTH}} (PDG: {str(_pdg_ratio_mu)}, Pull: {get_sigma_unc(_pdg_ratio_mu, _ratio_psi2s_jpsi_mu_2022)}σ)")
    print(f"{l_el_ratio:<{L_WIDTH}} : {str(_ratio_psi2s_jpsi_el_2022):<{V_WIDTH}} (PDG: {str(_pdg_ratio_ee)}, Pull: {get_sigma_unc(_pdg_ratio_ee, _ratio_psi2s_jpsi_el_2022)}σ)")
    print(f"{'='*80}")

if __name__ == "__main__":
    main()
