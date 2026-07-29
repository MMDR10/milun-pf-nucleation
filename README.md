# Milun Fault P_f Nucleation

**High-Frequency (200 Hz) Pore Fluid Pressure Phase-Space Singularity Prior to the 2024 M7.4 Milun Fault Rupture**

*Authors: DR (tygtDc), MM (nnRpMr) | License: CC BY 4.0*

## Abstract

This proof-of-concept case study applies the Ô-HAT phase-space topological framework to a 200 Hz FBG pore fluid pressure (P_f) record from four depths (170–464 m) intersecting the Milun Fault zone, Taiwan, spanning 28 minutes before and 2 hours after the 2024-04-03 M7.4 Hualien earthquake (M7.4).

**Key findings:**

| Finding | Value | Significance |
|---------|-------|-------------|
| P_f decline (P170) | −0.55 MPa/hr, R²=0.84 | Begins ~19 min pre-rupture |
| θ-space divergence | P424×P464: +72.5% (p=0.002) | Deep sensors decouple — fluid migration across ~440m gas zone |
| EffRank singularity | 2.977, θ₁=6.3° at T−27s | Phase-space critical transition before rupture |
| Coseismic asymmetry | 37× (23.74 vs 0.64 bar) | Rupture-induced pressure imbalance |
| Control days | 0% false positives | Specific to EQ day (preliminary, N=5) |

## Repository Structure

```
├── papers/
│   └── milun-pf-nucleation.md          # Preprint manuscript (full paper)
├── scripts/
│   ├── deep_dive.py                     # Main analysis: P_f decline, θ divergence, gradient map
│   ├── ohat_200hz_analysis.py           # Ô-HAT phase-space analysis (initial)
│   ├── ohat_200hz_detail.py             # Detailed phase-space analysis
│   ├── ohat_200hz_dirA.py               # Directory A analysis
│   ├── ohat_200hz_A2_deep.py            # A2 deep analysis (pair-wise θ, EffRank)
│   ├── ohat_200hz_A2b_lag.py            # Lag analysis
│   ├── run_analysis.py                  # Analysis runner
│   ├── run_theta.py                     # θ-space computation
│   ├── summary.py                       # Results summary
│   ├── generate_report.py               # Report generator
│   └── analyze_pf.py                    # P_f data inspection
├── tests/
│   ├── test_pf_decline.py               # Validate P170 slope (−0.55 MPa/hr)
│   ├── test_phase_space.py              # Validate θ divergence
│   └── test_sample_pipeline.py          # End-to-end pipeline test
├── data/
│   ├── sample_data/                     # Sample CSV subsets for quick testing
│   │   ├── pressure_60s.csv             # First 60s of P_f (12,000 rows)
│   │   ├── pressure_pre_eq_5min.csv     # Last 5 min before EQ (60,000 rows)
│   │   └── accel_30s.csv               # First 30s of accel data (6,000 rows)
│   ├── results/                         # Analysis output JSONs
│   │   ├── ohat_200hz_results.json      # Phase Space & H values
│   │   ├── A2c_full_matrix_results.json # Full covariance matrix
│   │   ├── A2d_late_trend_results.json  # Late trend analysis
│   │   ├── A2e_pair_ohat_results.json   # Pair-wise Ô-HAT
│   │   ├── A2f_phase_transition_results.json
│   │   ├── A3_3sensor_ohat_results.json
│   │   ├── A3b_effrank_timing_*.json    # EffRank timing at 5s/15s/30s
│   │   └── A4_4sensor_ohat_results.json
│   └── mendeley_dataset.zip             # Original dataset download
├── docs/
│   ├── MILUN_PF_SUMMARY_REPORT.md       # Summary report (extended)
│   └── OHAT_200Hz_FINAL_REPORT.md       # Final analysis report
├── figures/
│   ├── phase_trajectory.png             # θ-space trajectory plot
│   └── dashboard_dual.png               # Dual dashboard visualization
└── README.md
```

## How to Reproduce

### 1. Get the Data
Download the original dataset from Mendeley:
```
Marginingsih et al. (2026), doi: 10.17632/dn4338tv9g.2
```
Extract to `data/` or modify the file paths in the scripts.

### 2. Quick Test with Sample Data
```bash
# Run pipeline test (uses sample data in data/sample_data/)
python3 tests/test_sample_pipeline.py

# Validate P_f decline
python3 tests/test_pf_decline.py

# Validate phase space divergence
python3 tests/test_phase_space.py
```

### 3. Full Analysis
```bash
python3 scripts/deep_dive.py
python3 scripts/ohat_200hz_A2_deep.py
python3 scripts/ohat_200hz_dirA.py
```

## Dataset Citation

> Marginingsih, A. S., Wang, S.-J., Ma, K.-F., & Lin, Y.-Y. (2026). The 200 Hz Sampling Interval of Porewater Pressure from The Integrated Groundwater Observation System of Milun Fault, Taiwan [Data set]. Mendeley Data, V2. https://doi.org/10.17632/dn4338tv9g.2

## Zenodo

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21682478.svg)](https://doi.org/10.5281/zenodo.21682478)

## Related

- **Ô-HAT Framework:** [MMDR10/ohat-framework](https://github.com/MMDR10/ohat-framework)
- **ENSO Watch:** [MMDR10/enso-watch](https://github.com/MMDR10/enso-watch)
- **Typhoon dH_curl:** [MMDR10/typhoon-dh-curl](https://github.com/MMDR10/typhoon-dh-curl)
