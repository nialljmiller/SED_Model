# Extinction Prescription Comparison

`demos/extinction_prescription_comparison.py` generates synthetic Gaia-band observations of a star at 500 pc reddened with Fitzpatrick (1999) at Av = 0.8, then re-fits them once per prescription — no extinction, CCM89, O'Donnell94, Fitzpatrick99, FM07, Calzetti00, the three Gordon+2023 environments (MW/LMC/SMC), and Fitzpatrick99/CCM89 with non-standard R_V — and compares the recovered `Teff`, `logg`, `[M/H]` against truth.

This doubles as a sensitivity test: if all prescriptions agree, your photometry isn't constraining Av and the choice of law doesn't matter; if they disagree, the choice needs justifying for your sight line.

The output figure `extinction_comparison.png` shows median ± 1σ per prescription for each parameter, with the truth as a dashed line; per-parameter summary tables (including the bias of each prescription) are printed to the terminal.

Concepts demonstrated: [the extinction module](../guide/extinction.md), `make_extinction_model`, fixed-Av fitting under a deliberately wrong law.

## Source

```python
--8<-- "demos/extinction_prescription_comparison.py"
```
