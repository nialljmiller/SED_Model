"""
extinction_prescription_comparison.py
======================================
Demonstrates how different extinction laws affect the inverse fit.

Scenario
--------
We have a star at 500 pc with genuine Av=0.8 reddening.  We generate
synthetic observed magnitudes from the forward model using Fitzpatrick99
(the "truth"), then re-fit it six times — once per prescription — and
compare how well each one recovers the true Teff, logg, and [M/H].

This is also useful as a sensitivity test: if all six prescriptions
agree on the posterior, the photometry data is not constraining Av much
and the choice of law doesn't matter.  If they disagree, you need to
think about which is appropriate for your sight line.

Usage
-----
    python extinction_prescription_comparison.py

Requirements: sed_model built (make), emcee.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from sed_model import (
    load_grid, load_filters,
    run_forward, run_inverse,
    ExtinctionModel, make_extinction_model, AVAILABLE_LAWS,
    InverseResult,
)

import os

GRID_PATH    = os.path.expanduser("~/SED_Tools/data/stellar_models/Kurucz2003all")
FILTER_DIR   = os.path.expanduser("~/SED_Tools/data/filters/GAIA/GAIA")
VEGA_SED     = os.path.expanduser("~/SED_Tools/data/stellar_models/vega_flam.csv")

FILTER_PATHS = [
    os.path.join(FILTER_DIR, "G.dat"),
    os.path.join(FILTER_DIR, "Gbp.dat"),
    os.path.join(FILTER_DIR, "Grp.dat"),
]


# "True" stellar parameters used to generate synthetic photometry
TRUE_TEFF  = 5800.0    # K
TRUE_LOGG  = 4.40
TRUE_META  = 0.02


# Physical setup
R_SUN_CM   = 6.957e10              # cm
PC_TO_CM   = 3.085677581491367e18  # cm per parsec
TRUE_R     = 1.0 * R_SUN_CM       # 1 R_sun
TRUE_D     = 500.0 * PC_TO_CM     # 500 pc


# True extinction used to generate the synthetic observations
TRUE_AV    = 0.8
TRUE_RV    = 3.1
TRUE_LAW   = "gordon23"


# Noise added to synthetic magnitudes
OBS_SIGMA  = 0.02   # mag (1-sigma per band)

# MCMC settings — keep short for a quick demo; increase for publication
N_WALKERS  = 32
N_STEPS    = 1500
N_BURN     = 400
N_THIN     = 1


MAG_SYSTEM = "Vega"

ext_kwargs = dict(enabled=True, law=TRUE_LAW,      a_v=TRUE_AV, gordon23_env="mw")

# ---------------------------------------------------------------------------
# Helper: summary stats from an InverseResult
# ---------------------------------------------------------------------------

def posterior_stats(result: InverseResult):
    """Return (median, lo_1sigma, hi_1sigma) for each of Teff, logg, meta."""
    s = result.samples
    med  = np.percentile(s, 50,   axis=0)
    lo   = np.percentile(s, 15.87, axis=0)
    hi   = np.percentile(s, 84.13, axis=0)
    return med, lo, hi




def main():
    # -----------------------------------------------------------------------
    # Load grid and filters once
    # -----------------------------------------------------------------------
    print("Loading atmosphere grid...")
    grid = load_grid(GRID_PATH)
    print(f"  Grid: Teff {grid.teff_bounds}, logg {grid.logg_bounds}, "
          f"[M/H] {grid.meta_bounds}")

    print("Loading filters...")
    vega_path = VEGA_SED if os.path.isfile(VEGA_SED) else None
    filters = load_filters(FILTER_PATHS, vega_sed_path=vega_path)
    filter_names = [f.name for f in filters]
    print(f"  Filters: {filter_names}")

    # -----------------------------------------------------------------------
    # Generate synthetic observations using the true law (Fitzpatrick99)
    # -----------------------------------------------------------------------
    print(f"\nGenerating synthetic observations:")
    print(f"  Teff={TRUE_TEFF} K  logg={TRUE_LOGG}  [M/H]={TRUE_META}")
    print(f"  R={TRUE_R:.3e} cm  d={TRUE_D:.3e} cm  ({TRUE_D/PC_TO_CM:.0f} pc)")
    print(f"  True extinction: law={TRUE_LAW}  Av={TRUE_AV}  Rv={TRUE_RV}")

    true_ext = make_extinction_model(
        enabled=True, law=TRUE_LAW, a_v=TRUE_AV, r_v=TRUE_RV
    )
    truth = run_forward(
        teff=TRUE_TEFF, logg=TRUE_LOGG, meta=TRUE_META,
        R=TRUE_R, d=TRUE_D,
        grid=grid, filters=filters,
        mag_system=MAG_SYSTEM,
        extinction=true_ext,
    )

    rng = np.random.default_rng(42)
    obs_mags  = np.array([truth.magnitudes[n] for n in filter_names])
    obs_mags += rng.normal(0.0, OBS_SIGMA, size=len(filter_names))
    obs_errs  = np.full(len(filter_names), OBS_SIGMA)

    print("\n  Band         truth_mag   obs_mag")
    for name, tm, om in zip(filter_names, [truth.magnitudes[n] for n in filter_names], obs_mags):
        print(f"  {name:12s}  {tm:8.4f}    {om:8.4f}")

    # -----------------------------------------------------------------------
    # Run one inverse fit per prescription
    # -----------------------------------------------------------------------
    results = {}

    print(f"\n{'='*60}")
    print(f"{'='*60}")

    ext_model = make_extinction_model(**ext_kwargs)

    posterior = run_inverse(
        obs_magnitudes=obs_mags,
        obs_uncertainties=obs_errs,
        filter_names=filter_names,
        R=TRUE_R,
        d=TRUE_D,
        grid=grid,
        filters=filters,
        mag_system=MAG_SYSTEM,
        extinction=ext_model,
        n_walkers=N_WALKERS,
        n_steps=N_STEPS,
        n_burn=N_BURN,
        n_thin=N_THIN,
        p0_teff=TRUE_TEFF,
        p0_logg=TRUE_LOGG,
        p0_meta=TRUE_META,
        p0_scatter=0.02,
        progress=True,
    )

    med, lo, hi = posterior_stats(posterior)


    print(f"\n  Results:")
    print(f"  Teff  = {med[0]:.1f}  [{lo[0]:.1f}, {hi[0]:.1f}]  (true: {TRUE_TEFF})")
    print(f"  logg  = {med[1]:.3f}  [{lo[1]:.3f}, {hi[1]:.3f}]  (true: {TRUE_LOGG})")
    print(f"  [M/H] = {med[2]:.3f}  [{lo[2]:.3f}, {hi[2]:.3f}]  (true: {TRUE_META})")

# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
