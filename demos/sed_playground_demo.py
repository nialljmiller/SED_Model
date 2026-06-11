"""
sed_playground_demo.py
======================
One-file playground for SED_Model forward + inverse tests.

Change the two dictionaries below:

    TRUE_STAR   = what Nature/the fake observation really is
    FIT_CONFIG  = what the fitter is allowed to know or infer

For each fit parameter, use either:

    {"mode": "fixed", "value": ...}
    {"mode": "fit",   "bounds": (..., ...), "start": ...}

Distance is given in parsecs in the user config, but converted to cm
internally because SED_Model uses cgs.

Run:
    python sed_playground_demo.py
"""

from __future__ import annotations

import os
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt

from sed_model import (
    load_grid,
    load_filters,
    run_forward,
    run_inverse,
    make_extinction_model,
    FitParams,
    fixed,
    free,
    PC_TO_CM,
    RSUN_TO_CM,
)

# =============================================================================
# 1. Data paths: edit these once for your machine
# =============================================================================

GRID_PATH = Path(os.path.expanduser("~/SED_Tools/data/stellar_models/Kurucz2003all"))
FILTER_DIR = Path(os.path.expanduser("~/SED_Tools/data/filters/GAIA/GAIA"))
VEGA_SED = Path(os.path.expanduser("~/SED_Tools/data/stellar_models/vega_flam.csv"))

FILTER_PATHS = [
    FILTER_DIR / "G.dat",
    FILTER_DIR / "Gbp.dat",
    FILTER_DIR / "Grp.dat",
]

MAG_SYSTEM = "Vega"
INTERP_METHOD = "hermite"

# =============================================================================
# 2. Forward model playground: this is the fake real star
# =============================================================================

TRUE_STAR = {
    # Stellar atmosphere parameters
    "teff": 5800.0,       # K
    "logg": 4.40,         # log10(g / cm s^-2)
    "meta": 0.02,         # [M/H] dex

    # Flux scaling
    "radius_rsun": 1.0,   # stellar radius in solar radii
    "distance_pc": 500.0, # distance in parsecs

    # Dust extinction
    "a_v": 0.8,           # V-band extinction in magnitudes
    "r_v": 3.1,           # extinction-curve shape parameter
    "law": "gordon23",   # ccm89, odonnell94, fitzpatrick99, fm07, calzetti00, gordon23
    "gordon23_env": "mw",# mw, lmc, smc; only used by law='gordon23'
}

OBS_SIGMA = 0.02          # fake magnitude uncertainty per band
RANDOM_SEED = 42

# =============================================================================
# 3. Inverse model playground: decide what the fitter is allowed to fit
# =============================================================================

# Toggle any parameter between fixed and fit.
#
# Examples:
#   Fit Teff only:
#       teff = {"mode": "fit", "bounds": (5200, 6400), "start": 5800}
#       logg/meta/a_v/d_pc fixed
#
#   Fit Teff + Av:
#       teff = {"mode": "fit", "bounds": (5200, 6400), "start": 5800}
#       a_v  = {"mode": "fit", "bounds": (0.0, 1.5), "start": 0.5}
#
#   Fit everything:
#       set all five to mode='fit', but expect degeneracies, especially Av--distance.

FIT_CONFIG = {
    "teff": {"mode": "fit",   "bounds": (5200.0, 6400.0), "start": 5800.0},
    "logg": {"mode": "fixed", "value": TRUE_STAR["logg"]},
    "meta": {"mode": "fixed", "value": TRUE_STAR["meta"]},
    "a_v":  {"mode": "fit",   "bounds": (0.0, 1.5), "start": 0.5},
    "d_pc": {"mode": "fixed", "value": TRUE_STAR["distance_pc"]},
}

FIT_EXTINCTION = {
    "enabled": True,
    "law": "gordon23",
    "r_v": 3.1,
    "gordon23_env": "mw",
}

MCMC = {
    "n_walkers": 32,
    "n_steps": 1500,
    "n_burn": 400,
    "n_thin": 1,
    "p0_scatter": 0.02,
    "seed": 123,
    "progress": True,
}

# =============================================================================
# Helpers
# =============================================================================

DISPLAY_NAMES = {
    "teff": "Teff [K]",
    "logg": "logg",
    "meta": "[M/H]",
    "a_v": "Av [mag]",
    "d": "distance [pc]",
}


def _mode_to_spec(name: str, cfg: dict, grid=None):
    """Convert the user-friendly FIT_CONFIG entry into a ParamSpec."""
    mode = cfg["mode"].lower()
    if mode == "fixed":
        value = float(cfg["value"])
        if name == "d":
            value *= PC_TO_CM
        return fixed(value, name=name)

    if mode == "fit":
        lo, hi = cfg["bounds"]
        lo = float(lo)
        hi = float(hi)
        if name == "d":
            lo *= PC_TO_CM
            hi *= PC_TO_CM
        return free(lo, hi, name=name)

    raise ValueError(f"Unknown mode for {name!r}: {cfg['mode']!r}. Use fixed or fit.")


def build_fit_params(grid, fit_config: dict) -> FitParams:
    """Build FitParams from the nice human config block."""
    return FitParams(
        teff=_mode_to_spec("teff", fit_config["teff"], grid=grid),
        logg=_mode_to_spec("logg", fit_config["logg"], grid=grid),
        meta=_mode_to_spec("meta", fit_config["meta"], grid=grid),
        a_v=_mode_to_spec("a_v", fit_config["a_v"], grid=grid),
        d=_mode_to_spec("d", fit_config["d_pc"], grid=grid),
    )


def build_p0_centre(fit_config: dict) -> dict:
    """Starting values for only the parameters that are actually fitted."""
    centre = {}
    for user_name, model_name in [("teff", "teff"), ("logg", "logg"), ("meta", "meta"),
                                  ("a_v", "a_v"), ("d_pc", "d")]:
        cfg = fit_config[user_name]
        if cfg["mode"].lower() == "fit" and "start" in cfg:
            value = float(cfg["start"])
            if model_name == "d":
                value *= PC_TO_CM
            centre[model_name] = value
    return centre


def make_fit_extinction_model(fit_params: FitParams):
    """Build the extinction model used inside the inverse fit."""
    enabled = bool(FIT_EXTINCTION.get("enabled", False))

    # If Av is fixed, keep the ExtinctionModel and FitParams consistent.
    # If Av is free, run_forward will replace this placeholder value at each step.
    if fit_params.a_v.is_fixed:
        av_for_model = fit_params.a_v.value
    else:
        av_for_model = 0.0

    return make_extinction_model(
        enabled=enabled,
        law=FIT_EXTINCTION.get("law", "fitzpatrick99"),
        r_v=FIT_EXTINCTION.get("r_v", 3.1),
        a_v=av_for_model,
        gordon23_env=FIT_EXTINCTION.get("gordon23_env", "mw"),
    )


def print_observations(filter_names, truth_mags, obs_mags, obs_errs):
    print("\nSynthetic observations")
    print("  band          truth       observed     sigma")
    print("  ----       --------     --------     -----")
    for name, tm, om, oe in zip(filter_names, truth_mags, obs_mags, obs_errs):
        print(f"  {name:<8s}   {tm:9.4f}    {om:9.4f}   {oe:6.3f}")


def print_dynamic_summary(result, free_names, true_star):
    """Summarise arbitrary free-parameter posteriors.

    This avoids assuming the posterior columns are always Teff/logg/[M/H].
    """
    print("\nPosterior summary")
    print("  parameter          median        -1sigma       +1sigma        truth")
    print("  ---------       ----------    ----------    ----------    ----------")

    true_values = {
        "teff": true_star["teff"],
        "logg": true_star["logg"],
        "meta": true_star["meta"],
        "a_v": true_star["a_v"],
        "d": true_star["distance_pc"],
    }

    for i, name in enumerate(free_names):
        samples = result.samples[:, i]
        if name == "d":
            samples = samples / PC_TO_CM

        med = np.percentile(samples, 50.0)
        lo = np.percentile(samples, 15.865)
        hi = np.percentile(samples, 84.135)
        label = DISPLAY_NAMES.get(name, name)
        truth = true_values.get(name, np.nan)

        print(
            f"  {label:<13s}   {med:10.4f}    {med-lo:10.4f}    {hi-med:10.4f}    {truth:10.4f}"
        )


SHOW_PLOTS = False


def plot_sed(truth, path="sed_playground_sed.png"):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(truth.wavelengths, truth.observed_flux)
    ax.set_xlabel("Wavelength [Å]")
    ax.set_ylabel("Observed flux")
    ax.set_title("Forward-model observed SED")
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    print(f"\nSaved SED plot: {path}")
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig)



def _plot_samples_and_truths(posterior, fit_params, true_star):
    samples = np.asarray(posterior.samples, dtype=float)
    free_names = list(fit_params.free_names)

    label_map = {
        "teff": r"$T_{\rm eff}$ (K)",
        "logg": r"$\log g$",
        "meta": r"$[\rm M/H]$",
        "a_v":  r"$A_V$ (mag)",
        "d":    r"$d$ (pc)",
    }

    truth_map = {
        "teff": true_star["teff"],
        "logg": true_star["logg"],
        "meta": true_star["meta"],
        "a_v":  true_star["a_v"],
        "d":    true_star["distance_pc"],
    }

    plot_samples = samples.copy()
    for i, name in enumerate(free_names):
        if name == "d":
            plot_samples[:, i] = plot_samples[:, i] / PC_TO_CM

    labels = [label_map.get(name, name) for name in free_names]
    truths = [truth_map.get(name, None) for name in free_names]

    return plot_samples, free_names, labels, truths


def plot_posterior_diagnostics(
    posterior,
    fit_params,
    true_star,
    corner_path="inverse_corner.png",
    trace_path="inverse_chains.png",
):
    plot_samples, free_names, labels, truths = _plot_samples_and_truths(
        posterior, fit_params, true_star
    )

    npar = len(free_names)
    if npar == 0:
        print("No free parameters; skipping posterior plots.")
        return

    fig, axes = plt.subplots(npar, npar, figsize=(3.0 * npar, 3.0 * npar))

    if npar == 1:
        axes = np.array([[axes]])

    fig.suptitle("Posterior distributions", fontsize=13)

    for i in range(npar):
        for j in range(npar):
            ax = axes[i, j]

            if i == j:
                x = plot_samples[:, i]

                ax.hist(
                    x,
                    bins=40,
                    color="steelblue",
                    alpha=0.7,
                    density=True,
                    edgecolor="none",
                )

                med = float(np.percentile(x, 50.0))
                lo = float(np.percentile(x, 15.865))
                hi = float(np.percentile(x, 84.135))

                ax.axvline(med, color="navy", lw=1.5, label=f"median = {med:.4g}")
                ax.axvline(lo, color="navy", lw=0.8, ls="--")
                ax.axvline(hi, color="navy", lw=0.8, ls="--")

                if truths[i] is not None:
                    fmt = ".1f" if free_names[i] == "teff" else ".4g"
                    ax.axvline(
                        truths[i],
                        color="crimson",
                        lw=1.5,
                        label=f"true = {truths[i]:{fmt}}",
                    )

                ax.set_xlabel(labels[i], fontsize=10)
                ax.set_yticks([])
                ax.legend(fontsize=7, loc="best")

            elif i > j:
                ax.scatter(
                    plot_samples[:, j],
                    plot_samples[:, i],
                    s=2,
                    alpha=0.12,
                    color="steelblue",
                    rasterized=True,
                )

                if truths[j] is not None and truths[i] is not None:
                    ax.plot(truths[j], truths[i], "r+", ms=10, mew=1.5)

                ax.set_xlabel(labels[j], fontsize=10)
                ax.set_ylabel(labels[i], fontsize=10)

            else:
                ax.set_visible(False)

    fig.tight_layout()
    fig.savefig(corner_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {corner_path}")

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig)

    fig2, axes2 = plt.subplots(npar, 1, figsize=(10, 2.0 * npar), sharex=True)

    if npar == 1:
        axes2 = [axes2]

    for i, ax in enumerate(axes2):
        ax.plot(plot_samples[:, i], lw=0.4, alpha=0.6, color="steelblue")

        if truths[i] is not None:
            ax.axhline(truths[i], color="crimson", lw=1.2, ls="--")

        ax.set_ylabel(labels[i], fontsize=10)
        ax.grid(True, alpha=0.3)

    axes2[-1].set_xlabel("Flattened posterior sample index", fontsize=10)
    fig2.suptitle("Flattened posterior traces", fontsize=11)

    fig2.tight_layout()
    fig2.savefig(trace_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {trace_path}")

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig2)

# =============================================================================
# Main
# =============================================================================


def main():
    print("Loading grid...")
    grid = load_grid(GRID_PATH)
    print(grid)

    print("\nLoading filters...")
    vega_path = VEGA_SED if VEGA_SED.exists() else None
    filters = load_filters(FILTER_PATHS, vega_sed_path=vega_path)
    filter_names = [f.name for f in filters]
    print("Filters:", ", ".join(filter_names))

    print("\nTrue forward model")
    for key, value in TRUE_STAR.items():
        print(f"  {key:<13s} = {value}")

    true_ext = make_extinction_model(
        enabled=TRUE_STAR["a_v"] > 0.0,
        law=TRUE_STAR["law"],
        r_v=TRUE_STAR["r_v"],
        a_v=TRUE_STAR["a_v"],
        gordon23_env=TRUE_STAR["gordon23_env"],
    )

    truth = run_forward(
        teff=TRUE_STAR["teff"],
        logg=TRUE_STAR["logg"],
        meta=TRUE_STAR["meta"],
        R=TRUE_STAR["radius_rsun"] * RSUN_TO_CM,
        d=TRUE_STAR["distance_pc"] * PC_TO_CM,
        grid=grid,
        filters=filters,
        mag_system=MAG_SYSTEM,
        interp_method=INTERP_METHOD,
        extinction=true_ext,
    )

    rng = np.random.default_rng(RANDOM_SEED)
    truth_mags = np.array([truth.magnitudes[name] for name in filter_names])
    obs_mags = truth_mags + rng.normal(0.0, OBS_SIGMA, size=len(filter_names))
    obs_errs = np.full(len(filter_names), OBS_SIGMA)
    print_observations(filter_names, truth_mags, obs_mags, obs_errs)

    fit_params = build_fit_params(grid, FIT_CONFIG)
    p0_centre = build_p0_centre(FIT_CONFIG)
    fit_ext = make_fit_extinction_model(fit_params)

    print("\nFit configuration")
    print(fit_params.summary())
    print("Starting point:", p0_centre)

    posterior = run_inverse(
        obs_magnitudes=obs_mags,
        obs_uncertainties=obs_errs,
        filter_names=filter_names,
        R=TRUE_STAR["radius_rsun"] * RSUN_TO_CM,
        grid=grid,
        filters=filters,
        fit_params=fit_params,
        extinction=fit_ext,
        mag_system=MAG_SYSTEM,
        interp_method=INTERP_METHOD,
        p0_centre=p0_centre,
        **MCMC,
    )

    print_dynamic_summary(posterior, fit_params.free_names, TRUE_STAR)
    plot_sed(truth)

    plot_posterior_diagnostics(
        posterior=posterior,
        fit_params=fit_params,
        true_star=TRUE_STAR,
        corner_path="inverse_corner.png",
        trace_path="inverse_chains.png",
    )


if __name__ == "__main__":
    main()
