"""
demos/plot_utils.py
===================
Shared posterior visualisation utilities for SED_Model demos.

Import from a demo script with::

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent))
    from plot_utils import plot_posterior, plot_sed_posterior
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from sed_model.io import InverseResult
    from sed_model.params import FitParams
    from sed_model.grid import AtmosphereGrid
    from sed_model.filters import Filter
    from sed_model.sed_extinction import ExtinctionModel

from sed_model import PC_TO_CM

_LABELS: dict[str, str] = {
    "teff": r"$T_{\rm eff}$ (K)",
    "logg": r"$\log g$",
    "meta": r"$[\rm M/H]$",
    "a_v":  r"$A_V$ (mag)",
    "d":    r"$d$ (pc)",
}


# ---------------------------------------------------------------------------
# Corner / trace plot
# ---------------------------------------------------------------------------

def plot_posterior(
    result: "InverseResult",
    corner_path: str | Path,
    chains_path: str | Path,
    truths: dict[str, float] | None = None,
) -> None:
    """Save a corner plot and chain traces for an InverseResult.

    Parameters
    ----------
    result : InverseResult
    corner_path : path for the corner-plot PNG
    chains_path : path for the chain-trace PNG
    truths : optional {param_name: true_value} dict drawn as reference lines
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for posterior plots") from exc

    param_names = list(result.param_names)
    n = len(param_names)
    if n == 0:
        print("No free parameters; skipping posterior plots.")
        return

    # Convert distance column to parsecs for display
    samples = np.asarray(result.samples, dtype=float).copy()
    for i, name in enumerate(param_names):
        if name == "d":
            samples[:, i] /= PC_TO_CM

    labels     = [_LABELS.get(name, name) for name in param_names]
    truth_vals = [truths.get(name) if truths else None for name in param_names]

    # ------------------------------------------------------------------
    # Corner plot
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(n, n, figsize=(3.0 * n, 3.0 * n), squeeze=False)
    fig.suptitle("Posterior distributions", fontsize=13)

    for row in range(n):
        for col in range(n):
            ax = axes[row, col]
            if row == col:
                vals = samples[:, row]
                ax.hist(vals, bins=40, color="steelblue", alpha=0.7,
                        density=True, edgecolor="none")
                med, lo, hi = np.percentile(vals, [50.0, 15.865, 84.135])
                ax.axvline(med, color="navy", lw=1.5, label=f"median={med:.4g}")
                ax.axvline(lo,  color="navy", lw=0.8, ls="--")
                ax.axvline(hi,  color="navy", lw=0.8, ls="--")
                if truth_vals[row] is not None:
                    ax.axvline(truth_vals[row], color="crimson", lw=1.5,
                               label=f"true={truth_vals[row]:.4g}")
                ax.set_xlabel(labels[row], fontsize=9)
                ax.set_yticks([])
                ax.legend(fontsize=7, loc="best")
            elif row > col:
                ax.scatter(samples[:, col], samples[:, row],
                           s=2, alpha=0.12, color="steelblue", rasterized=True)
                if truth_vals[col] is not None and truth_vals[row] is not None:
                    ax.plot(truth_vals[col], truth_vals[row], "r+", ms=10, mew=1.5)
                ax.set_xlabel(labels[col], fontsize=9)
                ax.set_ylabel(labels[row], fontsize=9)
            else:
                ax.set_visible(False)

    corner_path = Path(corner_path)
    corner_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(corner_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {corner_path}")

    # ------------------------------------------------------------------
    # Chain traces
    # ------------------------------------------------------------------
    fig2, axes2 = plt.subplots(n, 1, figsize=(10, 2.0 * n), sharex=True, squeeze=False)
    for i, ax in enumerate(axes2[:, 0]):
        ax.plot(samples[:, i], lw=0.4, alpha=0.6, color="steelblue")
        if truth_vals[i] is not None:
            ax.axhline(truth_vals[i], color="crimson", lw=1.2, ls="--")
        ax.set_ylabel(labels[i], fontsize=9)
        ax.grid(True, alpha=0.3)
    axes2[-1, 0].set_xlabel("Flattened posterior sample index", fontsize=10)
    fig2.suptitle("Posterior traces (post burn-in)", fontsize=11)

    chains_path = Path(chains_path)
    chains_path.parent.mkdir(parents=True, exist_ok=True)
    fig2.tight_layout()
    fig2.savefig(chains_path, dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print(f"Saved: {chains_path}")


# ---------------------------------------------------------------------------
# SED posterior plot
# ---------------------------------------------------------------------------

def plot_sed_posterior(
    result: "InverseResult",
    fit_params: "FitParams",
    grid: "AtmosphereGrid",
    filters: list["Filter"],
    path: str | Path,
    n_samples: int = 100,
    interp_method: str = "hermite",
    extinction: "ExtinctionModel | None" = None,
    rng: "np.random.Generator | None" = None,
) -> None:
    """Plot the posterior SED alongside the observed photometry.

    Draws *n_samples* random posterior samples, reruns the forward model on
    each, and shows the resulting 16–84th percentile SED band with the
    median overlaid.  Observed photometric points (circles) and median
    predicted band fluxes (squares) are overplotted at filter pivot
    wavelengths.  A sub-panel below shows each filter's transmission curve
    in matching colours so the wavelength coverage is immediately visible.

    Parameters
    ----------
    result : InverseResult
    fit_params : FitParams
        The same FitParams used for the fit.
    grid : AtmosphereGrid
    filters : list of Filter
        Loaded filter objects in the same order as ``result.filter_names``.
    path : output PNG path
    n_samples : number of posterior samples to use for the SED ensemble
    interp_method : {'hermite', 'linear'}
    extinction : ExtinctionModel or None
        The same extinction model used for the fit.  Required when Av was
        free so each sample uses the correct per-sample A_V.
    rng : optional random generator for reproducible sample selection
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for SED posterior plots") from exc

    from sed_model import run_forward

    if rng is None:
        rng = np.random.default_rng()

    # ------------------------------------------------------------------
    # Sample the posterior and run forward models
    # ------------------------------------------------------------------
    n_avail = len(result.samples)
    chosen  = rng.choice(n_avail, size=min(n_samples, n_avail), replace=False)

    seds:       list[np.ndarray] = []
    band_fluxes: dict[str, list[float]] = {fname: [] for fname in result.filter_names}
    wavelengths: np.ndarray | None = None

    for idx in chosen:
        theta = result.samples[idx]
        try:
            fwd = run_forward(
                fit_params=fit_params,
                theta=theta,
                R=result.R,
                grid=grid,
                filters=filters,
                mag_system=result.mag_system,
                interp_method=interp_method,
                extinction=extinction,
            )
        except Exception:
            continue
        seds.append(fwd.observed_flux)
        if wavelengths is None:
            wavelengths = fwd.wavelengths
        for fname in result.filter_names:
            band_fluxes[fname].append(fwd.band_fluxes.get(fname, np.nan))

    if not seds:
        print("Warning: no valid forward models; skipping SED posterior plot.")
        return

    seds_arr   = np.array(seds)                # (n_valid, n_wave)
    wave_um    = wavelengths / 1e4             # Å → µm
    sed_lo, sed_med, sed_hi = np.percentile(seds_arr, [16, 50, 84], axis=0)

    # ------------------------------------------------------------------
    # Filter utilities
    # ------------------------------------------------------------------
    filter_map = {f.name: f for f in filters}

    def _pivot_um(filt: "Filter") -> float:
        num = np.trapezoid(filt.wavelengths * filt.transmission, filt.wavelengths)
        den = np.trapezoid(filt.transmission, filt.wavelengths)
        return (num / den / 1e4) if den > 0 else np.nan

    # Observed magnitudes → in-band fluxes (erg/s/cm²/Å)
    obs_flux:     dict[str, float] = {}
    obs_flux_err: dict[str, float] = {}
    pivot:        dict[str, float] = {}
    for fname, m, e in zip(result.filter_names, result.obs_magnitudes, result.obs_uncertainties):
        filt = filter_map.get(fname)
        if filt is None:
            continue
        zp         = filt.zero_point(result.mag_system)
        f          = zp * 10.0**(-0.4 * m)
        obs_flux[fname]     = f
        obs_flux_err[fname] = 0.4 * np.log(10.0) * f * e
        pivot[fname]        = _pivot_um(filt)

    # Median predicted band flux across the sampled posterior
    pred_flux_med: dict[str, float] = {}
    for fname, vals in band_fluxes.items():
        finite = [v for v in vals if np.isfinite(v)]
        if finite:
            pred_flux_med[fname] = float(np.median(finite))

    # ------------------------------------------------------------------
    # Assign one colour per filter, consistent across both panels
    # ------------------------------------------------------------------
    cmap   = plt.get_cmap("tab10")
    n_filt = len(result.filter_names)
    fcolors = {fname: cmap(i % 10) for i, fname in enumerate(result.filter_names)}

    # ------------------------------------------------------------------
    # Layout: SED panel (top) + transmission panel (bottom)
    # ------------------------------------------------------------------
    fig, (ax_sed, ax_t) = plt.subplots(
        2, 1, figsize=(11, 7),
        gridspec_kw={"height_ratios": [4, 1]},
        sharex=True,
    )

    # --- SED credible band ---
    ax_sed.fill_between(
        wave_um, sed_lo, sed_hi,
        color="steelblue", alpha=0.25, label="16–84th pctile",
    )
    ax_sed.plot(wave_um, sed_med, color="steelblue", lw=1.5, label="Median SED")

    # --- Photometric data points and model predictions ---
    for fname in result.filter_names:
        if fname not in obs_flux or fname not in pivot:
            continue
        col   = fcolors[fname]
        piv   = pivot[fname]
        short = fname.split("/")[-1]   # strip facility/instrument prefix

        # Observed (circle)
        ax_sed.errorbar(
            piv, obs_flux[fname], yerr=obs_flux_err[fname],
            fmt="o", color=col, ms=7, capsize=3, zorder=5, label=short,
        )
        # Predicted median (square, same colour, slightly smaller)
        if fname in pred_flux_med:
            ax_sed.plot(
                piv, pred_flux_med[fname],
                "s", color=col, ms=5, alpha=0.7, zorder=4,
            )

    ax_sed.set_yscale("log")
    ax_sed.set_ylabel(r"$F_\lambda$  (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)", fontsize=10)
    ax_sed.set_title("Posterior SED  (○ observed  □ predicted)", fontsize=11)
    ax_sed.legend(fontsize=8, loc="upper right", ncol=max(1, n_filt // 6 + 1))
    ax_sed.grid(True, alpha=0.2)

    # --- Filter transmission sub-panel ---
    for fname in result.filter_names:
        filt = filter_map.get(fname)
        if filt is None:
            continue
        col      = fcolors[fname]
        fw_um    = filt.wavelengths / 1e4
        t_max    = filt.transmission.max()
        t_norm   = filt.transmission / t_max if t_max > 0 else filt.transmission
        ax_t.fill_between(fw_um, 0, t_norm, color=col, alpha=0.35)
        ax_t.plot(fw_um, t_norm, color=col, lw=0.8)

    ax_t.set_ylim(0, 1.15)
    ax_t.set_ylabel("T (norm.)", fontsize=9)
    ax_t.set_xlabel("Wavelength (µm)", fontsize=10)
    ax_t.grid(True, alpha=0.2)

    # --- x-axis limits driven by the observed filters ---
    all_pivots = [v for v in pivot.values() if np.isfinite(v)]
    if all_pivots:
        x_lo = max(0.2, min(all_pivots) * 0.65)
        x_hi = min(6.0, max(all_pivots) * 1.55)
        ax_sed.set_xlim(x_lo, x_hi)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")
