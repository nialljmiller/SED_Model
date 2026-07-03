"""
sed_model.forward
======================
Forward model: stellar parameters → SED + synthetic photometry.

Pipeline
--------
  (Teff, logg, [M/H])  →  SED interpolation  (Fortran, Hermite or linear)
                       →  distance dilution   (Fortran, (R/d)²)
                       →  extinction          (Python, optional, sed_extinction)
                       →  bolometric          (Fortran)
                       →  filter convolution  (Fortran, per filter)
                       →  ForwardResult

All five physical parameters — Teff, logg, [M/H], Av, distance — can be
fixed or free, described by a FitParams object (sed_model.params).
This shared vocabulary is what makes the module bidirectional: the inverse
model unpacks a theta vector with FitParams.unpack and passes the result
directly to run_forward.

Calling conventions
-------------------
Classic (keyword-explicit)::

    run_forward(teff=5778, logg=4.44, meta=0.0,
                R=6.957e10, d=3.086e19,
                grid=grid, filters=filters)

FitParams (used by the inverse model)::

    run_forward(fit_params=params, theta=theta_vec,
                R=R_sun, grid=grid, filters=filters)

In the FitParams convention ``theta`` contains only the *free* parameters
in canonical order (teff, logg, meta, a_v, d — skipping fixed ones).
``fit_params.unpack(theta)`` fills in the fixed values.  If Av or
distance are free parameters they are taken from theta, not from any
keyword argument.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING

import numpy as np

from .grid import AtmosphereGrid
from .filters import Filter
from .params import FitParams
from ._cc_ext import get_cc_api

if TYPE_CHECKING:
    from sed_extinction import ExtinctionModel


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class ForwardResult:
    """Output of one forward-model evaluation.

    The five physical parameters that produced this result are all stored
    so the result is self-describing.  The inverse model uses
    ``ForwardResult.magnitudes`` to compute the likelihood; the forward
    model reads parameters from FitParams — both sides work with the same
    data structure.

    Attributes
    ----------
    nearest_grid_distance : float
        Euclidean distance in normalised parameter space from the nearest
        grid point (mirrors the MESA ``Interp_rad`` diagnostic).  Zero for
        an exact grid point, larger as the query moves away from the grid.
    magnitudes : dict[str, float]
        Synthetic magnitudes keyed by filter name, in the order the
        ``filters`` list was provided to ``run_forward``.
    """
    wavelengths:           np.ndarray
    surface_flux:          np.ndarray
    observed_flux:         np.ndarray
    magnitudes:            dict
    band_fluxes:           dict
    bol_flux:              float
    bol_mag:               float
    nearest_grid_distance: float
    clamped:               bool
    teff:                  float
    logg:                  float
    meta:                  float
    R:                     float
    d:                     float
    a_v:                   float = 0.0
    mag_system:            str   = "AB"
    extinction_applied:    bool  = False

    def __repr__(self) -> str:
        mags = ", ".join(f"{k}={v:.3f}" for k, v in self.magnitudes.items())
        ext  = f", Av={self.a_v:.3f}" if self.extinction_applied else ""
        return (
            f"ForwardResult("
            f"Teff={self.teff:.0f} K, logg={self.logg:.2f}, [M/H]={self.meta:.2f}"
            f"{ext}, bol_mag={self.bol_mag:.3f}, "
            f"mags=[{mags}], clamped={self.clamped})"
        )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_forward(
    # classic positional args
    teff:    Optional[float] = None,
    logg:    Optional[float] = None,
    meta:    Optional[float] = None,
    R:       Optional[float] = None,
    d:       Optional[float] = None,
    grid:    Optional[AtmosphereGrid] = None,
    filters: Optional[list] = None,
    # options
    mag_system:    str = "AB",
    interp_method: str = "hermite",
    extinction:    Optional["ExtinctionModel"] = None,
    # FitParams interface
    fit_params: Optional[FitParams] = None,
    theta:      Optional[np.ndarray] = None,
) -> ForwardResult:
    """Evaluate the forward model for one set of stellar parameters.

    Two calling conventions are supported.

    **Classic** (backward-compatible)::

        run_forward(teff, logg, meta, R, d, grid, filters)

    **FitParams** (used by the inverse model)::

        run_forward(fit_params=params, theta=theta_vec,
                    R=R_sun, grid=grid, filters=filters)

    In the FitParams convention ``theta`` contains only the *free* parameters
    in canonical order (teff, logg, meta, a_v, d — skipping fixed ones).
    ``fit_params.unpack(theta)`` fills in the fixed values.  If Av or
    distance are free parameters they are taken from theta, not from any
    keyword argument.

    Parameters
    ----------
    teff, logg, meta : float
        Atmospheric parameters.  Required in classic mode.
    R : float
        Stellar radius in cm.  Always required.
    d : float
        Distance in cm.  Required in classic mode; ignored if distance is
        a free parameter in fit_params.
    grid : AtmosphereGrid
    filters : list of Filter
    mag_system : {'AB', 'Vega', 'ST'}
        Case-insensitive.
    interp_method : {'hermite', 'linear'}
    extinction : ExtinctionModel or None
        Applied after dilution and before filter convolution.
        When Av is free in fit_params, the model's stored a_v is
        overridden by the value from theta at each call.
    fit_params : FitParams or None
    theta : array-like or None
        Required when fit_params is provided.
    """
    if R is None:
        raise ValueError("R (stellar radius in cm) is always required.")

    cc = get_cc_api(required=True)

    # ------------------------------------------------------------------
    # Resolve parameters from whichever calling convention is used
    # ------------------------------------------------------------------
    if fit_params is not None:
        if theta is None:
            raise ValueError("theta must be supplied when fit_params is used")
        p     = fit_params.unpack(np.asarray(theta, dtype=np.float64))
        teff  = float(p['teff'])
        logg  = float(p['logg'])
        meta  = float(p['meta'])
        d_use = float(p['d'])
        av    = float(p['a_v'])
    else:
        if any(x is None for x in (teff, logg, meta, d, grid, filters)):
            raise ValueError(
                "teff, logg, meta, d, grid, and filters must all be supplied "
                "when fit_params is not used."
            )
        d_use = float(d)
        av    = 0.0

    # ------------------------------------------------------------------
    # 1. SED interpolation (Fortran)
    # ------------------------------------------------------------------
    clamped  = not grid.in_bounds(teff, logg, meta)
    grid_dist = grid.nearest_grid_distance(teff, logg, meta)
    teff_q, logg_q, meta_q = grid.clamp(teff, logg, meta)

    # grid.flux_fortran is pre-computed at load time — zero-copy to Fortran.
    # The axis arrays are already float64 C-contiguous from load_grid.
    if interp_method == "hermite":
        surface_flux, ierr = cc.interp_sed_hermite(
            teff_q, logg_q, meta_q,
            grid.teff_grid, grid.logg_grid, grid.meta_grid,
            grid.flux_fortran,
        )
    elif interp_method == "linear":
        surface_flux, ierr = cc.interp_sed_linear(
            teff_q, logg_q, meta_q,
            grid.teff_grid, grid.logg_grid, grid.meta_grid,
            grid.flux_fortran,
        )
    else:
        raise ValueError(
            f"Unknown interp_method '{interp_method}'. "
            "Choose 'hermite' or 'linear'."
        )

    if ierr != 0:
        clamped = True

    # ------------------------------------------------------------------
    # 2. Distance dilution: F_obs = F_surface × (R/d)²  (Fortran)
    # ------------------------------------------------------------------
    observed_flux = cc.dilute_flux(surface_flux, float(R), d_use)

    # ------------------------------------------------------------------
    # 3. Extinction (Python, optional)
    #
    # Applied AFTER dilution, BEFORE filter convolution.
    # When Av is a free parameter, rebuild the ExtinctionModel with the
    # Av value from theta so every likelihood call uses the correct value.
    # ------------------------------------------------------------------
    extinction_applied = False
    if extinction is not None and getattr(extinction.config, 'enabled', False):
        if fit_params is not None:
            extinction = extinction.with_av(av)

        observed_flux      = extinction.apply(grid.wavelengths, observed_flux)
        extinction_applied = True
        av                 = extinction.config.a_v

    # ------------------------------------------------------------------
    # 4. Bolometric quantities (Fortran)
    # ------------------------------------------------------------------
    bol_flux, bol_mag, bol_ierr = cc.bolometric(grid.wavelengths, observed_flux)
    if bol_ierr != 0:
        clamped = True

    # ------------------------------------------------------------------
    # 5. Synthetic photometry per filter (Fortran)
    #
    # Filter arrays are float64 C-contiguous from load_filters — no copy needed.
    # ------------------------------------------------------------------
    magnitudes:  dict = {}
    band_fluxes: dict = {}

    for filt in filters:
        zp = filt.zero_point(mag_system)
        mag, band_flux, mag_ierr = cc.synthetic_magnitude(
            grid.wavelengths, observed_flux,
            filt.wavelengths, filt.transmission,
            zp,
        )
        magnitudes[filt.name]  = float(mag) if mag_ierr == 0 else float('nan')
        band_fluxes[filt.name] = float(band_flux)

    return ForwardResult(
        wavelengths=grid.wavelengths,
        surface_flux=surface_flux,
        observed_flux=observed_flux,
        magnitudes=magnitudes,
        band_fluxes=band_fluxes,
        bol_flux=float(bol_flux),
        bol_mag=float(bol_mag),
        nearest_grid_distance=float(grid_dist),
        clamped=clamped,
        teff=float(teff),
        logg=float(logg),
        meta=float(meta),
        R=float(R),
        d=d_use,
        a_v=av,
        mag_system=mag_system,
        extinction_applied=extinction_applied,
    )


# ---------------------------------------------------------------------------
# Batch wrapper (classic convention only)
# ---------------------------------------------------------------------------

def run_forward_batch(
    params: np.ndarray,
    R: float,
    d: float,
    grid: AtmosphereGrid,
    filters: list,
    mag_system: str = "AB",
    interp_method: str = "hermite",
    extinction: Optional["ExtinctionModel"] = None,
) -> list:
    """Run run_forward over an array of (teff, logg, meta) rows.

    Parameters
    ----------
    params : ndarray, shape (N, 3)
    Returns
    -------
    list of ForwardResult, length N.
    """
    params = np.atleast_2d(params)
    if params.shape[1] != 3:
        raise ValueError("params must have shape (N, 3): teff, logg, meta")
    return [
        run_forward(
            teff=float(r[0]), logg=float(r[1]), meta=float(r[2]),
            R=R, d=d, grid=grid, filters=filters,
            mag_system=mag_system, interp_method=interp_method,
            extinction=extinction,
        )
        for r in params
    ]
