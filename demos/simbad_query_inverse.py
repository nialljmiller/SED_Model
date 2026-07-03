#!/usr/bin/env python3
"""Resolve a SIMBAD object, fetch its filters, and run SED_Model inverse mode.

This is deliberately a standalone integration demo.  It does not add catalogue
code to SED_Model or SED_Tools; it only composes their current public APIs.

Example
-------
python demos/simbad_query_inverse.py "HD 209458" \
    --grid ~/SED_Tools/data/stellar_models/Kurucz2003all \
    --radius-rsun 1.20 --steps 4000 --burn 1000

The current SED_Model inverse API requires a fixed radius and one magnitude
system for the whole fit.  This script converts supported Vega magnitudes to
AB before fitting.  By default Teff/logg/[M/H] are free over the model grid,
Av is free in [0, 2] mag, and distance is fixed from SIMBAD's selected
parallax.  Each fit parameter can instead be fixed or given explicit limits.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from pathlib import Path
import sys
from typing import Any

import numpy as np

import sys
from sed_model import PC_TO_CM, RSUN_TO_CM
sys.path.insert(0, str(Path(__file__).parent))
from plot_utils import plot_posterior, plot_sed_posterior


@dataclass(frozen=True)
class BandSpec:
    facility: str
    instrument: str
    filename: str
    native_system: str
    note: str = ""


# SIMBAD's flux table uses short filter labels.  Keep this mapping deliberately
# small and auditable: an incorrect passband is worse than dropping a point.
SIMBAD_BANDS: dict[str, BandSpec] = {
    "U": BandSpec("Generic", "Johnson", "U.dat", "Vega"),
    "B": BandSpec("Generic", "Johnson", "B.dat", "Vega"),
    "V": BandSpec("Generic", "Johnson", "V.dat", "Vega"),
    "R": BandSpec("Generic", "Johnson", "R.dat", "Vega"),
    "I": BandSpec("Generic", "Johnson", "I.dat", "Vega"),
    "J": BandSpec("2MASS", "2MASS", "J.dat", "Vega"),
    "H": BandSpec("2MASS", "2MASS", "H.dat", "Vega"),
    # SIMBAD reports the 2MASS short-K measurement under the historical label K.
    "K": BandSpec("2MASS", "2MASS", "Ks.dat", "Vega", "SIMBAD K -> 2MASS Ks"),
    "G": BandSpec("GAIA", "GAIA", "G.dat", "Vega"),
}


@dataclass
class Observation:
    simbad_band: str
    magnitude: float
    uncertainty: float
    bibcode: str
    spec: BandSpec
    filter_name: str = ""
    ab_magnitude: float = np.nan


def _column(table: Any, *names: str) -> Any:
    """Return a column across astroquery/Astropy naming variations."""
    lookup = {str(name).lower(): name for name in table.colnames}
    for candidate in names:
        key = candidate.lower()
        if key in lookup:
            return table[lookup[key]]
    raise KeyError(f"None of {names!r} found; SIMBAD returned {table.colnames}")


def _value(value: Any, default: Any = None) -> Any:
    if value is None or np.ma.is_masked(value):
        return default
    if isinstance(value, bytes):
        return value.decode(errors="replace").strip()
    return value


def query_simbad(identifier: str, default_uncertainty: float | None) -> tuple[str, list[Observation], float]:
    try:
        from astroquery.simbad import Simbad
    except ImportError as exc:
        raise RuntimeError("astroquery is required (it is a SED_Tools dependency)") from exc

    flux_query = Simbad()
    flux_query.add_votable_fields("flux")
    fluxes = flux_query.query_object(identifier)
    if fluxes is None or len(fluxes) == 0:
        raise RuntimeError(f"SIMBAD did not resolve {identifier!r} or returned no fluxes")

    main_ids = _column(fluxes, "main_id")
    resolved = str(_value(main_ids[0], identifier)).strip()
    bands = _column(fluxes, "flux.filter", "flux_filter", "filter")
    mags = _column(fluxes, "flux", "flux.value", "flux_value")
    errors = _column(fluxes, "flux_err", "flux.err", "flux_error")
    try:
        bibcodes = _column(fluxes, "flux.bibcode", "flux_bibcode", "bibcode")
    except KeyError:
        bibcodes = [""] * len(fluxes)

    observations: list[Observation] = []
    for band_raw, mag_raw, err_raw, bib_raw in zip(bands, mags, errors, bibcodes):
        band = str(_value(band_raw, "")).strip()
        spec = SIMBAD_BANDS.get(band)
        mag = _value(mag_raw)
        err = _value(err_raw)
        if spec is None or mag is None:
            continue
        try:
            mag = float(mag)
            err = float(err) if err is not None else np.nan
        except (TypeError, ValueError):
            continue
        if not np.isfinite(mag):
            continue
        if not np.isfinite(err) or err <= 0:
            if default_uncertainty is None:
                continue
            err = default_uncertainty
        observations.append(
            Observation(band, mag, err, str(_value(bib_raw, "")).strip(), spec)
        )

    basic_query = Simbad()
    basic_query.add_votable_fields("parallax")
    basic = basic_query.query_object(identifier)
    if basic is None or len(basic) == 0:
        raise RuntimeError(f"Could not retrieve parallax for {resolved}")
    plx_col = _column(basic, "plx_value", "parallax", "plx")
    parallax_mas = float(_value(plx_col[0], np.nan))
    if not np.isfinite(parallax_mas) or parallax_mas <= 0:
        raise RuntimeError(
            f"{resolved} has no positive selected SIMBAD parallax; "
            "the current inverse API needs a fixed distance"
        )
    if not observations:
        raise RuntimeError(
            "No supported SIMBAD photometry with usable uncertainties was found. "
            "Use --default-uncertainty to admit measurements lacking errors."
        )
    return resolved, observations, parallax_mas


def acquire_filters(
    observations: list[Observation], filter_root: Path, vega_sed: Path
) -> tuple[list[Observation], list[Any]]:
    try:
        from sed_tools.api import Filters as ToolFilters
        from sed_model import load_filters
    except ImportError as exc:
        raise RuntimeError(
            "Both sed-tools and sed-model must be importable in this environment"
        ) from exc

    instrument_dirs: dict[tuple[str, str], Path] = {}
    for obs in observations:
        key = (obs.spec.facility, obs.spec.instrument)
        if key not in instrument_dirs:
            instrument_dirs[key] = Path(
                ToolFilters.fetch(*key, filter_root=filter_root)
            ).resolve()

    accepted: list[Observation] = []
    loaded_filters: list[Any] = []
    for obs in observations:
        directory = instrument_dirs[(obs.spec.facility, obs.spec.instrument)]
        path = directory / obs.spec.filename
        if not path.exists():
            print(f"[skip] {obs.simbad_band}: expected profile not found: {path}")
            continue
        filt = load_filters([path], vega_sed_path=vega_sed)[0]
        canonical = f"{obs.spec.facility}/{obs.spec.instrument}/{path.stem}"
        filt.name = canonical

        # m_AB = m_native - 2.5 log10(Fzp_native / Fzp_AB)
        if obs.spec.native_system.upper() == "VEGA":
            offset = -2.5 * np.log10(filt.vega_zero_point / filt.ab_zero_point)
            ab_mag = obs.magnitude + offset
        elif obs.spec.native_system.upper() == "AB":
            ab_mag = obs.magnitude
        else:
            print(f"[skip] {obs.simbad_band}: unknown system {obs.spec.native_system}")
            continue

        accepted.append(replace(obs, filter_name=canonical, ab_magnitude=float(ab_mag)))
        loaded_filters.append(filt)
    return accepted, loaded_filters




def parse_args() -> argparse.Namespace:
    default_tools = Path.home() / "SED_Tools" / "data"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("identifier", help="SIMBAD-resolvable identifier, e.g. 'HD 209458'")
    parser.add_argument("--grid", type=Path, required=True, help="SED_Tools atmosphere grid directory")
    parser.add_argument("--radius-rsun", type=float, required=True, help="Fixed stellar radius in solar radii")

    def add_parameter_options(
        name: str, label: str, unit: str = "", aliases: tuple[str, ...] = ()
    ) -> None:
        option_names = [f"--{name}", *(f"--{alias}" for alias in aliases)]
        limits_names = [
            f"--{name}-limits", *(f"--{alias}-limits" for alias in aliases)
        ]
        group = parser.add_mutually_exclusive_group()
        group.add_argument(
            *option_names,
            dest=name.replace("-", "_"),
            type=float,
            metavar="VALUE",
            help=f"Fix {label} at VALUE{f' {unit}' if unit else ''}",
        )
        group.add_argument(
            *limits_names,
            dest=f"{name.replace('-', '_')}_limits",
            type=float,
            nargs=2,
            metavar=("LOW", "HIGH"),
            help=f"Sample {label} between LOW and HIGH{f' {unit}' if unit else ''}",
        )

    add_parameter_options("teff", "Teff", "K")
    add_parameter_options("logg", "logg")
    add_parameter_options("meta", "[M/H]", "dex", aliases=("metallicity",))
    add_parameter_options("av", "Av", "mag")
    add_parameter_options("distance-pc", "distance", "pc", aliases=("distance",))
    parser.add_argument(
        "--filter-root", type=Path, default=default_tools / "filters",
        help="SED_Tools filter cache/download root",
    )
    parser.add_argument(
        "--vega-sed", type=Path, default=default_tools / "stellar_models" / "vega_flam.csv",
        help="Vega reference spectrum used for Vega-to-AB conversion",
    )
    parser.add_argument(
        "--default-uncertainty", type=float, default=None,
        help="Admit measurements without errors using this uncertainty in mag",
    )
    parser.add_argument("--walkers", type=int, default=32)
    parser.add_argument("--steps", type=int, default=2000)
    parser.add_argument("--burn", type=int, default=500)
    parser.add_argument("--thin", type=int, default=1)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--no-progress", action="store_true")
    parser.add_argument("--output", type=Path, help="Optional output .npz posterior")
    parser.add_argument(
        "--plot-prefix", type=Path, default=Path("inverse"),
        help="Plot filename prefix (default: inverse, producing inverse_corner.png and inverse_chains.png)",
    )
    parser.add_argument("--no-plots", action="store_true", help="Do not create posterior plots")
    parser.add_argument(
        "--extinction-law",
        default="fitzpatrick99",
        choices=["ccm89", "odonnell94", "fitzpatrick99", "fm07", "calzetti00", "gordon23"],
        help="Extinction law to use when Av is free (default: fitzpatrick99)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.radius_rsun <= 0:
        raise SystemExit("--radius-rsun must be positive")
    if args.default_uncertainty is not None and args.default_uncertainty <= 0:
        raise SystemExit("--default-uncertainty must be positive")
    if not args.vega_sed.expanduser().exists():
        raise SystemExit(f"Vega SED not found: {args.vega_sed.expanduser()}")

    resolved, observations, parallax_mas = query_simbad(
        args.identifier, args.default_uncertainty
    )
    distance_pc = 1000.0 / parallax_mas
    print(f"Resolved: {resolved}")
    print(f"SIMBAD selected parallax: {parallax_mas:.6g} mas ({distance_pc:.6g} pc)")

    accepted, filters = acquire_filters(
        observations, args.filter_root.expanduser(), args.vega_sed.expanduser()
    )
    if len(accepted) < 2:
        raise RuntimeError(f"Only {len(accepted)} usable band(s) remain; need at least 2")

    print("\nPhotometry admitted to the fit (normalized to AB):")
    for obs in accepted:
        note = f"; {obs.spec.note}" if obs.spec.note else ""
        print(
            f"  {obs.filter_name:28s} native={obs.magnitude:9.4f} "
            f"AB={obs.ab_magnitude:9.4f} err={obs.uncertainty:.4f} "
            f"bibcode={obs.bibcode or '-'}{note}"
        )

    from sed_model import load_grid, run_inverse

    grid = load_grid(args.grid.expanduser())

    def parameter_value(fixed_value: float | None, limits: list[float] | None, default: Any) -> Any:
        if fixed_value is not None:
            return fixed_value
        if limits is not None:
            return tuple(limits)
        return default

    from sed_model import fit_params_from_grid

    teff = parameter_value(args.teff, args.teff_limits, None)
    logg = parameter_value(args.logg, args.logg_limits, None)
    meta = parameter_value(args.meta, args.meta_limits, None)
    a_v = parameter_value(args.av, args.av_limits, (0.0, 2.0))
    distance_pc_spec = parameter_value(
        args.distance_pc, args.distance_pc_limits, distance_pc
    )
    if isinstance(distance_pc_spec, tuple):
        distance_spec = tuple(value * PC_TO_CM for value in distance_pc_spec)
    else:
        distance_spec = distance_pc_spec * PC_TO_CM

    fit_params = fit_params_from_grid(
        grid,
        teff=teff,
        logg=logg,
        meta=meta,
        a_v=a_v,
        d_cm=distance_spec,
    )
    if len(accepted) < fit_params.n_free:
        raise RuntimeError(
            f"Only {len(accepted)} usable bands but {fit_params.n_free} free parameters "
            f"({fit_params.free_names}); the fit is underdetermined"
        )

    # Build the extinction model explicitly so it can be reused for plotting.
    # Fitzpatrick (1999) is the default; pass --extinction-law to override.
    from sed_model.sed_extinction import make_extinction_model
    av_spec = fit_params.a_v
    if av_spec.is_free or (av_spec.is_fixed and av_spec.value > 0.0):
        ext_model = make_extinction_model(
            enabled=True,
            law=getattr(args, "extinction_law", "fitzpatrick99"),
            a_v=0.0 if av_spec.is_free else av_spec.value,
        )
    else:
        ext_model = None

    result = run_inverse(
        obs_magnitudes=[obs.ab_magnitude for obs in accepted],
        obs_uncertainties=[obs.uncertainty for obs in accepted],
        filter_names=[obs.filter_name for obs in accepted],
        R=args.radius_rsun * RSUN_TO_CM,
        grid=grid,
        filters=filters,
        fit_params=fit_params,
        extinction=ext_model,
        mag_system="AB",
        n_walkers=args.walkers,
        n_steps=args.steps,
        n_burn=args.burn,
        n_thin=args.thin,
        seed=args.seed,
        progress=not args.no_progress,
    )
    result.print_summary()
    if args.output:
        result.save(args.output.expanduser())
        print(f"Saved posterior: {args.output.expanduser().with_suffix('.npz')}")
    if not args.no_plots:
        plot_prefix = args.plot_prefix.expanduser()
        plot_posterior(
            result,
            plot_prefix.parent / f"{plot_prefix.name}_corner.png",
            plot_prefix.parent / f"{plot_prefix.name}_chains.png",
        )
        plot_sed_posterior(
            result,
            fit_params,
            grid,
            filters,
            plot_prefix.parent / f"{plot_prefix.name}_sed.png",
            extinction=ext_model,
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130) from None
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
