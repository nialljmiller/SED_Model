"""
sed_model.filters
======================
Loads photometric filter transmission curves from SED_Tools-prepared
``.dat`` files and precomputes Vega, AB, and ST zero-points using the
same photon-counting integrals as the MESA colors module
(colors/private/synthetic.f90).

Zero-point definitions
-----------------------
  Vega  :  F_zp = ∫ F_vega(λ) T(λ) λ dλ / ∫ T(λ) λ dλ
  AB    :  F_zp = ∫ F_AB(λ)   T(λ) λ dλ / ∫ T(λ) λ dλ
               where F_AB(λ) = 3.631e-20 × c / λ²  [erg/s/cm²/Å]
  ST    :  F_zp = 3.63e-9  (flat F_lambda, constant)

Zero-points are computed through the compiled cc_api Fortran extension so
the Python loader uses the same numerical kernels as run_forward() and the
current MESA Colors implementation.  A pure-Python trapezoid fallback is
kept for environments where the extension has not yet been built.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ._cc_ext import get_cc_api


# ---------------------------------------------------------------------------
# Constants  (cgs, wavelength in Angstroms)
# ---------------------------------------------------------------------------

_CLIGHT_CM_S = 2.99792458e10   # speed of light — only used in Python AB fallback
_AB_FNU_ZP   = 3.631e-20       # 3631 Jy in erg/s/cm^2/Hz
_ST_FLAM_ZP  = 3.63e-9         # flat f_lambda zero-point erg/s/cm^2/Å


# ---------------------------------------------------------------------------
# Public dataclass
# ---------------------------------------------------------------------------

@dataclass
class Filter:
    """A single photometric filter with precomputed zero-points.

    Attributes
    ----------
    name : str
        Short filter identifier (filename stem, e.g. ``'B'``, ``'Gbp'``).
    path : Path
        Absolute path to the source ``.dat`` file.
    wavelengths : ndarray, shape (n,)
        Filter wavelength grid in Angstroms.
    transmission : ndarray, shape (n,)
        Dimensionless transmission in [0, 1].
    vega_zero_point : float
        Photon-counting flux zero-point for the Vega magnitude system.
        ``-1.0`` if no Vega SED was supplied at load time.
    ab_zero_point : float
        Photon-counting flux zero-point for the AB magnitude system.
    st_zero_point : float
        Photon-counting flux zero-point for the ST magnitude system.
    """
    name:            str
    path:            Path
    wavelengths:     np.ndarray
    transmission:    np.ndarray
    vega_zero_point: float = field(default=-1.0)
    ab_zero_point:   float = field(default=-1.0)
    st_zero_point:   float = field(default=-1.0)

    def zero_point(self, system: str) -> float:
        """Return the zero-point for *system*.

        Parameters
        ----------
        system : str
            One of ``'Vega'``, ``'AB'``, ``'ST'`` (case-insensitive).
        """
        s = system.upper()
        if s == "VEGA":
            if self.vega_zero_point < 0:
                raise ValueError(
                    f"Filter '{self.name}': Vega zero-point not available. "
                    "Supply a Vega SED path to load_filters()."
                )
            return self.vega_zero_point
        if s == "AB":
            return self.ab_zero_point
        if s == "ST":
            return self.st_zero_point
        raise ValueError(f"Unknown magnitude system '{system}'. Choose Vega, AB, or ST.")

    def __repr__(self) -> str:
        return (
            f"Filter(name='{self.name}', "
            f"λ={self.wavelengths[0]:.0f}–{self.wavelengths[-1]:.0f} Å, "
            f"n={len(self.wavelengths)})"
        )


# ---------------------------------------------------------------------------
# Public loaders
# ---------------------------------------------------------------------------

def load_filters(
    filter_paths: list[str | Path],
    vega_sed_path: str | Path | None = None,
) -> list[Filter]:
    """Load a list of filter ``.dat`` files and return a list of
    :class:`Filter` objects with precomputed zero-points.

    Parameters
    ----------
    filter_paths:
        Paths to individual ``*.dat`` filter transmission files in the
        SED_Tools two-column format (wavelength Å, transmission 0–1).
        Comments starting with ``#`` are skipped.
    vega_sed_path:
        Optional path to a Vega reference SED CSV file
        (``wavelength,flux`` header, wavelength in Å,
        flux in erg/s/cm²/Å).  Required for Vega zero-points.

    Returns
    -------
    list of Filter, in the same order as *filter_paths*.
    """
    vega_wave: np.ndarray | None = None
    vega_flux: np.ndarray | None = None
    if vega_sed_path is not None:
        vega_wave, vega_flux = _load_vega_sed(Path(vega_sed_path))

    filters: list[Filter] = []
    for p in filter_paths:
        p = Path(p).resolve()
        wave, trans = _load_filter_dat(p)
        filt = Filter(
            name=p.stem,
            path=p,
            wavelengths=wave,
            transmission=trans,
        )
        filt.ab_zero_point = _compute_ab_zero_point(wave, trans)
        filt.st_zero_point = _compute_st_zero_point(wave, trans)
        if vega_wave is not None and vega_flux is not None:
            filt.vega_zero_point = _compute_vega_zero_point(
                vega_wave, vega_flux, wave, trans
            )
        filters.append(filt)

    return filters


def load_filters_from_instrument_dir(
    instrument_dir: str | Path,
    vega_sed_path: str | Path | None = None,
) -> list[Filter]:
    """Load all filters listed in a SED_Tools instrument directory.

    The directory must contain an index file whose name matches the
    last component of *instrument_dir* (e.g. ``Johnson/Johnson``) that
    lists one ``.dat`` filename per line.  This mirrors the structure
    expected by the MESA colors module.

    Filters are returned in index-file order, or alphabetically when
    the index file is absent.  Capture ``[f.name for f in filters]``
    immediately after loading to record the stable order for use with
    ``obs_magnitudes`` arrays.

    Parameters
    ----------
    instrument_dir:
        Path to a ``data/filters/<Facility>/<Instrument>/`` directory.
    vega_sed_path:
        Optional Vega reference SED for Vega zero-point computation.

    Returns
    -------
    list of Filter
    """
    instrument_dir = Path(instrument_dir).resolve()
    index_name = instrument_dir.name
    index_file = instrument_dir / index_name

    if not index_file.exists():
        dat_files = sorted(instrument_dir.glob("*.dat"))
        if not dat_files:
            raise FileNotFoundError(
                f"No index file '{index_name}' and no .dat files in {instrument_dir}"
            )
    else:
        dat_files = []
        for line in index_file.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            dat_files.append(instrument_dir / line)

    return load_filters(dat_files, vega_sed_path=vega_sed_path)


# ---------------------------------------------------------------------------
# Internal I/O helpers
# ---------------------------------------------------------------------------

def _sniff_delimiter(path: Path) -> str:
    """Return ',' if the file looks like CSV, else None (whitespace)."""
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            return "," if "," in line else None
    return None


def _load_filter_dat(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a two-column filter .dat file (wavelength Å, transmission).

    Handles both plain whitespace-separated files and CSV files with an
    optional header row (e.g. 'Wavelength,Transmission').
    """
    if not path.exists():
        raise FileNotFoundError(f"Filter file not found: {path}")

    delimiter = _sniff_delimiter(path)
    data = np.genfromtxt(path, comments="#", delimiter=delimiter,
                         invalid_raise=False)

    if data.ndim == 2:
        mask = np.isfinite(data[:, 0]) & np.isfinite(data[:, 1])
        data = data[mask]

    if data.ndim != 2 or data.shape[1] < 2 or len(data) == 0:
        raise ValueError(
            f"{path}: expected two-column file (wavelength, transmission)"
        )

    wave  = data[:, 0].astype(np.float64)
    trans = data[:, 1].astype(np.float64)

    if not np.all(np.diff(wave) > 0):
        order = np.argsort(wave)
        wave, trans = wave[order], trans[order]

    return wave, trans


def _load_vega_sed(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a Vega reference SED CSV (wavelength Å, flux erg/s/cm²/Å)."""
    if not path.exists():
        raise FileNotFoundError(f"Vega SED not found: {path}")

    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    wave = data[:, 0].astype(np.float64)
    flux = data[:, 1].astype(np.float64)

    if not np.all(np.diff(wave) > 0):
        order = np.argsort(wave)
        wave, flux = wave[order], flux[order]

    return wave, flux


# ---------------------------------------------------------------------------
# Zero-point computation  (mirrors synthetic.f90)
# ---------------------------------------------------------------------------

def _run_fortran_zero_point(fortran_method: str, *arrays: np.ndarray) -> float | None:
    """Attempt a Fortran zero-point kernel call.

    Returns the zero-point as a float, or None if the extension is absent or
    the kernel signals failure (ierr != 0).  Callers fall back to Python.
    """
    cc = get_cc_api()
    if cc is None:
        return None
    prepared = [np.ascontiguousarray(a, dtype=np.float64) for a in arrays]
    zp, ierr = getattr(cc, fortran_method)(*prepared)
    return float(zp) if int(ierr) == 0 else None


def _trapezoid(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.trapezoid(y, x))


def _compute_vega_zero_point(
    vega_wave: np.ndarray,
    vega_flux: np.ndarray,
    filt_wave: np.ndarray,
    filt_trans: np.ndarray,
) -> float:
    result = _run_fortran_zero_point(
        'vega_zero_point', vega_wave, vega_flux, filt_wave, filt_trans
    )
    if result is not None:
        return result

    trans_on_vega = np.interp(vega_wave, filt_wave, filt_trans, left=0.0, right=0.0)
    num = _trapezoid(vega_wave, vega_flux * trans_on_vega * vega_wave)
    den = _trapezoid(vega_wave, trans_on_vega * vega_wave)
    return num / den if den > 0.0 else -1.0


def _compute_ab_zero_point(
    filt_wave: np.ndarray,
    filt_trans: np.ndarray,
) -> float:
    result = _run_fortran_zero_point('ab_zero_point', filt_wave, filt_trans)
    if result is not None:
        return result

    f_ab = _AB_FNU_ZP * (_CLIGHT_CM_S * 1e8) / (filt_wave ** 2)
    num = _trapezoid(filt_wave, f_ab * filt_trans * filt_wave)
    den = _trapezoid(filt_wave, filt_trans * filt_wave)
    return num / den if den > 0.0 else -1.0


def _compute_st_zero_point(
    filt_wave: np.ndarray,
    filt_trans: np.ndarray,
) -> float:
    result = _run_fortran_zero_point('st_zero_point', filt_wave, filt_trans)
    if result is not None:
        return result

    f_st = np.full_like(filt_wave, _ST_FLAM_ZP)
    num = _trapezoid(filt_wave, f_st * filt_trans * filt_wave)
    den = _trapezoid(filt_wave, filt_trans * filt_wave)
    return num / den if den > 0.0 else -1.0
