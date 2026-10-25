"""
sed_model.grid
===================
Loads a SED_Tools-prepared stellar atmosphere grid (flux_cube.bin +
lookup_table.csv) into memory and exposes the axes and flux array
for downstream interpolation.

Binary format (flux_cube.bin)
------------------------------
  Header  : 4 × int32  -> (n_teff, n_logg, n_meta, n_lambda)
  Axes    : n_teff + n_logg + n_meta + n_lambda float64 values (written
             consecutively, no padding)
  Payload : n_teff × n_logg × n_meta × n_lambda float64 values stored
             in (W, M, L, T) order on disk, loaded into (T, L, M, W).

This matches the layout written by SED_Tools and read by the MESA
colors module (colors_utils.f90 :: load_flux_cube).
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Public dataclass
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AtmosphereGrid:
    """Immutable container for a loaded stellar atmosphere grid.

    Attributes
    ----------
    teff_grid : ndarray, shape (n_T,)
        Effective temperature grid points in Kelvin.
    logg_grid : ndarray, shape (n_L,)
        log10(g / cm s^-2) grid points.
    meta_grid : ndarray, shape (n_M,)
        Metallicity [M/H] grid points.
    wavelengths : ndarray, shape (n_W,)
        Wavelength grid in Angstroms.
    flux : ndarray, shape (n_T, n_L, n_M, n_W)
        Surface flux in erg/s/cm^2/Å at each (Teff, logg, [M/H]) node.
    flux_fortran : ndarray, shape (n_T, n_L, n_M, n_W)
        Fortran-contiguous copy of ``flux`` for zero-copy Fortran calls.
        Pre-computed at load time so the hot path avoids repeated copies.
    model_dir : Path
        Directory from which the grid was loaded.
    lookup_nodes : ndarray, shape (N, 3), optional
        (Teff, logg, [M/H]) of the tabulated models, i.e. the cube nodes
        that hold a spectrum.  MESA's ``Interp_rad`` is measured to the
        entries of ``lookup_table.csv``, which lists exactly these models
        (a rectangular cube also has empty nodes where no model exists).
        Computed from ``flux`` when not supplied.
    """
    teff_grid:    np.ndarray
    logg_grid:    np.ndarray
    meta_grid:    np.ndarray
    wavelengths:  np.ndarray
    flux:         np.ndarray
    flux_fortran: np.ndarray
    model_dir:    Path
    lookup_nodes: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    @property
    def teff_bounds(self) -> tuple[float, float]:
        return float(self.teff_grid[0]), float(self.teff_grid[-1])

    @property
    def logg_bounds(self) -> tuple[float, float]:
        return float(self.logg_grid[0]), float(self.logg_grid[-1])

    @property
    def meta_bounds(self) -> tuple[float, float]:
        return float(self.meta_grid[0]), float(self.meta_grid[-1])

    def in_bounds(self, teff: float, logg: float, meta: float) -> bool:
        """Return True if (teff, logg, meta) is within the grid axes."""
        return (
            self.teff_bounds[0] <= teff <= self.teff_bounds[1]
            and self.logg_bounds[0] <= logg <= self.logg_bounds[1]
            and self.meta_bounds[0] <= meta <= self.meta_bounds[1]
        )

    def clamp(
        self, teff: float, logg: float, meta: float
    ) -> tuple[float, float, float]:
        """Clamp (teff, logg, meta) to the grid boundary."""
        teff = float(np.clip(teff, *self.teff_bounds))
        logg = float(np.clip(logg, *self.logg_bounds))
        meta = float(np.clip(meta, *self.meta_bounds))
        return teff, logg, meta

    def _tabulated_nodes(self) -> np.ndarray:
        """(N, 3) array of (Teff, logg, [M/H]) of the nodes holding a model."""
        if self.lookup_nodes is not None:
            return self.lookup_nodes
        return _tabulated_nodes(self.teff_grid, self.logg_grid, self.meta_grid, self.flux)

    @cached_property
    def _node_has_model(self) -> np.ndarray:
        """Boolean (n_T, n_L, n_M) cube, True where a node holds a spectrum."""
        return self.flux.max(axis=3) > 0.0

    def cell_has_empty_nodes(self, teff: float, logg: float, meta: float) -> bool:
        """True if the cell containing the query has a corner with no model.

        The cell is the one MESA's interpolators use (out-of-range queries
        take the boundary cell, dummy and single-point axes collapse).  A
        corner counts only if it carries interpolation weight, so a query
        sitting exactly on a populated node is never reported.  Corners that
        hold no spectrum are all-zero in the cube, and interpolating through
        them gives unreliable fluxes.
        """
        sel = []
        for axis, val in ((self.teff_grid, teff), (self.logg_grid, logg),
                          (self.meta_grid, meta)):
            n = len(axis)
            if n < 2 or np.all(axis == 0.0) or np.all(axis == 999.0) \
                    or np.all(axis == -999.0):
                sel.append([0])
                continue
            if val <= axis[0]:
                sel.append([0])
            elif val >= axis[-1]:
                sel.append([n - 1])
            else:
                i = int(np.searchsorted(axis, val, side="right")) - 1
                sel.append([i] if axis[i] == val else [i, i + 1])
        return not bool(self._node_has_model[np.ix_(*sel)].all())

    def nearest_grid_distance(
        self, teff: float, logg: float, meta: float
    ) -> float:
        """Distance in normalised parameter space from the nearest tabulated
        model; identical to the ``Interp_rad`` diagnostic of the MESA colors
        module (``compute_interp_radius`` in colors/private/bolometric.f90).

        Each axis is normalised by the range spanned by the tabulated models
        (floored at 1e-12), and an axis on which every model has the value 0,
        999 or -999 (i.e. the axis does not exist for this grid) is ignored.
        The query is not clamped to the grid.
        """
        nodes = self._tabulated_nodes()
        q = (teff, logg, meta)
        total = np.zeros(len(nodes))
        for k in range(3):
            col = nodes[:, k]
            if np.all(col == 0.0) or np.all(col == 999.0) or np.all(col == -999.0):
                continue
            lo = col.min()
            rng = max(col.max() - lo, 1.0e-12)
            total += ((q[k] - lo) / rng - (col - lo) / rng) ** 2
        return float(np.sqrt(total.min()))

    def __repr__(self) -> str:
        return (
            f"AtmosphereGrid("
            f"Teff={self.teff_bounds[0]:.0f}–{self.teff_bounds[1]:.0f} K, "
            f"logg={self.logg_bounds[0]:.2f}–{self.logg_bounds[1]:.2f}, "
            f"[M/H]={self.meta_bounds[0]:.2f}–{self.meta_bounds[1]:.2f}, "
            f"n_wave={len(self.wavelengths)}, "
            f"flux={self.flux.shape})"
        )


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------

def load_grid(model_dir: str | Path) -> AtmosphereGrid:
    """Load a SED_Tools atmosphere grid from *model_dir*.

    The directory must contain:
      - ``flux_cube.bin``   — binary flux cube
      - ``lookup_table.csv`` — parameter index (used for validation only)

    Parameters
    ----------
    model_dir:
        Path to the directory produced by ``sed-tools rebuild``.

    Returns
    -------
    AtmosphereGrid

    Raises
    ------
    FileNotFoundError
        If ``flux_cube.bin`` is missing.
    ValueError
        If the binary header is invalid or the file is truncated.
    """
    model_dir = Path(model_dir).resolve()
    cube_path = model_dir / "flux_cube.bin"

    if not cube_path.exists():
        raise FileNotFoundError(
            f"flux_cube.bin not found in {model_dir}. "
            "Run 'sed-tools rebuild' to generate it."
        )

    teff_grid, logg_grid, meta_grid, wavelengths, flux = _read_flux_cube(cube_path)

    return AtmosphereGrid(
        teff_grid=teff_grid,
        logg_grid=logg_grid,
        meta_grid=meta_grid,
        wavelengths=wavelengths,
        flux=flux,
        flux_fortran=np.asfortranarray(flux, dtype=np.float64),
        model_dir=model_dir,
        lookup_nodes=_tabulated_nodes(teff_grid, logg_grid, meta_grid, flux),
    )


def _tabulated_nodes(teff_grid, logg_grid, meta_grid, flux) -> np.ndarray:
    """Coordinates of the cube nodes that carry a spectrum (flux > 0 somewhere)."""
    idx = np.argwhere(flux.max(axis=3) > 0.0)
    if len(idx) == 0:
        idx = np.argwhere(np.ones(flux.shape[:3], dtype=bool))
    return np.column_stack((teff_grid[idx[:, 0]], logg_grid[idx[:, 1]], meta_grid[idx[:, 2]]))


# ---------------------------------------------------------------------------
# Binary reader
# ---------------------------------------------------------------------------

def _read_flux_cube(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Parse a ``flux_cube.bin`` file.

    Returns (teff_grid, logg_grid, meta_grid, wavelengths, flux)
    where flux has shape (n_T, n_L, n_M, n_W).
    """
    with open(path, "rb") as fh:
        header_bytes = fh.read(16)
        if len(header_bytes) < 16:
            raise ValueError(f"{path}: file too short to contain a valid header")

        n_teff, n_logg, n_meta, n_lambda = struct.unpack("4i", header_bytes)

        if any(n <= 0 for n in (n_teff, n_logg, n_meta, n_lambda)):
            raise ValueError(
                f"{path}: invalid header dimensions "
                f"({n_teff}, {n_logg}, {n_meta}, {n_lambda})"
            )

        n_axis_total = n_teff + n_logg + n_meta + n_lambda
        axes_raw = np.frombuffer(fh.read(n_axis_total * 8), dtype=np.float64)
        if axes_raw.size != n_axis_total:
            raise ValueError(f"{path}: truncated axis data")

        offset = 0
        teff_grid   = axes_raw[offset : offset + n_teff].copy();  offset += n_teff
        logg_grid   = axes_raw[offset : offset + n_logg].copy();  offset += n_logg
        meta_grid   = axes_raw[offset : offset + n_meta].copy();  offset += n_meta
        wavelengths = axes_raw[offset : offset + n_lambda].copy()

        n_total = n_teff * n_logg * n_meta * n_lambda
        payload = np.frombuffer(fh.read(n_total * 8), dtype=np.float64)
        if payload.size != n_total:
            raise ValueError(
                f"{path}: flux payload truncated "
                f"(got {payload.size}, expected {n_total})"
            )

        # Disk layout: (n_lambda, n_meta, n_logg, n_teff) → transpose to (T, L, M, W)
        flux = (
            payload
            .reshape(n_lambda, n_meta, n_logg, n_teff)
            .transpose(3, 2, 1, 0)
            .copy()
        )

    return teff_grid, logg_grid, meta_grid, wavelengths, flux


# ---------------------------------------------------------------------------
# Optional validation helper
# ---------------------------------------------------------------------------

def validate_lookup_table(grid: AtmosphereGrid):
    """Read and return the ``lookup_table.csv`` from *grid.model_dir*.

    Useful for sanity-checking that the grid axes in the binary cube
    are consistent with the CSV metadata.  Not called at load time to
    avoid the pandas overhead on every initialisation.

    Returns
    -------
    pandas.DataFrame
    """
    import pandas as pd

    lookup_path = grid.model_dir / "lookup_table.csv"
    if not lookup_path.exists():
        raise FileNotFoundError(f"lookup_table.csv not found in {grid.model_dir}")

    df = pd.read_csv(lookup_path, dtype=str)
    df.columns = [c.strip().lstrip("#").strip() for c in df.columns]
    return df
