"""Empty-node detection on the real Kurucz2003all grid.

Of the 2,723 cells whose eight corners all hold a model, 455 still touch an
empty node in the Hermite derivative stencil.  ``stencil=False`` flags none
of them, and ``run_forward`` warns for the Hermite methods there while
returning finite, unchanged values.
"""

from __future__ import annotations

import os
import warnings

import numpy as np
import pytest

from _data_paths import KURUCZ_DIR, JOHNSON_DIR, VEGA_SED

pytest.importorskip("sed_model.cc_api", reason="cc_api Fortran extension not built")

if not os.path.isdir(KURUCZ_DIR):
    pytest.skip(f"Kurucz2003all grid not found at {KURUCZ_DIR}", allow_module_level=True)

N_CORNER_VALID = 2723
N_STENCIL_FLAGGED = 455


@pytest.fixture(scope="module")
def grid():
    from sed_model import load_grid
    return load_grid(KURUCZ_DIR)


@pytest.fixture(scope="module")
def valid_cells(grid):
    """Cell-centre queries for every cell with all eight corners populated."""
    has = grid._node_has_model
    t, g, m = grid.teff_grid, grid.logg_grid, grid.meta_grid
    corners = (has[:-1, :-1, :-1] & has[1:, :-1, :-1] & has[:-1, 1:, :-1]
               & has[1:, 1:, :-1] & has[:-1, :-1, 1:] & has[1:, :-1, 1:]
               & has[:-1, 1:, 1:] & has[1:, 1:, 1:])
    out = []
    for i, j, k in np.argwhere(corners):
        out.append((0.5 * (t[i] + t[i + 1]), 0.5 * (g[j] + g[j + 1]),
                    0.5 * (m[k] + m[k + 1])))
    return out


def test_corner_valid_cell_count(valid_cells):
    assert len(valid_cells) == N_CORNER_VALID


def test_stencil_flags_exactly_455(grid, valid_cells):
    flagged = [q for q in valid_cells if grid.cell_has_empty_nodes(*q, stencil=True)]
    assert len(flagged) == N_STENCIL_FLAGGED


def test_no_stencil_flags_none(grid, valid_cells):
    assert not any(grid.cell_has_empty_nodes(*q, stencil=False) for q in valid_cells)


@pytest.mark.skipif(not os.path.isdir(JOHNSON_DIR) or not os.path.exists(VEGA_SED),
                    reason="Johnson filters or Vega SED not found")
@pytest.mark.parametrize("method", ["hermite", "hermite_bounded"])
def test_run_forward_warns_and_stays_finite(grid, valid_cells, method):
    from sed_model import load_filters_from_instrument_dir, run_forward
    filters = load_filters_from_instrument_dir(JOHNSON_DIR, vega_sed_path=VEGA_SED)
    q = next(q for q in valid_cells if grid.cell_has_empty_nodes(*q, stencil=True))

    def fwd():
        return run_forward(teff=q[0], logg=q[1], meta=q[2], R=6.957e10, d=3.0857e19,
                           grid=grid, filters=filters, mag_system="Vega",
                           interp_method=method)

    with pytest.warns(UserWarning, match="empty nodes"):
        res = fwd()
    assert all(np.isfinite(v) for v in res.magnitudes.values())
    assert np.all(np.isfinite(res.observed_flux))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        again = fwd()
    assert res.magnitudes == again.magnitudes
