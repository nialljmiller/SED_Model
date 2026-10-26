"""
tests/test_inverse_guards.py
============================
Guard rails of the inverse fitter (0.4.4): Hermite-stencil rejection, the
non-convergence warning and the ``fit_jitter`` hint, on a small synthetic grid
with a hole (no data files needed), plus a seeded injection-recovery test on
the real Kurucz grid.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np
import pytest

from _data_paths import KURUCZ_DIR, JOHNSON_DIR, VEGA_SED

R_SUN = 6.957e10
D_10PC = 3.0857e19

try:
    from sed_model import cc_api  # noqa: F401
    HAVE_CC = True
except ImportError:
    HAVE_CC = False

pytestmark = pytest.mark.filterwarnings("ignore:The chain is not converged")

requires_fortran = pytest.mark.skipif(not HAVE_CC, reason="cc_api not built")
requires_data = pytest.mark.skipif(
    not (os.path.isdir(KURUCZ_DIR) and os.path.isdir(JOHNSON_DIR)),
    reason="Kurucz grid / Johnson filters not found",
)


def _hole_grid():
    """3x3x1 grid, all populated except (6000 K, 3.0).  Only the cell
    T<5000, logg<2 is clean for Hermite; linear only avoids cells whose
    corner is the hole."""
    from sed_model.grid import AtmosphereGrid, _tabulated_nodes

    T = np.array([4000.0, 5000.0, 6000.0])
    L = np.array([1.0, 2.0, 3.0])
    M = np.array([0.0])
    wave = np.linspace(1000.0, 9000.0, 8)
    flux = np.ones((3, 3, 1, 8))
    flux[2, 2, 0, :] = 0.0
    return AtmosphereGrid(T, L, M, wave, flux, np.asfortranarray(flux), Path("."),
                          lookup_nodes=_tabulated_nodes(T, L, M, flux))


def _tophat_filters():
    from sed_model.filters import Filter

    out = []
    for name, lo, hi in (("a", 2000.0, 3000.0), ("b", 4000.0, 5000.0),
                         ("c", 6000.0, 7000.0)):
        w = np.linspace(lo, hi, 20)
        out.append(Filter(name=name, path=Path("."), wavelengths=w,
                          transmission=np.ones_like(w), ab_zero_point=1.0,
                          vega_zero_point=1.0, st_zero_point=1.0))
    return out


def _hole_obs(offsets=(0.0, 0.0, 0.0)):
    """Magnitudes of the flat-flux grid at R_sun, 10 pc, shifted by offsets."""
    from sed_model import run_forward

    flt = _tophat_filters()
    r = run_forward(teff=4500.0, logg=1.5, meta=0.0, R=R_SUN, d=D_10PC,
                    grid=_hole_grid(), filters=flt, interp_method="linear")
    names = [f.name for f in flt]
    mags = np.array([r.magnitudes[n] for n in names]) + np.asarray(offsets)
    return mags, names, flt


def _fit(method="hermite_bounded", offsets=(0.0, 0.0, 0.0), **kw):
    from sed_model import run_inverse

    mags, names, flt = _hole_obs(offsets)
    from sed_model.params import FitParams, fixed, free

    g = _hole_grid()
    fp = FitParams(teff=free(4000.0, 6000.0, name="teff"),
                   logg=free(1.0, 3.0, name="logg"),
                   meta=fixed(0.0, name="meta"),
                   a_v=fixed(0.0, name="a_v"), d=fixed(D_10PC, name="d"))
    args = dict(R=R_SUN, fit_params=fp,
                grid=g, filters=flt,
                interp_method=method, n_walkers=16, n_steps=300, n_burn=100,
                progress=False, verbose=False, seed=3)
    args.update(kw)
    return run_inverse(mags, np.full(3, 0.05), names, **args)


# ---------------------------------------------------------------------------
# (a) stencil-aware rejection
# ---------------------------------------------------------------------------

@requires_fortran
def test_fastmodel_supported_depends_on_method():
    from sed_model._inv_core import _FastModel

    g, flt = _hole_grid(), _tophat_filters()
    names = [f.name for f in flt]
    herm = _FastModel(g, flt, names, R_SUN, "AB", "hermite_bounded", None)
    lin = _FastModel(g, flt, names, R_SUN, "AB", "linear", None)
    for q, h, li in (
        ((4500.0, 1.5, 0.0), True, True),     # clean cell
        ((4500.0, 2.5, 0.0), False, True),    # corners fine, stencil hits hole
        ((5500.0, 1.5, 0.0), False, True),    # stencil hit along Teff
        ((5500.0, 2.5, 0.0), False, False),   # empty corner
    ):
        assert herm.supported(*q) is h, q
        assert lin.supported(*q) is li, q
    # the cached mask agrees with the grid's own per-call check
    for q in ((4500.0, 1.5, 0.0), (4500.0, 2.5, 0.0), (5500.0, 1.5, 0.0)):
        assert herm.supported(*q) == (not g.cell_has_empty_nodes(*q, stencil=True))


@requires_fortran
def test_hermite_fit_stays_out_of_stencil_cells():
    res = _fit("hermite")
    t, l = res.samples[:, 0], res.samples[:, 1]
    assert np.all(t <= 5000.0 + 1e-9) and np.all(l <= 2.0 + 1e-9)
    assert np.all(np.isfinite(res.log_prob))


@requires_fortran
def test_linear_fit_may_use_stencil_cells():
    res = _fit("linear")
    assert np.any(res.samples[:, 0] > 5000.0) or np.any(res.samples[:, 1] > 2.0)


# ---------------------------------------------------------------------------
# (b) non-convergence warning, (c) fit_jitter hint
# ---------------------------------------------------------------------------

@requires_fortran
def test_warns_when_not_converged():
    with pytest.warns(UserWarning, match="not converged") as rec:
        res = _fit(n_steps=80, n_burn=20)
    assert res.converged is False
    assert any("auto_stop" in str(w.message) for w in rec)


@requires_fortran
def test_no_convergence_warning_when_converged():
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        res = _fit(n_steps=4000, n_burn=None)
    assert res.converged is True
    assert not any("not converged" in str(w.message) for w in rec)


@requires_fortran
def test_jitter_hint_on_poor_fit_only():
    bad = (0.0, 0.6, -0.6)     # inconsistent with a flat spectrum, err 0.05
    with pytest.warns(UserWarning, match="fit_jitter=True"):
        res = _fit(offsets=bad, n_steps=120, n_burn=30)
    assert res.chi2_best / res.n_data > 5.0

    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        _fit(offsets=bad, n_steps=120, n_burn=30, fit_jitter=True)
        _fit(n_steps=120, n_burn=30)                     # good fit
    assert not any("fit_jitter" in str(w.message) for w in rec)


# ---------------------------------------------------------------------------
# (d) seeded injection-recovery, loose bounds
# ---------------------------------------------------------------------------

@requires_fortran
@requires_data
def test_injection_recovery_seeded():
    from sed_model import (load_grid, load_filters_from_instrument_dir,
                           run_forward, run_inverse)

    grid = load_grid(KURUCZ_DIR)
    vega = VEGA_SED if os.path.isfile(VEGA_SED) else None
    flt = load_filters_from_instrument_dir(JOHNSON_DIR, vega_sed_path=vega)
    names = [f.name for f in flt]
    truth = dict(teff=7000.0, logg=4.0, meta=-0.5)
    r = run_forward(R=R_SUN, d=D_10PC, grid=grid, filters=flt, **truth)
    rng = np.random.default_rng(11)
    err = np.full(len(names), 0.02)
    obs = np.array([r.magnitudes[n] for n in names]) + rng.normal(0.0, err)

    res = run_inverse(obs, err, names, R=R_SUN, d=D_10PC, grid=grid, filters=flt,
                      n_walkers=16, n_steps=400, n_burn=150, progress=False,
                      verbose=False, seed=5)
    med = np.median(res.samples, axis=0)
    assert abs(med[0] - truth["teff"]) < 500.0
    assert abs(med[1] - truth["logg"]) < 0.75
    assert abs(med[2] - truth["meta"]) < 0.75
