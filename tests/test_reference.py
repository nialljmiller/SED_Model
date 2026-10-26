"""
tests/test_reference.py
=======================
Frozen cross-code reference tests: SED_Model against MESA ``colors`` and MSG.

The reference values live in ``tests/data/reference_points.npz`` (about 30 KB).
They were produced once, outside the repository, during the 0.4.x validation:

* **MESA** -- a standalone driver around the real MESA ``colors`` sources,
  Kurucz2003all grid, ``Generic/Johnson`` and ``GAIA/GAIA`` filters, Vega / AB /
  ST systems, all three interpolation methods.  33 queries cover on-node,
  off-node, node +- epsilon, cell-edge, cell-face, corner, grid-edge and
  out-of-range (clamped) positions.  Magnitudes, bolometric magnitude/flux and
  zero points must match to double-precision rounding (``TOL_MESA``).
* **MSG** -- v2.2.1 on a grid built from the same Kurucz cube, 24 MIST
  isochrone points spanning 3500-30000 K.  These are *not* bit-for-bit
  comparisons (MSG's quadrature and cubic differ), so the documented
  tolerances apply: a few mmag in V, R, I and Gaia, up to ~13 mmag in the
  coarsely sampled Johnson U and B (see docs/testing.md).

Data needs: the compiled extension and the Kurucz2003all grid (plus the
Johnson / Gaia filter directories and the Vega SED); each test skips with a
reason if any is missing.

Usage (run from outside the source tree against an installed build)
-----
    python -I -m pytest tests/test_reference.py -v
"""

from __future__ import annotations

import os
import warnings

import numpy as np
import pytest

from _data_paths import KURUCZ_DIR, JOHNSON_DIR, VEGA_SED

TOL_MESA = 1e-10          # mag; observed agreement is ~3e-13
TOL_MSG = 6e-3            # mag; V, R, I, Gaia (documented 1-4 mmag, cool stars reach ~5)
TOL_MSG_UB = 1.5e-2       # mag; Johnson U, B (documented up to 13 mmag)
PC_CM = 3.0856775814913673e19   # 10 pc

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "reference_points.npz")
GAIA_DIR = os.path.join(os.path.dirname(os.path.dirname(JOHNSON_DIR)), "GAIA", "GAIA")


def _cc_api_available():
    try:
        from sed_model import cc_api  # noqa: F401
        return True
    except ImportError:
        return False


def _skip_reason():
    if not _cc_api_available():
        return "cc_api Fortran extension not built"
    if not os.path.isdir(KURUCZ_DIR):
        return f"Kurucz2003all grid not found at {KURUCZ_DIR}"
    if not os.path.isfile(VEGA_SED):
        return f"Vega SED not found at {VEGA_SED}"
    if not os.path.isdir(JOHNSON_DIR):
        return f"Johnson filters not found at {JOHNSON_DIR}"
    if not os.path.isdir(GAIA_DIR):
        return f"Gaia filters not found at {GAIA_DIR}"
    if not os.path.isfile(DATA):
        return f"reference data not found at {DATA}"
    return None


_REASON = _skip_reason()
pytestmark = pytest.mark.skipif(_REASON is not None, reason=_REASON or "")

REF = np.load(DATA) if os.path.isfile(DATA) else None
N_MESA = len(REF["mesa_q"]) if REF is not None else 0
METHODS = list(REF["mesa_methods"]) if REF is not None else []
SYSTEMS = list(REF["mesa_systems"]) if REF is not None else []


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def grid():
    from sed_model import load_grid
    return load_grid(KURUCZ_DIR)


@pytest.fixture(scope="module")
def filter_sets():
    from sed_model import load_filters_from_instrument_dir
    return {
        "johnson": load_filters_from_instrument_dir(JOHNSON_DIR, vega_sed_path=VEGA_SED),
        "gaia": load_filters_from_instrument_dir(GAIA_DIR, vega_sed_path=VEGA_SED),
    }


def _forward(grid, filters, q, R, d, system, method):
    from sed_model import run_forward
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")   # empty-cell warnings are not under test here
        return run_forward(teff=float(q[0]), logg=float(q[1]), meta=float(q[2]),
                           R=float(R), d=float(d), grid=grid, filters=filters,
                           mag_system=system, interp_method=method)


# ---------------------------------------------------------------------------
# MESA colors: double-precision agreement
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tag", ["johnson", "gaia"])
def test_mesa_zero_points(filter_sets, tag):
    filters = filter_sets[tag]
    names = [str(n) for n in REF[f"{tag}_filters"]]
    assert sorted(f.name for f in filters) == sorted(names)
    byname = {f.name: f for f in filters}
    zp = REF[f"{tag}_zp"]
    for i, name in enumerate(names):
        f = byname[name]
        got = np.array([f.vega_zero_point, f.ab_zero_point, f.st_zero_point])
        np.testing.assert_allclose(got, zp[i], rtol=TOL_MESA, err_msg=f"{tag}/{name}")


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("i", range(N_MESA), ids=lambda i: f"{i:02d}-{REF['mesa_cat'][i]}")
def test_mesa_magnitudes(grid, filter_sets, i, method):
    mi = METHODS.index(method)
    q, R, d = REF["mesa_q"][i], REF["mesa_R"][i], REF["mesa_d"][i]
    for tag in ("johnson", "gaia"):
        names = [str(n) for n in REF[f"{tag}_filters"]]
        for si, system in enumerate(SYSTEMS):
            res = _forward(grid, filter_sets[tag], q, R, d, str(system), method)
            got = np.array([res.magnitudes[n] for n in names])
            ref = REF[f"{tag}_mags"][i, mi, si]
            np.testing.assert_allclose(
                got, ref, rtol=0, atol=TOL_MESA,
                err_msg=f"{tag}/{system}/{method} q={q} cat={REF['mesa_cat'][i]}")
            if si == 0:
                assert abs(res.bol_mag - REF[f"{tag}_bol_mag"][i, mi]) < TOL_MESA
                np.testing.assert_allclose(res.bol_flux, REF[f"{tag}_bol_flux"][i, mi], rtol=1e-10)


def test_mesa_reference_covers_categories():
    """Guard the data file itself: the advertised position classes are present."""
    cats = set(REF["mesa_cat"])
    for c in ("node_valid", "offnode", "edge_T", "face_T_node", "oob", "corner_exact"):
        assert c in cats, c
    assert N_MESA >= 30


# ---------------------------------------------------------------------------
# MSG: documented tolerances
# ---------------------------------------------------------------------------

MSG_BANDS = {"U": "U", "B": "B", "V": "V", "R": "R", "I": "I",
             "G": "G", "BP": "Gbp", "RP": "Grp"}
N_MSG = len(REF["msg_teff"]) if REF is not None else 0


@pytest.mark.parametrize("order,method,suffix", [(3, "hermite", ""), (1, "linear", "_o1")],
                         ids=["cubic-vs-hermite", "linear-vs-linear"])
@pytest.mark.parametrize("i", range(N_MSG), ids=lambda i: f"T{REF['msg_teff'][i]:.0f}")
def test_msg_photometry(grid, filter_sets, i, order, method, suffix):
    filters = filter_sets["johnson"] + filter_sets["gaia"]
    q = (REF["msg_teff"][i], REF["msg_logg"][i], REF["msg_feh"][i])
    res = _forward(grid, filters, q, REF["msg_R_cm"][i], PC_CM, "Vega", method)
    for msg_name, name in MSG_BANDS.items():
        tol = TOL_MSG_UB if msg_name in ("U", "B") else TOL_MSG
        ref = REF["msg_" + msg_name + suffix][i]
        assert abs(res.magnitudes[name] - ref) < tol, (
            f"{name}: SED_Model {res.magnitudes[name]:.5f} vs MSG {ref:.5f} (order {order}) at {q}")
