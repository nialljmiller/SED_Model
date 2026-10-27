"""
tests/test_extinction.py
========================
Extinction laws against a frozen reference.

``tests/data/extinction_reference.npz`` holds k = A(lambda)/A(V) from the MESA
colors extinction module (itself validated against the dust_extinction
package to ~1e-14) for all seven MESA laws at R_V = 3.1 and 5.0 on 100
log-spaced wavelengths from 500 A to 50000 A, so both the fitted range and
the out-of-range behaviour are covered.  No external data is needed except
for the forward-model tests at the end, which use the SED_Tools grid.
"""

from __future__ import annotations

import os
import warnings

import numpy as np
import pytest

from _data_paths import KURUCZ_DIR, JOHNSON_DIR, VEGA_SED
from sed_model import sed_extinction as se

REF = np.load(os.path.join(os.path.dirname(__file__), "data", "extinction_reference.npz"))
WAVE = REF["wave"]

LAWS = {
    "CCM89":      ("ccm89", True),
    "O94":        ("odonnell94", True),
    "F99":        ("fitzpatrick99", True),
    "FM07":       ("fm07", False),
    "G23":        ("gordon23", True),
    "G03_SMCBAR": ("g03_smcbar", False),
    "G03_LMCAVG": ("g03_lmcavg", False),
}
# blue limit of each law (A); shortward k is held constant
BLUE_LIMIT = {"CCM89": 1000.0, "O94": 1000.0, "F99": 1000.0, "FM07": 910.0,
              "G23": 912.0, "G03_SMCBAR": 1000.0, "G03_LMCAVG": 1000.0}


def _k(key, rv):
    law, _ = LAWS[key]
    model = se.ExtinctionModel(enabled=True, law=law, a_v=1.0, r_v=rv)
    return model.extinction_curve(WAVE)


@pytest.mark.parametrize("rv", ["3.10", "5.00"])
@pytest.mark.parametrize("key", list(LAWS))
def test_matches_mesa_reference(key, rv):
    k = _k(key, float(rv))
    ref = REF[f"{key}_{rv}"]
    np.testing.assert_allclose(k, ref, rtol=0, atol=1e-12)


@pytest.mark.parametrize("key", list(LAWS))
def test_reference_spans_in_and_out_of_range(key):
    assert (WAVE < BLUE_LIMIT[key]).any() and (WAVE > 33333.0).any()


@pytest.mark.parametrize("key", list(LAWS))
def test_blue_end_is_held_constant(key):
    law, _ = LAWS[key]
    w = np.array([BLUE_LIMIT[key], 0.8 * BLUE_LIMIT[key], 500.0, 100.0])
    k = se.ExtinctionModel(enabled=True, law=law, a_v=1.0).extinction_curve(w)
    np.testing.assert_allclose(k, k[0], rtol=0, atol=1e-12)


def test_non_positive_wavelength_sent_to_far_uv():
    k = se.ccm89(np.array([0.0, -5.0, 100.0]), 1.0)
    np.testing.assert_allclose(k, k[2], atol=1e-12)


def test_scalar_and_shape():
    assert np.ndim(se.fitzpatrick99(5500.0, 1.0)) == 0
    w = np.linspace(2000.0, 20000.0, 12).reshape(3, 4)
    assert se.ccm89(w, 1.0).shape == (3, 4)


def test_a_v_scales_linearly_and_v_band_is_one():
    for law in ("ccm89", "odonnell94", "fitzpatrick99", "fm07", "gordon23",
                "g03_smcbar", "g03_lmcavg"):
        m1 = se.ExtinctionModel(enabled=True, law=law, a_v=1.0)
        m2 = se.ExtinctionModel(enabled=True, law=law, a_v=0.37)
        w = np.array([5500.0])
        np.testing.assert_allclose(m2.extinction_curve(w), 0.37 * m1.extinction_curve(w))
        assert abs(m1.extinction_curve(w)[0] - 1.0) < 0.03


def test_g23_optical_and_ir_regression():
    # The pre-fix gordon23 gave 1.31 at 5500 A and 0.86 at 2.2 micron
    k = se.ExtinctionModel(enabled=True, law="gordon23", a_v=1.0).extinction_curve(
        np.array([5500.0, 22000.0]))
    assert abs(k[0] - 1.0) < 0.02
    assert abs(k[1] - 0.10) < 0.03


# ---------------------------------------------------------------------------
# R_V handling
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("law,lo,hi", [
    ("ccm89", 2.0, 6.0), ("odonnell94", 2.0, 6.0),
    ("fitzpatrick99", 2.0, 6.0), ("gordon23", 2.3, 5.6),
])
def test_rv_range_enforced(law, lo, hi):
    for rv in (lo, hi, 0.5 * (lo + hi)):
        se.ExtinctionConfig(law=law, r_v=rv)            # ok at the limits
    for rv in (lo - 0.01, hi + 0.01, 0.0, -3.0, float("nan")):
        with pytest.raises(ValueError, match="R_V"):
            se.ExtinctionConfig(law=law, r_v=rv)
        with pytest.raises(ValueError, match="R_V"):
            getattr(se, law)(np.array([5500.0]), 1.0, rv)
    with pytest.raises(ValueError, match=f"{lo:g} <= R_V <= {hi:g}"):
        se.ExtinctionConfig(law=law, r_v=hi + 1.0)


def test_rv_error_via_replace_and_factory():
    with pytest.raises(ValueError, match="R_V"):
        se.make_extinction_model(enabled=True, law="ccm89", a_v=0.1, r_v=7.0)
    ok = se.make_extinction_model(enabled=True, law="ccm89", a_v=0.1, r_v=3.1)
    assert ok.with_av(0.5).config.a_v == 0.5


@pytest.mark.parametrize("law", ["fm07", "g03_smcbar", "g03_lmcavg"])
def test_fixed_shape_laws_ignore_rv(law):
    w = np.geomspace(1000.0, 30000.0, 40)
    base = se.ExtinctionModel(enabled=True, law=law, a_v=1.0, r_v=3.1).extinction_curve(w)
    for rv in (1.0, 5.0, 50.0):
        m = se.ExtinctionModel(enabled=True, law=law, a_v=1.0, r_v=rv)
        np.testing.assert_array_equal(m.extinction_curve(w), base)


def test_negative_av_and_unknown_law():
    with pytest.raises(ValueError):
        se.ExtinctionConfig(law="ccm89", a_v=-0.1)
    with pytest.raises(ValueError):
        se.ExtinctionConfig(law="nope")


def test_new_laws_registered():
    for law in ("g03_smcbar", "g03_lmcavg"):
        assert law in se.AVAILABLE_LAWS and law in se.__all__


def test_calzetti_unchanged():
    k = se.calzetti00(np.array([5500.0]), 1.0)
    assert 0.9 < k[0] < 1.1


def test_gordon23_legacy_environments_warn():
    with pytest.warns(UserWarning, match="could not be verified"):
        k = se.gordon23(np.array([4000.0]), 1.0, environment="lmc")
    assert np.isfinite(k).all()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        se.gordon23(np.array([4000.0]), 1.0, environment="mw")
    with pytest.raises(ValueError):
        se.gordon23(np.array([4000.0]), 1.0, environment="xyz")


# ---------------------------------------------------------------------------
# apply / remove
# ---------------------------------------------------------------------------

def test_apply_remove_round_trip():
    w = np.geomspace(1000.0, 40000.0, 50)
    flux = np.linspace(1.0, 3.0, 50)
    m = se.ExtinctionModel(enabled=True, law="fitzpatrick99", a_v=0.8, r_v=3.1)
    red = m.apply(w, flux)
    assert np.all(red < flux)
    np.testing.assert_allclose(m.remove(w, red), flux, rtol=1e-13)
    a = m.extinction_curve(w)
    np.testing.assert_allclose(red, se.apply_extinction(a, flux))
    np.testing.assert_allclose(se.remove_extinction(a, red), flux, rtol=1e-13)


def test_disabled_model_is_identity():
    w = np.array([4000.0, 6000.0])
    f = np.array([1.0, 2.0])
    m = se.ExtinctionModel.disabled()
    np.testing.assert_array_equal(m.apply(w, f), f)
    np.testing.assert_array_equal(m.extinction_curve(w), 0.0)


# ---------------------------------------------------------------------------
# Forward model: bolometric quantities are intrinsic
# ---------------------------------------------------------------------------

def _cc_api_available():
    try:
        from sed_model import cc_api  # noqa: F401
        return True
    except ImportError:
        return False


requires_data = pytest.mark.skipif(
    not (_cc_api_available() and os.path.isdir(KURUCZ_DIR) and os.path.isdir(JOHNSON_DIR)),
    reason="needs the cc_api extension, Kurucz2003all grid and Johnson filters",
)


@requires_data
def test_forward_bolometric_intrinsic_bands_reddened():
    from sed_model import (load_grid, load_filters_from_instrument_dir,
                           run_forward, make_extinction_model)
    grid = load_grid(KURUCZ_DIR)
    vega = VEGA_SED if os.path.isfile(VEGA_SED) else None
    filters = load_filters_from_instrument_dir(JOHNSON_DIR, vega_sed_path=vega)
    args = dict(teff=5778.0, logg=4.44, meta=0.0, R=6.957e10, d=3.0857e19,
                grid=grid, filters=filters, mag_system="AB")

    plain = run_forward(**args)
    for law in ("fitzpatrick99", "gordon23", "g03_lmcavg"):
        red = run_forward(**args, extinction=make_extinction_model(
            enabled=True, law=law, a_v=1.0))
        assert red.extinction_applied
        assert red.bol_flux == plain.bol_flux
        assert red.bol_mag == plain.bol_mag
        np.testing.assert_array_equal(red.surface_flux, plain.surface_flux)
        assert np.all(red.observed_flux <= plain.observed_flux)
        for name in plain.magnitudes:
            assert red.magnitudes[name] > plain.magnitudes[name]
        # V-band dimming close to A_V = 1 for these laws
        assert 0.8 < red.magnitudes["V"] - plain.magnitudes["V"] < 1.2
