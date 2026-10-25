"""
tests/test_kernels.py
=====================
Data-free checks of the Fortran kernels against the MESA colors unit test
(colors/test/src/test_colors.f90).  Each case builds a tiny synthetic grid or
SED, so no SED_Tools data is required -- only the compiled cc_api extension.

Usage
-----
    pytest tests/test_kernels.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

try:
    from sed_model._cc_ext import get_cc_api
    cc = get_cc_api()
except Exception:  # pragma: no cover
    cc = None

pytestmark = pytest.mark.skipif(cc is None, reason="cc_api Fortran extension not built")

LSUN = 3.828e33
PC   = 3.0856775814913673e18


def _cube(axis_values, per_lambda):
    """(nx,1,1,nlam) Fortran-ordered cube from per-wavelength columns."""
    cube = np.stack(per_lambda, axis=-1)[:, None, None, :]
    return np.asfortranarray(cube, dtype=np.float64)


def test_filter_support():
    sed_wave = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    sed_flux = np.array([100.0, 2.0, 2.0, 2.0, 100.0])
    filt_wave = np.array([2.0, 3.0, 4.0])
    trans = np.ones(3)

    mag, band_flux, ierr = cc.synthetic_magnitude(
        sed_wave, sed_flux, filt_wave, trans, 1.0)
    vega_zp, zp_ierr = cc.vega_zero_point(sed_wave, sed_flux, filt_wave, trans)

    assert ierr == 0 and zp_ierr == 0
    assert mag == pytest.approx(-0.75257498915995, abs=1e-10)
    assert vega_zp == pytest.approx(2.0, abs=1e-12)


def test_nonuniform_hermite_derivatives():
    x = np.array([0.0, 1.0, 3.0, 6.0])
    wave = np.array([1.0, 2.0, 4.0])
    cube = _cube(x, [100.0 + x * x] * 3)
    one = np.array([0.0])

    flux, ierr = cc.interp_sed_hermite(1.5, 0.0, 0.0, x, one, one, cube)
    assert ierr == 0
    np.testing.assert_allclose(flux, 102.25, rtol=1e-12)

    flux, ierr = cc.interp_sed_hermite_bounded(1.5, 0.0, 0.0, x, one, one, cube, wave)
    assert ierr == 0
    np.testing.assert_allclose(flux, 102.25, rtol=1e-12)


def test_bounded_hermite_clips_small_negative():
    x = np.array([0.0, 1.0, 2.0])
    wave = np.array([1.0, 2.0, 4.0])
    cube = _cube(x, [
        np.full(3, 100.0),
        np.array([1.0, 1.0, 17.000016]),
        np.full(3, 100.0),
    ])
    one = np.array([0.0])

    raw, _ = cc.interp_sed_hermite(0.5, 0.0, 0.0, x, one, one, cube)
    assert raw[1] < 0.0          # unbounded Hermite undershoots

    flux, ierr = cc.interp_sed_hermite_bounded(0.5, 0.0, 0.0, x, one, one, cube, wave)
    assert ierr == 0
    np.testing.assert_allclose(flux, [100.0, 0.0, 100.0], atol=1e-12)


def test_bounded_hermite_falls_back_to_linear():
    x = np.array([0.0, 1.0, 2.0])
    wave = np.array([1.0, 2.0, 4.0])
    cube = _cube(x, [np.array([1.0, 1.0, 100.0])] * 3)
    one = np.array([0.0])

    flux, ierr = cc.interp_sed_hermite_bounded(0.5, 0.0, 0.0, x, one, one, cube, wave)
    assert ierr == 1
    np.testing.assert_allclose(flux, 1.0, rtol=1e-12)


def test_bolometric_solar_zero_point():
    solar_flux = LSUN / (4.0 * np.pi * (10.0 * PC) ** 2)
    wave = np.array([1.0, 2.0, 3.0])
    flux = np.full(3, solar_flux / 2.0)

    bol_flux, bol_mag, ierr = cc.bolometric(wave, flux)
    assert ierr == 0
    assert bol_mag == pytest.approx(4.74, abs=1e-10)


def test_bolometric_ignores_unphysical_samples():
    wave = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    good = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
    bad = good.copy()
    bad[2] = -5.0                # negative sample is zeroed, as in MESA
    f_good, _, _ = cc.bolometric(wave, good)
    f_bad, _, ierr = cc.bolometric(wave, bad)
    assert ierr == 0
    assert 0.0 < f_bad < f_good


def test_simpson_nonuniform_exact_for_quadratic():
    x = np.array([0.0, 1.0, 3.0, 4.0, 7.0])
    result = cc.simpson(x, x * x)
    assert result == pytest.approx(7.0 ** 3 / 3.0, rel=1e-12)


def test_bolometric_zeroes_nan_and_inf_samples():
    """NaN/Inf flux samples must be treated like MESA does (zeroed), under -O2."""
    wave = np.arange(1.0, 8.0)
    clean = np.ones(7)
    clean[[2, 4]] = 0.0
    dirty = np.ones(7)
    dirty[2] = np.nan
    dirty[4] = np.inf
    f_clean, _, ierr_clean = cc.bolometric(wave, clean)
    f_dirty, _, ierr_dirty = cc.bolometric(wave, dirty)
    assert ierr_clean == 0 and ierr_dirty == 0
    assert np.isfinite(f_dirty)
    assert f_dirty == pytest.approx(f_clean, rel=1e-12)


def test_linear_interp_handles_nan_cube_node():
    """A NaN node must trigger the nearest-node fallback, not propagate."""
    x = np.array([0.0, 1.0])
    one = np.array([0.0])
    cube = _cube(x, [np.array([1.0, np.nan]), np.array([2.0, 4.0])])
    flux, ierr = cc.interp_sed_linear(0.4, 0.0, 0.0, x, one, one, cube)
    assert np.all(np.isfinite(flux))
    assert ierr >= 1


def test_filter_header_row_with_metadata_field(tmp_path):
    """SED_Tools >= 0.4.0 appends ``;``-separated metadata as a third comma
    field of the column row.  The loader reads only the two data columns, for
    that layout, the bare row, and the older leading ``# Key: value`` lines."""
    from sed_model.filters import load_filters

    rows = "3200.0,2.4e-08\n3300.0,0.5\n3400.0,0.9\n3500.0,0.0\n"
    layouts = {
        "bare": "Wavelength,Transmission\n" + rows,
        "meta": ("Wavelength,Transmission,FilterID=GAIA/GAIA3.G;DetectorType=photon;"
                 "DetectorTypeSource=svo-curve-match;WavelengthUnit=Angstrom;"
                 "Components=filter+qe+instrument\n" + rows),
        "comments": "# FilterID: GAIA/GAIA3.G\n# DetectorType: photon\nWavelength,Transmission\n" + rows,
    }
    waves = {}
    for name, text in layouts.items():
        (tmp_path / f"{name}.dat").write_text(text)
        (f,) = load_filters([tmp_path / f"{name}.dat"])
        waves[name] = (f.wavelengths.tolist(), f.transmission.tolist())
    assert waves["meta"] == waves["bare"] == waves["comments"]
    assert waves["meta"][0] == [3200.0, 3300.0, 3400.0, 3500.0]
