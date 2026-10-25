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

# ---------------------------------------------------------------------------
# MESA-parity regressions (found by sweeping against the MESA colors module)
# ---------------------------------------------------------------------------

def test_bolometric_nonpositive_flux_uses_mesa_sentinel():
    """MESA's calculate_bolometric_phot returns Mag_bol = 99.0 for a SED that
    integrates to zero (e.g. an empty grid node)."""
    wave = np.array([1.0, 2.0, 3.0])
    bol_flux, bol_mag, ierr = cc.bolometric(wave, np.zeros(3))
    assert bol_flux == 0.0
    assert bol_mag == 99.0
    assert ierr != 0


def test_energy_counting_filter_matches_mesa():
    """colors/test check_tagged_energy_filter: a flat filter on [1, 3, 5] and
    F(lambda) = lambda gives <F> = 3 with energy weighting and a larger
    <F> = 3.44 with photon weighting (so a smaller magnitude)."""
    sed_wave = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    filt_wave = np.array([1.0, 3.0, 5.0])
    trans = np.ones(3)

    e_mag, e_flux, e_ierr = cc.synthetic_magnitude(
        sed_wave, sed_wave, filt_wave, trans, 1.0, photon_counting=0)
    p_mag, p_flux, p_ierr = cc.synthetic_magnitude(
        sed_wave, sed_wave, filt_wave, trans, 1.0, photon_counting=1)
    d_mag, _, _ = cc.synthetic_magnitude(sed_wave, sed_wave, filt_wave, trans, 1.0)

    assert e_ierr == 0 and p_ierr == 0
    assert e_mag == pytest.approx(-2.5 * np.log10(3.0), abs=1e-10)
    assert p_mag < e_mag
    assert d_mag == p_mag            # photon counting stays the default


def test_energy_counting_zero_points():
    """Zero points of an energy filter use the unweighted band average:
    ST is exactly the flat f_lambda constant, and AB differs from the
    photon-weighted value."""
    wave = np.array([1000.0, 2000.0, 3000.0, 4000.0])
    trans = np.array([0.0, 1.0, 1.0, 0.0])
    st_e, ierr = cc.st_zero_point(wave, trans, photon_counting=0)
    st_p, _ = cc.st_zero_point(wave, trans)
    ab_e, _ = cc.ab_zero_point(wave, trans, photon_counting=0)
    ab_p, _ = cc.ab_zero_point(wave, trans)
    assert ierr == 0
    assert st_e == pytest.approx(3.63e-9, rel=1e-13)
    assert st_p == pytest.approx(3.63e-9, rel=1e-13)
    f_ab = 3.631e-20 * 2.99792458e18 / wave ** 2
    expect = cc.simpson(wave, f_ab * trans) / cc.simpson(wave, trans)
    assert ab_e != ab_p
    assert ab_e == pytest.approx(expect, rel=1e-12)


_DET_BODY = "1000,0\n2000,1\n3000,0.5\n4000,0\n"


def _load_one(tmp_path, name, text, **kw):
    from sed_model.filters import load_filters
    (tmp_path / name).write_text(text)
    (f,) = load_filters([tmp_path / name], **kw)
    return f


def test_detector_type_read_from_column_row(tmp_path):
    """DetectorType on the column row (SED_Tools >= 0.4.0), case-insensitive,
    anywhere among the ';'-separated fields; a missing tag is photon."""
    import warnings
    warnings.simplefilter("error")
    plain = _load_one(tmp_path, "plain.dat", "Wavelength,Transmission\n" + _DET_BODY)
    photon = _load_one(tmp_path, "photon.dat",
        "Wavelength,Transmission,FilterID=A%2Fb;DetectorType=photon;WavelengthUnit=Angstrom\n" + _DET_BODY)
    energy = _load_one(tmp_path, "energy.dat",
        "Wavelength,Transmission,FilterID=A%2Fb;DetectorType=Energy;WavelengthUnit=Angstrom\n" + _DET_BODY)
    assert plain.detector_type == photon.detector_type == "photon"
    assert energy.detector_type == "energy"
    assert plain.transmission.tolist() == photon.transmission.tolist()
    assert plain.photon_counting and photon.photon_counting and energy.photon_counting


def test_energy_curve_is_converted_like_mesa(tmp_path):
    """An energy-counting curve becomes T/lambda, rescaled to the old peak
    (MESA photon_counting_filters load_filter), and its zero points equal the
    unweighted energy integrals of the file's own curve."""
    energy = _load_one(tmp_path, "e.dat",
        "Wavelength,Transmission,DetectorType=energy\n" + _DET_BODY)
    wave = np.array([1000.0, 2000.0, 3000.0, 4000.0])
    raw = np.array([0.0, 1.0, 0.5, 0.0])
    conv = raw / wave
    conv *= raw.max() / conv.max()
    np.testing.assert_allclose(energy.transmission, conv, rtol=1e-14)
    assert energy.transmission.max() == pytest.approx(raw.max(), rel=1e-14)

    ab_e, _ = cc.ab_zero_point(wave, raw, photon_counting=0)
    st_e, _ = cc.st_zero_point(wave, raw, photon_counting=0)
    assert energy.ab_zero_point == pytest.approx(ab_e, rel=1e-12)
    assert energy.st_zero_point == pytest.approx(st_e, rel=1e-12)
    photon = _load_one(tmp_path, "p.dat", "Wavelength,Transmission\n" + _DET_BODY)
    assert energy.ab_zero_point != photon.ab_zero_point


def test_detector_type_comment_lines_are_ignored(tmp_path):
    """Only the column row carries metadata; leading '# DetectorType:' lines
    (SED_Tools 0.3.0) are plain comments."""
    f = _load_one(tmp_path, "l.dat",
        "# FilterID: X/y\n# DetectorType: energy\nWavelength,Transmission\n" + _DET_BODY)
    assert f.detector_type == "photon"
    g = _load_one(tmp_path, "b.dat",
        "# DetectorType: energy\nWavelength,Transmission,DetectorType=photon\n" + _DET_BODY)
    assert g.detector_type == "photon"


def test_unrecognised_detector_type_warns(tmp_path):
    with pytest.warns(UserWarning, match="unknown DetectorType"):
        f = _load_one(tmp_path, "u.dat",
            "Wavelength,Transmission,DetectorType=bolometer\n" + _DET_BODY)
    assert f.detector_type == "photon"


def test_wavelength_unit_is_converted_to_angstrom(tmp_path):
    """WavelengthUnit names the unit of the file's own wavelength column."""
    ang = _load_one(tmp_path, "a.dat",
        "Wavelength,Transmission,WavelengthUnit=Angstrom\n" + _DET_BODY)
    nm = _load_one(tmp_path, "n.dat",
        "Wavelength,Transmission,WavelengthUnit=nm\n100,0\n200,1\n300,0.5\n400,0\n")
    assert nm.wavelengths.tolist() == ang.wavelengths.tolist()
    assert nm.ab_zero_point == pytest.approx(ang.ab_zero_point, rel=1e-13)
    with pytest.raises(ValueError, match="WavelengthUnit"):
        _load_one(tmp_path, "x.dat",
            "Wavelength,Transmission,WavelengthUnit=parsec\n" + _DET_BODY)


def test_nearest_grid_distance_matches_mesa_interp_rad():
    """Interp_rad is the normalised distance to the nearest *tabulated model*
    (rows of lookup_table.csv), with axes normalised by the range of the
    tabulated models; empty cube nodes do not count."""
    from pathlib import Path
    from sed_model.grid import AtmosphereGrid, _tabulated_nodes

    T = np.array([4000.0, 5000.0, 6000.0])
    L = np.array([1.0, 2.0])
    M = np.array([0.0])
    flux = np.ones((3, 2, 1, 4))
    flux[2, 1, 0, :] = 0.0                      # corner (6000, 2.0) has no model
    wave = np.arange(1.0, 5.0)
    g = AtmosphereGrid(T, L, M, wave, flux, np.asfortranarray(flux), Path("."),
                       lookup_nodes=_tabulated_nodes(T, L, M, flux))

    # [M/H] axis is degenerate (all zero): ignored.  Nearest model to the
    # empty corner is (5000, 2.0) or (6000, 1.0): both 0.5 away in Teff / 1.0 in logg.
    d_hole = g.nearest_grid_distance(6000.0, 2.0, 0.0)
    assert d_hole == pytest.approx(0.5, abs=1e-14)          # 1000 K / 2000 K range
    assert g.nearest_grid_distance(5000.0, 1.0, 0.0) == pytest.approx(0.0, abs=1e-14)
    # off-grid queries are not clamped
    assert g.nearest_grid_distance(7000.0, 1.0, 0.0) == pytest.approx(0.5, abs=1e-14)


def _hole_grid():
    from pathlib import Path
    from sed_model.grid import AtmosphereGrid, _tabulated_nodes

    T = np.array([4000.0, 5000.0, 6000.0])
    L = np.array([1.0, 2.0, 3.0])
    M = np.array([0.0])                      # degenerate axis
    wave = np.linspace(1000.0, 9000.0, 8)
    flux = np.ones((3, 3, 1, 8))
    flux[2, 2, 0, :] = 0.0                   # no model at (6000, 3.0)
    return AtmosphereGrid(T, L, M, wave, flux, np.asfortranarray(flux), Path("."),
                          lookup_nodes=_tabulated_nodes(T, L, M, flux))


def test_cell_has_empty_nodes():
    g = _hole_grid()
    assert g.cell_has_empty_nodes(5500.0, 2.5, 0.0)        # cell touches the hole
    assert g.cell_has_empty_nodes(6000.0, 3.0, 0.0)        # the empty node itself
    assert not g.cell_has_empty_nodes(4500.0, 1.5, 0.0)    # clean cell
    # on a populated node the empty neighbour carries no weight
    assert not g.cell_has_empty_nodes(5000.0, 2.0, 0.0)
    assert not g.cell_has_empty_nodes(6000.0, 2.0, 0.0)
    # out-of-range queries use the boundary node/cell
    assert g.cell_has_empty_nodes(9000.0, 9.0, 0.0)
    assert not g.cell_has_empty_nodes(1000.0, 0.0, 0.0)


def test_run_forward_warns_in_empty_cell():
    from sed_model import run_forward

    g = _hole_grid()
    kw = dict(R=6.957e10, d=3.0857e19, grid=g, filters=[])
    with pytest.warns(UserWarning, match="empty nodes"):
        run_forward(teff=5500.0, logg=2.5, meta=0.0, **kw)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        run_forward(teff=4500.0, logg=1.5, meta=0.0, **kw)


def test_linear_zero_bin_fallback_is_not_clamped():
    """A populated on-node, in-bounds query is not 'clamped' just because some
    wavelength bins hold zero flux and fell back to the nearest node."""
    from sed_model import run_forward

    g = _hole_grid()
    g.flux_fortran[1, 1, 0, 0] = 0.0         # one zero bin at a populated node
    r = run_forward(teff=5000.0, logg=2.0, meta=0.0, R=6.957e10, d=3.0857e19,
                    grid=g, filters=[], interp_method="linear")
    assert not r.clamped
