"""
tests/test_filters.py
=====================
Data-free checks of the filter reader (sed_model.filters): the SED_Tools
column-row metadata (``Wavelength,Transmission,Key=value;...``), DetectorType
handling, energy -> photon conversion, and WavelengthUnit conversion.  Every
case writes a tiny synthetic filter file into ``tmp_path``.

Usage
-----
    pytest tests/test_filters.py -v
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from sed_model.filters import load_filters

WAVE = np.array([4000.0, 4500.0, 5000.0, 5500.0, 6000.0])
TRANS = np.array([0.1, 0.5, 1.0, 0.6, 0.2])


def _write(path, column_row, wave=WAVE, trans=TRANS, comments=()):
    lines = list(comments)
    if column_row is not None:
        lines.append(column_row)
    lines += [f"{float(w)!r},{float(t)!r}" for w, t in zip(wave, trans)]
    path.write_text("\n".join(lines) + "\n")
    return path


def _load_one(path):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        return load_filters([path])[0]


def _expected_energy(wave, trans):
    out = trans / wave
    return out * (trans.max() / out.max())


# -- metadata / DetectorType ------------------------------------------------

@pytest.mark.parametrize("row", [
    None,
    "Wavelength,Transmission",
    "Wavelength,Transmission,FilterID=GAIA/GAIA3.G;Facility=GAIA",
    "Wavelength,Transmission,DetectorType=photon",
    "Wavelength,Transmission,FilterID=X;DetectorType=Photon;WavelengthUnit=Angstrom",
])
def test_photon_is_unchanged(tmp_path, row):
    f = _load_one(_write(tmp_path / "G.dat", row))
    assert f.name == "G"
    assert f.detector_type == "photon"
    assert f.photon_counting is True
    np.testing.assert_allclose(f.wavelengths, WAVE)
    np.testing.assert_allclose(f.transmission, TRANS)


def test_detector_type_not_first_field_and_case_insensitive(tmp_path):
    row = "Wavelength,Transmission,FilterID=A/B;DetectorType=ENERGY;Other=1"
    f = _load_one(_write(tmp_path / "E.dat", row))
    assert f.detector_type == "energy"


def test_unknown_detector_type_warns_and_is_photon(tmp_path):
    p = _write(tmp_path / "U.dat", "Wavelength,Transmission,DetectorType=bolometer")
    with pytest.warns(UserWarning, match="unknown DetectorType 'bolometer'"):
        f = load_filters([p])[0]
    assert f.detector_type == "photon"
    np.testing.assert_allclose(f.transmission, TRANS)


def test_percent_encoded_value_is_decoded(tmp_path):
    p = _write(tmp_path / "P.dat", "Wavelength,Transmission,DetectorType=en%65rgy")
    assert _load_one(p).detector_type == "energy"


@pytest.mark.parametrize("comment", [
    "# DetectorType: energy",
    "# DetectorType = energy",
    "# DetectorType: bogus",
    "#DetectorType=energy",
])
def test_legacy_comment_lines_are_ignored_silently(tmp_path, comment):
    # Only the column-row format is supported; comment tags are not read.
    p = _write(tmp_path / "L.dat", "Wavelength,Transmission", comments=[comment])
    f = _load_one(p)  # warnings are errors here
    assert f.detector_type == "photon"
    np.testing.assert_allclose(f.transmission, TRANS)


def test_legacy_comment_ignored_without_column_row(tmp_path):
    p = _write(tmp_path / "L2.dat", None, comments=["# DetectorType: energy"])
    f = _load_one(p)
    assert f.detector_type == "photon"
    np.testing.assert_allclose(f.transmission, TRANS)


def test_comment_does_not_override_column_row(tmp_path):
    p = _write(tmp_path / "L3.dat", "Wavelength,Transmission,DetectorType=energy",
               comments=["# DetectorType: photon"])
    assert _load_one(p).detector_type == "energy"


# -- energy -> photon conversion --------------------------------------------

def test_energy_curve_converted_and_peak_preserved(tmp_path):
    p = _write(tmp_path / "E.dat", "Wavelength,Transmission,DetectorType=energy")
    f = _load_one(p)
    assert f.detector_type == "energy"
    assert f.photon_counting is True
    np.testing.assert_allclose(f.wavelengths, WAVE)
    np.testing.assert_allclose(f.transmission, _expected_energy(WAVE, TRANS))
    assert f.transmission.max() == pytest.approx(TRANS.max())
    # T/lambda tilts the curve toward the blue relative to the input
    assert not np.allclose(f.transmission, TRANS)


def test_energy_zero_points_match_energy_weighted_band_average(tmp_path):
    # With the converted curve, the lambda-weighted integrals equal
    # int F T_E dl / int T_E dl.
    e = _load_one(_write(tmp_path / "E.dat",
                         "Wavelength,Transmission,DetectorType=energy"))
    w = np.linspace(WAVE[0], WAVE[-1], 4001)
    te = np.interp(w, WAVE, TRANS)
    f_ab = 3.631e-20 * 2.99792458e18 / w**2
    expected = np.trapezoid(f_ab * te, w) / np.trapezoid(te, w)
    assert e.ab_zero_point == pytest.approx(expected, rel=1e-2)
    assert e.st_zero_point == pytest.approx(3.63e-9, rel=1e-6)


def test_energy_and_photon_zero_points_differ(tmp_path):
    e = _load_one(_write(tmp_path / "E.dat",
                         "Wavelength,Transmission,DetectorType=energy"))
    p = _load_one(_write(tmp_path / "P.dat",
                         "Wavelength,Transmission,DetectorType=photon"))
    assert e.ab_zero_point != pytest.approx(p.ab_zero_point, rel=1e-4)


# -- WavelengthUnit ----------------------------------------------------------

@pytest.mark.parametrize("unit, factor", [
    ("Angstrom", 1.0), ("nm", 10.0), ("nanometer", 10.0),
    ("um", 1.0e4), ("micron", 1.0e4), ("µm", 1.0e4), ("mm", 1.0e7), ("m", 1.0e10),
])
def test_wavelength_unit_converted_to_angstrom(tmp_path, unit, factor):
    wave = WAVE / factor
    p = _write(tmp_path / "W.dat", f"Wavelength,Transmission,WavelengthUnit={unit}",
               wave=wave)
    f = _load_one(p)
    np.testing.assert_allclose(f.wavelengths, WAVE, rtol=1e-12)
    np.testing.assert_allclose(f.transmission, TRANS)


def test_wavelength_unit_applies_before_energy_conversion(tmp_path):
    p = _write(tmp_path / "W.dat",
               "Wavelength,Transmission,DetectorType=energy;WavelengthUnit=nm",
               wave=WAVE / 10.0)
    f = _load_one(p)
    np.testing.assert_allclose(f.wavelengths, WAVE)
    np.testing.assert_allclose(f.transmission, _expected_energy(WAVE, TRANS))


def test_unknown_wavelength_unit_raises(tmp_path):
    p = _write(tmp_path / "W.dat", "Wavelength,Transmission,WavelengthUnit=parsec")
    with pytest.raises(ValueError, match="unrecognised WavelengthUnit 'parsec'"):
        load_filters([p])
