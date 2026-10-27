"""
sed_extinction.py
=================
Interstellar dust extinction for SED fitting, in pure NumPy.

The laws are a NumPy port of the extinction module of the MESA ``colors``
code (same structure, same building blocks, same handling outside each law's
fitted range).  Coefficients follow the ``dust_extinction`` package (Gordon
2024, JOSS, 9, 7023) and, for FM07, the ``extinction`` package (Barbary 2016).
k = A(lambda)/A(V) agrees with the MESA module to about 1e-12.

Usage
-----
Apply extinction to a spectrum before filter convolution::

    from sed_model import ExtinctionModel

    ext = ExtinctionModel(law='fitzpatrick99', r_v=3.1, a_v=0.3, enabled=True)
    flux_obs = ext.apply(wavelength_aa, flux_intrinsic)

Supported laws (config key; calibrated R_V range)
-------------------------------------------------
'ccm89'          Cardelli, Clayton & Mathis (1989) ApJ 345 245;  2 <= R_V <= 6
'odonnell94'     O'Donnell (1994) ApJ 422 158 (CCM89, revised optical); 2 <= R_V <= 6
'fitzpatrick99'  Fitzpatrick (1999) PASP 111 63;  2 <= R_V <= 6
'fm07'           Fitzpatrick & Massa (2007) ApJ 663 320;  R_V = 3.1 fixed
'gordon23'       Gordon et al. (2023) ApJ 950 86 (Milky Way);  2.3 <= R_V <= 5.6
'g03_smcbar'     Gordon et al. (2003) ApJ 594 279, SMC bar;  R_V = 2.74 fixed
'g03_lmcavg'     Gordon et al. (2003) ApJ 594 279, LMC average;  R_V = 3.41 fixed
'calzetti00'     Calzetti et al. (2000) ApJ 533 682 (starburst attenuation);
                 not in MESA

Fixed-shape laws ignore ``r_v``.  For the R_V-dependent laws an R_V outside
the calibrated range raises ``ValueError``.

Outside a law's fitted range, longward the law's own functional form is
continued (each tends to zero as the wavelength goes to infinity); shortward,
k is held at its value at the law's blue limit (CCM89, O94, F99, G03: 1000 A;
FM07: 910 A; G23: 912 A), because the far-UV fits diverge beyond it.

Convention
----------
All wavelengths are in Angstroms.  Fluxes are in any linear unit.  The
extinction is applied as::

    F_obs(lambda) = F_intrinsic(lambda) * 10**(-0.4 * A_V * k(lambda))

with k = A(lambda)/A(V) the normalised curve of the chosen law.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np

__all__ = [
    "ExtinctionModel",
    "ExtinctionConfig",
    "make_extinction_model",
    "apply_extinction",
    "remove_extinction",
    # individual law functions
    "ccm89",
    "odonnell94",
    "fitzpatrick99",
    "fm07",
    "calzetti00",
    "gordon23",
    "g03_smcbar",
    "g03_lmcavg",
    "AVAILABLE_LAWS",
    "RV_RANGES",
]

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

AVAILABLE_LAWS = (
    "ccm89",
    "odonnell94",
    "fitzpatrick99",
    "fm07",
    "calzetti00",
    "gordon23",
    "g03_smcbar",
    "g03_lmcavg",
)

#: R_V range each R_V-dependent law was calibrated over.  Laws not listed
#: (fm07, g03_smcbar, g03_lmcavg, calzetti00) do not validate R_V.
RV_RANGES = {
    "ccm89": (2.0, 6.0),
    "odonnell94": (2.0, 6.0),
    "fitzpatrick99": (2.0, 6.0),
    "gordon23": (2.3, 5.6),
}


def _check_rv(law: str, r_v: float) -> None:
    """Raise ValueError when r_v is outside the calibrated range of ``law``."""
    rng = RV_RANGES.get(law)
    if rng is None:
        return
    lo, hi = rng
    if not (lo <= r_v <= hi):      # also rejects NaN
        raise ValueError(
            f"{law} is calibrated for {lo:g} <= R_V <= {hi:g}, got R_V = {r_v:g}"
        )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_X_FUV = 1.0e4 / 910.0


def _wave_to_x(wave) -> np.ndarray:
    """Angstroms -> wavenumber in inverse microns.

    A non-positive wavelength is sent to x = 1e4 (far-UV), where every law is
    held constant, as in the MESA module.
    """
    wave = np.atleast_1d(np.asarray(wave, dtype=np.float64))
    x = np.full(wave.shape, 1.0e4)
    pos = wave > 0.0
    x[pos] = 1.0e4 / wave[pos]
    return x


def _polynomial(c, y):
    """sum of c[i] * y**i (Horner)."""
    out = np.full_like(y, c[-1])
    for ci in c[-2::-1]:
        out = out * y + ci
    return out


def _smoothstep(v, lo, hi):
    t = np.clip((v - lo) / (hi - lo), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _drude(x, amp, x0, fwhm):
    return amp * (fwhm / x0) ** 2 / ((x / x0 - x0 / x) ** 2 + (fwhm / x0) ** 2)


def _modified_drude(x, amp, x0, fwhm, asym):
    return _drude(x, amp, x0, 2.0 * fwhm / (1.0 + np.exp(asym * (x - x0))))


def _fm90(x, c):
    """Fitzpatrick & Massa (1990) UV curve E(x-V)/E(B-V); c = [c1..c4, x0, gamma]."""
    x = np.asarray(x, dtype=np.float64)
    out = c[0] + c[1] * x + c[2] * x**2 / ((x**2 - c[4]**2) ** 2 + (x * c[5]) ** 2)
    y = np.maximum(x - 5.9, 0.0)
    return out + np.where(x >= 5.9, c[3] * (0.5392 * y**2 + 0.05644 * y**3), 0.0)


def _natural_spline(xs, ys):
    """Second derivatives of the natural cubic spline through (xs, ys)."""
    n = len(xs)
    m = np.zeros(n)
    c = np.zeros(n)
    r = np.zeros(n)
    for i in range(1, n - 1):
        h_lo = xs[i] - xs[i - 1]
        h_hi = xs[i + 1] - xs[i]
        w = 2.0 * (h_lo + h_hi) - h_lo * c[i - 1]
        c[i] = h_hi / w
        r[i] = (6.0 * ((ys[i + 1] - ys[i]) / h_hi - (ys[i] - ys[i - 1]) / h_lo)
                - h_lo * r[i - 1]) / w
    for i in range(n - 2, 0, -1):
        m[i] = r[i] - c[i] * m[i + 1]
    return m


def _spline_value(xs, ys, m, x):
    """Value at x of the cubic spline with knots (xs, ys) and second derivatives m.

    Outside the knots the end cubic is extended."""
    j = np.searchsorted(xs[1:-1], x, side="left")
    h = xs[j + 1] - xs[j]
    t_lo = xs[j + 1] - x
    t_hi = x - xs[j]
    return ((m[j] * t_lo**3 + m[j + 1] * t_hi**3) / (6.0 * h)
            + (ys[j] / h - m[j] * h / 6.0) * t_lo
            + (ys[j + 1] / h - m[j + 1] * h / 6.0) * t_hi)


# ---------------------------------------------------------------------------
# Building blocks shared by the laws.  All take x = 1/lambda [um^-1] and
# return k = A(lambda)/A(V).
# ---------------------------------------------------------------------------

def _cardelli_form(x, r_v, opt_a, opt_b):
    """k = a(x) + b(x)/R_V in the CCM89 form with the given optical polynomials."""
    xi = np.minimum(x, 10.0)
    a = np.empty_like(xi)
    b = np.empty_like(xi)

    m = xi < 1.1                                  # infrared power law
    p = xi[m] ** 1.61
    a[m] = 0.574 * p
    b[m] = -0.527 * p

    m = (xi >= 1.1) & (xi < 3.3)                  # optical
    y = xi[m] - 1.82
    a[m] = _polynomial(opt_a, y)
    b[m] = _polynomial(opt_b, y)

    m = (xi >= 3.3) & (xi <= 8.0)                 # near-UV
    xm = xi[m]
    am = 1.752 - 0.316 * xm - 0.104 / ((xm - 4.67) ** 2 + 0.341)
    bm = -3.090 + 1.825 * xm + 1.206 / ((xm - 4.62) ** 2 + 0.263)
    y = np.maximum(xm - 5.9, 0.0)                 # extra curvature above x = 5.9
    hi = xm >= 5.9
    am = am + np.where(hi, -0.04473 * y**2 - 0.009779 * y**3, 0.0)
    bm = bm + np.where(hi, 0.2130 * y**2 + 0.1207 * y**3, 0.0)
    a[m] = am
    b[m] = bm

    m = xi > 8.0                                  # far-UV
    y = xi[m] - 8.0
    a[m] = -1.073 - 0.628 * y + 0.137 * y**2 - 0.070 * y**3
    b[m] = 13.670 + 4.257 * y - 0.420 * y**2 + 0.374 * y**3

    return a + b / r_v


def _f99_form(x, r_v, fm90_c, knot_x, knot_k):
    """The F99 construction: FM90 shortward of 2700 A, longward a natural cubic
    spline through k = 0 at x = 0, the optical/IR anchors (knot_x, knot_k) and
    the FM90 curve at 2700 and 2600 A.  FM90 gives E(x-V)/E(B-V), so
    k = 1 + FM90/R_V."""
    x_uv = 1.0e4 / 2700.0
    x_uv2 = 1.0e4 / 2600.0
    sx = np.concatenate(([0.0], knot_x, [x_uv, x_uv2]))
    sk = np.concatenate(([0.0], knot_k,
                         [1.0 + _fm90(x_uv, fm90_c) / r_v, 1.0 + _fm90(x_uv2, fm90_c) / r_v]))
    m = _natural_spline(sx, sk)

    xi = np.minimum(x, 10.0)
    k = np.empty_like(xi)
    uv = xi >= x_uv
    k[uv] = 1.0 + _fm90(xi[uv], fm90_c) / r_v
    k[~uv] = _spline_value(sx, sk, m, xi[~uv])
    return k


# ---------------------------------------------------------------------------
# Laws: k(x, ...) for the Milky Way and Magellanic Cloud curves
# ---------------------------------------------------------------------------

_CCM89_A = [1.0, 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530, 0.32999]
_CCM89_B = [0.0, 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260, -2.09002]

_O94_A = [1.0, 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505]
_O94_B = [0.0, 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347]


def _k_ccm89(x, r_v):
    return _cardelli_form(x, r_v, _CCM89_A, _CCM89_B)


def _k_o94(x, r_v):
    return _cardelli_form(x, r_v, _O94_A, _O94_B)


_F99_KNOT_LAMBDA = np.array([26500.0, 12200.0, 6000.0, 5470.0, 4670.0, 4110.0])


def _k_f99(x, r_v):
    # A(lambda)/E(B-V) at the anchors (the last optical anchor uses +1.208,
    # which reproduces F99 table 3; table 4 prints -1.208)
    knot_k = np.empty(6)
    knot_k[0:2] = np.array([0.265, 0.829]) * r_v / 3.1
    knot_k[2] = -0.426 + 1.0044 * r_v
    knot_k[3] = -0.050 + 1.0016 * r_v
    knot_k[4] = 0.701 + 1.0016 * r_v
    knot_k[5] = 1.208 + 1.0032 * r_v - 0.00033 * r_v**2
    c2 = -0.824 + 4.717 / r_v
    return _f99_form(x, r_v, [2.030 - 3.007 * c2, c2, 3.23, 0.41, 4.596, 0.99],
                     1.0e4 / _F99_KNOT_LAMBDA, knot_k / r_v)


def _fm07_uv(xu):
    xu = np.asarray(xu, dtype=np.float64)
    out = (-0.175 + 0.807 * xu
           + 2.991 * xu**2 / ((xu**2 - 4.592**2) ** 2 + (0.922 * xu) ** 2))
    return out + np.where(xu > 6.097, 0.319 * (xu - 6.097) ** 2, 0.0)


_FM07_R_V = 3.1
_FM07_X_UV = 1.0e4 / 2700.0
_FM07_SX = np.array([0.0, 0.25, 0.50, 0.75, 1.0, 1.0e4 / 5530.0, 1.0e4 / 4000.0,
                     1.0e4 / 3300.0, 1.0e4 / 2700.0, 1.0e4 / 2600.0])
_FM07_SK = np.concatenate((
    (0.63 * _FM07_R_V - 0.83) * _FM07_SX[:5] ** 1.84 - _FM07_R_V,
    [0.0, 1.322, 2.055],
    [float(_fm07_uv(_FM07_SX[8])), float(_fm07_uv(_FM07_SX[9]))],
))
_FM07_M = _natural_spline(_FM07_SX, _FM07_SK)


def _k_fm07(x):
    """FM07 average Milky Way curve, fixed at R_V = 3.1; FM90 in the UV with a
    quadratic far-UV rise above x = 6.097, a natural cubic spline in
    E(x-V)/E(B-V) longward of 2700 A."""
    xi = np.minimum(x, _X_FUV)
    k = np.empty_like(xi)
    uv = xi >= _FM07_X_UV
    k[uv] = 1.0 + _fm07_uv(xi[uv]) / _FM07_R_V
    k[~uv] = 1.0 + _spline_value(_FM07_SX, _FM07_SK, _FM07_M, xi[~uv]) / _FM07_R_V
    return k


def _g23_ir(x):
    """G23 near/mid-IR: a is a power law with a smooth break near 4.3 micron
    plus the two silicate features, b a single power law."""
    scale, alpha1, alpha2 = 0.38526, 1.68467, 0.78791
    break_lam, break_width = 4.30578, 4.78338
    lam = 1.0 / x
    w = _smoothstep(lam, break_lam - 0.5 * break_width, break_lam + 0.5 * break_width)
    a = (scale * ((1.0 - w) * lam ** (-alpha1)
                  + w * break_lam ** (alpha2 - alpha1) * lam ** (-alpha2))
         + _modified_drude(lam, 0.06652, 9.8434, 2.21205, -0.24703)
         + _modified_drude(lam, 0.0267, 19.58294, 17.0, -0.27))
    b = -1.01251 * x**1.06099
    return a, b


def _g23_optical(x):
    """G23 optical: quartic in x plus three intermediate-scale (Drude) features."""
    a = (_polynomial([-0.35848, 0.7122, 0.08746, -0.05403, 0.00674], x)
         + _drude(x, 0.03893, 2.288, 0.243)
         + _drude(x, 0.02965, 2.054, 0.179)
         + _drude(x, 0.01747, 1.587, 0.243))
    b = (_polynomial([0.12354, -2.68335, 2.01901, -0.39299, 0.03355], x)
         + _drude(x, 0.18453, 2.288, 0.243)
         + _drude(x, 0.19728, 2.054, 0.179)
         + _drude(x, 0.1713, 1.587, 0.243))
    return a, b


def _g23_uv(x):
    """G23 ultraviolet: FM90 curves for a and b."""
    a = _fm90(x, [0.81297, 0.2775, 1.06295, 0.11303, 4.60, 0.99])
    b = _fm90(x, [-2.97868, 1.89808, 3.10334, 0.65484, 4.60, 0.99])
    return a, b


def _k_g23(x, r_v):
    """G23 (Milky Way, R_V dependent); fitted for 0.0912 <= lambda <= 32 micron.
    Optical blends into the IR over 0.9-1.1 micron and the UV into the optical
    over 0.30-0.33 micron.  The 1/0.9854 renormalisation corrects the F19
    photometry to the HST/STIS spectrophotometry."""
    xi = np.minimum(x, 1.0 / 0.0912)
    lam = 1.0 / xi
    a = np.empty_like(xi)
    b = np.empty_like(xi)

    m = lam > 1.1
    a[m], b[m] = _g23_ir(xi[m])

    m = (lam >= 0.9) & (lam <= 1.1)
    a1, b1 = _g23_optical(xi[m])
    a2, b2 = _g23_ir(xi[m])
    w = _smoothstep(lam[m], 0.9, 1.1)
    a[m] = (1.0 - w) * a1 + w * a2
    b[m] = (1.0 - w) * b1 + w * b2

    m = (lam > 0.33) & (lam < 0.9)
    a[m], b[m] = _g23_optical(xi[m])

    m = (lam >= 0.3) & (lam <= 0.33)
    a1, b1 = _g23_uv(xi[m])
    a2, b2 = _g23_optical(xi[m])
    w = _smoothstep(lam[m], 0.3, 0.33)
    a[m] = (1.0 - w) * a1 + w * a2
    b[m] = (1.0 - w) * b1 + w * b2

    m = lam < 0.3
    a[m], b[m] = _g23_uv(xi[m])

    return (a + b * (1.0 / r_v - 1.0 / 3.1)) / 0.9854


def _k_g03_smcbar(x):
    """Gordon et al. (2003) SMC bar average, R_V = 2.74."""
    return _f99_form(x, 2.74, [-4.959, 2.264, 0.389, 0.461, 4.6, 1.0],
                     1.0 / np.array([2.198, 1.65, 1.25, 0.81, 0.65, 0.55, 0.44, 0.37]),
                     np.array([0.110, 0.169, 0.250, 0.567, 0.801, 1.000, 1.374, 1.672]))


def _k_g03_lmcavg(x):
    """Gordon et al. (2003) LMC average, R_V = 3.41."""
    return _f99_form(x, 3.41, [-0.890, 0.998, 2.719, 0.400, 4.579, 0.934],
                     1.0 / np.array([2.198, 1.65, 1.25, 0.55, 0.44, 0.37]),
                     np.array([0.100, 0.186, 0.257, 1.000, 1.293, 1.518]))


# ---------------------------------------------------------------------------
# Public law functions: A(lambda) in magnitudes
# ---------------------------------------------------------------------------

def _shape(wave, k):
    return k if np.ndim(wave) else k[0]


def ccm89(wave: np.ndarray, a_v: float, r_v: float = 3.1) -> np.ndarray:
    """Cardelli, Clayton & Mathis (1989) extinction in magnitudes.

    Fitted for 1000 A to 3.3 micron (0.3 <= x <= 10 um^-1).  Shortward of
    1000 A k is held at its 1000 A value; longward the infrared power law
    continues.  Requires 2 <= R_V <= 6 (``ValueError`` otherwise).

    Parameters
    ----------
    wave : array_like
        Wavelengths in Angstroms.
    a_v : float
        V-band extinction in magnitudes.
    r_v : float
        Total-to-selective extinction ratio (default 3.1).

    Returns
    -------
    A_lambda : ndarray
        Extinction in magnitudes at each wavelength.
    """
    _check_rv("ccm89", r_v)
    return _shape(wave, a_v * _k_ccm89(_wave_to_x(wave), r_v))


def odonnell94(wave: np.ndarray, a_v: float, r_v: float = 3.1) -> np.ndarray:
    """O'Donnell (1994) extinction: CCM89 with a refitted optical section.

    Same fitted range and out-of-range behaviour as :func:`ccm89`.
    Requires 2 <= R_V <= 6 (``ValueError`` otherwise).

    Parameters
    ----------
    wave : array_like
        Wavelengths in Angstroms.
    a_v : float
        V-band extinction in magnitudes.
    r_v : float
        Total-to-selective extinction ratio (default 3.1).

    Returns
    -------
    A_lambda : ndarray
        Extinction in magnitudes at each wavelength.
    """
    _check_rv("odonnell94", r_v)
    return _shape(wave, a_v * _k_o94(_wave_to_x(wave), r_v))


def fitzpatrick99(wave: np.ndarray, a_v: float, r_v: float = 3.1) -> np.ndarray:
    """Fitzpatrick (1999) R_V-dependent extinction.

    FM90 in the UV (shortward of 2700 A); longward a natural cubic spline
    through the optical/IR anchors used by the ``dust_extinction`` package.
    Fitted for 1000 A to 3.3 micron; shortward of 1000 A k is held at its
    1000 A value, longward the spline continues.  Requires 2 <= R_V <= 6
    (``ValueError`` otherwise).

    Parameters
    ----------
    wave : array_like
        Wavelengths in Angstroms.
    a_v : float
        V-band extinction in magnitudes.
    r_v : float
        Total-to-selective extinction ratio (default 3.1).

    Returns
    -------
    A_lambda : ndarray
        Extinction in magnitudes at each wavelength.
    """
    _check_rv("fitzpatrick99", r_v)
    return _shape(wave, a_v * _k_f99(_wave_to_x(wave), r_v))


def fm07(wave: np.ndarray, a_v: float) -> np.ndarray:
    """Fitzpatrick & Massa (2007) average Milky Way extinction (R_V = 3.1 fixed).

    Fitted from 910 A to 6 micron; shortward of 910 A k is held at its 910 A
    value, longward the spline continues.

    Parameters
    ----------
    wave : array_like
        Wavelengths in Angstroms.
    a_v : float
        V-band extinction in magnitudes.

    Returns
    -------
    A_lambda : ndarray
        Extinction in magnitudes at each wavelength.
    """
    return _shape(wave, a_v * _k_fm07(_wave_to_x(wave)))


def g03_smcbar(wave: np.ndarray, a_v: float) -> np.ndarray:
    """Gordon et al. (2003) SMC bar average extinction (R_V = 2.74 fixed).

    Fitted for 1000 A to 3.3 micron (0.3 <= x <= 10 um^-1); shortward of
    1000 A k is held at its 1000 A value.

    Parameters
    ----------
    wave : array_like
        Wavelengths in Angstroms.
    a_v : float
        V-band extinction in magnitudes.
    """
    return _shape(wave, a_v * _k_g03_smcbar(_wave_to_x(wave)))


def g03_lmcavg(wave: np.ndarray, a_v: float) -> np.ndarray:
    """Gordon et al. (2003) LMC average extinction (R_V = 3.41 fixed).

    Fitted for 1000 A to 3.3 micron (0.3 <= x <= 10 um^-1); shortward of
    1000 A k is held at its 1000 A value.

    Parameters
    ----------
    wave : array_like
        Wavelengths in Angstroms.
    a_v : float
        V-band extinction in magnitudes.
    """
    return _shape(wave, a_v * _k_g03_lmcavg(_wave_to_x(wave)))


# ---------------------------------------------------------------------------
# Calzetti et al. (2000)  — starburst attenuation law
# ---------------------------------------------------------------------------

def calzetti00(wave: np.ndarray, a_v: float, r_v: float = 4.05) -> np.ndarray:
    """Calzetti et al. (2000) starburst galaxy attenuation law.

    Valid range: 1200–22 000 Å.

    Parameters
    ----------
    wave : array_like
        Wavelengths in Angstroms.
    a_v : float
        V-band attenuation in magnitudes (a_v = r_v * E(B−V)_s).
    r_v : float
        R_V for the starburst law (Calzetti default: 4.05).

    Returns
    -------
    A_lambda : ndarray
        Attenuation in magnitudes at each wavelength.
    """
    wave  = np.asarray(wave, dtype=np.float64)
    w_um  = wave * 1e-4

    k = np.zeros_like(w_um)

    m = (w_um >= 0.12) & (w_um < 0.63)
    if m.any():
        w = w_um[m]
        k[m] = (2.659*(-2.156 + 1.509/w - 0.198/w**2 + 0.011/w**3) + r_v)

    m = (w_um >= 0.63) & (w_um <= 2.2)
    if m.any():
        w = w_um[m]
        k[m] = 2.659*(-1.857 + 1.040/w) + r_v

    ebv_s = a_v / r_v
    return ebv_s * k


# ---------------------------------------------------------------------------
# Gordon et al. (2023)
# ---------------------------------------------------------------------------

# Parameters of the earlier-version 'lmc' and 'smc' environments.  They are an
# FM90-style curve in the UV with a c6 * x**1.84 - R_V power law longward of
# 3 micron; they are not taken from G23 (which publishes only the Milky Way
# R_V-dependent relation) and are kept unchanged for compatibility.  They
# could not be checked against any reference implementation.
_G23_LEGACY_PARAMS = {
    'lmc': dict(
        c1=-1.475, c2=1.132, c3=1.463, c4=0.294, c5=5.9,
        x0=4.558, gamma=0.945, r_v=3.41,
        c6=0.389, alpha=1.84,
    ),
    'smc': dict(
        c1=-4.959, c2=2.264, c3=0.389, c4=0.461, c5=5.9,
        x0=4.600, gamma=1.000, r_v=2.74,
        c6=1.057, alpha=1.84,
    ),
}

_G23_ENVIRONMENTS = ('mw', 'lmc', 'smc')


def _gordon23_legacy(wave, a_v, rv, env):
    """Pre-0.5 'lmc'/'smc' curves (see ``_G23_LEGACY_PARAMS``); returns A in mag."""
    p = _G23_LEGACY_PARAMS[env]
    c1, c2, c3, c4, c5 = p['c1'], p['c2'], p['c3'], p['c4'], p['c5']
    x0, gamma = p['x0'], p['gamma']
    c6, alpha = p['c6'], p['alpha']
    if rv is None:
        rv = p['r_v']

    x = 1e4 / np.asarray(wave, dtype=np.float64)
    k = np.zeros_like(x)

    ir_mask = x < (1.0 / 3.0)
    if ir_mask.any():
        k[ir_mask] = c6 * x[ir_mask]**alpha - rv

    opt_uv_mask = ~ir_mask
    if opt_uv_mask.any():
        xm = x[opt_uv_mask]
        xm2 = xm**2
        d = xm2 / ((xm2 - x0**2)**2 + xm2 * gamma**2)
        k_uv = c1 + c2*xm + c3*d
        fuv = xm > c5
        if fuv.any():
            y = xm[fuv] - c5
            k_uv[fuv] += c4*(0.5392*y**2 + 0.05644*y**3)
        k[opt_uv_mask] = k_uv

    return a_v / rv * (k + rv)


def gordon23(wave: np.ndarray, a_v: float, r_v: Optional[float] = None,
             environment: str = 'mw') -> np.ndarray:
    """Gordon et al. (2023) dust extinction law.

    ``environment='mw'`` (default) is the Milky Way R_V-dependent relation of
    Gordon, K. D. et al. 2023, ApJ, 950, 86, ported from the MESA ``colors``
    module.  It is fitted for 912 A to 32 micron, including the 10 and 20
    micron silicate features; shortward of 912 A k is held at its 912 A value,
    longward of 32 micron the power laws continue.  Requires
    2.3 <= R_V <= 5.6 (``ValueError`` otherwise).

    ``environment='lmc'`` and ``'smc'`` select the older curves this package
    used in earlier versions; they are not part of G23, not in MESA, and could not be
    verified against any reference, so they emit a ``UserWarning`` and R_V is
    not range-checked.  For the Magellanic Clouds use :func:`g03_lmcavg` and
    :func:`g03_smcbar`, which follow Gordon et al. (2003) and are validated.

    Parameters
    ----------
    wave : array_like
        Wavelengths in Angstroms.
    a_v : float
        V-band extinction in magnitudes.
    r_v : float or None
        Total-to-selective extinction ratio.  None means 3.1 for 'mw' and the
        environment default (LMC 3.41, SMC 2.74) for the older curves.
    environment : {'mw', 'lmc', 'smc'}
        Which curve to use (default 'mw').

    Returns
    -------
    A_lambda : ndarray
        Extinction in magnitudes at each wavelength.
    """
    env = environment.lower()
    if env not in _G23_ENVIRONMENTS:
        raise ValueError(
            f"environment must be one of {list(_G23_ENVIRONMENTS)}; got '{env}'")

    if env == 'mw':
        rv = 3.1 if r_v is None else float(r_v)
        _check_rv("gordon23", rv)
        return _shape(wave, a_v * _k_g23(_wave_to_x(wave), rv))

    warnings.warn(
        f"gordon23 environment '{env}' is the earlier-version curve, which is not part of "
        "Gordon et al. (2023) and could not be verified against a reference; use "
        f"law='{'g03_lmcavg' if env == 'lmc' else 'g03_smcbar'}' for the Gordon "
        "et al. (2003) curve.",
        UserWarning, stacklevel=2,
    )
    wave = np.asarray(wave, dtype=np.float64)
    return _gordon23_legacy(wave, a_v, r_v, env)


# ---------------------------------------------------------------------------
# apply / remove helpers
# ---------------------------------------------------------------------------

def apply_extinction(a_lambda: np.ndarray, flux: np.ndarray) -> np.ndarray:
    """Apply extinction magnitudes to a flux array.

    F_obs = F_intrinsic × 10^(−0.4 × A_lambda)

    Parameters
    ----------
    a_lambda : ndarray
        Extinction in magnitudes at each wavelength (same shape as flux).
    flux : ndarray
        Intrinsic flux array.

    Returns
    -------
    flux_obs : ndarray
        Observed (extincted) flux.
    """
    return flux * 10.0**(-0.4 * a_lambda)


def remove_extinction(a_lambda: np.ndarray, flux: np.ndarray) -> np.ndarray:
    """Remove extinction from an observed flux (de-redden).

    F_intrinsic = F_obs × 10^(+0.4 × A_lambda)

    Parameters
    ----------
    a_lambda : ndarray
        Extinction in magnitudes at each wavelength.
    flux : ndarray
        Observed (extincted) flux.

    Returns
    -------
    flux_intrinsic : ndarray
        De-reddened flux.
    """
    return flux * 10.0**(0.4 * a_lambda)


# ---------------------------------------------------------------------------
# High-level config dataclass
# ---------------------------------------------------------------------------

@dataclass
class ExtinctionConfig:
    """Configuration container for extinction settings.

    Parameters
    ----------
    enabled : bool
        If False, ``ExtinctionModel.apply()`` is a no-op.  Default False.
    law : str
        Name of the extinction law.  One of ``AVAILABLE_LAWS``.
    r_v : float
        R_V value.  Ignored by the fixed-shape laws 'fm07' (3.1),
        'g03_smcbar' (2.74) and 'g03_lmcavg' (3.41).  For 'ccm89',
        'odonnell94', 'fitzpatrick99' (2 to 6) and 'gordon23' with
        environment 'mw' (2.3 to 5.6) a value outside the calibrated range
        raises ``ValueError``.  'calzetti00' is not range-checked.
    a_v : float
        V-band extinction in magnitudes.  Default 0.0.
    gordon23_env : str
        Only used when ``law='gordon23'``.  One of 'mw', 'lmc', 'smc'
        ('lmc' and 'smc' are unverified older curves, see :func:`gordon23`).
    """
    enabled: bool      = False
    law: str           = 'fitzpatrick99'
    r_v: float         = 3.1
    a_v: float         = 0.0
    gordon23_env: str  = 'mw'

    def __post_init__(self) -> None:
        if self.law not in AVAILABLE_LAWS:
            raise ValueError(
                f"law '{self.law}' not recognised.  "
                f"Choose from {AVAILABLE_LAWS}."
            )
        if self.a_v < 0.0:
            raise ValueError("a_v must be >= 0")
        if self.law == 'gordon23' and str(self.gordon23_env).lower() != 'mw':
            return
        _check_rv(self.law, self.r_v)


# ---------------------------------------------------------------------------
# Law dispatch table — unified signature (wave, av, rv, env)
# ---------------------------------------------------------------------------

_LAW_DISPATCH = {
    'ccm89':         lambda w, av, rv, env: ccm89(w, av, rv),
    'odonnell94':    lambda w, av, rv, env: odonnell94(w, av, rv),
    'fitzpatrick99': lambda w, av, rv, env: fitzpatrick99(w, av, rv),
    'fm07':          lambda w, av, rv, env: fm07(w, av),
    'calzetti00':    lambda w, av, rv, env: calzetti00(w, av, rv),
    'g03_smcbar':    lambda w, av, rv, env: g03_smcbar(w, av),
    'g03_lmcavg':    lambda w, av, rv, env: g03_lmcavg(w, av),
    # For the older 'lmc'/'smc' curves, pass rv=None when r_v was left at the
    # generic default (3.1) so they use their own R_V, as before.
    'gordon23':      lambda w, av, rv, env: gordon23(
        w, av, rv if str(env).lower() == 'mw' or rv != 3.1 else None, env),
}


# ---------------------------------------------------------------------------
# ExtinctionModel — the main user-facing class
# ---------------------------------------------------------------------------

class ExtinctionModel:
    """Apply extinction to a spectrum.

    Parameters
    ----------
    config : ExtinctionConfig, optional
        Full configuration object.  If None, a default (disabled) config is used.
    **kwargs
        Convenience keyword arguments forwarded to ``ExtinctionConfig``.
        These override fields in *config* if both are supplied.

    Examples
    --------
    Disabled (default) — passes flux through unchanged::

        ext = ExtinctionModel()
        flux_out = ext.apply(wave, flux)

    Enabled, Fitzpatrick99, Av=0.5::

        ext = ExtinctionModel(enabled=True, law='fitzpatrick99', a_v=0.5, r_v=3.1)
        flux_obs = ext.apply(wave, flux)

    Gordon+2023 SMC curve::

        ext = ExtinctionModel(enabled=True, law='gordon23', a_v=0.8, gordon23_env='smc')
        flux_obs = ext.apply(wave, flux)
    """

    def __init__(self,
                 config: Optional[ExtinctionConfig] = None,
                 **kwargs) -> None:
        if config is None:
            config = ExtinctionConfig(**kwargs)
        elif kwargs:
            import dataclasses
            d = dataclasses.asdict(config)
            d.update(kwargs)
            config = ExtinctionConfig(**d)
        self.config = config

    # ------------------------------------------------------------------
    # Core methods
    # ------------------------------------------------------------------

    def extinction_curve(self, wave: np.ndarray) -> np.ndarray:
        """Return A(λ) in magnitudes for the configured law and A_V.

        Parameters
        ----------
        wave : array_like
            Wavelengths in Angstroms.

        Returns
        -------
        A_lambda : ndarray
            Extinction in magnitudes.  Returns zeros if ``enabled`` is False.
        """
        wave = np.asarray(wave, dtype=np.float64)
        if not self.config.enabled or self.config.a_v == 0.0:
            return np.zeros(len(wave))
        fn = _LAW_DISPATCH[self.config.law]
        return fn(wave, self.config.a_v, self.config.r_v, self.config.gordon23_env)

    def apply(self, wave: np.ndarray, flux: np.ndarray) -> np.ndarray:
        """Apply extinction to a flux array.

        Parameters
        ----------
        wave : array_like
            Wavelengths in Angstroms.
        flux : array_like
            Intrinsic flux in any linear unit.

        Returns
        -------
        flux_out : ndarray
            Extincted flux.  If ``enabled`` is False this is identical to
            the input flux (no copy is made).
        """
        wave = np.asarray(wave, dtype=np.float64)
        flux = np.asarray(flux, dtype=np.float64)
        if not self.config.enabled:
            return flux
        a_lam = self.extinction_curve(wave)
        return apply_extinction(a_lam, flux)

    def remove(self, wave: np.ndarray, flux: np.ndarray) -> np.ndarray:
        """Remove (de-redden) extinction from an observed flux.

        Parameters
        ----------
        wave : array_like
            Wavelengths in Angstroms.
        flux : array_like
            Observed flux.

        Returns
        -------
        flux_intrinsic : ndarray
            De-reddened flux.
        """
        wave = np.asarray(wave, dtype=np.float64)
        flux = np.asarray(flux, dtype=np.float64)
        if not self.config.enabled:
            return flux
        a_lam = self.extinction_curve(wave)
        return remove_extinction(a_lam, flux)

    def with_av(self, a_v: float) -> "ExtinctionModel":
        """Return a new ExtinctionModel identical to this one but with a different A_V."""
        import dataclasses
        return ExtinctionModel(dataclasses.replace(self.config, a_v=a_v))

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        c = self.config
        if not c.enabled:
            return "ExtinctionModel(enabled=False)"
        return (f"ExtinctionModel(law='{c.law}', a_v={c.a_v}, r_v={c.r_v})")

    @classmethod
    def disabled(cls) -> "ExtinctionModel":
        """Return an extinction model that does nothing (convenience factory)."""
        return cls(ExtinctionConfig(enabled=False))

    @classmethod
    def from_dict(cls, d: dict) -> "ExtinctionModel":
        """Construct from a plain dict (e.g. loaded from a config file)."""
        return cls(ExtinctionConfig(**d))


# ---------------------------------------------------------------------------
# Integration helper for forward models
# ---------------------------------------------------------------------------

def make_extinction_model(
    enabled: bool = False,
    law: str = 'fitzpatrick99',
    r_v: float = 3.1,
    a_v: float = 0.0,
    gordon23_env: str = 'mw',
) -> ExtinctionModel:
    """Convenience factory — mirrors the ExtinctionConfig keyword signature.

    This is the recommended way to build an extinction model inside a fitter's
    forward-model function, e.g.::

        ext = make_extinction_model(enabled=True, law='fitzpatrick99',
                                    a_v=theta['a_v'])
        flux_obs = ext.apply(wave, flux_intrinsic)

    Parameters
    ----------
    enabled : bool
        Master switch.  Default False — extinction is skipped entirely.
    law : str
        Extinction law name.  See ``AVAILABLE_LAWS``.
    r_v : float
        R_V (shape parameter).  Default 3.1.
    a_v : float
        V-band extinction.  Default 0.0.
    gordon23_env : str
        Gordon+2023 environment preset: 'mw', 'lmc', or 'smc'.

    Returns
    -------
    ExtinctionModel
    """
    return ExtinctionModel(ExtinctionConfig(
        enabled=enabled,
        law=law,
        r_v=r_v,
        a_v=a_v,
        gordon23_env=gordon23_env,
    ))
