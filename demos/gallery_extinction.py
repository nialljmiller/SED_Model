"""
gallery_extinction.py
=====================
Interstellar extinction showcase for the gallery.

Writes to ``docs/assets/gallery/``:

* ``dust.mp4`` / ``dust.gif``   dust piling up in front of a 7000 K star
* ``extinction_laws.png``       A(lambda)/A(V) for the supported laws
* ``reddening_vectors.png``     colour-magnitude diagram with A_V = 1 vectors

Usage::

    python demos/gallery_extinction.py            # full quality
    python demos/gallery_extinction.py --quick    # fast test, written to a temp dir

Needs the SED_Tools Kurucz2003all grid, the Johnson filters and the Vega SED
(see ``viz_common.data_paths``) and ffmpeg for the movie.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import viz_common as vc  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

import sed_model as sm  # noqa: E402

RSUN = sm.RSUN_TO_CM
PC = sm.PC_TO_CM
BANDS = "UBVRI"
D_CM = 200.0 * PC


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def spectrum_color(wave_A, flux):
    """Apparent sRGB chromaticity of a spectrum (same pipeline as star_color)."""
    nm = np.linspace(380.0, 780.0, 401)
    f = np.interp(nm * 10.0, wave_A, flux)
    x, y, z = vc._cie_xyz(nm)
    xyz = np.array([np.trapezoid(f * x, nm), np.trapezoid(f * y, nm), np.trapezoid(f * z, nm)])
    xyz = xyz / xyz.sum()
    X, Y, Z = xyz[0] / xyz[1], 1.0, xyz[2] / xyz[1]
    rgb = np.array([3.2406 * X - 1.5372 * Y - 0.4986 * Z,
                    -0.9689 * X + 1.8758 * Y + 0.0415 * Z,
                    0.0557 * X - 0.2040 * Y + 1.0570 * Z])
    rgb = np.clip(rgb, 0.0, None)
    rgb = rgb / rgb.max()
    rgb = np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb ** (1 / 2.4) - 0.055)
    return np.clip(rgb, 0.0, 1.0)


def load_all(paths):
    grid = sm.load_grid(paths["grid"])
    jf = [f"{paths['johnson']}/{b}.dat" for b in BANDS]
    filters = sm.load_filters(jf, paths["vega"])
    return grid, filters


def forward(grid, filters, teff, logg, r_rsun, ext=None):
    if grid.cell_has_empty_nodes(teff, logg, 0.0, stencil=True):
        raise ValueError(f"({teff}, {logg}) is next to empty grid nodes")
    with warnings.catch_warnings():
        warnings.simplefilter("error")           # any grid warning is a bug here
        return sm.run_forward(teff=teff, logg=logg, meta=0.0, R=r_rsun * RSUN, d=D_CM,
                              grid=grid, filters=filters, mag_system="Vega", extinction=ext)


def ext_model(law, a_v, r_v=3.1):
    return sm.make_extinction_model(enabled=a_v > 0, law=law, r_v=r_v, a_v=a_v)


def mix(c, bg, t):
    """Fade colour c toward the background (t=1 keeps it)."""
    c, bg = np.asarray(c, float), np.asarray(matplotlib_rgb(bg), float)
    return tuple(bg + (c - bg) * t)


def matplotlib_rgb(c):
    import matplotlib.colors as mc
    return mc.to_rgb(c)


def ms_radius(teff):
    return (teff / 5778.0) ** 1.6      # crude main-sequence radius, solar units


# ---------------------------------------------------------------------------
# 1. the movie
# ---------------------------------------------------------------------------

def make_movie(grid, filters, quick):
    teff, logg, rad = 7000.0, 4.3, 1.6
    law, r_v = "fitzpatrick99", 3.1
    n = 36 if quick else 330
    a_peak = 3.0
    u = 1 - abs(2 * np.arange(n) / n - 1)                       # triangle 0 -> 1 -> 0, loops
    a_v = a_peak * (0.65 * u + 0.35 * u * u * (3 - 2 * u))      # short dwell at both ends

    base = forward(grid, filters, teff, logg, rad)
    wave = base.wavelengths
    sel = (wave > 2800) & (wave < 10500)
    w_um = wave[sel] / 1e4
    ymax = 1.12 * base.observed_flux[sel].max()
    ghost = base.observed_flux[sel] / ymax

    # filter strips: edges where transmission falls to half maximum
    edges = []
    for f in filters:
        t = np.asarray(f.transmission)
        w = np.asarray(f.wavelengths)
        hi = w[t >= 0.5 * t.max()]
        edges.append((hi.min() / 1e4, hi.max() / 1e4))

    frames = []
    for a in a_v:
        r = forward(grid, filters, teff, logg, rad, ext_model(law, a, r_v))
        frames.append((r.observed_flux[sel] / ymax, spectrum_color(wave, r.observed_flux),
                       r.magnitudes))
    base_mag = base.magnitudes
    band_dim = np.array([[10 ** (-0.4 * (fr[2][b] - base_mag[b])) for b in BANDS] for fr in frames])
    track = np.array([[fr[2]["B"] - fr[2]["V"], fr[2]["V"]] for fr in frames])

    # backdrop of unreddened main-sequence stars for the companion panel
    bg = []
    for t in np.arange(4500, 9600, 250):
        try:
            r = forward(grid, filters, float(t), 4.4, ms_radius(t))
            bg.append((r.magnitudes["B"] - r.magnitudes["V"], r.magnitudes["V"]))
        except ValueError:
            pass
    bg = np.array(bg)

    # --- figure ----------------------------------------------------------
    vc.apply_style()
    fig = vc.new_figure(1280, 720)
    ax = fig.add_axes([0.07, 0.20, 0.57, 0.72])
    axs = fig.add_axes([0.07, 0.10, 0.57, 0.055])
    axs.set_xlim(0.28, 1.05)
    axc = fig.add_axes([0.73, 0.20, 0.24, 0.72])

    ax.set_xlim(0.28, 1.05)
    ax.set_yscale('log')
    ax.set_ylim(4e-4, 1.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(["0.4", "0.6", "0.8", "1.0 µm"])
    ax.plot(w_um, ghost, color=vc.INK_3, lw=1.1, zorder=2)
    ax.fill_between(w_um, 1e-5, ghost, color=vc.INK_3, alpha=0.12, lw=0, zorder=1)
    ax.set_ylabel("flux", labelpad=4)
    counter = ax.text(0.02, 0.95, "", transform=ax.transAxes, color=vc.INK_2, fontsize=17,
                      va="top", ha="left")
    disc = ax.scatter([0.93], [0.88], s=520, zorder=5, lw=0)
    halo = ax.scatter([0.93], [0.88], s=1800, zorder=4, lw=0, alpha=0.12)
    (glow,) = ax.plot(w_um, ghost, lw=8, alpha=0.18, zorder=3, solid_capstyle='round')
    (line,) = ax.plot(w_um, ghost, lw=2.8, zorder=4)
    state = {"fill": None}

    axs.set_ylim(0, 1)
    axs.set_yticks([])
    for sp in ("left",):
        axs.spines[sp].set_visible(False)
    axs.spines["bottom"].set_visible(False)
    axs.tick_params(axis="x", length=0, labelsize=12)
    centers = [0.5 * (a + b) for a, b in edges]
    axs.set_xticks(centers)
    axs.set_xticklabels(list(BANDS))
    axs.set_xlabel("wavelength (µm)", labelpad=0, color=vc.INK_3, fontsize=10)
    axs.xaxis.set_label_position("bottom")
    strips = []
    for (lo, hi_) in edges:
        axs.add_patch(Rectangle((lo, 0.1), hi_ - lo, 0.9, fill=False, ec=vc.INK_3, lw=0.6, zorder=1))
        p = Rectangle((lo, 0.1), hi_ - lo, 0.9, lw=2.5, ec=vc.SURFACE, zorder=2)
        axs.add_patch(p)
        strips.append(p)
    ax.set_xlabel("")
    axs.set_xlabel("")

    # companion: colour-magnitude diagram
    axc.set_xlim(-0.15, 1.45)
    v0 = track[0, 1]
    axc.set_ylim(max(track[:, 1].max(), bg[:, 1].max()) + 0.4, min(track[:, 1].min(), bg[:, 1].min()) - 0.4)
    axc.scatter(bg[:, 0], bg[:, 1], s=16, color=vc.INK_3, lw=0, zorder=1)
    axc.plot(track[:, 0], track[:, 1], color=vc.INK_3, lw=0.9, ls=(0, (2, 3)), zorder=2)
    for target in (1.0, 2.0, 3.0):
        j = int(np.argmin(abs(a_v - target)))
        axc.plot(*track[j], marker="o", ms=3, color=vc.INK_2, zorder=3)
        axc.text(track[j, 0] + 0.05, track[j, 1], f"{target:.0f}", color=vc.INK_2, fontsize=10,
                 va="center")
    axc.set_xlabel("B − V", labelpad=4)
    axc.set_ylabel("V (mag)", labelpad=4)
    dot = axc.scatter([track[0, 0]], [track[0, 1]], s=90, zorder=5, lw=0)
    trail, = axc.plot([], [], lw=1.6, zorder=4)

    def update(i):
        flux, col, _ = frames[i]
        dim = band_dim[i]
        a = a_v[i]
        line.set_color(col)
        line.set_ydata(flux)
        glow.set_ydata(flux)
        glow.set_color(col)
        if state["fill"] is not None:
            state["fill"].remove()
        state["fill"] = ax.fill_between(w_um, 1e-5, flux, color=col, alpha=0.28, lw=0, zorder=3)
        for p, d in zip(strips, dim):
            p.set_facecolor(mix(col, vc.SURFACE, 0.12 + 0.88 * d ** 0.6))
        counter.set_text(f"A$_V$  {a:.1f}")
        vis = 0.2 + 0.8 * 10 ** (-0.4 * a * 0.6)
        disc.set_facecolor([mix(col, vc.SURFACE, vis)])
        halo.set_facecolor([mix(col, vc.SURFACE, vis * 0.5)])
        dot.set_offsets([track[i]])
        dot.set_facecolor([col])
        k = i + 1 if a_v[i] >= a_v[i - 1] or i == 0 else i + 1
        trail.set_data(track[:k, 0], track[:k, 1])
        trail.set_color(col)

    update(0)
    out = vc.save_movie(fig, update, n, "dust")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# 2. A(lambda)/A(V) for every law
# ---------------------------------------------------------------------------

def make_laws():
    vc.apply_style()
    wave_um = np.logspace(np.log10(0.1), np.log10(3.3), 600)
    wave = wave_um * 1e4

    def k(law, r_v=3.1, lo=0.1, hi=3.3):
        out = np.asarray(ext_model(law, 1.0, r_v).extinction_curve(wave), float)
        out = np.where((wave_um >= lo) & (wave_um <= hi), out, np.nan)
        return out

    fig = plt.figure(figsize=(18, 5.8), dpi=100)
    axes = fig.subplots(1, 3, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.045, right=0.985, top=0.88, bottom=0.14, wspace=0.05)
    G = "#6a6962"                                     # neutral grey for context lines

    def draw(ax, spec, sub):
        for (y, color, label, dx, dy, lw) in spec:
            ax.plot(wave_um, y, color=color, lw=lw, solid_capstyle="round", zorder=3 if color != G else 2)
            ok = np.where(np.isfinite(y))[0]
            j = ok[int(len(ok) * 0.0)] if label.startswith("@L") else None
        ax.set_xscale("log")
        ax.set_xlim(0.1, 3.9)
        ax.set_ylim(-0.1, 6.4)
        ax.set_xticks([0.1, 0.2, 0.3, 0.5, 1, 2, 3])
        ax.set_xticklabels(["0.1", "0.2", "0.3", "0.5", "1", "2", "3"])
        ax.minorticks_off()
        ax.axvline(0.2175, color=vc.INK_3, lw=0.8, ls=(0, (2, 3)), zorder=1)
        ax.text(0.5, 1.03, sub, transform=ax.transAxes, color=vc.INK_2, fontsize=13, ha="center")
        ax.grid(axis="y", color=vc.INK_3, lw=0.4, alpha=0.5)
        ax.set_axisbelow(True)

    def note(ax, text, color, at, wl, y_fn, leader=True):
        """Direct label at ``at`` with a hairline leader to the curve at wavelength ``wl``."""
        j = int(np.argmin(abs(wave_um - wl)))
        ax.annotate(text, xy=(wave_um[j], y_fn[j]), xytext=at, color=color, fontsize=12, va="center",
                    ha="left", zorder=6,
                    arrowprops=dict(arrowstyle="-", color=color, lw=0.7, shrinkA=3, shrinkB=2,
                                    alpha=0.8) if leader else None)

    B, O, A = vc.BLUE, vc.ORANGE, vc.AQUA

    # panel 1: R_V dependence in the Milky Way (Fitzpatrick 1999)
    ax = axes[0]
    ys = {rv: k("fitzpatrick99", rv) for rv in (2.5, 3.1, 5.0)}
    draw(ax, [(ys[2.5], B, "", 0, 0, 2.2), (ys[3.1], O, "", 0, 0, 2.2), (ys[5.0], A, "", 0, 0, 2.2)],
         "R$_V$ dependence")
    note(ax, "R$_V$ 2.5", B, (0.245, 3.95), 0.19, ys[2.5], leader=False)
    note(ax, "R$_V$ 3.1", O, (0.245, 3.2), 0.19, ys[3.1], leader=False)
    note(ax, "R$_V$ 5.0", A, (0.115, 1.0), 0.16, ys[5.0])
    ax.text(0.2175 * 1.04, 5.75, "2175 Å bump", color=vc.INK_2, fontsize=11, va="center")
    ax.text(3.6, 0.35, "IR tail", color=vc.INK_2, fontsize=11, ha="right", va="bottom")
    ax.set_ylabel("A(λ) / A(V)")

    # panel 2: parametrisations at R_V = 3.1
    ax = axes[1]
    y_f99, y_ccm, y_g23 = k("fitzpatrick99"), k("ccm89"), k("gordon23")
    y_o94, y_fm = k("odonnell94"), k("fm07", lo=0.0911, hi=6.0)
    draw(ax, [(y_o94, G, "", 0, 0, 1.4), (y_fm, G, "", 0, 0, 1.4),
              (y_f99, B, "", 0, 0, 2.2), (y_ccm, O, "", 0, 0, 2.2), (y_g23, A, "", 0, 0, 2.2)],
         "Milky Way laws at R$_V$ 3.1")
    note(ax, "Gordon 23", A, (0.30, 4.6), 0.2175, y_g23)
    note(ax, "CCM 89", O, (0.30, 4.1), 0.2175, y_ccm)
    note(ax, "Fitzpatrick 99", B, (0.30, 3.6), 0.2175, y_f99)
    ax.text(0.86, 5.05, "O'Donnell 94,\nFM07 (grey)", color=G, fontsize=11, ha="left", va="center")

    # panel 3: other galaxies
    ax = axes[2]
    y_smc, y_lmc = k("g03_smcbar", 2.74), k("g03_lmcavg", 3.41)
    y_cal = k("calzetti00", 4.05, lo=0.12, hi=2.2)
    draw(ax, [(y_f99, G, "", 0, 0, 1.4), (y_smc, B, "", 0, 0, 2.2), (y_lmc, O, "", 0, 0, 2.2),
              (y_cal, A, "", 0, 0, 2.2)],
         "Magellanic and starburst")
    note(ax, "SMC bar", B, (0.17, 5.4), 0.145, y_smc, leader=False)
    note(ax, "LMC average", O, (0.27, 3.85), 0.2175, y_lmc)
    note(ax, "Calzetti 00", A, (0.6, 2.3), 1.0, y_cal)
    ax.text(0.86, 5.05, "Milky Way F99\n(grey, for scale)", color=G, fontsize=11, ha="left", va="center")

    for a in axes:
        a.set_xlabel("wavelength (µm)")
    path = vc.save_still(fig, "extinction_laws")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# 3. reddening vectors on a CMD
# ---------------------------------------------------------------------------

def make_vectors(grid, filters, quick):
    """Two-colour diagram: a main-sequence sweep and A_V = 1 vectors for three laws."""
    vc.apply_style()
    fig = vc.new_figure(1000, 800)
    ax = fig.add_axes([0.11, 0.11, 0.85, 0.85])

    def colours(t, ext=None):
        m = forward(grid, filters, float(t), 4.4, ms_radius(t), ext).magnitudes
        return m["B"] - m["V"], m["U"] - m["B"]

    pts = []
    for t in np.arange(5000, 12001, 500 if quick else 250):
        try:
            pts.append((t, *colours(t)))
        except ValueError:
            pass
    pts = np.array(pts)
    ax.plot(pts[:, 1], pts[:, 2], color=vc.INK_3, lw=0.8, zorder=1)
    for t, x, y in pts:
        ax.scatter([x], [y], s=70, color=vc.star_color(float(t)), lw=0, zorder=3)

    t0 = 9000.0
    x0, y0 = colours(t0)
    ax.scatter([x0], [y0], s=250, facecolor="none", edgecolor=vc.INK_2, lw=0.8, zorder=4)

    laws = [("fitzpatrick99", 3.1, vc.BLUE, "Milky Way\n(Fitzpatrick 99)", "left", (0.03, 0.0)),
            ("g03_smcbar", 2.74, vc.ORANGE, "SMC bar", "left", (0.03, 0.04)),
            ("calzetti00", 4.05, vc.AQUA, "Calzetti 00\n(starburst)", "right", (-0.02, 0.13))]
    for law, r_v, col, name, ha, (dx, dy) in laws:
        x1, y1 = colours(t0, ext_model(law, 1.0, r_v))
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=6,
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=2.2, shrinkA=9, shrinkB=0,
                                    mutation_scale=17))
        ax.text(x1 + dx, y1 + dy, name, color=col, fontsize=12.5, va="center", ha=ha, linespacing=1.2)
        print(f"  {law:14s} dBV={x1 - x0:.3f}  dUB={y1 - y0:.3f}")

    ax.text(x0 - 0.04, y0, "A$_V$ = 1", color=vc.INK_2, fontsize=12, ha="right", va="center")
    ax.set_xlabel("B − V")
    ax.set_ylabel("U − B")
    ax.set_xlim(-0.15, 1.0)
    ax.set_ylim(0.85, -0.45)
    path = vc.save_still(fig, "reddening_vectors")
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="short movie, sparse star sweep")
    args = ap.parse_args()
    if args.quick:
        print("quick run, written to", vc.set_output_dir(tempfile.mkdtemp(prefix="gallery_quick_")))
    paths = vc.require_data()
    grid, filters = load_all(paths)
    print(make_laws())
    print(make_vectors(grid, filters, args.quick))
    print(make_movie(grid, filters, args.quick))


if __name__ == "__main__":
    main()
