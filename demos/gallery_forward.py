"""
demos/gallery_forward.py
========================
Forward-model showcase for the gallery (uses the released ``sed_model``).

Writes to ``docs/assets/gallery/``:

* ``star_heating.mp4`` / ``.gif``  a star heated from 3500 K to 30000 K at fixed
  radius: its SED, the Johnson bands it falls through, and its colour-magnitude track.
* ``cmd_grid.png``      constant-gravity lines (dwarfs to giants) in the B-V / M_V plane.
* ``between_nodes.png`` Hermite interpolation between two neighbouring grid nodes.

Every parameter used lies in a grid cell whose nodes (and Hermite stencil) are all
populated, so the forward model never emits the empty-cell warning; warnings are
turned into errors so a regression would be loud rather than silent.

    python demos/gallery_forward.py            # full render
    python demos/gallery_forward.py --quick    # fewer frames, written to a temp dir
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import viz_common as vc  # noqa: E402

RSUN = 6.957e10          # cm
D10PC = 3.0857e19        # cm
MSUN_G = 1.989e33
G_CGS = 6.674e-8

BANDS = ["U", "B", "V", "R", "I", "J"]
T_LO, T_HI = 3500.0, 30000.0


def load_all():
    from sed_model import load_grid, load_filters_from_instrument_dir
    p = vc.require_data()
    grid = load_grid(p["grid"])
    filters = load_filters_from_instrument_dir(p["johnson"], vega_sed_path=p["vega"])
    byname = {f.name: f for f in filters}
    return grid, [byname[b] for b in BANDS]


def safe(grid, teff, logg, meta=0.0):
    return not grid.cell_has_empty_nodes(teff, logg, meta, stencil=True)


def fwd(grid, filters, teff, logg=4.5, meta=0.0, R=RSUN, method="hermite_bounded"):
    from sed_model import run_forward
    assert safe(grid, teff, logg, meta), (teff, logg, meta)
    return run_forward(teff=teff, logg=logg, meta=meta, R=R, d=D10PC, grid=grid,
                       filters=filters, mag_system="Vega", interp_method=method)


def radius_for_logg(logg, mass=1.0):
    return float(np.sqrt(G_CGS * mass * MSUN_G / 10.0 ** logg))


def lamflam(r):
    w = r.wavelengths
    return w / 10.0, w * r.observed_flux          # nm, lambda*F_lambda


def glow(ax, x, y, color, lw=1.8, zorder=5):
    """A bright core with soft halos; returns the artists (widest first)."""
    arts = []
    for w, a in ((11, 0.05), (6, 0.10), (3, 0.20)):
        arts.append(ax.plot(x, y, color=color, lw=w, alpha=a, zorder=zorder - 1,
                            solid_capstyle="round")[0])
    arts.append(ax.plot(x, y, color=color, lw=lw, alpha=1.0, zorder=zorder,
                        solid_capstyle="round")[0])
    return arts


def lighten(c, f=0.35):
    return tuple(np.array(c) * (1 - f) + f)


# ---------------------------------------------------------------------------
# 1. movie
# ---------------------------------------------------------------------------

def make_movie(grid, filters, quick):
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    n_frames = 60 if quick else 330
    lg_lo, lg_hi = np.log10(T_LO), np.log10(T_HI)
    # cosine-eased out-and-back: frame 0 == frame n (seamless loop)
    ph = np.arange(n_frames) / n_frames
    temps = 10 ** (lg_lo + (lg_hi - lg_lo) * (0.5 - 0.5 * np.cos(2 * np.pi * ph)))

    # whole-sweep reference for scales
    sweep = np.geomspace(T_LO, T_HI, 90)
    res = [fwd(grid, filters, t) for t in sweep]
    bmag = {b: np.array([r.magnitudes[b] for r in res]) for b in BANDS}
    bv = bmag["B"] - bmag["V"]
    mv = bmag["V"]
    peak = np.array([lamflam(r)[1][(lamflam(r)[0] > 100) & (lamflam(r)[0] < 2500)].max() for r in res])

    vc.apply_style()
    fig = vc.new_figure(1280, 720)
    ax = fig.add_axes([0.07, 0.13, 0.62, 0.78])
    axc = fig.add_axes([0.765, 0.17, 0.215, 0.66])

    # --- main panel ---------------------------------------------------------
    ax.set_xscale("log")
    ax.set_xlim(110, 2200)
    ax.set_ylim(-0.2, 1.15)          # each SED is normalised to its own peak; zero sits above the band lane
    ax.set_xticks([200, 400, 700, 1000, 2000])
    ax.set_xticklabels(["200", "400", "700", "1000", "2000"])
    ax.minorticks_off()
    ax.set_yticks([])
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("relative flux")

    bandfill = {}
    ref = grid.wavelength if hasattr(grid, "wavelength") else None
    for f in filters:
        wv = f.wavelengths / 10.0
        tr = f.transmission / f.transmission.max()
        pv = np.sum(wv * tr) / np.sum(tr)
        poly = ax.fill_between(wv, 0, 0.12 * tr, transform=ax.get_xaxis_transform(),
                               color=vc.INK, alpha=0.05, lw=0, zorder=1)
        ax.plot(wv, 0.12 * tr, transform=ax.get_xaxis_transform(), color=vc.INK_3, lw=0.6, zorder=1)
        ax.text(pv, 0.012, f.name, transform=ax.get_xaxis_transform(), color=vc.INK_2,
                ha="center", va="bottom", fontsize=10)
        bandfill[f.name] = poly

    w0, y0 = lamflam(res[0])
    y0 = y0 / y0[(w0 > 110) & (w0 < 2200)].max()
    arts = glow(ax, w0, y0, vc.star_color(T_LO))
    under = [ax.fill_between(w0, 0, y0, color="w", alpha=0.0, lw=0)]
    counter = ax.text(0.02, 0.95, "", transform=ax.transAxes, color=vc.INK_2, fontsize=14,
                      va="top", ha="left")

    # --- companion CMD ------------------------------------------------------
    pts = np.column_stack([bv, mv]).reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    cols = [vc.star_color(t, 1.0) + (0.40,) for t in sweep[:-1]]
    axc.add_collection(LineCollection(segs, colors=cols, linewidths=2.0, capstyle="round"))
    axc.set_xlim(bv.min() - 0.15, bv.max() + 0.15)
    axc.set_ylim(mv.max() + 1.0, mv.min() - 1.0)
    axc.set_xlabel("B – V")
    axc.set_ylabel("M$_V$")
    axc.tick_params(labelsize=10)
    halo = [axc.plot([], [], "o", color="w", ms=m, alpha=a, mew=0, zorder=4)[0]
            for m, a in ((22, 0.06), (14, 0.12), (8, 0.25))]
    dot, = axc.plot([], [], "o", color="w", ms=5, mew=0, zorder=5)

    bmin = {b: bmag[b].min() for b in BANDS}
    bmax = {b: bmag[b].max() for b in BANDS}
    cache = {}

    def update(i):
        t = float(temps[i])
        r = fwd(grid, filters, t)
        c = vc.star_color(t)
        x, y = lamflam(r)
        y = y / y[(x > 110) & (x < 2200)].max()
        under[0].remove()
        under[0] = ax.fill_between(x, 0, y, color=c, alpha=0.10, lw=0, zorder=2)
        for k, a in enumerate(arts):
            a.set_data(x, y)
            a.set_color(lighten(c, 0.25) if k == 3 else c)
        for b in BANDS:
            frac = (bmax[b] - r.magnitudes[b]) / (bmax[b] - bmin[b])
            bandfill[b].set_alpha(0.04 + 0.30 * float(np.clip(frac, 0, 1)) ** 1.2)
            bandfill[b].set_facecolor(c)
        cx, cy = r.magnitudes["B"] - r.magnitudes["V"], r.magnitudes["V"]
        for h in halo:
            h.set_data([cx], [cy]); h.set_color(c)
        dot.set_data([cx], [cy]); dot.set_color(lighten(c, 0.3))
        counter.set_text(f"Teff {int(round(t, -1))} K")

    paths = vc.save_movie(fig, update, n_frames, "star_heating", fps=30)
    plt.close(fig)
    return paths


# ---------------------------------------------------------------------------
# 2. CMD from a grid sweep
# ---------------------------------------------------------------------------

def make_cmd(grid, filters):
    import matplotlib.pyplot as plt

    vc.apply_style()
    fig = vc.new_figure(1280, 720)
    ax = fig.add_axes([0.08, 0.12, 0.88, 0.82])
    loggs = [4.5, 4.0, 3.0, 2.0, 1.0]
    temps = np.arange(3500, 30001, 100.0)
    ends = {}
    for lg in loggs:
        R = radius_for_logg(lg)
        pts = []
        for t in temps:
            if not safe(grid, t, lg):
                continue
            r = fwd(grid, filters, t, lg, R=R)
            pts.append((t, r.magnitudes["B"] - r.magnitudes["V"], r.magnitudes["V"]))
        a = np.array(pts)
        ax.plot(a[:, 1], a[:, 2], color=vc.INK_3, lw=0.7, zorder=1)
        d = a[np.isclose(a[:, 0] % np.where(a[:, 0] < 13000, 500.0, 1000.0), 0.0)]
        ax.scatter(d[:, 1], d[:, 2], s=14, c=[vc.star_color(t) for t in d[:, 0]], lw=0, zorder=3)
        ends[lg] = a

    ax.invert_yaxis()
    ax.set_xlabel("B – V")
    ax.set_ylabel("M$_V$  (1 M$_\\odot$ at each log g)")
    # labels at the hot end of every line
    for lg, a in ends.items():
        ax.annotate(f"log g {lg:g}", (a[-1, 1], a[-1, 2]), xytext=(-8, 0), textcoords="offset points",
                    color=vc.INK_2, fontsize=10, va="center", ha="right")
    # a few anchors on the log g 4.5 line, labelled underneath it
    a = ends[4.5]
    for tt in (5000, 10000, 20000):
        k = np.argmin(abs(a[:, 0] - tt))
        ax.scatter(a[k, 1], a[k, 2], s=40, facecolor="none", edgecolor=vc.INK_2, lw=0.6, zorder=4)
        ax.annotate(f"{tt:,} K".replace(",", " "), (a[k, 1], a[k, 2]), xytext=(0, -14),
                    textcoords="offset points", color=vc.INK, fontsize=11, ha="center", va="top")
    ax.margins(x=0.06, y=0.06)
    return vc.save_still(fig, "cmd_grid")


# ---------------------------------------------------------------------------
# 3. interpolation between two nodes
# ---------------------------------------------------------------------------

def make_between(grid, filters):
    import matplotlib.pyplot as plt

    def lf(r):
        x, y = lamflam(r)
        return x, y * 1e5

    lg = 4.5
    t0, t1 = 13000.0, 14000.0                    # adjacent Teff nodes, populated cell
    assert t0 in grid.teff_grid and t1 in grid.teff_grid
    A, B = fwd(grid, filters, t0, lg), fwd(grid, filters, t1, lg)
    mids = [fwd(grid, filters, t, lg) for t in (13250.0, 13500.0, 13750.0)]

    vc.apply_style()
    fig = vc.new_figure(1280, 720)
    ax = fig.add_axes([0.07, 0.12, 0.55, 0.82])
    zx = fig.add_axes([0.70, 0.12, 0.17, 0.82])
    win = (340, 420)

    def draw(a):
        for r in (A, B):
            x, y = lf(r)
            a.plot(x, y, color=vc.INK_2, lw=0.9, alpha=0.5, zorder=2)
        for r, main in zip(mids, (False, True, False)):
            x, y = lf(r)
            c = vc.star_color(r.teff)
            if main:
                for w, al in ((9, 0.06), (5, 0.12), (2.6, 0.25)):
                    a.plot(x, y, color=c, lw=w, alpha=al, zorder=3.5)
                a.plot(x, y, color=lighten(c, 0.25), lw=1.5, zorder=4)
            else:
                a.plot(x, y, color=c, lw=0.9, alpha=0.55, zorder=3)

    draw(ax)
    ax.set_xscale("log")
    ax.set_xlim(110, 1100)
    ax.set_ylim(0, None)
    ax.set_xticks([200, 400, 700, 1000]); ax.set_xticklabels(["200", "400", "700", "1000"])
    ax.minorticks_off()
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("flux at 10 pc  ($\\lambda F_\\lambda$, $10^{-5}$ cgs)")
    ax.axvspan(*win, color=vc.INK, alpha=0.06, lw=0, zorder=0)

    draw(zx)
    sl = [lf(r) for r in (A, B, *mids)]
    inwin = [y[(x > win[0]) & (x < win[1])] for x, y in sl]
    lo, hi = min(v.min() for v in inwin), max(v.max() for v in inwin)
    zx.set_xlim(*win)
    zx.set_ylim(lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))
    zx.set_xlabel("nm  (Balmer jump)")
    zx.set_yticks([])
    # direct labels at the right edge of the zoom, no legend box
    for r, txt, col in ((A, "13 000 K node", vc.INK_2), (mids[1], "13 500 K  Hermite", lighten(vc.star_color(13500), 0.25)),
                        (B, "14 000 K node", vc.INK_2)):
        x, y = lf(r)
        zx.annotate(txt, (win[1], y[np.argmin(abs(x - win[1]))]), xytext=(6, 0),
                    textcoords="offset points", color=col, fontsize=10, va="center", annotation_clip=False)
    return vc.save_still(fig, "between_nodes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", choices=["movie", "cmd", "between"])
    args = ap.parse_args()
    if args.quick:
        print("quick run, written to", vc.set_output_dir(tempfile.mkdtemp(prefix="gallery_quick_")))
    t0 = time.time()
    warnings.simplefilter("error", UserWarning)      # an empty-cell warning must be loud
    grid, filters = load_all()
    for name, fn in (("between", lambda: make_between(grid, filters)),
                     ("cmd", lambda: make_cmd(grid, filters)),
                     ("movie", lambda: make_movie(grid, filters, args.quick))):
        if args.only in (None, name):
            out = fn()
            print(name, "->", out)
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
