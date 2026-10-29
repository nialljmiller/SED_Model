"""
demos/gallery_grid.py
=====================
Gallery pieces about the atmosphere grid and the posterior.

    python demos/gallery_grid.py            # all three assets
    python demos/gallery_grid.py --quick    # fast smoke test (written to a temp dir)

Writes to ``docs/assets/gallery``:

* ``grid_walk.mp4`` / ``.gif``  a walk through the 3-D Kurucz grid: a (Teff, logg)
  map with populated and empty nodes, and the interpolated SED at the query point
  with the surrounding node SEDs as ghosts.
* ``posterior_corner.png``      a corner plot from a real ``run_inverse`` fit of
  synthetic photometry (Teff, logg, [M/H], A_V).
* ``hermite_vs_linear.png``     'linear' against 'hermite_bounded' across grid cells.

Needs the released ``sed-model`` wheel, matplotlib, ffmpeg, and the SED_Tools data
(see ``viz_common.require_data``).
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

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, to_rgb  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

import sed_model as sm  # noqa: E402

RSUN = 6.957e10
PC = 3.0857e18


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _mix(c1, c2, f):
    a, b = np.array(to_rgb(c1)), np.array(to_rgb(c2))
    return tuple(a * (1 - f) + b * f)


def _populated(grid) -> np.ndarray:
    """(n_teff, n_logg, n_meta) bool: node carries a model spectrum."""
    return grid.flux.max(axis=-1) > 0


def _surface_flux(grid, filters, teff, logg, meta, method="hermite_bounded"):
    r = sm.run_forward(teff=teff, logg=logg, meta=meta, R=RSUN, d=10 * PC, grid=grid,
                       filters=filters, mag_system="AB", interp_method=method)
    return r


# ---------------------------------------------------------------------------
# 1. grid walk
# ---------------------------------------------------------------------------

def _spline_path(way: np.ndarray, n: int) -> np.ndarray:
    """Natural cubic spline through ``way`` (k, 3), sampled at constant speed."""
    k = len(way)
    # chord-length parameter in a normalised space (log Teff range ~1, logg ~5, meta ~3)
    scale = np.array([1.0 / 0.35, 1.0 / 2.5, 1.0 / 1.5])
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(way * scale, axis=0), axis=1))]
    # tridiagonal solve for second derivatives (natural)
    M = np.zeros((k, 3))
    h = np.diff(d)
    A = np.zeros((k, k))
    rhs = np.zeros((k, 3))
    A[0, 0] = A[-1, -1] = 1.0
    for i in range(1, k - 1):
        A[i, i - 1], A[i, i], A[i, i + 1] = h[i - 1], 2 * (h[i - 1] + h[i]), h[i]
        rhs[i] = 6 * ((way[i + 1] - way[i]) / h[i] - (way[i] - way[i - 1]) / h[i - 1])
    M = np.linalg.solve(A, rhs)

    def ev(s):
        i = np.clip(np.searchsorted(d, s) - 1, 0, k - 2)
        hh = h[i][:, None]
        a = (d[i + 1] - s)[:, None] / hh
        b = (s - d[i])[:, None] / hh
        return (a * way[i] + b * way[i + 1]
                + ((a ** 3 - a) * M[i] + (b ** 3 - b) * M[i + 1]) * hh ** 2 / 6)

    # constant speed (arc length), then ease in/out in time
    fine = np.linspace(0, d[-1], 4000)
    pts = ev(fine)
    seg = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts * scale, axis=0), axis=1))]
    u = np.linspace(0, 1, n)
    u = u * u * (3 - 2 * u) * 0.35 + u * 0.65            # gentle ease, never stalls
    s_of = np.interp(u * seg[-1], seg, fine)
    return ev(s_of)


def path_through_grid(grid, n_frames: int) -> np.ndarray:
    """Waypoints: cool giant -> Sun -> A star -> hot dwarf, all in populated cells."""
    way = np.array([
        [np.log10(4100.0), 1.9, -1.6],
        [np.log10(4900.0), 3.2, -0.9],
        [np.log10(5778.0), 4.44, 0.0],
        [np.log10(7300.0), 4.05, 0.2],
        [np.log10(10500.0), 3.95, -0.6],
        [np.log10(16500.0), 4.25, -0.2],
        [np.log10(23000.0), 4.2, -1.2],
        [np.log10(30000.0), 4.55, -0.9],
        [np.log10(36000.0), 4.75, -0.4],
    ])
    p = _spline_path(way, n_frames)
    teff = 10 ** p[:, 0]
    bad = [i for i in range(n_frames)
           if grid.cell_has_empty_nodes(teff[i], p[i, 1], p[i, 2], stencil=True)]
    if bad:
        pts = [(round(float(teff[i])), round(float(p[i, 1]), 2), round(float(p[i, 2]), 2)) for i in bad[::max(1, len(bad) // 6)]]
        raise RuntimeError(f"path enters empty cells at {len(bad)} frames, e.g. {pts}")
    return np.column_stack([teff, p[:, 1], p[:, 2]])


def make_grid_walk(paths, n_frames: int, outdir: Path | None, formats=("mp4", "gif")):
    vc.apply_style()
    grid = sm.load_grid(paths["grid"])
    filters = sm.load_filters([paths["johnson"] + "/V.dat"], paths["vega"])
    path = path_through_grid(grid, n_frames)
    pop = _populated(grid)
    tg, lg, mg = grid.teff_grid, grid.logg_grid, grid.meta_grid
    full = pop.all(axis=2)                                   # (T, g) same for every [M/H]
    TT, LL = np.meshgrid(tg, lg, indexing="ij")

    fig = vc.new_figure(1280, 720)
    axm = fig.add_axes([0.075, 0.14, 0.32, 0.74])
    axs = fig.add_axes([0.475, 0.14, 0.50, 0.74])

    # ---- map (static layer) ----
    axm.set_xlim(np.log10(52000), np.log10(3200))
    axm.set_ylim(5.25, -0.25)
    axm.scatter(np.log10(TT[~full]), LL[~full], s=9, facecolors="none", edgecolors=vc.INK_3,
                linewidths=0.5, alpha=0.9, zorder=1)
    axm.scatter(np.log10(TT[full]), LL[full], s=9, color=vc.INK_2, alpha=0.9, linewidths=0, zorder=2)
    xt = [4000, 6000, 10000, 20000, 40000]
    axm.set_xticks(np.log10(xt), [f"{t // 1000}k" for t in xt])
    axm.set_yticks([0, 1, 2, 3, 4, 5])
    axm.set_xlabel("Teff  (K)", labelpad=6)
    axm.set_ylabel("log g", labelpad=6)
    axm.tick_params(labelsize=10)
    # whole route, faint
    axm.plot(np.log10(path[:, 0]), path[:, 1], color=vc.INK_3, lw=0.7, ls=(0, (2, 4)), zorder=1, alpha=0.9)

    # ---- sed axes (static layer) ----
    lam = grid.wavelengths * 1e-4                             # micron
    axs.set_xscale("log")
    axs.set_xlim(0.05, 4.0)
    axs.set_ylim(0, 1.22)
    axs.set_xticks([0.1, 0.3, 1, 3], ["0.1", "0.3", "1", "3"])
    axs.minorticks_off()
    axs.set_yticks([])
    axs.spines["left"].set_visible(False)
    axs.set_xlabel("wavelength  (µm)", labelpad=6)
    axs.tick_params(labelsize=10)

    dyn = {}                                                  # artists removed each frame
    n_trail = 70
    counter = fig.text(0.475, 0.915, "", fontsize=13, color=vc.INK_2, family="monospace", va="bottom")

    def lamF(flux):
        return lam * flux

    def smooth(y, k=5):
        return np.convolve(np.pad(y, k // 2, mode="edge"), np.ones(k) / k, mode="valid")

    def update(i):
        for a in dyn.values():
            for x in (a if isinstance(a, list) else [a]):
                x.remove()
        dyn.clear()
        T, g, m = path[i]
        col = vc.star_color(T)
        # trail
        j0 = max(0, i - n_trail)
        if i - j0 >= 2:
            pts = np.column_stack([np.log10(path[j0:i + 1, 0]), path[j0:i + 1, 1]])
            segs = np.stack([pts[:-1], pts[1:]], axis=1)
            alpha = np.linspace(0.0, 0.8, len(segs))
            rgba = [(*col, a) for a in alpha]
            lc = LineCollection(segs, colors=rgba, linewidths=1.6, zorder=3, capstyle="round")
            axm.add_collection(lc)
            dyn["trail"] = lc
        # current cell and its corners
        it = int(np.clip(np.searchsorted(tg, T) - 1, 0, len(tg) - 2))
        ig = int(np.clip(np.searchsorted(lg, g) - 1, 0, len(lg) - 2))
        x0, x1 = np.log10(tg[it]), np.log10(tg[it + 1])
        rect = Rectangle((x0, lg[ig]), x1 - x0, lg[ig + 1] - lg[ig], fill=False,
                         edgecolor=col, lw=0.8, alpha=0.7, zorder=3)
        axm.add_patch(rect)
        dyn["cell"] = rect
        cx = [x0, x1, x0, x1]
        cy = [lg[ig], lg[ig], lg[ig + 1], lg[ig + 1]]
        dyn["corners"] = axm.scatter(cx, cy, s=26, facecolors="none", edgecolors=col, linewidths=1.0, zorder=4)
        # glowing query point
        xq = np.log10(T)
        glow = []
        for s, a in ((1500, 0.05), (700, 0.09), (300, 0.16), (110, 0.3)):
            glow.append(axm.scatter([xq], [g], s=s, color=col, alpha=a, linewidths=0, zorder=5))
        glow.append(axm.scatter([xq], [g], s=26, color=_mix(col, "#ffffff", 0.5), linewidths=0, zorder=6))
        dyn["glow"] = glow

        # ---- SED ----
        r = _surface_flux(grid, filters, T, g, m)
        y = lamF(r.surface_flux)
        win = (lam > 0.05) & (lam < 4.0)
        norm = y[win].max()
        ghosts = []
        # corner nodes of the 3-D cell: 2 x 2 x 2 neighbours
        im = int(np.clip(np.searchsorted(mg, m) - 1, 0, len(mg) - 2))
        for a in (it, it + 1):
            for b in (ig, ig + 1):
                for c in (im, im + 1):
                    yn = lamF(grid.flux[a, b, c])
                    ghosts += axs.plot(lam, smooth(yn / norm), color=col, lw=0.7, alpha=0.22, zorder=2)
        dyn["ghosts"] = ghosts
        ys = smooth(y / norm)
        dyn["fill"] = axs.fill_between(lam, 0, ys, color=col, alpha=0.10, linewidth=0, zorder=3)
        dyn["halo"] = axs.plot(lam, ys, color=col, lw=7, alpha=0.10, zorder=4, solid_capstyle="round")
        dyn["line"] = axs.plot(lam, ys, color=col, lw=2.0, zorder=5)
        counter.set_text(f"{T:7.0f} K  ·  {g:4.2f}  ·  {m:+4.2f}")

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)           # a warning means a bad path
        written = vc.save_movie(fig, update, n_frames, "grid_walk", outdir=outdir, formats=formats)
    plt.close(fig)
    return written


# ---------------------------------------------------------------------------
# 2. posterior corner
# ---------------------------------------------------------------------------

def _blur(img, sigma):
    k = int(4 * sigma) + 1
    x = np.arange(-k, k + 1)
    w = np.exp(-0.5 * (x / sigma) ** 2)
    w /= w.sum()
    out = np.apply_along_axis(lambda v: np.convolve(v, w, mode="same"), 0, img)
    return np.apply_along_axis(lambda v: np.convolve(v, w, mode="same"), 1, out)


def run_fit(paths, quick: bool):
    grid = sm.load_grid(paths["grid"])
    j = sm.load_filters_from_instrument_dir(paths["johnson"], vega_sed_path=paths["vega"])
    ga = sm.load_filters_from_instrument_dir(paths["gaia"], vega_sed_path=paths["vega"])
    tm = sm.load_filters_from_instrument_dir(paths["2mass"], vega_sed_path=paths["vega"])
    pick = ([f for f in j if f.name in ("U", "B", "V", "R", "I")]
            + [f for f in ga if f.name in ("G", "Gbp", "Grp")]
            + [f for f in tm if f.name in ("J", "H", "Ks")])
    names = [f.name for f in pick]

    truth = dict(teff=6150.0, logg=4.15, meta=-0.30, a_v=0.55)
    R, d = 1.35 * RSUN, 85.0 * PC
    ext = sm.make_extinction_model(enabled=True, law="fitzpatrick99", r_v=3.1, a_v=truth["a_v"])
    fwd = sm.run_forward(teff=truth["teff"], logg=truth["logg"], meta=truth["meta"], R=R, d=d,
                         grid=grid, filters=pick, mag_system="AB", extinction=ext)
    sig = np.full(len(pick), 0.03)
    rng = np.random.default_rng(7)
    obs = np.array([fwd.magnitudes[n] for n in names]) + rng.normal(0, sig)

    fp = sm.fit_params_from_grid(grid, teff=(5000.0, 8000.0), logg=(3.0, 5.0), meta=(-1.2, 0.3),
                                 a_v=(0.0, 1.5), d_cm=d)
    res = sm.run_inverse(obs, sig, names, R=R, grid=grid, filters=pick, fit_params=fp,
                         extinction=ext, mag_system="AB",
                         n_walkers=24 if quick else 48, n_steps=500 if quick else 3000,
                         n_burn=150 if quick else 800, n_thin=1 if quick else 2,
                         seed=11, progress=False, verbose=False)
    return res, truth


def make_corner(paths, quick: bool, outdir: Path | None):
    vc.apply_style()
    res, truth = run_fit(paths, quick)
    keys = list(res.param_names)
    S = np.asarray(res.samples)
    labels = {"teff": "Teff (K)", "logg": "log g", "meta": "[M/H]", "a_v": "A$_V$"}
    n = len(keys)
    cmap = LinearSegmentedColormap.from_list(
        "post", [vc.SURFACE, _mix(vc.SURFACE, vc.BLUE, 0.35), vc.BLUE, _mix(vc.BLUE, "#ffffff", 0.55)])

    fig = vc.new_figure(1200, 1200)
    gs = fig.add_gridspec(n, n, left=0.105, right=0.975, bottom=0.095, top=0.975, wspace=0.07, hspace=0.07)
    lim = {}
    for k in keys:
        x = S[:, keys.index(k)]
        lo, hi = np.percentile(x, [0.15, 99.85])
        pad = 0.12 * (hi - lo)
        lim[k] = (min(lo - pad, truth[k]), max(hi + pad, truth[k]))
        if truth[k] < lo - pad or truth[k] > hi + pad:
            lim[k] = (min(lim[k][0], truth[k] - pad), max(lim[k][1], truth[k] + pad))

    for r in range(n):
        for c in range(r + 1):
            ax = fig.add_subplot(gs[r, c])
            kx, ky = keys[c], keys[r]
            ax.set_facecolor(vc.SURFACE)
            if r == c:
                x = S[:, c]
                h, e = np.histogram(x, bins=70, range=lim[kx], density=True)
                h = np.convolve(np.pad(h, 3, mode="edge"), np.ones(7) / 7, mode="valid")
                xc = 0.5 * (e[1:] + e[:-1])
                ax.fill_between(xc, 0, h, color=vc.BLUE, alpha=0.14, linewidth=0)
                ax.plot(xc, h, color=vc.BLUE, lw=1.0)
                ax.axvline(truth[kx], color=vc.ORANGE, lw=0.8)
                ax.set_xlim(*lim[kx])
                ax.set_ylim(0, h.max() * 1.1)
                ax.set_yticks([])
                ax.spines["left"].set_visible(False)
            else:
                H, xe, ye = np.histogram2d(S[:, c], S[:, r], bins=90, range=[lim[kx], lim[ky]])
                H = _blur(H, 1.8).T
                srt = np.sort(H.ravel())[::-1]
                cum = np.cumsum(srt) / srt.sum()
                lev = [srt[np.searchsorted(cum, q)] for q in (0.99, 0.95, 0.86, 0.68, 0.40)]
                lev = sorted(set(lev))
                xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
                ax.contourf(xc, yc, H, levels=lev + [H.max() * 1.001], cmap=cmap, antialiased=True)
                ax.contour(xc, yc, H, levels=[lev[1], lev[3]], colors=[_mix(vc.BLUE, "#ffffff", 0.35)],
                           linewidths=0.5, alpha=0.7)
                ax.axvline(truth[kx], color=vc.ORANGE, lw=0.8, alpha=0.9)
                ax.axhline(truth[ky], color=vc.ORANGE, lw=0.8, alpha=0.9)
                ax.set_xlim(*lim[kx])
                ax.set_ylim(*lim[ky])
            ax.xaxis.set_major_locator(plt.MaxNLocator(3, prune="both"))
            ax.tick_params(labelsize=10, length=2.5)
            if r == n - 1:
                ax.set_xlabel(labels[kx], labelpad=6, fontsize=13)
            else:
                ax.set_xticklabels([])
            if c == 0 and r > 0:
                ax.yaxis.set_major_locator(plt.MaxNLocator(3, prune="both"))
                ax.set_ylabel(labels[ky], labelpad=6, fontsize=13)
            elif r != c:
                ax.set_yticklabels([])
    path = vc.save_still(fig, "posterior_corner", outdir)
    plt.close(fig)
    return path, res, truth


# ---------------------------------------------------------------------------
# 3. hermite vs linear
# ---------------------------------------------------------------------------

def make_interp_compare(paths, outdir: Path | None):
    vc.apply_style()
    grid = sm.load_grid(paths["grid"])
    V = sm.load_filters([paths["johnson"] + "/V.dat"], paths["vega"])
    logg, meta = 4.0, 0.0
    lo, hi = 5000.0, 6500.0
    T = np.linspace(lo, hi, 1201)
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        mv = {m: np.array([_surface_flux(grid, V, t, logg, meta, m).magnitudes["V"] for t in T])
              for m in ("linear", "hermite_bounded")}
    diff = (mv["hermite_bounded"] - mv["linear"]) * 1000.0                   # mmag
    sl = {m: np.gradient(mv[m], T) * 1000.0 for m in mv}                     # mmag / K
    nodes = grid.teff_grid[(grid.teff_grid >= lo) & (grid.teff_grid <= hi)]

    fig = vc.new_figure(1200, 900)
    a1 = fig.add_axes([0.11, 0.55, 0.86, 0.40])
    a2 = fig.add_axes([0.11, 0.10, 0.86, 0.38], sharex=a1)
    for a in (a1, a2):
        for t in nodes:
            a.axvline(t, color=vc.INK_3, lw=0.5, alpha=0.8, zorder=0)
        a.set_xlim(lo, hi)
    a1.plot(T, diff, color=vc.BLUE, lw=1.6)
    a1.axhline(0, color=vc.INK_2, lw=0.6)
    a1.set_ylabel("hermite − linear  (mmag)", labelpad=8)
    plt.setp(a1.get_xticklabels(), visible=False)
    a2.plot(T, sl["linear"], color=vc.ORANGE, lw=1.4, drawstyle="default")
    a2.plot(T, sl["hermite_bounded"], color=vc.BLUE, lw=1.6)
    a2.set_ylabel("dV / dTeff  (mmag / K)", labelpad=8)
    a2.set_xlabel("Teff  (K)", labelpad=6)
    path = vc.save_still(fig, "hermite_vs_linear", outdir)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="fast smoke test, written to a temp dir")
    ap.add_argument("--only", choices=["walk", "corner", "interp"], default=None)
    a = ap.parse_args()
    paths = vc.require_data()
    outdir = Path(tempfile.mkdtemp(prefix="gallery_grid_")) if a.quick else None
    t0 = time.time()
    todo = [a.only] if a.only else ["walk", "corner", "interp"]
    out = []
    if "walk" in todo:
        out += make_grid_walk(paths, 45 if a.quick else 330, outdir)
    if "corner" in todo:
        out.append(make_corner(paths, a.quick, outdir)[0])
    if "interp" in todo:
        out.append(make_interp_compare(paths, outdir))
    for p in out:
        print(f"{p}  {p.stat().st_size / 1e6:.2f} MB")
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
