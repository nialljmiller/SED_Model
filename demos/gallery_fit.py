"""
demos/gallery_fit.py
====================
Hero movie: a real ``run_inverse`` fit, watched converging.

Synthetic photometry (Johnson UBVRI + 2MASS JHKs) is forward-modelled for a
known F/G star, 0.03 mag of noise is added, and ``run_inverse`` recovers
(Teff, logg, [M/H]).  The walkers start scattered over a wide part of the
prior box (``init='ball'`` with a broad scatter, ``n_burn=0`` so nothing is
discarded) and the movie replays the retained chain step by step:

* big panel   one model SED per walker (coloured by the physical colour of
              that walker's Teff) against the observed photometry (orange);
* small panel the walkers in the Teff-[M/H] plane, with the history they
              leave behind.

Every frame is a state of the real chain; nothing is interpolated.  The
fixed radius and distance are the "known" inputs of the demo.

Usage
-----
    python demos/gallery_fit.py            # full movie (1280x720, ~12 s)
    python demos/gallery_fit.py --quick    # fast low-res test (written to a temp dir)

Needs the Kurucz2003all grid, the Johnson and 2MASS filters and the Vega SED
(see ``viz_common.data_paths``).
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import viz_common as vc  # noqa: E402

from sed_model import (  # noqa: E402
    load_grid, load_filters_from_instrument_dir, run_forward, run_inverse,
)
from sed_model.params import fit_params_from_grid  # noqa: E402

# --- the known star and the observing setup -------------------------------
TRUE = dict(teff=6150.0, logg=4.30, meta=-0.20)
R = 6.957e10                    # 1 R_sun (cm)
D = 100.0 * 3.0857e18           # 100 pc (cm)
SIGMA = 0.03                    # mag
SYSTEM = "Vega"
BOX = dict(teff=(4200.0, 9000.0), logg=(2.0, 5.0), meta=(-1.5, 0.5))   # prior box
START = dict(teff=6900.0, logg=3.5, meta=-0.6)                         # centre of the starting cloud
SCATTER = 0.30                  # starting scatter, fraction of each range
N_WALKERS = 48
N_STEPS = 160                   # chain length kept (n_burn=0)
SHOW_STEPS = 100                # steps replayed at easing speed (the rest plays slowly while the panel zooms)


def build_problem():
    p = vc.require_data()
    grid = load_grid(p["grid"])
    vega = p["vega"]
    joh = {f.name: f for f in load_filters_from_instrument_dir(p["johnson"], vega_sed_path=vega)}
    two = load_filters_from_instrument_dir(p["2mass"], vega_sed_path=vega)
    filters = [joh[n] for n in "UBVRI"] + sorted(two, key=lambda f: f.name != "J")  # J, H, Ks
    return grid, filters


def band_lambda(f) -> float:
    """Transmission-weighted effective wavelength (Angstrom)."""
    w, t = np.asarray(f.wavelengths, float), np.asarray(f.transmission, float)
    return float(np.trapezoid(w * t, w) / np.trapezoid(t, w))


def run_fit(grid, filters, seed=7):
    names = [f.name for f in filters]
    fwd = run_forward(teff=TRUE["teff"], logg=TRUE["logg"], meta=TRUE["meta"], R=R, d=D,
                      grid=grid, filters=filters, mag_system=SYSTEM)
    rng = np.random.default_rng(3)
    noise = rng.normal(0.0, SIGMA, len(names))
    mags = np.array([fwd.magnitudes[n] for n in names]) + noise
    sig = np.full(len(names), SIGMA)
    fp = fit_params_from_grid(grid, teff=BOX["teff"], logg=BOX["logg"], meta=BOX["meta"], d_cm=D)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")       # a deliberately short chain
        res = run_inverse(list(mags), list(sig), names, R, grid, filters, fit_params=fp,
                          mag_system=SYSTEM, n_walkers=N_WALKERS, n_steps=N_STEPS, n_burn=0,
                          n_thin=1, init="ball", p0_centre=START, p0_scatter=SCATTER,
                          seed=seed, progress=False, verbose=False)
    chain = res.samples.reshape(N_STEPS, N_WALKERS, -1)       # (step, walker, param): flat order is step-major
    # observed flux points: the true SED at each band's effective wavelength, scaled by the noise drawn
    lam = np.array([band_lambda(f) for f in filters])
    f_true = np.interp(lam, fwd.wavelengths, fwd.observed_flux)
    obs = lam * f_true * 10.0 ** (-0.4 * noise)
    err = obs * SIGMA * 0.4 * np.log(10.0)
    return res, chain, lam, obs, err


def step_schedule(n_frames, n_lead, n_hold, n_show):
    """Chain step shown on each frame: hold, ease in over the run, hold the end."""
    n_move = n_frames - n_lead - n_hold
    t = np.linspace(0.0, 1.0, n_move)
    s = np.round((n_show - 1) * t ** 1.3).astype(int)
    tail = np.round(np.linspace(n_show - 1, N_STEPS - 1, n_hold)).astype(int)   # converged chain keeps running
    return np.concatenate([np.zeros(n_lead, int), s, tail])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--quick", action="store_true", help="low-res short test, written to a temp dir")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    grid, filters = build_problem()
    res, chain, lam, obs, err = run_fit(grid, filters)
    print(f"chain {chain.shape}: Teff {res.samples[-N_WALKERS:, 0].mean():.0f} K at the end "
          f"(truth {TRUE['teff']:.0f})")

    if args.quick:
        W, H, fps, secs = 640, 360, 15, 5
        outdir = Path(tempfile.mkdtemp(prefix="gallery_fit_"))
    else:
        W, H, fps, secs = 1280, 720, 30, 12.5
        outdir = vc.ASSETS
    n_frames = int(round(fps * secs))
    n_lead, n_hold = int(0.5 * fps), int(2.4 * fps)
    fade = int(0.35 * fps)
    sched = step_schedule(n_frames, n_lead, n_hold, SHOW_STEPS)

    # model SEDs for every step that will be shown (one forward model per walker)
    sig_t = chain[SHOW_STEPS:, :, 0].std()
    mu_t = chain[SHOW_STEPS:, :, 0].mean()
    sig_m = chain[SHOW_STEPS:, :, 2].std()
    mu_m = chain[SHOW_STEPS:, :, 2].mean()
    zoom_t = (mu_t + 5 * sig_t, mu_t - 5 * sig_t)          # final window of the small panel
    zoom_m = (max(BOX["meta"][0], mu_m - 4 * sig_m), min(BOX["meta"][1], mu_m + 4 * sig_m))
    cache: dict[int, np.ndarray] = {}
    wave = None

    def seds(step):
        nonlocal wave
        if step not in cache:
            rows = []
            for th in chain[step]:
                r = run_forward(teff=th[0], logg=th[1], meta=th[2], R=R, d=D, grid=grid,
                                filters=filters[:1], mag_system=SYSTEM)
                wave = r.wavelengths
                rows.append(r.wavelengths * r.observed_flux)
            cache[step] = np.array(rows)
        return cache[step]

    seds(0)
    keep = (wave > 2500) & (wave < 30000)
    micron = wave[keep] / 1e4

    vc.apply_style()
    fig = vc.new_figure(W, H)
    s = W / 1280.0
    ax = fig.add_axes([0.055, 0.115, 0.625, 0.80])
    ax2 = fig.add_axes([0.745, 0.115, 0.215, 0.80])

    ymid = np.log10(obs.mean())
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.30, 2.6)
    ax.set_ylim(10 ** (ymid - 0.95), 10 ** (ymid + 0.55))
    ax.set_xticks([0.4, 0.6, 1.0, 2.0])
    ax.set_xticklabels(["0.4", "0.6", "1", "2"])
    ax.minorticks_off()
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.tick_params(labelsize=11 * s, length=3)

    ax2.set_xlim(*BOX["teff"][::-1])           # hot on the left, as in an HR diagram
    ax2.set_ylim(*BOX["meta"])
    ax2.set_xticks([])
    ax2.set_yticks([])
    for sp in ("left", "bottom"):
        ax2.spines[sp].set_color(vc.INK_3)

    # walker SED bundle
    cols0 = [vc.star_color(t) for t in chain[0][:, 0]]
    lines = [ax.plot(micron, np.full(keep.sum(), np.nan), lw=0.9 * s, color=c, alpha=0.30,
                     solid_capstyle="round")[0] for c in cols0]
    # observed photometry (drawn last, on top)
    ax.errorbar(lam / 1e4, obs, yerr=err, fmt="o", ms=7.5 * s, mfc=vc.ORANGE, mec=vc.INK,
                mew=0.8 * s, ecolor=vc.ORANGE, elinewidth=1.4 * s, capsize=3 * s, zorder=5)

    # Teff-[M/H] panel
    hist = ax2.scatter([], [], s=5 * s * s, c=vc.INK_2, alpha=0.10, linewidths=0, zorder=1)
    trail = ax2.scatter([], [], s=9 * s * s, linewidths=0, zorder=2)
    now = ax2.scatter([], [], s=30 * s * s, edgecolors=vc.SURFACE, linewidths=0.5 * s, zorder=3)

    counter = ax.text(0.985, 0.955, "", transform=ax.transAxes, ha="right", va="top",
                      fontsize=19 * s, color=vc.INK)

    def draw(i):
        step = int(sched[i])
        th = chain[step]
        cols = np.array([vc.star_color(t) for t in th[:, 0]])
        a_vis = min(1.0, (i + 1) / fade, (n_frames - i) / fade)      # fade in/out for a clean loop
        y = seds(step)[:, keep]
        for ln, row, c in zip(lines, y, cols):
            ln.set_ydata(row)
            ln.set_color(c)
            ln.set_alpha(0.30 * a_vis)
        # history cloud, recent trail, current walkers
        hist.set_offsets(chain[:step + 1].reshape(-1, 3)[:, [0, 2]])
        recent = chain[max(0, step - 4):step]
        if len(recent):
            trail.set_offsets(recent.reshape(-1, 3)[:, [0, 2]])
            trail.set_facecolor(np.tile(np.c_[cols, np.full(len(cols), 0.35)], (len(recent), 1)))
        else:
            trail.set_offsets(np.empty((0, 2)))
        now.set_offsets(th[:, [0, 2]])
        now.set_facecolor(np.c_[cols, np.full(len(cols), 0.95 * a_vis)])
        now.set_edgecolor(vc.SURFACE)
        z = np.clip((i - (n_frames - n_hold) - 0.08 * n_hold) / (0.65 * n_hold), 0.0, 1.0)
        z = z * z * (3 - 2 * z)                                 # smoothstep zoom of the small panel
        ax2.set_xlim(BOX["teff"][1] + z * (zoom_t[0] - BOX["teff"][1]),
                     BOX["teff"][0] + z * (zoom_t[1] - BOX["teff"][0]))
        ax2.set_ylim(BOX["meta"][0] + z * (zoom_m[0] - BOX["meta"][0]),
                     BOX["meta"][1] + z * (zoom_m[1] - BOX["meta"][1]))
        counter.set_text(f"Teff {np.median(th[:, 0]):.0f} K")
        counter.set_alpha(a_vis)

    paths = vc.save_movie(fig, draw, n_frames, "sed_fit", fps=fps, outdir=outdir,
                          gif_width=720 if not args.quick else 480)

    # still: mid-convergence, the messy bundle already pulling onto the data
    draw(int(np.argmax(sched >= 9)))
    still = vc.save_still(fig, "sed_fit_still", outdir=outdir)
    for p in paths + [still]:
        print(f"wrote {p}  ({p.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
