"""
demos/viz_common.py
===================
Shared look and helpers for the gallery demos (``gallery_*.py``).

Everything here is presentation: a dark surface, a small validated palette,
the physical colour of a star from its temperature, data-path lookup, and
helpers that write stills and movies.  Import from a demo with::

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent))
    import viz_common as vc

Palette
-------
Three categorical hues on the dark surface (they pass the all-pairs colour
checks of the data-viz palette validator; a fourth hue does not, so a panel
never uses more than three).  Stars are drawn in their blackbody colour,
which is a physical encoding and not a categorical one.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

# --- surface and ink -------------------------------------------------------
SURFACE = "#121214"
INK = "#f0efec"        # primary text
INK_2 = "#a8a79e"      # secondary text, axes
INK_3 = "#55554f"      # grid, quiet marks

# --- categorical (fixed order) ---------------------------------------------
BLUE = "#3987e5"       # model
ORANGE = "#d95926"     # data / observed
AQUA = "#199e70"       # third series

ASSETS = Path(__file__).resolve().parent.parent / "docs" / "assets" / "gallery"


# ---------------------------------------------------------------------------
# Matplotlib style
# ---------------------------------------------------------------------------

def apply_style() -> None:
    """Set matplotlib rcParams for the gallery look (call once per script)."""
    import matplotlib as mpl

    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": INK_3,
        "axes.labelcolor": INK_2,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": INK_2,
        "ytick.color": INK_2,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "text.color": INK,
        "font.family": "sans-serif",
        "font.size": 12,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "legend.frameon": False,
        "lines.solid_capstyle": "round",
        "figure.dpi": 100,
        "savefig.dpi": 100,
    })


def new_figure(width_px: int = 1280, height_px: int = 720):
    """A figure of an exact pixel size (default 16:9, 1280x720)."""
    import matplotlib.pyplot as plt

    return plt.figure(figsize=(width_px / 100, height_px / 100), dpi=100)


# ---------------------------------------------------------------------------
# The colour of a star
# ---------------------------------------------------------------------------

def _cie_xyz(wave_nm: np.ndarray):
    """CIE 1931 colour-matching functions (Wyman, Sloan & Shirley 2013)."""
    w = wave_nm

    def g(x, mu, s1, s2):
        s = np.where(x < mu, s1, s2)
        return np.exp(-0.5 * ((x - mu) / s) ** 2)

    x = 1.056 * g(w, 599.8, 37.9, 31.0) + 0.362 * g(w, 442.0, 16.0, 26.7) - 0.065 * g(w, 501.1, 20.4, 26.2)
    y = 0.821 * g(w, 568.8, 46.9, 40.5) + 0.286 * g(w, 530.9, 16.3, 31.1)
    z = 1.217 * g(w, 437.0, 11.8, 36.0) + 0.681 * g(w, 459.0, 26.0, 13.8)
    return x, y, z


def star_color(teff: float, brightness: float = 1.0) -> tuple[float, float, float]:
    """sRGB colour (0-1) of a blackbody at ``teff`` kelvin.

    The chromaticity is normalised so the brightest channel is 1, then
    ``brightness`` scales it.  Cool stars come out orange-red, the Sun
    near white, hot stars blue-white.
    """
    wave = np.linspace(380.0, 780.0, 401)                    # nm
    lam = wave * 1e-9
    h, c, k = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    b = 1.0 / (lam ** 5 * np.expm1(h * c / (lam * k * teff)))
    x, y, z = _cie_xyz(wave)
    xyz = np.array([np.trapezoid(b * x, wave), np.trapezoid(b * y, wave), np.trapezoid(b * z, wave)])
    xyz = xyz / xyz.sum()
    X, Y, Z = xyz[0] / xyz[1], 1.0, xyz[2] / xyz[1]
    rgb = np.array([3.2406 * X - 1.5372 * Y - 0.4986 * Z,
                    -0.9689 * X + 1.8758 * Y + 0.0415 * Z,
                    0.0557 * X - 0.2040 * Y + 1.0570 * Z])
    rgb = np.clip(rgb, 0.0, None)
    rgb = rgb / rgb.max()
    rgb = np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb ** (1 / 2.4) - 0.055)
    return tuple(np.clip(rgb * brightness, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Data lookup (same order as tests/_data_paths.py)
# ---------------------------------------------------------------------------

_ROOTS = [os.path.expanduser("~/.sed_tools/data"), os.path.expanduser("~/SED_Tools/data")]


def _find(env: str, *parts: str) -> str:
    if os.environ.get(env):
        return os.environ[env]
    cands = [os.path.join(r, *parts) for r in _ROOTS]
    return next((c for c in cands if os.path.exists(c)), cands[0])


def data_paths() -> dict[str, str]:
    """Locations of the SED_Tools products the gallery uses.

    Overridable with ``CC_KURUCZ_DIR``, ``CC_JOHNSON_DIR``, ``CC_GAIA_DIR``,
    ``CC_2MASS_DIR`` and ``CC_VEGA_SED``.
    """
    return {
        "grid": _find("CC_KURUCZ_DIR", "stellar_models", "Kurucz2003all"),
        "johnson": _find("CC_JOHNSON_DIR", "filters", "Generic", "Johnson"),
        "gaia": _find("CC_GAIA_DIR", "filters", "GAIA", "GAIA"),
        "2mass": _find("CC_2MASS_DIR", "filters", "2MASS", "2MASS"),
        "vega": _find("CC_VEGA_SED", "stellar_models", "vega_flam.csv"),
    }


def require_data() -> dict[str, str]:
    """``data_paths()`` or exit with a message naming what is missing."""
    p = data_paths()
    missing = [f"{k}: {v}" for k, v in p.items() if k in ("grid", "johnson", "vega") and not os.path.exists(v)]
    if missing:
        raise SystemExit("Missing SED_Tools data:\n  " + "\n  ".join(missing)
                         + "\nDownload it with SED_Tools or set the CC_* variables (see docs/testing.md).")
    return p


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def save_still(fig, name: str, outdir: Path | None = None) -> Path:
    """Save ``fig`` as ``<outdir>/<name>.png`` (default ``docs/assets/gallery``)."""
    outdir = Path(outdir or ASSETS)
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{name}.png"
    fig.savefig(path, facecolor=SURFACE)
    return path


def _ffmpeg() -> str:
    exe = shutil.which("ffmpeg")
    if exe is None:
        raise SystemExit("ffmpeg is needed to write movies (not found on PATH).")
    return exe


def save_movie(fig, update, n_frames: int, name: str, fps: int = 30,
               outdir: Path | None = None, gif_width: int = 720, gif_fps: int = 15,
               formats: tuple[str, ...] = ("mp4", "gif")) -> list[Path]:
    """Render ``update(i)`` for ``i in range(n_frames)`` and encode the movie.

    ``update`` redraws ``fig`` for frame ``i``.  Frames are written as PNGs to a
    temporary directory and encoded with ffmpeg: H.264 MP4 (yuv420p, plays
    everywhere) and a palette-optimised GIF at ``gif_width`` pixels for READMEs.
    """
    outdir = Path(outdir or ASSETS)
    outdir.mkdir(parents=True, exist_ok=True)
    ff = _ffmpeg()
    written: list[Path] = []
    with tempfile.TemporaryDirectory() as tmp:
        for i in range(n_frames):
            update(i)
            fig.savefig(Path(tmp) / f"f{i:05d}.png", facecolor=SURFACE)
        src = str(Path(tmp) / "f%05d.png")
        if "mp4" in formats:
            out = outdir / f"{name}.mp4"
            subprocess.run([ff, "-y", "-loglevel", "error", "-framerate", str(fps), "-i", src,
                            "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
                            "-movflags", "+faststart", str(out)], check=True)
            written.append(out)
        if "gif" in formats:
            out = outdir / f"{name}.gif"
            vf = (f"fps={gif_fps},scale={gif_width}:-1:flags=lanczos,"
                  "split[a][b];[a]palettegen=max_colors=96:stats_mode=diff[p];"
                  "[b][p]paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle")
            subprocess.run([ff, "-y", "-loglevel", "error", "-framerate", str(fps), "-i", src,
                            "-vf", vf, str(out)], check=True)
            written.append(out)
    return written
