"""Smoke tests for the gallery scripts in demos/ (quick mode).

Each ``demos/gallery_*.py --quick`` renders a low-resolution version into a
temporary directory, so these tests never touch the committed assets.  They
skip when matplotlib, ffmpeg, the compiled extension or the SED_Tools data are
missing.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("matplotlib")
if shutil.which("ffmpeg") is None:
    pytest.skip("ffmpeg not on PATH", allow_module_level=True)

from _data_paths import JOHNSON_DIR, KURUCZ_DIR, VEGA_SED  # noqa: E402

for _p, _what in ((KURUCZ_DIR, "Kurucz2003all grid"), (JOHNSON_DIR, "Johnson filters"),
                  (VEGA_SED, "Vega SED")):
    if not os.path.exists(_p):
        pytest.skip(f"{_what} not found at {_p}", allow_module_level=True)

try:
    from sed_model import cc_api  # noqa: F401
except Exception:                                              # pragma: no cover
    pytest.skip("Fortran extension not built", allow_module_level=True)

DEMOS = Path(__file__).resolve().parent.parent / "demos"
sys.path.insert(0, str(DEMOS))


def test_star_colour_runs_from_orange_to_blue():
    import viz_common as vc

    cool, sun, hot = (vc.star_color(t) for t in (3000.0, 5772.0, 20000.0))
    assert cool[0] > cool[2]          # red beats blue when cool
    assert hot[2] > hot[0]            # blue beats red when hot
    assert min(sun) > 0.85            # the Sun is close to white
    assert all(0.0 <= c <= 1.0 for c in cool + sun + hot)


@pytest.mark.parametrize("script", ["gallery_forward", "gallery_extinction", "gallery_grid", "gallery_fit"])
def test_gallery_script_quick(script):
    env = dict(os.environ, MPLBACKEND="Agg", CC_KURUCZ_DIR=KURUCZ_DIR,
               CC_JOHNSON_DIR=JOHNSON_DIR, CC_VEGA_SED=VEGA_SED)
    out = subprocess.run([sys.executable, str(DEMOS / f"{script}.py"), "--quick"],
                         env=env, capture_output=True, text=True, timeout=600)
    assert out.returncode == 0, out.stderr[-2000:]
