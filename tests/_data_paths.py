"""Locations of the SED_Tools data products the tests use.

Each path can be set with an environment variable; otherwise the first
existing candidate is used (SED_Tools' own cache, the older ``~/SED_Tools``
layout).
"""

from __future__ import annotations

import os

_ROOTS = [
    os.path.expanduser("~/.sed_tools/data"),
    os.path.expanduser("~/SED_Tools/data"),
]


def _find(env: str, *parts: str) -> str:
    if os.environ.get(env):
        return os.environ[env]
    candidates = [os.path.join(root, *parts) for root in _ROOTS]
    return next((c for c in candidates if os.path.exists(c)), candidates[0])


KURUCZ_DIR  = _find("CC_KURUCZ_DIR", "stellar_models", "Kurucz2003all")
JOHNSON_DIR = _find("CC_JOHNSON_DIR", "filters", "Generic", "Johnson")
VEGA_SED    = _find("CC_VEGA_SED", "stellar_models", "vega_flam.csv")
