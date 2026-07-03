"""
sed_model._cc_ext
=================
Single source-of-truth accessor for the compiled Fortran extension.

Importing from this module instead of calling the extension directly lets
filters.py and forward.py share one import path with one documented failure
policy, rather than each defining their own _get_cc_api() with divergent
error handling.
"""

from __future__ import annotations


def get_cc_api(required: bool = False):
    """Return the f2py-wrapped cc_api module object, or None if unavailable.

    Parameters
    ----------
    required : bool
        If True, raise ImportError when the extension is absent.
        If False (default), return None so the caller can fall back to a
        pure-Python implementation.
    """
    try:
        from . import cc_api as _ext
        return _ext.cc_api
    except Exception as exc:
        if required:
            raise ImportError(
                "The Fortran extension 'cc_api' is not built. "
                "Run 'make' in the SED_Model root directory."
            ) from exc
        return None
