"""
Lazy import cache

This is a pure Python file with no dependencies so it can be used in setup.py.
"""
import os
import hashlib

from ..env import SAGE_LIB
from sage.misc.dot_sage import dot_sage


def get_cache_file():
    """
    Return the canonical filename for caching names of lazily imported
    modules.

    EXAMPLES::

        sage: from sage.misc.dot_sage
        sage: from sage.misc.lazy_import_cache import get_cache_file
        sage: str(get_cache_file())
        '...-lazy_import_cache.pickle'
        sage: get_cache_file().parent == dot_sage() / 'cache'
        True
    """
    mangled = hashlib.sha256(os.path.realpath(SAGE_LIB).encode('utf-8')).hexdigest()
    return dot_sage() / 'cache' / ("%s-lazy_import_cache.pickle" % mangled)
