"""
Operations for LibGAP Elements

GAP functions for which several methods can be available are called
operations, so GAP ``Size`` is an example of an operation. This module
is for inspecting GAP operations from Python. In particular, it can
list the operations that take a particular LibGAP element as first
argument. This is used in tab completion, where Python ``x.[TAB]``
lists all GAP operations for which ``Operation(x, ...)`` is defined.
"""

from sage.misc.superseded import deprecation
from gappy.operations import *
deprecation(31297,
    'this module has been subsumed by the gappy package and will be removed')
