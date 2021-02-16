"""
GAP element wrapper

This document describes the individual wrappers for various GAP elements. For
general information about GAP, you should read the :mod:`~sage.libs.gap.libgap`
module documentation.

This module and the classes in it are deprecated as it has been superseded by
the gappy package; see the documentation for :mod:`gappy.gapobj` for the
documentation on the equivalent classes in gappy.
"""

# ****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.superseded import deprecation

from gappy import gapobj


deprecation(31404,
    'this module has been subsumed by the gappy package and will be removed '
    'in a future version; all the GapElement classes are still available but '
    'as aliases for the equivalent gappy.gapobj.GapObj classes, and are no '
    'longer subclasses of Element')


def import_gapobjs_as_gapelements():
    """
    Creates aliases for the various `~gappy.gapobj.GapObj` subclasses under
    the names of the original ``GapElement`` subclasses to provide a bit of
    backwards compatibility for ``isinstance(..., GapElement)`` checks in
    existing code.

    EXAMPLES::

        sage: from gappy.gapobj import GapObj, GapList
        sage: from sage.libs.gap.element import GapElement, GapElement_List
        sage: GapObj is GapElement
        True
        sage: GapList is GapElement_List
        True
        sage: isinstance(libgap.eval('[]'), GapElement_List)
        True
    """

    for name, cls in vars(gapobj).items():
        if name[0] == '_':
            continue

        if not (isinstance(cls, type) and issubclass(cls, gapobj.GapObj)):
            continue

        elem_name = 'GapElement'
        sub = name[3:]  # strip 'Gap' prefix
        if sub and sub != 'Obj':
            elem_name += '_' + sub

        locals()[elem_name] = cls


import_gapobjs_as_gapelements()
