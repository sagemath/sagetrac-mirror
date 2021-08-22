"""
Substitution Maps

This object wraps Pynac ``exmap`` objects. These encode substitutions
of symbolic expressions. The main use of this module is to hook into
Pynac's ``subs()`` methods and pass a wrapper for the substitution map
back to Python.
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Re-export
from .expression import SubstitutionMap as SubstitutionMap
from .expression import make_map as make_map
