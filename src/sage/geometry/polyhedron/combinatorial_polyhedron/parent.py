r"""
Parent class for combinatorial polyhedra
"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

class CombinatorialPolyhedra(UniqueRepresentation, Parent):

    r"""
    Parent class for combinatorial polyhedra
    """

    def __init__(self):
        """
        Construct the parent of class ``CombinatorialPolyhedra``.

        """
        from sage.categories.combinatorial_polyhedral_sets import CombinatorialPolyhedralSets
        Parent.__init__(self, category=CombinatorialPolyhedralSets())

    def _repr_(self):
        """
        Return a string representation.
        """
        return "Combinatorial polyhedra"
