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
from .base import CombinatorialPolyhedron_class

class CombinatorialPolyhedra(UniqueRepresentation, Parent):

    r"""
    Parent class for combinatorial polyhedra
    """

    Element = CombinatorialPolyhedron_class

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

    def _Hom_(self, other, category=None):
        """
        EXAMPLES::

        sage: C = polytopes.cube().combinatorial_polyhedron()
        sage: Hom(C.parent(), C.parent()) # indirect doctest
        Set of Morphisms
         from Combinatorial polyhedra
         to Combinatorial polyhedra
         in Category of combinatorial polyhedral sets
        """
        from sage.categories.combinatorial_polyhedral_sets import CombinatorialPolyhedralSets
        if category is not None and not category.is_subcategory(CombinatorialPolyhedralSets()):
            raise TypeError(f"{category} is not a subcategory of CombinatorialPolyhedralSets()")
        from .homset import CombinatorialPolyhedraHomset
        return CombinatorialPolyhedraHomset(self, other)

    def hom(self, Vrep_dict, codomain=None, check=True, category=None):
        """
        Create a morphism from ``self`` to another set of combinatorial polyhedra.

        EXAMPLES::

            sage: C = polytopes.cube().combinatorial_polyhedron()
            sage: C.parent().hom({1: 2, 2: 2}, C.parent())
            Combinatorial polyhedral set endomorphism of Combinatorial polyhedra

        """
        if codomain is None:
            codomain = self
        homset = self.Hom(codomain)
        return homset(Vrep_dict, check=check)
