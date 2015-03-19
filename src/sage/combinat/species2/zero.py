# -*- coding: utf-8 -*-
"""
The species `0`

References
----------

 _[BBL] Combinatorial species and tree-like structures,
 François Bergeron, Gilbert Labelle and Pierre Leroux,
 1998, Cambridge University Press

"""
#*******************************************************************************
#       Copyright (C) 2015 Jean-Baptiste Priez <jbp@kerios.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*******************************************************************************
from sage.combinat.species2 import SpeciesDesign

class ZeroSpecies(SpeciesDesign):
    """
    The species `0` defined by

    MATH::

        0[U] := \emptyset

    for any finite set `U`.

    It is the neutral element for addition:

    MATH::

        F + 0 = 0 + F = F\,,

    and it is the absorbing element for the multiplication:

    MATH::

        F \cdot 0 = 0 \cdot F = 0\,,

    for all species `F`.
    (section 1.3, _[BBL])

    EXAMPLES::

        sage: from sage.combinat.species2.zero import ZeroSpecies
        sage: Z = ZeroSpecies()
        sage: for i in range(5): Z.structures(Set(range(1,i+1))).list()
        []
        []
        []
        []
        []

    TESTS::

        sage: from sage.combinat.species2.zero import ZeroSpecies
        sage: Z0 = ZeroSpecies()
        sage: TestSuite(Z0).run(skip=["_test_an_element", "_test_elements"])

    """
    # TODO there is a correct way to fix the TestSuite??
    def grading(self, I):
        raise ValueError("We never should be here, there is no structure associated to this species")

    def _repr_(self):
        return "`0`"

    def transport(self, sigma):
        def trsprt(struct):
            raise ValueError("Should never be here")
        return trsprt

    def some_elements(self):
        return iter([])

    def __iter__(self):
        return iter([])

    class Structures(SpeciesDesign.Structures):

        def cardinality(self):
            return 0

        def __iter__(self):
            return iter([])