# -*- coding: utf-8 -*-
"""
The species `E`, of *sets*.

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
from sage.sets.set import Set_object_enumerated, Set
from sage.structure.element_wrapper import ElementWrapper

class SetStructure(ElementWrapper):
    # FIXME it could be convenient to use `Set(U)` with Set the current sage class (not this one)

    wrapped_class = Set_object_enumerated

    def cardinality(self):
        return self.value.cardinality()

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self.value == other.value
        return False

    def __lt__(self, other):
        if isinstance(other, type(self)):
            return self.value < other.value
        return other >= self.value


class SetsSpecies(SpeciesDesign):
    """
    The species `E`, of *sets*, defined by

    MATH::

        E[U] = \{U\}

    for any finite set `U`.

    (section 1.1, _[BLL])

    EXAMPLES::

        sage: from sage.combinat.species2.sets import SetsSpecies
        sage: E = SetsSpecies()
        sage: for i in range(5): E.structures(Set(range(1,i+1))).list()
        [{}]
        [{1}]
        [{1, 2}]
        [{1, 2, 3}]
        [{1, 2, 3, 4}]

    TESTS::

        sage: from sage.combinat.species2.sets import SetsSpecies
        sage: E = SetsSpecies()
        sage: TestSuite(E).run()

    """

    def _repr_(self):
        return "`E`"
    
    def transport(self, sigma):

        def Fsigma(U):
            #assert(U == sigma.domain())
            #print U, "-->", Set(map(sigma, U.value))
            return self.element_class(self, Set(map(sigma, U.value))) # ??

        return Fsigma

    def grading(self, U):
        return U.cardinality()

    class Structures(SpeciesDesign.Structures):

        def cardinality(self):
            return 1

        def __iter__(self):
            F = self.ambient()
            yield F.element_class(F, self.finite_set())

    Element = SetStructure