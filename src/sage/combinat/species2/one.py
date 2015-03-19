# -*- coding: utf-8 -*-
"""
The species `1`

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
from sage.combinat.structures import Structure
from sage.misc.lazy_attribute import lazy_class_attribute, lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation


class One(UniqueRepresentation, Structure):

    def __init__(self, parent):
        Structure.__init__(self, parent)

    def _repr_(self):
        return r"·"

    @lazy_class_attribute
    def _auto_parent_(self):
        return OneSpecies()


class OneSpecies(SpeciesDesign):
    """
    The species `1`, characteristic of the *empty set*, defined by

    MATH::

        1[U] := \begin{dcases*}
            \{U\} & if `U = \emptyset`,\\
            \emptyset & otherwise,
        \end{dcases*}

    for any finite set `U`.

    It is the neutral element for the multiplication.
    (section 1.1 and 1.3, _[BBL])

    EXAMPLES::

        sage: from sage.combinat.species2.one import OneSpecies
        sage: O = OneSpecies()
        sage: for i in range(5): O.structures(Set(range(1,i+1))).list()
        [·]
        []
        []
        []
        []

    TESTS::

        sage: from sage.combinat.species2.one import OneSpecies
        sage: O = OneSpecies()
        sage: TestSuite(O).run()

    """

    def grading(self, I):
        assert(isinstance(I, One))
        return 0

    def _repr_(self):
        return "`1`"

    def transport(self, sigma):
        return lambda one: self._element_constructor_(sigma(one._label_))

    def some_elements(self):
        return iter([self.first()])

    def __iter__(self):
        yield self.graded_component(0).first()

    class Structures(SpeciesDesign.Structures):

        def cardinality(self):
            return 0 if self.grading() != 1 else 1

        def __iter__(self):
            if self.grading() == 0:
                yield self._element_constructor_()
            return

    Element = One