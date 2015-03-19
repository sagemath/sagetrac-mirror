# -*- coding: utf-8 -*-
"""
The species `X`, characteristic of *singletons*.

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
from sage.misc.lazy_attribute import lazy_attribute

class Singleton(Structure):

    @lazy_attribute
    def _auto_parent_(self):
        return SingletonsSpecies()

    def __init__(self, parent, label):
        Structure.__init__(self, parent)
        self._label_ = label

    def _repr_(self):
        return repr(self._label_)

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return other._label_ == self._label_
        return False


class SingletonsSpecies(SpeciesDesign):
    """
    The species `X`, characteristic of *singletons*, defined by

    MATH::

        X[U] = \begin{dcases*}
            \{U\} & if `|U| = 1`,\\
            \emptyset & otherwise,
        \end{dcases*}

    for any finite set `U`.

    (section 1.1, _[BLL])

    EXAMPLES::

        sage: from sage.combinat.species2.singletons import SingletonsSpecies
        sage: X = SingletonsSpecies()
        sage: for i in range(5): X.structures(Set(range(1,i+1))).list()
        []
        [1]
        []
        []
        []

    TESTS::

        sage: from sage.combinat.species2.singletons import SingletonsSpecies
        sage: X = SingletonsSpecies()
        sage: TestSuite(X).run()

    """

    def _repr_(self):
        return "`X`"
    
    def transport(self, sigma):

        def Fsigma(x):
            return self._element_constructor_(sigma.codomain())

        return Fsigma

    def some_elements(self):
        return iter([self.first()])

    def grading(self, U):
        assert(isinstance(U, Singleton))
        return 1

    def __iter__(self):
        yield self.graded_component(1).first()

    class Structures(SpeciesDesign.Structures):

        def cardinality(self):
            return 1 if self.grading() == 1 else 0

        def __iter__(self):
            if self.grading() == 1:
                yield self._element_constructor_(iter(self.finite_set()).next())

    Element = Singleton