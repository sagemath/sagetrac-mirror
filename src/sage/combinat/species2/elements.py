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
from sage.combinat.structures import Structure
from sage.misc.lazy_attribute import lazy_class_attribute


class Elem(Structure):

    @lazy_class_attribute
    def _auto_parent_(self):
        return ElementsSpecies()

    def __init__(self, parent, label, grading=None):
        Structure.__init__(self, parent)
        self._label_ = label
        self._grading_ = grading

    def _repr_(self):
        return repr(self._label_)

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return other._label_ == self._label_
        return False

    def __ne__(self, other):
        if isinstance(other, type(self)):
            return other._label_ != self._label_
        return True


class ElementsSpecies(SpeciesDesign):
    """
    The species `\varepsilon`, of *elements*, defined by

    MATH::

        E[U] = U

    for any finite set `U`.

    (section 1.1, _[BLL])

    EXAMPLES::

        sage: from sage.combinat.species2.elements import ElementsSpecies
        sage: e = ElementsSpecies()
        sage: for i in range(5): e.graded_component(i).list()
        []
        [1]
        [1, 2]
        [1, 2, 3]
        [1, 2, 3, 4]

    TESTS::

        sage: from sage.combinat.species2.elements import ElementsSpecies
        sage: e = ElementsSpecies()
        sage: TestSuite(e).run()

    """

    def _repr_(self):
        return "`e`"
    
    def transport(self, sigma):

        def Fsigma(u):
            #assert(u in sigma.domain())
            return self._element_constructor_(sigma(u._label_))

        return Fsigma

    def grading(self, u):
        assert(isinstance(u, self.element_class))
        return u._grading_

    class Structures(SpeciesDesign.Structures):

        def cardinality(self):
            return self.grading()

        def __iter__(self):
            for u in self.finite_set():
                yield self._element_constructor_(u, grading=self.grading())

    Element = Elem