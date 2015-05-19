# -*- coding: utf-8 -*-
"""
Recursive Species

References
----------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux,
  1998, Cambridge University Press

"""
# *******************************************************************************
#       Copyright (C) 2015 Jean-Baptiste Priez <jbp@kerios.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
# *******************************************************************************
from sage.categories.species import Species
from sage.combinat.species2 import SpeciesDesign
from sage.combinat.species2.cycle_index_series.operations.recursive_cis import RecursiveCIS
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.set import Set
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent


class RecursiveSpecies(Parent):
    """
    Species defined recursively

    EXAMPLES::

        sage: Sp = Species()
        sage: B = Sp.recursive_species()
        sage: X = Sp.singletons()
        sage: O = Sp.one()
        sage: B.define(O + B*X*B) # species of binary trees
        sage: list(B.type_generating_series().coefficients(10))
        [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796]

        sage: R = Sp.recursive_species()
        sage: E = Sp.sets()
        sage: R.define(X*E.composite(R)) # species of rooted trees
        sage: list(R.generating_series().coefficients(6))
        [0, 1, 2, 9, 64, 625, 7776]

    TESTS::

        sage: Sp = Species()
        sage: E = Sp.sets()
        sage: X = Sp.singletons()
        sage:
        sage: T = Sp.recursive_species()
        sage: T.define(E.functorial_composite(X*T))
        sage: T._valuation_()
        0
        sage: T = Sp.recursive_species()
        sage: T.define(X*E.functorial_composite(T))
        sage: T._valuation_()
        1

    """

    def __init__(self, name="F"):
        Parent.__init__(self, category=Species())
        self._def_ = None
        self._active_ = False
        self._name_ = name
        self._cis_ = RecursiveCIS(name="Z"+self._name_)
        self._structures = {}

    def define(self, definition):
        if self._def_:
            raise ValueError("This species is already defined: %s" % repr(self))
        self._def_ = definition
        self._cis_.define(self._def_.cycle_index_series())

    def _repr_(self):
        if self._active_:
            return self._name_
        self._active_ = True
        s = self._name_ + " := " + repr(self._def_)
        self._active_ = False

        return s

    def grading(self, s):
        pass

    def structures(self, U):
        if U in self._structures:
            return self._structures[U]

        FU = self.Structures(self, U)
        self._structures[U] = FU
        FU._ambient_ = self
        return FU

    def graded_component(self, k):
        return self.structures(Set(range(1, k+1)))

    def cycle_index_series(self):
        return self._cis_

    def transport(self, sigma):
        def ssigma(f):
            return self._element_constructor_(self._def_.transport(sigma)(f.value))
        return ssigma

    def _element_constructor_(self, *args, **options):
        return self.element_class(self, *args, **options)

    Element = ElementWrapper

    class Structures(SpeciesDesign.Structures):

        def __init__(self, ambient, U):
            SpeciesDesign.Structures.__init__(self, ambient, U)
            self._active_ = False

        @lazy_attribute
        def _element_constructor_(self):
            return self.ambient()._element_constructor_

        def __iter__(self):
            """
            TESTS::

                sage: Sp = Species()
                sage: B = Sp.recursive_species(name="B")
                sage: X = Sp.singletons()
                sage: O = Sp.one()
                sage: B.define(O + B*X*B)
                sage: B.structures(Set([1,2])).list()
                [[·, 1, [·, 2, ·]],
                 [·, 2, [·, 1, ·]],
                 [[·, 1, ·], 2, ·],
                 [[·, 2, ·], 1, ·]]

            """
            if self.ambient()._valuation_() > self.finite_set().cardinality():
                return

            for s in self.ambient()._def_.structures(self.finite_set()):
                yield self._element_constructor_(s)