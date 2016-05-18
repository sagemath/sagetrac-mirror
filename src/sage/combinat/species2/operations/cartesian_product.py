# -*- coding: utf-8 -*-
"""
Cartesian product of species

Reference
---------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux
  Cambridge University Press, 1998

"""
# *****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ******************************************************************************
from itertools import product
from sage.categories.cycle_index_series import CycleIndexSeries
from sage.categories.species import Species
from sage.combinat.species2 import SpeciesDesign
from sage.combinat.species2.operations.add import Add
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.list_clone import ClonableArray


class CartProdStructure(ClonableArray):

    def check(self):
        assert(all(s in self.parent()._species_[i] for i, s in enumerate(self)))

class CartesianProduct(SpeciesDesign):
    """
    The *cartesian product* of species

    The species `F \times G`, called *Cartesian product* of `F` and `G`, is defined as follows:
    An `(F \times G)`-structure on a finite set `U` is a pair `s = (f, g)`, where

     - `f` is an `F`-structure on `U`,
     - `g` is a `G`-structure on `U`.

    In other words, for all finite sets `U`, one has

    MATH::

        (F \times G)[U] = F[U] \times G[U]\,.

    The transport along a bijection `\sigma : U \to V` is carried out by setting

    MATH::

        (F \times G)[\sigma](s) = (F[\sigma](f), G[\sigma](g))\,,

    for any `(F \times G)`-structure `s = (f, g)` on `U`.

    Properties:
    -----------

     - Associativity: `(F \times G) \times H = F \times (G \times H)`,
     - Commutativity: `F \times G = G \times F`,
     - Neutral element: `F \times E = E \times F = F`,
     - Distributivity: `F \times (G + H) = F \times G + F \times H`,
     - `(F \times G)^\bullet = F^\bullet \times G = F \times G^\bullet`.

    (section 2.1, _[BBL)

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *species):

        list_of_species = []
        for i, F in enumerate(species):
            # neutral element
            if F == Species().sets():
                continue
            # associativity
            if isinstance(F, CartesianProduct):
                list_of_species.extend(F._species_)
            # distributivity
            if isinstance(F, Add):
                return Add(*[CartesianProduct(G, *(tuple(list_of_species) + species[i+1:]))
                           for G in F._species_])
            # otherwise
            list_of_species.append(F)

        if len(list_of_species) == 0:
            return Species().sets()
        elif len(list_of_species) == 1:
            return list_of_species[0]

        return super(CartesianProduct, cls).__classcall__(cls, *sorted(list_of_species))

    def __init__(self, *species):
        SpeciesDesign.__init__(self)
        self._species_ = species

    def _repr_(self):
        r = ""
        Flast = self._species_[0]
        acc = 1
        for F in self._species_[1:]:
            if F == Flast:
                acc += 1
            else:
                r += ("" if r == "" else "x") + repr(Flast) + ("^(x%d)" % acc if acc > 1 else "")
                acc = 1
            Flast = F
        r += ("" if r == "" else "x") + repr(Flast) + ("^(x%d)" % acc if acc > 1 else "")
        return r

    def transport(self, sigma):
        def Fsigma(u):
            # TODO
            pass
        return Fsigma

    def cycle_index_series(self):
        return reduce(lambda ZF, G: ZF.Hadamard_product(G.cycle_index_series()),
                      self._species_, CycleIndexSeries().sets())

    def grading(self, s):
        return self._species_[0].grading(s[0])

    def _element_constructor_(self, *args, **options):
        return self.element_class(self, *args, **options)

    class Structures(SpeciesDesign.Structures):

        def __iter__(self):
            """
            TESTS::

                sage: P = Permutations()
                sage: PP = P.cartesian_product(P)
                sage: PP.graded_component(2).list()
                [[[1, 2], [1, 2]], [[1, 2], [2, 1]], [[2, 1], [1, 2]], [[2, 1], [2, 1]]]

            """
            if self.ambient().valuation() > self.finite_set().cardinality(): return
            for s in product(*map(lambda F: F.structures(self.finite_set()), self.ambient()._species_)):
                yield self._element_constructor_(s)

    Element = CartProdStructure
