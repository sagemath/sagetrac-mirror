# -*- coding: utf-8 -*-
"""
Product of species.

References
----------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux,
  1998, Cambridge University Press

"""
# *****************************************************************************
#       Copyright (C) 2015 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
# *****************************************************************************
from sage.categories.species import Species
from sage.combinat.species2 import SpeciesDesign
from sage.combinat.species2.operations.add import Add
from sage.combinat.species2.operations.cartesian_product import CartesianProduct
from sage.misc.classcall_metaclass import ClasscallMetaclass


class FunctorialComposite(SpeciesDesign):
    """
    The *functorial composite* of species

    The species `F \Box G` (also denoted `F[G]`) is the *functorial composite* of `F` and `G`. It is defined as follows:
    An `(F \Box G)`-structures on `U` is an `F`-structure placed on the set `G[U]` of all the `G`-structures on `U`.

    In other words, for any finite set `U`,

    MATH::

        (F \Box G)[U] = F[G[U]]\,.

    The transport along a bijection `\sigma : U \to V` is carried out by setting

    MATH::

        (F \Box G)[\sigma] = F[G[\sigma]]\,.

    Properties:
    -----------

     - Associative,
     - Neutral element `E^\bullet`,
     - Right distributive by the cartesian product of species: `(F \times G) \Box H = (F \Box H) \times (G \Box H)`,
     - `F \times F = ((X + X^2) \cdot E) \Box F`.
     - Distributive on the right by the sum of species: `(F + G) \Bod H = F \Box H + G \Box H`.

    (section 2.2, _[BBL])

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, F, G):
        # neutral
        if F == Species().sets().pointing():
            return G
        if G == Species().sets().pointing():
            return F
        # associative
        if isinstance(G, FunctorialComposite):
            return FunctorialComposite(FunctorialComposite(F, G._F_), G._G_)
        # right distributive
        if isinstance(F, CartesianProduct):
            return CartesianProduct(*[FunctorialComposite(H, G) for H in F._species_])
        if isinstance(F, Add):
            return Add(*[FunctorialComposite(H, G) for H in F._species_])
        # otherwise
        return super(FunctorialComposite, cls).__classcall__(cls, F, G)

    def __init__(self, F, G):
        SpeciesDesign.__init__(self)
        self._F_, self._G_ = F, G

        self.Element = self._F_.Element

    def transport(self, sigma):
        def FBoxGsigma(u):
            s = self._F_.transport(self._G_.transport(sigma))(u)
            s._set_parent(self)
            return s
        return FBoxGsigma

    def grading(self, s):
        # FIXME
        F = self._F_
        s._set_parent(F)
        n = F.grading(s)
        s._set_parent(self)
        return n

    def cycle_index_series(self):
        ZF = self._F_.cycle_index_series()
        ZG = self._G_.cycle_index_series()
        return ZF.functorial_composite(ZG)

    class Structures(SpeciesDesign.Structures):

        def cardinality(self):
            pass

        def __iter__(self):
            F, G = self.ambient()._F_, self.ambient()._G_
            for s in F.structures(G.structures(self.finite_set())):
                s._set_parent(self.ambient())
                yield s
