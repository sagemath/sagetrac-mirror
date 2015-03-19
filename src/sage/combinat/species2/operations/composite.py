# -*- coding: utf-8 -*-
"""
(Partitional) Composite of species.

References
----------

 _[BBL] Combinatorial species and tree-like structures,
 François Bergeron, Gilbert Labelle and Pierre Leroux,
 1998, Cambridge University Press

AUTHOR:

- Jean-Baptiste Priez (2015)
"""
#*****************************************************************************
#       Copyright (C) 2015 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.species import Species
from sage.combinat.set_partition import SetPartitions
from sage.combinat.species2 import SpeciesDesign
from sage.combinat.species2.operations.add import Add
from sage.combinat.species2.operations.product import Prod
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.sets.set import Set
from sage.structure.list_clone import ClonableArray


class CompositeStructure(ClonableArray):

    def check(self):
        pass


class Composite(SpeciesDesign):
    """
    (Partitional) Composite of species.

    Let `F` and `G` be two species such that `G[\emptyset] = \emptyset`.
    The species `F \circ G`, also denoted `F(G)`, the (partitional) composite of `G` in `F`, is defined as follows:
    An `(F \circ G)`-structure on `U` is a triplet `s = (\pi, \varphi, \gamma)`, where

     - `\pi` is a partition of `U`,
     - `\varphi` is an `F`-structure on the set of classes of `\pi`,
     - `\gamma = (\gamma_p)_{p \in \pi}`, where for each class `p` of `\pi`, `\gamma_p` is a `G`-structure on `p`.

    In other words, for any finite set `U`, one has

    MATH::

        (F \circ G)[U] = \sum_{\pi \text{ partition of } U} F[\pi] \times \prod_{p \in \pi} G[p]\,,

    the (disjoint) sum being taken over the set of partitions `\pi` of `U` (*i.e.*, `\pi \in Par[U]`).

    The transport along a bijection `\sigma : U \to V` is carried out by setting, for any `(F \circ G)`-structure
    `s = (\pi, \varphi, (\gamma_p)_{p \in \pi})` on `U`,

    MATH::

        (F \circ G)[\sigma](s) = (\bar\pi, \bar\varphi, (\bar\gamma_{\bar{p}})_{\bar{p} \in \bar\pi})\,,

    where

     - `\bar\pi` is the partition of `V` obtained by the transport of `\pi` along `\sigma`,
     - for each `\bar{p} = \sigma(p) \in \bar\pi`, the structure `\bar\gamma_{\bar{p}}` is obtained from the structure
     `gamma_p` by `G`-transport along `\sigma_{|p}`,
     - the structure `\bar\varphi` is obtained from the transport `\varphi` by `F`-transport along the bijection
     `\bar\sigma` induced on `\pi` by `\sigma`.

    The partitional composition has the following properties:

     - associativity: `(F \circ G) \circ H = F \circ (G \circ H)`,
     - neutral element: `F \circ X = X \circ F = F`,
     - distributivity: `(F + G) \circ H = F \circ H + G \circ H`,
     - distributivity: `(F \cdot G) \circ H = (F \circ H) \cdot (G \circ H)`,
     - `F_0 = F \circ 0` and `F[\emptyset] = \emptyset \Longleftrightarrow F(0) = 0`.

    (section 1.4, _[BBL])
    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        F = args[0]
        G = args[1]
        ## Neutral element
        if F == Species().singleton():
            return G
        if G == Species().singleton():
            return F

        ## Associativity re-orientation F ○ (G ○ H) --> (F ○ G) ○ H
        if isinstance(G, Composite):
            G, H = G._F_, G._G_
            return Composite(Composite(F, G), H)

        ## Distrib +
        if isinstance(F, Add):
            return Add(*[Composite(H, G) for H in F._species_])

        ## Distrib ·
        if isinstance(F, Prod):
            return Prod(*[Composite(H, G) for H in F._species_])

        return super(Composite, cls).__classcall__(cls, *args, **opts)

    def __init__(self, F, G):
        assert(G.graded_component(0).cardinality() == 0), r"%s[∅] must be empty"%repr(G)
        SpeciesDesign.__init__(self)
        self._F_ = F
        self._G_ = G

    def _repr_(self):
        G = self._G_
        rG = ("(" + repr(G) + ")") if isinstance(G, (Add, Prod)) else repr(G)
        return "(" + repr(self._F_) + r"○" + rG + ")"

    def transport(self, sigma):

        def Fsigma(s):
            phi, gamma = s[0], s[1]
            ##############
            barsigma = lambda part: Set(map(sigma, part)) # ???
            barphi   = self._F_.transport(barsigma)(phi)
            bargamma = Set(self._G_.transport(sigma)(t) for t in gamma[:])
            return self._element_constructor_((barphi, bargamma))

        return Fsigma

    def _element_constructor_(self, *args, **options):
        return self.element_class(self, *args, **options)

    class Structures(SpeciesDesign.Structures):

        def __iter__(self):
            def rec_prod(pi):
                if len(pi) == 0:
                    yield ()
                    return
                p = pi[0]
                for s in G.structures(p):
                    for tup in rec_prod(pi[1:]):
                        yield (s,) + tup

            F = self.ambient()._F_
            G = self.ambient()._G_
            for pi in SetPartitions(self.finite_set()):
                for phi in F.structures(Set(pi[:])):
                    for tup in rec_prod(list(pi)):
                        yield self._element_constructor_((phi, Set(tup)))

    Element = CompositeStructure


