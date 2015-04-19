# -*- coding: utf-8 -*-
"""
Derivative of species.

References
----------

 _[BBL] Combinatorial species and tree-like structures,
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
from sage.combinat.species2.operations.composite import Composite
from sage.combinat.species2.operations.product import Prod
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.sets.set import Set
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation


class Ast(UniqueRepresentation, SageObject):

    def _repr_(self):
        return "*"


class Derivative(SpeciesDesign):
    """
    The *derivative* of species

    The *derivative* of species `F`, noted `F'` is defined as follows:
    An `F'`-structure on `U` is an `F`-structure on `U^+ = U \cup \{\ast\}` where `\ast` is an element chosen outside of
    `U`. In other words, one sets `F'[U] = F[U^+]`.
    The transport along a bijection `\sigma : U \to V` is carried out by setting, by setting

    MATH::

        F'[\sigma](s) = F[\sigma^+](s)

    where `\sigma^+ : U \sqcup \{\ast\} \to V \sqcup \{\ast\}` is the canonical extension of `\sigma` obtained by
    setting `\sigma^+(u) = \sigma(u)` if `u \in U` and `\sigma^+(\ast) = \ast`.

    Properties:

     - Additivity: `(F + G)' = F' + G'`,
     - Product rule: `(F \cdot G)' = F' \cdot G + F \cdot G'`,
     - Chain rule: `(F \circ G)' = (F' \circ G) \cdot G'`,
     - Neutral:  `E' = E`.

    (section 1.4, _[BBL])
    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, F):
        # neutral
        if F == Species().sets():
            return F
        # product rule
        if isinstance(F, Prod):
            return Add(*[Prod(*(F._species_[:i] + [Derivative(F[i])] + F._species_[i+1:]))
                         for i in enumerate(F._species_)])
        # additivity
        if isinstance(F, Add):
            return Add(*[Derivative(G) for G in F._species_])
        # chain rule
        if isinstance(F, Composite):
            return Prod(Composite(Derivative(F._F_), F._G_), Derivative(F._G_))
        # otherwise
        return super(Derivative, cls).__classcall__(cls, F)

    def __init__(self, F):
        SpeciesDesign.__init__(self)
        self._F_ = F

    def _repr_(self):
        return repr(self._F_) + "'"

    def transport(self, sigma):
        """
        TESTS::

            sage: SP = SetPartitions()
            sage: SPd = SP.derivative()
            sage: SPd.isomorphism_types(2).list() # indirect doctest
            [[{{1, 2}, {*}}],
             [{{1}, {2, *}}, {{1, *}, {2}}],
             [{{1}, {2}, {*}}],
             [{{1, 2, *}}]]
        :param sigma:
        :return:
        """
        def nsigma(u):
            return u if u == Ast() else sigma(u)

        def Fpsigma(u):
            s = self._F_.transport(nsigma)(u)
            s._set_parent(self)
            return s
        return Fpsigma

    def grading(self, s):
        return self._F_.grading(s)

    def cycle_index_series(self):
        """
        The derivative of cycle index series

        MATH::

            Z_{F'}(p_1, p_2, \cdots) = \left(\frac{\partial}{\partial p_1} Z_F\right)(p_1, p_2, \cdots)

        (section 1.4, _[BBL])
        """
        return self._F_.cycle_index_series().derivative()

    class Structures(SpeciesDesign.Structures):

        def cardinality(self):
            pass

        def __iter__(self):
            """
            TESTS::

                sage: SP = SetPartitions()
                sage: SPd = SP.derivative()
                sage: SPd.graded_component(3).list()
                [{{1, 2, 3, *}},
                 {{1}, {2, 3, *}},
                 {{1, 3, *}, {2}},
                 {{1, 2, *}, {3}},
                 {{1, 2, 3}, {*}},
                 {{1, 2}, {3, *}},
                 {{1, 3}, {2, *}},
                 {{1, *}, {2, 3}},
                 {{1}, {2}, {3, *}},
                 {{1}, {2, *}, {3}},
                 {{1}, {2, 3}, {*}},
                 {{1, *}, {2}, {3}},
                 {{1, 3}, {2}, {*}},
                 {{1, 2}, {3}, {*}},
                 {{1}, {2}, {3}, {*}}]

            """
            if self.ambient()._valuation_() > self.finite_set().cardinality(): return
            for s in self.ambient()._F_.structures(Set([Ast()]).union(self.finite_set())):
                s._set_parent(self.ambient())
                yield s
