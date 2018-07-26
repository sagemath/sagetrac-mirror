# -*- coding: utf-8 -*-
r"""
The fundamental dual basis of PBT Hopf algebra.
"""
#*****************************************************************************
#       Copyright (C) 2012 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.cha.pbt import PlanarBinaryTreeFunctions
from sage.combinat.permutation import Permutation


class FundamentalDual(PlanarBinaryTreeFunctions.Bases.Base):
        '''
        This is the fundamental basis of ``PBT``::

            sage: P = PBT(QQ).P()
            The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field in the realization Fundamental

        .. MATH::

            \mathbb{F}_{\sigma} =
            \sum_{w\in \mathfrak{A}^*; std(w) = \sigma} w

        where `\mathfrak A` is an infinite and totally ordered alphabet such
        that `\mathfrak A^*` is free monoid. (``FQSym`` defined as a
        sub-algebra of the free algebra `\mathbb{K} \\langle \mathfrak{A}
        \\rangle`.)

        The product of `(\mathbb{F}_\sigma)_{\sigma\in\mathfrak{G}}` is
        described by

        .. MATH::

            \mathbb{F}_{\sigma} \\times \mathbb{F}_{\mu} =
            \sum_{\gamma \in \sigma \Cup \mu} \mathbb{F}_{\gamma}

        where `\sigma \Cup \mu` is the shifted shuffle of `\sigma` and `\mu`
        with. We denote by `std` the standardization map which  associate to
        a word `w` a permutation `\sigma` (see ``Permutations``).

        EXAMPLES::

            sage: F[1,2] * F[2,1]
            F[1, 2, 4, 3] + F[1, 4, 2, 3] + F[1, 4, 3, 2] + F[4, 1, 2, 3] + F[4, 1, 3, 2] + F[4, 3, 1, 2]

        And the coproduct is described by

        .. MATH::

            \Delta(\mathbb{F}_{\sigma}) =
            \sum_{i\in [n]} \mathbb{F}_{\sigma_{\mid i[}} \otimes 
            \mathbb{F}_{\sigma_{\mid [i}}

        where `[n] := [1,..,n]`, `i[` is the sub-interval of `[n]` defined by
        `[1,..,i-1]` and `[i` the sub-interval defined by `[i,..,n]`.

        EXAMPLES::

            sage: F[3,1,2,4].coproduct()
            F[] # F[3, 1, 2, 4] + F[1] # F[1, 2, 3] + F[2, 1] # F[1, 2] + F[3, 1, 2] # F[1] + F[3, 1, 2, 4] # F[]

        (See [MalReut]_ and [NCSF-VI]_.)

        EXAMPLES::

            sage: F().antipode()
            0
            sage: F[[]].antipode()
            F[]
            sage: F[1].antipode()
            -F[1]
            sage: F[1, 2].antipode()
            F[2, 1]
            sage: F[2, 1].antipode()
            F[1, 2]

        TESTS::

            sage: F = FQSym(QQ).F(); F
            Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field in the Fundamental basis
            sage: TestSuite(F).run()
        '''
        _prefix = "Q"

        def build_morphisms(self):
            self.morph_G_to_Q_FQSym()
            self.morph_Q_P()

        def morph_G_to_Q_FQSym(self):
            '''
            TESTS::


            '''
            from sage.combinat.cha.all import FQSym
            G = FQSym(self.base()).G()
            G.module_morphism(
                on_basis=lambda sigma: self(self.BST(sigma)),
                codomain=self
            ).register_as_coercion()

        def morph_Q_P(self):
            '''
            TESTS::
                sage: Q = PBT(QQ).Q(); P = PBT(QQ).P()
                sage: Q(P[2,3,1])
                Q[1, 3, 2]
                sage: Q(P[3,1,2])
                Q[2, 3, 1] + Q[1, 3, 2]
                sage: Q(P[4,2,1,3])
                Q[3, 2, 4, 1] + Q[1, 3, 4, 2] + Q[2, 1, 4, 3]
                sage: Q(P[5,4,2,1,3])
                Q[4, 3, 5, 2, 1] + Q[2, 4, 5, 3, 1] + Q[3, 2, 5, 4, 1] + Q[1, 4, 5, 3, 2] + Q[1, 3, 5, 4, 2] + Q[2, 1, 5, 4, 3]
                sage: Q(P[5,3,4,1,2])
                Q[2, 4, 5, 3, 1] + Q[3, 2, 5, 4, 1] + Q[2, 3, 5, 4, 1] + Q[1, 4, 5, 3, 2] + 2*Q[1, 3, 5, 4, 2] + Q[2, 1, 5, 4, 3] + Q[1, 2, 5, 4, 3]
            '''
            P = self.realization_of().P()
            P.module_morphism(
                on_basis=lambda bt: self.sum_of_monomials(map(
                    lambda sigma: self.BST(sigma.inverse()),
                    self.realization_of().sylvester_class(bt)
                )), codomain=self,
                #### FIXME:: must be invertible ####
                ## use a good comparaison method
                #triangular = ...,
                #cmp = ...
            ).register_as_coercion()

        def dual_basis(self):
            return self.realization_of().P()

        def product_on_basis(self, c1, c2):
            r'''
            The product of the F basis : concatenation of tree

            EXAMPLES::

                sage: Q = PBT(QQ).Q()
                sage: Q[2,1]*Q[3,1,2]
                Q[1, 5, 4, 3, 2] + Q[1, 3, 5, 4, 2] + Q[1, 4, 3, 5, 2] + Q[2, 1, 5, 4, 3] + Q[2, 1, 4, 5, 3] + Q[3, 2, 1, 5, 4] + Q[1, 2, 5, 4, 3] + Q[1, 2, 4, 5, 3] + Q[1, 3, 2, 5, 4] + Q[2, 1, 3, 5, 4]

            .. FIXME::

                improve the product : don't use permutation, compute directly
                on tree
                ... may be false ...
            '''
            from sage.combinat.words.word import Word
            c1 = Word(self.canonical_permutation(c1).inverse())
            c2 = Word(self.canonical_permutation(c2).inverse())
            ens = c1.shifted_shuffle(c2)
            return self.sum_of_monomials(
                [self.BST(Permutation(s).inverse()) for s in ens])

        def coproduct_on_basis(self, p):
            r"""
            The coproduct of the Q basis : de-shuffle of tree

            TESTS::

                sage: Q = PBT(QQ).Q()
                sage: Q[[]].coproduct()
                Q[] # Q[]
                sage: Q[1].coproduct()
                Q[] # Q[1] + Q[1] # Q[]
                sage: Q[1,2].coproduct()
                Q[] # Q[1, 2] + Q[1] # Q[1] + Q[1, 2] # Q[]
                sage: Q[1,2,3].coproduct()
                Q[] # Q[1, 2, 3] + Q[1] # Q[1, 2] + Q[1, 2] # Q[1] + Q[1, 2, 3] # Q[]

            EXAMPLES::

                sage: Q = PBT(QQ).Q()
                sage: Q[6,4,5,2,1,3].coproduct()
                Q[] # Q[2, 1, 4, 6, 5, 3] + Q[1] # Q[1, 3, 5, 4, 2] + Q[2, 1] # Q[2, 4, 3, 1] + Q[2, 1, 3] # Q[1, 3, 2] + Q[2, 1, 4, 3] # Q[2, 1] + Q[2, 1, 4, 5, 3] # Q[1] + Q[2, 1, 4, 6, 5, 3] # Q[]
            """
            def restrict(bt):
                if bt.is_empty():
                    yield (bt, bt)
                    return

                for bm, bM in restrict(bt[1]):
                    yield (self.indices()([bt[0], bm]), bM)
                for bm, bM in restrict(bt[0]):
                    yield (bm, self.indices()([bM, bt[1]]))
            return self.tensor_square().sum_of_monomials(restrict(p))
