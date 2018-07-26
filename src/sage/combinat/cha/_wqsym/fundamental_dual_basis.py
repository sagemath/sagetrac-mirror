# -*- coding: utf-8 -*-
r"""
The fundamental dual basis of WQSym Hopf algebra.
"""
#*****************************************************************************
#       Copyright (C) 2012 Jean-Baptiste Priez <jbp@kerios.fr>,
#                          Rémi Maurice <maurice@univ-mlv.fr>.
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
from sage.combinat.cha.wqsym import WordQuasiSymmetricFunctions


class FundamentalDual(WordQuasiSymmetricFunctions.Bases.Base):
        '''

        '''
        _prefix = "M"

        def build_morphisms(self):
            self.M_to_Monomial_of_QSym()
            self.G_of_FQSym_to_M()

        def G_of_FQSym_to_M(self):
            '''
            TESTS::

                sage: G = FQSym(QQ).G()
                sage: M = WQSym(QQ).M()
                sage: M(G[1,2])
                M[1, 1] + M[1, 2]
                sage: M(G[3,1,2])
                M[2, 1, 1] + M[3, 1, 2]
            '''
            from sage.combinat.cha.fqsym import FreeQuasiSymmetricFunctions
            G = FreeQuasiSymmetricFunctions(self.base()).G()
            G.module_morphism(
                on_basis=lambda sigma: self.sum_of_monomials(
                    self.basis().keys().permutation_to_packed_words(sigma)
                ), codomain=self
            ).register_as_coercion()

        def M_to_Monomial_of_QSym(self):
            '''

            TESTS::

                sage: M = WQSym(QQ).M()
                sage: Mon = QuasiSymmetricFunctions(QQ).M()
                sage: Mon(M[2,1,2])
                M[1, 2]
                sage: Mon(M[1]**2)
                2*M[1, 1] + M[2]
                sage: Mon(M[2,1,2]**2)
                4*M[1, 1, 2, 2] + 2*M[1, 1, 4] + 2*M[1, 2, 1, 2] + 2*M[1, 3, 2] + 2*M[2, 2, 2] + M[2, 4]
            '''
            from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
            Mo = QuasiSymmetricFunctions(self.base()).M()
            self.module_morphism(
                on_basis=lambda pw: Mo(pw.to_composition()),
                codomain=Mo
            ).register_as_coercion()

        def dual_basis(self):
            return self.realization_of().S()

        def product_on_basis(self, pw1, pw2):
            '''
            TESTS::

                sage: M = WQSym(QQ).M()
                sage: M[1,1]*M[1]
                M[1, 1, 1] + M[1, 1, 2] + M[2, 2, 1]
            '''
            from quasi_shuffle_product import Shifted_QuasiShuffleProduct
            from sage.combinat.packed_words import \
                ordered_partition_sets_to_packed_word
            return self.sum_of_monomials(
                ordered_partition_sets_to_packed_word(osp)
                for osp in Shifted_QuasiShuffleProduct(
                        pw1.to_ordered_partition_sets(),
                        pw2.to_ordered_partition_sets()))

        def coproduct_on_basis(self, e):
            '''
            TESTS::

                sage: M = WQSym(QQ).M()
                sage: M[1,1].coproduct()
                M[] # M[1, 1] + M[1, 1] # M[]
                sage: M[[]].coproduct()
                M[] # M[]
                sage: M[2,1,2,1,4,3].coproduct()
                M[] # M[2, 1, 2, 1, 4, 3] + M[1, 1] # M[1, 1, 3, 2] + M[2, 1, 2, 1] # M[2, 1] + M[2, 1, 2, 1, 3] # M[1] + M[2, 1, 2, 1, 4, 3] # M[]
            '''
            def restriction(w, i):
                r''' Restriction of intervalle [min, i] and [i+1, max]

                '''
                from sage.combinat.packed_words import to_pack
                left = []
                right = []
                for l in w:
                    if l >= i + 1:
                        right.append(l)
                    else:
                        left.append(l)
                return (self.basis().keys()(left), to_pack(right))

            return self.tensor_square().sum_of_monomials(
                restriction(e, i)
                for i in set([-1] + list(e)))
