# -*- coding: utf-8 -*-
r"""
The fundamental basis of WQSym Hopf algebra.
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


class Fundamental(WordQuasiSymmetricFunctions.Bases.Base):
        '''

        '''
        _prefix = "S"

        # def build_morphisms(self):
        #     '''
        #     '''
        #     self.F_of_FQSym_to_S()

        # def F_of_FQSym_to_S(self):
        #     from sage.combinat.cha.fqsym import FreeQuasiSymmetricFunctions
        #     from sage.combinat.permutation import to_standard
        #     F = FreeQuasiSymmetricFunctions(self.base()).F()
        #     self.module_morphism(
        #         on_basis=lambda pw: F(to_standard(pw)),
        #         codomain=F
        #     ).register_as_coercion()

            ##################################################################
            ##################################################################
            ##################################################################
                
            
        def dual_basis(self):
            return self.realization_of().M()

        def product_on_basis(self, sigma, mu):
            '''
            TESTS::

                sage: S = WQSym(QQ).S()
                sage: sage: S[1,1]* S[1,2]
                S[1, 1, 2, 3] + S[1, 2, 1, 3] + S[1, 2, 3, 1] + S[2, 1, 1, 3] + S[2, 1, 3, 1] + S[2, 3, 1, 1]
                sage: S[[]] * S[1]
                S[1]
                sage: S[1] * S[[]]
                S[1]
            '''
            from sage.combinat.words.word import Word
            Keys = self.basis().keys()
            return self.sum_of_monomials(map(
                lambda w: Keys(w),
                Word(sigma).shuffle(Word(
                    [l + max(list(sigma) + [0]) for l in mu]))
            ))

        def coproduct_on_basis(self, sigma):
            '''
            TESTS::

                sage: S = WQSym(QQ).S()
                sage: S[1].coproduct()
                S[] # S[1] + S[1] # S[]
                sage: S[1,1].coproduct()
                S[] # S[1, 1] + S[1] # S[1] + S[1, 1] # S[]
                sage: S[2,1,2,1,3].coproduct()
                S[] # S[2, 1, 2, 1, 3] + S[1] # S[1, 2, 1, 3] + S[2, 1] # S[2, 1, 3] + S[2, 1, 2] # S[1, 2] + S[2, 1, 2, 1] # S[1] + S[2, 1, 2, 1, 3] # S[]

            '''
            from sage.combinat.packed_words import to_pack
            return self.tensor_square().sum_of_monomials(
                (to_pack(sigma[:i]), to_pack(sigma[i:]))
                    for i in range(len(sigma) + 1)
            )

        def internal_product_on_basis(self,sigma,mu):
                '''
                '''
                from sage.combinat.packed_words import to_pack

                if len(sigma)!=len(mu):
                        raise IndexError
                n=len(mu)
                print sigma, mu
                return self.monomial(to_pack(
                        [str(sigma[i])+str(mu[i]) for i in range(n)]))
                        
                
