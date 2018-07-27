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
from sage.combinat.set_partition_ordered import OrderedSetPartition


class FundamentalDual(WordQuasiSymmetricFunctions.Bases.Base):
        '''

        '''
        _prefix = "M"

        def build_morphisms(self):
#            self.M_to_Monomial_of_QSym()
#            self.G_of_FQSym_to_M()
            self.L_to_M()
#            self.GL_to_M()
#            self.SR_to_M()
#            self.R_to_M()
            
                
        # def G_of_FQSym_to_M(self):
        #     '''
        #     TESTS::

        #         sage: G = FQSym(QQ).G()
        #         sage: M = WQSym(QQ).M()
        #         sage: M(G[1,2])
        #         M[1, 1] + M[1, 2]
        #         sage: M(G[3,1,2])
        #         M[2, 1, 1] + M[3, 1, 2]
        #     '''
        #     from sage.combinat.cha.fqsym import FreeQuasiSymmetricFunctions
        #     G = FreeQuasiSymmetricFunctions(self.base()).G()
        #     G.module_morphism(
        #         on_basis=lambda sigma: self.sum_of_monomials(
        #             self.basis().keys().permutation_to_packed_words(sigma)
        #         ), codomain=self
        #     ).register_as_coercion()

        # def M_to_Monomial_of_QSym(self):
        #     '''

        #     TESTS::

        #         sage: M = WQSym(QQ).M()
        #         sage: Mon = QuasiSymmetricFunctions(QQ).M()
        #         sage: Mon(M[2,1,2])
        #         M[1, 2]
        #         sage: Mon(M[1]**2)
        #         2*M[1, 1] + M[2]
        #         sage: Mon(M[2,1,2]**2)
        #         4*M[1, 1, 2, 2] + 2*M[1, 1, 4] + 2*M[1, 2, 1, 2] + 2*M[1, 3, 2] + 2*M[2, 2, 2] + M[2, 4]
        #     '''
        #     from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
        #     Mo = QuasiSymmetricFunctions(self.base()).M()
        #     self.module_morphism(
        #         on_basis=lambda pw: Mo(pw.to_composition()),
        #         codomain=Mo
        #     ).register_as_coercion()


        def _calcul_coef(self,e):
                if len(e)==0:
                        return [1]
                i=min(e)
                res=[]
                while i<=max(e):
                        k=1
                        while i in e:
                                i+=1
                                k+=1
                        res+=[k]
                        i+=1
                return res

        def _func_L_to_M(self,sigma):
                from sage.combinat.subset import Subsets
                from sage.combinat.packed_words import to_pack
                from sage.combinat.packed_words import PackedWord
                from sage.functions.other import factorial
                from sage.rings.all import ZZ
                from sage.misc.misc_c import prod
                
                if len(sigma)==0:
                        return self.one()
                res = []
                m=max(sigma)
                for e in Subsets(m-1):
                        r0=self._calcul_coef(e)
                        r1=[]
                        for x in sigma:
                                k=1
                                while x-k in e:
                                        k+=1
                                r1+=[x-k+1]
                        res+=[(ZZ(1)/prod(factorial(k) for k in r0),to_pack(r1))]
                return sum(c*self.monomial(PackedWord(pw)) for (c,pw) in res)
        
        def L_to_M(self):
                morph = lambda X, Y, func, tri = None, comp = None: (
                        X._module_morphism(func, codomain=Y, triangular=tri, cmp=comp)
                        if comp is not None else
                        X._module_morphism(func, codomain=Y, triangular=tri))
                L = self.realization_of().L()
                M = self
                
                # L to M and back
                L_to_M = morph(L, M,
                               lambda sigma: self._func_L_to_M(sigma),
                               tri="upper")
                L_to_M.register_as_coercion()
                (~L_to_M).register_as_coercion()
                
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
            return self.sum_of_monomials(OrderedSetPartition(osp).to_packed_word()
                for osp in Shifted_QuasiShuffleProduct(
                        pw1.to_ordered_set_partition(),
                        pw2.to_ordered_set_partition()))

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
