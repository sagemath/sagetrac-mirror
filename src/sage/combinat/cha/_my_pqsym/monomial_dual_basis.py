# -*- coding: utf-8 -*-
r"""
The monomial dual basis of PQSym Hopf algebra.
"""


from my_pqsym import ParkingQuasiSymmetricFunctions


class MonomialDual(ParkingQuasiSymmetricFunctions.Bases.Base):
        '''

        '''
        _prefix = "S"

        
        def dual_basis(self):
                return self.realization_of().M()

            
        def product_on_basis(self, pk1, pk2):
            '''
            TESTS::

                sage: S = PQSym(QQ).S()
                sage: S[1,1]*S[1]
                S[1, 1, 2] + S[1, 2, 1] + S[2, 1, 1]
            '''
            from sage.combinat.words.word import Word
            from sage.combinat.words.shuffle_product import \
                ShuffleProduct_w1w2
            if len(pk1)==0:
                return self.monomial(pk2)
            n=len(pk1)
            Keys = self.basis().keys()
            return self.sum_of_monomials(
                map(lambda w: Keys(list(w)), ShuffleProduct_w1w2(
                        Word(pk1),
                        Word([x+n for x in pk2])
            )))
        
        # def coproduct_on_basis(self, e):
        #     '''
        #     TESTS::

        #         sage: M = WQSym(QQ).M()
        #         sage: M[1,1].coproduct()
        #         M[] # M[1, 1] + M[1, 1] # M[]
        #         sage: M[[]].coproduct()
        #         M[] # M[]
        #         sage: M[2,1,2,1,4,3].coproduct()
        #         M[] # M[2, 1, 2, 1, 4, 3] + M[1, 1] # M[1, 1, 3, 2] + M[2, 1, 2, 1] # M[2, 1] + M[2, 1, 2, 1, 3] # M[1] + M[2, 1, 2, 1, 4, 3] # M[]
        #     '''
        #     def restriction(w, i):
        #         r''' Restriction of intervalle [min, i] and [i+1, max]

        #         '''
        #         from sage.combinat.packed_words import to_pack
        #         left = []
        #         right = []
        #         for l in w:
        #             if l >= i + 1:
        #                 right.append(l)
        #             else:
        #                 left.append(l)
        #         return (self.basis().keys()(left), to_pack(right))

        #     return self.tensor_square().sum_of_monomials(
        #         restriction(e, i)
        #         for i in set([-1] + list(e)))
