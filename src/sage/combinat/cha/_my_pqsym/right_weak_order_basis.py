# -*- coding: utf-8 -*-
r"""
The Left Weak Order basis of PQSym Hopf algebra.
"""


from my_pqsym import ParkingQuasiSymmetricFunctions


class RightWeakOrder(ParkingQuasiSymmetricFunctions.Bases.Base):
        '''
        
        '''
        _prefix = "R"
        
        
        def dual_basis(self):
                return self.realization_of().L()
                    
        def product_on_basis(self, pk1, pk2):
            '''
            TESTS::

                sage: R = PQSym(QQ).R()
                sage: R[1,1]*R[1]
                R[1, 1, 2] + R[1, 2, 1] + R[2, 1, 1]
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
