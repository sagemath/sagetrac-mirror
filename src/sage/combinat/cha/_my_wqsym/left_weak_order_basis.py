# -*- coding: utf-8 -*-
r"""
La base des G du papier de Vargas
"""



from sage.combinat.cha.wqsym import WordQuasiSymmetricFunctions


class LeftWeakOrder(WordQuasiSymmetricFunctions.Bases.Base):
    '''

    '''
    _prefix = "L"
    
    # def build_morphisms(self):
    #     from sage.combinat.permutation import to_standard
    #     F = FQSym(self.base()).F()
    #     self.module_morphism(
    #         on_basis=lambda pw: F(to_standard(pw)),
    #         codomain=F
    #     ).register_as_coercion()
    
    # def dual_basis(self):
    #     return self.realization_of().M()

    def product_on_basis(self, sigma, mu):
        '''
        TESTS::
        
        sage: L = WQSym(QQ).L()
        sage: sage: L[1,2,2] * L[1]
        L[1, 2, 2, 3] + L[1, 3, 3, 2] + L[2, 3, 3, 1]
        sage: L[[]] * L[1]
        L[1]
        sage: L[1] * L[[]]
        L[1]
        '''
        if len(sigma)==0 :
            return self.monomial(mu)
        if len(mu)==0 :
            return self.monomial(sigma)
        from sage.combinat.combination import Combinations
        Keys = self.basis().keys()
        m1=max(sigma)
        m2=max(mu)
        
        r1=list(Combinations(range(1,m1+m2+1),m1))
        r2=list(Combinations(range(1,m1+m2+1),m2))
        r2.reverse()
        res=[]
        
        for i in range(len(r1)):
            res_i=[r1[i][k-1] for k in sigma]+[r2[i][k-1] for k in mu]
            res.append(res_i)
        return self.sum_of_monomials(map(lambda w: Keys(w), res))


    def coproduct_on_basis(self, sigma):
        '''
        TESTS::
        
        sage: L = WQSym(QQ).L()
        sage: L[1].coproduct()
        L[] # L[1] + L[1] # L[]
        sage: L[1,1].coproduct()
        L[] # L[1, 1] + L[1] # L[1] + L[1, 1] # L[]
        sage: L[2,1,2,1,3].coproduct()
        L[] # L[2, 1, 2, 1, 3] + L[1, 1] # L[1, 1, 2] + L[2, 1, 2, 1] # L[1] + L[2, 1, 2, 1, 3] # L[]
        
        '''
        from sage.combinat.packed_words import to_pack
        return self.tensor_square().sum_of_monomials(
            (to_pack([x for x in sigma if x<=i]),
             to_pack([x for x in sigma if x>i]))
            for i in range(max(list(sigma)+[0]) + 1)
        )

