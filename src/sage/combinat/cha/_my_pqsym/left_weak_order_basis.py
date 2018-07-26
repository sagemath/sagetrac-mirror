# -*- coding: utf-8 -*-
r"""
The Left Weak Order basis of PQSym Hopf algebra.
"""


from my_pqsym import ParkingQuasiSymmetricFunctions


class LeftWeakOrder(ParkingQuasiSymmetricFunctions.Bases.Base):
        '''
        
        '''
        _prefix = "L"
        
        
        def dual_basis(self):
                return self.realization_of().R()
        
        def product_on_basis(self, sigma, mu):
                '''
                TESTS::
                
                sage: L = PQSym(QQ).L()
                sage: sage: L[1,1]* L[1,2]
                L[1, 1, 2, 3] + L[1, 2, 1, 3] + L[1, 2, 3, 1] + L[2, 1, 1, 3] + L[2, 1, 3, 1] + L[2, 3, 1, 1]
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
                M1=len(sigma)
                m1=M1
                while m1 in sigma:
                        m1=m1-1
                M2=1
                while M2 in mu:
                        M2=M2+1
                
                r1=list(Combinations(range(m1+1,M1+M2),M1-m1))
                r2=list(Combinations(range(m1+1,M1+M2),M2-1))
                r2.reverse()
                res=[]
                
                for i in range(len(r1)):
                        res_i=[]
                        for k in sigma:
                                if k<m1:
                                        res_i.append(k)
                                else:
                                        res_i.append(r1[i][k-m1-1])
                        for k in mu:
                                if k>M2:
                                        res_i.append(k+M1)
                                else:
                                        res_i.append(r2[i][k-1])
                        res.append(res_i)
                return self.sum_of_monomials(map(lambda w: Keys(w), res))
