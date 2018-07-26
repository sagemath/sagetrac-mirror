# -*- coding: utf-8 -*-
r"""
The monomial basis of PQSym Hopf algebra.
"""


from my_pqsym import ParkingQuasiSymmetricFunctions


class Monomial(ParkingQuasiSymmetricFunctions.Bases.Base):
        '''

        '''
        _prefix = "M"

        #TODO
        def build_morphisms(self):
            self.L_to_M()


        def _calcul_coef(self,sigma,groups,t):
                # TODO c'est pour l'instant ici que je trouve le plus simple de placer le 0 pour les cas où il faut (cf my_pk.py lignes 60 et autour)
                # si sigma n'est pas de la forme [grands, petits] alors il faut annuler le res pour ceux qui sautent des valeurs dans t
                # si sigma n'est pas de la forme [petits, trou, grands] alors il faut annuler le res pour ceux qui ne sautent pas de valeurs dans t
                #print sigma, groups, t
                from sage.functions.other import factorial
                from sage.rings.all import ZZ
                from sage.misc.misc_c import prod

                # ICI je fais du cas par cas pour comprendre!!!
                # if (sigma==[1,1,2] or sigma==[1,2,1]) and t==((1,1,3),):
                #         #print 0
                #         return 0
                # if (sigma==[3,1,1] or sigma==[1,3,1]) and t==((1,1,1,1),(1,1,1,2)):
                #         #print 0
                #         return 0
                # if (sigma==[1,2,2,3] or sigma==[2,1,1,3]) and \
                #    (t==((1,1,1,3),) or t[0][3]==4 or t[0][2]==3):
                #         #print 0
                #         return 0
                # if (sigma==[3,1,1,2]) and t[0][2]==3:
                #         #print 0
                #         return 0
                # if sigma==[2,1,1,3] and t==((1,1,2,3),):
                #         return 0
                
                
                res=1
                for g in groups:
                        res_g=[1]
                        mg=min(g)
                        Mg=max(g)
                        for x in range(mg,Mg):
                                k=0
                                mgk=min(groups[k])
                                Mgk=max(groups[k])
                                while Mgk < x :
                                        k+=1
                                        mgk=min(groups[k])
                                        Mgk=max(groups[k])
                                #print t,k,x
                                if t[k][x] == t[k][x+1]:
                                        res_g[-1]+=1
                                else:
                                        res_g+=[1]
                        #print res_g
                        res=res*(ZZ(1)/prod(factorial(k) for k in res_g))
                #print res
                return res

        
        def reeval(self,m,mg,Mg,d,i,strict_p1,res_g):
                res_g+=[tuple(d)]#[x] for x in range(mg,Mg+1))]
                for j in xrange(i+1,Mg+1):
#                        if strict_p1 <= i:
#                                d[j]=d[j-1]
#                                res_g=self.reeval(m,mg,Mg,d,j,strict_p1,res_g)
#                                if d[j-1] < m:
#                                        d[j]=d[j-1]+1
#                                        res_g=self.reeval(m,mg,Mg,d,j,strict_p1,res_g)
#                        else:
                        for k in xrange(d[j-1],m+1):
                                d[j]=k
                                res_g=self.reeval(m,mg,Mg,d,j,strict_p1,res_g)
                return res_g

                
        def _func_L_to_M(self,sigma):
                from sage.combinat.subset import Subsets
                from my_pk import to_park
                from sage.combinat.parking_functions import ParkingFunction
                from sage.combinat.parking_functions import is_a
                from sage.misc.misc_c import prod
                from sage.misc.misc import uniq
                from sage.categories.cartesian_product import cartesian_product
                
                
                if len(sigma)==0:
                        return self.one()
                m=max(sigma)
                n=len(sigma)
                d=dict()
                l=sorted(sigma)
                groups=[[]]

                j=l[0]
                for i in l:
                        if i-j>1:
                                groups+=[[]]
                        groups[-1]+=[i]
                        j=i

                res_g=[]    
                for g in groups:
                        mg=min(g)
                        Mg=max(g)
                        strict_p1_g = mg
                        while g.count(strict_p1_g)==1:
                                strict_p1_g+=1
                        d=[1]*(m+1)
                        #print "res_g avant", res_g
                        res_g+=[tuple(uniq(self.reeval(n,mg,Mg,d,mg-1,strict_p1_g,[])))]
                        #print "res_g apres", res_g

                res=[]
#                res_g=uniq(res_g)
                #print sigma
                #print groups
                #print res_g
                for t in cartesian_product(res_g):
                        new_pk=[]
                        for x in sigma:
                                k=0
                                mg=min(groups[k])
                                Mg=max(groups[k])
                                while Mg < x :
                                        k+=1
                                        mg=min(groups[k])
                                        Mg=max(groups[k])
                                #print t,k,x
                                new_pk+=[t[k][x]]
                        #print new_pk
                        if is_a(new_pk):
                                c=self._calcul_coef(sigma,groups,t)
                                res+=[(c,new_pk)]


                return sum(c*self.monomial(ParkingFunction(pk)) for (c,pk) in res \
                           if is_a(pk))
        
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

            
        def product_on_basis(self, pk1, pk2):
            '''
            TESTS::

                sage: M = PQSym(QQ).M()
                sage: M[1,1]*M[1]
                M[1, 1, 1] + M[1, 1, 2] + M[1, 1, 3] + M[2, 2, 1]
                sage: M[1,1,3]*M[1]
                M[1, 1, 3, 1] + M[1, 1, 3, 2] + M[1, 1, 3, 3] + M[1, 1, 3, 4] + \
                M[1, 1, 4, 1] + M[1, 1, 4, 2] + M[1, 1, 4, 3] + M[2, 2, 4, 1]
            '''
            from my_pk import to_park
            from sage.combinat.parking_functions import ParkingFunctions

            n=len(pk1)
            m=len(pk2)
            
            return self.sum_of_monomials(p for p in ParkingFunctions(n+m)
                                         if to_park(p[:n])==pk1 and to_park(p[n:])==pk2)

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
