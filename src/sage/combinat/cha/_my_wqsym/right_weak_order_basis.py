# -*- coding: utf-8 -*-
r"""
La base des F du papier de Vargas
"""



from sage.combinat.cha.wqsym import WordQuasiSymmetricFunctions


class RightWeakOrder(WordQuasiSymmetricFunctions.Bases.Base):
        '''

        '''
        _prefix = "R"

        # def build_morphisms(self):
        #     from sage.combinat.permutation import to_standard
        #     F = FQSym(self.base()).F()
        #     self.module_morphism(
        #         on_basis=lambda pw: F(to_standard(pw)),
        #         codomain=F
        #     ).register_as_coercion()

#        def dual_basis(self):
#            return self.realization_of().M()

        def product_on_basis(self, sigma, mu):
            '''
            TESTS::

                sage: R = WQSym(QQ).R()
                sage: R[1,2,2]* R[1]
                R[1, 2, 2, 3] + R[1, 2, 3, 2] + R[1, 3, 2, 2] + R[3, 1, 2, 2]
            '''
            from sage.combinat.words.word import Word
            from sage.combinat.words.shuffle_product import \
                ShuffleProduct_w1w2
            if len(sigma)==0:
                return self.monomial(mu)
            m=max(sigma)
            Keys = self.basis().keys()
            return self.sum_of_monomials(
                map(lambda w: Keys(list(w)), ShuffleProduct_w1w2(
                        Word(sigma),
                        Word([x+m for x in mu])
            )))
            
        def coproduct_on_basis(self, sigma):
            '''
            TESTS::

                sage: R = WQSym(QQ).R()
                sage: R[1].coproduct()
                R[] # R[1] + R[1] # R[]
                sage: R[1,1].coproduct()
                R[] # R[1, 1] + R[1, 1] # R[]
                sage: R[2,1,2,1,3].coproduct()
                R[] # R[2, 1, 2, 1, 3] + R[2, 1, 2, 1] # R[1] + R[2, 1, 2, 1, 3] # R[]

            '''
            from sage.combinat.packed_words import to_pack
            return self.tensor_square().sum_of_monomials(
                (to_pack(sigma[:i]), to_pack(sigma[i:]))
                    for i in range(len(sigma) + 1)
                    if [sigma[:i].count(x)
                        for x in sigma[i:]].count(0)==len(sigma[i:])
            )


        def build_morphisms(self):
                self.S_to_R()
                

        def sub_composition(self,l):
                if len(l)==0:
                        return []
                if len(l)==1:
                        return [[i] for i in range(l[0])]
                res=self.sub_composition(l[1::])
                return [[i]+ll for ll in res for i in range(l[0])]

        def packedwords_of_alphabet(self,n,al):
                from sage.combinat.set_partition_ordered import OrderedSetPartitions
                from sage.combinat.packed_words import ordered_partition_sets_to_packed_word 
                
                l=[ordered_partition_sets_to_packed_word(e)
                   for e in OrderedSetPartitions(range(n),len(al))]
                return [[al[i-1] for i in ll] for ll in l]
    
        def untass_pw(self,sigma,ll):
                from sage.misc.misc_c import prod
                
                m=max(sigma)
                l=sigma.to_composition()
                r=[[i+1+sum(ll[:i])]*c for i,c in enumerate(l)]
                rr=[self.packedwords_of_alphabet(len(li),range(li[0],li[0]+ll[i]+1))
                    for i,li in enumerate(r)]
                res=[]
                for k in range(prod(len(l) for l in rr)):
                        ll=[0]*m
                        r0=[]
                        for x in sigma:
                                r0+=[rr[x-1]
                                     [k/prod(len(rr[x-i-1])
                                             for i in range(1,x))%(len(rr[x-1]))]
                                     [ll[x-1]]]
                                ll[x-1]+=1
                        res+=[r0]
                return res
    
        def _func_S_to_R(self,sigma):
                from sage.combinat.subset import Subsets
                from sage.combinat.packed_words import to_pack
                from sage.combinat.packed_words import PackedWord
                from sage.functions.other import factorial
                from sage.rings.all import ZZ
                from sage.misc.misc_c import prod
                
                if len(sigma)==0:
                        return self.one()
                res = []
                l=sigma.to_composition()
                for ll in self.sub_composition(l): # l=[3,2] truc(l) = [[0,0],[0,1],[1,0],[1,1],[2,0],[2,1]]
                        r0=ZZ(1)/prod(factorial(k+1) for k in ll)
                        r1=self.untass_pw(sigma,ll)
                        res+=[(r0,to_pack(rr1)) for rr1 in r1]
                return sum(c*self.monomial(PackedWord(pw)) for (c,pw) in res)

        def S_to_R(self):
                morph = lambda X, Y, func, tri = None, comp = None: (
                        X._module_morphism(func, codomain=Y, triangular=tri, cmp=comp)
                        if comp is not None else
                        X._module_morphism(func, codomain=Y, triangular=tri))
                R = self
                S = self.realization_of().S()
                
                # S to R and back
                S_to_R = morph(S, R,
                               lambda sigma: self._func_S_to_R(sigma),
                               tri="lower")
                S_to_R.register_as_coercion()
                (~S_to_R).register_as_coercion()
    


        def internal_product_on_basis(self,sigma,mu):
                R=self
                S=self.realization_of().S()
                return R(S(R(sigma)).internal_product(S(R(mu))))
