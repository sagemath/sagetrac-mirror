# -*- coding: utf-8 -*-
r"""
The elementary dual basis of WQSym Hopf algebra.
"""

from sage.combinat.cha.wqsym import WordQuasiSymmetricFunctions


class GreaterLeftDual(WordQuasiSymmetricFunctions.Bases.Base):
    '''
    '''
    _prefix = "GLD"

    def dual_basis(self):
        return self.realization_of().WE()

    def build_morphisms(self):
        '''
        GLD(u):=somme pour u <=L v de muL(u,v)Rv
        TESTS::

            sage: WQS = WQSym(QQ); R = WQS.R(); GLD = WQS.GLD()
            sage: R(GLD[2,1,3])
faire des tests et voir si c'est la bonne formule!!
        '''
        morph = lambda X, Y, func, tri = None, comp = None: (
            X._module_morphism(func, codomain=Y, triangular=tri, cmp=comp)
            if comp is not None else
            X._module_morphism(func, codomain=Y, triangular=tri))

        R = self.realization_of().R()
        GLD = self

        # GLD to R and back
        R_to_GLD = morph(R, GLD,
            lambda sig: GLD.sum_of_monomials(
                sig.left_weak_order_greater()),
            tri="lower")
        R_to_GLD.register_as_coercion()
        (~R_to_GLD).register_as_coercion()


    def coproduct_on_basis(self, sigma):
        '''
        '''
        from sage.combinat.packed_words import to_pack
        return self.tensor_square().sum_of_monomials(
            (to_pack(sigma[:i]), to_pack(sigma[i:]))
            for i in range(len(sigma) + 1)
            if i==0 or i==len(sigma) or
            min(sigma[:i])>max(sigma[i:])
        )

        
