# -*- coding: utf-8 -*-
r"""
The homogene dual basis of WQSym Hopf algebra.
"""

from sage.combinat.cha.wqsym import WordQuasiSymmetricFunctions

class SmallerRightDual (WordQuasiSymmetricFunctions.Bases.Base):
    '''
    '''
    _prefix = "SRD"

    def dual_basis(self):
        return self.realization_of().WH()

    def build_morphisms(self):
        '''
        TESTS::

            sage: WQS = WQSym(QQ); R = WQS.R(); SRD = WQS.SRD()
            sage: R(SRD[2,1,3])
faire des tests et voir si c'est la bonne formule!!
        '''
        morph = lambda X, Y, func, tri = None, comp = None: (
            X._module_morphism(func, codomain=Y, triangular=tri, cmp=comp)
            if comp is not None else
            X._module_morphism(func, codomain=Y, triangular=tri))

        L = self.realization_of().L()
        SRD = self

        # L to SRD and back
        L_to_SRD = morph(L, SRD,
            lambda sig: SRD.sum_of_monomials(
                sig.right_weak_order_greater()),
            tri="lower")
        L_to_SRD.register_as_coercion()
        (~L_to_SRD).register_as_coercion()

        
    def coproduct_on_basis(self, sigma):
        '''
        '''
        from sage.combinat.packed_words import to_pack
        return self.tensor_square().sum_of_monomials(
            (to_pack(sigma[i:]), to_pack(sigma[:i]))
            for i in range(len(sigma) + 1)
            if i==0 or i==len(sigma) or
            min(sigma[:i])>max(sigma[i:])
        )

        
