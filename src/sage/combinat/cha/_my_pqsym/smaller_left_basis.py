# -*- coding: utf-8 -*-
r"""
The homogene basis of WQSym for Weak Order
"""

from my_pqsym import ParkingQuasiSymmetricFunctions


class SmallerLeft(ParkingQuasiSymmetricFunctions.Bases.Base):

    _prefix = "SL"

    def dual_basis(self):
        return self.realisation_of().SLD()
        
    def product_on_basis(self, e1, e2):
        '''
        TESTS::
        
        sage: SL = WQSym(QQ).SL()
        sage: SL[1,1]*SL[1]
        SL[2, 2, 1]
        '''
        if len(e2) == 0:
            return self.monomial(e1)
        m2=max(e2)
        return self.monomial(
                self.basis().keys()([i + m2 for i in e1] + list(e2)))
        
    def build_morphisms(self):
        '''

        '''
        self.morph_L_SL()

    def morph_L_SL(self):
        '''
        '''
        from sage.combinat.parking_functions import ParkingFunction_class
        from sage.combinat.parking_functions import ParkingFunction
        
        morph = lambda L, T, func, tri = None, comp = None: (
            L._module_morphism(func, codomain=T, triangular=tri, cmp=comp)
            if comp is not None else
            L._module_morphism(func, codomain=T, triangular=tri))

        L = self.realization_of().L()
        SL = self

        # L <-> SL
        SL_to_L = morph(SL, L,
            lambda sigma: L.sum_of_monomials(
                ParkingFunction_class(list(sigma)).multi_left_permutohedron_smaller()),
            tri="upper")
        SL_to_L.register_as_coercion()
        (~SL_to_L).register_as_coercion()
