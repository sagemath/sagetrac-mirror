# -*- coding: utf-8 -*-
r"""
The left elementary basis of WQSym for Weak Order
"""

from my_pqsym import ParkingQuasiSymmetricFunctions

class GreaterRight(ParkingQuasiSymmetricFunctions.Bases.Base):
    '''

    '''
    _prefix = "GR"

    def dual_basis(self):
        return self.realisation_of().GRD()

    def product_on_basis(self, e1, e2):
        '''
            TESTS::

                sage: GR = WQSym(QQ).GR()
                sage: GR[1,1]*GR[1]
                GR[1, 1, 2]
        '''
        if len(e1) == 0:
            return self.monomial(e2)
        m1=max(e1)
        return self.monomial(
            self.basis().keys()(list(e1) + [i + m1 for i in e2])
        )

    def build_morphisms(self):
        '''
        TESTS::

        '''
        self.morph_R_GR()
#        self.morph_SL_GR()

    def morph_R_GR(self):
        '''
        '''
        from sage.combinat.parking_functions import ParkingFunction_class
        from sage.combinat.parking_functions import ParkingFunction
        
        morph = lambda R, T, func, tri = None, comp = None: (
            R._module_morphism(func, codomain=T, triangular=tri, cmp=comp)
            if comp is not None else
            R._module_morphism(func, codomain=T, triangular=tri))

        R = self.realization_of().R()
        GR = self

        # R <-> GR
        GR_to_R = morph(GR, R,
            lambda sigma: R.sum_of_monomials(
                ParkingFunction_class(list(sigma)).multi_right_permutohedron_greater()),
            tri="lower")
        GR_to_R.register_as_coercion()
        (~GR_to_R).register_as_coercion()
            
    # def morph_SL_GR(self):
    #     '''
    #     '''
    #     # morph = lambda X, T, func, tri = None, comp = None: (
    #     #     X._module_morphism(func, codomain=T, triangular=tri, cmp=comp)
    #     #     if comp is not None else
    #     #     X._module_morphism(func, codomain=T, triangular=tri))

    #     SL = self.realization_of().SL()
    #     GR = self

    #     # SL <-> GR
    #     # GR_to_SL = morph(GR, SL,
    #     #     lambda sigma: SL.monomial(PackedWord(sigma[::-1])),
    #     #     tri="lower")
    #     # GR_to_SL.register_as_coercion()
    #     # (~GR_to_SL).register_as_coercion()

    #     GR._module_morphism(
    #         lambda sigma:
    #         SL(PackedWord([max(sigma)-i+1 for i in sigma])),
    #         codomain=SL
    #     ).register_as_coercion()

    #     SL._module_morphism(
    #         lambda sigma:
    #         GR(PackedWord([max(sigma)-i+1 for i in sigma])),
    #         codomain=GR
    #     ).register_as_coercion()
