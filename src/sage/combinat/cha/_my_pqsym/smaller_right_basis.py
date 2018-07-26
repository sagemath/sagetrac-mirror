# -*- coding: utf-8 -*-
r"""
The smaller right basis of PQSym for Weak Order
"""

from my_pqsym import ParkingQuasiSymmetricFunctions

class SmallerRight(ParkingQuasiSymmetricFunctions.Bases.Base):

    _prefix = "SR"

#    def dual_basis(self):
#        return self.realisation_of().SRD()
        
    def product_on_basis(self, e1, e2):
        '''
        TESTS::
        
        sage: SR = WQSym(QQ).SR()
        sage: SR[1,1]*SR[1]
        SR[2, 1, 1]
        '''
        if len(e1) == 0:
            return self.monomial(e2)
        n1=len(e1)
        return self.monomial(
                self.basis().keys()([i + n1 for i in e2] + list(e1)))
        
    def build_morphisms(self):
        '''

        '''
        self.morph_R_SR()

    def morph_R_SR(self):
        '''
        '''
        from sage.combinat.parking_functions import ParkingFunction_class
        from sage.combinat.parking_functions import ParkingFunction
        
        morph = lambda R, T, func, tri = None, comp = None: (
            R._module_morphism(func, codomain=T, triangular=tri, cmp=comp)
            if comp is not None else
            R._module_morphism(func, codomain=T, triangular=tri))

        R = self.realization_of().R()
        SR = self

        # R <-> SR
        SR_to_R = morph(SR, R,
            lambda sigma: R.sum_of_monomials(
                ParkingFunction_class(list(sigma)).multi_right_permutohedron_smaller()),
            tri="upper")
        SR_to_R.register_as_coercion()
        (~SR_to_R).register_as_coercion()
            
        # morph = lambda R, T, func, tri = None, comp = None: (
        #     R._module_morphism(func, codomain=T, triangular=tri, cmp=comp)
        #     if comp is not None else
        #     R._module_morphism(func, codomain=T, triangular=tri))

        # R = self.realization_of().R()
        # SR = self

        # # SR <-> R
        # SR_to_R = morph(SR, R,
        #     lambda sigma: R.sum_of_monomials(sigma.multi_right_permutohedron_smaller()),
        #     tri="lower")
        # SR_to_R.register_as_coercion()
        # (~SR_to_R).register_as_coercion()

