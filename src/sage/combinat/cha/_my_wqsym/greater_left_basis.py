# -*- coding: utf-8 -*-
r"""
The left elementary basis of WQSym for Weak Order
"""

from sage.combinat.cha.wqsym import WordQuasiSymmetricFunctions
from sage.combinat.packed_words import PackedWord

class GreaterLeft(WordQuasiSymmetricFunctions.Bases.Base):
    '''

    '''
    _prefix = "GL"

    def dual_basis(self):
        return self.realisation_of().GLD()

    def product_on_basis(self, e1, e2):
        '''
            TESTS::

                sage: GL = WQSym(QQ).GL()
                sage: GL[1,1]*GL[1]
                GL[1, 1, 2]
        '''
        if len(e1) == 0:
            return self.monomial(e2)
        m1=max(e1)
        return self.monomial(
            self.basis().keys()(list(e1) + [i + m1 for i in e2])
        )
    
    # def build_morphisms(self):
    #     '''
    #     TESTS::

    #     '''
        
    #     morph = lambda R, T, func, tri = None, comp = None: (
    #         R._module_morphism(func, codomain=T, triangular=tri, cmp=comp)
    #         if comp is not None else
    #         R._module_morphism(func, codomain=T, triangular=tri))

    #     R = self.realization_of().R()
    #     GL = self

    #     # R <-> GL
    #     GL_to_R = morph(GL, R,
    #         lambda sigma: R.sum_of_monomials(sigma.multi_right_permutohedron_greater()),
    #         tri="lower")
    #     GL_to_R.register_as_coercion()
    #     (~GL_to_R).register_as_coercion()

    def build_morphisms(self):
        '''
        TESTS::

        '''
        self.morph_L_GL()
        self.morph_SR_GL()


    def morph_L_GL(self):
        '''
        '''
        morph = lambda L, T, func, tri = None, comp = None: (
            L._module_morphism(func, codomain=T, triangular=tri, cmp=comp)
            if comp is not None else
            L._module_morphism(func, codomain=T, triangular=tri))

        L = self.realization_of().L()
        GL = self
        
        # L <-> GL
        GL_to_L = morph(GL, L,
            lambda sigma: L.sum_of_monomials(sigma.left_weak_order_greater()),
            tri="lower")
        GL_to_L.register_as_coercion()
        (~GL_to_L).register_as_coercion()
            
    def morph_SR_GL(self):
        '''
        '''
        # morph = lambda X, T, func, tri = None, comp = None: (
        #     X._module_morphism(func, codomain=T, triangular=tri, cmp=comp)
        #     if comp is not None else
        #     X._module_morphism(func, codomain=T, triangular=tri))

        SR = self.realization_of().SR()
        GL = self

        # SR <-> GL
        # GL_to_SR = morph(GL, SR,
        #     lambda sigma: SR.monomial(PackedWord(sigma[::-1])),
        #     tri="lower")
        # GL_to_SR.register_as_coercion()
        # (~GL_to_SR).register_as_coercion()

        GL._module_morphism(
            lambda sigma: SR(PackedWord(sigma[::-1])),
            codomain=SR
        ).register_as_coercion()

        SR._module_morphism(
            lambda sigma: GL(PackedWord(sigma[::-1])),
            codomain=GL
        ).register_as_coercion()
