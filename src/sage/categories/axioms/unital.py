"""
Unital Axiom
"""

from sage.misc.cachefunc import cached_method

from axiom import Axiom


class Unital(Axiom):
    
    @cached_method
    def Inverse(self):
        r"""
        Returns the full subcategory of the inverse objects of ``self``.

        An inverse (multiplicative) magma is a unital magma
        such that every element admits both an inverse on the
        left and on the right. Such a magma is also called a
        *loop*.

        .. SEEALSO:: :wikipedia:`Inverse_element`, :wikipedia:`Quasigroup`

        EXAMPLES::

            sage: Magmas().Unital().Inverse()
            Category of inverse unital magmas
            sage: Monoids().Inverse()
            Category of groups

        TESTS::

            sage: TestSuite(Magmas().Unital().Inverse()).run()
            sage: Algebras(QQ).Inverse.__module__
            'sage.categories.magmas'
        """
        from inverse import Inverse
        return self._with_axiom(Inverse)


class AdditiveUnital(Axiom): 
    pass
