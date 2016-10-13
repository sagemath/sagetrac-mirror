r"""
Super algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.super_modules import SuperModulesCategory
from sage.categories.algebras import Algebras
from sage.categories.modules import Modules
from sage.categories.tensor import TensorProductsCategory, tensor

class SuperAlgebrasWithBasis(SuperModulesCategory):
    """
    The category of super algebras with a distinguished basis

    EXAMPLES::

        sage: C = Algebras(ZZ).WithBasis().Super(); C
        Category of super algebras with basis over Integer Ring

    TESTS::

        sage: TestSuite(C).run()
    """
    def extra_super_categories(self):
        """
        EXAMPLES::

            sage: C = Algebras(ZZ).WithBasis().Super()
            sage: sorted(C.super_categories(), key=str) # indirect doctest
            [Category of graded algebras with basis over Integer Ring,
             Category of super algebras over Integer Ring,
             Category of super modules with basis over Integer Ring]
        """
        return [self.base_category().Graded()]

    class ParentMethods:
        def graded_algebra(self):
            r"""
            Return the associated graded module to ``self``.

            See :class:`~sage.algebras.associated_graded.AssociatedGradedAlgebra`
            for the definition and the properties of this.

            .. SEEALSO::

                :meth:`~sage.categories.filtered_modules_with_basis.ParentMethods.graded_algebra`

            EXAMPLES::

                sage: W.<x,y> = algebras.DifferentialWeyl(QQ)
                sage: W.graded_algebra()
                Graded Algebra of Differential Weyl algebra of
                 polynomials in x, y over Rational Field
            """
            from sage.algebras.associated_graded import AssociatedGradedAlgebra
            return AssociatedGradedAlgebra(self)


    class TensorProducts(TensorProductsCategory):
        """
        The category of super algebras with basis constructed by tensor product of superalgebras with basis
        """
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: SuperAlgebrasWithBasis(QQ).TensorProducts().extra_super_categories()
                [Category of super algebras with basis over Rational Field]
                sage: SuperAlgebrasWithBasis(QQ).TensorProducts().super_categories()
                [Category of super algebras with basis over Rational Field,
                 Category of tensor products of algebras with basis over Rational Field]
            """
            return [self.base_category()]

        class ParentMethods:
            """
            implements operations on tensor products of superalgebras with basis
            """

            def product_on_basis(self, t1, t2):
                """
                The product of the superalgebra on the basis, as per
                ``SuperalgebrasWithBasis.ParentMethods.product_on_basis``.
                """

                l = (module.monomial(x1)*module.monomial(x2) for (module, x1, x2) in zip(self._sets, t1, t2))

                return (-1)^(x1.degree()*x2.degree())*tensor(l)
