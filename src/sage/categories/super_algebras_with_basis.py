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
        The category of algebras with basis constructed by tensor product of super algebras with basis
        """

        class ParentMethods:
            """
            implements operations on tensor products of algebras with basis
            """

            def product_on_basis(self, t0, t1):
                """
                The product of the algebra on the basis, as per
                ``AlgebrasWithBasis.ParentMethods.product_on_basis``.

                EXAMPLES:

                Test the sign in the super tensor product -- see :trac:`25603`::

                    sage: A = SteenrodAlgebra(3)
                    sage: x = A.Q(0)
                    sage: y = x.coproduct()
                    sage: y^2
                    0

                TODO: optimize this implementation!
                """
                basic = tensor((module.monomial(x0)*module.monomial(x1)
                                for (module, x0, x1) in zip(self._sets, t0, t1)))
                n = len(self._sets)
                parity1 = [self._sets[0].degree_on_basis(x0) for x0 in t0]
                parity2 = [self._sets[1].degree_on_basis(x1) for x1 in t1]
                parity = sum(parity1[i] * parity2[j]
                           for j in range(n) for i in range(j+1,n))
                return (-1)**parity * basic
