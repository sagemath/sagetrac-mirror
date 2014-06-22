"""
Differential algebras with basis

AUTHORS:

- Miguel Marco, John Palmieri, Travis Scrimshaw (2014-06-21): Initial version
"""
#*****************************************************************************
#  Copyright (C) 2014 Miguel Marco, John Palmieri, Travis Scrimshaw
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.differential_algebras import DifferentialAlgebrasCategory
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.tensor import TensorProductsCategory
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute

class DifferentialAlgebrasWithBasis(DifferentialAlgebrasCategory):
    """
    The category of differential algebras with a distinguished basis.

    EXAMPLES::

        sage: Algebras(QQ).Differential().WithBasis()
        Category of differential algebras with basis over Rational Field
        sage: Algebras(QQ).Differential().WithBasis().super_categories()
        [Category of algebras with basis over Rational Field,
         Category of chain complexes over Rational Field]

    TESTS::

        sage: TestSuite(Algebras(QQ).Differential().WithBasis()).run()
    """
    class ParentMethods:
        @abstract_method(optional=True)
        def differential_on_basis(self, t):
            """
            The differential of the algebra on the basis (optional).

            INPUT:

            - ``t`` -- the indices of an element of the basis of ``self``

            Returns the differential of the corresponding basis element.
            If implemented, the differential of the algebra is defined
            from it by linearity.

            EXAMPLES::

                sage: from sage.algebras.differential_algebras.commutative_dga import CommutativeDGA
                sage: D = CommutativeDGA(degrees=(1,2), differential=(None, (((1,1), 1),)))
                sage: D.gens()
                (a_1, b_2)
                sage: D.differential_on_basis((0,1))
                a_1 * b_2
            """
            pass

        @lazy_attribute
        def differential(self):
            """
            The differential of this algebra.

            If :meth:`.differential_basis` is available, this
            constructs the differential morphism from ``self``
            to ``self`` by extending it by linearity. Otherwise,
            :meth:`self.differential_by_coercion` is used, if
            available.

            EXAMPLES::

                sage: from sage.algebras.differential_algebras.commutative_dga import CommutativeDGA
                sage: D = CommutativeDGA(degrees=(1,2), differential=(None, (((1,1), 1),)))
                sage: D.gens()
                (a_1, b_2)
                sage: D.differential(D.gen(1))
                a_1 * b_2
            """
            if self.differential_on_basis is not NotImplemented:
                return self._module_morphism(self.differential_on_basis, codomain = self)
            elif hasattr(self, "differential_by_coercion"):
                return self.differential_by_coercion

    class ElementMethods:
        pass

    class CartesianProducts(CartesianProductsCategory):
        def extra_super_categories(self):
            """
            A cartesian product of differential algebras is endowed with a
            natural algebra structure and differential.

            EXAMPLES::

                sage: C = Algebras(QQ).WithBasis().Differential().CartesianProducts()
                sage: C.extra_super_categories()
                [Category of algebras over Rational Field]
                sage: sorted(C.super_categories(), key=str)
                [Category of Cartesian products of commutative additive groups,
                 Category of Cartesian products of distributive magmas and additive magmas,
                 Category of Cartesian products of semigroups,
                 Category of Cartesian products of unital magmas,
                 Category of algebras over Rational Field]
            """
            return [self.base_category()]

    class TensorProducts(TensorProductsCategory):
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: C = Algebras(QQ).WithBasis().Differential().TensorProducts()
                sage: C.extra_super_categories()
                [Category of algebras over Rational Field]
                sage: C.super_categories()
                [Category of algebras over Rational Field]

            Meaning: a tensor product of differential algebras is a
            differential algebra.
            """
            return [self.base_category()]

