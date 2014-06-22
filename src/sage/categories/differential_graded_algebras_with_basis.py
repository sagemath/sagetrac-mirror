"""
Differential graded algebras with basis

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
from sage.misc.cachefunc import cached_method

class DifferentialGradedAlgebrasWithBasis(DifferentialAlgebrasCategory):
    """
    The category of differential graded algebras with a distinguished basis.

    EXAMPLES::

        sage: Algebras(QQ).Differential().Graded().WithBasis()
        Category of differential graded algebras with basis over Rational Field
        sage: Algebras(QQ).Differential().Graded().WithBasis().super_categories()
        [Category of graded algebras with basis over Rational Field,
         Category of chain complexes over Rational Field]

    TESTS::

        sage: TestSuite(Algebras(QQ).Differential().Graded().WithBasis()).run()
    """
    class CartesianProducts(CartesianProductsCategory):
        def extra_super_categories(self):
            """
            A cartesian product of differential graded algebras with a
            distinguished basis is endowed with a natural graded algebra
            structure, a distinguished basis, and differential.

            EXAMPLES::

                sage: C = Algebras(QQ).Differential().Graded().WithBasis().CartesianProducts()
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

                sage: C = Algebras(QQ).Differential().Graded().WithBasis().TensorProducts()
                sage: C.extra_super_categories()
                [Category of algebras over Rational Field]
                sage: C.super_categories()
                [Category of algebras over Rational Field]

            Meaning: a tensor product of differential graded algebras with
            basis is a differential graded algebra with basis.
            """
            return [self.base_category()]

