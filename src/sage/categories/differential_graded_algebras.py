"""
Differential graded algebras

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
from sage.misc.lazy_import import LazyImport
from sage.misc.cachefunc import cached_method

class DifferentialGradedAlgebras(DifferentialAlgebrasCategory):
    """
    The category of differential graded algebras.

    EXAMPLES::

        sage: from sage.categories.differential_graded_algebras import DifferentialGradedAlgebras
        sage: DifferentialGradedAlgebras(QQ)
        Category of differential graded algebras with basis over Rational Field
        sage: DifferentialGradedAlgebras(QQ).super_categories()
        [Category of graded algebras with basis over Rational Field,
         Category of chain complexes over Rational Field]

    TESTS::

        sage: TestSuite(DifferentialGradedAlgebras(QQ)).run()
    """
    class CartesianProducts(CartesianProductsCategory):
        def extra_super_categories(self):
            """
            A cartesian product of differential graded algebras is endowed
            with a natural graded algebra structure and differential.

            EXAMPLES::

                sage: C = Algebras(QQ).Differential().Graded().CartesianProducts()
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

                sage: Algebras(QQ).Differential().Graded().TensorProducts().extra_super_categories()
                [Category of algebras over Rational Field]
                sage: Algebras(QQ).Differential().Graded().TensorProducts().super_categories()
                [Category of algebras over Rational Field]

            Meaning: a tensor product of differential algebras is a
            differential graded algebra.
            """
            return [self.base_category()]

