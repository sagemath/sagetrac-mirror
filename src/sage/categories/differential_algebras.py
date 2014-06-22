"""
Differential algebras

AUTHORS:

- Miguel Marco, John Palmieri, Travis Scrimshaw (2014-06-21): Initial version
"""
#*****************************************************************************
#  Copyright (C) 2014 Miguel Marco, John Palmieri, Travis Scrimshaw
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.algebras import Algebras
from sage.categories.category_types import ChainComplexes
from sage.categories.tensor import TensorProductsCategory
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.misc.lazy_import import LazyImport
from sage.misc.cachefunc import cached_method

class DifferentialAlgebras(CategoryWithAxiom_over_base_ring):
    """
    The category of differential algebras.

    EXAMPLES::

        sage: Algebras(QQ).Differential()
        Category of differential algebras with basis over Rational Field
        sage: Algebras(QQ).Differential().super_categories()
        [Category of algebras with basis over Rational Field,
         Category of chain complexes over Rational Field]

    TESTS::

        sage: TestSuite(Algebras((QQ).Differential()).run()
    """
    _base_category_class_and_axiom = (Algebras, "Differential")

    def extra_super_categories(self):
        r"""
        Return the :class:`ChainComplexes` category.

        This method specifies that a differential algebra with
        a basis is a chain complex.

        .. SEEALSO::

            The :ref:`axioms-deduction-rules` section in the
            documentation of axioms

        EXAMPLES::

            sage: C = Algebras(QQ).Differential()
            sage: C.extra_super_categories()
            (Category of chain complexes over Rational Field,)
        """
        return (ChainComplexes(self.base_ring()),)

    def example(self):
        """
        An example of differential algebra.

        EXAMPLES::

            sage: from sage.categories.differential_algebras import DifferentialAlgebras
            sage: DifferentialAlgebras(QQ).example()
            Free commutative differential algebra over Rational Field on generators in degrees 1, 1, 1, 1, 1
        """
        from sage.algebras.differential_algebras.commutative_dga import CommutativeDGA
        return CommutativeDGA(degrees=(1,1,1,1,1), differential=(None, None, None, (((1,1,0,0,0), 1),), (((0,1,1,0,0), 1),)))

    WithBasis = LazyImport("sage.categories.differential_algebras_with_basis",
                           "DifferentialAlgebrasWithBasis", as_name="WithBasis")

    Graded = LazyImport("sage.categories.differential_graded_algebras",
                        "DifferentialGradedAlgebras", as_name="Graded")

    class CartesianProducts(CartesianProductsCategory):
        """
        The category of differential algebras constructed as Cartesian
        products of differential algebras.

        This construction gives the direct product of differential algebras.
        See discussion on:

        - http://en.wikipedia.org/wiki/Direct_product
        """
        def extra_super_categories(self):
            """
            A cartesian product of differential algebras is endowed with a
            natural algebra structure and differential.

            EXAMPLES::

                sage: C = Algebras(QQ).Differential().CartesianProducts()
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

                sage: Algebras(QQ).Differential().TensorProducts().extra_super_categories()
                [Category of algebras over Rational Field]
                sage: Algebras(QQ).Differential().TensorProducts().super_categories()
                [Category of algebras over Rational Field]

            Meaning: a tensor product of differential algebras is a
            differential algebra.
            """
            return [self.base_category()]

