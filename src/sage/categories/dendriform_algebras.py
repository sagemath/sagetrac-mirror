# -*- coding: utf-8 -*-
r"""
Dendriform algebras

A dendriform algebra `A` over a field `K` is a (non-unitary) associative
algebra, such that the product `\times` can be split into two part `\prec` and
`\succ` (called left and right products) with good compatibilities (which mean
that `(A, \prec, \succ)` is a bimodule over it self). (see _[Loday], _[Aguiar],
_[Arith], _[PermBT] and _[Foissy])


A dendriform algebra `(A, \prec, \succ)` is:
 - A is `K`-vector space and

    MATH::

         \prec : \left\{\begin{array}{l}
             A \otimes A \to \A \\
             x \otimes y \to x \prec y
         \end{array}\right. \text{ and }
         \succ : \left\{\begin{array}{l}
             A \otimes A \to \A \\
             x \otimes y \to x \succ y
         \end{array}\right.

 such that for all `x, y, z \in A`:

    MATH::

        \begin{align}
            ( x \prec y ) \prec z &= x \prec ( y \prec z )\\
            ( x \succ y ) \prec z &= x \succ ( y \prec  z )\\
            ( x \times y ) \succ z &= x \succ ( y \succ z )
        \end{align}

And if we put `\times := \prec + \succ` the structure `(A, m)` is a non-unitary
associative algebra.

References:

.. [Loday] Dialgebras,
    Jean-Louis Loday

.. [Aguiar] Infinitesimal bialgebras, pre-lie and dendriform algebras,
    Marcelo Aguiar

.. [Arith] Arithmetree,
    J. Algebra (258) - 2002

.. [PermBT] Order structure on the algebra of permutations and of planar binary
    trees,
    J. Algebraic Combin. (15) - 2002

.. [Foissy] Bidendriform bialgebras, trees, and free quasi-symmetric functions,
    Loïc Foissy

AUTHORS:

 - Jean-Baptiste Priez
"""
#*****************************************************************************
#  Copyright (C) 2013   Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.category_types import Category_over_base_ring
from sage.categories.modules import Modules
from sage.misc.abstract_method import abstract_method
from sage.categories.dual import DualObjectsCategory
from sage.categories.realizations import RealizationsCategory

from sage.categories.category_with_axiom import \
    CategoryWithAxiom_over_base_ring
from sage.categories.with_realizations import WithRealizationsCategory
from sage.categories.associative_algebras import AssociativeAlgebras
from sage.misc.cachefunc import cached_method


class DendriformAlgebras(Category_over_base_ring):
    """
    The category of dendriform algebras over a given base ring.

    A dendriform algebra over a ring `R` is an algebra with `\times` is a
    product defined by `\times := \prec + \succ` two operator such that respect
    following properties, for `x,y,z` three basis element:

    MATH::

        \begin{align}
            ( x \prec y ) \prec z &= x \prec ( y \prec z )\\
            ( x \succ y ) \prec z &= x \succ ( y \prec  z )\\
            ( x \times y ) \succ z &= x \succ ( y \succ z )
        \end{align}

    TODO: should `R` be a commutative ring?

    EXAMPLES::

        sage: from sage.combinat.hopf_algebras.categories.dendriform_algebras import \
        ....:     DendriformAlgebras
        sage: DendriformAlgebras(ZZ)
        Category of dendriform algebras over Integer Ring

    TESTS::

        sage: TestSuite(DendriformAlgebras(ZZ)).run()

    """

    @cached_method
    def super_categories(self):
        """
        TESTS::

            sage: from sage.categories.dendriform_algebras import DendriformAlgebras
            sage: DendriformAlgebras(QQ)
            Category of dendriform algebras over Rational Field
        """
        R = self.base_ring()
        return [Modules(R)]

    class ParentMethods:

        def _test_dendriform(self, **options):
            r"""
            Run generic tests on the methods :meth:`.left_product` and
            `.right_product'.

            See also: :class:`TestSuite`.

            """
            from sage.categories.sets_cat import EmptySetError
            tester = self._tester(**options)
            try:
                S = self.some_elements()
            except EmptySetError:
                return
            if AssociativeAlgebras(self.base_ring()).Unital() in self.categories():
                minus = lambda x: x - x[self.one().support()[0]] * self.one()
            else:
                minus = lambda x: x
            for x in S:
                x = minus(x)
                for y in S:
                    y = minus(y)
                    for z in S:
                        z = minus(z)
                        tester.assertEqual(
                            (x << y) << z,
                            x << ((y << z) + (y >> z)),
                            "(x << y) << z != x << ((y << z) + (y >> z))"
                        )
                        tester.assertEqual(
                            (x >> y) << z,
                            x >> (y << z),
                            "(x >> y) << z != x >> (y << z)"
                        )
                        tester.assertEqual(
                            ((x << y) + (x >> y)) >> z,
                            x >> (y >> z),
                            "((x << y) + (x >> y)) >> z != x >> (y >> z)"
                        )

        @abstract_method
        def left_product(self, l_elem, r_elem):
            r"""
            The left product (`\prec`) of two elements.

            ##TESTS::

             ##   sage: F = FQSym(QQ).F()
             ##   sage: F[1] << F[1] == F.left_product(F[1], F[1])
             ##  True
            """

        @abstract_method
        def right_product(self, l_elem, r_elem):
            r"""
            The right product (`\succ`) of two elements.

            ##TESTS::

              #  sage: F = FQSym(QQ).F()
              #  sage: F[1] >> F[1] == F.right_product(F[1], F[1])
              #  True
            """

        def m_product(self, l, r):
            """
            By default a dendriform algebra is a non-unital algebra
            with the product `\times` is defined for any `x,y` by:

            MATH::

                x \times y := x \prec y + x \succ y
            """
            return self.right_product(l, r) + self.left_product(l, r)

    class ElementMethods:

        def __rshift__(self, right):
            """
            Shorthand for the right product

            ##TESTS::

              #  sage: F = FQSym(QQ).F()
              #  sage: F[1] >> F[1] == F[1].right_product(F[1])
              #  True
            """
            return self.parent().right_product(self, right)

        def __lshift__(self, right):
            """
            Shorthand for the left product

            ##TESTS::

              #  sage: F = FQSym(QQ).F()
              #  sage: F[1] << F[1] == F[1].left_product(F[1])
              #  True
            """
            return self.parent().left_product(self, right)

    class Realizations(RealizationsCategory):

        class ParentMethods:

            def m_product(self, l, r):
                """
                (See :meth:`DendriformAlgebras.ParentMethods.product`.)
                """
                return self.left_product(l, r) + self.right_product(l, r)

    class WithRealizations(WithRealizationsCategory):

        class ParentMethods:
            pass

    class WithBasis(CategoryWithAxiom_over_base_ring):

        class ParentMethods:
            pass

    class DualObjects(DualObjectsCategory):

        def extra_super_categories(self):
            r"""
            Returns the dual category

            EXAMPLES:

            The category of algebras over the Rational Field is dual
            to the category of coalgebras over the same field::

                sage: from sage.combinat.hopf_algebras.categories.dendriform_algebras \
                ....:     import DendriformAlgebras
                sage: DendriformAlgebras(QQ).dual()
                Category of duals of dendriform algebras over Rational Field
                sage: DendriformAlgebras(QQ).dual().extra_super_categories()
                [Category of dendriform coalgebras over Rational Field]
            """
            from sage.categories.dendriform_coalgebras import DendriformCoalgebras
            return [DendriformCoalgebras(self.base_category().base_ring())]

from sage.categories.operators_tools import expand_binary_operator

expand_binary_operator(DendriformAlgebras, "left_product")
expand_binary_operator(DendriformAlgebras, "right_product")
