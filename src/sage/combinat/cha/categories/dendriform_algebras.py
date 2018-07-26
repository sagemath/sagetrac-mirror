# -*- coding: utf-8 -*-
r"""
Dendriform algebras

AUTHORS:

 - Rémi Maurice & Jean-Baptiste Priez (first version)
"""
#*****************************************************************************
#  Copyright (C) 2013      Rémi Maurice <maurice@univ-mlv.fr>
#                          Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_types import Category_over_base_ring
from sage.categories.dual import DualObjectsCategory
from sage.categories.algebras import Algebras
from sage.structure.element import have_same_parent
from sage.categories.realizations import RealizationsCategory
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.misc.lazy_attribute import lazy_attribute


class DendriformAlgebras(Category_over_base_ring):
    """
    The category of dendriform algebras over a given base ring.

    An dendriform algebra over a ring `R` is an algebra with `\times` is a
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

        sage: from sage.combinat.cha.tools.dendriform_algebras import \
        ....:     DendriformAlgebras
        sage: DendriformAlgebras(ZZ)
        Category of dendriform algebras over Integer Ring

    TESTS::

        sage: TestSuite(DendriformAlgebras(ZZ)).run()
    """

    def super_categories(self):
        """
        TESTS::

            sage: from sage.combinat.cha.tools.dendriform_algebras import \
            ....:     DendriformAlgebras
            sage: DendriformAlgebras(ZZ).super_categories()
            [Category of algebras over Integer Ring]
        """
        R = self.base_ring()
        return [Algebras(R)]

    class ParentMethods:

        def _test_dendriform(self, **options):
            """
            TODO:: how to do test of test_methods
            """
            tester = self._tester(**options)
            S = tester.some_elements()
            for x in S.some_elements():
                for y in S.some_elements():
                    for z in S.some_elements():
                        tester.assertEqual(
                            (x << y) << z,
                            x << (y * z)
                        )
                        tester.assertEqual(
                            (x >> y) << z,
                            x >> (y << z)
                        )
                        tester.assertEqual(
                            (x * y) >> z,
                            x >> (y >> z)
                        )

        def left_product(self, l_elem, r_elem):
            return l_elem.left_product(r_elem)

        def right_product(self, l_elem, r_elem):
            return l_elem.right_product(r_elem)

        def product(self, l_elem, r_elem):
            return self.left_product(l_elem, r_elem) + \
                   self.right_product(l_elem, r_elem)

    class ElementMethods:

        def right_product(self, Nelt):
            return self.__rshift__(Nelt)

        def left_product(self, Nelt):
            return self.__lshift__(Nelt)

        def __rshift__(self, right):
            '''
            .. FIXME:: must be done
            attention.... c + b << a == (c+b)<<a
            '''
            if have_same_parent(self, right) and\
               hasattr(self.parent(), "right_product"):
                return self.parent().right_product(self, right)
            from sage.structure.element import get_coercion_model
            import operator
            return get_coercion_model().bin_op(self, right, operator.rshift)

        def __lshift__(self, right):
            '''
            .. FIXME:: must be done
            '''
            if have_same_parent(self, right) and \
               hasattr(self.parent(), "left_product"):
                return self.parent().left_product(self, right)
            from sage.structure.element import get_coercion_model
            import operator
            return get_coercion_model().bin_op(self, right, operator.lshift)

    class Realizations(RealizationsCategory):

        class ParentMethods:

            def left_product_by_coercion(self, left, right):
                r"""
                TODO::
                """
                for R in self.realization_of().realizations():
                    self_to_R = R.coerce_map_from(self)
                    R_to_self = self.coerce_map_from(R)
                    if self_to_R is not None and R_to_self is not None and \
                       R.left_product != R.left_product_by_coercion:
                        return self(R(self(left)) << R(self(right)))
                return NotImplementedError

            def right_product_by_coercion(self, left, right):
                r"""
                TODO::
                """
                for R in self.realization_of().realizations():
                    self_to_R = R.coerce_map_from(self)
                    R_to_self = self.coerce_map_from(R)
                    if self_to_R is not None and R_to_self is not None and \
                       R.right_product != R.right_product_by_coercion:
                        return self(R(self(left)) >> R(self(right)))
                return NotImplementedError

    class WithBasis(Category_over_base_ring):

        @cached_method
        def super_categories(self):
            R = self.base_ring()
            return [DendriformAlgebras(R), AlgebrasWithBasis(R)]

        def __repr__(self):
            return "Category of dendriform algebras with basis over %s" % \
                    self.base_ring()

        class ParentMethods:

            @abstract_method(optional=True)
            def left_product_on_basis(self, left, right):
                """
                TODO::
                """

            @abstract_method(optional=True)
            def right_product_on_basis(self, left, right):
                """
                TODO::
                """

            # .. TODO:: add a lazy product which is left + right product?

            @lazy_attribute
            def left_product(self):
                if self.left_product_on_basis is NotImplemented:
                    if hasattr(self, "left_product_by_coercion"):
                        return self.left_product_by_coercion
                    else:
                        return NotImplemented
                return self._left_product_from_product_on_basis_multiply

            @lazy_attribute
            def right_product(self):
                if self.right_product_on_basis is NotImplemented:
                    if hasattr(self, "right_product_by_coercion"):
                        return self.right_product_by_coercion
                    else:
                        return NotImplemented
                return self._right_product_from_product_on_basis_multiply

            def _left_product_from_product_on_basis_multiply(
                    self, left, right):
                r"""
                """
                return self.linear_combination(
                    (self.left_product_on_basis(mon_left, mon_right),
                     coeff_left * coeff_right)
                    for (mon_left, coeff_left) in left.\
                        monomial_coefficients().iteritems()
                    for (mon_right, coeff_right) in right.\
                        monomial_coefficients().iteritems()
                )

            def _right_product_from_product_on_basis_multiply(
                    self, left, right):
                r"""

                """
                return self.linear_combination(
                    (self.right_product_on_basis(mon_left, mon_right),
                     coeff_left * coeff_right)
                    for (mon_left, coeff_left) in left.\
                        monomial_coefficients().iteritems()
                    for (mon_right, coeff_right) in right.\
                        monomial_coefficients().iteritems()
                )

    class DualObjects(DualObjectsCategory):

        def extra_super_categories(self):
            r"""
            Returns the dual category

            EXAMPLES:

            The category of algebras over the Rational Field is dual
            to the category of coalgebras over the same field::

                sage: C = Algebras(QQ)
                sage: C.dual()
                Category of duals of algebras over Rational Field
                sage: C.dual().extra_super_categories()
                [Category of coalgebras over Rational Field]
            """
            from dendriform_coalgebras import DendriformCoalgebras
            return [DendriformCoalgebras(self.base_category().base_ring())]
