# -*- coding: utf-8 -*-
r"""
Dialgebras

A dialgebra _[Loday] over a ring `R` is a vector space D equipped with two
associative operation `\left` and `\right`, called respectively left and right
product, statisfying 3 more axioms:

MATH::

    \begin{align}
        x \left ( y \left z ) &= x \left ( y \right z )\\
        ( x \right y ) \left z &= x \right ( y \left  z )\\
        ( x \left y ) \right z &= ( x \right y ) \right z
    \end{align}

References:

.. [Loday] Dialgebras,
    Jean-Louis Loday

AUTHORS:

 - Jean-Baptiste Priez (first version)
"""
#*****************************************************************************
#  Copyright (C) 2013      Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.modules import Modules
from sage.misc.abstract_method import abstract_method
from sage.categories.realizations import RealizationsCategory
from sage.categories.category_with_axiom import \
    CategoryWithAxiom_over_base_ring
from sage.categories.category_types import Category_over_base_ring

#from sage.categories.algebras_over_operads import AlgebrasOverOperad


class Dialgebras(Category_over_base_ring):
    """
    The category of dialgebras over a given base ring.

    A dialgebra _[Loday] over a ring `R` is a vector space D equipped with two
    associative operation `\left` and `\right`, called respectively left and
    right product, statisfying 3 more axioms:

    MATH::

        \begin{align}
            x \left ( y \left z ) &= x \left ( y \right z )\\
            ( x \right y ) \left z &= x \right ( y \left  z )\\
            ( x \left y ) \right z &= ( x \right y ) \right z
        \end{align}

    EXAMPLES::

        sage: from sage.combinat.hopf_algebras.categories.dialgebras import Dialgebras
        sage: Dialgebras(ZZ)
        Category of dialgebras over Integer Ring

    TESTS::

        sage: from sage.combinat.hopf_algebras.categories.dialgebras import Dialgebras
        sage: TestSuite(Dialgebras(ZZ)).run()
    """

    def super_categories(self):
        """
        TESTS::

            sage: from sage.combinat.hopf_algebras.categories.dialgebras import Dialgebras
            sage: Dialgebras(ZZ).super_categories()
            [Category of graded modules over Integer Ring]
        """
        R = self.base_ring()
        return [Modules(R).Graded()] #AlgebrasOverOperad(R)]

    class ParentMethods:

        def _test_dialgebras(self, **options):
            r"""
            Run generic tests on the methods :meth:`.left_product` and
            `.right_product'.

            See also: :class:`TestSuite`.

            """
            # TODO: tests
            from sage.categories.sets_cat import EmptySetError
            tester = self._tester(**options)
            try:
                S = self.some_elements()
            except EmptySetError:
                return
            for x in S:
                for y in S:
                    for z in S:
                        tester.assertEqual(
                            x << (y << z),
                            x << (y >> z)
                        )
                        tester.assertEqual(
                            (x >> y) << z,
                            x >> (y << z)
                        )
                        tester.assertEqual(
                            (x << y) >> z,
                            (x >> y) >> z
                        )

        @abstract_method
        def left_product(self, l_elem, r_elem):
            r"""
            The left product of two elements.
            """

        @abstract_method
        def right_product(self, l_elem, r_elem):
            r"""
            The right product of two elements.
            """

    class ElementMethods:

        def __rshift__(self, right):
            """
            Shorthand for the right product
            """
            # TODO: tests
            self.parent().right_product(self, right)

        def __lshift__(self, right):
            """
            Shorthand for the left product
            """
            # TODO: tests
            self.parent().left_product(self, right)

    class Realizations(RealizationsCategory):

        class ParentMethods:

            pass

    class WithBasis(CategoryWithAxiom_over_base_ring):

        class ParentMethods:

            pass

from sage.categories.operators_tools import expand_binary_operator


expand_binary_operator(Dialgebras, "left_product")
expand_binary_operator(Dialgebras, "right_product")
