# -*- coding: utf-8 -*-
r"""
 Category of algebras with polynomial realization

AUTHORS:

     - Jean-Baptiste Priez
"""
#*****************************************************************************
#  Copyright (C) 2013      Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.algebras import Algebras
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.realizations import RealizationsCategory
from sage.categories.with_realizations import WithRealizationsCategory
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method


class PolynomialRealizationAlgebras(Category_over_base_ring):
    """
    TESTS::

        sage: TestSuite(PolynomialRealizationAlgebras(ZZ)).run()

    """

    @cached_method
    def super_categories(self):
        return [Algebras(self.base_ring())]

    class ParentMethods:

        @abstract_method
        def expand_to_polynomial(self, elem, k):
            """
            Give a polynomial realization of *elem* over *k* variables.

            EXAMPLES::

                sage: G = FQSym(QQ).G()
                sage: G[3,1,2].expand_to_polynomial(4)
                a2*a1^2 + a3*a1^2 + a3*a1*a2 + a3*a2^2 + a4*a1^2 + a4*a1*a2 + a4*a1*a3 + a4*a2^2 + a4*a2*a3 + a4*a3^2
                sage: (G[1].expand_to_polynomial(2))^2 == (G[1]^2).expand_to_polynomial(2)
                True
            """

    class ElementMethods:
        pass

    class Realizations(RealizationsCategory):
        class ParentMethods:

            def expand_to_polynomial_by_coercion(self, x, k):
                for R in self.realization_of().realizations():
                    self_to_R = R.coerce_map_from(self)
                    R_to_self = self.coerce_map_from(R)
                    prod_meth = R.expand_to_polynomial
                    if self_to_R is not None and R_to_self is not None and \
                       prod_meth != R.expand_to_polynomial_by_coercion:
                        return prod_meth(R(self(x)), k)
                return NotImplementedError

    class WithRealizations(WithRealizationsCategory):
        class ParentMethods:
            pass

    class WithBasis(CategoryWithAxiom_over_base_ring):
        class ParentMethods:

            def expand_to_polynomial_by_linearity(self, Nelt, k):
                return sum(coef * self.expand_to_polynomial_on_basis(mon, k)
                    for (mon, coef) in Nelt.monomial_coefficients().iteritems())

from sage.categories.operators_tools import expand_binary_operator

expand_binary_operator(PolynomialRealizationAlgebras, "expand_to_polynomial")