# -*- coding: utf-8 -*-
r"""
 Category of algebras with scalar product

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


class ScalarProductAlgebras(Category_over_base_ring):

    @cached_method
    def super_categories(self):
        return [Algebras(self.base_ring())]

    class ParentMethods:

        @abstract_method
        def scalar_product(self, left, right):
            """
            Scalar product of <left, right>.

            EXAMPLES::

                sage: F = FQSym(QQ).F()
                sage: sigma = Permutation([1,4,2,3])
                sage: F.scalar_product_on_basis(sigma, sigma.inverse())
                1
                sage: F.scalar_product_on_basis(sigma, sigma)
                0
            """

    class ElementMethods:
        pass

    class Realizations(RealizationsCategory):
        class ParentMethods:

            def scalar_product_by_coercion(self, left, right):
                for R in self.realization_of().realizations():
                    self_to_R = R.coerce_map_from(self)
                    R_to_self = self.coerce_map_from(R)
                    prod_meth = R.scalar_product
                    if self_to_R is not None and R_to_self is not None and \
                       prod_meth != R.scalar_product_by_coercion:
                        return prod_meth(R(left), R(right))
                return NotImplementedError

    class WithRealizations(WithRealizationsCategory):
        class ParentMethods:
            pass

    class WithBasis(CategoryWithAxiom_over_base_ring):
        class ParentMethods:

            def scalar_product_by_linearity(self, left, right):
                if left.parent() != right.parent():
                    return self.scalar_product_by_coercion(left, right)
                
                return sum(coefL * coefR * self.scalar_product_on_basis(monL, monR)
                    for (monL, coefL) in left.monomial_coefficients().iteritems()
                    for (monR, coefR) in right.monomial_coefficients().iteritems())

from sage.categories.operators_tools import expand_binary_operator

expand_binary_operator(ScalarProductAlgebras, "scalar_product")