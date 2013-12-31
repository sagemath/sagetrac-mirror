# -*- coding: utf-8 -*-
r"""
The #-product algebra category

AUTHORS:

     - Jean-Baptiste Priez (first version)

References:

.. [AvaVien] The product of trees in the Loday-Ronco algebra through Catalan
    alternative tableaux,
    Jean-Christophe Aval and
    Xavier Viennot

.. [AvNoThi] The # product in combinatorial Hopf algebras,
    Jean-Christophe Aval,
    Jean-Christophe Novelli and
    Jean-Yves Thibon

.. [Chapo] Some dendriform functors,
    Frédéric Chapoton
"""
#*****************************************************************************
#  Copyright (C) 2013      Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.modules import Modules
from sage.categories.with_realizations import WithRealizationsCategory
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.realizations import RealizationsCategory


class DieseProductAlgebras(Category_over_base_ring):
    r"""
    Constructs the class of algebras with #-product. This is used to give an
    #-product structure to several combinatorial Hopf algebras.

    TESTS::

        sage: from sage.combinat.hopf_algebras.diese_product_category import \
        ....:   DieseProductAlgebras
        sage: TestSuite(DieseProductAlgebras(ZZ)).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.combinat.hopf_algebras.categories.diese_product_category import \
                    DieseProductAlgebras
            sage: DieseProductAlgebras(ZZ).super_categories()
            [Category of graded modules over Integer Ring]
        """
        R = self.base_ring()
        return [Modules(R).Graded()]

    class ParentMethods:

        @abstract_method(optional=True)
        def diese_product(self, left, right):
            r"""
            #-product as an endomorphism of ``self``.

            This is constructed by extending the method
            :meth:`diese_product_on_basis` bilinearly, if available,
            or :meth:`diese_linear_operator_on_basis` if available, or
            using the method :meth:`diese_product_by_coercion`.

            OUTPUT:

            - The #-product endomorphism of the module elements.

            EXAMPLES::

                sage: F = FQSym(QQ).F()
                sage: G = FQSym(QQ).G()
                sage: F.diese_product == F.diese_product_by_coercion
                True
                sage: F[2,1].diese_product(F[1])
                F[2, 1]
                sage: F[2,1].diese_product(F[1,2])
                F[2, 1, 3] + F[2, 3, 1]
                sage: G.diese_product == G._default_diese_product_from_diese_linear_operator_on_basis
                True
                sage: G[3,1,2].diese_product(G[3,1,2])
                G[5, 1, 4, 2, 3] + G[5, 2, 4, 1, 3] + G[5, 3, 4, 1, 2]

            """

        @abstract_method(optional=True)
        def diese_linear_operator_on_basis(self, k, sigma):
            """
            Let `d_k` be the linear operator on the free associative algebra
            `\mathbb{K}\langle A\rangle` (over some field `\mathbb{K}`) defined
            by

            MATH::

                d_k(w) = \left{\begin{array}{cc}
                    uav & \text{if } w = uaav \text{ for some } a, \text{with}
                    \mid u \mid = k - 1,\\
                    0 & \text{otherwise.}
                \end{array}\right.

            The ``diese_linear_operator_on_basis`` is the projection of `d_k`
            to the graded module.

            INPUT:

                - `k` an indice on `sigma` with `sigma` an element of the basis

            If this method is implemented and :meth:`diese_product_on_basis` is
            not, the diese product is defined as the linearity application of
            this method on the result of current product.
            (`d_k(m_\sigma m_\mu)`).
            """

        def _default_diese_product_from_diese_linear_operator_on_basis(self, x, y):
            res = self.zero()
            for (mon1, coef1) in x.monomial_coefficients().iteritems():
                for (mon2, coef2) in y.monomial_coefficients().iteritems():
                    for (mon, coef) in (self(mon1) * self(mon2)).\
                            monomial_coefficients().iteritems():
                        res += coef1 * coef2 * coef * self.\
                            diese_linear_operator_on_basis(mon1.size(), mon)
            return res

        def _default_diese_product_from_diese_linear_operator_on_basis(self, x, y):
            res = self.zero()
            for (mon1, coef1) in x.monomial_coefficients().iteritems():
                for (mon2, coef2) in y.monomial_coefficients().iteritems():
                    for (mon, coef) in (self(mon1) * self(mon2)).\
                            monomial_coefficients().iteritems():
                        res += coef1 * coef2 * coef * self.\
                            diese_linear_operator_on_basis(mon1.size(), mon)
            return res

    class ElementMethods:
        pass

    class Realizations(RealizationsCategory):
        class ParentMethods:
            pass

    class WithRealizations(WithRealizationsCategory):
        class ParentMethods:
            pass

    class WithBasis(CategoryWithAxiom_over_base_ring):
        class ParentMethods:

            @lazy_attribute
            def diese_product(self):
                r"""
                #-product as an endomorphism of ``self``.

                This is constructed by extending the method
                :meth:`diese_product_on_basis` bilinearly, if available,
                or :meth:`diese_linear_operator_on_basis` if available, or
                using the method :meth:`diese_product_by_coercion`.

                OUTPUT:

                - The #-product endomorphism of the module elements.

                EXAMPLES::

                    sage: F = FQSym(QQ).F()
                    sage: G = FQSym(QQ).G()
                    sage: F.diese_product == F.diese_product_by_coercion
                    True
                    sage: F[2,1].diese_product(F[1])
                    F[2, 1]
                    sage: F[2,1].diese_product(F[1,2])
                    F[2, 1, 3] + F[2, 3, 1]
                    sage: G.diese_product == G._default_diese_product_from_diese_linear_operator_on_basis
                    True
                    sage: G[3,1,2].diese_product(G[3,1,2])
                    G[5, 1, 4, 2, 3] + G[5, 2, 4, 1, 3] + G[5, 3, 4, 1, 2]

                """
                if self.diese_product_on_basis is NotImplemented:
                    if self.diese_linear_operator_on_basis is NotImplemented:
                        return self.diese_product_by_coercion
                    return self.\
                        _default_diese_product_from_diese_linear_operator_on_basis
                else:
                    return self.module_morphism(
                    self.module_morphism(self.diese_product_on_basis,
                         position=0,
                         codomain=self),
                    position=1)

from sage.categories.operators_tools import expand_binary_operator

expand_binary_operator(DieseProductAlgebras, "diese_product")