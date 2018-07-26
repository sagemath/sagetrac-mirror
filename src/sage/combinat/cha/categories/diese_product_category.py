# -*- coding: utf-8 -*-
r"""
The #-product algebra category

AUTHORS:

     - Jean-Baptiste Priez (first version)

References:

.. [AvaVien] The product of trees in the Loday-Ronco algebra through Catalan
    alternative tableaux,
    Jean-Christophe Aval and
    Xavier Viennot,

.. [AvNoThi] The # product in combinatorial Hopf algebras,
    Jean-Christophe Aval,
    Jean-Christophe Novelli and
    Jean-Yves Thibon

.. [Chapo] Some dendriform functors,
    Frédéric Chapoton
"""
#*****************************************************************************
#  Copyright (C) 2013      Rémi Maurice <maurice@univ-mlv.fr>
#                          Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.category_types import Category_over_base_ring
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.realizations import RealizationsCategory


class GradedAlgebrasWithDieseProduct(Category_over_base_ring):
    r"""
    Constructs the class of algebras with #-product. This is used to give an
    #-product structure to several combinatorial Hopf algebras.

    EXAMPLES::

        sage: from sage.combinat.ncsf_qsym.generic_basis_code import\
        ....:     GradedModulesWithInternalProduct
        sage: N = NonCommutativeSymmetricFunctions(QQ)
        sage: R = N.ribbon()
        sage: R in GradedModulesWithInternalProduct(QQ)
        True
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.combinat.ncsf_qsym.generic_basis_code import\
            ....:     GradedModulesWithInternalProduct
            sage: GradedModulesWithInternalProduct(ZZ).super_categories()
            [Category of graded modules over Integer Ring]
        """
        from sage.categories.graded_modules import GradedModules
        R = self.base_ring()
        return [GradedModules(R)]

    class ParentMethods:

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

        @abstract_method(optional=True)
        def diese_product_on_basis(self, sigma, mu):
            """
            The diese product of the two basis elements indexed by ``sigma`` and
            ``mu``.

            INPUT:

                - ``sigma``, ``mu`` -- two elements of the basis of self

            Returns the diese product of the corresponding basis elements.
            If this method is implemented, the diese product is defined from
            it by linearity.

            TESTS::

                sage: G = FQSym(QQ).G()
                sage: G.diese_product == G._default_diese_product_from_diese_linear_operator_on_basis
                True
                sage: G[3,1,2].diese_product(G[3,1,2])
                G[5, 1, 4, 2, 3] + G[5, 2, 4, 1, 3] + G[5, 3, 4, 1, 2]
            """

        def _default_diese_product_from_diese_linear_operator_on_basis(
                    self, x, y):
            res = self.zero()
            for (mon1, coef1) in x.monomial_coefficients().iteritems():
                for (mon2, coef2) in y.monomial_coefficients().iteritems():
                    for (mon, coef) in (self(mon1) * self(mon2)).\
                            monomial_coefficients().iteritems():
                        res += coef1 * coef2 * coef * self.\
                            diese_linear_operator_on_basis(mon1.size(), mon)
            return res


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

    class ElementMethods:
        def diese_product(self, other):
            r"""
            Returns the #-product of two elements of the module.

            INPUT:

            - ``other`` -- another element of the module

            OUTPUT:

            - Returns the result of taking the #-product of ``self`` with
              ``other``.

            EXAMPLES::

                sage: G = FQSym(QQ).G()
                sage: F = FQSym(QQ).F()
                sage: F[2,1].diese_product(G[2,1])
                F[3, 2, 1]
            """
            return self.parent().diese_product(self, other)

    class Realizations(RealizationsCategory):
        class ParentMethods:
            def diese_product_by_coercion(self, left, right):
                r"""
                #-product of ``left`` and ``right``.

                This is a default implementation that try to compute
                the #-product in an other realization specified which
                implement :meth:``diese_product``.

                INPUT:

                - ``left`` and ``right`` -- two elements of the module

                OUTPUT:

                - The #-product of ``left`` and ``right``.

                TESTS::

                    sage: E = FQSym(QQ).E()
                    sage: E.diese_product == E.diese_product_by_coercion
                    True
                    sage: E[1,3,2].diese_product(E[2,3,1])
                    E[1, 4, 5, 3, 2]
                """
                for R in self.realization_of().realizations():
                    self_to_R = R.coerce_map_from(self)
                    R_to_self = self.coerce_map_from(R)
                    if self_to_R is not None and R_to_self is not None and \
                       R.diese_product != R.diese_product_by_coercion:
                        return self(R(left).diese_product(R(right)))
                return NotImplementedError
