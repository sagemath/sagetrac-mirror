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
from sage.categories.tensor import TensorProductsCategory
from sage.categories.realizations import RealizationsCategory
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.coalgebras import Coalgebras
from sage.categories.coalgebras_with_basis import CoalgebrasWithBasis
from sage.categories.all import ModulesWithBasis, tensor, Hom


class DendriformCoalgebras(Category_over_base_ring):
    """
    The category of dendriform coalgebras over a given base ring.

    An dendriform coalgebra over a ring `R` is an coalgebra with `\delta` is a
    coproduct defined by `\delta := \delta_\prec + \delta_\succ + \bar\delta`
    three operators with `\bar\delta(x) = 1\otimes x + x\otimes 1`
    such that respect following properties, for `x,y,z` three basis element:

    MATH::

        \begin{align}
            ( \delta_\prec \otimes Id ) \circ \delta_\prec( x ) &=
            ( Id \otimes \bar\delta )\circ \delta_\prec(x)\\
            ( \delta_\succ \otimes Id ) \circ \delta_\prec( x ) &=
            ( Id \otimes \delta_\prec )\circ \delta_\succ(x)\\
            ( \bar\delta \otimes Id ) \circ \delta_\succ( y ) &=
            ( Id \otimes \delta_\succ )\circ \delta_\succ(x)\\
        \end{align}

    TODO: should `R` be a commutative ring?

    EXAMPLES::

        sage: DendriformAlgebras(ZZ)
        Category of dendriform algebras over Integer Ring
        sage: Algebras(ZZ).super_categories()
        [Category of non unital algebras over Integer Ring, Category of rings]

    TESTS::

        sage: TestSuite(Algebras(ZZ)).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: Algebras(ZZ).super_categories()
            [Category of non unital algebras over Integer Ring, Category of rings]
        """
        R = self.base_ring()
        return [Coalgebras(R)]

    class ParentMethods:

        def _test_codendriform(self, **options):
            tester = self._tester(**options)
            S = tester.some_elements()
            for e in self.some_elements():
# ( \delta_\prec \otimes Id ) \circ \delta_\prec( x ) &= ( Id \otimes \bar\delta )\circ \delta_\prec(x)\\
                print tensor( (self,self,self) ).sum( [
                        c * tensor( (self(l).left_coproduct(), self(r))) 
                        for ((l,r), c) in self.left_coproduct( e ).monomial_coefficients().iteritems()
                ] )
                tester.assertEqual( 
                    tensor( (self,self,self) ).sum( [
                        c * tensor( (self(l).left_coproduct(), self(r))) 
                        for ((l,r), c) in self.left_coproduct( e ).monomial_coefficients().iteritems()
                    ] ), tensor( (self,self,self) ).sum( [
                        c * tensor( (
                                self(l), 
                                tensor( (
                                    self(r), 
                                    self.one_basis() 
                                ) ) + tensor( (
                                    self.one_basis(),
                                    self(r) 
                                ) ) 
                        ) ) 
                    for ((l,r), c) in self.left_coproduct( e ).monomial_coefficients().iteritems()
                ] ) )
# ( \delta_\succ \otimes Id ) \circ \delta_\prec( x ) &= ( Id \otimes \delta_\prec )\circ \delta_\succ(x)\\
# ( \bar\delta \otimes Id ) \circ \delta_\succ( y ) &= ( Id \otimes \delta_\succ )\circ \delta_\succ(x)\\

        @abstract_method
        def left_coproduct(self, elem):
            '''
            '''

        @abstract_method
        def right_coproduct(self, elem):
            '''
            '''

#         def coproduct(self, elem):
#             return self.tensor_square().sum_of_monomials(
#                 [(self.one_basis(), elem), (elem, self.one_basis())]
#             ) + self.left_coproduct(elem) + \
#                 self.right_coproduct(elem)

    class ElementMethods:

        def right_coproduct(self):
            return self.parent().right_coproduct(self)

        def left_coproduct(self):
            return self.parent().left_coproduct(self)

    class Realizations(RealizationsCategory):

        class ParentMethods:

            def left_coproduct_by_coercion(self, elem):
                r"""
                TODO::
                """
                for R in self.realization_of().realizations():
                    self_to_R = R.coerce_map_from(self)
                    R_to_self = self.coerce_map_from(R)
                    if self_to_R is not None and R_to_self is not None and \
                       R.left_product != R.left_product_by_coercion:
                        return self.tensor_square().sum(
                            coeff * tensor([self(R(I)), self(R(J))])
                            for ((I, J), coeff) in R(self(elem)).coproduct()
                        )
                return NotImplementedError

            def right_coproduct_by_coercion(self, elem):
                r"""
                TODO::
                """
                for R in self.realization_of().realizations():
                    self_to_R = R.coerce_map_from(self)
                    R_to_self = self.coerce_map_from(R)
                    if self_to_R is not None and R_to_self is not None and \
                       R.right_product != R.right_product_by_coercion:
                        return self.tensor_square().sum(
                            coeff * tensor([self(R(I)), self(R(J))])
                            for ((I, J), coeff) in R(self(elem)).\
                                    right_coproduct()
                        )
                return NotImplementedError

    class WithBasis(Category_over_base_ring):

        @cached_method
        def super_categories(self):
            R = self.base_ring()
            return [DendriformCoalgebras(R), CoalgebrasWithBasis(R)]

        def __repr__(self):
            return "Category of dendriform algebras with basis over %s" % (
                    self.base_ring())

        class ParentMethods:

#             @lazy_attribute
#             def coproduct_on_basis(self):
#                 if self.left_coproduct_on_basis is NotImplemented and \
#                    self.right_coproduct_on_basis is NotImplemented:
#                     return self.coproduct_by_coercion
#                 else:
#                     return self.default_coproduct_on_basis
#
#             def default_coproduct_on_basis(self, elem):
#                 DxD = lambda x, y: self.tensor_square()(
#                     tensor((self(x), self(y)))
#                 )
#                 return DxD(self.one_basis(), elem) + \
#                        DxD(elem, self.one_basis()) + \
#                        self.left_coproduct_on_basis(elem) + \
#                        self.right_coproduct_on_basis(elem)

            @abstract_method(optional=True)
            def left_coproduct_on_basis(self, elem):
                """
                TODO::
                """

            @abstract_method(optional=True)
            def right_coproduct_on_basis(self, elem):
                """
                TODO::
                """

            @lazy_attribute
            def left_coproduct(self):
                if self.left_coproduct_on_basis is not NotImplemented:
                    return Hom(
                        self, tensor([self, self]),
                        ModulesWithBasis(self.base_ring())
                    )(
                      on_basis=self.left_coproduct_on_basis
                    )
                elif hasattr(self, "left_coproduct_by_coercion"):
                    return self.left_coproduct_by_coercion

            @lazy_attribute
            def right_coproduct(self):
                if self.right_coproduct_on_basis is not NotImplemented:
                    return Hom(
                        self, tensor([self, self]),
                        ModulesWithBasis(self.base_ring())
                    )(
                      on_basis=self.right_coproduct_on_basis
                    )
                elif hasattr(self, "right_coproduct_by_coercion"):
                    return self.right_coproduct_by_coercion

#             @lazy_attribute
#             def coproduct(self):
#                 if self.coproduct_on_basis is not NotImplemented:
#                     return Hom(
#                         self, tensor([self, self]),
#                         ModulesWithBasis(self.base_ring())
#                     )(
#                       on_basis=self.coproduct_on_basis
#                     )
#                 if self.right_coproduct_on_basis is not NotImplemented and \
#                    self.left_coproduct_on_basis is not NotImplemented:
#                     return Hom(
#                         self, tensor([self, self]),
#                         ModulesWithBasis(self.base_ring())
#                     )(
#                       on_basis=self.default_coproduct_on_basis
#                     )
#                 elif hasattr(self, "coproduct_by_coercion"):
#                     return self.coproduct_by_coercion

    class TensorProducts(TensorProductsCategory):
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Algebras(QQ).TensorProducts().extra_super_categories()
                [Category of algebras over Rational Field]
                sage: Algebras(QQ).TensorProducts().super_categories()
                [Category of algebras over Rational Field]

            Meaning: a tensor product of algebras is an algebra
            """
            return [self.base_category()]

        class ParentMethods:
            #def coproduct(self):
            #    tensor products of morphisms are not yet implemented
            #    return tensor(module.coproduct for module in self.modules)
            pass

        class ElementMethods:
            pass

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
            from dendriform_algebras import DendriformAlgebras
            return [DendriformAlgebras(self.base_category().base_ring())]
