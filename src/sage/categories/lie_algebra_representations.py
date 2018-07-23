r"""
Lie Algebra representations

AUTHORS:

- Michael Walter (2018-07-22): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Michael Walter <m.walter@uva.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import absolute_import

from sage.categories.category_types import Category_module
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.modules import Modules
from sage.categories.tensor import tensor, TensorProductsCategory
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method


class LieAlgebraRepresentations(Category_module):
    r"""
    The category of representations of Lie algebras over a field ``base``.

    INPUT:

    - ``base`` -- the base field

    EXAMPLES::

        sage: LieAlgebraRepresentations(QQ)
        Category of lie algebra representations over Rational Field

    TESTS::

        sage: TestSuite(LieAlgebraRepresentations(QQ)).run()
        sage: TestSuite(LieAlgebraRepresentations(CC)).run()
    """

    @cached_method
    def super_categories(self):
        """
        Return the super categories of ``self``.

        EXAMPLES::

            sage: LieAlgebraRepresentations(CC).super_categories()
            [Category of vector spaces over Complex Field with 53 bits of precision]
        """
        return [Modules(self.base_ring())]

    def example(self):
        """
        Return an example of a Lie algebra representation.

        EXAMPLES::

            sage: LieAlgebraRepresentations(QQ).example()
            Lie algebra representation of ['A', 2] with highest weight (2, 1, 0)
        """
        from sage.combinat.root_system.root_system import RootSystem
        from sage.algebras.lie_algebras.representations import (
            lie_algebra_representation
        )

        return lie_algebra_representation([2, 1, 0], "A2", self.base_ring())

    class ParentMethods:
        @abstract_method
        def cartan_type(self):
            r"""
            Return the Cartan type of ``self``.

            EXAMPLES::

                sage: V = lie_algebra_representation(Partition([2, 1]), "A3")
                sage: V.cartan_type()
                ['A', 3]
            """

        @abstract_method
        def action(self, X):
            r"""
            Return the action of a Lie algebra element.

            INPUT:

            - ``X`` -- a Lie algebra element 

            OUTPUT:

            An endomorphism of ``self``.

            EXAMPLES::
            
                sage: V = lie_algebra_representation([2, 0], "A1")
                sage: V.action(matrix([[1, 0], [0, 0]])).matrix()
                [0 0 0]
                [0 1 0]
                [0 0 2]
            """

    class TensorProducts(TensorProductsCategory):
        """
        The category of Lie algebra representations constructed by tensor product of Lie algebra representations.
        """

        @cached_method
        def extra_super_categories(self):
            r"""
            EXAMPLES::

                sage: LieAlgebraRepresentations(QQ).TensorProducts().extra_super_categories()
                [Category of lie algebra representations over Rational Field]
            """
            return [self.base_category()]

        class ParentMethods:
            @cached_method
            def cartan_type(self):
                r"""
                Return the Cartan type of ``self``.

                EXAMPLES::

                    sage: V = lie_algebra_representation([2, 1], "A1")
                    sage: W = lie_algebra_representation([7, 3], "A1")
                    sage: T = tensor([V, W])
                    sage: T.cartan_type()
                    ['A', 1]
                """
                return self._sets[0].cartan_type()

    class WithBasis(CategoryWithAxiom_over_base_ring):
        """
        The category of Lie algebra representations with a distinguished basis.
        """

        class ParentMethods:
            @cached_method
            def tensor(*factors):
                r"""
                Return the tensor product of the given ``factors``.

                If $V_1, \dots, V_n$ are representations of a Lie algebra $\mathfrak g$
                then the tensor product $V_1 \otimes \dots \otimes V_n$ is naturally
                a representation of $\mathfrak g$. The action of an element $X$ is given
                by

                .. MATH::

                    X_1 \otimes \id \otimes \dots \otimes \id + \dots +
                    \id \otimes \dots \otimes \id \otimes X_n,

                where $X_k$ denotes the action of $X$ on $V_k$.

                INPUT:

                - ``factors`` -- a list (or other iterable) of factors

                OUTPUT:

                The Lie algebra representation on the tensor product.

                EXAMPLES::

                    sage: V = lie_algebra_representation([2, 1], "A1")
                    sage: W = lie_algebra_representation([7, 3], "A1")
                    sage: T = tensor([V, W])
                    sage: TestSuite(T).run()

                    sage: X = lie_algebra_representation([7, 3, 1], "A2")
                    sage: tensor([V, X])
                    Traceback (most recent call last):
                    ...
                    ValueError: all factors must have same Cartan type
                """
                cartan_type = factors[0].cartan_type()
                if any(V.cartan_type() != cartan_type for V in factors):
                    raise ValueError("all factors must have same Cartan type")

                cat = tensor.category_from_parents(factors)
                return factors[0].__class__.Tensor(factors, category=cat)

        class TensorProducts(TensorProductsCategory):
            """
            The category of Lie algebra representations constructed by tensor product of Lie algebra representations with a distinguished basis.
            """

            @cached_method
            def extra_super_categories(self):
                r"""
                EXAMPLES::

                    sage: LieAlgebraRepresentations(QQ).WithBasis().TensorProducts().extra_super_categories()
                    [Category of lie algebra representations with basis over Rational Field]
                """
                return [self.base_category()]

            class ParentMethods:
                def action(self, X):
                    r"""
                    Return the action of a Lie algebra element.

                    INPUT:

                    - ``X`` -- a Lie algebra element 

                    OUTPUT:

                    An endomorphism of ``self``.

                    EXAMPLES::

                        sage: V = lie_algebra_representation([1, 0], "A1")
                        sage: T = tensor([V, V])
                        sage: T.action(matrix([[1, 0], [0, 1]])).matrix()
                        [2 0 0 0]
                        [0 2 0 0]
                        [0 0 2 0]
                        [0 0 0 2]
                        sage: T.action(matrix([[1, 0], [0, -1]])).matrix()
                        [-2  0  0  0]
                        [ 0  0  0  0]
                        [ 0  0  0  0]
                        [ 0  0  0  2]
                        sage: T.action(matrix([[0, 1], [0, 0]])).matrix()
                        [0 0 0 0]
                        [1 0 0 0]
                        [1 0 0 0]
                        [0 1 1 0]
                        sage: T.action(matrix([[0, 0], [1, 0]])).matrix()
                        [0 1 1 0]
                        [0 0 0 1]
                        [0 0 0 1]
                        [0 0 0 0]

                        sage: V = lie_algebra_representation([3, 1], "A1")
                        sage: W = lie_algebra_representation([5, 2], "A1")
                        sage: T = tensor([V, W])
                        sage: T.action(matrix([[1, 0], [0, -1]])).matrix()
                        [-5  0  0  0  0  0  0  0  0  0  0  0]
                        [ 0 -3  0  0  0  0  0  0  0  0  0  0]
                        [ 0  0 -1  0  0  0  0  0  0  0  0  0]
                        [ 0  0  0  1  0  0  0  0  0  0  0  0]
                        [ 0  0  0  0 -3  0  0  0  0  0  0  0]
                        [ 0  0  0  0  0 -1  0  0  0  0  0  0]
                        [ 0  0  0  0  0  0  1  0  0  0  0  0]
                        [ 0  0  0  0  0  0  0  3  0  0  0  0]
                        [ 0  0  0  0  0  0  0  0 -1  0  0  0]
                        [ 0  0  0  0  0  0  0  0  0  1  0  0]
                        [ 0  0  0  0  0  0  0  0  0  0  3  0]
                        [ 0  0  0  0  0  0  0  0  0  0  0  5]
                    """

                    tc = self.tensor_constructor(self._sets)
                    actions = [V.action(X) for V in self._sets]

                    def on_basis(bs):
                        # XXX: there must be a better way
                        result = self.zero()
                        for k in range(len(self._sets)):
                            x = [V(b) for (V, b) in zip(self._sets, bs)]
                            x[k] = actions[k](x[k])
                            result += tc(*x)
                        return result

                    return self.module_morphism(on_basis=on_basis, codomain=self)
