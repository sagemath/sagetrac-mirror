r"""
Additive semigroups
"""
#*****************************************************************************
#  Copyright (C) 2013-2014 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

import operator
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom_singleton
from sage.categories.algebra_functor import AlgebrasCategory
from sage.categories.additive_magmas import AdditiveMagmas
from sage.categories.action import Action
from sage.structure.element import generic_power

class AdditiveSemigroups(CategoryWithAxiom_singleton):
    """
    The category of additive semigroups.

    An *additive semigroup* is an associative :class:`additive magma
    <AdditiveMagmas>`, that is a set endowed with an operation `+`
    which is associative.

    EXAMPLES::

        sage: from sage.categories.additive_semigroups import AdditiveSemigroups
        sage: C = AdditiveSemigroups(); C
        Category of additive semigroups
        sage: C.super_categories()
        [Category of additive magmas]
        sage: C.all_super_categories()
        [Category of additive semigroups,
         Category of additive magmas,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

        sage: C.axioms()
        frozenset(['AdditiveAssociative'])
        sage: C is AdditiveMagmas().AdditiveAssociative()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (AdditiveMagmas, "AdditiveAssociative")

    AdditiveCommutative = LazyImport('sage.categories.commutative_additive_semigroups', 'CommutativeAdditiveSemigroups', at_startup=True)
    AdditiveUnital = LazyImport('sage.categories.additive_monoids', 'AdditiveMonoids', at_startup=True)

    class ParentMethods:
        def _test_additive_associativity(self, **options):
            r"""
            Test associativity for (not necessarily all) elements of this
            additive semigroup.

            INPUT:

            - ``options`` -- any keyword arguments accepted by :meth:`_tester`

            EXAMPLES:

            By default, this method tests only the elements returned by
            ``self.some_elements()``::

                sage: S = CommutativeAdditiveSemigroups().example()
                sage: S._test_additive_associativity()

            However, the elements tested can be customized with the
            ``elements`` keyword argument::

                sage: (a,b,c,d) = S.additive_semigroup_generators()
                sage: S._test_additive_associativity(elements = (a, b+c, d))

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)
            S = tester.some_elements()
            from sage.combinat.cartesian_product import CartesianProduct
            for x,y,z in tester.some_elements(CartesianProduct(S,S,S)):
                tester.assert_((x + y) + z == x + (y + z))

        def __init_extra__(self):
            """
            Register the (partial) action of `\ZZ` by multiplication on the left and on the right.

            .. WARNING:: This is actually only an action of `\NN`.

            EXAMPLES::

                sage: E = CommutativeAdditiveMonoids().example()
                sage: e = E.an_element(); e
                a + 3*c + 2*b + 4*d
                sage: 2*e
                2*a + 6*c + 4*b + 8*d
                sage: e*3
                3*a + 9*c + 6*b + 12*d
            """
            from sage.categories.modules import Modules
            if self in Modules:
                 # In this case multiplication by integers is already
                 # defined via coercion into the base ring. Since this
                 # is likely to be faster, we don't redefine it here.
                 return
            # TODO: it ought to be possible to define such an action
            # from just a method (or method name), without having to
            # create a class just for that.
            from sage.rings.integer_ring import ZZ
            left_action  = IntMultAction(self, ZZ, is_left=1, op=operator.mul)
            right_action = IntMultAction(self, ZZ, is_left=0, op=operator.mul)
            self.register_action(left_action)
            self.register_action(right_action)

    class ElementMethods:
        def __mul__(left, right):
            """
            Return ``left*right``.

            The calculation is delegated to the coercion model.

            EXAMPLES::

                sage: E = CommutativeAdditiveMonoids().example()
                sage: e = E.an_element(); e
                a + 3*c + 2*b + 4*d
                sage: e * 3
                3*a + 9*c + 6*b + 12*d

            TESTS::

                sage: F = CombinatorialFreeModule(QQ, ["a", "b"])
                sage: x = F.monomial("a")
                sage: x * int(2)
                2*B['a']
            """
            from sage.structure.element import get_coercion_model
            import operator
            return get_coercion_model().bin_op(left, right, operator.mul)

        def __rmul__(right, left):
            """
            Return ``left*right`` with ``right`` in this class.

            The calculation is delegated to the coercion model.

            EXAMPLES::

                sage: E = CommutativeAdditiveMonoids().example()
                sage: e = E.an_element(); e
                a + 3*c + 2*b + 4*d
                sage: 3 * e
                3*a + 9*c + 6*b + 12*d

            TESTS::

                sage: F = CombinatorialFreeModule(QQ, ["a", "b"])
                sage: x = F.monomial("a")
                sage: int(2) * x
                2*B['a']

            .. TODO::

                Add an example where multiplication on the left and on
                the right differ.
            """
            from sage.structure.element import get_coercion_model
            import operator
            return get_coercion_model().bin_op(left, right, operator.mul)

        def _intmul_(self, n):
            """
            Return ``self`` multiplied by ``n``.

            INPUT:

            - ``n`` -- a positive integer

            EXAMPLES::

                sage: E = CommutativeAdditiveMonoids().example()
                sage: e = E.an_element(); e
                a + 3*c + 2*b + 4*d
                sage: e._intmul_(3)
                3*a + 9*c + 6*b + 12*d
                sage: e._intmul_(1)
                a + 3*c + 2*b + 4*d

                sage: e._intmul_(-1)
                Traceback (most recent call last):
                ...
                ValueError: n should be positive

            Multiplication by `0` could be implemented in
            :class:`.additive_monoids.AdditiveMonoids`::

                sage: e._intmul_(0)
                Traceback (most recent call last):
                ...
                ValueError: n should be positive
                sage: e._intmul_(0)  # todo: not implemented
                0
            """
            n = int(n)
            if n <= 0:
                 raise ValueError("n should be positive")
            # This should use binary powering and share the code with generic_power
            return sum([self for i in range(n-1)], self)

    class Algebras(AlgebrasCategory):

        def extra_super_categories(self):
            """
            Implement the fact that the algebra of a semigroup is an
            associative (but not necessarily unital) algebra.

            EXAMPLES::

                sage: from sage.categories.additive_semigroups import AdditiveSemigroups
                sage: AdditiveSemigroups().Algebras(QQ).extra_super_categories()
                [Category of semigroups]
                sage: CommutativeAdditiveSemigroups().Algebras(QQ).super_categories()
                [Category of additive semigroup algebras over Rational Field,
                 Category of additive commutative additive magma algebras over Rational Field]
            """
            from sage.categories.semigroups import Semigroups
            return [Semigroups()]

        class ParentMethods:

            @cached_method
            def algebra_generators(self):
                r"""
                Return the generators of this algebra, as per
                :meth:`MagmaticAlgebras.ParentMethods.algebra_generators()
                <.magmatic_algebras.MagmaticAlgebras.ParentMethods.algebra_generators>`.

                They correspond to the generators of the additive semigroup.

                EXAMPLES::

                    sage: S = CommutativeAdditiveSemigroups().example(); S
                    An example of a commutative monoid: the free commutative monoid generated by ('a', 'b', 'c', 'd')
                    sage: A = S.algebra(QQ)
                    sage: A.algebra_generators()
                    Finite family {0: B[a], 1: B[b], 2: B[c], 3: B[d]}
                """
                return self.basis().keys().additive_semigroup_generators().map(self.monomial)

            def product_on_basis(self, g1, g2):
                r"""
                Product, on basis elements, as per
                :meth:`MagmaticAlgebras.WithBasis.ParentMethods.product_on_basis()
                <sage.categories.magmatic_algebras.MagmaticAlgebras.WithBasis.ParentMethods.product_on_basis>`.

                The product of two basis elements is induced by the
                addition of the corresponding elements of the group.

                EXAMPLES::

                    sage: S = CommutativeAdditiveSemigroups().example(); S
                    An example of a commutative monoid: the free commutative monoid generated by ('a', 'b', 'c', 'd')
                    sage: A = S.algebra(QQ)
                    sage: a,b,c,d = A.algebra_generators()
                    sage: a * b + b * d * c
                    B[c + b + d] + B[a + b]
                """
                return self.monomial(g1 + g2)

class IntMultAction(Action):
     r"""
     Action of integers on an additive semigroup by multiplication.
     """

     def _call_(self, x,  n):
         r"""
         Return ``self*n`` or ``n*self``.

         INPUT::

         - ``n`` -- an integer

         This assumes that ``self*n == n*self``, and calls
         ``self._int_mul_(n)``.

         EXAMPLES::

             sage: from sage.categories.additive_semigroups import IntMultAction
             sage: E = CommutativeAdditiveMonoids().example()
             sage: left_action  = IntMultAction(E, ZZ, is_left=True)
             sage: right_action = IntMultAction(E, ZZ, is_left=False)
             sage: e = E.an_element(); e
             a + 3*c + 2*b + 4*d
             sage: left_action._call_(e, 2)
             2*a + 6*c + 4*b + 8*d
             sage: right_action._call_(3, e)
             3*a + 9*c + 6*b + 12*d
         """
         if not self.is_left():
             x,n = n,x
         return x._intmul_(n)

