r"""
Complex reflection groups
"""
#*****************************************************************************
#       Copyright (C) 2011-2015 Christian Stump <christian.stump at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.complex_reflection_or_generalized_coxeter_groups import ComplexReflectionOrGeneralizedCoxeterGroups

class ComplexReflectionGroups(Category_singleton):
    r"""
    The category of complex reflection groups.

    Let `V` be a complex vector space. A *complex reflections* is an
    element of `\operatorname{GL}(V)` fixing an hyperplane pointwise
    and acting by multiplication by a root of unity on a complementary
    line.

    A *complex reflection group* is a group `W` that is (isomorphic
    to) a subgroup of some general linear group `\operatorname{GL}(V)`
    generated by a distinguished set of complex reflections.

    The dimension of `V` is the *rank* of `W`.

    For a comprehensive treatment of complex reflection groups and
    many definitions and theorems used here, we refer to [LT2009]_.
    See also :wikipedia:`Reflection_group`.

    .. SEEALSO::

        :func:`ReflectionGroup` for usage examples of this category.

    REFERENCES:

    .. [LT2009] G.I. Lehrer and D.E. Taylor. *Unitary reflection groups*.
       Australian Mathematical Society Lecture Series, 2009.

    EXAMPLES::

        sage: from sage.categories.complex_reflection_groups import ComplexReflectionGroups
        sage: ComplexReflectionGroups()
        Category of complex reflection groups
        sage: ComplexReflectionGroups().super_categories()
        [Category of groups]
        sage: ComplexReflectionGroups().all_super_categories()
        [Category of complex reflection groups,
         Category of groups,
         Category of monoids,
         Category of semigroups,
         Category of inverse unital magmas,
         Category of unital magmas,
         Category of magmas,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

    An example of a reflection group::

        sage: W = ComplexReflectionGroups().example(); W
        5-colored permutations of size 3

    ``W`` is in the category of complex reflection groups::

        sage: W in ComplexReflectionGroups()
        True

    TESTS::

        sage: TestSuite(W).run()
        sage: TestSuite(ComplexReflectionGroups()).run()
    """

    @cached_method
    def super_categories(self):
        r"""
        Return the super categories of ``self``.

        EXAMPLES::

            sage: from sage.categories.complex_reflection_groups import ComplexReflectionGroups
            sage: ComplexReflectionGroups().super_categories()
            [Category of groups]
        """
        return [ComplexReflectionOrGeneralizedCoxeterGroups()]

    def example(self):
        r"""
        Return an example of a complex reflection group.

        EXAMPLES::

            sage: from sage.categories.complex_reflection_groups import ComplexReflectionGroups
            sage: ComplexReflectionGroups().example()
            5-colored permutations of size 3
        """
        from sage.combinat.colored_permutations import ColoredPermutations
        return ColoredPermutations(5, 3)

    class ParentMethods:

        @abstract_method(optional=True)
        def hyperplane_index_set(self):
            r"""
            Return the index set of the reflection hyperplanes of
            ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,4))
                sage: W.hyperplane_index_set()
                (1, 2, 3, 4, 5, 6)
                sage: W = ReflectionGroup((1,1,4), hyperplane_index_set=[1,3,'asdf',7,9,11])
                sage: W.hyperplane_index_set()
                (1, 3, 'asdf', 7, 9, 11)
                sage: W = ReflectionGroup((1,1,4), hyperplane_index_set=('a','b','c','d','e','f'))
                sage: W.hyperplane_index_set()
                ('a', 'b', 'c', 'd', 'e', 'f')
            """

        # TODO TODO This only makes sense in the finite case; even then, one
        # probably wants not to 
        def some_elements_disabled(self):
            r"""
            Return a list of typical elements of ``self``.

            Implements :meth:`Sets.ParentMethods.some_elements` by
            returning some typical element of ``self``.

            EXAMPLES::

                sage: W = ColoredPermutations(1,4)
                sage: W.some_elements()
                [[[0, 0, 0, 0], [2, 1, 3, 4]],
                 [[0, 0, 0, 0], [1, 3, 2, 4]],
                 [[0, 0, 0, 0], [1, 2, 4, 3]],
                 [[0, 0, 0, 0], [1, 2, 3, 4]],
                 [[0, 0, 0, 0], [4, 1, 2, 3]]]
                sage: W.order()
                24
            """
            prod_ref = self.prod(self.reflection(i) for i in self.index_set())
            return list(self.simple_reflections()) + [self.one(), self.an_element(), prod_ref]

        @cached_method
        def rank(self):
            r"""
            Return the rank of ``self``.

            The rank of ``self`` is the dimension of the smallest
            faithfull reflection representation of ``self``.

            EXAMPLES::

                sage: W = CoxeterGroups().example(); W
                sage: W.rank()
                3
            """

    class ElementMethods:

        def apply_distinguished_reflection(self, i, side='right'):
            r"""
            Return the result of the (left/right) multiplication of
            the ``i``-th distingiushed reflection to ``self``.

            INPUT:

            - ``i`` -- an index of a distinguished reflection
            - ``side`` -- (default: ``'right'``) multiplying from left/right

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: W.one().apply_distinguished_reflection(1)
                (1,4)(2,3)(5,6)
                sage: W.one().apply_distinguished_reflection(2)
                (1,3)(2,5)(4,6)
                sage: W.one().apply_distinguished_reflection(3)
                (1,5)(2,4)(3,6)

                sage: W = ReflectionGroup((1,1,3), hyperplane_index_set=['A','B','C']); W
                Irreducible real reflection group of rank 2 and type A2
                sage: W.one().apply_distinguished_reflection('A')
                (1,4)(2,3)(5,6)
                sage: W.one().apply_distinguished_reflection('B')
                (1,3)(2,5)(4,6)
                sage: W.one().apply_distinguished_reflection('C')
                (1,5)(2,4)(3,6)
            """
            G = self.parent()
            if not i in G.hyperplane_index_set():
                raise ValueError("the given index %s is not an index of a hyperplane"%i)
            if side == 'right':
                return self * G.distinguished_reflection(i)
            else:
                return self.parent().reflection(i) * self

        def apply_distinguished_reflections(self, word, side='right'):
            r"""
            Return the result of the (left/right) multiplication of the
            distinguished reflections indexed by the elements in
            ``word`` to ``self``.

            INPUT:

             - ``word`` -- iterable of distinguished reflections indices
             - ``side`` -- (default: ``'right'``) multiplying from left/right

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: W.one().apply_distinguished_reflections([1])
                (1,4)(2,3)(5,6)
                sage: W.one().apply_distinguished_reflections([2])
                (1,3)(2,5)(4,6)
                sage: W.one().apply_distinguished_reflections([2,1])
                (1,2,6)(3,4,5)
            """
            for i in word:
                self = self.apply_distinguished_reflection(i, side=side)
            return self

        def apply_reflection(self, i, side='right'):
            r"""
            Return the result of the (left/right) multiplication of
            the ``i``-th reflection to ``self``.

            INPUT:

             - ``i`` -- an index of a reflection
             - ``side`` -- (default: ``'right'``) multiplying from left/right

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: W.one().apply_reflection(1)
                (1,4)(2,3)(5,6)
                sage: W.one().apply_reflection(2)
                (1,3)(2,5)(4,6)
                sage: W.one().apply_reflection(3)
                (1,5)(2,4)(3,6)

                sage: W = ReflectionGroup((1,1,3), reflection_index_set=['A','B','C']); W
                Irreducible real reflection group of rank 2 and type A2
                sage: W.one().apply_reflection('A')
                (1,4)(2,3)(5,6)
                sage: W.one().apply_reflection('B')
                (1,3)(2,5)(4,6)
                sage: W.one().apply_reflection('C')
                (1,5)(2,4)(3,6)
            """
            W = self.parent()
            if i not in W.reflection_index_set():
                raise ValueError("the given index %s is not an index of a reflection"%i)
            if side == 'right':
                return self * W.reflection(i)
            else:
                return W.reflection(i) * self

        def apply_reflections(self, word, side='right'):
            r"""
            Return the result of the (left/right) multiplication of the
            reflections indexed by the elements in ``word`` to ``self``.

            INPUT:

             - ``word`` -- iterable of reflections indices
             - ``side`` -- (default: ``'right'``) multiplying from left/right

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: W.one().apply_reflections([1])
                (1,4)(2,3)(5,6)
                sage: W.one().apply_reflections([2])
                (1,3)(2,5)(4,6)
                sage: W.one().apply_reflections([2,1])
                (1,2,6)(3,4,5)
            """
            for i in word:
                self = self.apply_reflection(i, side=side)
            return self

    Finite = LazyImport('sage.categories.finite_complex_reflection_groups', 'FiniteComplexReflectionGroups', as_name='Finite')
