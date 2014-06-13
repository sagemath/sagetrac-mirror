# -*- coding: utf-8 -*-
"""
The category of classes of combinatorial structures

AUTHORS:

- Jean-Baptiste Priez (2014)
"""
from sage.categories.category import Category
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.misc.abstract_method import abstract_method


class CombinatorialStructures(Category):
    """
    "In mathematics, a combinatorial class is a countable set of mathematical objects."
    (Wikipedia)
    """

    def super_categories(self):
        return [SetsWithGrading()]#, InfiniteEnumeratedSets()] that

    class ParentMethods:

        def product(F, G):
            """
            @param F, G: both are classes of combinatorial structures
            @return: the product of combinatorial structures *FG*

            TESTS::

                sage: B = BinaryTrees()
                sage: C = Compositions()
                sage: CB = C*B; CB
                Product of structures : `Compositions of non-negative integers`, `Binary trees`
                sage: ascii_art(CB.graded_component(2).list())
                [[[1, 1, 1], .],
                 [[1, 2], .],
                 [[2, 1], .],
                 [[3], .],
                 [[1, 1], [., .]],
                 [[2], [., .]],
                 [[1], [., [., .]]],
                 [[1], [[., .], .]],
                 [[], [., [., [., .]]]],
                 [[], [., [[., .], .]]],
                 [[], [[., .], [., .]]],
                 [[], [[., [., .]], .]],
                 [[], [[[., .], .], .]]]
            """
            from sage.combinat.structures.operations.product import CauchyProduct
            return CauchyProduct(F, G)
        _mul_ = product

        def cartesian_product(F, G):
            """
            @param F, G: both are classes of combinatorial structures
            @return: the cartesian product of combinatorial structures *F \square G*

            TESTS::

                sage: B = BinaryTrees()
                sage: C = Compositions()
                sage: BcC = B.cartesian_product(C); BcC
                Cartesian product of structures : 'Binary trees, Compositions of non-negative integers'
                sage: ascii_art(BcC.graded_component(2).list())
                [ [      * ]               [      * ]              ]
                [ [ o    * ]  [ o    ** ]  [   o  * ]  [   o  ** ] ]
                [ [  \     ]  [  \      ]  [  /     ]  [  /      ] ]
                [ [   o,   ], [   o,    ], [ o  ,   ], [ o  ,    ] ]
            """
            from sage.combinat.structures.operations.cartesian_product import CartesianProduct
            return CartesianProduct(F, G)

        def sum(F, G):
            """
            @param F, G: both are classes of combinatorial structures
            @return: the sum of combinatorial structures *F + G*

            TESTS::

                sage: B = BinaryTrees()
                sage: C = Compositions()
                sage: BpC = B + C; BpC
                Sum of structures : `Binary trees`, `Compositions of non-negative integers`
                sage: ascii_art(BpC.graded_component(3).list())
                [                                *              ]
                [                                *  **   *      ]
                [ o      o      o      o      o  *  *   **  *** ]
                [  \      \    / \    /      /                  ]
                [   o      o  o   o  o      o                   ]
                [    \    /           \    /                    ]
                [     o, o  ,      ,   o, o    ,  ,   ,   ,     ]
            """
            from sage.combinat.structures.operations.sum import Sum
            return Sum(F, G)

        __add__ = sum

        def restricted_structures(self, min, max):
            """
            @param min, max: both define an interval of graduation of the class of combinatorial structures
            @return: the restriction of the combinatorial structure to the homogeneous components between *min* and
                     *max*.

            TESTS::

                sage: B = BinaryTrees()
                sage: RB = B.restricted_structures(min=3, max=4); RB
                Binary trees with min=`3`, max=`4`
                sage: RB.graded_component(2).list()
                []
                sage: RB.graded_component(35).list()
                []
                sage: RB.graded_component(3).list()
                [[., [., [., .]]],
                 [., [[., .], .]],
                 [[., .], [., .]],
                 [[., [., .]], .],
                 [[[., .], .], .]]
                sage: _[2].parent()
                Binary trees with min=`3`, max=`4`
            """
            from sage.combinat.structures.operations import RestrictedStructures
            return RestrictedStructures(self, min=min, max=max)

        def __iter__(self):
            """

            """
            for grade in self.grading_set():
                for obj in self.graded_component(grade):
                    yield obj

    class GradedComponents(Category):

        def super_categories(self):
            return [FiniteEnumeratedSets()]

        class ParentMethods:

            def _repr_(self):
                """
                TESTS::

                    sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
                    sage: Compositions(4)
                    Compositions of integers of degree 4
                """
                return repr(self.ambient()) + " of degree " + repr(self.grading())

            @abstract_method(optional=False)
            def ambient(self):
                """
                @return: the combinatorial class of the homogeneous component

                TESTS::

                    sage: BinaryTrees(4).ambient()
                    Binary trees
                """

            @abstract_method(optional=False)
            def grading(self):
                """
                @return: the grading of the homogeneous component

                TESTS::

                    sage: BinaryTrees(4).grading()
                    4
                """