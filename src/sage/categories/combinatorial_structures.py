# -*- coding: utf-8 -*-
"""
The category of classes of combinatorial structures

AUTHOR:

- Jean-Baptiste Priez (2014)
"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.category import Category
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.misc.abstract_method import abstract_method


class CombinatorialStructures(Category):
    """
    "In mathematics, a combinatorial class is a countable set of mathematical
    objects." (Wikipedia)
    """

    def super_categories(self):
        return [SetsWithGrading(), InfiniteEnumeratedSets()]

    class ParentMethods:

        def cycle(F, **options):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: BC = B.cycle(); BC
                Cycle of `Binary trees`
                sage: ascii_art(BC.graded_component(3,2).list())
                [                         (         )  (       )               (       )  (
                [                         ( , o     )  ( , o   )               ( ,   o )  ( ,
                [ (        )  (        )  (    \    )  (    \  )  (         )  (    /  )  (
                [ ( o, o   )  ( o,   o )  (     o   )  (     o )  ( ,   o   )  (   o   )  (
                [ (     \  )  (     /  )  (      \  )  (    /  )  (    / \  )  (    \  )  (
                [ (      o ), (    o   ), (       o ), (   o   ), (   o   o ), (     o ), (
                <BLANKLINE>
                      ) ]
                    o ) ]
                   /  ) ]
                  o   ) ]
                 /    ) ]
                o     ) ]

                sage: BC2 = B.restricted_structures(min=1).cycle(grading_set="sum"); BC2
                Cycle of `Binary trees with grading min=`1``
                sage: ascii_art(BC2.graded_component(3).list())
                [ ( o     )  ( o   )             (   o )  (     o )
                [ (  \    )  (  \  )             (  /  )  (    /  )  (        )  (        )
                [ (   o   )  (   o )  (   o   )  ( o   )  (   o   )  ( o, o   )  ( o,   o )
                [ (    \  )  (  /  )  (  / \  )  (  \  )  (  /    )  (     \  )  (     /  )  (
                [ (     o ), ( o   ), ( o   o ), (   o ), ( o     ), (      o ), (    o   ), (
                <BLANKLINE>
                          ]
                          ]
                          ]
                        ) ]
                o, o, o ) ]
            """
            from sage.combinat.structures.operations.cycle import Cycle
            return Cycle(F, **options)

        def multiset(F, **options):
            """
            TESTS::

                sage: C = Compositions()
                sage: MC = C.multiset(); MC
                Multi-set of `Compositions of non-negative integers`
                sage: MC32 = MC.graded_component(3,2); MC32
                Multi-set of `Compositions of non-negative integers` of degree
                (3, 2)
                sage: MC32.list()
                [{[1], [1, 1]}, {[1], [2]}, {[], [1, 1, 1]}, {[], [1, 2]}, {[],
                [2, 1]}, {[], [3]}]
            """
            from sage.combinat.structures.operations.multi_sets import MultiSet
            return MultiSet(F, **options)

        def sequence(F, **options):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: BS = B.sequence(); BS
                Sequence of `Binary trees`
                sage: ascii_art(BS.graded_component(3,2).list())
                [
                [ [ o       ]  [ o     ]               [   o   ]  [     o   ]
                [ [  \      ]  [  \    ]               [  /    ]  [    /    ]
                [ [   o     ]  [   o   ]  [   o     ]  [ o     ]  [   o     ]  [ o    o ]  [
                [ [    \    ]  [  /    ]  [  / \    ]  [  \    ]  [  /      ]  [  \     ]  [
                [ [     o,  ], [ o  ,  ], [ o   o,  ], [   o,  ], [ o    ,  ], [   o,   ], [
                <BLANKLINE>
                                                  [         ]  [       ]               [
                                                  [ , o     ]  [ , o   ]               [ ,   o
                          [        ]  [        ]  [    \    ]  [    \  ]  [         ]  [    /
                  o  o ]  [ o, o   ]  [ o,   o ]  [     o   ]  [     o ]  [ ,   o   ]  [   o
                 /     ]  [     \  ]  [     /  ]  [      \  ]  [    /  ]  [    / \  ]  [    \
                o  ,   ], [      o ], [    o   ], [       o ], [   o   ], [   o   o ], [     o
                <BLANKLINE>
                 ]  [         ] ]
                 ]  [ ,     o ] ]
                 ]  [      /  ] ]
                 ]  [     o   ] ]
                 ]  [    /    ] ]
                 ], [   o     ] ]
                sage: BS2 = B.restricted_structures(min=1).sequence(grading_set="sum"); BS2
                Sequence of `Binary trees with grading min=`1``
                sage: ascii_art(BS2.graded_component(3).list())
                [ [ o     ]  [ o   ]             [   o ]  [     o ]
                [ [  \    ]  [  \  ]             [  /  ]  [    /  ]                          [
                [ [   o   ]  [   o ]  [   o   ]  [ o   ]  [   o   ]  [ o    o ]  [   o  o ]  [
                [ [    \  ]  [  /  ]  [  / \  ]  [  \  ]  [  /    ]  [  \     ]  [  /     ]  [
                [ [     o ], [ o   ], [ o   o ], [   o ], [ o     ], [   o,   ], [ o  ,   ], [
                <BLANKLINE>
                                                  ]
                       ]  [        ]              ]
                o, o   ]  [ o,   o ]              ]
                    \  ]  [     /  ]  [         ] ]
                     o ], [    o   ], [ o, o, o ] ]

            """
            from sage.combinat.structures.operations.sequences import Sequence
            return Sequence(F, **options)

        def product(F, G):
            """
            Return the Cauchy product of two structure.

            INPUT:

            -``F``, ``G`` -- two combinatorial structures.

            OUTPUT:

            The Cauchy product of ``F`` and ``G``.

            TESTS::

                sage: B = BinaryTrees()
                sage: C = Compositions()
                sage: CB = C*B; CB
                Product of structures : `Compositions of non-negative integers`,
                `Binary trees`
                sage: ascii_art(CB.graded_component(2).list())
                [                              [       ]  [       ] ]
                [                              [ , o   ]  [ ,   o ] ]
                [ [ *   ]  [      ]  [      ]  [    \  ]  [    /  ] ]
                [ [ *,  ], [ **,  ], [ *, o ], [     o ], [   o   ] ]
            """
            from sage.combinat.structures.operations.product import CauchyProduct
            return CauchyProduct(F, G)
        _mul_ = product

        def cartesian_product(F, G):
            """
            Return the cartesian products of two structures.

            INPUT:

            -``F``, ``G`` -- two combinatorial structures.

            OUTPUT:

            The cartesian product of ``F`` and ``G``.

            TESTS::

                sage: B = BinaryTrees()
                sage: C = Compositions()
                sage: BcC = B.cartesian_product(C); BcC
                Cartesian product of structures : 'Binary trees, Compositions of
                non-negative integers'
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
            Return the sum of two structures.

            INPUT:

            -``F``, ``G`` -- two combinatorial structures.

            OUTPUT:

            The sum of ``F`` and ``G``.

            TESTS::

                sage: B = BinaryTrees()
                sage: C = Compositions()
                sage: BpC = B + C; BpC
                Sum of structures : `Binary trees`, `Compositions of
                non-negative integers`
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

        def restricted_structures(self, min=None, max=None):
            """
            Return the restriction of the combinatorial structure to the
            homogeneous component between ``min`` and ``max``.

            Return the sum of two structures.

            INPUT:

            -``min`` -- integer (default: None) lower bound of the interval.
            -``max`` -- integer (default: None) upper bound of the interval.

            OUTPUT:

            The restricted structure of ``self`` on the interval bounded by
            ``min`` and ``max``.

            TESTS::

                sage: B = BinaryTrees()
                sage: RB = B.restricted_structures(min=3, max=4); RB
                Binary trees with grading min=`3`, max=`4`
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
                Binary trees with grading min=`3`, max=`4`
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
                Return the combinatorial class of the homogeneous component.

                TESTS::

                    sage: BinaryTrees(4).ambient()
                    Binary trees
                """

            @abstract_method(optional=False)
            def grading(self):
                """
                Return the grading of the homogenous component.

                TESTS::

                    sage: BinaryTrees(4).grading()
                    4
                """
