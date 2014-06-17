# -*- coding: utf-8 -*-
"""
Sequence class of combinatorial structures.

Let `F` be class of combinatorial structures, the *sequence* class `\mathbf{Seq}(F)` is
defined the infinite sum:

MATH::

    \mathbf{Seq}(F) := \{\epsilon\} + F + (F\times F) + (F\times F \times F) + \cdots

with `\epsilon` being a neutral element _[FS].


References:
-----------

.. [FS] Analytic combinatorics
  Philippe Flajolet and Robert Sedgewick

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
from itertools import product, chain
from sage.categories.combinatorial_structures import CombinatorialStructures
from sage.combinat.integer_vector import IntegerVectors
from sage.misc.ascii_art import ascii_art_list
from sage.misc.misc_c import prod
from sage.rings.integer import Integer
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.list_clone import ClonableArray
from sage.combinat.structures import Structures, Structure
from sage.combinat.structures.operations import _Operations


class Sequence(_Operations):
    """
    Sequence class of combinatorial structures.

    Let `F` be class of combinatorial structures, the *sequence* class `\mathbf{Seq}(F)` is
    defined the infinite sum:

    MATH::

        \mathbf{Seq}(F) := \{\epsilon\} + F + (F\times F) + (F\times F \times F) + \cdots

    with `\epsilon` being a neutral element _[FS].

    TESTS::

        sage: B = BinaryTrees()
        sage: BSeq = B.sequence(); BSeq
        Sequence of `Binary trees`
        sage: BSeq23 = BSeq.graded_component((2,3))
        sage: ascii_art(BSeq23.list())
        [                                                   [         ]  [         ]
        [ [ o       ]  [   o     ]                          [ , o     ]  [ ,   o   ]
        [ [  \      ]  [  /      ]  [        ]  [        ]  [    \    ]  [    /    ]
        [ [   o, ,  ], [ o  , ,  ], [ o, o,  ], [ o, , o ], [     o,  ], [   o  ,  ],
        <BLANKLINE>
                    [         ]  [         ] ]
                    [ , , o   ]  [ , ,   o ] ]
        [        ]  [      \  ]  [      /  ] ]
        [ , o, o ], [       o ], [     o   ] ]
        sage: BSeq23.first().parent()
        Sequence of `Binary trees`
    """

    class GradedComponentByLengthSum(Structures.GradedComponent):
        """
        `\mathbf{Seq}(F).graded_component((i, j))` with:
            - `i` the sum of the grading element in the sequence and
            - `j` the length of the sequence
        """

        def cardinality(self):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: BSeq = B.sequence(); BSeq
                Sequence of `Binary trees`
                sage: BSeq.graded_component((2,3)).cardinality()
                9
            """
            k, length = self.grading()
            F = self.ambient()._F
            return sum(map(
                lambda vect: prod(F.graded_component(i).cardinality() for i in vect),
                IntegerVectors(k, length)
            ))

        def __iter__(self):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: BSeq = B.sequence()
                sage: ascii_art(BSeq.graded_component((2,3)).list())
                [                                                   [         ]  [         ]
                [ [ o       ]  [   o     ]                          [ , o     ]  [ ,   o   ]
                [ [  \      ]  [  /      ]  [        ]  [        ]  [    \    ]  [    /    ]
                [ [   o, ,  ], [ o  , ,  ], [ o, o,  ], [ o, , o ], [     o,  ], [   o  ,  ],
                <BLANKLINE>
                            [         ]  [         ] ]
                            [ , , o   ]  [ , ,   o ] ]
                [        ]  [      \  ]  [      /  ] ]
                [ , o, o ], [       o ], [     o   ] ]

                sage: ascii_art(BSeq.graded_component((1,5)).list())
                [ [           ]  [           ]  [           ]  [           ]  [           ] ]
                [ [ o, , , ,  ], [ , o, , ,  ], [ , , o, ,  ], [ , , , o,  ], [ , , , , o ] ]
            """
            k, length = self.grading()
            F = self.ambient()._F

            for vect in IntegerVectors(k, length):
                for obj in product(*(F.graded_component(i) for i in vect)):
                    yield self._element_constructor_(obj)

    class GradedComponentBySum(Structures.GradedComponent):
        """
        `\mathbf{Seq}(F).graded_component(i)` with `i` the sum of the grading elements in the sequence
        """

        def cardinality(self):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: B1Seq = B.restricted_structures(min=1).sequence(grading_set="sum")
                sage: B1Seq3 = B1Seq.graded_component(3)
                sage: B1Seq3.cardinality()
                10

                sage: BSeq = B.sequence(grading_set="sum"); BSeq
                Sequence of `Binary trees`
                sage: BSeq.graded_component(4).cardinality()
                Traceback (most recent call last):
                ...
                Exception: `Binary trees` contains an element of grading 0, so the cardinality of the sequence is not defined.

            """
            if self.ambient()._F.graded_component(0).cardinality() != 0:
                raise Exception("`%s` contains an element of grading 0, so the cardinality"%self.ambient()._F + \
                                " of the sequence is not defined.")
            else: return sum(self.ambient().GradedComponentByLengthSum(self.ambient(),
                                                                       (self.grading(), i)).cardinality()
                             for i in range(self.grading() + 1)
                         )
            # FIXME: that suppose *self.grading* is an integer but... we want (may be) the sum or not...
            # what append with the sequence of sequence??

        def __iter__(self):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: B1Seq = B.restricted_structures(min=1).sequence(grading_set="sum")
                sage: B1Seq3 = B1Seq.graded_component(3); B1Seq3
                Sequence of `Binary trees with grading min=`1`` of degree 3
                sage: ascii_art(B1Seq3.list())
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

                sage: BSeq = B.sequence(grading_set="sum")
                sage: B4 = BSeq.graded_component(4); B4
                Sequence of `Binary trees` of degree 4
                sage: B4.list()
                Traceback (most recent call last):
                ...
                Exception: `Binary trees` contains an element of grading 0, so the cardinality of the sequence is not defined.
            """
            if self.ambient()._F.graded_component(0).cardinality() != 0:
                raise Exception("`%s` contains an element of grading 0, so the cardinality"%self.ambient()._F + \
                                " of the sequence is not defined.")
            return chain(*(self.ambient().GradedComponentByLengthSum(self.ambient(), (self.grading(), i))
                           for i in range(self.grading() + 1)))

    _graded_component = {"both" : (GradedComponentByLengthSum,
                                   lambda F: NonNegativeIntegers().cartesian_product(F.grading_set())),

                         "sum"  : (GradedComponentBySum,
                                   lambda F: F.grading_set())}

    def __init__(self, F, grading_set="both"):
        """
        @param F: a class of combinatorial structures
        @param grading_set: a string which could be:
           - "sum" : that means `\mathbf{Seq}(F)` is only graded by the grading sum of the sequence elements
        (with that grading set, one has to be sure `[0]F(x) = 0`).
           - "both": `\mathbf{Seq}(F)` is graded by the sum of the sequence elements and by the length of the sequence
        (with that constraint, the enumeration of the sequence is always well defined).

        TESTS::

            sage: B = BinaryTrees()
            sage: BSeq = B.sequence(); BSeq
            Sequence of `Binary trees`
            sage: #TestSuite(BSeq).run() # FIXME: NOT WORKING

        """
        # FIXME: that is correct??
        super(Sequence, self).__init__(self)

        assert(F in CombinatorialStructures())
        self._F = F

        self.GradedComponent = self._graded_component[grading_set][0]
        self._grading_set = self._graded_component[grading_set][1]

    def grading_set(self):
        """
        TESTS::

            sage: B = BinaryTrees()
            sage: BS = B.sequence()
            sage: BS.grading_set()
            The cartesian product of (Non negative integers, Non negative integers)
            sage: BS2 = B.sequence(grading_set="sum")
            sage: BS2.grading_set()
            Non negative integers
        """
        return self._grading_set(self._F)

    def graded_component(self, k, length=None):
        """
        TESTS::

            sage: B = BinaryTrees()
            sage: BSeq = B.sequence()
            sage: BSeq.graded_component((2,4)) == BSeq.graded_component(2,4)
            True
            sage: BSeq.graded_component((2,4))
            Sequence of `Binary trees` of degree (2, 4)
        """

        if length: k = (k, length)

        return _Operations.graded_component(self, k)

    def _repr_(self):
        """
        TESTS::

            sage: B = BinaryTrees()
            sage: B.sequence()
            Sequence of `Binary trees`

        """

        return "Sequence of `" + repr(self._F) + "`"

    def generating_series(self):
        """
        The generating series the sequence of `F` is given by:

        MATH::

            Seq_F(t) = \frac{1}{1 - f(t)}

        with `f(t)` the generating series of `F` and `[0]f(t) = 0`.

        """
        return Integer(1) / (Integer(1) - self._structures[0].generating_series())
        
    class Element(Structure, ClonableArray):
        """
        TESTS::

            sage: B = BinaryTrees()
            sage: BSeq = B.sequence()
            sage: BSeq.graded_component((2,3)).first().parent() is BSeq
            True

        """

        def check(self):
            pass

        def _ascii_art_(self):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: BSeq = B.sequence()
                sage: ascii_art(BSeq.graded_component((2,3)).list())
                [                                                   [         ]  [         ]
                [ [ o       ]  [   o     ]                          [ , o     ]  [ ,   o   ]
                [ [  \      ]  [  /      ]  [        ]  [        ]  [    \    ]  [    /    ]
                [ [   o, ,  ], [ o  , ,  ], [ o, o,  ], [ o, , o ], [     o,  ], [   o  ,  ],
                <BLANKLINE>
                            [         ]  [         ] ]
                            [ , , o   ]  [ , ,   o ] ]
                [        ]  [      \  ]  [      /  ] ]
                [ , o, o ], [       o ], [     o   ] ]

                sage: ascii_art(BSeq.graded_component((1,5)).list())
                [ [           ]  [           ]  [           ]  [           ]  [           ] ]
                [ [ o, , , ,  ], [ , o, , ,  ], [ , , o, ,  ], [ , , , o,  ], [ , , , , o ] ]
            """
            return ascii_art_list(self)