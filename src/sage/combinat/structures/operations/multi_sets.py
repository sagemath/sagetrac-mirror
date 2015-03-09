# -*- coding: utf-8 -*-
"""
Class of multi set of combinatorial structures.

Let `F` be a class of combinatorial structures, the *multi set* class
`\mathbf{MSet}(F)` is defined to be the infinite sum:

MATH::

    \mathbf{MSet}(F) := \mathbf{Seq}(F)/\mathcal{R}

with `\mathcal{R}` being the equivalence relation on sequences defined by `s_1
\mathcal{R} s_2` if and only if it exists some permutation `\sigma` such that
`s_1 \sigma = s_2` _[FS].


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
from itertools import product, combinations_with_replacement
from sage.combinat.partition import Partitions
from sage.combinat.structures.operations.sequences import Sequence
from sage.misc.ascii_art import ascii_art_set
from sage.misc.misc_c import prod
from sage.rings.arith import multinomial, factorial
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.list_clone import ClonableArray
from sage.combinat.structures import Structures, Structure


class MultiSet(Sequence):
    r"""
    Class of multi set of combinatorial structures.

    Let `F` be a class of combinatorial structures, the *multi set* class
    `\mathbf{MSet}(F)` is defined to be the infinite sum:

    MATH::

        \mathbf{MSet}(F) := \mathbf{Seq}(F)/\mathcal{R}

    with `\mathcal{R}` being the equivalence relation on sequences defined by
    `s_1 \mathcal{R} s_2` if and only if it exists some permutation `\sigma`
    such that `s_1 \sigma = s_2` _[FS].

    TESTS::

        sage: C = Compositions()
        sage: MC = C.restricted_structures(min=1).multiset(grading_set="sum")
        sage: MC.graded_component(4).list()
        [{[1, 1, 1, 1]},
         {[1, 1, 2]},
         {[1, 2, 1]},
         {[1, 3]},
         {[2, 1, 1]},
         {[2, 2]},
         {[3, 1]},
         {[4]},
         {[1], [1, 1, 1]},
         {[1], [1, 2]},
         {[1], [2, 1]},
         {[1], [3]},
         {[1, 1], [1, 1]},
         {[1, 1], [2]},
         {[2], [2]},
         {[1], [1], [1, 1]},
         {[1], [1], [2]},
         {[1], [1], [1], [1]}]

    """

    class GradedComponentByLengthSum(Structures.GradedComponent):
        r"""
        `\mathbf{MSet}(F).graded_component((i, j))` with:
            - `i` the sum of the grading element in the sequence and
            - `j` the length of the sequence
        """

        def cardinality(self):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: BSet = B.multiset()
                sage: BSet.graded_component(2,3).cardinality()
                3
                sage: BSet2 = B.restricted_structures(min=1).multiset(grading_set="sum")
                sage: #[BSet2.graded_component(i).cardinality() for i in range(10)]

                #[1, 1, 3, 8, 25, 77, 256, 854, 2940, 10229]
            """
            # FIXME: ...
            k, length = self.grading()
            if length == 0 and k == 0:
                return 1

            F = self.ambient()._F.graded_component

            hasZero = F(0).cardinality() > 0
            acc = 0
            for nbZero in range(length):
                if nbZero > 0 and not hasZero:
                    continue

                for I in Partitions(k, length=length - nbZero):
                    I = tuple(I) + (0,) * nbZero

                    acc += prod(F(i).cardinality() for i in I)

            return acc

        def __iter__(self):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: BSet = B.multiset()
                sage: ascii_art(BSet.graded_component(2,3).list())
                [             {         }  {         } ]
                [             { , , o   }  { , ,   o } ]
                [ {        }  {      \  }  {      /  } ]
                [ { , o, o }, {       o }, {     o   } ]

                sage: BSet2 = B.restricted_structures(min=1).multiset(grading_set="sum")
                sage: [len(BSet2.graded_component(i).list()) for i in range(7)]
                [1, 1, 3, 8, 25, 77, 256]
            """
            k, length = self.grading()
            if length == 0 and k == 0:
                yield self._element_constructor_(())

            F = self.ambient()._F.graded_component

            hasZero = F(0).cardinality() > 0

            for nbZero in range(length):
                if nbZero > 0 and not hasZero:
                    continue

                for I in Partitions(k, length=length - nbZero):

                    I = list(I)
                    iterators = []
                    ### These lines are use to do not repeat tuple of elements of     ###
                    ### the same graded component.                                    ###
                    for i in set(I):
                        iterators.append(combinations_with_replacement(F(i), I.count(i)))

                    iterators.append(combinations_with_replacement(F(0), nbZero))
                    #####################################################################

                    for tup in product(*iterators):
                        yield self._element_constructor_(reduce(lambda a, b: a+b, tup, ()))

    _graded_component = {"both": (GradedComponentByLengthSum,
                                   lambda F: NonNegativeIntegers().cartesian_product(F.grading_set())),

                         "sum": (Sequence.GradedComponentBySum,
                                   lambda F: F.grading_set())}

    def _repr_(self):
        """
        TESTS::

            sage: B = BinaryTrees()
            sage: B.multiset()
            Multi-set of `Binary trees`

        """

        return "Multi-set of `" + repr(self._F) + "`"

    def generating_series(self):
        r"""
        Return the generating serie of ``self``.

        The generating serie is given by:

        MATH::

            Seq_F(t) = ...

        with `f(t)` the generating serie of `F` and `[0]f(t) = 0`.

        """
        ##return Integer(1) / (Integer(1) - self._structures[0].generating_series())

    class Element(Structure, ClonableArray):
        """
        TESTS::

            sage: B = BinaryTrees()
            sage: BSet = B.multiset()
            sage: BSet.graded_component((4,3)).first().parent() is BSet
            True

        """

        def __init__(self, parent, li):
            ClonableArray.__init__(self, parent, sorted(li))

        def _repr_(self):
            """
            TESTS::

                sage: C = Compositions()
                sage: MC = C.restricted_structures(min=1).multiset()
                sage: MC.graded_component(3,2).list() # order not determinist
                [{[1], [1, 1]}, {[1], [2]}]

            """
            return "{" + ClonableArray._repr_(self)[1:-1] + "}"

        def check(self):
            pass

        def _ascii_art_(self):
            """
            TESTS::

                sage: C = Compositions()
                sage: MC = C.restricted_structures(min=1).multiset()
                sage: ascii_art(MC.graded_component(3,2).list()) # order not determinist
                [ {    * }  {       } ]
                [ { *, * }, { *, ** } ]
            """
            return ascii_art_set(self)
