# -*- coding: utf-8 -*-
r"""
Cycle class of combinatorial structures.

Let `F` a class of combinatorial structures, the *cycle* class `\mathbf{Cyc}(F)`
is defined to be the infinite sum:

MATH::

    \mathbf{Cyc}(F) := \mathbf{Seq}(F)/\mathcal{R}

with `\mathcal{R}` being the equivalence relation on sequences being defined by
`s_1 \mathcal{R} s_2` if and only if it exists some cyclic permutation
`\sigma` such that `s_1 \sigma = s_2` _[FS].


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
from sage.combinat.permutation import CyclicPermutations
from sage.combinat.structures.operations.multi_sets import MultiSet
from sage.combinat.structures.operations.sequences import Sequence
from sage.misc.ascii_art import ascii_art_tuple
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.list_clone import ClonableArray
from sage.combinat.structures import Structures, Structure


class Cycle(Sequence):
    r"""
    Cycle class of combinatorial structures.

    Let `F` be class of combinatorial structures, the *cycle* class
    `\mathbf{Cyc}(F)` is defined to be the infinite sum:

    MATH::

        \mathbf{Cyc}(F) := \mathbf{Seq}(F)/\mathcal{R}

    with `\mathcal{R}` being the equivalence relation on sequences being defined
    by `s_1 \mathcal{R} s_2` if and only if it exists some cyclic permutation
    `\sigma` such that `s_1 \sigma = s_2` _[FS].

    TESTS::

        sage: C = Compositions()
        sage: CC = C.cycle(); CC
        Cycle of `Compositions of non-negative integers`

    """

    class GradedComponentByLengthSum(Structures.GradedComponent):
        r"""
        `\mathbf{Cyc}(F).graded_component((i, j))` with:
            - `i` the sum of the grading element in the sequence and
            - `j` the length of the sequence
        """

        def cardinality(self):
            """
            TESTS::

                sage: C = Compositions()
                sage: CC = C.cycle()
                sage: len(CC.graded_component(3,2).list())
                6
                sage: CC.graded_component(3,2).cardinality()
                6

                sage: B = BinaryTrees()
                sage: BC = B.restricted_structures(min=1).cycle(grading_set="sum")
                sage: len(BC.graded_component(3).list())
                8
                sage: BC.graded_component(3).cardinality()
                8
            """
            # FIXME
            F = self.ambient()._F
            acc = 0
            for mset in MultiSet(F).graded_component(self.grading()):
                acc += CyclicPermutations(mset).cardinality()
            return acc

        def __iter__(self):
            """
            TESTS::

                sage: C = Compositions()
                sage: CC = C.cycle()
                sage: CC.graded_component(3,2).list()
                [{[1], [1, 1]},
                 {[1], [2]},
                 {[], [1, 1, 1]},
                 {[], [1, 2]},
                 {[], [2, 1]},
                 {[], [3]}]

                sage: B = BinaryTrees()
                sage: BC = B.restricted_structures(min=1).cycle(grading_set="sum"); BC
                Cycle of `Binary trees with grading min=`1``
                sage: BC.graded_component(3).list()
                [{[., [., [., .]]]},
                 {[., [[., .], .]]},
                 {[[., .], [., .]]},
                 {[[., [., .]], .]},
                 {[[[., .], .], .]},
                 {[., .], [., [., .]]},
                 {[., .], [[., .], .]},
                 {[., .], [., .], [., .]}]
            """
            F = self.ambient()._F

            for mset in MultiSet(F).graded_component(self.grading()):
                for cyc in CyclicPermutations(mset):
                    yield self._element_constructor_(cyc)


    _graded_component = {"both": (GradedComponentByLengthSum,
                                   lambda F: NonNegativeIntegers().cartesian_product(F.grading_set())),

                         "sum": (Sequence.GradedComponentBySum,
                                   lambda F: F.grading_set())}

    def _repr_(self):
        """
        TESTS::

            sage: B = BinaryTrees()
            sage: B.cycle()
            Cycle of `Binary trees`

        """

        return "Cycle of `" + repr(self._F) + "`"

    def generating_series(self):
        r"""
        Return the generating serie of ``self`` which is given by:

        MATH::

            Cyc_F(t) = ...

        with `f(t)` the generating series of `F` and `[0]f(t) = 0`.

        """
        ##return Integer(1) / (Integer(1) - self._structures[0].generating_series())

    class Element(Structure, ClonableArray):
        """
        TESTS::

            sage: B = BinaryTrees()
            sage: BC = B.cycle()
            sage: BC.graded_component((4,3)).first().parent() is BC
            True

        """

        def __init__(self, parent, li):
            ClonableArray.__init__(self, parent, sorted(li))

        def _repr_(self):
            """
            TESTS::

                sage: C = Compositions()
                sage: CC = C.cycle()
                sage: CC.graded_component(3,2).list()
                [{[1], [1, 1]},
                 {[1], [2]},
                 {[], [1, 1, 1]},
                 {[], [1, 2]},
                 {[], [2, 1]},
                 {[], [3]}]

            """
            return "{" + ClonableArray._repr_(self)[1:-1] + "}"

        def check(self):
            pass

        def _ascii_art_(self):
            """
            TESTS::

                sage: C = Compositions()
                sage: CC = C.cycle()
                sage: ascii_art(CC.graded_component(3,2).list())
                [                      (   * )                                ]
                [ (    * )  (       )  (   * )  (   ** )  (    * )  (       ) ]
                [ ( *, * ), ( *, ** ), ( , * ), ( , *  ), ( , ** ), ( , *** ) ]
            """
            return ascii_art_tuple(self)
