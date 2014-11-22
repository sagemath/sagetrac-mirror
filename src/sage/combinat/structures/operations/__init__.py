# -*- coding: utf-8 -*-
"""
That module implements several operations on class of combinatorial structures

    - Sum
    - Cauchy product
    - Hadamard product
    - Differentiation
    - ...

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
from itertools import imap
from sage.categories.combinatorial_structures import CombinatorialStructures
from sage.combinat.structures import Structures, StructuresWithArguments

class _Operations(Structures):
    """
    That class is an interface of operation like Sum, Cartesian Product...

    EXAMPLES::

        sage: B, C = BinaryTrees(), Compositions()
        sage: B + C
        Sum of structures : `Binary trees`, `Compositions of non-negative integers`

    """

    @staticmethod
    def __classcall__(cls, *structures, **options):
        """
        Method to overload the *Structures* classcall method

        TESTS::

            sage: B, C = BinaryTrees(), Compositions()
            sage: B + C
            Sum of structures : `Binary trees`, `Compositions of non-negative integers`

        """
        return super(Structures, cls).__classcall__(cls, *structures, **options)

    def __init__(self, *structures):
        """

        @param structures: a sequences of *Structures*

        TESTS::

            sage: B, C = BinaryTrees(), Compositions()
            sage: BpC = B * C
            sage: BpC._structures
            (Binary trees, Compositions of non-negative integers)

        """
        Structures.__init__(self)
        self._structures = structures

class RestrictedStructures(StructuresWithArguments):
    """
    That class is a wrapper use to restrict a class of combinatorial structures
    to some graded components.

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

    def __init__(self, F, min=None, max=None):
        """
        @param F: a class of combinatorial structures
        @param min, max: both element of the grading set


        TESTS::

            sage: B = BinaryTrees()
            sage: RB = B.restricted_structures(min=3, max=4); RB
            Binary trees with grading min=`3`, max=`4`

        """
        Structures.__init__(self)
        assert(F in CombinatorialStructures()), "`%s` has to be a class of combinatorial structures"
        self._F = F

        assert((not min or min in self.grading_set()) and (not max or max in self.grading_set())), \
                "min and max have to be element of the grading set"
        self._min = min
        self._max = max

    def _repr_(self):
        """
        TESTS::

            sage: B = BinaryTrees()
            sage: RB = B.restricted_structures(min=3, max=4); RB
            Binary trees with grading min=`3`, max=`4`

        """
        # TODO improve the representation
        s = ""
        if self._min:
            s += "min=`" + repr(self._min) + "`"
            if self._max:
                s += ", max=`" + repr(self._max) + "`"
        else:
            s += "max=`" + repr(self._max) + "`"
        return repr(self._F) + " with grading " + s

    def generating_series(self, variable="x"):
        """
        TESTS::

        """
        # TODO: find a canonical way to define the generating series!!!
        max = (self._max + 1) if self._max is not None else None
        return self._F.generating_series(variable)[self._min:max]

    def _element_constructor_(self, *args, **options):
        """
        The element constructor use the original one and change the parent

        TESTS::

            sage: B = BinaryTrees()
            sage: RB = B.restricted_structures(min=3, max=4); RB
            Binary trees with grading min=`3`, max=`4`
            sage: RB.graded_component(3).list()[0].parent() is RB
            True

        """
        obj = self._F._element_constructor_(*args, **options)
        obj._set_parent(self)
        return obj


    class GradedComponent(Structures.GradedComponent):

        def __iter__(self):
            """
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

            """
            ambient = self.ambient()
            F = ambient._F
            k = self.grading()
            element_constructor = ambient._element_constructor_

            if (ambient._min and k < ambient._min) or \
               (ambient._max and k > ambient._max):
                return iter([])

            return imap(
                element_constructor,
                F.graded_component(k)
            )
