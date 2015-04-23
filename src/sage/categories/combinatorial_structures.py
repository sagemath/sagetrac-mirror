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
    The category of classes of combinatorial structures

    "In mathematics, a combinatorial class is a countable set of mathematical
    objects." (Wikipedia)

    A class of combinatorial structures is a graded set `C = \bigsqcup_{n \geqslant 0} C_n` such that each graded
    component `C_n` is a finite set.

    EXAMPLES::

        sage: Partitions().category()
        Category of combinatorial structures

    The generating series of `C` is a *formal power series*

    MATH::

        C(t) = \sum_{n \geqslant 0} c_n t^n\,,

    where `c_n` is the cardinality of the graded component `C_n`.

    .. TODO waiting for TICKET species (u/elixyre/species)

    """

    def super_categories(self):
        return [SetsWithGrading(), InfiniteEnumeratedSets()]

    class ParentMethods:

        @abstract_method(optional=True)
        def generating_series(self):
            """
            The generating series of the class of combinatorial structure.

            MATH::

                C(t) = \sum_{n \geqslant 0} c_n t^n\,.

            """
            # TODO waiting for TICKET species (u/elixyre/species)

        def gs(self): return self.generating_series()

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
                    sage: Compositions().graded_component(4)
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
