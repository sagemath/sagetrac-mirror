# -*- coding: utf-8 -*-
"""
The category of classes of combinatorial structures
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
        return [SetsWithGrading(), InfiniteEnumeratedSets()]

    class ParentMethods:

        @staticmethod
        def product(F, G):
            """
            @param F, G: both are classes of combinatorial structures
            @return: the product of combinatorial structures *FG*
            """

        @staticmethod
        def cartesian_product(F, G):
            """
            @param F, G: both are classes of combinatorial structures
            @return: the cartesian product of combinatorial structures *F \oblong G*
            """

        @staticmethod
        def sum(F, G):
            """
            @param F, G: both are classes of combinatorial structures
            @return: the sum of combinatorial structures *F + G*
            """
        __add__ = sum

        def restrict_class_of_structures(self, min, max):
            """
            @param min, max: both define an interval of graduation of the class of combinatorial structures
            @return: the restriction of the combinatorial structure to the homogeneous components between *min* and
                     *max*.
            """

    class GradedComponents(Category):

        def super_categories(self):
            return [FiniteEnumeratedSets()]

        class ParentMethods:

            def _repr_(self):
                return repr(self.ambient()) + " of degree " + repr(self.grade())

            @abstract_method(optional=False)
            def ambient(self):
                """
                @return: the combinatorial class of the homogeneous component
                """

            @abstract_method(optional=False)
            def grading(self):
                """
                @return: the grading of the homogeneous component
                """