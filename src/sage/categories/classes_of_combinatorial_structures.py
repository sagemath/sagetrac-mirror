# -*- coding: utf-8 -*-
"""
Classes of combinatorial structures

A class of combinatorial structures is a *denumerable set* of discrete objects
(structures) on which a *degree* function is defined, satisfying the following
conditions:

   (i)  the degree of an element is a non-negative integer;
   (ii) the number of elements of any given degree is finite.


REFERENCES:
-----------

.. [FlaSed] Analytic combinatorics,
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
from sage.categories.category import Category
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.misc.abstract_method import abstract_method


class ClassesOfCombinatorialStructures(Category):
    """
    A class of combinatorial structures is a *denumerable set* of discrete
    objects (structures) on which a *degree* function is defined, satisfying the
    following conditions:

       (i)  the degree of a structure is a non-negative integer;
       (ii) the number of structures of any given degree is finite.

    To implement this category, a parent must implement methods excepted in
    the category ``SetsWithGrading`` and
    """

    def super_categories(self):
        return [SetsWithGrading()]

    class ParentMethods:

        def _test_graded_components_2(self, **options):
            """
            Test that some graded components of ``self`` are that the parent
            has a properly implemented ``ambient`` and ``grade`` methods.

            EXAMPLE::

                sage: ClassesOfCombinatorialStructures().example().\
                      _test_graded_components_2()

            """
            tester = self._tester(**options)
            for grade in self.grading_set().some_elements():
                G = self.graded_component(grade)

                for elt in G.some_elements():
                    tester.assertEqual(G.grade(), self.grading(elt))

                tester.assertEqual(G.ambient().graded_component(grade), G)

    class GradedComponents(Category):
        """
        A graded component of a class of combinatorial structures `C` is an
        finite set of structures which all have same degree `i`:
        the graded component of degree `i`.
        """

        def super_categories(self):
            return [FiniteEnumeratedSets()]

        class ParentMethods:

            @abstract_method(optional=False)
            def ambient(self):
                """
                The underlying class of the graded component.

                EXAMPLE::

                    sage: C = ClassesOfCombinatorialStructures().example(); C
                    Compositions of integers
                    sage: C.graded_component(3).ambient()
                    Compositions of integers

                """

            @abstract_method(optional=False)
            def grade(self):
                """
                Return the degree associated to  all structures of the component.

                EXEMPLE::

                    sage: ClassesOfCombinatorialStructures().example().\
                          graded_component(42).grade()
                    42

                """
