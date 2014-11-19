# -*- coding: utf-8 -*-
"""
Framework of class of combinatorial structures

A class of combinatorial structures is a *denumerable set* of discrete objects
(structures) on which a *degree* function is defined, satisfying the following
conditions:

   (i)  the degree of an element is a non-negative integer;
   (ii) the number of elements of any given degree is finite.


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
from sage.structure.element import Element
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.structure.parent import Parent
from sage.categories.classes_of_combinatorial_structures import \
    ClassesOfCombinatorialStructures
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation


class Structure(Element):
    """
    A *structure* is a discrete object with a degree notion. It is an element of
    a combinatorial class of structures.

    For example, a (combinatorial) structure is a *composition of integer*.
    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall__(cls, *args, **options):
        """
        The constructor of an Element class must have explicit argument *parent*
        but most of the time, the user does not want explicit this argument::

            Permutation([1,3,2])

        This method is called before the *__init__* method and provides a
        compromise.
        If user does not want specify a parent, this method looks for a method
        *_auto_parent_*. Otherwise, the parent could be specified as an optional
        argument "parent"::

            Composition([2,1,2], parent=Parent())

        TESTS::

            sage: from sage.categories.examples.\
                  classes_of_combinatorial_structures import Composition
            sage: Composition([2,1,2], parent=Parent()).parent()
            <type 'sage.structure.parent.Parent'>
            sage: Composition([2,1,2]).parent()
            Compositions of integers

        """
        par = "parent"
        if options.has_key(par):
            parent = options[par]
            del options[par]
        else:
            try:
                parent = cls._auto_parent_
            except AttributeError:
                raise NotImplementedError("A *Structure* (%s) should"%cls +
                                          " implement *_auto_parent_*")

        return typecall(cls, parent, *args, **options)


class StructuresClass(UniqueRepresentation, Parent):
    """
    *Structures* is the python class of the classes of combinatorial structures.

    For example, the set of all *compositions of integers*.

    TESTS::

        sage: from sage.categories.examples.\
                   classes_of_combinatorial_structures import Compositions
        sage: TestSuite(Compositions()).run()

    """

    def __init__(self, category=ClassesOfCombinatorialStructures(),
                 *args, **options):
        """

        TESTS::

            sage: from sage.categories.examples.\
                  classes_of_combinatorial_structures import Compositions
            sage: C = Compositions(); C
            Compositions of integers
            sage: C.graded_component(2)
            Compositions of integers of degree 2
            sage: C.graded_component(18).ambient()
            Compositions of integers
            sage: C.graded_component(18).first().parent()
            Compositions of integers

        """
        Parent.__init__(self, category=category)

    def _element_constructor_(self, *args, **options):
        """
        Redefinition of that method to be coherent with the *_classcall_* of
        *Structure*.

        TESTS::

            sage: from sage.categories.examples.\
                  classes_of_combinatorial_structures import Compositions
            sage: C = Compositions()
            sage: C([3,2,2])
            [3, 2, 2]
            sage: C._element_constructor_([2,2,2])
            [2, 2, 2]

        """
        return self.element_class(parent=self, *args, **options)

    def graded_component(self, grading):
        """
        Return the graded component of degree *grading*.

        TESTS::

            sage: from sage.categories.examples.\
                  classes_of_combinatorial_structures import Compositions
            sage: Compositions().graded_component(4)
            Compositions of integers of degree 4

        """
        return self.GradedComponent(self, grading)

    def _an_element_(self):
        """
        Default implementation

        TESTS::

            sage: from sage.categories.examples.\
                  classes_of_combinatorial_structures import Compositions
            sage: I = Compositions()._an_element_(); I
            []
            sage: I.parent()
            Compositions of integers

        """
        for n in self.grading_set():
            try:
                return self.graded_component(n).first()
            except:
                pass

    def some_elements(self):
        """
        TESTS::

            sage: from sage.categories.examples.\
                  classes_of_combinatorial_structures import Compositions
            sage: list(Compositions().some_elements())
            [[1], [2], [3], [4]]

        """
        for i in range(1,5):
            if self.graded_component(i).cardinality() > 0:
                yield self.graded_component(i).first()

    class GradedComponent(UniqueRepresentation, Parent):
        """
        A *graded component* (of degree `n`) of a class of combinatorial structures
        is the finite set of all structures of degree `n`.

        TESTS::

            sage: from sage.categories.examples.\
                  classes_of_combinatorial_structures import Compositions
            sage: TestSuite(Compositions().graded_component(4)).run()

        """

        def __init__(self, ambient, grading,
                     category=ClassesOfCombinatorialStructures.GradedComponents()):
            """

            """
            Parent.__init__(self, category=category)
            self._ambient_ = ambient
            self._grading_ = grading

        def _repr_(self):
            return repr(self.ambient()) + " of degree " + repr(self.grade())

        def ambient(self):
            return self._ambient_

        def grade(self):
            return self._grading_

        @lazy_attribute
        def _element_constructor_(self, *args, **opts):
            return self.ambient()._element_constructor_
