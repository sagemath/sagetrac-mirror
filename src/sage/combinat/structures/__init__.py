# -*- coding: utf-8 -*-
"""
Framework of class of combinatorial structures

Some notes extracts from _[BLL] and _[FL] and implementation interpretation.
----------------------------------------------------------------------------

"A *structure* `s` is a *construction* `\gamma` which one performs on a set `U` of data." and
"The ''general form'' that isomorphic structures have in common is their *isomorphism type*. [...]
the elements of the underlying set are represented by ''indistinguishable'' points."_[BLL]


In that implementation, we distinct *species* (or more commonly *labeled structures*) and *unlabeled structures*.


References:
-----------

.. [BLL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle, and Pierre Leroux

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
from sage.rings.integer import Integer
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.unique_representation import UniqueRepresentation

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.categories.combinatorial_structures import CombinatorialStructures


class Structure(Element):
    """
    A *structure* describe the structure of a combinatorial object.

    TESTS::

        sage: from sage.categories.examples.combinatorial_structures_compositions import Composition
        sage: Composition([3,1,2,2])
        [3, 1, 2, 2]
    """

    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall__(cls, *args, **options):
        """
        This method is defined to allow the user to not give explicitly a parent in argument.

        If the user want use a specific parent then he has to pass the option *parent*.

        TESTS::

            sage: from sage.categories.examples.combinatorial_structures_compositions import Composition
            sage: I = Composition([3,1,2,2]); I
            [3, 1, 2, 2]
            sage: I.parent()
            Compositions of integers

            sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
            sage: J = Composition([2,2,2], parent=Compositions()); J
            [2, 2, 2]
            sage: J.parent()
            Compositions of integers

            sage: P = Parent()
            sage: K = Composition([2,2,2],parent=P); K
            [2, 2, 2]
            sage: K.parent()
            <type 'sage.structure.parent.Parent'>

        """
        if options.has_key("parent"):
            parent = options["parent"]
            del options["parent"]
        else:
            try:
                parent = cls._auto_parent_
            except AttributeError:
                # see:: for example *BinaryTree*
                raise AttributeError("A *Structure* (%s) must have an attribute *_auto_parent_*" %cls)
        return typecall(cls, parent, *args, **options)


class Structures(UniqueRepresentation, Parent):
    """
    *Structures* is the python class of the classes of combinatorial structures.

    TESTS::

        sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
        sage: TestSuite(Compositions()).run()
    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall__(cls, grade=None, *args, **options):
        """
        This method is a default shortcut of each class of combinatorial structure to obtain
        directly a graded component

        For example::

            sage: BinaryTrees()
            Binary trees
            sage: BinaryTrees(4)
            Binary trees of size 4

        The first call returns the class of all binary trees and
        the second returns only binary trees of size 4; that one is equivalent to::

            sage: BinaryTrees().graded_component(4)
            Binary trees of size 4

        .. WARNING:: if the user define is own method *_init_* with arguments then he has to overload that *_classcall_*
        method. --> see class *StructuresWithArguments*

        TESTS::

            sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
            sage: Compositions()
            Compositions of integers
            sage: Compositions(4)
            Compositions of integers of degree 4
            sage: Compositions(0)
            Compositions of integers of degree 0
        """
        self = super(Structures, cls).__classcall__(cls, *args, **options)

        if grade is not None and grade in cls.grading_set():
            return self.graded_component(grade)
        else: return self

    def __init__(self, category=CombinatorialStructures(), *args, **options):
        """
        This constructor is use to:
          - memoize the graded component
          - define the "ambient" as self of each graded component
          - define an auto parent of the elements


        NOTE:: in the context, it could be (may be) usefull to define a wrapper for each parent...

        @param category: by default the category is a (simple) combinatorial structures
                        but it is possible to use enriched category (like the category of species).

        TESTS::

            sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
            sage: C = Compositions(); C
            Compositions of integers
            sage: C.graded_component(2)
            Compositions of integers of degree 2
            sage: len(C._graded_components) > 0
            True

            sage: C.graded_component(18).ambient()
            Compositions of integers

            sage: C.graded_component(18).first().parent()
            Compositions of integers
        """
        Parent.__init__(self, category=category)
        self._graded_components = {}

    @staticmethod
    def grading_set():
        """
        This method should static to be call in classcall... (so that can not be defined by a category).

        TESTS::

            sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
            sage: Compositions.grading_set()
            Non negative integers
        """
        return NonNegativeIntegers()

    def _element_constructor_(self, *args, **options):
        """
        Redefinition of that method to be coherent with the *_classcall_* of *Structure*.

        TESTS::

            sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
            sage: C = Compositions()
            sage: C([3,2,2])
            [3, 2, 2]
            sage: C._element_constructor_([2,2,2])
            [2, 2, 2]
        """
        return self.element_class(parent=self, *args, **options)

    def graded_component(self, k):
        """
        TESTS::

            sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
            sage: C = Compositions()
            sage: C.graded_component(2)
            Compositions of integers of degree 2
            sage: len(C._graded_components) > 0
            True
        """
        # memoization of graded component: that is less fat than UniqueRepresentation...
        if self._graded_components.has_key(k):
            return self._graded_components[k]

        component = self.GradedComponent(self, k)
        self._graded_components[k] = component
        return component

    # FIXME:: should be implemented somewhere else no??
    def _an_element_(self):
        """
        Default implementation

        TESTS::

            sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
            sage: I = Compositions()._an_element_(); I
            []
            sage: I.parent()
            Compositions of integers

        """
        return self.first()

    def some_elements(self):
        """
        TESTS::

            sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
            sage: list(Compositions().some_elements())
            [[1], [2], [3], [4]]
        """
        for i in range(1,5):
            if self.graded_component(i).cardinality() > 0:
                yield self.graded_component(i).first()


    class GradedComponent(Parent):

        def __init__(self, ambient, grade, category=CombinatorialStructures.GradedComponents()):
            """
            TESTS::

                sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
                sage: Compositions(4)
                Compositions of integers of degree 4
                sage: Compositions().graded_component(4)
                Compositions of integers of degree 4
                sage: Compositions(3).list()
                [[3], [1, 2], [2, 1], [1, 1, 1]]
            """
            Parent.__init__(self, category=category)
            self._grade_ = grade
            self._ambient_ = ambient

        def __contains__(self, obj):
            """
            TESTS::

                sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions, \
                        Composition
                sage: C4 = Compositions(4)
                sage: I = Composition([2,2,1])
                sage: I in C4
                False
                sage: J = Composition([3,1])
                sage: J in C4
                True

            """
            if not isinstance(obj, self.ambient().element_class):
                try:
                    obj = self._element_constructor_(obj)
                except: return False
            return obj in self.ambient() and obj.grade() == self.grading()

        def ambient(self):
            """
            TESTS::

                sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
                sage: Compositions(4).ambient()
                Compositions of integers
            """
            return self._ambient_

        def grading(self):
            """
            TESTS::

                sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
                sage: Compositions(4).grading()
                4
            """
            return self._grade_

        @lazy_attribute
        def _element_constructor_(self, *args, **options):
            """
            TESTS::

                sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
                sage: Compositions(4)([1,2,1])
                [1, 2, 1]
            """
            return self.ambient()._element_constructor_

    Element = Structure

class StructuresWithArguments(Structures):
    """
    That class propose a generic support to implement a class of structures with arguments
    (such that the first argument is not an integer ~ a graduation)
    """

    @staticmethod
    def __classcall__(cls, grade=None, *args, **options):
        """
        This method is a default shortcut of each class of combinatorial structure to obtain
        directly a graded component
        """
        if len(args) == 0 or grade not in cls.grading_set():
            self = super(Structures, cls).__classcall__(cls, grade, *args, **options)
            return self
        else:
            self = super(Structures, cls).__classcall__(cls, *args, **options)
            return self.graded_component(grade)


