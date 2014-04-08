r"""
Set Operads
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category
from sage.categories.sets_cat import Sets
from sage.misc.cachefunc import cached_method
from sage.categories.cartesian_product import cartesian_product
from sage.misc.abstract_method import abstract_method

class SetOperads(Category):
    """
    The category of set operads

    EXAMPLES::

      sage: SetOperads()
      Category of set operads
      sage: SetOperads().super_categories()
      [Category of sets]

    TESTS::

        sage: C = SetOperads()
        sage: TestSuite(C).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: SetOperads().super_categories()
            [Category of sets]
        """
        return [Sets()]

    def example(self):
        """
        Returns an example of set operad ::

            sage: SetOperads().example()
            An example of a set operad: the Associative operad
        """
        from sage.categories.examples.set_operads import Example
        return Example()

    class ParentMethods:

        @abstract_method(optional = True)
        def composition(self, left, right, index):
            """
            returns the composition of left with right at position index
            """

        @abstract_method(optional = True)
        def composition_with_numbers(self):
            """
            This is a variant of composition, where one assumes that
            the objects are labelled by integers from 1 to n. The
            result is labelled in the same way.
            """

        def global_composition(self, left, list_right):
            r"""
            returns the global composition of left with a list of elements
            """
            if self.composition is not NotImplemented:
                assert left.degree() == len(list_right), "degree of x is not equal to the length of list_right"
                res = left
                for i in xrange(left.degree(), 0, -1):
                    res = res.compose(list_right[i - 1], i)
                    return res
            else:
                return NotImplemented

        def global_composition_with_numbers(self, left, list_right):
            r"""
            returns the global composition of left with a list of elements
            """
            if self.composition_with_numbers is not NotImplemented:
                assert left.degree() == len(list_right), "degree of x is not equal to the length of list_right"
                res = left
                for i in xrange(left.degree(), 0, -1):
                    res = res.compose_with_numbers(list_right[i - 1], i)
                    return res
            else:
                return NotImplemented

        @abstract_method(optional = True)
        def operad_morphism(self, arg, codomain):
            """
            returns the image of arg by a morphism from self to codomain
            """

        @abstract_method(optional = True)
        def one(self,letter):
            """
            returns the one of the operad
            """

        @abstract_method(optional = True)
        def is_symmetric(self):
            r"""
            returns `True` if the operad is symmetric
            """
            pass

        @abstract_method(optional = True)
        def elements(self, n):
            """
            returns the set of elements in degree `n`
            """
            pass

        def cardinality(self, n):
            """
            returns the cardinality in degree `n`
            """
            return len(self.elements(n))


    class ElementMethods:

        def compose(self, other, index):
            """
            returns the composition of self with other at position index

            EXAMPLES::
            """
            return self.parent().composition(self, other, index)

        def compose_with_numbers(self, other, index):
            """
            returns the composition of self with other at position index

            EXAMPLES::
            """
            return self.parent().composition_with_numbers(self, other, index)

        @abstract_method(optional = True)
        def degree(self, x):
            """
            returns the degree of an element
            """
            pass

        @abstract_method(optional = True)
        def map_labels(self, f):
            """
            applies the function `f` to the labels of an element
            """
            pass
