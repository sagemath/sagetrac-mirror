r"""
Example of a set with grading
"""

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_with_grading import SetsWithGrading

from sage.rings.integer_ring import IntegerRing
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet

class NonNegativeIntegers(UniqueRepresentation, Parent):
    r"""
    Non negative integers graded by themselves.

    EXAMPLES::

        sage: E = SetsWithGrading().example()
        sage: E
        Non negative integers
        sage: E.graded_component(0)
        {0}
        sage: E.graded_component(100)
        {100}
    """
    def __init__(self):
        r"""
        TESTS::

            sage: TestSuite(SetsWithGrading().example()).run()
        """
        Parent.__init__(self, category=SetsWithGrading(), facade=IntegerRing())

    def an_element(self):
        r"""
        Returns 0.

        EXAMPLES::

            sage: SetsWithGrading().example().an_element()
            0
        """
        return 0

    def _repr_(self):
        r"""
        TESTS::

            sage: SetsWithGrading().example() # indirect example
            Non negative integers
        """
        return "Non negative integers"

    def graded_component(self, grade):
        r"""
        Returns the component with grade ``grade``.

        EXAMPLES::

            sage: N = SetsWithGrading().example()
            sage: N.graded_component(65)
            {65}
        """
        return FiniteEnumeratedSet([grade])

    def grading(self, elt):
        r"""
        Returns the grade of ``elt``.

        EXAMPLES::

            sage: N = SetsWithGrading().example()
            sage: N.grading(10)
            10
        """
        return elt

Example = NonNegativeIntegers
