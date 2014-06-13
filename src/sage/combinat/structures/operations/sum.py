# -*- coding: utf-8 -*-
"""
The sum of classes of combinatorial structures.

Let `F` and `G` be both disjoint classes.

We denote `H = F + G` the addition of `F` and `G`.
"""
from sage.combinat.structures.operations import _Operations
from sage.combinat.structures import Structures


class Sum(_Operations):
    """
    The sum of classes of combinatorial structures.

    Let `F` and `G` be both disjoint classes.

    We denote `H = F + G` the addition of `F` and `G`.
    """

    def _repr_(self):
        """
        TESTS::

            sage: B = BinaryTrees()
            sage: C = Compositions()
            sage: B + C
            Sum of structures : `Binary trees`, `Compositions of non-negative integers`
        """
        return "Sum of structures : `" + "`, `".join(map(repr, self._structures)) + "`"

    def generating_series(self, variable="x"):
        """
        The generating series `h` of `H`: the sum of `F` and `G` is
        defined by the sum of its generating series:

        MATH::

            h(t) = f(t) + g(t)

        """
        return sum(map(lambda F: F.generating_series(variable), self._structures))

    def _element_constructor_(self, F, *args, **options):
        """
        The element constructor use the original one and change the parent

        TESTS::

            sage: B = BinaryTrees()
            sage: C = Compositions()
            sage: BpC = B + C
            sage: BpC.graded_component(2).list()[0].parent() is BpC
            True
            sage: BpC.graded_component(2).list()[-1].parent() is BpC
            True
        """
        obj = F._element_constructor_(*args, **options)
        obj._set_parent(self)
        return obj

    class GradedComponent(Structures.GradedComponent):

        def __iter__(self):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: C = Compositions()
                sage: BpC = B + C
                sage: BpC.graded_component(2).list()
                [[., [., .]], [[., .], .], [1, 1], [2]]

            """
            k = self.grading()
            for F in self.ambient()._structures:
                for obj in F.graded_component(k):
                    yield self._element_constructor_(F, obj)