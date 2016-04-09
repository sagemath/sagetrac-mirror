r"""
Examples of finite Weyl groups
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.all import FiniteWeylGroups
from sage.structure.unique_representation import UniqueRepresentation

class SymmetricGroup(UniqueRepresentation, Parent):
    r"""
    An example of finite Weyl group: the symmetric group, with
    elements in list notation.

    The purpose of this class is to provide a minimal template for
    implementing finite Weyl groups. See
    :class:`~sage.groups.perm_gps.permgroup_named.SymmetricGroup` for
    a full featured and optimized implementation.

    EXAMPLES::

        sage: S = FiniteWeylGroups().example()
        sage: S
        The symmetric group on {0, ..., 3}
        sage: S.category()
        Category of finite weyl groups

    The elements of this group are permutations of the set `\{0,\ldots,3\}`::

        sage: S.one()
        (0, 1, 2, 3)
        sage: S.an_element()
        (1, 2, 3, 0)

    The group itself is generated by the elementary transpositions::

        sage: S.simple_reflections()
        Finite family {0: (1, 0, 2, 3), 1: (0, 2, 1, 3), 2: (0, 1, 3, 2)}

    Only the following basic operations are implemented:

    - :meth:`.one`
    - :meth:`.product`
    - :meth:`.simple_reflection`
    - :meth:`.Element.has_right_descent`.

    All the other usual Weyl group operations are inherited from the
    categories::

        sage: S.cardinality()
        24
        sage: S.long_element()
        (3, 2, 1, 0)
        sage: S.cayley_graph(side = "left").plot()
        Graphics object consisting of 120 graphics primitives

    Alternatively, one could have implemented
    :meth:`sage.categories.coxeter_groups.CoxeterGroups.ElementMethods.apply_simple_reflection`
    instead of :meth:`.simple_reflection` and :meth:`.product`. See
    ``CoxeterGroups().example()``.

    TESTS::

        sage: TestSuite(S).run(verbose=True)
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_has_descent() . . . pass
        running ._test_inverse() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_reduced_word() . . . pass
        running ._test_simple_projections() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_well_generated() . . . pass
    """

    def __init__(self, n = 4):
        """
        EXAMPLES::

            sage: S = sage.categories.examples.finite_weyl_groups.SymmetricGroup(4)
            sage: S == FiniteWeylGroups().example(4)
            True
        """
        Parent.__init__(self, category = FiniteWeylGroups())
        self.n = n

    def _repr_(self):
        """
        EXAMPLES::

            sage: FiniteWeylGroups().example()
            The symmetric group on {0, ..., 3}

        """
        return "The symmetric group on {0, ..., %s}"%(self.n-1)

    @cached_method
    def one(self):
        """
        Implements :meth:`Monoids.ParentMethods.one`.

        EXAMPLES::

            sage: FiniteWeylGroups().example().one()
            (0, 1, 2, 3)
        """
        return self(tuple(range(self.n)))

    def index_set(self):
        """
        Implements :meth:`CoxeterGroups.ParentMethods.index_set`.

        EXAMPLES::

            sage: FiniteWeylGroups().example().index_set()
            [0, 1, 2]
        """
        return range(self.n-1)

    def simple_reflection(self, i):
        """
        Implements :meth:`CoxeterGroups.ParentMethods.simple_reflection`
        by returning the transposition `(i, i+1)`.

        EXAMPLES::

            sage: FiniteWeylGroups().example().simple_reflection(2)
            (0, 1, 3, 2)
        """
        assert i in self.index_set()
        return self(tuple(range(i)+[i+1,i]+range(i+2,self.n)))

    def product(self, x, y):
        """
        Implements :meth:`Semigroups.ParentMethods.product`.

        EXAMPLES::

            sage: s = FiniteWeylGroups().example().simple_reflections()
            sage: s[1] * s[2]
            (0, 2, 3, 1)
            """
        assert x in self
        assert y in self
        return self(tuple(x.value[i] for i in y.value))

    def degrees(self):
        """
        Return the degrees of ``self``.

        EXAMPLES::

            sage: W = FiniteWeylGroups().example()
            sage: W.degrees()
            (2, 3, 4)

        TESTS::

            sage: W = FiniteWeylGroups().example()
            sage: prod(W.degrees()) == len(W)
            True
        """
        return tuple(range(2, self.n + 1))

    class Element(ElementWrapper):

        def has_right_descent(self, i):
            """
            Implements :meth:`CoxeterGroups.ElementMethods.has_right_descent`.

            EXAMPLES::

                sage: S = FiniteWeylGroups().example()
                sage: s = S.simple_reflections()
                sage: (s[1] * s[2]).has_descent(2)
                True
                sage: S._test_has_descent()
            """
            return (self.value[i] > self.value[i+1])


Example = SymmetricGroup
