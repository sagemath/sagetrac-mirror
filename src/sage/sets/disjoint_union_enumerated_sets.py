"""
Disjoint union of enumerated sets

AUTHORS:

- Florent Hivert (2009-07/09): initial implementation.
- Florent Hivert (2010-03): classcall related stuff.
- Florent Hivert (2010-12): fixed facade element construction.
"""
#****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#****************************************************************************

# from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.sets.family import Family
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.rings.infinity import Infinity
from sage.misc.all import cached_method, lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation

class DisjointUnionEnumeratedSets(UniqueRepresentation, Parent):
    """
    A class for disjoint unions of enumerated sets.

    INPUT:

     - ``family``  -- a list (or iterable or family) of enumerated sets
     - ``keepkey`` -- a boolean
     - ``facade``  -- a boolean

    This models the enumerated set obtained by concatenating together
    the specified ordered sets. The later are supposed to be pairwise
    disjoint; otherwise, a multiset is created.

    The argument ``family`` can be a list, a tuple, a dictionary, or a
    family. If it is not a family it is first converted into a family
    (see :func:`sage.sets.family.Family`).

    Experimental options:

    By default, there is no way to tell from which set of the union an
    element is generated. The option ``keepkey=True`` keeps track of
    those by returning pairs ``(key, el)`` where ``key`` is the index
    of the set to which ``el`` belongs. When this option is specified,
    the enumerated sets need not be disjoint anymore.

    With the option ``facade=False`` the elements are wrapped in an
    object whose parent is the disjoint union itself. The wrapped
    object can then be recovered using the 'value' attribute.

    The two options can be combined.

    The names of those options is imperfect, and subject to change in
    future versions. Feedback welcome.

    EXAMPLES:

    The input can be a list or a tuple of FiniteEnumeratedSets::

        sage: U1 = DisjointUnionEnumeratedSets((
        ....:       FiniteEnumeratedSet([1,2,3]),
        ....:       FiniteEnumeratedSet([4,5,6])))
        sage: U1
        Disjoint union of Family ({1, 2, 3}, {4, 5, 6})
        sage: U1.list()
        [1, 2, 3, 4, 5, 6]
        sage: U1.cardinality()
        6

    The input can also be a dictionary::

        sage: U2 = DisjointUnionEnumeratedSets({1: FiniteEnumeratedSet([1,2,3]),
        ....:                                   2: FiniteEnumeratedSet([4,5,6])})
        sage: U2
        Disjoint union of Finite family {1: {1, 2, 3}, 2: {4, 5, 6}}
        sage: U2.list()
        [1, 2, 3, 4, 5, 6]
        sage: U2.cardinality()
        6

    However in that case the enumeration order is not specified.

    In general the input can be any family::

        sage: U3 = DisjointUnionEnumeratedSets(
        ....:     Family([2,3,4], Permutations, lazy=True))
        sage: U3
        Disjoint union of Lazy family (<class 'sage.combinat.permutation.Permutations'>(i))_{i in [2, 3, 4]}
        sage: U3.cardinality()
        32
        sage: it = iter(U3)
        sage: [next(it), next(it), next(it), next(it), next(it), next(it)]
        [[1, 2], [2, 1], [1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1]]
        sage: U3.unrank(18)
        [2, 4, 1, 3]

    This allows for infinite unions::

        sage: U4 = DisjointUnionEnumeratedSets(
        ....:     Family(NonNegativeIntegers(), Permutations))
        sage: U4
        Disjoint union of Lazy family (<class 'sage.combinat.permutation.Permutations'>(i))_{i in Non negative integers}
        sage: U4.cardinality()
        +Infinity
        sage: it = iter(U4)
        sage: [next(it), next(it), next(it), next(it), next(it), next(it)]
        [[], [1], [1, 2], [2, 1], [1, 2, 3], [1, 3, 2]]
        sage: U4.unrank(18)
        [2, 3, 1, 4]

    .. warning ::

        Beware that some of the operations assume in that case that infinitely
        many of the enumerated sets are non empty.


    .. rubric:: Experimental options

    We demonstrate the ``keepkey`` option::

        sage: Ukeep = DisjointUnionEnumeratedSets(
        ....:            Family(range(4), Permutations), keepkey=True)
        sage: it = iter(Ukeep)
        sage: [next(it) for i in range(6)]
        [(0, []), (1, [1]), (2, [1, 2]), (2, [2, 1]), (3, [1, 2, 3]), (3, [1, 3, 2])]
        sage: type(next(it)[1])
        <class 'sage.combinat.permutation.StandardPermutations_n_with_category.element_class'>

    We now demonstrate the ``facade`` option::

        sage: UNoFacade = DisjointUnionEnumeratedSets(
        ....:                Family(range(4), Permutations), facade=False)
        sage: it = iter(UNoFacade)
        sage: [next(it) for i in range(6)]
        [[], [1], [1, 2], [2, 1], [1, 2, 3], [1, 3, 2]]
        sage: el = next(it); el
        [2, 1, 3]
        sage: type(el)
        <type 'sage.structure.element_wrapper.ElementWrapper'>
        sage: el.parent() == UNoFacade
        True
        sage: elv = el.value; elv
        [2, 1, 3]
        sage: type(elv)
        <class 'sage.combinat.permutation.StandardPermutations_n_with_category.element_class'>

    Possible extensions: the current enumeration order is not suitable
    for unions of infinite enumerated sets (except possibly for the
    last one). One could add options to specify alternative enumeration
    orders (anti-diagonal, round robin, ...) to handle this case.


    .. rubric:: Inheriting from ``DisjointUnionEnumeratedSets``

    There are two different use cases for inheriting from
    :class:`DisjointUnionEnumeratedSets`: writing a parent which
    happens to be a disjoint union of some known parents, or writing
    generic disjoint unions for some particular classes of
    :class:`sage.categories.enumerated_sets.EnumeratedSets`.

      - In the first use case, the input of the ``__init__`` method is
        most likely different from that of
        :class:`DisjointUnionEnumeratedSets`. Then, one simply
        writes the ``__init__`` method as usual::

            sage: class MyUnion(DisjointUnionEnumeratedSets):
            ....:   def __init__(self):
            ....:       DisjointUnionEnumeratedSets.__init__(self,
            ....:            Family([1,2], Permutations))
            sage: pp = MyUnion()
            sage: pp.list()
            [[1], [1, 2], [2, 1]]

        In case the :meth:`__init__` method takes optional arguments,
        or does some normalization on them, a specific method
        ``__classcall_private__`` is required (see the
        documentation of :class:`UniqueRepresentation`).

      - In the second use case, the input of the ``__init__`` method
        is the same as that of :class:`DisjointUnionEnumeratedSets`;
        one therefore wants to inherit the :meth:`__classcall_private__`
        method as well, which can be achieved as follows::

            sage: class UnionOfSpecialSets(DisjointUnionEnumeratedSets):
            ....:  __classcall_private__ = staticmethod(DisjointUnionEnumeratedSets.__classcall_private__)
            ...
            sage: psp = UnionOfSpecialSets(Family([1,2], Permutations))
            sage: psp.list()
            [[1], [1, 2], [2, 1]]

    TESTS::

        sage: TestSuite(U1).run()
        sage: TestSuite(U2).run()
        sage: TestSuite(U3).run()
        sage: TestSuite(U4).run()
        doctest:...: UserWarning: Disjoint union of Lazy family (<class 'sage.combinat.permutation.Permutations'>(i))_{i in Non negative integers} is an infinite union
        The default implementation of __contains__ can loop forever. Please overload it.
        sage: TestSuite(Ukeep).run()
        sage: TestSuite(UNoFacade).run()

    The following three lines are required for the pickling tests,
    because the classes ``MyUnion`` and ``UnionOfSpecialSets`` have
    been defined interactively::

        sage: import __main__
        sage: __main__.MyUnion = MyUnion
        sage: __main__.UnionOfSpecialSets = UnionOfSpecialSets

        sage: TestSuite(pp).run()
        sage: TestSuite(psp).run()

    """

    @staticmethod
    def __classcall_private__(cls, fam, facade=True,
                              keepkey=False, category=None):
        """
        Normalization of arguments; see :class:`UniqueRepresentation`.

        TESTS:

        We check that disjoint unions have unique representation::

            sage: U1 = DisjointUnionEnumeratedSets({1: FiniteEnumeratedSet([1,2,3]),
            ....:                                   2: FiniteEnumeratedSet([4,5,6])})
            sage: U2 = DisjointUnionEnumeratedSets({1: FiniteEnumeratedSet([1,2,3]),
            ....:                                   2: FiniteEnumeratedSet([4,5,6])})
            sage: U1 == U2
            True
            sage: U1 is U2        # indirect doctest
            True
            sage: U3 = DisjointUnionEnumeratedSets({1: FiniteEnumeratedSet([1,2,3]),
            ....:                                   2: FiniteEnumeratedSet([4,5])})
            sage: U1 == U3
            False
        """
        # facade  = options.pop('facade', True);
        # keepkey = options.pop('keepkey', False);
        assert(isinstance(facade,  bool))
        assert(isinstance(keepkey, bool))
        return super(DisjointUnionEnumeratedSets, cls).__classcall__(
            cls, Family(fam),
            facade = facade, keepkey = keepkey, category=category)


    def __init__(self, family, facade=True, keepkey=False, category = None):
        """
        TESTS::

            sage: U = DisjointUnionEnumeratedSets({1: FiniteEnumeratedSet([1,2,3]),
            ....:                                  2: FiniteEnumeratedSet([4,5,6])})
            sage: TestSuite(U).run()
        """
        self._family = family
        self._facade  = facade
        self._keepkey = keepkey
        if self._is_category_initialized():
            return
        if category is None:
            # try to guess if the result is infinite or not.
            if self._family in InfiniteEnumeratedSets():
                category = InfiniteEnumeratedSets()
            elif self._family.last().cardinality() == Infinity:
                category = InfiniteEnumeratedSets()
            else:
                category = FiniteEnumeratedSets()
        Parent.__init__(self, category = category)

    def _repr_(self):
        """
        TESTS::

            sage: U = DisjointUnionEnumeratedSets({1: FiniteEnumeratedSet([1,2,3]),
            ....:                                  2: FiniteEnumeratedSet([4,5,6])})
            sage: U  # indirect doctest
            Disjoint union of Finite family {1: {1, 2, 3}, 2: {4, 5, 6}}
        """
        return "Disjoint union of %s"%self._family


    def _is_a(self, x):
        """
        Check if a Sage object belongs to self. This methods is a helper for
        :meth:`__contains__` and the constructor :meth:`_element_constructor_`.

        EXAMPLES::

            sage: U4 = DisjointUnionEnumeratedSets(
            ....:          Family(NonNegativeIntegers(), Compositions))
            sage: U4._is_a(Composition([3,2,1,1]))
            doctest:...: UserWarning: Disjoint union of Lazy family (<class 'sage.combinat.composition.Compositions'>(i))_{i in Non negative integers} is an infinite union
            The default implementation of __contains__ can loop forever. Please overload it.
            True
        """
        if self._keepkey:
            return (isinstance(x, tuple) and
                    x[0] in self._family.keys() and
                    x[1] in self._family[x[0]])
        else:
            from warnings import warn
            if self._family.cardinality() == Infinity:
                warn("%s is an infinite union\nThe default implementation of __contains__ can loop forever. Please overload it."%(self))
            return any(x in a for a in self._family)


    def __contains__(self, x):
        """
        .. warning::

        If ``self`` is an infinite union and if the answer is
        logically False, this will loop forever and never answer
        ``False``. Therefore, a warning is issued.

        EXAMPLES::

            sage: U4 = DisjointUnionEnumeratedSets(
            ....:          Family(NonNegativeIntegers(), Partitions))
            sage: Partition([]) in U4
            doctest:...: UserWarning: Disjoint union of Lazy family (<class 'sage.combinat.partition.Partitions'>(i))_{i in Non negative integers} is an infinite union
            The default implementation of __contains__ can loop forever. Please overload it.
            True

        Note: one has to use a different family from the previous one in this
        file otherwise the warning is not re-issued::

            sage: Partition([3,2,1,1]) in U4
            True

        The following call will loop forever::

            sage: 2 in U4 # not tested, loop forever
        """
        if self._facade:
            return self._is_a(x)
        else:
            if isinstance(x, self.element_class):
                return True
            else:
                return self._is_a(x)

    def __iter__(self):
        """
        TESTS::

            sage: U4 = DisjointUnionEnumeratedSets(
            ....:          Family(NonNegativeIntegers(), Permutations))
            sage: it = iter(U4)
            sage: [next(it), next(it), next(it), next(it), next(it), next(it)]
            [[], [1], [1, 2], [2, 1], [1, 2, 3], [1, 3, 2]]

            sage: U4 = DisjointUnionEnumeratedSets(
            ....:          Family(NonNegativeIntegers(), Permutations),
            ....:          keepkey=True, facade=False)
            sage: it = iter(U4)
            sage: [next(it), next(it), next(it), next(it), next(it), next(it)]
            [(0, []), (1, [1]), (2, [1, 2]), (2, [2, 1]), (3, [1, 2, 3]), (3, [1, 3, 2])]
            sage: el = next(it); el.parent() == U4
            True
            sage: el.value == (3, Permutation([2,1,3]))
            True
        """
        for k in self._family.keys():
            for el in self._family[k]:
                if self._keepkey:
                    el = (k, el)
                if self._facade:
                    yield el
                else:
                    yield self.element_class(self, el) # Bypass correctness tests

    def an_element(self):
        """
        Returns an element of this disjoint union, as per :meth:`Sets.ParentMethods.an_element`.

        EXAMPLES::

            sage: U4 = DisjointUnionEnumeratedSets(
            ....:          Family([3, 5, 7], Permutations))
            sage: U4.an_element()
            [1, 2, 3]
        """
        return self._an_element_from_iterator()

    @cached_method
    def cardinality(self):
        """
        Returns the cardinality of this disjoint union.

        EXAMPLES:

        For finite disjoint unions, the cardinality is computed by
        summing the cardinalities of the enumerated sets::

            sage: U = DisjointUnionEnumeratedSets(Family([0,1,2,3], Permutations))
            sage: U.cardinality()
            10

        For infinite disjoint unions, this makes the assumption that
        the result is infinite::

            sage: U = DisjointUnionEnumeratedSets(
            ....:         Family(NonNegativeIntegers(), Permutations))
            sage: U.cardinality()
            +Infinity

        .. warning::

            as pointed out in the main documentation, it is
            possible to construct examples where this is incorrect::

                sage: U = DisjointUnionEnumeratedSets(
                ....:         Family(NonNegativeIntegers(), lambda x: []))
                sage: U.cardinality()  # Should be 0!
                +Infinity

        """
        if self._family.cardinality() == Infinity:
            return Infinity
        return sum(set.cardinality() for set in self._family)

    @lazy_attribute
    def _element_constructor_(self):
        """
        TESTS::

            sage: U = DisjointUnionEnumeratedSets(
            ....:          Family([1,2,3], Partitions), facade=False)
            sage: U._element_constructor_
            <bound method DisjointUnionEnumeratedSets_with_category._element_constructor_default of Disjoint union of Finite family {1: Partitions of the integer 1, 2: Partitions of the integer 2, 3: Partitions of the integer 3}>
            sage: U = DisjointUnionEnumeratedSets(
            ....:          Family([1,2,3], Partitions), facade=True)
            sage: U._element_constructor_
            Traceback (most recent call last):
            ...
            AttributeError: 'DisjointUnionEnumeratedSets_with_category' object has no attribute '_element_constructor_'
        """
        if not self._facade:
            return self._element_constructor_default
        else:
            return NotImplemented

    def _element_constructor_default(self, el):
        r"""
        TESTS::

            sage: U = DisjointUnionEnumeratedSets(
            ....:         Family([1,2,3], Partitions), facade=False)
            sage: U([1])       # indirect doctest
            [1]
            sage: U([2,1])     # indirect doctest
            [2, 1]
            sage: U([1,3,2])   # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Value [1, 3, 2] does not belong to Disjoint union of Finite family {1: Partitions of the integer 1, 2: Partitions of the integer 2, 3: Partitions of the integer 3}

            sage: U = DisjointUnionEnumeratedSets(
            ....:          Family([1,2,3], Partitions), keepkey=True, facade=False)
            sage: U((1, [1]))    # indirect doctest
            (1, [1])
            sage: U((3,[2,1]))   # indirect doctest
            (3, [2, 1])
            sage: U((4,[2,1]))   # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Value (4, [2, 1]) does not belong to Disjoint union of Finite family {1: Partitions of the integer 1, 2: Partitions of the integer 2, 3: Partitions of the integer 3}
        """
        if isinstance(el, self.element_class):
            el = el.value
        if self._is_a(el):
            return self.element_class(self, el)
        else:
            raise ValueError("Value %s does not belong to %s"%(el, self))

    @lazy_attribute
    def Element(self):
        """
        TESTS::

            sage: U = DisjointUnionEnumeratedSets(
            ....:          Family([1,2,3], Partitions), facade=False)
            sage: U.Element
            <type 'sage.structure.element_wrapper.ElementWrapper'>
            sage: U = DisjointUnionEnumeratedSets(
            ....:          Family([1,2,3], Partitions), facade=True)
            sage: U.Element
            Traceback (most recent call last):
            ...
            AttributeError: 'DisjointUnionEnumeratedSets_with_category' object has no attribute 'Element'
        """
        if not self._facade:
            return ElementWrapper
        else:
            return NotImplemented
