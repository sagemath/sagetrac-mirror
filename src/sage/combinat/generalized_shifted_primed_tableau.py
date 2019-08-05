# -*- coding: utf-8 -*-
r"""
Generalized shifted primed tableaux

AUTHORS:

- Kirill Paramonov (2017-08-18): initial implementation
- Chaman Agrawal (2017-07-24): Modified shifted primed tableaux to include
  primed entries on main diagonal.
"""

from __future__ import print_function, absolute_import, division
from six import add_metaclass

from sage.combinat.partition import Partition, Partitions, _Partitions, OrderedPartitions
from sage.combinat.partitions import ZS1_iterator
from sage.combinat.tableau import Tableaux
from sage.combinat.skew_partition import SkewPartition
from sage.combinat.integer_vector import IntegerVectors
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute

from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.regular_crystals import RegularCrystals
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.sets_cat import Sets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.shifted_primed_tableau import ShiftedPrimedTableau, _add_strip
from sage.combinat.shifted_primed_tableau import PrimedEntry


@add_metaclass(InheritComparisonClasscallMetaclass)
class GeneralizedShiftedPrimedTableau(ShiftedPrimedTableau):
    r"""
    A generalized shifted primed tableau.
    
    A generalized primed tableau is a tableau of shifted shape in the alphabet
    `X' = \{1' < 1 < 2' < 2 < \cdots < n' < n\}` such that:

    1. the entries are weakly increasing along rows and columns;
    2. a row cannot have two repeated primed elements, and a column
       cannot have two repeated non-primed elements.

    ..NOTE::

        This class is a variation of the class
        :class:`~sage.combinat.shifted_primed_tableau.ShiftedPrimedTableau`.
        The main difference between
        :class:`~sage.combinat.shifted_primed_tableau.ShiftedPrimedTableau` and
        :class:`~sage.combinat.generalized_shifted_primed_tableau.
        GeneralizedShiftedPrimedTableau` is that the first one does not allow
        primed entries on the main diagonal while second does allow.

    Skew shape of the generalized shifted primed tableaux is specified either
    with an optional argument ``skew`` or with ``None`` entries.

    EXAMPLES::

        sage: T = GeneralizedShiftedPrimedTableaux([4,2])
        sage: T([[1,"2'","3'",3], [2,"3'"]])[1]
        (2, 3')
        sage: t = GeneralizedShiftedPrimedTableau([[1,"2p",2.5,3], [2,2.5]])
        sage: t[1]
        (2, 3')
        sage: GeneralizedShiftedPrimedTableau([["2p",2,3], ["2p","3p"], [2]], 
        ....:                                 skew=[2,1])
        [(None, None, 2', 2, 3), (None, 2', 3'), (2,)]
        sage: GeneralizedShiftedPrimedTableau([[None,None,"2p"], [None,"2p"]])
        [(None, None, 2'), (None, 2')]

    TESTS:

        sage: t = GeneralizedShiftedPrimedTableau([[1,2,2.5,3], [2,2.5]])
        Traceback (most recent call last):
        ...
        ValueError: [[1, 2, 2.50000000000000, 3], [2, 2.50000000000000]]
         is not an element of Generalized Shifted Primed Tableaux
    """
    @staticmethod
    def __classcall_private__(cls, T, skew=None):
        r"""
        Ensure that a generalized shifted tableau is only ever constructed as
        an ``element_class`` call of an appropriate parent.

        EXAMPLES::

            sage: data = [[1,"2'","2",3], [2,"3'"]]
            sage: t = GeneralizedShiftedPrimedTableau(data)
            sage: T = GeneralizedShiftedPrimedTableaux(shape=[4,2],
            ....:                                      weight=(1,3,2))
            sage: t == T(data)
            True
            sage: S = GeneralizedShiftedPrimedTableaux(shape=[4,2])
            sage: t == S(data)
            True
            sage: t = GeneralizedShiftedPrimedTableau([["2p",2,3],["2p"]],
            ....:                                     skew=[2,1])
            sage: t.parent()
            Generalized Shifted Primed Tableaux skewed by [2, 1]
            sage: s = GeneralizedShiftedPrimedTableau([[None, None,"2p",2,3],
            ....:                                     [None,"2p"]])
            sage: s.parent()
            Generalized Shifted Primed Tableaux skewed by [2, 1]

        TESTS:

            sage: GeneralizedShiftedPrimedTableau([])
            []
            sage: GeneralizedShiftedPrimedTableau([tuple([])])
            []
        """
        if isinstance(T, GeneralizedShiftedPrimedTableau) and T._skew == skew:
            return T

        skew_ = Partition([row.count(None) for row in T])
        if skew_:
            if skew and Partition(skew) != skew_:
                raise ValueError("skew shape does not agree with None entries")
            skew = skew_
        return GeneralizedShiftedPrimedTableaux(skew=skew)(T)

    @staticmethod
    def _preprocess(T, skew=None):
        r"""
        Preprocessing list ``T`` to initialize the tableau.
        The output is a list of rows as tuples, with explicit
        ``None``'s to indicate the skew shape, and entries being
        ``PrimedEntry``s.

        Trailing empty rows are removed.

        TESTS:

            sage: GeneralizedShiftedPrimedTableau._preprocess([["2'","3p",3.5]]
            ....:                                             , skew=[1])
            [(None, 2', 3', 4')]
            sage: GeneralizedShiftedPrimedTableau._preprocess([[None]],
            ....:                                             skew=[1])
            [(None,)]
            sage: GeneralizedShiftedPrimedTableau._preprocess([], skew=[2,1])
            [(None, None), (None,)]
            sage: GeneralizedShiftedPrimedTableau._preprocess([], skew=[])
            []
        """
        if isinstance(T, GeneralizedShiftedPrimedTableau):
            return T
        # Preprocessing list t for primes and other symbols
        T = [[PrimedEntry(entry) for entry in row if entry is not None]
             for row in T]
        while len(T) > 0 and len(T[-1]) == 0:
            T = T[:-1]
        row_min = min(len(skew), len(T)) if skew else 0
        T_ = [(None,)*skew[i] + tuple(T[i]) for i in range(row_min)]

        if row_min < len(T):
            T_ += [tuple(T[i]) for i in range(row_min, len(T))]
        elif skew:
            T_ += [(None,)*skew[i] for i in range(row_min, len(skew))]
        return T_

    def check(self):
        r"""
        Check that ``self`` is a valid primed tableau.

        EXAMPLES::

            sage: T = GeneralizedShiftedPrimedTableaux([4,2])
            sage: t = T([[1,'2p',2,2],[2,'3p']])
            sage: t.check()
            sage: s = GeneralizedShiftedPrimedTableau([["2p",2,3],["2p"],[2]],
            ....:                                     skew=[2,1])
            sage: s.check()
            sage: t = T([['1p','1p',2,2],[2,'3p']])
            Traceback (most recent call last):
            ...
            ValueError: [['1p', '1p', 2, 2], [2, '3p']] is not an element of
            Generalized Shifted Primed Tableaux of shape [4, 2] and maximum
            entry 6
        """
        if not self.parent()._contains_tableau(self):
            raise ValueError("{} is not an element of Generalized Shifted "
                             "Primed Tableaux".format(self))

    def __eq__(self, other):
        r"""
        Check whether ``self`` is equal to ``other``.

        INPUT:

        - ``other`` -- the element that ``self`` is compared to

        EXAMPLES::

            sage: t = GeneralizedShiftedPrimedTableau([[1,"2p"]])
            sage: t == GeneralizedShiftedPrimedTableaux([2])([[1,3/2]])
            True
            sage: s = GeneralizedShiftedPrimedTableau([["2p",3]], skew=[1])
            sage: s == [[None, "2p", 3]]
            True
        """
        if isinstance(other, GeneralizedShiftedPrimedTableau):
            return self._skew == other._skew and list(self) == list(other)
        try:
            Tab = GeneralizedShiftedPrimedTableau(other)
        except (ValueError, TypeError):
            return False
        return self._skew == Tab._skew and list(self) == list(Tab)

    def restrict(self, n):
        r"""
        Return the restriction of the shifted tableau to all
        the numbers less than or equal to ``n``.

        EXAMPLES::

            sage: t = GeneralizedShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.restrict(2).pp()
            1  2' 2  2
               2
            sage: t.restrict("2p").pp()
            1  2'
            sage: s = GeneralizedShiftedPrimedTableau([["2p",2,3],["2p"]],
            ....:                                     skew=[2,1])
            sage: s.restrict(2).pp()
            .  .  2' 2
               .  2'
            sage: s.restrict(1.5).pp()
            .  .  2'
               .  2'
        """
        t = self[:]
        n = PrimedEntry(n)
        return GeneralizedShiftedPrimedTableau([z for z in [[y for y in x if y is not None and y <= n]
                                                 for x in t] if z], skew=self._skew)


class GeneralizedShiftedPrimedTableaux(UniqueRepresentation, Parent):
    r"""
    Return the combinatorial class of generalized shifted primed tableaux
    subject to the constraints given by the arguments.

    A generalized primed tableau is a tableau of shifted shape on the alphabet
    `X' = \{1' < 1 < 2' < 2 < \cdots < n' < n\}` such that:

    1. the entries are weakly increasing along rows and columns

    2. a row cannot have two repeated primed entries, and a column
       cannot have two repeated non-primed entries.

    INPUT:

    Valid optional keywords:

    - ``shape`` -- the (outer skew) shape of tableaux

    - ``weight`` -- the weight of tableaux

    - ``max_entry`` -- the maximum entry of tableaux

    - ``skew`` -- the inner skew shape of tableaux

    The weight of a tableau is defined to be the vector with `i`-th
    component equal to the number of entries `i` and `i'` in the tableau.
    The sum of the coordinates in the weight vector must be equal to the
    number of entries in the partition.

    The ``shape`` and ``skew`` must be strictly decreasing partitions.

    EXAMPLES::

        sage: SPT = GeneralizedShiftedPrimedTableaux(weight=(1,2,2),
        ....:                                        shape=[3,2]); SPT
        Generalized Shifted Primed Tableaux of weight (1, 2, 2) and shape [3, 2]
        sage: SPT.list()
        [[(1, 2, 2), (3, 3)],
         [(1, 2, 2), (3', 3)],
         [(1, 2', 3'), (2, 3)],
         [(1, 2', 3'), (2, 3')],
         [(1, 2', 3'), (2', 3)],
         [(1, 2', 3'), (2', 3')],
         [(1, 2', 2), (3, 3)],
         [(1, 2', 2), (3', 3)],
         [(1', 2, 2), (3, 3)],
         [(1', 2, 2), (3', 3)],
         [(1', 2', 3'), (2, 3)],
         [(1', 2', 3'), (2, 3')],
         [(1', 2', 3'), (2', 3)],
         [(1', 2', 3'), (2', 3')],
         [(1', 2', 2), (3, 3)],
         [(1', 2', 2), (3', 3)]]
        sage: SPT = GeneralizedShiftedPrimedTableaux(weight=(1,2)); SPT
        Generalized Shifted Primed Tableaux of weight (1, 2)
        sage: list(SPT)
        [[(1, 2, 2)],
         [(1, 2', 2)],
         [(1', 2, 2)],
         [(1', 2', 2)],
         [(1, 2'), (2,)],
         [(1, 2'), (2',)],
         [(1', 2'), (2,)],
         [(1', 2'), (2',)]]
        sage: SPT = GeneralizedShiftedPrimedTableaux([3,2], max_entry = 2); SPT
        Generalized Shifted Primed Tableaux of shape [3, 2] and maximum entry 2
        sage: list(SPT)
        [[(1, 1, 1), (2, 2)],
         [(1, 1, 1), (2', 2)],
         [(1', 1, 1), (2, 2)],
         [(1', 1, 1), (2', 2)],
         [(1, 1, 2'), (2, 2)],
         [(1, 1, 2'), (2', 2)],
         [(1', 1, 2'), (2, 2)],
         [(1', 1, 2'), (2', 2)]]

    TESTS:

        sage: [(1,'2p',2,2),(2,'3p')] in GeneralizedShiftedPrimedTableaux()
        True
        sage: [(1,1),(2,2)] in GeneralizedShiftedPrimedTableaux()
        False
        sage: [] in GeneralizedShiftedPrimedTableaux()
        True
    """
    Element = GeneralizedShiftedPrimedTableau
    options = Tableaux.options

    @staticmethod
    def __classcall_private__(cls, shape=None, weight=None,
                              max_entry=None, skew=None):
        r"""
        Normalize and process input to return the correct parent and
        ensure a unique representation.

        TESTS:

            sage: GeneralizedShiftedPrimedTableaux([])
            Generalized Shifted Primed Tableaux of shape [] and maximum entry 0
            sage: GeneralizedShiftedPrimedTableaux(3)
            Traceback (most recent call last):
            ...
            ValueError: invalid shape argument
            sage: GeneralizedShiftedPrimedTableaux(weight=(2,2,2), shape=[3,2])
            Traceback (most recent call last):
            ...
            ValueError: weight and shape are incompatible
            sage: GeneralizedShiftedPrimedTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: invalid shape argument
            sage: GeneralizedShiftedPrimedTableaux(weight=(2,2,2), max_entry=2)
            Traceback (most recent call last):
            ...
            ValueError: maximum entry is incompatible with the weight
            sage: GeneralizedShiftedPrimedTableaux(shape=[4,1],skew=[3,2])
            Traceback (most recent call last):
            ...
            ValueError: skew shape must be inside the given tableau shape

            sage: SPT1 = GeneralizedShiftedPrimedTableaux(weight=())
            sage: SPT2 = GeneralizedShiftedPrimedTableaux(weight=(0,0,0))
            sage: SPT1 is SPT2
            True
        """
        if skew is not None:
            try:
                skew = Partition(skew)
            except ValueError:
                raise ValueError("invalid skew argument")
            if not all(skew[i] > skew[i+1] for i in range(len(skew)-1)):
                raise ValueError("skew shape must be a strict partition")

        if weight is not None:
            weight = tuple(weight)

        if shape is not None:
            if isinstance(shape, SkewPartition):
                skew = shape.inner()
                shape = shape.outer()
            try:
                shape = Partition(shape)
            except (ValueError, TypeError):
                raise ValueError("invalid shape argument")

            if not all(shape[i] > shape[i+1] for i in range(len(shape)-1)):
                raise ValueError("shape {} is not a strict partition".format(shape))

            if (skew is not None and not all(skew[i] <= shape[i]
                                             for i in range(len(skew)))):
                raise ValueError("skew shape must be inside the given tableau shape")

        if weight is not None:
            while weight and weight[-1] == 0:
                weight = weight[:-1]

        if max_entry is not None and weight is not None:
            if len(weight) > max_entry:
                raise ValueError("maximum entry is incompatible with the weight")

        if shape is None:
            if weight is None:
                if max_entry is not None:
                    raise ValueError("specify shape or weight argument")
                return GeneralizedShiftedPrimedTableaux_all(skew=skew)
            else:
                return GeneralizedShiftedPrimedTableaux_weight(weight, skew=skew)
        else:
            if weight is None:
                return GeneralizedShiftedPrimedTableaux_shape(shape, max_entry=max_entry, skew=skew)

            if (skew is not None and sum(shape) - sum(skew) != sum(weight)
                    or skew is None and sum(shape) != sum(weight)):
                raise ValueError("weight and shape are incompatible")
            return GeneralizedShiftedPrimedTableaux_weight_shape(weight, shape, skew=skew)

    def __init__(self, skew=None):
        r"""
        Initialize the parent class with given skew shape.

        TESTS:

            sage: SPT = GeneralizedShiftedPrimedTableaux(skew=[1])
            sage: TestSuite(SPT).run()  # known bug
        """
        self._skew = skew

    def _element_constructor_(self, T):
        r"""
        Construct an object from ``T`` as an element of generalized shifted
        primed tableaux, if possible.

        INPUT:

        - ``T`` -- data which can be interpreted as a primed tableau

        OUTPUT:

        - the corresponding primed tableau object

        EXAMPLES::

            sage: SPT = GeneralizedShiftedPrimedTableaux()
            sage: tab = SPT([[1,1,"2p"]]); tab
            [(1, 1, 2')]
            sage: tab.parent() is SPT
            True
            sage: tab = SPT([[1,1,2],[2,2]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1, 2], [2, 2]] is not an element of Generalized
            Shifted Primed Tableaux
            sage: SPT([[1,"2p","2p"]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, '2p', '2p']] is not an element of Generalized
            Shifted Primed Tableaux

            sage: SPT = GeneralizedShiftedPrimedTableaux(skew=[1])
            sage: SPT([["2p",2]])
            [(None, 2', 2)]

            sage: SPT = GeneralizedShiftedPrimedTableaux(weight=(2,1))
            sage: tab = SPT([[1,1,1.5]]); tab
            [(1, 1, 2')]
            sage: tab.parent() is SPT
            True

            sage: SPT = GeneralizedShiftedPrimedTableaux([3])
            sage: tab = SPT([[1,1,1.5]]); tab
            [(1, 1, 2')]
            sage: tab.parent() is SPT
            True
            sage: SPT([[1,1]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1]] is not an element of Generalized Shifted 
            Primed Tableaux of shape [3] and maximum entry 3
            sage: SPT = GeneralizedShiftedPrimedTableaux([3], weight=(2,1))
            sage: tab = SPT([[1,1,1.5]]); tab
            [(1, 1, 2')]
            sage: tab.parent() is SPT
            True
            sage: SPT([[1,1]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1]] is not an element of Generalized Shifted
            Primed Tableaux of weight (2, 1) and shape [3]
        """
        try:
            return self.element_class(self, T, skew=self._skew)
        except ValueError:
            raise ValueError("{} is not an element of {}".format(T, self))

    def _contains_tableau(self, T):
        r"""
        Check if ``self`` contains preprocessed tableau ``T``.

        TESTS:

            sage: Tabs = GeneralizedShiftedPrimedTableaux()
            sage: tab = GeneralizedShiftedPrimedTableau._preprocess(
            ....: [[1,"2p","3p","3p"]])
            sage: tab
            [(1, 2', 3', 3')]
            sage: Tabs._contains_tableau(tab)
            False
            sage: Tabs = GeneralizedShiftedPrimedTableaux(skew=[1])
            sage: tab = GeneralizedShiftedPrimedTableau._preprocess(
            ....: [["2p","3p",3]], skew=[1])
            sage: tab
            [(None, 2', 3', 3)]
            sage: Tabs._contains_tableau(tab)
            True
        """
        if not all(len(T[i]) > len(T[i+1]) for i in range(len(T)-1)):
            return False
        if self._skew is not None:
            skew = self._skew + [0]*(len(T)-len(self._skew))
        else:
            skew = [0] * len(T)
        for i, row in enumerate(T):
            if i > 0:
                if not all(val > T[i-1][j+1]
                           for j, val in enumerate(row)
                           if j+1 >= skew[i-1] and val.is_unprimed()):
                    return False
                if not all(val >= T[i-1][j+1]
                           for j, val in enumerate(row)
                           if j+1 >= skew[i-1] and val.is_primed()):
                    return False
            if not all(row[j] <= row[j+1]
                       for j in range(skew[i], len(row)-1)
                       if row[j].is_unprimed()):
                return False
            if not all(row[j] < row[j+1]
                       for j in range(skew[i], len(row)-1)
                       if row[j].is_primed()):
                return False
        return True


class GeneralizedShiftedPrimedTableaux_all(GeneralizedShiftedPrimedTableaux):
    r"""
    The class of all generalized shifted primed tableaux.
    """
    def __init__(self, skew=None):
        r"""
        Initialize the class of all generalized shifted tableaux.

        TESTS:

            sage: SPT = GeneralizedShiftedPrimedTableaux()
            sage: [[1,1.5],[2]] in SPT
            True
            sage: [[1,1.5],[1.5]] in SPT
            True
            sage: [[1,1],[1]] in SPT
            False
            sage: [[1,1],[2,2]] in SPT
            False
            sage: TestSuite(SPT).run()  # long time
        """
        if skew is None:
            Parent.__init__(self, category=InfiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=Sets().Infinite())
        GeneralizedShiftedPrimedTableaux.__init__(self, skew=skew)
        self._skew = skew

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS:

            sage: GeneralizedShiftedPrimedTableaux()
            Generalized Shifted Primed Tableaux
        """
        if self._skew is None:
            return "Generalized Shifted Primed Tableaux"
        return "Generalized Shifted Primed Tableaux skewed by {}".format(self._skew)

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: Tabs = GeneralizedShiftedPrimedTableaux()
            sage: Tabs[:5]
            [[], [(1,)], [(1',)], [(2,)], [(2',)]]
        """
        if self._skew is not None:
            raise NotImplementedError("skew tableau must be empty")
        yield self.element_class(self, [])  

        max_entry = 1
        while True:
            for size in range(1, max_entry+1):
                for shape in Partitions(size, max_slope=-1):
                    for weight in OrderedPartitions(size+max_entry-1,
                                                    k=max_entry):
                        weight = [weight[i]-1 for i in range(max_entry)]
                        weight[-1] += 1
                        for tab in GeneralizedShiftedPrimedTableaux(shape=shape,
                                                         weight=weight):
                            yield self.element_class(self, tab, check=False,
                                                     preprocessed=True)
            max_entry += 1


class GeneralizedShiftedPrimedTableaux_shape(GeneralizedShiftedPrimedTableaux):
    r"""
    Generalized Shifted primed tableaux of a fixed shape.

    EXAMPLES::

        sage: GeneralizedShiftedPrimedTableaux([4,3,1], max_entry=4)
        Generalized Shifted Primed Tableaux of shape [4, 3, 1] and
        maximum entry 4
        sage: GeneralizedShiftedPrimedTableaux([4,3,1], 
        ....:                                  max_entry=4).cardinality()
        2688
        sage: SPTC = GeneralizedShiftedPrimedTableaux([3,2], max_entry=3)
        sage: T = SPTC[-1]
        sage: T
        [(1', 2', 2), (3', 3)]
        sage: SPTC[0]
        [(1, 1, 1), (2, 2)]
        sage: SPTC.cardinality()
        88
    """
    @staticmethod
    def __classcall_private__(cls, shape, max_entry=None, skew=None):
        r"""
        Normalize the attributes for the class.

        TESTS:

            sage: SPT = GeneralizedShiftedPrimedTableaux(shape=[2,1])
            sage: SPT._shape.parent()
            Partitions
            sage: SPT1 = GeneralizedShiftedPrimedTableaux(shape=(2,1), max_entry=3)
            sage: SPT2 = GeneralizedShiftedPrimedTableaux(shape=[2,1], max_entry=3)
            sage: SPT1 is SPT2
            True
        """
        shape = _Partitions(shape)
        return super(GeneralizedShiftedPrimedTableaux_shape, cls).__classcall__(cls,
                     shape=shape, max_entry=max_entry, skew=skew)

    def __init__(self, shape, max_entry, skew):
        r"""
        Initialize the class of generalized shifted primed tableaux
        of a given shape.

        TESTS:

            sage: SPT = GeneralizedShiftedPrimedTableaux([4,2,1], max_entry=4)
            sage: TestSuite(SPT).run()  # long time
        """
        if max_entry is None:
            max_entry = sum(shape)

        GeneralizedShiftedPrimedTableaux.__init__(self, skew=skew)
        if skew is None:
            Parent.__init__(self, category=FiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=Sets().Finite())
        self._max_entry = max_entry
        self._skew = skew
        if skew is None:
            self._shape = Partition(shape)
        else:
            self._shape = SkewPartition((shape, skew))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS:

            sage: GeneralizedShiftedPrimedTableaux([3,2,1])
            Generalized Shifted Primed Tableaux of shape [3, 2, 1] and maximum
            entry 6
        """
        base = "Generalized Shifted Primed Tableaux of shape " + self._shape._repr_()
        if self._max_entry is not None:
            base += " and maximum entry {}".format(self._max_entry)
        return base

    def _contains_tableau(self, T):
        r"""
        Check if ``self`` contains preprocessed tableau ``T``.

        TESTS:

            sage: t = GeneralizedShiftedPrimedTableau._preprocess([[1,'2p',2,2],[2,'3p']])
            sage: GeneralizedShiftedPrimedTableaux([4,2],max_entry=4)._contains_tableau(t)
            True
            sage: s = GeneralizedShiftedPrimedTableau._preprocess([[1,'2p',2],[2,'3p']])
            sage: GeneralizedShiftedPrimedTableaux([4,2])._contains_tableau(s)
            False
        """
        if not super(GeneralizedShiftedPrimedTableaux_shape, self)._contains_tableau(T):
            return False

        shape = [len(row) for row in T]
        skew = [row.count(None) for row in T]
        if sum(skew) == 0:
            shape = Partition(shape)
        else:
            shape = SkewPartition((shape, skew))
        if self._shape != shape:
            return False

        if self._max_entry is not None:
            flat = [item.integer() for sublist in T for item in sublist]
            if flat == []:
                max_entry = 0
            else:
                max_entry = int(max(flat))
            if max_entry > self._max_entry:
                return False

        return True

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: Tabs = GeneralizedShiftedPrimedTableaux(shape=(3,2))
            sage: Tabs[:4]
            [[(1, 1, 1), (2, 2)],
             [(1, 1, 1), (2', 2)],
             [(1', 1, 1), (2, 2)],
             [(1', 1, 1), (2', 2)]]
            sage: len(list(Tabs))
            568
        """
        from sage.combinat.permutation import Permutations
        list_weights = []
        for partition in Partitions(sum(self._shape)):
            if len(partition) <= self._max_entry:
                for p in Permutations(partition):
                    for padding in range(self._max_entry-len(p)+1):
                        list_weights.append([0]*padding + list(p))
        for weight in list_weights:
            for T in GeneralizedShiftedPrimedTableaux(weight=tuple(weight), shape=self._shape):
                yield T


class GeneralizedShiftedPrimedTableaux_weight(GeneralizedShiftedPrimedTableaux):
    r"""
    Generalized Shifted primed tableaux of fixed weight.

    EXAMPLES::

        sage: GeneralizedShiftedPrimedTableaux(weight=(2,3,1))
        Generalized Shifted Primed Tableaux of weight (2, 3, 1)
        sage: GeneralizedShiftedPrimedTableaux(weight=(2,3,1)).cardinality()
        64
        sage: T = GeneralizedShiftedPrimedTableaux(weight=(3,2))
        sage: T[:5]
        [[(1, 1, 1, 2, 2)],
         [(1, 1, 1, 2', 2)],
         [(1', 1, 1, 2, 2)],
         [(1', 1, 1, 2', 2)],
         [(1, 1, 1, 2), (2,)]]
        sage: T.cardinality()
        16
    """
    def __init__(self, weight, skew=None):
        r"""
        Initialize the class of generalized shifted primed tableaux of a
        given weight.

        TESTS:

            sage: TestSuite( GeneralizedShiftedPrimedTableaux(
            ....:                                   weight=(3,2,1))).run()
        """
        GeneralizedShiftedPrimedTableaux.__init__(self, skew=skew)
        if skew is None:
            Parent.__init__(self, category=FiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=Sets().Finite())
        self._weight = weight
        self._skew = skew

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS:

            sage: GeneralizedShiftedPrimedTableaux(weight=(3,2,1))
            Generalized Shifted Primed Tableaux of weight (3, 2, 1)
        """
        if self._skew is None:
            return "Generalized Shifted Primed Tableaux of weight {}".format(self._weight)
        return "Generalized Shifted Primed Tableaux of weight {} skewed by {}".format(self._weight, self._skew)

    def _contains_tableau(self, T):
        r"""
        Check if ``self`` contains preprocessed tableau ``T``.

        TESTS:

            sage: t = GeneralizedShiftedPrimedTableau._preprocess([[1,1.5],[2]])
            sage: GeneralizedShiftedPrimedTableaux(weight=(1,2))._contains_tableau(t)
            True
            sage: s = GeneralizedShiftedPrimedTableau._preprocess([[1,1.5],[3]])
            sage: GeneralizedShiftedPrimedTableaux(weight=(1,2))._contains_tableau(s)
            False

            sage: u = GeneralizedShiftedPrimedTableau._preprocess([])
            sage: GeneralizedShiftedPrimedTableaux(weight=())._contains_tableau(u)
            True
            sage: GeneralizedShiftedPrimedTableaux(weight=(1,2))._contains_tableau(u)
            False
        """
        if not super(GeneralizedShiftedPrimedTableaux_weight, self)._contains_tableau(T):
            return False
        flat = [item.integer() for sublist in T for item in sublist]
        if not flat:
            return not self._weight
        max_ind = max(flat)
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        return self._weight == weight

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: Tabs = GeneralizedShiftedPrimedTableaux(weight=(2,3))
            sage: Tabs[:4]
            [[(1, 1, 2, 2, 2)],
             [(1, 1, 2', 2, 2)],
             [(1', 1, 2, 2, 2)],
             [(1', 1, 2', 2, 2)]]
            sage: len(list(Tabs))
            16
        """
        for shape_ in ZS1_iterator(sum(self._weight)):
            if all(shape_[i] > shape_[i+1] for i in range(len(shape_)-1)):
                for tab in GeneralizedShiftedPrimedTableaux(shape=shape_, weight=self._weight,
                                                 skew=self._skew):
                    yield self.element_class(self, tab, check=False,
                                             preprocessed=True)


class GeneralizedShiftedPrimedTableaux_weight_shape(GeneralizedShiftedPrimedTableaux):
    r"""
    Generalized Shifted primed tableaux of the fixed weight and shape.

    EXAMPLES::

        sage: GeneralizedShiftedPrimedTableaux([4,2,1], weight=(2,3,2))
        Generalized Shifted Primed Tableaux of weight (2, 3, 2) and
        shape [4, 2, 1]
        sage: T = GeneralizedShiftedPrimedTableaux([4,2,1], weight=(2,3,2))
        sage: T[:6]
        [[(1, 1, 2, 2), (2, 3'), (3,)],
         [(1, 1, 2, 2), (2, 3'), (3',)],
         [(1, 1, 2, 2), (2', 3'), (3,)],
         [(1, 1, 2, 2), (2', 3'), (3',)],
         [(1, 1, 2', 3), (2, 2), (3,)],
         [(1, 1, 2', 3), (2, 2), (3',)]]
        sage: T.cardinality()
        32
    """
    def __init__(self, weight, shape, skew=None):
        r"""
        Initialize the class of generalized shifted primed tableaux of the
        given weight and shape.

        TESTS:

            sage: TestSuite( GeneralizedShiftedPrimedTableaux([4,2,1],
            ....:                                     weight=(3,2,2))).run()
        """
        GeneralizedShiftedPrimedTableaux.__init__(self, skew=skew)
        if skew is None:
            Parent.__init__(self, category=FiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=Sets().Finite())
        self._weight = weight
        self._skew = skew
        if skew is None:
            self._shape = _Partitions(shape)
        else:
            self._shape = SkewPartition((shape, skew))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS:

            sage: GeneralizedShiftedPrimedTableaux([3,2,1], weight=(4,2))
            Generalized Shifted Primed Tableaux of weight (4, 2) and
            shape [3, 2, 1]
        """
        return ("Generalized Shifted Primed Tableaux of weight {} and shape {}"
                .format(self._weight, self._shape))

    def _contains_tableau(self, T):
        r"""
        Check if ``self`` contains preprocessed tableau ``T``.

        TESTS:

            sage: t = GeneralizedShiftedPrimedTableau._preprocess([[1,1.5],[2]])
            sage: GeneralizedShiftedPrimedTableaux([2,1],
            ....:                            weight=(1,2))._contains_tableau(t)
            True
            sage: GeneralizedShiftedPrimedTableaux([2,1],
            ....:                            weight=(2,1))._contains_tableau(t)
            False
            sage: s = GeneralizedShiftedPrimedTableau._preprocess([[1,1.5,2,3],
            ....:                                                 [3]])
            sage: GeneralizedShiftedPrimedTableaux([3,2],
            ....:                          weight=(1,2,2))._contains_tableau(s)
            False
            sage: u = GeneralizedShiftedPrimedTableau._preprocess([])
            sage: GeneralizedShiftedPrimedTableaux([3,2],
            ....:                          weight=(1,2,2))._contains_tableau(u)
            False
            sage: GeneralizedShiftedPrimedTableaux([],
            ....:                               weight=())._contains_tableau(u)
            True
        """
        if not super(GeneralizedShiftedPrimedTableaux_weight_shape, self)._contains_tableau(T):
            return False

        flat = [item.integer() for sublist in T for item in sublist]
        if not flat:
            # It is sufficient only to check this because the weight
            #   and shape must be compatible
            return not self._weight

        max_ind = max(flat)
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        if self._weight != weight:
            return False

        shape = [len(row) for row in T]
        skew = [row.count(None) for row in T]
        if sum(skew) == 0:
            shape = _Partitions(shape)
        else:
            shape = SkewPartition((shape, skew))
        return self._shape == shape

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: Tabs = GeneralizedShiftedPrimedTableaux([3,2],weight=(1,2,2))
            sage: Tabs[:4]
            [[(1, 2, 2), (3, 3)],
             [(1, 2, 2), (3', 3)],
             [(1, 2', 3'), (2, 3)],
             [(1, 2', 3'), (2, 3')]]
            sage: len(list(Tabs))
            16

        TESTS:

            sage: Tabs = GeneralizedShiftedPrimedTableaux([3,2], weight=(1,4))
            sage: list(Tabs)
            []
        """
        if self._skew is not None:
            raise NotImplementedError("skew tableau must be empty")

        if not self._shape.dominates(sorted(self._weight, reverse=True)):
            return
        full_shape = self._shape
        sub_tab = []
        tab_list_new = [[]]
        half = ~QQ(2)
        for i, w in enumerate(self._weight):
            tab_list_old = tab_list_new
            tab_list_new = []
            for sub_tab in tab_list_old:
                sub_shape = [len(row) for row in sub_tab]
                for strip in _add_strip(sub_shape, full_shape, w):
                    l = len(strip) // 2
                    # new row can be added
                    new_tab = []
                    new_tab1 = None
                    if len(sub_shape) < len(full_shape):
                        new_tab = [sub_tab[r] + [i+half]*strip[r] + [i+1]*strip[-r-1]
                                   for r in range(l-1)]
                        if strip[l] != 0:
                            new_tab1 = new_tab[:]
                            new_tab.append([i+1] * strip[l])
                            new_tab1.append([i+half] + [i+1] * (strip[l]-1))
                    else:
                        new_tab = [sub_tab[r] + [i+half]*strip[r] + [i+1]*strip[-r-1]
                                   for r in range(l)]
                    tab_list_new.append(new_tab)
                    if new_tab1:
                        tab_list_new.append(new_tab1)
        for tab in tab_list_new:
            yield self.element_class(self, tab)
