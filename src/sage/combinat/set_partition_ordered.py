r"""
Ordered Set Partitions

AUTHORS:

- Mike Hansen

- MuPAD-Combinat developers (for algorithms and design inspiration)

- Travis Scrimshaw (2013-02-28): Removed ``CombinatorialClass`` and added
  entry point through :class:`OrderedSetPartition`.

- Aaron Lauve (2015-06-10): Added :class:`OrderedSetPartitions_all` and
  refactored some `__iter__` calls.
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import itertools
from sage.arith.all import factorial
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.sets.set import Set, is_Set
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.all import prod
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.combinat import stirling_number2
from sage.combinat.composition import Composition, Compositions
from sage.combinat.words.word import Word
import sage.combinat.permutation as permutation
from functools import reduce


class OrderedSetPartition(ClonableArray):
    """
    An ordered partition of a set.

    An ordered set partition `p` of a set `s` is a list of pairwise
    disjoint nonempty subsets of `s` such that the union of these
    subsets is `s`. These subsets are called the parts of the partition.
    We represent an ordered set partition as a list of sets. By
    extension, an ordered set partition of a nonnegative integer `n` is
    the set partition of the integers from `1` to `n`. The number of
    ordered set partitions of `n` is called the `n`-th ordered Bell
    number.

    There is a natural integer composition associated with an ordered
    set partition, that is the sequence of sizes of all its parts in
    order.

    The number `T_n` of ordered set partitions of
    `\{ 1, 2, ..., n \}` is the so-called `n`-th *Fubini number* (also
    known as the `n`-th ordered Bell number; see
    :wikipedia:`Ordered Bell number`). Its exponential generating
    function is

    .. MATH::

        \sum_n {T_n \over n!} x^n = {1 \over 2-e^x}.

    (See sequence A000670 in OEIS.)

    EXAMPLES:

    There are 13 ordered set partitions of `\{1,2,3\}`::

        sage: OrderedSetPartitions(3).cardinality()
        13

    Here is the list of them::

        sage: OrderedSetPartitions(3).list()
        [[{1}, {2}, {3}],
         [{1}, {3}, {2}],
         [{2}, {1}, {3}],
         [{3}, {1}, {2}],
         [{2}, {3}, {1}],
         [{3}, {2}, {1}],
         [{1}, {2, 3}],
         [{2}, {1, 3}],
         [{3}, {1, 2}],
         [{1, 2}, {3}],
         [{1, 3}, {2}],
         [{2, 3}, {1}],
         [{1, 2, 3}]]

    There are 12 ordered set partitions of `\{1,2,3,4\}` whose underlying
    composition is `[1,2,1]`::

        sage: OrderedSetPartitions(4,[1,2,1]).list()
        [[{1}, {2, 3}, {4}],
         [{1}, {2, 4}, {3}],
         [{1}, {3, 4}, {2}],
         [{2}, {1, 3}, {4}],
         [{2}, {1, 4}, {3}],
         [{3}, {1, 2}, {4}],
         [{4}, {1, 2}, {3}],
         [{3}, {1, 4}, {2}],
         [{4}, {1, 3}, {2}],
         [{2}, {3, 4}, {1}],
         [{3}, {2, 4}, {1}],
         [{4}, {2, 3}, {1}]]

    Since :trac:`14140`, we can create an ordered set partition directly by
    :class:`OrderedSetPartition` which creates the parent object by taking the
    union of the partitions passed in. However it is recommended and
    (marginally) faster to create the parent first and then create the ordered
    set partition from that. ::

        sage: s = OrderedSetPartition([[1,3],[2,4]]); s
        [{1, 3}, {2, 4}]
        sage: s.parent()
        Ordered set partitions of {1, 2, 3, 4}

    REFERENCES:

    :wikipedia:`Ordered_partition_of_a_set`
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, parts):
        """
        Create a set partition from ``parts`` with the appropriate parent.

        EXAMPLES::

            sage: s = OrderedSetPartition([[1,3],[2,4]]); s
            [{1, 3}, {2, 4}]
            sage: s.parent()
            Ordered set partitions of {1, 2, 3, 4}
            sage: t = OrderedSetPartition([[2,4],[1,3]]); t
            [{2, 4}, {1, 3}]
            sage: s != t
            True
            sage: OrderedSetPartition([])
            []
        """
        OP = OrderedSetPartitions()
        return OP.element_class(OP, parts)

    def __init__(self, parent, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: s = OS([[1, 3], [2, 4]])
            sage: TestSuite(s).run()
        """
        ClonableArray.__init__(self, parent, map(Set, s))

    def check(self):
        """
        Check that we are a valid ordered set partition.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: s = OS([[1, 3], [2, 4]])
            sage: s.check()
        """
        assert self in self.parent()

    @combinatorial_map(name='to packed word') # should this line be here?
    def to_packed_word(self):
        # TODO: create PackedWords class.
        r"""
        Return the word (a packed word) via standard bijection between set compositions and packed words.
        """
        # TODO: add a check about base_set (this bijection assumes base_set has a natural ordering)
        w = {}
        for i in range(len(self)):
            for j in self[i]:
                w[j] = i+1
        keys = sorted(w.keys())
        return Word([w[k] for k in keys])

    @combinatorial_map(name='to composition') # should this line be here?
    def to_composition(self):
        r"""
        Return the integer composition whose parts are the sizes of the sets
        in ``self``.

        EXAMPLES::

            sage: S = OrderedSetPartitions(5)
            sage: x = S([[3,5,4], [1, 2]])
            sage: x.to_composition()
            [3, 2]
            sage: y = S([[3,1], [2], [5,4]])
            sage: y.to_composition()
            [2, 1, 2]
        """
        return Composition(map(len, self))

    def base_set(self):
        """
        Return the base set of ``self``, which is the union of all parts
        of ``self``.

        EXAMPLES::

            sage: OrderedSetPartition([[1], [2,3], [4]]).base_set()
            {1, 2, 3, 4}
            sage: OrderedSetPartition([[1,2,3,4]]).base_set()
            {1, 2, 3, 4}
            sage: OrderedSetPartition([]).base_set()
            {}
        """
        return reduce(lambda x,y: x.union(y), self, Set([]))

    def base_set_cardinality(self):
        """
        Return the cardinality of the base set of ``self``, which is the sum
        of the sizes of the parts of ``self``.

        This is also known as the *size* (sometimes the *weight*) of
        a set partition.

        EXAMPLES::

            sage: OrderedSetPartition([[1], [2,3], [4]]).base_set_cardinality()
            4
            sage: OrderedSetPartition([[1,2,3,4]]).base_set_cardinality()
            4
        """
        return sum(len(x) for x in self)

    size = base_set_cardinality

    def apply_permutation(self, p):
        r"""
        Apply ``p`` to the underlying set of ``self``.

        INPUT:

        - ``p`` -- A permutation

        EXAMPLES::

            sage: x = OrderedSetPartition([[3,5,4], [1,2]])
            sage: p = Permutation([2,1,4,5,3])
            sage: x.apply_permutation(p)
            [{3, 4, 5}, {1, 2}]
            sage: q = Permutation([3,5,1,2,4])
            sage: x.apply_permutation(q)
            [{1, 2, 4}, {3, 5}]
        """
        return self.__class__(self.parent(), [Set(map(p, bloc)) for bloc in self])

    def restriction(self, I):
        """
        Return the restriction of ``self`` to a subset ``I``
        (which is given as a set or list or any other iterable).

        EXAMPLES::

            sage: A = OrderedSetPartition([[1], [2,3]])
            sage: A.restriction([1,2])
            {{1}, {2}}
            sage: A.restriction([2,3])
            {{2, 3}}
            sage: A.restriction([])
            {}
            sage: A.restriction([4])
            {}
        """
        assert Set(I).issubset(self.base_set()), "%s is not a subset of %s" %(I,self.base_set())
        ret = []
        for bloc in self:
            newbloc = [i for i in bloc if i in I]
            if len(newbloc) != 0:
                ret.append(newbloc)
        return OrderedSetPartition(ret)

class OrderedSetPartitions(UniqueRepresentation, Parent):
    """
    Return the combinatorial class of ordered set partitions of ``s``.

    EXAMPLES::

        sage: OS = OrderedSetPartitions([1,2,3,4]); OS
        Ordered set partitions of {1, 2, 3, 4}
        sage: OS.cardinality()
        75
        sage: OS.first()
        [{1}, {2}, {3}, {4}]
        sage: OS.last()
        [{1, 2, 3, 4}]
        sage: OS.random_element()
        [{3}, {1}, {2}, {4}]

    ::

        sage: OS = OrderedSetPartitions([1,2,3,4], [2,2]); OS
        Ordered set partitions of {1, 2, 3, 4} into parts of size [2, 2]
        sage: OS.cardinality()
        6
        sage: OS.first()
        [{1, 2}, {3, 4}]
        sage: OS.last()
        [{3, 4}, {1, 2}]
        sage: OS.list()
        [[{1, 2}, {3, 4}],
         [{1, 3}, {2, 4}],
         [{1, 4}, {2, 3}],
         [{2, 3}, {1, 4}],
         [{2, 4}, {1, 3}],
         [{3, 4}, {1, 2}]]

    ::

        sage: OS = OrderedSetPartitions("cat"); OS
        Ordered set partitions of {'a', 'c', 't'}
        sage: OS.list()
        [[{'a'}, {'c'}, {'t'}],
         [{'a'}, {'t'}, {'c'}],
         [{'c'}, {'a'}, {'t'}],
         [{'t'}, {'a'}, {'c'}],
         [{'c'}, {'t'}, {'a'}],
         [{'t'}, {'c'}, {'a'}],
         [{'a'}, {'c', 't'}],
         [{'c'}, {'a', 't'}],
         [{'t'}, {'a', 'c'}],
         [{'a', 'c'}, {'t'}],
         [{'a', 't'}, {'c'}],
         [{'c', 't'}, {'a'}],
         [{'a', 'c', 't'}]]
    """
    @staticmethod
    def __classcall_private__(cls, s=None, comp=None):
        """
        Choose the correct parent based upon input.

        EXAMPLES::

            #TODO: add another example...
            #sage: OrderedSetPartitions()
            #Ordered set partitions of integers

            sage: OrderedSetPartitions(4)
            Ordered set partitions of {1, 2, 3, 4}
            sage: OrderedSetPartitions(4, [1, 2, 1])
            Ordered set partitions of {1, 2, 3, 4} into parts of size [1, 2, 1]
        """
        if s is None:
            return OrderedSetPartitions_all()
        if isinstance(s, (int, sage.rings.integer.Integer)):
            if s < 0:
                raise ValueError("s must be non-negative")
            s = frozenset(range(1, s+1))
        else:
            try:
                if s.cardinality() == infinity:
                    raise ValueError("The set must be finite")
            except AttributeError:
                pass
            s = frozenset(s)

        if comp is None:
            return OrderedSetPartitions_s(s)

        if isinstance(comp, (int, sage.rings.integer.Integer)):
            return OrderedSetPartitions_sn(s, comp)
        if comp not in Compositions(len(s)):
            raise ValueError("comp must be a composition of %s"%len(s))
        return OrderedSetPartitions_scomp(s, Composition(comp))

    def __contains__(self, x):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4])
            sage: all([sp in OS for sp in OS])
            True
        """
        # x must be a list/tuple
        if not isinstance(x, (OrderedSetPartition, list, tuple)):
            return False

        # Check that all parts are disjoint
        base_set = reduce( lambda x,y: x.union(y), map(Set, x), Set([]) )
        if len(base_set) != sum(map(len, x)):
            return False

        # Check to make sure each element of x is a set
        for s in x:
            if not (isinstance(s, (set, frozenset)) or is_Set(s)):
                return False

        return True

    def _element_constructor_(self, s):
        """
        Construct an element of ``self`` from ``s``.

        INPUT:

        - ``s`` -- A set of sets

        EXAMPLES::

            sage: S = OrderedSetPartitions(4)
            sage: elt = S([[1,3],[2,4]]); elt
            {{1, 3}, {2, 4}}
            sage: P = OrderedSetPartitions()
            sage: P(elt).parent() is P
            True
            sage: S = OrderedSetPartitions([])
            sage: S([])
            {}
        """
        if isinstance(s, OrderedSetPartition):
            if s.parent() is self:
                return s
            if isinstance(s.parent(), OrderedSetPartitions):
                return self.element_class(self, list(s))
            raise ValueError("cannot convert %s into an element of %s"%(s, self))
        return self.element_class(self, s)

    Element = OrderedSetPartition

    def _iterator_comp(self, comp):
        """
        Return an iterator for the ordered set partitions with block sizes
        corresponding to the composition ``comp``.

        INPUT:

        -  ``comp`` -- a :class:`Composition` object

        EXAMPLES::

            sage: OP = OrderedSetPartitions(3)
            sage: it = OP._iterator_comp(Composition([1,1,1]))
            sage: list(sorted(map(list, next(it))))
            [[1], [2], [3]]
            sage: OP12 = OrderedSetPartitions(3,Composition([1,2]))
            sage: len(list(OP._iterator_comp(Composition([1,2])))) == OP12.cardinality()
            True
        """
        #comp = Composition(comp)
        lset = [x for x in self._set]
        n = len(comp)
        dcomp = [-1] + comp.descents(final_descent=True)

        p = []
        for j in range(n):
            p += [j + 1] * comp[j]

        for x in permutation.Permutations(p):
            res = permutation.to_standard(x).inverse()
            res = [lset[x - 1] for x in res]
            yield self.element_class(self, [Set(res[dcomp[i]+1:dcomp[i+1]+1])
                                            for i in range(n)])

class OrderedSetPartitions_all(OrderedSetPartitions):
    r"""
    All ordered set partitions.
    """
    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: OP = OrderedSetPartitions()
            sage: TestSuite(OP).run()
        """
        OrderedSetPartitions.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: OrderedSetPartitions()
            Ordered set partitions
        """
        return "Ordered set partitions"

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: it = SetPartitions().__iter__()
            sage: [next(it) for x in range(10)]
            [{}, {{1}}, {{1, 2}}, {{1}, {2}}, {{1, 2, 3}}, {{1}, {2, 3}},
             {{1, 3}, {2}}, {{1, 2}, {3}}, {{1}, {2}, {3}}, {{1, 2, 3, 4}}]
        """
        n = 0
        while True:
            for x in OrderedSetPartitions_s(frozenset(range(1, n+1))):
                yield self.element_class(self, list(x))
            n += 1

class OrderedSetPartitions_s(OrderedSetPartitions):
    """
    Class of ordered partitions of a set `S`.
    """
    @staticmethod
    def __classcall_private__(cls, s):
        """
        Normalize ``s`` to ensure a unique representation.

        EXAMPLES::

            sage: OP1 = OrderedSetPartitions(set([2,1,4]))
            sage: OP2 = OrderedSetPartitions([4,1,2])
            sage: OP3 = OrderedSetPartitions((1,2,4))
            sage: OP1 is OP2, OP1 is OP3
            (True, True)
        """
        return super(OrderedSetPartitions_s, cls).__classcall__(cls, frozenset(s))

    def __init__(self, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: OP = OrderedSetPartitions(3)
            sage: TestSuite(OP).run()
            sage: OrderedSetPartitions(0).list()
            [{}]
            sage: OrderedSetPartitions([]).list()
            [{}]
        """
        self._set = s
        #self._size = len(s)
        OrderedSetPartitions.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: OrderedSetPartitions([1,2,3,4])
            Ordered set partitions of {1, 2, 3, 4}
        """
        return "Ordered set partitions of %s"%Set(self._set)

    def cardinality(self):
        """
        EXAMPLES::

            sage: OrderedSetPartitions(0).cardinality()
            1
            sage: OrderedSetPartitions(1).cardinality()
            1
            sage: OrderedSetPartitions(2).cardinality()
            3
            sage: OrderedSetPartitions(3).cardinality()
            13
            sage: OrderedSetPartitions([1,2,3]).cardinality()
            13
            sage: OrderedSetPartitions(4).cardinality()
            75
            sage: OrderedSetPartitions(5).cardinality()
            541
        """
        return sum([factorial(k)*stirling_number2(len(self._set),k) for k in range(len(self._set)+1)])

    def __iter__(self):
        """
        EXAMPLES::

            sage: [ p for p in OrderedSetPartitions([1,2,3]) ]
            [[{1}, {2}, {3}],
             [{1}, {3}, {2}],
             [{2}, {1}, {3}],
             [{3}, {1}, {2}],
             [{2}, {3}, {1}],
             [{3}, {2}, {1}],
             [{1}, {2, 3}],
             [{2}, {1, 3}],
             [{3}, {1, 2}],
             [{1, 2}, {3}],
             [{1, 3}, {2}],
             [{2, 3}, {1}],compsn
             [{1, 2, 3}]]
        """
        for cc in Compositions(len(self._set)):
            for sc in self._iterator_comp(cc):
                yield self.element_class(self,sc)

class OrderedSetPartitions_sn(OrderedSetPartitions):
    def __init__(self, s, n):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4], 2)
            sage: OS == loads(dumps(OS))
            True
        """
        OrderedSetPartitions.__init__(self, s)
        self.n = n

    def __repr__(self):
        """
        TESTS::

            sage: OrderedSetPartitions([1,2,3,4], 2)
            Ordered set partitions of {1, 2, 3, 4} into 2 parts
        """
        return "Ordered set partitions of %s into %s parts"%(Set(self._set),self.n)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        The number of ordered partitions of a set of size `n` into `k`
        parts is equal to `k! S(n,k)` where `S(n,k)` denotes the Stirling
        number of the second kind.

        EXAMPLES::

            sage: OrderedSetPartitions(4,2).cardinality()
            14
            sage: OrderedSetPartitions(4,1).cardinality()
            1
        """
        return factorial(self.n)*stirling_number2(len(self._set), self.n)

    def __iter__(self):
        """
        EXAMPLES::

            sage: [ p for p in OrderedSetPartitions([1,2,3,4], 2) ]
            [[{1, 2, 3}, {4}],
             [{1, 2, 4}, {3}],
             [{1, 3, 4}, {2}],
             [{2, 3, 4}, {1}],
             [{1, 2}, {3, 4}],
             [{1, 3}, {2, 4}],
             [{1, 4}, {2, 3}],
             [{2, 3}, {1, 4}],
             [{2, 4}, {1, 3}],
             [{3, 4}, {1, 2}],
             [{1}, {2, 3, 4}],
             [{2}, {1, 3, 4}],
             [{3}, {1, 2, 4}],
             [{4}, {1, 2, 3}]]
        """
        for x in Compositions(len(self._set),length=self.n):
            for z in OrderedSetPartitions_scomp(self._set,x):
                yield self.element_class(self, z)

    def __contains__(self, x):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4], 2)
            sage: all([sp in OS for sp in OS])
            True
            sage: OS.cardinality()
            14
            sage: len(filter(lambda x: x in OS, OrderedSetPartitions([1,2,3,4])))
            14
        """
        if not OrderedSetPartitions_s.__contains__(self, x):
            return False
        return len(x) == self.n

class OrderedSetPartitions_scomp(OrderedSetPartitions):
    def __init__(self, s, comp):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4], [2,1,1])
            sage: OS == loads(dumps(OS))
            True
        """
        self._set = s
        OrderedSetPartitions.__init__(self, s)
        self.c = Composition(comp)
        self.n = self.c.size()

    def __repr__(self):
        """
        TESTS::

            sage: OrderedSetPartitions([1,2,3,4], [2,1,1])
            Ordered set partitions of {1, 2, 3, 4} into parts of size [2, 1, 1]
        """
        return "Ordered set partitions of %s into parts of size %s"%(Set(self._set), self.c)

    def __contains__(self, x):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4], [2,1,1])
            sage: all(sp in OS for sp in OS)
            True
            sage: OS.cardinality()
            12
            sage: len(filter(lambda x: x in OS, OrderedSetPartitions([1,2,3,4])))
            12
        """
        return OrderedSetPartitions.__contains__(self, x) and map(len, x) == self.c

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        The number of ordered set partitions of a set of length `k` with
        composition shape `\mu` is equal to

        .. MATH::

            \frac{k!}{\prod_{\mu_i \neq 0} \mu_i!}.

        EXAMPLES::

            sage: OrderedSetPartitions(5,[2,3]).cardinality()
            10
            sage: OrderedSetPartitions(0, []).cardinality()
            1
            sage: OrderedSetPartitions(0, [0]).cardinality()
            1
            sage: OrderedSetPartitions(0, [0,0]).cardinality()
            1
            sage: OrderedSetPartitions(5, [2,0,3]).cardinality()
            10
        """
        return factorial(len(self._set))/prod([factorial(i) for i in self.c])

    def __iter__(self):
        """
        TESTS::

            sage: [ p for p in OrderedSetPartitions([1,2,3,4], [2,1,1]) ]
            [[{1, 2}, {3}, {4}],
             [{1, 2}, {4}, {3}],
             [{1, 3}, {2}, {4}],
             [{1, 4}, {2}, {3}],
             [{1, 3}, {4}, {2}],
             [{1, 4}, {3}, {2}],
             [{2, 3}, {1}, {4}],
             [{2, 4}, {1}, {3}],
             [{3, 4}, {1}, {2}],
             [{2, 3}, {4}, {1}],
             [{2, 4}, {3}, {1}],
             [{3, 4}, {2}, {1}]]

            sage: len(OrderedSetPartitions([1,2,3,4], [1,1,1,1]))
            24

            sage: [ x for x in OrderedSetPartitions([1,4,7], [3]) ]
            [[{1, 4, 7}]]

            sage: [ x for x in OrderedSetPartitions([1,4,7], [1,2]) ]
            [[{1}, {4, 7}], [{4}, {1, 7}], [{7}, {1, 4}]]

            sage: [ p for p in OrderedSetPartitions([], []) ]
            [[]]

            sage: [ p for p in OrderedSetPartitions([1], [1]) ]
            [[{1}]]

        Let us check that it works for large size (:trac:`16646`)::

            sage: OrderedSetPartitions(42).first()
            [{1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {12},
            {13}, {14}, {15}, {16}, {17}, {18}, {19}, {20}, {21}, {22}, {23},
            {24}, {25}, {26}, {27}, {28}, {29}, {30}, {31}, {32}, {33}, {34},
            {35}, {36}, {37}, {38}, {39}, {40}, {41}, {42}]
        """
        comp = self.c
        for sc in self._iterator_comp(comp):
            yield self.element_class(self, sc)

##########################################################
# Deprecations


class SplitNK(OrderedSetPartitions_scomp):
    def __setstate__(self, state):
        r"""
        For unpickling old ``SplitNK`` objects.

        TESTS::

            sage: loads("x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1+.\xc8\xc9,"
            ....:   "\x89\xcf\xcb\xe6\n\x061\xfc\xbcA\xccBF\xcd\xc6B\xa6\xda"
            ....:   "Bf\x8dP\xa6\xf8\xbcB\x16\x88\x96\xa2\xcc\xbc\xf4b\xbd\xcc"
            ....:   "\xbc\x92\xd4\xf4\xd4\"\xae\xdc\xc4\xec\xd4x\x18\xa7\x905"
            ....:   "\x94\xd1\xb45\xa8\x90\r\xa8>\xbb\x90=\x03\xc85\x02r9J\x93"
            ....:   "\xf4\x00\xb4\xc6%f")
            Ordered set partitions of {0, 1, 2, 3, 4} into parts of size [2, 3]
        """
        self.__class__ = OrderedSetPartitions_scomp
        n = state['_n']
        k = state['_k']
        OrderedSetPartitions_scomp.__init__(self, range(state['_n']), (k,n-k))

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override("sage.combinat.split_nk", "SplitNK_nk", SplitNK)