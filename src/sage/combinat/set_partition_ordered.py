r"""
Ordered Set Partitions

AUTHORS:

- Mike Hansen

- MuPAD-Combinat developers (for algorithms and design inspiration)

- Travis Scrimshaw (2013-02-28): Removed ``CombinatorialClass`` and added
  entry point through :class:`OrderedSetPartition`.
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
from six import add_metaclass

from sage.arith.all import factorial
from sage.sets.set import Set, Set_generic
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.all import prod
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.structure.richcmp import richcmp
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.combinat import stirling_number2
from sage.combinat.composition import Composition, Compositions
from sage.combinat.words.words import Words
from sage.combinat.words.finite_word import FiniteWord_class
import sage.combinat.permutation as permutation
from functools import reduce

@add_metaclass(InheritComparisonClasscallMetaclass)
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
    `\{ 1, 2, \ldots, n \}` is the so-called `n`-th *Fubini number*
    (also known as the `n`-th ordered Bell number; see
    :wikipedia:`Ordered Bell number`). Its exponential generating
    function is

    .. MATH::

        \sum_n {T_n \over n!} x^n = {1 \over 2-e^x}.

    (See sequence A000670 in OEIS.)

    INPUT:

    - ``parts`` -- an object or iterable that defines an ordered set partition
      (e.g., a list of pairwise disjoint sets) or a packed word (e.g., a list
      of letters on some alphabet). If there is ambiguity and if the input should
      be treated as a packed word, the keyword ``from_word`` should be used.

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

    We can construct the ordered set partition from a word,
    which we consider as packed::

        sage: OrderedSetPartition([2,4,1,2])
        [{3}, {1, 4}, {2}]
        sage: OrderedSetPartition(from_word=[2,4,1,2])
        [{3}, {1, 4}, {2}]
        sage: OrderedSetPartition(from_word='bdab')
        [{3}, {1, 4}, {2}]

    REFERENCES:

    :wikipedia:`Ordered_partition_of_a_set`
    """
    @staticmethod
    def __classcall_private__(cls, parts=None, from_word=None):
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
            sage: OrderedSetPartition()
            []
            sage: OrderedSetPartition([])
            []
            sage: OrderedSetPartition('')
            []
            sage: OrderedSetPartition('bdab') == OrderedSetPartition(from_word='bdab')
            True
            sage: OrderedSetPartition('bdab') == OrderedSetPartition(Word('bdab'))
            True
        """
        if parts is None and from_word is None:
            P = OrderedSetPartitions([])
            return P.element_class(P, [])
        if from_word:
            return OrderedSetPartitions().from_finite_word(Words()(from_word))
        # if `parts` looks like a sequence of "letters" then treat it like a word.
        if parts in Words() or (len(parts) > 0 and (parts[0] in ZZ or isinstance(parts[0], str))):
            return OrderedSetPartitions().from_finite_word(Words()(parts))
        else:
            P = OrderedSetPartitions( reduce(lambda x,y: x.union(y), map(Set, parts), Set([])) )
            return P.element_class(P, parts)

    def __init__(self, parent, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: s = OS([[1, 3], [2, 4]])
            sage: TestSuite(s).run()
        """
        self._base_set = reduce(lambda x,y: x.union(y), map(Set, s), Set([]))
        ClonableArray.__init__(self, parent, [Set(_) for _ in s])


    def _repr_(self):
        """
        Return a string representation of ``self``.

        .. TODO::

            Sort the repr output of Sage's :class:`Set` and remove
            this method.

        EXAMPLES::

            sage: OrderedSetPartition([[1,3],[2,4]])
            [{1, 3}, {2, 4}]
        """
        return '[' + ', '.join(('{' + repr(sorted(x))[1:-1] + '}' for x in self)) + ']'

    def check(self):
        """
        Check that we are a valid ordered set partition.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: s = OS([[1, 3], [2, 4]])
            sage: s.check()
        """
        assert self in self.parent()

    def _hash_(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: s = OS([[1, 3], [2, 4]])
            sage: OSP = OrderedSetPartitions()
            sage: hash(s) == hash(OSP(s))
            True
        """
        return hash(tuple(self))

    def size(self):
        r"""
        Return the cardinality of the base set of ``self``.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: s = OS([[1, 3], [2, 4]])
            sage: s.size()
            4
        """
        return sum(map(len, self))

    def length(self):
        r"""
        Return the number of parts of ``self``.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: s = OS([[1, 3], [2, 4]])
            sage: s.length()
            2
        """
        return len(self)

    @combinatorial_map(name='to composition')
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
        return Composition([len(_) for _ in self])

    @combinatorial_map(name='to packed word')
    def to_packed_word(self):
        r"""
        Return the packed word on alphabet `\{1,2,3,\ldots}`
        corresponding to ``self``.

        A *packed word* on alphabet `\{1,2,3,\ldots\}` is any word whose
        maximum letter is the same as its total number of distinct letters.
        Let `P` be an ordered set partition of a set `X`.
        The corresponding packed word `w_1 w_2 \cdots w_n` is constructed
        by having letter `w_i = j` if the `i`-th smallest entry in `X`
        occurs in the `j`-th block of `P`.

        .. SEEALSO::

            :meth:`Word.to_ordered_set_partition`

        .. WARNING::

            This assumes there is a total order on the underlying
            set (``self._base_set``).

        EXAMPLES::

            sage: S = OrderedSetPartitions()
            sage: x = S([[3,5], [2], [1,4,6]])
            sage: x.to_packed_word()
            word: 321313
            sage: x = S([['a', 'c', 'e'], ['b', 'd']])
            sage: x.to_packed_word()
            word: 12121
        """
        X = sorted(self._base_set)
        out = {}
        for i in range(len(self)):
            for letter in self[i]:
                out[letter] = i
        return Words()([out[letter]+1 for letter in X])

class OrderedSetPartitions(UniqueRepresentation, Parent):
    """
    Return the combinatorial class of ordered set partitions of ``s``.

    The optional argument ``c``, if specified, restricts the parts of
    the partition to have certain sizes (the entries of ``c``).

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
    def __classcall_private__(cls, s=None, c=None):
        """
        Choose the correct parent based upon input.

        EXAMPLES::

            sage: OrderedSetPartitions(4)
            Ordered set partitions of {1, 2, 3, 4}
            sage: OrderedSetPartitions(4, [1, 2, 1])
            Ordered set partitions of {1, 2, 3, 4} into parts of size [1, 2, 1]
        """
        if s is None:
            if c is not None:
                raise NotImplementedError("cannot specify 'c' without specifying 's'")
            return OrderedSetPartitions_all()
        if isinstance(s, (int, Integer)):
            if s < 0:
                raise ValueError("s must be non-negative")
            s = frozenset(range(1, s+1))
        else:
            s = frozenset(s)

        if c is None:
            return OrderedSetPartitions_s(s)

        if isinstance(c, (int, Integer)):
            return OrderedSetPartitions_sn(s, c)
        if c not in Compositions(len(s)):
            raise ValueError("c must be a composition of %s"%len(s))
        return OrderedSetPartitions_scomp(s, Composition(c))

    def __init__(self, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: TestSuite(OS).run()
        """
        self._set = s
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _element_constructor_(self, s):
        """
        Construct an element of ``self`` from ``s``.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: OS([[1,3],[2,4]])
            [{1, 3}, {2, 4}]
        """
        if isinstance(s, OrderedSetPartition):
            raise ValueError("cannot convert %s into an element of %s"%(s, self))
        return self.element_class(self, list(s))

    Element = OrderedSetPartition

    def __contains__(self, x):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4])
            sage: all(sp in OS for sp in OS)
            True
        """
        #x must be a list
        if not isinstance(x, (OrderedSetPartition, list, tuple)):
            return False

        #The total number of elements in the list
        #should be the same as the number is self._set
        if sum(map(len, x)) != len(self._set):
            return False

        #Check to make sure each element of the list
        #is a set
        u = Set([])
        for s in x:
            if not isinstance(s, (set, frozenset, Set_generic)):
                return False
            u = u.union(s)

        #Make sure that the union of all the
        #sets is the original set
        if u != Set(self._set):
            return False

        return True

    def from_finite_word(self, w):
        r"""
        Return the unique ordered set partition of `{1, 2, ..., n}` corresponding to a word
        of length `n`.'

        .. SEEALSO::

            :meth:`Word.to_ordered_set_partition`

        EXAMPLES::

            sage: A = OrderedSetPartitions().from_finite_word('abcabcabd'); A
            [{1, 4, 7}, {2, 5, 8}, {3, 6}, {9}]
            sage: B = OrderedSetPartitions().from_finite_word([1,2,3,1,2,3,1,2,4])
            sage: A == B
            True
        """
        # TODO: fix this if statement.
        #       In fact, what we need is for the underlying alphabet to be sortable.
        if isinstance(w, (list, tuple, str, FiniteWord_class)):
            return self.element_class(self, Words()(w).to_ordered_set_partition())
        else:
            raise ValueError("Something is wrong: `from_finite_word` expects an object of type list/tuple/str/Word representing a finite word, received {}.".format(str(w)))


class OrderedSetPartitions_s(OrderedSetPartitions):
    """
    Class of ordered partitions of a set `S`.
    """
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
             [{2, 3}, {1}],
             [{1, 2, 3}]]
        """
        for x in Compositions(len(self._set)):
            for z in OrderedSetPartitions(self._set, x):
                yield self.element_class(self, z)

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

    def __contains__(self, x):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4], 2)
            sage: all(sp in OS for sp in OS)
            True
            sage: OS.cardinality()
            14
            sage: len(filter(lambda x: x in OS, OrderedSetPartitions([1,2,3,4])))
            14
        """
        return OrderedSetPartitions.__contains__(self, x) and len(x) == self.n

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

class OrderedSetPartitions_scomp(OrderedSetPartitions):
    def __init__(self, s, comp):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4], [2,1,1])
            sage: OS == loads(dumps(OS))
            True
        """
        OrderedSetPartitions.__init__(self, s)
        self.c = Composition(comp)

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
        return OrderedSetPartitions.__contains__(self, x) and [len(_) for _ in x] == self.c

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
        lset = [x for x in self._set]
        l = len(self.c)
        dcomp = [-1] + comp.descents(final_descent=True)

        p = []
        for j in range(l):
            p += [j + 1] * comp[j]

        from sage.combinat.permutation import Permutations_mset
        for x in Permutations_mset(p):
            res = permutation.to_standard(x).inverse()
            res = [lset[x - 1] for x in res]
            yield self.element_class(self, [Set(res[dcomp[i]+1:dcomp[i+1]+1])
                                            for i in range(l)])

class OrderedSetPartitions_all(OrderedSetPartitions):
    r"""
    Ordered set partitions of `\{1, \ldots, n\}` for all
    `n \in \ZZ_{\geq 0}`.
    """
    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: OS = OrderedSetPartitions()
            sage: TestSuite(OS).run()  # long time
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets())

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: it = iter(OrderedSetPartitions())
            sage: [next(it) for _ in range(10)]
            [[], [{1}], [{1}, {2}], [{2}, {1}], [{1, 2}],
             [{1}, {2}, {3}], [{1}, {3}, {2}], [{2}, {1}, {3}],
             [{3}, {1}, {2}], [{2}, {3}, {1}]]
        """
        n = 0
        while True:
            for X in OrderedSetPartitions(n):
                yield self.element_class(self, list(X))
            n += 1

    def _element_constructor_(self, s):
        """
        Construct an element of ``self`` from ``s``.

        EXAMPLES::

            sage: OS = OrderedSetPartitions()
            sage: OS([[1,3],[2,4]])
            [{1, 3}, {2, 4}]
        """
        if isinstance(s, OrderedSetPartition):
            gset = s.parent()._set
            if gset == frozenset(range(1,len(gset)+1)):
                return self.element_class(self, list(s))
            raise ValueError("cannot convert %s into an element of %s"%(s, self))
        return self.element_class(self, list(s))

    def __contains__(self, x):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4])
            sage: AOS = OrderedSetPartitions()
            sage: all(sp in AOS for sp in OS)
            True
            sage: [[1,3],[4],[5,2]] in AOS
            True
            sage: [[1,4],[3]] in AOS
            False
            sage: [[1,3],[4,2],[2,5]] in AOS
            False
        """
        if isinstance(x, OrderedSetPartition):
            if x.parent() is self:
                return True
            gset = x.parent()._set
            return gset == frozenset(range(1, len(gset)+1))

        # x must be a list or a tuple
        if not isinstance(x, (list, tuple)):
            return False

        # Check to make sure each element of the list is a set
        if any(not isinstance(s, (set, frozenset, list, tuple, Set_generic))
               for s in x):
            return False
        if not all(isinstance(s, (set, frozenset, Set_generic)) or len(s) == len(set(s)) for s in x):
            return False
        X = set(reduce(lambda A,B: A.union(B), x, set()))
        return len(X) == sum(len(s) for s in x) and X == set(range(1,len(X)+1))

    def _coerce_map_from_(self, X):
        """
        Return ``True`` if there is a coercion map from ``X``.

        EXAMPLES::

            sage: OSP = OrderedSetPartitions()
            sage: OSP._coerce_map_from_(OrderedSetPartitions(3))
            True
            sage: OSP._coerce_map_from_(OrderedSetPartitions(['a','b']))
            False
        """
        if X is self:
            return True
        if isinstance(X, OrderedSetPartitions):
            return X._set == frozenset(range(1,len(X._set)+1))
        return super(OrderedSetPartitions_all, self)._coerce_map_from_(X)

    def _repr_(self):
        """
        TESTS::

            sage: OrderedSetPartitions()
            Ordered set partitions
        """
        return "Ordered set partitions"

    class Element(OrderedSetPartition):
        def _richcmp_(left, right, op):
            """
            TESTS::

                sage: OSP = OrderedSetPartitions()
                sage: el1 = OSP([[1,3], [4], [2]])
                sage: el2 = OSP([[3,1], [2], [4]])
                sage: el1 == el1, el2 == el2, el1 == el2    # indirect doctest
                (True, True, False)
                sage: el1 <= el2, el1 >= el2, el2 <= el1    # indirect doctest
                (False, True, True)
            """
            return richcmp([sorted(s) for s in left], [sorted(s) for s in right], op)

##########################################################
# Deprecations


class SplitNK(OrderedSetPartitions_scomp):
    def __setstate__(self, state):
        r"""
        For unpickling old ``SplitNK`` objects.

        TESTS::

            sage: loads(b"x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1+.\xc8\xc9,"
            ....:   b"\x89\xcf\xcb\xe6\n\x061\xfc\xbcA\xccBF\xcd\xc6B\xa6\xda"
            ....:   b"Bf\x8dP\xa6\xf8\xbcB\x16\x88\x96\xa2\xcc\xbc\xf4b\xbd\xcc"
            ....:   b"\xbc\x92\xd4\xf4\xd4\"\xae\xdc\xc4\xec\xd4x\x18\xa7\x905"
            ....:   b"\x94\xd1\xb45\xa8\x90\r\xa8>\xbb\x90=\x03\xc85\x02r9J\x93"
            ....:   b"\xf4\x00\xb4\xc6%f")
            Ordered set partitions of {0, 1, 2, 3, 4} into parts of size [2, 3]
        """
        self.__class__ = OrderedSetPartitions_scomp
        n = state['_n']
        k = state['_k']
        OrderedSetPartitions_scomp.__init__(self, range(state['_n']), (k,n-k))

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override("sage.combinat.split_nk", "SplitNK_nk", SplitNK)

