# -*- coding: utf-8 -*-
"""
Packed Words

References:
-----------

.. [NoTh06] Polynomial realizations of some trialgebras,
    J.-C. Novelli and J.-Y. Thibon.

AUTHOR:

- Jean-Baptiste Priez
"""
#*****************************************************************************
#       Copyright (C) 2013 Jean-Baptiste Priez <jbp@kerios.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import itertools
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.shuffle import ShuffleProduct
from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute
from sage.structure.parent import Parent
from sage.rings.integer import Integer
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.disjoint_union_enumerated_sets import \
    DisjointUnionEnumeratedSets
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.list_clone import ClonableIntArray
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.family import Family
from sage.sets.set import Set
from sage.combinat.set_partition_ordered import OrderedSetPartitions, OrderedSetPartition
from sage.misc.misc import uniq
from sage.combinat.composition import Composition


def to_pack(li):
    """
    The analogue map of the *standardization* (..see
    :func:`sage.combinat.permutation.to_standard`) for *packed words*.

    .. see _[NoTh06] §2. The Hopf algebra WQSym


    TESTS::

        sage: from sage.combinat.packed_word import to_pack
        sage: to_pack([])
        []
        sage: to_pack([3,1])
        [2, 1]
        sage: to_pack([1, 0, 0])
        [2, 1, 1]
        sage: to_pack([3,1,55])
        [2, 1, 3]
        sage: to_pack([11,4,1,55])
        [3, 2, 1, 4]
        sage: to_pack([11,4,1,11,4])
        [3, 2, 1, 3, 2]
    """
    l = uniq(li)
    return PackedWord([l.index(i) + 1 for i in li])


class PackedWord(ClonableIntArray):
    """
    The class of packed words

    TESTS::

        sage: PackedWord()
        []
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        """
        Ensure that packed words created by the enumerated sets and directly
        are the same and that they are instances of :class:`PackedWord`

        TESTS::

            sage: from sage.combinat.packed_word import PackedWords_all
            sage: issubclass(PackedWords_all().element_class, PackedWord)
            True
            sage: w0 = PackedWord([4,2,3,1,2])
            sage: w0.parent()
            Packed words
            sage: type(w0)
            <class 'sage.combinat.packed_word.PackedWords_all_with_category.element_class'>
            sage: w1 = PackedWords()([4,2,3,1,2])
            sage: w1.parent() is w0.parent()
            True
            sage: type(w1) is type(w0)
            True
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        return PackedWords_all()

    def __init__(self, parent, li=None, check=True):
        """

        TESTS::

            sage: PackedWord([]).parent()
            Packed words
        """
        if li is None:
            li = []
        ClonableIntArray.__init__(self, parent, li, check=check)

    def check(self):
        """
        Checks that ``self`` is a packed word

        TESTS::

            sage: PackedWord([3,3,2,1])
            [3, 3, 2, 1]

            #sage: PackedWord([2,2,1,0,4])
            #Traceback (most recent call last)
            #...
            #AssertionError: This is not a packed word
        """
        s = uniq(self)
        assert(len(s) == 0 or (max(s) == len(s) and min(s) == 1)
            ), "This is not a packed word %s" % str(self)

    def to_ordered_set_partition(self):
        """
        This method build an *ordered partition sets* associated to *self*.

        TESTS::

            sage: pw = PackedWords(6).random_element()
            sage: pw.to_ordered_set_partition().to_packed_word() == pw
            True
            sage: PackedWord([1,2,3,1,1,3]).to_ordered_set_partition()
            [{1, 4, 5}, {2}, {3, 6}]
        """
        import collections
        d = collections.defaultdict(list)
        for i in range(len(self)):
            d[self[i]].append(i + 1)
        return OrderedSetPartition([Set(d[k]) for k in sorted(d.keys())])

    @combinatorial_map(name='to composition')
    def to_composition(self):
        """
        Compute a *composition* associated to the parikh vector of *self*.

        TESTS::

            sage: PackedWord([1,2,3,1,1,3]).to_composition()
            [3, 1, 2]
            sage: PackedWord([1,2,3,1,1,3]).to_ordered_set_partition().to_composition()
            [3, 1, 2]
            sage: for pw in PackedWords(4):
            ....:     assert(pw.to_composition() == pw.to_ordered_set_partition().to_composition())
        """

        if len(self) == 0:
            return Composition([])
        li = list(self)
        return Composition([li.count(i) for i in set(self)])

    def shifted_shuffle(self, other):
        """
        The analogue map of the *shifted_shuffle* (..see
        :meth:`sage.combinat.permutation.Permutation_class.shifted_shuffle`)
        for *packed words*:

        MATH::

            p_1\dots p_k \Cup q_1\dots q_l[m] := p_1 \dot (p_2 \dots p_k \Cup (q_1 \dots q_l)[m])
                + (q_1 + m) \dot (p_1 \dots p_k \Cup (q_2 \dots q_l)[m])\,.

        with `m := \max(p_1 \dots p_k)`.

        .. see _[NoTh06] §2. The Hopf algebra WQSym

        TESTS::

            sage: PW = PackedWord
            sage: list(PW([1,1]).shifted_shuffle(PW([1,2])))
            [[1, 1, 2, 3],
             [1, 2, 1, 3],
             [1, 2, 3, 1],
             [2, 1, 1, 3],
             [2, 1, 3, 1],
             [2, 3, 1, 1]]
        """
        assert(other in self.parent())
        shift = max(self)
        return iter(ShuffleProduct(self, [i + shift for i in other]))

    def shifted_concatenation(self, other, side='right'):
        """
        The analogue map of the *shifted_concatenation* (..see
        :meth:`sage.combinat.permutation.Permutation_class.shifted_concatenation`)
        for *packed words*:

        MATH::

            p_1\dots p_k \bullet q_1\dots q_l[m] :=
                p_1 \dots p_k \dot (q_1 + m) \dots (q_l + m)\,.

        with `m := \max(p_1 \dots p_k)`.

        TESTS::

            sage: PackedWord([1,1,2]).shifted_concatenation([1,3,2,2])
            [1, 1, 2, 3, 5, 4, 4]
            sage: PackedWord([1,1,2]).shifted_concatenation([1,3,2,2], side='right')
            [1, 1, 2, 3, 5, 4, 4]
            sage: PackedWord([1,1,2]).shifted_concatenation([1,3,2,2], side='left')
            [3, 5, 4, 4, 1, 1, 2]
            sage: PackedWord([1,1,2]).shifted_concatenation([1,3,2,2], side='toto')
            Traceback (most recent call last):
            ...
            ValueError: toto must be "left" or "right"
        """
        assert(other in self.parent())
        PW = self.parent()._element_constructor
        shift = max(self)
        if side == "right" :
            return PW(list(self) + [a + shift for a in other])
        elif side == "left" :
            return PW([a + shift for a in other] + list(self))
        else :
            raise ValueError, "%s must be \"left\" or \"right\"" %(side)

    def is_empty(self):
        """
        Returns whether ``self`` is the empty word.

        EXAMPLES::

            sage: PackedWord().is_empty()
            True
            sage: PackedWord([]).is_empty()
            True
            sage: PackedWord([2,1,2]).is_empty()
            False
        """
        return not self

    def _latex_(self):
        """
        TESTS::

            sage: latex(PackedWord([1,2,3,1,1,3]))
            123113
            sage: latex(PackedWord(range(1,11)))
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        """
        if max(self) >= 10:
            return str(list(self))
        return str(list(self)).replace(
            '[', '').replace(
            ']', '').replace(
            ', ', '')

    def size(self):
        """
        EXAMPLES::

            sage: PackedWord().size()
            0
            sage: PackedWord([2,1,1]).size()
            3
        """
        return len(self)

    def pseudo_permutohedron_succ(self):
        r"""
        Iterate the successor of the packed word ``self``.

        ..see _[NoTh06] §2.6 The pseudo-permutohedron

        TESTS::

            sage: PW = PackedWord
            sage: list(PW([1,2,3]).pseudo_permutohedron_succ())
            [[1, 1, 2], [1, 2, 2]]
            sage: list(_[1].pseudo_permutohedron_succ())
            [[1, 1, 1], [1, 3, 2]]

        """
        return itertools.imap(
            lambda osp: osp.to_packed_word(),
            self.to_ordered_set_partition().pseudo_permutohedron_succ()
        )

    def pseudo_permutohedron_pred(self):
        r"""
        Iterate the predecessor of the packed word ``self``.

        ..see _[NoTh06] §2.6 The pseudo-permutohedron

        TESTS::

            sage: PW = PackedWord
            sage: list(PW([3,2,1]).pseudo_permutohedron_pred())
            [[2, 1, 1], [2, 2, 1]]
            sage: list(_[1].pseudo_permutohedron_pred())
            [[2, 3, 1], [1, 1, 1]]
        """
        return itertools.imap(
            lambda osp: osp.to_packed_word(),
            self.to_ordered_set_partition().pseudo_permutohedron_pred()
        )

    def pseudo_permutohedron_smaller(self):
        """
        Iterate through a list of packed words smaller than or equal to ``p``
        in the pseudo-permutohedron order.

        ..see _[NoTh06] §2.6 The pseudo-permutohedron

        TESTS::

            sage: set(PackedWord([3,2,1]).pseudo_permutohedron_smaller()) == \
                    set(PackedWords(3))
            True
        """
        return itertools.imap(
            lambda osp: osp.to_packed_word(),
            self.to_ordered_set_partition().pseudo_permutohedron_smaller()
        )

    def pseudo_permutohedron_greater(self):
        """
        Iterate through a list of packed words greater than or equal to ``p``
        in the pseudo-permutohedron order.

        ..see _[NoTh06] §2.6 The pseudo-permutohedron

        TESTS::

            sage: set(PackedWord([1,2,3]).pseudo_permutohedron_greater()) == \
                    set(PackedWords(3))
            True
        """
        return itertools.imap(
            lambda osp: osp.to_packed_word(),
            self.to_ordered_set_partition().pseudo_permutohedron_greater()
        )

    def half_inversions(self):
        """
        Return a list of the half inversions of ``self``.

        ..see: :meth:`sage.combinat.set_partition_ordered.OrderedSetPartition.half_inversions`.

        ..see _[NoTh06] §2.6 The pseudo-permutohedron

        TESTS::

            sage: PW = PackedWord
            sage: PW([1,1,2]).half_inversions()
            [(1, 2)]
            sage: PW([1,2,1]).half_inversions()
            [(1, 3)]
        """
        return self.to_ordered_set_partition().half_inversions()

    def inversions(self):
        """
        Return a list of the inversions of ``self``.

        An inversion of a packed word `p` is a pair `(i, j)` such that
        `i < j` and `p(i) > p(j)`.

        (The definition is same to inversions for permutations)

        ..see _[NoTh06] §2.6 The pseudo-permutohedron

        TESTS::

            sage: PW = PackedWord
            sage: PW([1,1,2]).inversions()
            []
            sage: PW([2,1,2]).inversions()
            [(1, 2)]
            sage: PW([2, 3, 2, 1, 1, 3, 3, 4]).inversions()
            [(1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5)]
        """
        n = len(self)
        return [(i + 1, j + 1) for i in range(n - 1)
                               for j in range(i + 1, n)
                    if self[i] > self[j]]

    def is_smaller_than(self, other):
        """
        ``self`` is smaller than ``other`` if the value of the
        inversion (i, j) in the table of inversions of ``self``
        is smaller than or equal to its value in the table of inversions
        of ``other``, for all (i, j).

        The table of inversion of ``p`` is given by the set of
        inversion ``p.inversions()`` with weight `1` and
        the set of half inversion with weight `1/2`.

        ..see _[NoTh06] §2.6 The pseudo-permutohedron

        TESTS::

            sage: PW = PackedWord
            sage: pw = PW([1,2,3])
            sage: forall(PackedWords(3), pw.is_smaller_than)
            (True, None)
            sage: PW([4,4,2,5,3,3,1,3]).is_smaller_than(PW([3,3,2,4,2,2,1,2]))
            True
        """
        assert(other in self.parent())
        invO = Set(other.inversions())
        for inv in self.inversions():
            if inv not in invO:
                return False
        halfO = Set(other.half_inversions())
        for half in self.half_inversions():
            if half not in invO and half not in halfO:
                return False
        return True


#==============================================================================
# Abstract class to serve as a Factory no instance are created
#==============================================================================
class PackedWords(UniqueRepresentation, Parent):
    """
    Factory class for packed words.

    INPUT:

    - ``size`` -- (optional) an integer

    OUTPUT:

    - the set of all packed words (of ``size`` (if specified))

    TESTS::

        sage: TestSuite(PackedWords()).run()

    EXAMPLES::

        sage: PackedWords()
        Packed words
        sage: PackedWords(4)
        Packed words of size 4
    """
    @staticmethod
    def __classcall_private__(cls, n=None):
        """
        TESTS::

            sage: from sage.combinat.packed_word import PackedWords_size, \
            ....:     PackedWords_all
            sage: isinstance(PackedWords(2), PackedWords)
            True
            sage: isinstance(PackedWords(), PackedWords)
            True
            sage: PackedWords(2) is PackedWords_size(2)
            True
            sage: PackedWords(5).cardinality()
            541
            sage: PackedWords() is PackedWords_all()
            True
        """
        if n is None:
            return PackedWords_all()
        else:
            assert(isinstance(n, (Integer, int)) and n >= 0), \
                "n must be a non negative integer"
            return PackedWords_size(Integer(n))

    def __init__(self, is_infinite=False):
        """
        Initialize ``self``

        TESTS::

            sage: TestSuite(PackedWords()).run()

        """
        if is_infinite:
            Parent.__init__(self, category=InfiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=FiniteEnumeratedSets())

    Element = PackedWord

    def __call__(self, x=None, *args, **keywords):
        """
        Ensure that ``None`` instead of ``0`` is passed by default.

        TESTS::

            sage: P = PackedWords()
            sage: P()
            []
        """
        return super(PackedWords, self).__call__(x, *args, **keywords)

    def __contains__(self, item):
        """
        TESTS::

            sage: [1,1,1] in PackedWords()
            True
            sage: [2,1,1] in PackedWords()
            True
            sage: [2,'2',1] in PackedWords()
            False
        """
        if isinstance(item, self.element_class):
            return True
        try:
            self._element_constructor(item)
            return True
        except TypeError:
            return False

    def permutation_to_packed_word(self, sigma):
        """
        TESTS::

            sage: PW = PackedWords()
            sage: PW.permutation_to_packed_word(Permutation([3,1,2,4]))
            [[2, 1, 1, 2], [2, 1, 1, 3], [3, 1, 2, 3], [3, 1, 2, 4]]
            sage: PW.permutation_to_packed_word(Permutation([1,2,3]))
            [[1, 1, 1], [1, 1, 2], [1, 2, 2], [1, 2, 3]]
        """
        return PackedWords_size(len(sigma)).permutation_to_packed_word(sigma)


#==============================================================================
# Enumerated set of all packed words
#==============================================================================
class PackedWords_all(DisjointUnionEnumeratedSets, PackedWords):

    def __init__(self):
        """
        TESTS::

            sage: from sage.combinat.packed_word import PackedWords_all
            sage: P = PackedWords_all()
            sage: P.cardinality()
            +Infinity
            sage: it = iter(P)
            sage: (it.next(), it.next(), it.next(), it.next(), it.next())
            ([], [1], [1, 2], [2, 1], [1, 1])
            sage: it.next().parent()
            Packed words
            sage: P([])
            []
            sage: P is PackedWords_all()
            True
            sage: TestSuite(P).run()
        """
        PackedWords.__init__(self, True)
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), PackedWords_size),
            facade=True, keepkey=False)

    def _repr_(self):
        """
        TESTS::

            sage: PackedWords()
            Packed words
        """
        return "Packed words"

    def combinatorial_class_of_size(self, size):
        return PackedWords(size)

    def __contains__(self, item):
        """
        TESTS::

            sage: [1,1,1] in PackedWords()
            True
            sage: [2,1,1] in PackedWords()
            True
            sage: [2,'2',1] in PackedWords()
            False
        """
        if isinstance(item, self.element_class):
            return True
        try:
            self._element_constructor(item)
            return True
        except TypeError:
            return False

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: P = PackedWords(0)
            sage: P([])
            []

            # sage: P([1])
            # Traceback (most recent call last)
            # ...
            # ValueError: Wrong size of word
        """
        return self.element_class(self, *args, **keywords)


#==============================================================================
# Enumerated set of packed words of a given size
#==============================================================================
class PackedWords_size(PackedWords):
    """
    TESTS::

        sage: from sage.combinat.packed_word import PackedWords_size
        sage: TestSuite(PackedWords_size(5)).run()

    """

    def __init__(self, size):
        super(PackedWords_size, self).__init__(False)
        assert(size >= 0), "Argument size (%d) must be a positive integer" % size
        self._size = size

    @lazy_attribute
    def _parent_for(self):
        return PackedWords_all()

    @lazy_attribute
    def element_class(self):
        """
        TESTS::

            sage: P4 = PackedWords(4)
            sage: type(P4.an_element())
            <class 'sage.combinat.packed_word.PackedWords_all_with_category.element_class'>
            sage: type(P4([1,2,3,4]))
            <class 'sage.combinat.packed_word.PackedWords_all_with_category.element_class'>
        """
        return self._parent_for.element_class

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: P = PackedWords(0)
            sage: P([])
            []

            # sage: P([1])
            # Traceback (most recent call last)
            # ...
            # ValueError: Wrong size of word
        """
        res = self.element_class(self._parent_for, *args, **keywords)
        if res.size() != self._size:
            raise ValueError, "Wrong size of word"
        return res

    def _repr_(self):
        """
        TESTS::

            sage: PackedWords(4)
            Packed words of size 4
        """
        return "Packed words of size %s" % (self._size)

    def __contains__(self, x):
        """
        TESTS::

            sage: P = PackedWords(4)
            sage: 1 in P
            False
            sage: PackedWord([]) in P
            False
            sage: PackedWord([1,2,1,3,1,4]) in P
            False
            sage: PackedWord([1,2,1,3]) in P
            True
            sage: [1,2,1,3] in P
            True

        """
        if isinstance(x, self.element_class):
            return self._size == x.size()
        else:
            return super(PackedWords_size, self).__contains__(x) and len(x) == self._size

    def _an_element_(self):
        """
        TESTS::

            sage: PackedWords(6).an_element()
            [1, 2, 3, 4, 5, 6]
        """
        return self.first()

    def cardinality(self):
        """
        Stirling number ???

        TESTS::

            sage: from sage.combinat.packed_word import PackedWords_size
            sage: PackedWords_size(0).cardinality()
            1
            sage: PackedWords_size(1).cardinality()
            1
            sage: PackedWords_size(2).cardinality()
            3
            sage: PackedWords_size(3).cardinality()
            13
        """
        return OrderedSetPartitions(self._size).cardinality()

    def __iter__(self):
        """
        TESTS::

            sage: from sage.combinat.packed_word import PackedWords_size
            sage: list(PackedWords_size(0))
            [[]]
            sage: list(PackedWords_size(1))
            [[1]]
            sage: list(PackedWords_size(2))
            [[1, 2], [2, 1], [1, 1]]
            sage: list(PackedWords_size(3))
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1],
             [1, 2, 2], [2, 1, 2], [2, 2, 1], [1, 1, 2], [1, 2, 1], [2, 1, 1],
             [1, 1, 1]]
        """
        if self._size == 0:
            yield self._element_constructor()
        else:
            for osp in OrderedSetPartitions(self._size):
                yield osp.to_packed_word()

    def permutation_to_packed_word(self, sigma):
        """
        Compute all packed words which give *sigma* by standardization.

        TESTS::

            sage: PW = PackedWords()
            sage: PW.permutation_to_packed_word(Permutation([3,1,2,4]))
            [[2, 1, 1, 2], [2, 1, 1, 3], [3, 1, 2, 3], [3, 1, 2, 4]]
            sage: PW.permutation_to_packed_word(Permutation([1,2,3]))
            [[1, 1, 1], [1, 1, 2], [1, 2, 2], [1, 2, 3]]
        """
        if self._size <= 1:
            if self._size == 0:
                return [self._element_constructor([])]
            if self._size == 1:
                return [self._element_constructor([1])]
        li = [({sigma.index(1):1}, sigma.index(1))]
        for i in range(2, self._size):
            index_i = sigma.index(i)
            tmp = []
            for (pw, l_index) in li:
                if l_index < index_i:
                    pw[index_i] = pw[l_index]
                    tmp.append((dict(pw), index_i))
                pw[index_i] = pw[l_index] + 1
                tmp.append((dict(pw), index_i))
            li = tmp
        index_i = sigma.index(self._size)
        res = []
        for (pw, l_index) in li:
            if l_index < index_i:
                pw[index_i] = pw[l_index]
                res.append(self._element_constructor(pw.values()))
            pw[index_i] = pw[l_index] + 1
            res.append(self._element_constructor(pw.values()))
        return res
