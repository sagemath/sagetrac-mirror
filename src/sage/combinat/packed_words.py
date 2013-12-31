# -*- coding: utf-8 -*-
"""
Packed Words

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
from sage.structure.parent import Parent
from sage.rings.integer import Integer
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.disjoint_union_enumerated_sets import \
    DisjointUnionEnumeratedSets
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute
from sage.structure.list_clone import ClonableIntArray
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.family import Family
from sage.sets.set import Set
from sage.combinat.set_partition_ordered import OrderedSetPartitions
from sage.combinat.tools import transitive_ideal
from sage.misc.misc import uniq
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets


def to_pack(li):
    """
    The analogue map of the *standardization* (..see
    :func:`sage.combinat.permutation.to_standard`) for *packed words*.

    TESTS::

        sage: from sage.combinat.packed_words import to_pack
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


def ordered_partition_sets_to_packed_word(li):
    """
    This map build the associated *ordered partition sets* of a *packed word*.

    TESTS::

        sage: from sage.combinat.packed_words import \
        ....:     ordered_partition_sets_to_packed_word
        sage: ordered_partition_sets_to_packed_word(
        ....:     OrderedSetPartitions(5)[34]
        ....: )
        [2, 3, 5, 1, 4]
        sage: ordered_partition_sets_to_packed_word([Set([2,3]), Set([1,4])])
        [2, 1, 1, 2]
    """
    dic = {}
    i = 1
    for set in li:
        for p in set:
            dic[p] = i
        i += 1
    return PackedWord(list(dic.values()))


def quasi_inv(x):
    """ Compte le nombre de quasi inversions
    """
    assert (x in PackedWords())
    invX = []
    for i in range(len(x)):
        ai = x[i]
        for j in range(i + 1, len(x)):
            aj = x[j]
            if ai == aj:
                invX.append((1. / 2, (i, j)))
            elif ai > aj:
                invX.append((1, (i, j)))
    return invX


def quasi_cmp(x, y):
    assert (x in PackedWords() and y in PackedWords())
    if x.size() == y.size():
        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                invXij = 0 if x[i] == x[j] else 1 if x[i] > x[j] else -1
                invYij = 0 if y[i] == y[j] else 1 if y[i] > y[j] else -1
                if invXij == invYij:
                    continue
                else:
                    return cmp(invXij, invYij)
        return 0
    elif x.size() > y.size():
        return 1
    else:
        return -1


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

            sage: from sage.combinat.packed_words import PackedWords_all
            sage: issubclass(PackedWords_all().element_class, PackedWord)
            True
            sage: w0 = PackedWord([4,2,3,1,2])
            sage: w0.parent()
            Packed words
            sage: type(w0)
            <class 'sage.combinat.packed_words.PackedWords_all_with_category.element_class'>
            sage: w1 = PackedWords()([4,2,3,1,2])
            sage: w1.parent() is w0.parent()
            True
            sage: type(w1) is type(w0)
            True
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        """
        The automatic parent of the element of this class

        When calling the constructor of an element of this class, one needs a
        parent. This class attribute specifies which parent is used.

        EXAMPLES::

            sage: PackedWord._auto_parent
            Packed words
            sage: PackedWord().parent()
            Packed words
         """
        return PackedWords_all()

    def __init__(self, parent, li=None, check=True):
        """
        TESTS::

            sage: PackedWord([]).parent()
            Packed words
        """
        if li is None:
            li = []
#        if li in OrderedSetPartitions(len(li)):
#            li = to_pack(li)
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

    def to_ordered_partition_sets(self):
        """
        This method build an *ordered partition sets* associated to *self*.

        TESTS::

            sage: from sage.combinat.packed_words import \
            ....:     ordered_partition_sets_to_packed_word
            sage: opspw = ordered_partition_sets_to_packed_word
            sage: pw = PackedWords(6).random_element()
            sage: opspw(pw.to_ordered_partition_sets()) == pw
            True
            sage: PackedWord([1,2,3,1,1,3]).to_ordered_partition_sets()
            [{1, 4, 5}, {2}, {3, 6}]
        """
        import collections
        d = collections.defaultdict(list)
        for i in range(len(self)):
            d[self[i]].append(i + 1)
        return [Set(d[k]) for k in sorted(d.keys())]

    def to_composition(self):
        """
        Compute a *composition* associated to the parikh vector of *self*.

        TESTS::

            sage: PackedWord([1,2,3,1,1,3]).to_composition()
            [3, 1, 2]
        """
        from sage.combinat.composition import Composition
        if len(self) == 0:
            return Composition([])
        from sage.combinat.words.word import Word
        return Composition(Word(self).parikh_vector(range(1, max(self) + 1)))

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

#    def __repr__( self ):
#        return str( list( self ) ).replace( '[', '(' ).replace( ']', ')' ).replace( ',', '' )

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

    def quasi_permutohedron_succ(self):
        r"""
        .. TODO:: should be in OrderedPartitionSets?
        implement
        TESTS::

            sage: PackedWord([4,4,2,5,3,3,1,3]).quasi_permutohedron_succ()
            [[3, 3, 2, 4, 2, 2, 1, 2], [4, 4, 2, 4, 3, 3, 1, 3], [5, 5, 2, 6, 4, 3, 1, 3], [5, 5, 2, 6, 4, 4, 1, 3], [5, 4, 2, 6, 3, 3, 1, 3]]
        """
        osp = self.to_ordered_partition_sets()
        succ = []
        #### we use the bijection with ordered sets partitions ####
        #####
        # first operation M_i
        # # "the operator m_i acts on the j-th and the j+1th
        # # parentheses of a 'quasi-permutation' as follows:
        # # if each element of the j-th parenthese is smaller
        # # than all the elements of the (j+1) th, then one
        # # can merge these two parentheses into one single
        # # parenthese which contains the union of the elements
        # # of these two parentheses."
        # ## [(1), (2,3,4), (5,7), (6)]
        # ##  -- > 1 < 2 == true for i = 0 :
        # ##      [(1,2,3,4), (5,7), (6)]
        # ##  -- > 4 < 5 == true for i = 1 :
        # ##      [(1), (2,3,4,5,7), (6)]
        # ##  -- > 7 < 6 == false for i = 2.
        for i_part in range(1, len(osp)):
            if max (osp[i_part - 1]) < min (osp[i_part]):
                succ.append(osp[:i_part - 1] + [osp[i_part - 1].union(osp[i_part])] + osp[i_part + 1:])
        #####
        # second operation S_i,j
        # # "the operator S_i,j acts on the j-th parenthese of a
        # # 'quasi-permutation' as follows : it splits this parentheses
        # # into two parentheses, the second one containing the j
        # # smallest elements of the initial parenthese and the first
        # # one containing the others
        # ## [(1), (2,3,4), (5,7), (6)]
        # ## -- > for i = 0 : too short
        # ## -- > for i = 1 : len (2,3,4) = 3 OK :
        # ##      [(1), (3,4), (2), ...]
        # ##      [(1), (4), (2,3), ...]
        # ## -- > for i = 2 : len (5,7) = 2 OK :
        # ##      [..., (7), (5), (6)]
        # ## -- > for i = 3 : too short
        for i_part in range(len(osp)):
            part = sorted(osp[i_part])
            for j in range(1, len(part)):
                succ.append(
                    osp[:i_part] +
                    [uniq(part[j:]), uniq(part[:j])] +
                    osp[i_part + 1:])
        return [ordered_partition_sets_to_packed_word(li) for li in succ]

    def quasi_permutohedron_pred(self):
        r"""
        .. TODO:: should be in OrderedPartitionSets?
        """
        osp = self.to_ordered_partition_sets()
        pred = []
        #### we use the bijection with ordered sets partitions ####
        #####
        # first operation 'M_i^{-1}'
        for i_part in range(len(osp)):
            part = sorted(osp[i_part])
            for j in range(1, len(part)):
                pred.append (
                    osp[:i_part] +
                    [uniq(part[:j])] + [uniq(part[j:])] +
                    osp[i_part + 1:])
        #####
        # second operation "S_i,j^{-1}"
        for i_part in range(1, len(osp)):
            if min(osp[i_part - 1]) > max(osp[i_part]):
                pred.append(
                    osp[:i_part - 1] +
                    [osp[i_part - 1] + osp[i_part]] +
                    osp[i_part + 1:])
        return [ordered_partition_sets_to_packed_word(li) for li in pred]

    def quasi_permutohedron_smaller(self):
        """
        .. TODO:: examples
        """
        return transitive_ideal(lambda x: x.quasi_permutohedron_pred(), self)

    def quasi_permutohedron_greater(self):
        """
        .. TODO:: examples
        """
        return transitive_ideal(lambda x: x.quasi_permutohedron_succ(), self)


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

            sage: from sage.combinat.packed_words import PackedWords_size, \
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


#==============================================================================
# Enumerated set of all packed words
#==============================================================================
class PackedWords_all(DisjointUnionEnumeratedSets, PackedWords):

    def __init__(self):
        """
        TESTS::

            sage: from sage.combinat.packed_words import PackedWords_all
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
        Parent.__init__(self, category=InfiniteEnumeratedSets())
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

    def __contains__(self, w):
        """
        TESTS::

            sage: P = PackedWords()
            sage: 1 in P
            False
            sage: PackedWord([]) in P
            True
            sage: [1,1,4,2,3] in P
            True

        """
        if isinstance(w, self.element_class):
            return True
        try:
            self(w)
            return True
        except:
            return False

    def __call__(self, x=None, *args, **keywords):
        """
        Ensure that ``None`` instead of ``0`` is passed by default.

        TESTS::

            sage: P = PackedWords()
            sage: P()
            []
        """
        return super(PackedWords, self).__call__(x, *args, **keywords)

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: P = PackedWords()
            sage: P._element_constructor_()
            []
            sage: P()
            []
            sage: P([1,1,1])
            [1, 1, 1]
        """
        return self.element_class(self, *args, **keywords)

    def combinatorial_class_of_size(self, size):
        assert(size >= 0), "size (%d) must be a positive integer" % size
        return PackedWords(size)

    def permutation_to_packed_words(self, sigma):
        """
        TESTS::

            sage: PW = PackedWords()
            sage: PW.permutation_to_packed_words(Permutation([3,1,2,4]))
            [[2, 1, 1, 2], [2, 1, 1, 3], [3, 1, 2, 3], [3, 1, 2, 4]]
            sage: PW.permutation_to_packed_words(Permutation([1,2,3]))
            [[1, 1, 1], [1, 1, 2], [1, 2, 2], [1, 2, 3]]
        """
        return PackedWords_size(len(sigma)).permutation_to_packed_words(sigma)

    Element = PackedWord


#==============================================================================
# Enumerated set of packed words of a given size
#==============================================================================
class PackedWords_size(PackedWords):
    """
    TESTS::

        sage: from sage.combinat.packed_words import PackedWords_size
        sage: for i in range(6): TestSuite(PackedWords_size(i)).run()

    """

    def __init__(self, size):
        super(PackedWords_size, self).__init__(category=FiniteEnumeratedSets())
        self._size = size

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
        """
        return isinstance(x, self.element_class) and x.size() == self._size

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

            sage: from sage.combinat.packed_words import PackedWords_size
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

            sage: from sage.combinat.packed_words import PackedWords_size
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
            yield self._element_constructor_()
        else:
            osp = OrderedSetPartitions(self._size)
            for part in osp:
                yield ordered_partition_sets_to_packed_word(part)

    @lazy_attribute
    def _parent_for(self):
        """
        The parent of the element generated by ``self`

        TESTS::

            sage: P = PackedWords(4)
            sage: P._parent_for
            Packed words
        """
        return PackedWords_all()

    @lazy_attribute
    def element_class(self):
        """
        TESTS::

            sage: P = PackedWords(4)
            sage: P.element_class
            <class 'sage.combinat.packed_words.PackedWords_all_with_category.element_class'>
            sage: P.first().__class__ == PackedWords().first().__class__
            True
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

    def permutation_to_packed_words(self, sigma):
        """
        Compute all packed words which give *sigma* by standardization.

        TESTS::

            sage: PW = PackedWords()
            sage: PW.permutation_to_packed_words(Permutation([3,1,2,4]))
            [[2, 1, 1, 2], [2, 1, 1, 3], [3, 1, 2, 3], [3, 1, 2, 4]]
            sage: PW.permutation_to_packed_words(Permutation([1,2,3]))
            [[1, 1, 1], [1, 1, 2], [1, 2, 2], [1, 2, 3]]
        """
        if self._size <= 1:
            if self._size == 0:
                return [self._element_constructor_([])]
            if self._size == 1:
                return [self._element_constructor_([1])]
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
                res.append(self._element_constructor_(pw.values()))
            pw[index_i] = pw[l_index] + 1
            res.append(self._element_constructor_(pw.values()))
        return res
