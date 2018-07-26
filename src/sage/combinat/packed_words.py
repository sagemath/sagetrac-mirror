# -*- coding: utf-8 -*-
r"""
Packed words are a way to represent ordered set partitions.
A word `w` with letters in `\{1,...,n\}` is a packed word if 
for each number `k > 1` appearing in `w`, the number `k - 1` appears in `w` too.

Thus, `w` can be obtained from an ordered set partition by setting 
`w_i=j` if `i` belongs to the `j`-th block.

Here are the Packed Words of size 0 to 3:
\epsilon
1
11, 12, 21
111, 112, 121, 211, 122, 212, 221, 123, 132, 213, 231, 312, 321


AUTHORS:

- Jean-Baptiste Priez
- Hugo Mlodecki
"""
#*****************************************************************************
#       Copyright (C) 2012 Jean-Baptiste Priez <jbp@kerios.fr>,
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
from sage.structure.parent import Parent
from sage.rings.integer import Integer
from six import add_metaclass
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute
from sage.structure.list_clone import ClonableIntArray
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.family import Family
from sage.combinat.set_partition_ordered import OrderedSetPartitions
from sage.combinat.set_partition_ordered import OrderedSetPartition
from sage.combinat.tools import transitive_ideal
from sage.misc.misc import uniq
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from collections import defaultdict
from sage.combinat.composition import Composition
from sage.combinat.words.word import Word
from sage.combinat.words.words import Words
from sage.combinat.combinatorial_map import combinatorial_map



@add_metaclass(InheritComparisonClasscallMetaclass)
class PackedWord(ClonableIntArray):
    r"""
    The class of packed words

    TESTS::

        sage: PackedWord()
        []
    """

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        r"""
        Ensure that packed words created by the enumerated sets and directly
        are the same and that they are instances of :class:`PackedWord`

        TESTS::

            sage: from sage.combinat.packed_words import PackedWords_all
            sage: issubclass(PackedWords_all().element_class, PackedWord)
            True
            sage: w0 = PackedWord([4, 2, 3, 1, 2])
            sage: w0.parent()
            Packed words
            sage: type(w0)
            <class 'sage.combinat.packed_words.PackedWords_all_with_category.element_class'>
            sage: w1 = PackedWords()([4, 2, 3, 1, 2])
            sage: w1.parent() is w0.parent()
            True
            sage: type(w1) is type(w0)
            True
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        r"""
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
        r"""
        TESTS::

            sage: PackedWord([]).parent()
            Packed words
        """
        if li is None:
            li=[]
        ClonableIntArray.__init__(self, parent, li, check=check)

    def check(self):
        r"""
        Checks that ``self`` is a packed word

        TESTS::

            sage: PackedWord([3, 3, 2, 1])
            [3, 3, 2, 1]

            sage: PackedWord([2, 2, 1, 0, 4])
            Traceback (most recent call last):
            ...
            AssertionError: This is not a packed word [2, 2, 1, 0, 4]
        """
        s = uniq(self)
        assert(len(s) == 0 or (max(s) == len(s) and min(s) == 1)
            ), "This is not a packed word %s" % self

    @combinatorial_map(name='to ordered set partition')
    def to_ordered_set_partition(self):
        r"""
        This method build an *ordered set partition* associated to *self*.

        TESTS::

            sage: pw = PackedWords(6).random_element()
            sage: pw in PackedWords(6)
            True
            sage: PackedWord([1, 2, 3, 1, 1, 3]).to_ordered_set_partition()
            [{1, 4, 5}, {2}, {3, 6}]
        """
        d = defaultdict(list)
        for i in range(len(self)):
            d[self[i]].append(i + 1)
        return [set(d[k]) for k in sorted(d)]

    def to_composition(self):
        r"""
        See http://trac.sagemath.org/17058 for details.

        Compute a *composition* associated to the parikh vector of *self*.

        TESTS::
            sage: PackedWord([]).to_composition()
            []
            sage: PackedWord([1, 1, 1]).to_composition()
            [3]
            sage: PackedWord([1, 2, 1]).to_composition()
            [2, 1]
            sage: PackedWord([1, 2, 3, 1, 1, 3]).to_composition()
            [3, 1, 2]
        """
        
        if len(self) == 0:
            return Composition([])
        W = Words(range(1, max(self) + 1))
        return Composition([Word(self).evaluation_dict()[i + 1] \
                            for i in range(max(self))])

    def is_empty(self):
        r"""
        Returns whether ``self`` is the empty word.

        EXAMPLES::

            sage: PackedWord().is_empty()
            True
            sage: PackedWord([]).is_empty()
            True
            sage: PackedWord([2, 1, 2]).is_empty()
            False
        """
        return not self


    def _latex_(self):
        r"""
        TESTS::

            sage: latex(PackedWord([1, 2, 3, 1, 1, 3]))
            123113
            sage: latex(PackedWord(range(1, 11)))
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        """
        if max(self) >= 10:
            return str(list(self))
        return str(list(self)).replace(
            '[', '').replace(
            ']', '').replace(
            ', ', '')

    def size(self):
        r"""
        EXAMPLES::

            sage: PackedWord().size()
            0
            sage: PackedWord([2, 1, 1]).size()
            3
        """
        return len(self)

    
# #################     Right Weak Order     ##################################

    def inversions_right(pw):
        r""" 
        Return the set of right weak order inversions with the definition :
        inversions_right(pw) := {(i, j) in [1..n]² : i < j and pw[i] > pw[j]}

        EXAMPLES::

            sage: PackedWord([]).inversions_right()
            set()
            sage: PackedWord([1, 2, 3]).inversions_right()
            set()
            sage: PackedWord([1, 2, 1]).inversions_right()
            {(2, 3)}
            sage: PackedWord([3, 2, 1]).inversions_right()
            {(1, 2), (1, 3), (2, 3)}
            sage: PackedWord([2, 3, 4, 1, 2, 4, 3]).inversions_right()
            {(1, 4), (2, 4), (2, 5), (3, 4), (3, 5), (3, 7), (6, 7)}
        """
        PackedWord.check(pw)
        n = len(pw)
        return set((i + 1, j + 1)
               for i in range(n - 1)
               for j in range(i + 1, n)
               if pw[i] > pw[j])

    def coinversions_right(pw):
        r""" 
        Return the set of right weak order coinversions with the definition :
        coinversions_right(pw) := {(pw[i], pw[j]) in [1..m]² : i < j and pw[i] > pw[j]}

        EXAMPLES::

            sage: PackedWord([]).coinversions_right()
            set()
            sage: PackedWord([1, 2, 3]).coinversions_right()
            set()
            sage: PackedWord([1, 2, 1]).coinversions_right()
            {(2, 1)}
            sage: PackedWord([3, 2, 1]).coinversions_right()
            {(2, 1), (3, 1), (3, 2)}
            sage: PackedWord([2, 3, 4, 1, 2, 4, 3]).coinversions_right()
            {(2, 1), (3, 1), (3, 2), (4, 1), (4, 2), (4, 3)}
        """
        PackedWord.check(pw)
        n = len(pw)
        return set((pw[i], pw[j])
               for i in range(n - 1)
               for j in range(i + 1, n)
               if pw[i] > pw[j])

    def right_weak_order_succ(self):
        r"""
        Return the list of successor for right weak order with the definition :

        v is a right successor of u if there exist i < n - 1 such that 
        v is equal to u where the u[i] and u[i + 1] are inversed 
        and len(inversions_right(v)) = len(inversions_right(u)) + 1

        EXAMPLES::

            sage: PackedWord([]).right_weak_order_succ()
            []
            sage: PackedWord([1, 1, 1]).right_weak_order_succ()
            []
            sage: PackedWord([1, 2, 1]).right_weak_order_succ()
            [[2, 1, 1]]
            sage: PackedWord([3, 1, 2]).right_weak_order_succ()
            [[3, 2, 1]]
            sage: PackedWord([3, 1, 2, 1, 1, 2, 4]).right_weak_order_succ()
            [[3, 2, 1, 1, 1, 2, 4], [3, 1, 2, 1, 2, 1, 4], [3, 1, 2, 1, 1, 4, 2]]
        """
        succ = []
        n = len(self)
        for i in range(n - 1):
            if self[i]<self[i + 1]:
                succ.append(self[:i] + [self[i + 1], self[i]] + self[i + 2:])
        return [PackedWord(p) for p in succ]

    def right_weak_order_pred(self):
        r"""
        Return the list of predecessor for right weak order with the definition :

        u is a predecessor of v if there exist i < n - 1 such that 
        v is equal to u where the u[i] and u[i + 1] are inversed 
        and len(inversions_right(v)) = len(inversions_right(u)) + 1

        EXAMPLES::

            sage: PackedWord([]).right_weak_order_pred()
            []
            sage: PackedWord([1, 1, 1]).right_weak_order_pred()
            []
            sage: PackedWord([1, 2, 1]).right_weak_order_pred()
            [[1, 1, 2]]
            sage: PackedWord([3, 1, 2]).right_weak_order_pred()
            [[1, 3, 2]]
            sage: PackedWord([3, 1, 2, 1, 1, 2, 4]).right_weak_order_pred()
            [[1, 3, 2, 1, 1, 2, 4], [3, 1, 1, 2, 1, 2, 4]]
        """
        pred = []
        n = len(self)
        for i in range(n - 1):
            if self[i] > self[i + 1]:
                pred.append(self[:i] + [self[i + 1], self[i]] + self[i + 2:])
        return [PackedWord(p) for p in pred]

    def right_weak_order_smaller(self):
        r"""
        Return the list of smaller or equal packed words for the right weak order.
        (..see :func:`sage.combinat.packed_words.right_weak_order_pred`) 
        for more informations.

        EXAMPLES::

            sage: PackedWord([]).right_weak_order_smaller()
            [[]]
            sage: PackedWord([1, 1, 1]).right_weak_order_smaller()
            [[1, 1, 1]]
            sage: PackedWord([1, 2, 1]).right_weak_order_smaller()
            [[1, 1, 2], [1, 2, 1]]
            sage: PackedWord([3, 1, 2]).right_weak_order_smaller()
            [[1, 2, 3], [1, 3, 2], [3, 1, 2]]
            sage: PackedWord([3, 1, 2, 2, 4]).right_weak_order_smaller()
            [[1, 2, 2, 3, 4], [1, 2, 3, 2, 4], [1, 3, 2, 2, 4], [3, 1, 2, 2, 4]]
        """
        return transitive_ideal(lambda x: x.right_weak_order_pred(), self)

    def right_weak_order_greater(self):
        r"""
        Return the list of greater or equal packed words for the right weak order.
        (..see :func:`sage.combinat.packed_words.right_weak_order_succ`) 
        for more informations.

        EXAMPLES::
        
            sage: PackedWord([]).right_weak_order_greater()
            [[]]
            sage: PackedWord([1, 1, 1]).right_weak_order_greater()
            [[1, 1, 1]]
            sage: PackedWord([1, 2, 1]).right_weak_order_greater()
            [[1, 2, 1], [2, 1, 1]]
            sage: PackedWord([3, 2, 1]).right_weak_order_greater()
            [[3, 2, 1]]
            sage: PackedWord([3, 1, 2, 2, 4]).right_weak_order_greater()
            [[3, 1, 2, 2, 4],
             [3, 1, 2, 4, 2],
             [3, 1, 4, 2, 2],
             [3, 2, 1, 2, 4],
             [3, 2, 1, 4, 2],
             [3, 2, 2, 1, 4],
             [3, 2, 2, 4, 1],
             [3, 2, 4, 1, 2],
             [3, 2, 4, 2, 1],
             [3, 4, 1, 2, 2],
             [3, 4, 2, 1, 2],
             [3, 4, 2, 2, 1],
             [4, 3, 1, 2, 2],
             [4, 3, 2, 1, 2],
             [4, 3, 2, 2, 1]]
        """
        return transitive_ideal(lambda x: x.right_weak_order_succ(), self)

# #################     Left Weak Order     ###################################

    def inversions_left(pw):
        r"""
        Return the set of left weak order inversions with the definition :
        inversions_left(pw) := {(pw[i], pw[j]) in [1..m]² : pw[i] < pw[j] 
            and the first occurence of 'pw[i]' in pw is after
                the last occurence of 'pw[j]' in pw}

        EXAMPLES::

            sage: PackedWord([]).inversions_left()
            set()
            sage: PackedWord([1, 2, 3]).inversions_left()
            set()
            sage: PackedWord([1, 2, 1]).inversions_left()
            set()
            sage: PackedWord([3, 1, 2]).inversions_left()
            {(1, 3), (2, 3)}
            sage: PackedWord([3, 1, 4, 1, 2]).inversions_left()
            {(1, 3), (2, 3), (2, 4)}

        """
        if len(pw) == 0:
            return set()
        PackedWord.check(pw)
        n = len(pw)
        m = max(pw)
        return set((i, j)
               for i in range(1, m)
               for j in range(i + 1, m + 1)
               if pw.index(i)>n - pw[::-1].index(j) - 1)

    def coinversions_left(pw):
        r"""
        Return the set of left weak order inversions with the definition :
        inversions_left(pw) := {(i, j) in [1..n]² : pw[i] < pw[j] 
            and the first occurence of 'pw[i]' in pw is after
                the last occurence of 'pw[j]' in pw}

        EXAMPLES::

            sage: PackedWord([]).inversions_right()
            set()
            sage: PackedWord([1, 2, 3]).coinversions_right()
            set()
            sage: PackedWord([1, 2, 1]).coinversions_right()
            {(2, 1)}
            sage: PackedWord([3, 1, 2]).coinversions_right()
            {(3, 1), (3, 2)}
            sage: PackedWord([3, 1, 4, 1, 2]).coinversions_right()
            {(3, 1), (3, 2), (4, 1), (4, 2)}

        """
        if len(pw) == 0:
            return set()
        PackedWord.check(pw)
        m = max(pw)
        return set((pw.index(i), n - pw[::-1].index(j) - 1)
               for i in range(1, n)
               for j in range(i + 1,n + 1)
                   if p.index(i) > n - pw[::-1].index(j) - 1)

    def left_weak_order_succ(self):
        r"""
        Return the list of successor for left weak order with the definition :

        v is a left successor of u if there exist i < n - 1 such that 
        v is equal to u where the i and i + 1 are inversed 
        and len(inversions_left(v)) = len(inversions_left(u)) + 1

        EXAMPLES::

            sage: PackedWord([]).left_weak_order_succ()
            []
            sage: PackedWord([1, 1, 1]).left_weak_order_succ()
            []
            sage: PackedWord([1, 2, 1]).left_weak_order_succ()
            []
            sage: PackedWord([3, 1, 2]).left_weak_order_succ()
            [[3, 2, 1]]
            sage: PackedWord([1, 2, 4, 3, 3, 2]).left_weak_order_succ()
            [[2, 1, 4, 3, 3, 1]]
            sage: PackedWord([1, 2, 4, 3, 3]).left_weak_order_succ()
            [[2, 1, 4, 3, 3], [1, 3, 4, 2, 2]]
        """
        succ = []
        if len(self) > 0:
            m = max(self)
            for i in range(1, m):
                if len(self) - 1 - self[::-1].index(i) < self.index(i + 1):
                    l = []
                    for x in self:
                        if x == i:
                            l.append(i + 1)
                        elif x == i + 1:
                            l.append(i)
                        else: l.append(x)
                    succ.append(l)
        return [PackedWord(p) for p in succ]

    def left_weak_order_pred(self):
        r"""
        Return the list of successor for left weak order with the definition :

        u is a left predecessor of v if there exist i < n - 1 such that 
        v is equal to u where the i and i + 1 are inversed 
        and len(inversions_left(v)) = len(inversions_left(u)) + 1

        EXAMPLES::

            sage: PackedWord([]).left_weak_order_pred()
            []
            sage: PackedWord([1, 1, 1]).left_weak_order_pred()
            []
            sage: PackedWord([1, 2, 1]).left_weak_order_pred()
            []
            sage: PackedWord([3, 1, 2]).left_weak_order_pred()
            [[2, 1, 3]]
            sage: PackedWord([3, 1, 2, 4, 4]).left_weak_order_pred()
            [[2, 1, 3, 4, 4]]
        """
        pred = []
        if len(self) > 0:
            m = max(self)
            for i in range(1, m):
                if self.index(i) > len(self) - 1 - self[::-1].index(i + 1):
                    l = []
                    for x in self:
                        if x == i:
                            l.append(i + 1)
                        elif x == i + 1:
                            l.append(i)
                        else: l.append(x)
                    pred.append(l)
        return [PackedWord(p) for p in pred]

    def left_weak_order_smaller(self):
        r"""
        Return the list of smaller or equal packed words for the left weak order.
        (..see :func:`sage.combinat.packed_words.left_weak_order_pred`) 
        for more informations.

        EXAMPLES::

            sage: PackedWord([]).left_weak_order_smaller()
            [[]]
            sage: PackedWord([1, 1, 1]).left_weak_order_smaller()
            [[1, 1, 1]]
            sage: PackedWord([1, 2, 1]).left_weak_order_smaller()
            [[1, 2, 1]]
            sage: PackedWord([3, 1, 2]).left_weak_order_smaller()
            [[1, 2, 3], [2, 1, 3], [3, 1, 2]]
            sage: PackedWord([3, 1, 2, 4, 1, 2]).left_weak_order_smaller()
            [[1, 2, 3, 4, 2, 3], [2, 1, 3, 4, 1, 3], [3, 1, 2, 4, 1, 2]]
        """
        return transitive_ideal(lambda x: x.left_weak_order_pred(), self)

    def left_weak_order_greater(self):
        r"""
        Return the list of greater or equal packed words for the left weak order.
        (..see :func:`sage.combinat.packed_words.left_weak_order_succ`) 
        for more informations.

        EXAMPLES::
        
            sage: PackedWord([]).left_weak_order_greater()
            [[]]
            sage: PackedWord([1, 1, 1]).left_weak_order_greater()
            [[1, 1, 1]]
            sage: PackedWord([1, 2, 1]).left_weak_order_greater()
            [[1, 2, 1]]
            sage: PackedWord([3, 1, 2]).left_weak_order_greater()
            [[3, 1, 2], [3, 2, 1]]
            sage: PackedWord([3, 1, 2, 4, 1, 2]).left_weak_order_greater()
            [[3, 1, 2, 4, 1, 2], [4, 1, 2, 3, 1, 2]]
        """
        return transitive_ideal(lambda x: x.left_weak_order_succ(), self)

#==============================================================================
# Abstract class to serve as a Factory no instance are created
#==============================================================================
class PackedWords(UniqueRepresentation, Parent):
    r"""
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
        r"""
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

        
    @staticmethod
    def to_pack(li):
        r"""
        The analogue map of the *standardization* (..see
        :func:`sage.combinat.permutation.to_standard`) for *packed words*.

        EXAMPLES::

            sage: PackedWords.to_pack([])
            []
            sage: PackedWords.to_pack([3, 1])
            [2, 1]
            sage: PackedWords.to_pack([1, 0, 0])
            [2, 1, 1]
            sage: PackedWords.to_pack([3, 1, 55])
            [2, 1, 3]
            sage: PackedWords.to_pack([11, 4, 1, 55])
            [3, 2, 1, 4]
            sage: PackedWords.to_pack([11, 4, 1, 11, 4])
            [3, 2, 1, 3, 2]
        """
        
        l = uniq(list(li))
        return PackedWord([l.index(i) + 1 for i in li])


#==============================================================================
# Enumerated set of all packed words
#==============================================================================
class PackedWords_all(DisjointUnionEnumeratedSets, PackedWords):

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.combinat.packed_words import PackedWords_all
            sage: P = PackedWords_all()
            sage: P.cardinality()
            +Infinity
            sage: it = iter(P)
            sage: (next(it), next(it), next(it), next(it), next(it))
            ([], [1], [1, 2], [2, 1], [1, 1])
            sage: next(it).parent()
            Packed words
            sage: P([])
            []
            sage: P is PackedWords_all()
            True
            sage: TestSuite(P).run()
        """
        Parent.__init__(self, category = InfiniteEnumeratedSets())
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), PackedWords_size),
            facade = True, keepkey = False)

    def _repr_(self):
        r"""
        TESTS::

            sage: PackedWords()
            Packed words
        """
        return "Packed words"

    def __contains__(self, w):
        r"""
        TESTS::

            sage: P = PackedWords()
            sage: 1 in P
            False
            sage: PackedWord([]) in P
            True
            sage: [1, 1, 4, 2, 3] in P
            True

        """
        if isinstance(w, self.element_class):
            return True
        try:
            self(w)
            return True
        except:
            return False

    def __call__(self, x = None, *args, **keywords):
        r"""
        Ensure that ``None`` instead of ``0`` is passed by default.

        TESTS::

            sage: P = PackedWords()
            sage: P([])
            []
        """
        return super(PackedWords, self).__call__(x, *args, **keywords)

    def _element_constructor_(self, *args, **keywords):
        r"""
        EXAMPLES::

            sage: P = PackedWords()
            sage: P._element_constructor_()
            []
            sage: P([])
            []
            sage: P([1, 1, 1])
            [1, 1, 1]
        """
        return self.element_class(self, *args, **keywords)

    def combinatorial_class_of_size(self, size):
        r"""
        EXAMPLES::
        
            sage: P = PackedWords()
            sage: P.combinatorial_class_of_size(6) is PackedWords(6)
            True
            sage: P.combinatorial_class_of_size(-1)
            Traceback (most recent call last):
            ...
            AssertionError: size (-1) must be a positive integer
        """
        assert(size >= 0), "size (%d) must be a positive integer" % size
        return PackedWords(size)

    def permutation_to_packed_words(self, sigma):
        r"""
        Compute all packed words which give *sigma* by standardization.

        TESTS::

            sage: PW = PackedWords()
            sage: PW.permutation_to_packed_words(Permutation([3, 1, 2, 4]))
            [[2, 1, 1, 2], [2, 1, 1, 3], [3, 1, 2, 3], [3, 1, 2, 4]]
            sage: PW.permutation_to_packed_words(Permutation([1, 2, 3]))
            [[1, 1, 1], [1, 1, 2], [1, 2, 2], [1, 2, 3]]
        """
        return PackedWords_size(len(sigma)).permutation_to_packed_words(sigma)

    Element = PackedWord


#==============================================================================
# Enumerated set of packed words of a given size
#==============================================================================
class PackedWords_size(PackedWords):
    r"""
    TESTS::

        sage: from sage.combinat.packed_words import PackedWords_size
        sage: for i in range(6): TestSuite(PackedWords_size(i)).run()

    """

    def __init__(self, size):
        super(PackedWords_size, self).__init__(category=FiniteEnumeratedSets())
        self._size = size

    def _repr_(self):
        r"""
        TESTS::

            sage: PackedWords(4)
            Packed words of size 4
        """
        return "Packed words of size %s" % (self._size)

    def __contains__(self, x):
        r"""
        TESTS::

            sage: P = PackedWords(4)
            sage: 1 in P
            False
            sage: PackedWord([]) in P
            False
            sage: PackedWord([1, 2, 1, 3, 1, 4]) in P
            False
            sage: PackedWord([1, 2, 1, 3]) in P
            True
        """
        return isinstance(x, self.element_class) and x.size() == self._size

    def _an_element_(self):
        r"""
        TESTS::

            sage: PackedWords(6).an_element()
            [1, 2, 3, 4, 5, 6]
        """
        return self.first()

    def cardinality(self):
        r"""
        Cardinality of Packed Words of size n

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
        r"""
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
                yield PackedWord(OrderedSetPartition.to_packed_word(part))

    @lazy_attribute
    def _parent_for(self):
        r"""
        The parent of the element generated by ``self`

        TESTS::

            sage: P = PackedWords(4)
            sage: P._parent_for
            Packed words
        """
        return PackedWords_all()

    @lazy_attribute
    def element_class(self):
        r"""
        TESTS::

            sage: P = PackedWords(4)
            sage: P.element_class
            <class 'sage.combinat.packed_words.PackedWords_all_with_category.element_class'>
            sage: P.first().__class__ == PackedWords().first().__class__
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, *args, **keywords):
        r"""
        EXAMPLES::

            sage: P = PackedWords(0)
            sage: P([])
            []

            sage: P([1])
            Traceback (most recent call last):
            ...
            ValueError: Wrong size of word
        """
        res = self.element_class(self._parent_for, *args, **keywords)
        if res.size() != self._size:
            raise ValueError("Wrong size of word")
        return res

    def permutation_to_packed_words(self, sigma):
        r"""
        Compute all packed words which give *sigma* by standardization.

        TESTS::

            sage: PW = PackedWords()
            sage: PW.permutation_to_packed_words(Permutation([3, 1, 2, 4]))
            [[2, 1, 1, 2], [2, 1, 1, 3], [3, 1, 2, 3], [3, 1, 2, 4]]
            sage: PW.permutation_to_packed_words(Permutation([1, 2, 3]))
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
