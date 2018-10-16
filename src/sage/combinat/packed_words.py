# -*- coding: utf-8 -*-
r"""
Packed Words

AUTHORS:

- Jean-Baptiste Priez
- Hugo Mlodecki
"""
#*****************************************************************************
#       Copyright (C) 2012 Jean-Baptiste Priez <jbp@kerios.fr>,
#                     2018 Hugo Mlodecki
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from six import add_metaclass

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import parent
from sage.structure.list_clone import ClonableIntArray
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.misc import uniq

from sage.rings.integer_ring import ZZ
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.family import Family
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
from sage.combinat.set_partition_ordered import OrderedSetPartition, OrderedSetPartitions
from sage.combinat.tools import transitive_ideal
from sage.combinat.composition import Composition
from sage.combinat.permutation import Permutations
from sage.combinat.permutation import to_standard
from sage.combinat.words.finite_word import evaluation_dict
from sage.combinat.combinatorial_map import combinatorial_map


@add_metaclass(InheritComparisonClasscallMetaclass)
class PackedWord(ClonableIntArray):
    r"""
    A packed word.

    A word `u` over the positive integers is a *packed word* if for each
    number `k > 1` appearing in `u`, the number `k - 1` also appears in `u`.

    .. SEEALSO:: :class:`PackedWords`

    EXAMPLES::

        sage: PackedWord([3, 4, 2, 2, 3, 5, 4, 2, 2, 3, 1, 5, 4, 1, 2])
        [3, 4, 2, 2, 3, 5, 4, 2, 2, 3, 1, 5, 4, 1, 2]
        sage: PackedWord([])
        []
        sage: PackedWord()
        []
        sage: PackedWord([1])
        [1]
        sage: PackedWord([2])
        Traceback (most recent call last):
        ...
        ValueError: [2] not in Packed words

    TESTS::

        sage: w = PackedWord()
        sage: TestSuite(w).run()
        sage: w = PackedWord([1, 3, 3, 2, 4, 1, 3])
        sage: TestSuite(w).run()
    """
    @staticmethod
    def __classcall_private__(cls, lst=[], check=True):
        r"""
        Ensure that packed words created by the enumerated sets and directly
        are the same and that they are instances of :class:`PackedWord`.

        TESTS::

            sage: from sage.combinat.packed_words import PackedWords_all
            sage: issubclass(PackedWords_all().element_class, PackedWord)
            True
            sage: w0 = PackedWord([4, 2, 3, 1, 2])
            sage: w0.parent()
            Packed words

            sage: w1 = PackedWords()([4, 2, 3, 1, 2])
            sage: w1.parent() is w0.parent()
            True
            sage: type(w1) is type(w0)
            True
        """
        P = PackedWords_all()
        return P(lst, check=check)


    def __init__(self, parent, lst, check=True):
        r"""
        Initialize ``self``.

        TESTS::

            sage: TestSuite(PackedWords()).run()  # long time
        """
        ClonableIntArray.__init__(self, parent, lst, check=check)
        self._max = 0 if not lst else max(lst)
        self._size = len(lst)

    def check(self):
        r"""
        Check that ``self`` is a packed word.

        TESTS::

            sage: PackedWord([3, 3, 2, 1])  # indirect doctest
            [3, 3, 2, 1]

            sage: PackedWord([2, 2, 1, 0, 4])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: [2, 2, 1, 0, 4] not in Packed words

            sage: PackedWords(3)([1,2])
            Traceback (most recent call last):
            ...
            ValueError: [1, 2] is not a packed word of size 3

            sage: PackedWords(3)([])
            Traceback (most recent call last):
            ...
            ValueError: [] is not a packed word of size 3
        """
        if self not in self.parent():
            raise ValueError("%s not in %s"%(self, self.parent()))

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        TESTS::

            sage: latex(PackedWord([1, 2, 3, 1, 1, 3]))
            123113
            sage: latex(PackedWord([1, 5, 10, 2, 3, 11, 4, 7, 6, 7, 2, 9, 8]))
            [1, 5, 10, 2, 3, 11, 4, 7, 6, 7, 2, 9, 8]
            sage: latex(PackedWord([]))
            \emptyset
        """
        if not self:
            return "\\emptyset"
        if self._max >= 10:
            return repr(self)
        return ''.join(repr(val) for val in self)

    @combinatorial_map(name='to ordered set partition')
    def to_ordered_set_partition(self):
        r"""
        Build an ordered set partition corresponding to ``self``.

        EXAMPLES::

            sage: PackedWord().to_ordered_set_partition()
            []
            sage: PackedWord([1]).to_ordered_set_partition()
            [{1}]
            sage: PackedWord([1, 2, 3, 1, 1, 3]).to_ordered_set_partition()
            [{1, 4, 5}, {2}, {3, 6}]
        """
        lst = [[]] * self._max
        for i, val in enumerate(self):
            lst[val-1] = lst[val-1] + [i+1]
        return OrderedSetPartition(lst)

    @combinatorial_map(name='to composition')
    def to_composition(self):
        r"""
        Return the compostion associated to ``self``.

        Given a packed word `u` with greatest letter `\ell`, we construct a
        composition `\alpha = (\alpha_1, \alpha_2, \dotsc, \alpha_{\ell})`
        by `\alpha_i` being the multiplicity of `i` in `u`.

        EXAMPLES::

            sage: PackedWord([]).to_composition()
            []
            sage: PackedWord([1, 1, 1]).to_composition()
            [3]
            sage: PackedWord([1, 2, 1]).to_composition()
            [2, 1]
            sage: PackedWord([1, 2, 3, 1, 1, 3]).to_composition()
            [3, 1, 2]
        """
        if not self:
            return Composition([])
        d = evaluation_dict(self)
        return Composition([d[i+1] for i in range(self._max)])

    def is_empty(self):
        r"""
        Return whether ``self`` is the empty word.

        EXAMPLES::

            sage: PackedWord().is_empty()
            True
            sage: PackedWord([]).is_empty()
            True
            sage: PackedWord([2, 1, 2]).is_empty()
            False
        """
        return not self

    def __add__(self, pw):
        r"""
        Return the concatenation of the two packed words.

        EXAMPLES::

            sage: pw1 = PackedWord([1, 3, 2, 1])
            sage: pw2 = PackedWord([1, 2, 1])
            sage: pw1 + pw2
            [1, 3, 2, 1, 1, 2, 1]
            sage: pw2 + pw1
            [1, 2, 1, 1, 3, 2, 1]
        """
        return PackedWord([i for i in self] + [i for i in pw])

    def __mul__(self, pw):
        r"""
        Return the usual composition of surjections if ``size(pw) == max(self)``.

        The right element can be a permutation.

        EXAMPLES::

            sage: pw1 = PackedWord([2, 1, 2])
            sage: pw2 = PackedWord([2, 1, 3, 2, 1])
            sage: p = Permutation([1, 3, 2])
            sage: pw2 * p
            Traceback (most recent call last):
            ...
            ValueError: the maximum value of the right packed word must be equal to the size of the left packed word

            sage: pw1 * PackedWord([1, 2, 3])
            [2, 1, 2]
            sage: pw1 * PackedWord([3, 2, 1])
            [2, 1, 2]
            sage: pw1 * p
            [2, 2, 1]
            sage: pw1 * pw2
            [1, 2, 2, 1, 2]
            sage: pw2 * PackedWord([1, 2, 3, 4, 5, 4])
            [2, 1, 3, 2, 1, 2]
            sage: pw2 * PackedWord([5, 4, 3, 2, 2, 1])
            [1, 2, 3, 1, 1, 2]
        """
        if pw in Permutations():
            m = pw.size()
        elif pw in PackedWords():
            m = pw._max
        if self.size() == m:
            return PackedWord([self[i-1] for i in pw])

        raise ValueError("the maximum value of the right packed word must be equal to the size of the left packed word")

    def size(self):
        r"""
        Return the size of ``self``.

        EXAMPLES::

            sage: PackedWord().size()
            0
            sage: PackedWord([2, 1, 1]).size()
            3
        """
        return self._size

    def max(self):
        r"""
        Return the maximum value of ``self``.

        It is also the length of the corresponding ordered set partition.

        EXAMPLES::

            sage: pw=PackedWord([2, 4, 3, 1, 2, 1])
            sage: pw.max()
            4
            sage: pw.to_ordered_set_partition().length()
            4
        """
        return self._max

    def reverse(self):
        """
        Return the packed word obtained by reversing the list ``self``.

        EXAMPLES::

            sage: PackedWord([]).reverse()
            []
            sage: PackedWord([3, 2, 1, 2]).reverse()
            [2, 1, 2, 3]
            sage: PackedWord([1, 2, 3, 4, 5]).reverse()
            [5, 4, 3, 2, 1]
        """
        return self.__class__(self.parent(), [i for i in reversed(self)])

    @combinatorial_map(name='complement')
    def complement(self):
        r"""
        Return the image of ``self`` under the complement involution
        on packed words.

        Given a packed word `u` of size `n` and greatest letter `\ell`,
        its *complement* is the word `[\ell+1-u_1, \ldots, \ell+1-u_n]`.

        EXAMPLES::

            sage: PackedWord([]).complement()
            []
            sage: PackedWord([1, 2, 1]).complement()
            [2, 1, 2]
            sage: PackedWord([1, 2, 3, 1, 1, 4]).complement()
            [4, 3, 2, 4, 4, 1]
        """
        if not self:
            return self
        return self.__class__(self.parent(), [self._max + 1 - i for i in self])

    def global_descents(self, final_descent=False, from_zero=False):
        r"""
        Return the list of global descents of ``self``.

        A *global descent* of a packed word `u` is a position (index) `d`
        of `u` such that `u_i > u_j` for all `i \leq d < j.`

        .. WARNING::

            By default, the positions of ``self`` start at `1`. If you want
            them to start at `0`, set the keyword ``from_zero`` to ``True``.

        .. SEEALSO:: :meth:`global_descents_factorization`

        INPUT:

        - ``final_descent`` -- boolean (default ``False``);
          if ``True``, the last position of a non-empty
          packed word is also considered as a global descent

        - ``from_zero`` -- boolean (default ``False``);
          if ``True``, return the positions starting from `0`

        EXAMPLES::

            sage: PackedWord([]).global_descents()
            []
            sage: PackedWord([1, 2, 1]).global_descents()
            []
            sage: PackedWord([2, 1]).global_descents()
            [1]
            sage: PackedWord([7, 5, 4, 6, 3, 1, 2, 1]).global_descents()
            [1, 4, 5]
            sage: PackedWord([7, 5, 4, 6, 3, 1, 2, 1]).global_descents(from_zero=True)
            [0, 3, 4]
            sage: PackedWord([7, 5, 4, 6, 3, 1, 2, 1]).global_descents(final_descent=True)
            [1, 4, 5, 8]
        """
        if not self:
            return []
        g_descents = []
        local_left_min = self._max
        for i in range(self._size - 1):
            local_left_min = min(local_left_min,self[i])
            if local_left_min > max(self[i+1::] + [0]):
                g_descents.append(i + 1)

        if final_descent:
            g_descents.append(self._size)

        if from_zero:
            return [i - 1 for i in g_descents]

        return g_descents

    def global_descents_factorization(self):
        r"""
        Return the list of packed words resulting from the factorization
        of ``self`` at global descents.

        .. SEEALSO:: :meth:`global_descents`

        EXAMPLES::

            sage: PackedWord([]).global_descents_factorization()
            [[]]
            sage: PackedWord([1, 2, 1]).global_descents_factorization()
            [[1, 2, 1]]
            sage: PackedWord([5, 4, 3, 4, 1, 2, 1]).global_descents_factorization()
            [[1], [2, 1, 2], [1, 2, 1]]
            sage: PackedWord([5, 4, 3, 4, 1, 2, 1, 4]).global_descents_factorization()
            [[1], [4, 3, 4, 1, 2, 1, 4]]
        """
        g_descents = self.global_descents(final_descent=True)
        if not g_descents:
            return [self]
        i = g_descents[0]
        g_d_f = [PackedWords.pack(self[:i])]
        for j in g_descents[1:]:
            g_d_f.append(PackedWords.pack(self[i:j]))
            i=j
        return g_d_f

    def global_ascents(self, initial_ascent=False, from_zero=False):
        r"""
        Return the list of the global ascents of ``self``.

        A *global ascent* of a packed word `u` is a position (index) `d`
        of `u` such that `u_i < u_j` for all `i < d \leq j`.

        .. WARNING::

            By default, the positions of ``self`` start at `1`. If you want
            them to start at `0`, set the keyword ``from_zero`` to ``True``.

        .. SEEALSO:: :meth:`global_ascents_factorization`

        INPUT:

        - ``initial_ascent`` -- boolean (default ``False``);
          if ``True``, the first position of a non-empty
          packed word is also considered as a global ascent

        - ``from_zero`` -- boolean (default ``False``);
          if ``True``, return the positions starting from `0`

        EXAMPLES::

            sage: PackedWord([]).global_ascents()
            []
            sage: PackedWord([1, 2, 1]).global_ascents()
            []
            sage: PackedWord([2, 1]).global_ascents()
            []
            sage: PackedWord([1, 2]).global_ascents()
            [2]
            sage: PackedWord([3, 1, 2, 1, 4, 6, 6, 5, 7, 8, 7]).global_ascents()
            [5, 6, 9]
            sage: PackedWord([3, 1, 2, 1, 4, 6, 6, 5, 7, 8, 7]).global_ascents(from_zero=True)
            [4, 5, 8]
            sage: PackedWord([3, 1, 2, 1, 4, 6, 6, 5, 7, 8, 7]).global_ascents(initial_ascent=True)
            [1, 5, 6, 9]
        """
        if not self:
            return []
        g_ascents = []
        local_left_max = 0
        for i in range(1, self._size):
            local_left_max = max(local_left_max, self[i-1])
            if local_left_max < min(self[i:] + [self._max]):
                g_ascents.append(i + 1)

        if initial_ascent:
            g_ascents.insert(0, 1)

        if from_zero:
            return [i - 1 for i in g_ascents]

        return g_ascents

    def global_ascents_factorization(self):
        r"""
        Return the list of packed words resulting from the factorization
        of ``self`` at global ascents.

        .. SEEALSO:: :meth:`global_ascents`

        EXAMPLES::

            sage: PackedWord([]).global_ascents_factorization()
            [[]]
            sage: PackedWord([1, 2, 1]).global_ascents_factorization()
            [[1, 2, 1]]
            sage: PackedWord([3, 1, 2, 1, 4]).global_ascents_factorization()
            [[3, 1, 2, 1], [1]]
            sage: PackedWord([3, 1, 2, 1, 4, 6, 6, 5, 7, 8, 7]).global_ascents_factorization()
            [[3, 1, 2, 1], [1], [2, 2, 1], [1, 2, 1]]
            sage: PackedWord([3, 1, 2, 1, 4, 6, 6, 5, 7, 8, 7, 4]).global_ascents_factorization()
            [[3, 1, 2, 1], [1, 3, 3, 2, 4, 5, 4, 1]]
        """
        g_a = self.global_ascents(from_zero=True, initial_ascent=True)
        if not g_a:
            return [self]
        g_a.append(self._size)
        return [PackedWords.pack(self[g_a[i]:g_a[i+1]]) for i in range(len(g_a)-1)]


    def inversions(self, side="right", support=None, from_zero=False):
        r"""
        Return the (``side``) weak order inversions (on ``support``)  of ``self``.

        .. RUBRIC:: Right inversions on positions

        The behavior when ``side`` is "right" and ``support`` is "position".

        This is the default behavior when no keywords are specified, and also
        when only the ``side`` keyword is specified.

        Let `u` be a packed word of size `n`. The *right weak order
        inversions* on *positions* of `u` are the pairs `(i, j)` such that
        `1 \leq i < j \leq n` and `u_i > u_j`.

        .. RUBRIC:: Right inversions on values

        The behavior when ``side`` is "right" and ``support`` is "value".

        The *right weak order inversions* on *values* of a packed word `u`
        are the pairs `(u_i, u_j)` such that `u_i > u_j` for some `i < j`.
   
        .. RUBRIC:: Left inversions on values

        The behavior when ``side`` is "left" and ``support`` is "value".

        This is the default behavior when only the ``side`` keyword is specified.

        The *left weak order inversions* on *values* of a packed word `u`
        are the pairs `(b, a)` such that `a < b` and the first occurence of
        `a` in `u` is after the last occurrence of `b` in `u`.

        .. RUBRIC:: Left inversions on positions

        The behavior when ``side`` is "left" and ``support`` is "position".

        The *left weak order inversions* on *positions* of a packed word `u`
        are the pairs `(i, j)` such that `i < j` and the first occurence of
        `u_j` in `u` is after the last occurrence of `u_i` in `u`.

        .. WARNING::

            By default, the inversions of ``self`` start at `1`. If you want
            them to start at `0` for support "position", set the keyword
            ``from_zero`` to ``True``. If the support is "value", then the
            keyword ``from_zero`` does not affect the result.

        INPUT:

        - ``side`` -- ``"left"`` or ``"right"``; the side of the weak order of inversions
        - ``support`` -- ``"position"`` or ``"value"``; the support of the result

        The ``right`` inversions are taken on ``positions`` by default,
        whereas the ``left`` inversions are taken by default on ``values``.

        Exemples when ``side`` is "right" and ``support`` is "position".
        
        EXAMPLES::

            sage: PackedWord([]).inversions()
            set()
            sage: PackedWord([1, 2, 3]).inversions()
            set()
            sage: PackedWord([1, 2, 1]).inversions()
            {(2, 3)}
            sage: PackedWord([3, 2, 1]).inversions()
            {(1, 2), (1, 3), (2, 3)}
            sage: PackedWord([3, 2, 1]).inversions(from_zero = True)
            {(0, 1), (0, 2), (1, 2)}
            sage: PackedWord([2, 3, 4, 1, 2, 4, 3]).inversions()
            {(1, 4), (2, 4), (2, 5), (3, 4), (3, 5), (3, 7), (6, 7)}

        Exemples when ``side`` is "right" and ``support`` is "value".

        EXAMPLES::

            sage: PackedWord([]).inversions(support = "value")
            set()
            sage: PackedWord([1, 2, 3]).inversions(support = "value")
            set()
            sage: PackedWord([1, 2, 1]).inversions(support = "value")
            {(2, 1)}
            sage: PackedWord([3, 2, 1]).inversions(support = "value")
            {(2, 1), (3, 1), (3, 2)}
            sage: PackedWord([3, 2, 1]).inversions(support = "value", from_zero = True)
            {(2, 1), (3, 1), (3, 2)}
            sage: PackedWord([2, 3, 4, 1, 2, 4, 3]).inversions(support = "value")
            {(2, 1), (3, 1), (3, 2), (4, 1), (4, 2), (4, 3)}

        Exemples when ``side`` is "left" and ``support`` is "value".

        EXAMPLES::

            sage: PackedWord([]).inversions(side = "left")
            set()
            sage: PackedWord([1, 2, 3]).inversions(side = "left")
            set()
            sage: PackedWord([1, 2, 1]).inversions(side = "left")
            set()
            sage: PackedWord([3, 1, 2]).inversions(side = "left")
            {(3, 1), (3, 2)}
            sage: PackedWord([3, 1, 2]).inversions(side = "left", from_zero = True)
            {(3, 1), (3, 2)}
            sage: PackedWord([3, 1, 4, 1, 2]).inversions(side = "left")
            {(3, 1), (3, 2), (4, 2)}

        Exemples when ``side`` is "left" and ``support`` is "position".

        EXAMPLES::

            sage: PackedWord([]).inversions(side = "left", support = "position")
            set()
            sage: PackedWord([1, 2, 3]).inversions(side = "left", support = "position")
            set()
            sage: PackedWord([1, 2, 1]).inversions(side = "left", support = "position")
            set()
            sage: PackedWord([3, 1, 2]).inversions(side = "left", support = "position")
            {(1, 2), (1, 3)}
            sage: PackedWord([3, 1, 2]).inversions(side = "left", support = "position", from_zero = True)
            {(0, 1), (0, 2)}
            sage: PackedWord([3, 1, 4, 1, 2]).inversions(side = "left",support = "position")
            {(1, 2), (1, 5), (3, 5)}
        """
        if not side in ["left", "right"]:
            raise ValueError("option 'side' must be 'left' or 'right'")
        if support == None:
            support = "position" if side == "right" else "value"
        if not support in ["position", "value"]:
            raise ValueError("option 'support' must be 'position' or 'value'")

        if not self:
            return set()

        n = self._size
        m = self._max

        if side == "right":

            if support == "position":
                return set((i + 1 - from_zero, j + 1 - from_zero)
                           for i in range(n - 1)
                           for j in range(i + 1, n)
                           if self[i] > self[j])

            if support == "value":
                return set((self[i], self[j])
                           for i in range(n - 1)
                           for j in range(i + 1, n)
                           if self[i] > self[j])

        if side == "left":
            rev = self[::-1]

            if support == "value":
                return set((j, i)
                           for i in range(1, m)
                           for j in range(i + 1, m + 1)
                           if self.index(i) > n - rev.index(j) - 1)

            if support == "position":
                return set((n - rev.index(j) - from_zero, self.index(i) + 1 - from_zero)
                           for i in range(1, m)
                           for j in range(i + 1, m + 1)
                           if self.index(i) > n - rev.index(j) - 1)

    ###################     Right Weak Order     ###############################

    def right_weak_order_succ(self):
        r"""
        Return the list of successors of ``self`` under the right weak order.

        For the right weak order on packed words of size `n`, we say `u` is a
        *right successor* of `v` if: (1) there exists a position `i < n - 1` such
        that `u` and `v` agree at all positions except `i` and `i + 1`, where the
        values are switched; and (2) `u` has one more right inversion on positions
        than `v` has.

        EXAMPLES::

            sage: PackedWord([]).right_weak_order_succ()
            []
            sage: PackedWord([1, 1, 1]).right_weak_order_succ()
            []

            sage: v = PackedWord([1, 2, 1])
            sage: u, = v.right_weak_order_succ(); u
            [2, 1, 1]
            sage: v.inversions()
            {(2, 3)}
            sage: u.inversions()
            {(1, 2), (1, 3)}

            sage: PackedWord([3, 1, 2]).right_weak_order_succ()
            [[3, 2, 1]]
            sage: PackedWord([3, 1, 2, 1, 1, 2, 4]).right_weak_order_succ()
            [[3, 2, 1, 1, 1, 2, 4], [3, 1, 2, 1, 2, 1, 4], [3, 1, 2, 1, 1, 4, 2]]
        """
        succ = []
        n = self._size
        P = parent(self)
        for i in range(n - 1):
            if self[i] < self[i+1]:
                p = self[:i] + [self[i+1], self[i]] + self[i+2:]
                succ.append(P.element_class(P, p, check=False))
        return succ

    def right_weak_order_pred(self):
        r"""
        Return the list of predecessors of ``self`` under the right weak order.

        For the right weak order on packed words of size `n`, we say `v` is a
        *right predecessor* of `u` if: (1) there exists a position `i < n - 1` such
        that `u` and `v` agree at all positions except `i` and `i + 1`, where the
        values are switched; and (2) `v` has one fewer right inversion on positions
        than `u` has.

        EXAMPLES::

            sage: PackedWord([]).right_weak_order_pred()
            []
            sage: PackedWord([1, 1, 1]).right_weak_order_pred()
            []

            sage: u = PackedWord([1, 2, 1])
            sage: v, = u.right_weak_order_pred(); v
            [1, 1, 2]
            sage: u.inversions()
            {(2, 3)}
            sage: v.inversions()
            set()

            sage: PackedWord([3, 1, 2]).right_weak_order_pred()
            [[1, 3, 2]]
            sage: PackedWord([3, 1, 2, 1, 1, 2, 4]).right_weak_order_pred()
            [[1, 3, 2, 1, 1, 2, 4], [3, 1, 1, 2, 1, 2, 4]]
        """
        pred = []
        n = self._size
        P = parent(self)
        for i in range(n - 1):
            if self[i] > self[i+1]:
                p = self[:i] + [self[i+1], self[i]] + self[i+2:]
                pred.append(P.element_class(P, p, check=False))
        return pred

    def right_weak_order_smaller(self):
        r"""
        Return the list of smaller or equal packed words of ``self``
        under the right weak order.

        .. SEEALSO:: :meth:`right_weak_order_pred`

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
        Return the list of greater or equal packed words of ``self``
        under the right weak order.

        .. SEEALSO:: :meth:`right_weak_order_succ`

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

    def right_weak_order_interval(self, other):
        r"""
        Return the list of packed words in the interval
        between `self` and `pw` under the right weak order.

        EXAMPLES::

            sage: P = PackedWord
            sage: P([1, 2, 1]).right_weak_order_interval(P([1, 2, 1, 3]))
            []
            sage: P([1, 2, 1]).right_weak_order_interval(P([1, 2, 1]))
            [[1, 2, 1]]
            sage: P([1, 2, 1]).right_weak_order_interval(P([2, 1, 1]))
            [[1, 2, 1], [2, 1, 1]]
            sage: P([2, 1, 1]).right_weak_order_interval(P([1, 2, 1]))
            [[1, 2, 1], [2, 1, 1]]
            sage: P([1, 2, 1]).right_weak_order_interval(P([2, 1, 2]))
            []
            sage: P([1, 2, 1, 3, 4, 1]).right_weak_order_interval(P([3, 1, 2, 1, 1, 4]))
            []
            sage: P([1, 2, 1, 3, 4, 1]).right_weak_order_interval(P([3, 1, 2, 4, 1, 1]))
            [[1, 2, 1, 3, 4, 1],
             [1, 2, 3, 1, 4, 1],
             [1, 3, 2, 1, 4, 1],
             [1, 2, 3, 4, 1, 1],
             [3, 1, 2, 1, 4, 1],
             [1, 3, 2, 4, 1, 1],
             [3, 1, 2, 4, 1, 1]]
        """
        if other in self.parent():
            other = self.parent()(other)
        else:
            raise ValueError("{0} is not a member of {1}".format(other, self.parent()))

        res = []
        up = lambda a: a.right_weak_order_succ()
        if self.is_gequal(other):
            UP_other = RecursivelyEnumeratedSet([other], up, structure="graded")
            for x in UP_other:
                if x.is_lequal(self):
                    res.append(x)
                elif len(x.inversions()) > len(self.inversions()):
                    break
                    
        elif other.is_gequal(self):
            UP_self = RecursivelyEnumeratedSet([self], up, structure="graded")
            for x in UP_self:
                if x.is_lequal(other):
                    res.append(x)
                elif len(x.inversions()) > len(other.inversions()):
                    break

        return res

    ###################     Left Weak Order     ################################

    def left_weak_order_succ(self):
        r"""
        Return the list of successors of ``self`` under the left weak order.

        For the left weak order on packed words, we say `u` is a *left successor*
        of `v` if: (1) `u` is built from `v` by switching all occurrences of the
        values `i` and `i+1`, for some `i` less than the maximum letter in `v`;
        and (2) `u` has one more left inversion on values than `v` has.

        EXAMPLES::

            sage: PackedWord([]).left_weak_order_succ()
            []
            sage: PackedWord([1, 1, 1]).left_weak_order_succ()
            []
            sage: PackedWord([1, 2, 1]).left_weak_order_succ()
            []

            sage: v = PackedWord([3, 1, 2])
            sage: u, = v.left_weak_order_succ(); u
            [3, 2, 1]
            sage: v.inversions(side="left")
            {(3, 1), (3, 2)}
            sage: u.inversions(side="left")
            {(2, 1), (3, 1), (3, 2)}

            sage: PackedWord([1, 2, 4, 3, 3, 2]).left_weak_order_succ()
            [[2, 1, 4, 3, 3, 1]]
            sage: PackedWord([1, 2, 4, 3, 3]).left_weak_order_succ()
            [[2, 1, 4, 3, 3], [1, 3, 4, 2, 2]]
        """
        if not self:
            return []

        succ = []
        m = self._max
        for i in range(1, m):
            if self._size - 1 - self[::-1].index(i) < self.index(i + 1):
                l = []
                for x in self:
                    if x == i:
                        l.append(i + 1)
                    elif x == i + 1:
                        l.append(i)
                    else:
                        l.append(x)
                succ.append(l)
        P = parent(self)
        return [P.element_class(P, p, check=False) for p in succ]

    def left_weak_order_pred(self):
        r"""
        Return the list of predecessors of ``self`` under the left weak order.

        For the left weak order on packed words, we say `v` is a *left predecessor*
        of `u` if: (1) `v` is built from `u` by switching all occurrences of the
        values `i` and `i+1`, for some `i` less than the maximum letter in `u`;
        and (2) `v` has one fewer left inversion on values than `u` has.

        EXAMPLES::

            sage: PackedWord([]).left_weak_order_pred()
            []
            sage: PackedWord([1, 1, 1]).left_weak_order_pred()
            []
            sage: PackedWord([1, 2, 1]).left_weak_order_pred()
            []

            sage: u = PackedWord([3, 1, 2])
            sage: v, = u.left_weak_order_pred(); v
            [2, 1, 3]
            sage: u.inversions(side="left")
            {(3, 1), (3, 2)}
            sage: v.inversions(side="left")
            {(2, 1)}

            sage: PackedWord([3, 1, 2, 4, 4]).left_weak_order_pred()
            [[2, 1, 3, 4, 4]]
            sage: PackedWord([3, 1, 3, 1, 2, 2, 2]).left_weak_order_pred()
            [[2, 1, 2, 1, 3, 3, 3]]
        """
        if not self:
            return []

        pred = []
        m = self._max
        for i in range(1, m):
            if self.index(i) > self._size - 1 - self[::-1].index(i + 1):
                l = []
                for x in self:
                    if x == i:
                        l.append(i + 1)
                    elif x == i + 1:
                        l.append(i)
                    else:
                        l.append(x)
                pred.append(l)
        P = parent(self)
        return [P.element_class(P, p, check=False) for p in pred]

    def left_weak_order_smaller(self):
        r"""
        Return the list of smaller or equal packed words of ``self``
        under the left weak order.

        .. SEEALSO:: :meth:`left_weak_order_pred`

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
        Return the list of greater or equal packed words of ``self``
        under the left weak order.

        .. SEEALSO:: :meth:`left_weak_order_succ`

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

    def left_weak_order_interval(self, other):
        r"""
        Return the list of packed words in the interval
        between `self` and `pw` under the let weak order.

        EXAMPLES::

            sage: P = PackedWord
            sage: P([1, 2, 1]).left_weak_order_interval(P([1, 2, 1, 3]))
            []
            sage: P([1, 2, 1]).left_weak_order_interval(P([1, 2, 1]))
            [[1, 2, 1]]
            sage: P([1, 2, 1]).left_weak_order_interval(P([2, 1, 1]))
            []
            sage: P([1, 2, 1]).left_weak_order_interval(P([2, 1, 2]))
            []
            sage: P([1, 1, 2]).left_weak_order_interval(P([2, 2, 1]))
            [[1, 1, 2], [2, 2, 1]]
            sage: P([2, 2, 1]).left_weak_order_interval(P([1, 1, 2]))
            [[1, 1, 2], [2, 2, 1]]
            sage: P([1, 4, 2, 3, 2, 3]).left_weak_order_interval(P([4, 3, 1, 2, 1, 2]))
            [[1, 4, 2, 3, 2, 3],
             [2, 4, 1, 3, 1, 3],
             [3, 4, 1, 2, 1, 2],
             [4, 3, 1, 2, 1, 2]]
        """
        if other in self.parent():
            other = self.parent()(other)
        else:
            raise ValueError("{0} is not a member of {1}".format(other, self.parent()))

        res = []
        up = lambda a: a.left_weak_order_succ()
        if self.is_gequal(other, side="left"):
            UP_other = RecursivelyEnumeratedSet([other], up, structure="graded")
            for x in UP_other:
                if x.is_lequal(self, side="left"):
                    res.append(x)
                elif len(x.inversions(side="left")) > len(self.inversions(side="left")):
                    break
                    
        elif other.is_gequal(self, side="left"):
            UP_self = RecursivelyEnumeratedSet([self], up, structure="graded")
            for x in UP_self:
                if x.is_lequal(other, side="left"):
                    res.append(x)
                elif len(x.inversions(side="left")) > len(other.inversions(side="left")):
                    break

        return res

    def is_gequal(self, pw, side="right"):
        r"""
        Return whether ``pw`` is greater than

            sage: u = PackedWord([3, 1, 1, 2])
            sage: v = PackedWord([3, 1, 2, 1])
            sage: w = PackedWord([3, 2, 2, 1])
            sage: v.is_gequal(v)
            True
            sage: v.is_gequal(u)
            True
            sage: v.is_gequal(w)
            False
            sage: w.is_gequal(u, side="left")
            True
            
            sage: u = PackedWords(3).random_element()
            sage: v = PackedWords(4).random_element()
            sage: v.is_gequal(u)
            False
            sage: v.is_gequal(u, side="left")
            False
            
            sage: u = PackedWords(5).random_element()
            sage: v = PackedWords(5).random_element()
            sage: v.is_gequal(u) == (v in u.right_weak_order_greater())
            True
            sage: v.is_gequal(u, side="left") == (v in u.left_weak_order_greater())
            True
        """
        if not side in ["left", "right"]:
            raise ValueError("option 'side' must be 'left' or 'right'")
        
        if side == "right" and self.to_composition() == pw.to_composition():
            std_self = PackedWord(to_standard(self))
            std_pw = PackedWord(to_standard(pw))
            return std_pw.inversions(side="left").issubset(std_self.inversions(side="left"))
        
        s1=set(self.to_ordered_set_partition())
        s2=set(pw.to_ordered_set_partition())
        if side == "left" and s1 == s2:
            return pw.inversions().issubset(self.inversions())
        
        return False
    
    def is_lequal(self, pw, side="right"):
        r"""
        Return whether ``pw`` is greater than

            sage: u = PackedWord([3, 1, 1, 2])
            sage: v = PackedWord([3, 1, 2, 1])
            sage: w = PackedWord([3, 2, 2, 1])
            sage: v.is_lequal(v)
            True
            sage: v.is_lequal(u)
            False
            sage: v.is_lequal(w)
            False
            sage: u.is_lequal(w, side="left")
            True
            
            sage: u = PackedWords(3).random_element()
            sage: v = PackedWords(4).random_element()
            sage: v.is_lequal(u)
            False
            sage: v.is_lequal(u, side="left")
            False
            
            sage: u = PackedWords(5).random_element()
            sage: v = PackedWords(5).random_element()
            sage: v.is_lequal(u) == (v in u.right_weak_order_greater())
            True
            sage: v.is_lequal(u, side="left") == (v in u.left_weak_order_greater())
            True
        """
        return pw.is_gequal(self, side)

#==============================================================================
# Parent classes
#==============================================================================

class PackedWords(UniqueRepresentation, Parent):
    r"""
    Packed words.

    A word `u` over the positive integers is a *packed word* if for each
    number `k > 1` appearing in `u`, the number `k - 1` also appears in `u`.

    TODO:

        - perhaps reformat the display of packed words of sizes 0..3.
          (I need to look at the documentation to see.)

    Packed words are in bijection with ordered set partitions as follows.
    Given an ordered set partition `O`, the corresponding packed word `u`
    is obtained by setting `u_i = j` if `i` belongs to the `j`-th block of `O`.

    Here are the Packed Words of sizes 0 to 3::

        \emptyset
        1
        11, 12, 21
        111, 112, 121, 211, 122, 212, 221, 123, 132, 213, 231, 312, 321

    INPUT:

    - ``size`` -- (optional) an integer

    EXAMPLES::

        sage: P = PackedWords(); P
        Packed words
        sage: P([])
        []
        sage: P([6, 2, 3, 3, 1, 2, 4, 2, 5, 2, 6])
        [6, 2, 3, 3, 1, 2, 4, 2, 5, 2, 6]

        sage: P = PackedWords(3); P
        Packed words of size 3
        sage: P.list()
        [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2],
         [3, 2, 1], [1, 2, 2], [2, 1, 2], [2, 2, 1], [1, 1, 2],
         [1, 2, 1], [2, 1, 1], [1, 1, 1]]
        sage: P.cardinality()
        13
        sage: O = OrderedSetPartitions(3)
        sage: O.cardinality()
        13

        sage: PackedWords(4)
        Packed words of size 4
    """
    @staticmethod
    def __classcall_private__(cls, n=None):
        r"""
        Construct the correct parent based upon input ``n``.

        TESTS::

            sage: from sage.combinat.packed_words import PackedWords_size, PackedWords_all
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

            sage: PackedWords(3/2)
            Traceback (most recent call last):
            ...
            ValueError: n must be a non-negative integer
        """
        if n is None:
            return PackedWords_all()
        if n not in ZZ or n < 0:
            raise ValueError("n must be a non-negative integer")
        return PackedWords_size(ZZ(n))

    @staticmethod
    def from_ordered_set_partition(osp):
        r"""
        Build a packed word corresponding to the ordered set partition ``osp``.

        INPUT::

        - ``osp`` -- an object or iterable (e.g., a list of lists)
          that defines an ordered set partition.

        EXAMPLES::

            sage: PackedWords.from_ordered_set_partition([])
            []
            sage: osp=[[1, 3],[2]]
            sage: PackedWords.from_ordered_set_partition(osp)
            [1, 2, 1]

            sage: pw=PackedWord([1, 4, 1, 3, 2])
            sage: PackedWords.from_ordered_set_partition(pw.to_ordered_set_partition()) == pw
            True
        """
        OSP = OrderedSetPartition
        return PackedWord(OSP(osp).to_packed_word())

    @staticmethod
    def pack(li):
        r"""
        Return the packed word associated to ``li``.

        This map is the analog for packed words of the standardization map
        for permutations.

        .. SEEALSO:: :meth:`sage.combinat.permutation.Permutation.to_standard`

        EXAMPLES::

            sage: PackedWords.pack([])
            []
            sage: PackedWords.pack([3, 1])
            [2, 1]
            sage: PackedWords.pack([1, 0, 0])
            [2, 1, 1]
            sage: PackedWords.pack([3, 1, 55])
            [2, 1, 3]
            sage: PackedWords.pack([11, 4, 1, 55])
            [3, 2, 1, 4]
            sage: PackedWords.pack([11, 4, 1, 11, 4])
            [3, 2, 1, 3, 2]
        """
        l = uniq(li)
        return PackedWord([l.index(i) + 1 for i in li])

    def permutation_to_packed_words(self, sigma):
        r"""
        Compute all packed words whose standardization is ``sigma``.

        If ``self._size`` is not None, first check that this size agrees
        with the size of ``sigma``.

        INPUT:

        - ``sigma`` -- a permutation of `n`

        EXAMPLES::

            sage: PW = PackedWords(4)
            sage: PW.permutation_to_packed_words(Permutation([3, 1, 2, 4]))
            [[2, 1, 1, 2], [2, 1, 1, 3], [3, 1, 2, 3], [3, 1, 2, 4]]

            sage: PW = PackedWords(3)
            sage: PW.permutation_to_packed_words(Permutation([1, 2, 3]))
            [[1, 1, 1], [1, 1, 2], [1, 2, 2], [1, 2, 3]]

        TESTS::

            sage: PackedWords(0).permutation_to_packed_words([])
            [[]]
            sage: PackedWords(1).permutation_to_packed_words([1])
            [[1]]

            sage: PW = PackedWords()
            sage: PW.permutation_to_packed_words([1,4,2,3])
            [[1, 2, 1, 1], [1, 3, 1, 2], [1, 3, 2, 2], [1, 4, 2, 3]]

            sage: PW = PackedWords(4)
            sage: PW.permutation_to_packed_words([1, 2, 3])
            Traceback (most recent call last):
            ...
            ValueError: [1, 2, 3] does not represent a packed word of size 4

            sage: PW.permutation_to_packed_words([1, 2, 4, 4])
            Traceback (most recent call last):
            ...
            ValueError: [1, 2, 4, 4] does not represent a permutation of 4

        """
        try:
            n = Permutations()(sigma).size()
        except ValueError:
            raise ValueError("{0} does not represent a permutation of {1}".format(sigma, len(sigma)))

        try:
            size = self._size
            if size != n:
                raise ValueError("{0} does not represent a packed word of size {1}".format(sigma, self._size))
        except AttributeError:
            size = None

        if size == 0:
            return [self.element_class(self, [], check=False)]
        if size == 1:
            return [self.element_class(self, [1], check=False)]

        li = [({sigma.index(1): 1}, sigma.index(1))]
        for i in range(2, n):
            index_i = sigma.index(i)
            tmp = []
            for (pw, l_index) in li:
                if l_index < index_i:
                    pw[index_i] = pw[l_index]
                    tmp.append((dict(pw), index_i))
                pw[index_i] = pw[l_index] + 1
                tmp.append((dict(pw), index_i))
            li = tmp
        index_i = sigma.index(n)
        res = []
        for (pw, l_index) in li:
            if l_index < index_i:
                pw[index_i] = pw[l_index]
                res.append(self.element_class(self, list(pw.values()), check=False))
            pw[index_i] = pw[l_index] + 1
            res.append(self.element_class(self, list(pw.values()), check=False))
        return res

    def __contains__(self, w):
        r"""
        Determine if ``w`` may be viewed as belonging to ``self``.

        TESTS::

            sage: P = PackedWords()
            sage: P([2,3,1]) in P
            True
            sage: PackedWord([]) in P
            True
            sage: P([]) in PackedWords(0)
            True

            sage: 1 in P
            False
            sage: [1, 1, 4, 2, 3] in P
            True
            sage: {2,3,4,1} in P  # perhaps undesireable answer
            True
        """
        try:
            # check for a list of integers
            if all(a in ZZ and ZZ(a) > 0 for a in w):
                s = sorted(set(int(a) for a in w))
                n = len(w)
            else:
                #raise ValueError("{} is not a packed word".format(w))
                return False
        except (ValueError, TypeError):
            # Elements not integers or w not iterable
            #raise ValueError("{} is not a packed word".format(w))
            return False

        # if parent has a size, it should agree with that of w
        size = getattr(self, "_size", None)
        if size and n != size:
            return False

        # the distinct letters in w should be {1, 2, ..., max(w)}
        if s != [i for i in range(1, max([0]+s[-1:])+1)]:
            return False

        return True

    Element = PackedWord

    def _element_constructor_(self, lst=[], check=True):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: P = PackedWords()
            sage: P4 = PackedWords(4)
            sage: P([])
            []
            sage: p = P([1, 3, 3, 2]); p
            [1, 3, 3, 2]
            sage: p.parent()
            Packed words
            sage: p4 = P4([1, 3, 3, 2]); p4
            [1, 3, 3, 2]
            sage: p4.parent()
            Packed words of size 4
            sage: P4([1, 3, 3, 2, 2])
            Traceback (most recent call last):
            ...
            ValueError: [1, 3, 3, 2, 2] is not a packed word of size 4
            sage: P4([1, 4, 4, 2])
            Traceback (most recent call last):
            ...
            ValueError: [1, 4, 4, 2] not in Packed words of size 4
            sage: P([1, 4, 4, 2])
            Traceback (most recent call last):
            ...
            ValueError: [1, 4, 4, 2] not in Packed words
        """
        # if parent has a size, it should agree with the size of lst
        size = getattr(self, "_size", None)
        if size and len(lst) != size:
            raise ValueError("{0} is not a packed word of size {1}".format(lst, size))

        if isinstance(parent(lst), self.__class__):
            return lst
        if not lst:
            return self.element_class(self, [], check=False)

        return self.element_class(self, lst, check=check)


#==============================================================================
# Enumerated set of all packed words
#==============================================================================

class PackedWords_all(PackedWords, DisjointUnionEnumeratedSets):
    """
    The set of all packed words of any size.

    EXAMPLES::

        sage: P = PackedWords()
        sage: P.cardinality()
        +Infinity
        sage: it = iter(P)
        sage: [next(it) for dummy in range(5)]
        [[], [1], [1, 2], [2, 1], [1, 1]]
        sage: P.an_element()
        []
        sage: P([])
        []
    """
    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: TestSuite(PackedWords()).run()  # long time
        """
        fam = Family(NonNegativeIntegers(), PackedWords_size)
        DisjointUnionEnumeratedSets.__init__(self, fam, facade=True, keepkey=False)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: PackedWords()
            Packed words
        """
        return "Packed words"

    def subset(self, size=None):
        r"""
        Return the set of packed words of size ``size``.

        EXAMPLES::

            sage: P = PackedWords()
            sage: P.subset(6) is PackedWords(6)
            True
            sage: P.subset(0)
            Packed words of size 0

        TESTS::

            sage: P.subset(-1)
            Traceback (most recent call last):
            ...
            ValueError: n must be a non-negative integer
            sage: P.subset(3/2)
            Traceback (most recent call last):
            ...
            ValueError: n must be a non-negative integer
        """
        if size is None:
            return self
        return PackedWords(size)


#==============================================================================
# Enumerated set of packed words of a given size
#==============================================================================

class PackedWords_size(PackedWords):
    r"""
    Packed words of a fixed size (or length) `n`.

    INPUT:

    - ``size`` -- integer; the size

    EXAMPLES:

    We create all packed words of size at most `3`::

        sage: P = PackedWords(0)
        sage: list(P)
        [[]]

        sage: P = PackedWords(1)
        sage: list(P)
        [[1]]

        sage: P = PackedWords(2)
        sage: list(P)
        [[1, 2], [2, 1], [1, 1]]

        sage: P = PackedWords(3)
        sage: list(P)
        [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2],
         [3, 2, 1], [1, 2, 2], [2, 1, 2], [2, 2, 1], [1, 1, 2],
         [1, 2, 1], [2, 1, 1], [1, 1, 1]]

    We can create specific packed words of a fixed size::

        sage: P = PackedWords(13)
        sage: P([2, 7, 5, 2, 3, 3, 4, 7, 1, 3, 6, 2, 5])
        [2, 7, 5, 2, 3, 3, 4, 7, 1, 3, 6, 2, 5]

    TESTS::

        sage: P = PackedWords(0)
        sage: P([])
        []

        sage: P([1])
        Traceback (most recent call last):
        ...
        ValueError: [1]  not in Packed words of size 0
    """
    def __init__(self, size):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: TestSuite(PackedWords(0)).run()
            sage: TestSuite(PackedWords(1)).run()
            sage: TestSuite(PackedWords(2)).run()
            sage: TestSuite(PackedWords(5)).run()  # long time
        """
        super(PackedWords_size, self).__init__(category=FiniteEnumeratedSets())
        self._size = size

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: PackedWords(4)
            Packed words of size 4
        """
        return "Packed words of size %s" % (self._size)

    def __contains__(self, x):
        r"""
        Determine if ``x`` may be viewed as belonging to ``self``.

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
        return PackedWords.__contains__(self, x) and len(x) == self._size

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        TESTS::

            sage: PackedWords(0).cardinality()
            1
            sage: PackedWords(1).cardinality()
            1
            sage: PackedWords(2).cardinality()
            3
            sage: PackedWords(3).cardinality()
            13
        """
        return OrderedSetPartitions(self._size).cardinality()

    def __iter__(self):
        r"""
        Iterate over ``self``.

        TESTS::

            sage: list(PackedWords(0))
            [[]]
            sage: list(PackedWords(1))
            [[1]]
            sage: list(PackedWords(2))
            [[1, 2], [2, 1], [1, 1]]
            sage: list(PackedWords(3))
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2],
             [3, 2, 1], [1, 2, 2], [2, 1, 2], [2, 2, 1], [1, 1, 2],
             [1, 2, 1], [2, 1, 1], [1, 1, 1]]
        """
        if self._size == 0:
            yield self.element_class(self, [], check=False)
            return

        OSP = OrderedSetPartitions(self._size)
        for osp in OSP:
            yield self.element_class(self, osp.to_packed_word(), check=False)

    def random_element(self):
        r"""
        Return a random element of ``self``.

        The packed words of size ``self._size`` are not covered with
        uniform probability, but they are all covered.

        TESTS::

            sage: PackedWords(0).random_element()
            []
            sage: PackedWords(1).random_element()
            [1]
            sage: PackedWords(5).random_element()  # random
            [1, 3, 3, 2, 4]
        """
        osp = OrderedSetPartitions(self._size).random_element()
        return self.element_class(self, osp.to_packed_word(), check=False)
