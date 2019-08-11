# -*- coding: utf-8 -*-
r"""
Semistandard set valued tableaux

AUTHORS:

- Jeremy Meza, Oliver Pechenik, Wencin Poh (2019): initial version
"""

#*****************************************************************************
#       Copyright (C) Jeremy Meza jdmeza@berkeley.edu
#                     Oliver Pechenik pechenik@umich.edu
#                     Wencin Poh wpoh@ucdavis.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import
from six.moves import range, zip, map
from six import add_metaclass
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.classical_crystals import ClassicalCrystals
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.partition import Partition, Partitions, _Partitions, Partitions_n
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.tableau import Tableau, Tableaux, SemistandardTableau, SemistandardTableaux
from sage.combinat.skew_tableau import SkewTableau, SemistandardSkewTableaux
from sage.rings.infinity import PlusInfinity
from sage.rings.integer import Integer
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc import powerset

#####################################
# Semistandard Set-Valued Tableaux  #
#####################################

class SemistandardSetValuedTableau(Tableau):
    r"""
    A semistandard set-valued tableau.

    A semistandard set-valued tableau is a tableau whose cells are filled with
    nonempty sets of integers such that:

    1. If cell `A` is to the left of cell `B` in the tableau, then 
       max(A) <= \min(B). 
    2. If cell `A` is above cell `B` in the tableau, then max(A) < min(B).

    EXAMPLES::

        sage: SSVT = SemistandardSetValuedTableaux([3,2])
        sage: SSVT([[[2,1], [2], [3,4]], [[3], [4,3]]])[1]
        ((3,), (3, 4))
        sage: T = SemistandardSetValuedTableau([[[2,1], [2], [3,4]], [[4,3], [4]]])
        sage: T[0]
        ((1, 2), (2,), (3, 4))
        sage: SemistandardSetValuedTableau([[[1], [1,3,2], [3], [4,3,5]], [[4], [5,4]], [[5]]])
        [[[1], [1, 2, 3], [3], [3, 4, 5]], [[4], [4, 5]], [[5]]]

    TESTS::

        sage: T = SemistandardSetValuedTableau([[[1,2], [1,3]], [[3]]])
        Traceback (most recent call last):
        ...
        ValueError: entries in each row of semistandard set-valued tableau must be weakly increasing
        sage: T = SemistandardSetValuedTableau([[[1,2], [2,3]], [[2,3]]])
        Traceback (most recent call last):
        ...
        ValueError: entries in each column of semistandard set-valued tableau must be strictly increasing
    """
    @staticmethod
    def __classcall_private__(self, t):
        """
        Ensure that a :class:`SemistandardSetValuedTableau` is only
        constructed as an ``element_class`` call of an appropriate parent.

        EXAMPLES::
        
            sage: t = SemistandardSetValuedTableau([[[1],[1],[8]],[[2],[6,2]],[[3,7,4]]])
            sage: t.shape()
            [3, 2, 1]

        TESTS::

            sage: SemistandardSetValuedTableau([[[1,2],[2,3]],[[4,6]]])
            [[[1, 2], [2, 3]], [[4, 6]]]
        """

        if isinstance(t, SemistandardSetValuedTableau):
            return t

        T = list([list(row) for row in t])
        for i in range(len(T)):
            for j in range(len(T[i])):
                T[i][j] =  tuple(sorted(T[i][j]))

        SSVT = SemistandardSetValuedTableaux(Tableau(T).shape())
        return SSVT.element_class(SSVT, T)

    def __init__(self, parent, t, check=True, preprocessed=False):
        """
        Initialize a semistandard set-valued tableau.

        TESTS::

            sage: T1 = SemistandardSetValuedTableau([[(1,2), [3,2]], [(4,3,5), [5]]])
            sage: T2 = SemistandardSetValuedTableaux([2,2])([[(1,2), [3,2]], [(4,3,5), [5]]])
            sage: T1 == T2
            True
            sage: T1.parent()
            Semistandard set-valued tableaux of shape [2, 2]
            sage: T2.parent()
            Semistandard set-valued tableaux of shape [2, 2]
            sage: T1 is T2  # identical shifted tableaux are distinct objects
            False

        A semistandard set-valued tableau is deeply immutable as the cells are
        stored as tuples::

            sage: T = SemistandardSetValuedTableau([[[1], [1,2]], [[2]]])
            sage: T[0][1] = (1,2,3)
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment

        """
        if not preprocessed:
            T = self._preprocess(t)
        Tableau.__init__(self, parent, T, check=check)

    @staticmethod
    def _preprocess(t):
        """
        Preprocess list ``t`` to initialize the tableau.

        The output is a list of list of sets (each set is represented by
        a sorted tuple).

        TESTS::

            sage: SemistandardSetValuedTableau._preprocess([[(1,2),[3,2]],[(4,6,5)]])
            [[(1, 2), (2, 3)], [(4, 5, 6)]]
            sage: SemistandardSetValuedTableau._preprocess([])
            []
        """
        if isinstance(t, SemistandardSetValuedTableau):
            return t
        # Preprocess list to represent sets as tuples.
        T = list([list(row) for row in t])
        for i in range(len(T)):
            for j in range(len(T[i])):
                T[i][j] = tuple(sorted(T[i][j]))
        return T

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: T = SemistandardSetValuedTableau([[(1,2),[3,2]],[(4,6,5)]])
            sage: T
            [[[1, 2], [2, 3]], [[4, 5, 6]]]
            sage: U = SemistandardSetValuedTableau([])
            sage: U
            []
        """
        L = [[list(cell) for cell in row] for row in self.to_list()]
        return repr(L)

    def _repr_tab(self):
        """
        Return a nested list of string representation of ``self``.

        EXAMPLES::

            sage: T1 = SemistandardSetValuedTableau([[[1], [1, 2]], [[2, 4], [4]]])
            sage: T1._repr_tab()
            [['1', '1, 2'], ['2, 4', '4']]
            sage: T2 = SemistandardSetValuedTableau([[[1], [1], [1], [1]], [[2], [2, 4, 5]], [[3, 4]]])
            sage: T2._repr_tab()
            [['1', '1', '1', '1'], ['2', '2, 4, 5'], ['3, 4']]
        """
        return [[repr(list(cell))[1:-1] for cell in row] for row in self]

    def _latex_(self):
        r"""
        Return LaTeX code for ``self``.

        EXAMPLES::

            sage: T = SemistandardSetValuedTableau([[[1], [1, 2]], [[2, 4], [4]]])
            sage: latex(T)
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{1}&\lr{1, 2}\\\cline{1-2}
            \lr{2, 4}&\lr{4}\\\cline{1-2}
            \end{array}$}
            }
        """
        from sage.combinat.output import tex_from_array
        L = [row for i, row in enumerate(self._repr_tab())]
        return tex_from_array(L)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: T1 = SemistandardSetValuedTableau([[[1], [1, 2]], [[2, 4], [4]]])
            sage: SSVT = SemistandardSetValuedTableaux([2,2],max_entry=4)
            sage: T2 = SSVT([[[1], [2, 1]], [[4, 2], [4]]])
            sage: hash(T1) == hash(T2)
            True
        """
        return hash(tuple(tuple(tuple(cell) for cell in row) for row in self))

    def check(self):
        """
        Check that ``self`` is a valid semistandard set-valued tableau.
        """
        super(SemistandardSetValuedTableau, self).check()

        # Tableau() has checked that t is tableau, so it remains to check that
        # the entries of t are positive integers which are weakly increasing
        # along rows
        for row in self:
            if any(max(row[c]) > min(row[c+1]) for c in range(len(row)-1)):
                raise ValueError("entries in each row of semistandard set-valued tableau must be weakly increasing")
        # and strictly increasing down columns
        if self:
            for row, next in zip(self, self[1:]):
                if not all(max(row[c]) < min(next[c]) for c in range(len(next))):
                    raise ValueError("entries in each column of semistandard set-valued tableau must be strictly increasing")

    def excess(self):
        """
        Return the excess statistic for ``self``.

        The excess of a semistandard set-valued tableau `T` is the total
        number of integers in `T` minus the size of `T`.

        EXAMPLES::

            sage: T = SemistandardSetValuedTableau([[[1,2],[2,3]],[[3,4,5]]])
            sage: T.excess()
            4

            sage: S = SemistandardSetValuedTableau([[[1],[3]],[[5]]])
            sage: S.excess()
            0
        """
        tot = sum([len(cells) for cells in self.entries()])
        return tot - self.size()

    def weight(self):
        r"""
        Return the weight of the set-valued tableau ``self``. 

        Trailing zeroes are omitted when returning the weight.

        The weight of a tableau `T` is the sequence `(a_1, a_2, a_3, \ldots )`,
        where `a_k` is the number of entries of `T` equal to `k`. This
        sequence contains only finitely many nonzero entries.

        .. WARNING::

            If tableau is considered as a crystal element, the weight returned
            will ignore all trailing zeroes.

            sage: SSVT = SemistandardSetValuedTableaux([2,1], max_entry=5)
            sage: T = SemistandardSetValuedTableau([[[1], [1,3]], [[3]]])
            sage: T.weight() # (2, 0, 2, 0, 0) is expected
            [2, 0, 2]

        EXAMPLES::

            sage: SemistandardSetValuedTableau([[[1],[1],[8]],[[2],[6,2]],[[3,7,4]]]).weight()
            [2, 2, 1, 1, 0, 1, 1, 1]
            sage: SemistandardSetValuedTableau([]).weight()
            []
        """
        if len(self) == 0:
            return []
        vec = sum([sum(row,()) for row in self],())
        m = max(vec)
        wt = [0] * m
        for i in vec:
            if i > 0:
                wt[i-1] += 1
        return wt

    @combinatorial_map(order=2, name='Bender-Knuth involution')
    def bender_knuth_involution(self, k, rows=None):
        r"""
        Return the image of ``self`` under the `k`-th Bender--Knuth involution,
        assuming ``self`` is a semistandard set-valued tableau. This function
        was introduced by Ikeda and Shimazaki [IS2014]_ in this context.

        Let `T` be a tableau and fix `k`. Then a free `k` in `T` means a cell of
        `T` that contains `k` and whose direct lower neighbor does not contain
        `k + 1` (in particular, this lower neighbor might not exist). A free
        `k + 1` in `T` is a cell of `T` that contains `k + 1` and whose direct
        upper neighbor does not contain `k` (in particular, this neighbor might
        not exist). Note that a cell that contains both `k` and `k + 1` is both
        a free `k` and a free `k + 1`. It is clear that for any row `r` of `T`,
        the free `k`'s and free `k + 1`'s in `r` together form a contiguous
        interval of `r`.

        The *`k`-th Bender--Knuth involution at row `i`* changes the entries of
        the cells in this interval in such a way that if it used to have
        `a` entries of `k` and `b` entries of `k + 1`, it will now have `b`
        entries of `k` and `a` entries of `k + 1`. For fixed `k`, the
        `k`-th Bender--Knuth switches for different `i` commute. The
        composition of the `k`-th Bender--Knuth switches for all rows is
        called the *`k`-th Bender--Knuth involution*. This is used to show that
        the symmetric Grothendieck polynomials defined as generating functions
        for semistandard set-valued tableaux are in fact symmetric polynomials.

        INPUT:

        - ``k`` -- an integer

        - ``rows`` -- (Default ``None``) When set to ``None``, the method
          computes the `k`-th Bender--Knuth involution as defined above.
          When an iterable, this computes the composition of the `k`-th
          Bender--Knuth switches at row `i` over all `i` in ``rows``. When set
          to an integer `i`, the method computes the `k`-th Bender--Knuth
          switch at row `i`. Note the indexing of the rows starts with `1`.

        OUTPUT:

        The image of ``self`` under either the `k`-th Bender--Knuth
        involution, the `k`-th Bender--Knuth switch at a certain row, or
        the composition of such switches, as detailed in the INPUT section.

        EXAMPLES::

            sage: t = SemistandardSetValuedTableau([[[1],[1,2,3]],[[4,6]]])
            sage: t
            [[[1], [1, 2, 3]], [[4, 6]]]
            sage: t.bender_knuth_involution(1)
            [[[1, 2], [2, 3]], [[4, 6]]]
            sage: t.bender_knuth_involution(1).bender_knuth_involution(2)
            [[[1, 2, 3], [3]], [[4, 6]]]

        The Bender--Knuth involution is an involution::

            sage: t = SemistandardSetValuedTableau([[[1],[1],[8]],[[2],[6,2]],[[3,7,4]]])
            sage: all(t.bender_knuth_involution(k).bender_knuth_involution(k) == t for k in range(1,8))
            True
        """
        ell = len(self)    # ell is the number of rows of self.
        # Sanitizing the rows input so that it always becomes a list of
        # nonnegative integers. We also subtract 1 from these integers
        # because the i-th row of a tableau T is T[i - 1].
        def rem(tup, num):
            mylist = list(tup)
            mylist.remove(num)
            return tuple(mylist)

        def app(tup, num):
            mylist = list(tup)
            mylist.append(num)
            return tuple(mylist)

        if rows is None:
            rows = list(range(ell))
        elif rows in ZZ:
            rows = [rows-1]
        else:
            rows = [i-1 for i in rows]
        # Now, rows should be iterable.

        # result_tab is going to be the result tableau (as a list of lists);
        # we will build it up step by step, starting with a deep copy of self.
        result_tab = self.to_list()
        for i in rows:
            if i >= ell:
                continue
            # Setup the previous and next rows
            if i == 0:
                prev_row = [[None]] * len(result_tab[i])
            else:
                prev_row = result_tab[i-1]
            if i == ell - 1:
                next_row = [[None]] * len(result_tab[i])
            else:
                next_row = result_tab[i+1] + [[None]] * (len(result_tab[i]) - len(result_tab[i+1]))
            a = 0 #counter for free k
            b = 0 #counter for free k+1
            c = 0 #counter for free box with both k and k+1
            sk = None # The column number of the first free k
            sk1 = None # The column number of the first free k+1
            skboth = None #The column number of the (unique) free box with both k and k+1
            for j, val in enumerate(result_tab[i]):
                if k in val and k+1 not in next_row[j]:
                    if k+1 in val:
                        c += 1
                        skboth = j
                    else:
                        if sk is None:
                            sk = j
                        a += 1
                elif k+1 in val and k not in prev_row[j]:
                    if sk1 is None:
                        sk1 = j
                    b += 1
            if skboth is None:
                if sk1 is not None:
                    if a > b:
                        for j in range(sk1-(a-b), sk1):
                            result_tab[i][j] = rem(result_tab[i][j],k)
                            result_tab[i][j] = app(result_tab[i][j],k+1)
                    elif a < b:
                        for j in range(sk1, sk1+b-a):
                            result_tab[i][j] = rem(result_tab[i][j],k+1)
                            result_tab[i][j] = app(result_tab[i][j],k)
                elif sk is not None:
                    for j in range(sk, sk+a):
                        result_tab[i][j] = rem(result_tab[i][j],k)
                        result_tab[i][j] = app(result_tab[i][j],k+1)
            else:
                if sk1 is not None:
                    if a > b:
                        result_tab[i][sk1-(a-b)-1] = app(result_tab[i][sk1-(a-b)-1],k+1)
                        for j in range(sk1-(a-b), sk1):
                            result_tab[i][j] = rem(result_tab[i][j],k)
                            result_tab[i][j] = app(result_tab[i][j],k+1)
                        result_tab[i][sk1] = rem(result_tab[i][sk1],k)
                    elif a < b:
                        result_tab[i][sk1-1] = rem(result_tab[i][sk1-1],k+1)
                        for j in range(sk1, sk1+b-a-1):
                            result_tab[i][j] = rem(result_tab[i][j],k+1)
                            result_tab[i][j] = app(result_tab[i][j],k)
                        result_tab[i][sk1+b-a-1] = app(result_tab[i][sk1+b-a-1],k)
                elif sk is not None:
                    result_tab[i][sk] = app(result_tab[i][sk],k+1)
                    for j in range(sk+1, sk+a):
                        result_tab[i][j] = rem(result_tab[i][j],k)
                        result_tab[i][j] = app(result_tab[i][j],k+1)
                    result_tab[i][sk+a] = rem(result_tab[i][sk+a],k)
        return SemistandardSetValuedTableau(result_tab)

    def pp(self):
        """
        Return the pretty print of ``self``.

        EXAMPLES::

            sage: T = SemistandardSetValuedTableau([[[1],[1],[8]],[[2],[6,2]],[[3,7,4]]])
            sage: T.pp()
            [   1   ][  1  ][ 8 ]
            [   2   ][ 2,6 ]
            [ 3,4,7 ]

            sage: U = [[[1,2,3,4,6],[6],[6,7],[8,9],[9,11,12],[12]], [[7],[7],[8,9,10],[10,11,13,14],[14]], [[8,9],[9,10],[11,13],[16,17,18]]]
            sage: SemistandardSetValuedTableau(U).pp()
            [ 1,2,3,4,6 ][  6   ][  6,7   ][     8,9     ][ 9,11,12 ][ 12 ]
            [     7     ][  7   ][ 8,9,10 ][ 10,11,13,14 ][   14    ]
            [    8,9    ][ 9,10 ][ 11,13  ][  16,17,18   ]
        """
        max_len = max(len(row) for row in self)
        str_len = [[sum([len(str(elt)) for elt in cell]) + len(cell) - 1 for cell in row] for row in self]
        col_max = [max(row[j] for row in str_len if j < len(row)) for j in range(max_len)]
        S = ""
        for row in self:
            for j in range(len(row)):
                s = ""
                for k in range(len(row[j])):
                    s += str(row[j][k])
                    if k < len(row[j])-1:
                        s+=","
                S+="[ "+'{st:{c}^{n}}'.format(st=s, c=" ", n=col_max[j])+" ]"
            S += "\n"
        print(S)

    def uncrowding(self):
        r"""
        Return the image of ``self`` under the uncrowding map.

        INPUT:

        - ``self`` -- semistandard set-valued tableau

        OUTPUT: a tuple of

        - ``P`` -- semistandard Young tableau

        - ``Q`` -- a flagged increasing tableau or skew tableau ``Q`` with
          same shape as ``P``

        EXAMPLES::
            sage: T = SemistandardSetValuedTableau([[[1], [1,2,3]], [[2,3]]])
            sage: T.uncrowding()
            ([[1, 1], [2, 2], [3, 3]], [['X', 'X'], ['X', 1], [1, 2]])

            sage: T = SemistandardSetValuedTableau([[[1], [1,2], [2]], [[2,3], [3,4,5]], [[4]]])
            sage: T.uncrowding()
            ([[1, 1, 2], [2, 2], [3, 3], [4, 4], [5]], [['X', 'X', 'X'], ['X', 'X'], ['X', 1], [2, 3], [3]])
        """
        P = SemistandardTableau([])
        Q = Tableau([])
        sequences = _insertion_sequence(self.to_list())
        for seq in sequences:
            P,Q = _uncrowding_insertion(seq, P, Q)
        return P,Q


class CrystalElementSemistandardSetValuedTableau(SemistandardSetValuedTableau):
    r"""
    Class for elements of ``crystals.SemistandardSetValuedTableaux``
    """
    def _get_signs(self, i):
        """
        Auxiliary function for `e_i` and `f_i` methods.

        Assign each column of ``self`` a +1, -1 or 0 according to
        +1 if there is an unmatched `i+1` aka left parenthesis '('
        -1 if there is an unmatched `i` aka right parenthesis ')'
        0 if all the `i`s and `i+1`s are matched
        
        Return list of +1, -1, 0 with length equal to number of columns of ``self``.
        """
        signs = []
        for col in self.conjugate():
            word = sum(col, ())
            if i in word and i+1 in word:
                signs += [0]
            elif i in word: # i in word, i+1 not in word
                signs += [-1]
            elif i+1 in word: # i not in word, i+1 is
                signs += [+1]
            else: # neither i nor i+1 in word
                signs += [0]
        return signs

    def _bracket(self, i, right=True):
        r"""
        Auxiliary function for `e_i` and `f_i` methods.

        Return index of column in self with rightmost `i` to be changed to `i+1` 
        or leftmost `i+1` to be changed to `i`.

        If right is True (default), then return index of rightmost `i`.
        If right is False, then return index of leftmost `i+1`.

        If no `i` can be changed to `i+1` or vice versa, return -1.
        """
        x = self._get_signs(i) # x is a list of +1, -1, 0
        if not right:
            x = [-j for j in x][::-1]
        count = 0
        index = -1
        for j in range(len(x)):
            if x[j] == -1:
                if count == 0:
                    index = j
                else:
                    count -= 1
            if x[j] == +1:
                count += 1
        if right:
            return index
        else:
            return -1 if index < 0 else len(x)-1-index

    def e(self, i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: SSVT = SemistandardSetValuedTableaux([2,1],max_entry=3)
            sage: T1 = SSVT([[[1, 2], [2]], [[3]]])
            sage: T1.e(1)
            [[[1], [1, 2]], [[3]]]

            sage: SSVT = SemistandardSetValuedTableaux([2,1],max_entry=3)
            sage: T2 = SSVT([[[1, 2], [3]], [[3]]])
            sage: T2.e(2)
            [[[1, 2], [2]], [[3]]]
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        col = self._bracket(i,right=False)
        if col == -1:
            return None
        import copy
        T = copy.deepcopy(self)
        column = [T(cell) for cell in T.cells() if cell[1]==col]
        row = min([ j for j in range(len(column)) if i+1 in column[j] ]) 
        # checks that there is a cell to the left and that the cell contains i and i+1
        if col>0 and all(x in T(row,col-1) for x in [i,i+1]):
            entry = list(T(row,col-1))
            entry.remove(i+1)
            T = T.add_entry((row,col-1),tuple(entry))
        else:
            entry = list(T(row,col))
            entry.remove(i+1)
            T = T.add_entry((row,col),tuple(entry))
        entry = sorted(list(T(row,col))+[i])
        T = T.add_entry((row,col),tuple(entry))
        return self.parent()(T)

    def f(self, i):
        r"""
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: SSVT = SemistandardSetValuedTableaux([2,1],max_entry=3)
            sage: T1 = SSVT([[[1, 2], [2]], [[3]]])
            sage: T1.f(2)
            [[[1, 2], [3]], [[3]]]

            sage: SSVT = SemistandardSetValuedTableaux([2,1],max_entry=3)
            sage: T2 = SSVT([[[1], [1, 2]], [[3]]])
            sage: T2.f(1)
            [[[1, 2], [2]], [[3]]]
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        col = self._bracket(i,right=True)
        if col == -1:
            return None
        import copy
        T = copy.deepcopy(self)
        column = [T(cell) for cell in T.cells() if cell[1]==col]
        row = min([ j for j in range(len(column)) if i in column[j] ])
        # checks that there is a cell to the right and that the cell contains i and i+1
        if col<len(T[row])-1 and all(x in T(row,col+1) for x in [i,i+1]):
            entry = list(T(row,col+1))
            entry.remove(i)
            T = T.add_entry((row,col+1),tuple(entry))
        else:
            entry = list(T(row,col))
            entry.remove(i)
            T = T.add_entry((row,col),tuple(entry))
        entry = sorted(list(T(row,col))+[i+1])
        T = T.add_entry((row,col),tuple(entry))
        return self.parent()(T)

    # def weight(self):
    #     r"""
    #     Return the weight for ``self``.

    #     EXAMPLES::
        
    #         sage: SSVT = SemistandardSetValuedTableaux([2,1,1],max_entry=5)
    #         sage: T = SSVT([[[1], [2, 3]], [[2]], [[3]]])
    #         sage: T.weight()
    #         (1, 2, 2, 0, 0)
    #     """
    #     WLR = self.parent().weight_lattice_realization()
    #     return sum(sum(WLR.monomial(elt-1) for elt in self(cell)) for cell in self.cells())

    def reading_word(self):
        r"""

        Return the reading word of ``self``.

        EXAMPLES::

            sage: SSVT = SemistandardSetValuedTableaux([5,3,1,1],max_entry=6)
            sage: T = SSVT([[[1,2],[2,3],[3],[3,4,5],[5,6]], [[3],[4,6],[6]], [[4,5]],[[6]]])
            sage: T.reading_word()
            [6, 5, 4, 3, 6, 6, 4, 2, 3, 3, 5, 6, 5, 4, 3, 2, 1]
        """
        R = []
        for s in _insertion_sequence(self):
            R += s
        return R

    # def to_tableau(self):
    #     SSVT = SemistandardSetValuedTableaux(self.size(),max_entry=self.parent().max_entry)
    #     return SSVT.element_class(SSVT,self)


class SemistandardSetValuedTableaux(Tableaux):
    """

    Return the class of semistandard set-valued tableaux.

    INPUT:
    
    Positional arguments:

    - ``p`` -- first argument is either a non-negative integer or a partition
    
    Keyword arguments:

    - ``max_entry`` -- positive integer that is the maximum allowed entry in the tableau

    OUTPUT:

    - With no shape or size, the class of all semistandard set-valued tableaux.

    - With no ``max_entry``, the class of all semistandard set-valued tableaux of 
      specified size or shape.

    - With a ``max_entry``, one of the following:
        - With a non-negative integer ``p``, the class of all semistandard 
          set-valued tableaux of size ``p`` and maximum integer ``max_entry``.

        - With a partition ``p``, the class of all semistandard 
          set-valued tableaux of shape ``p`` and maximum integer ``max_entry``.

    A semistandard set-valued tableau is a tableau whose cells are filled with 
    nonempty sets of integers such that the tableau is semistandard with respect 
    to the ordering on sets `A`, `B` given by `A <= B` if `max(A) < min(B)`.

    EXAMPLES::

        sage: SSVT = SemistandardSetValuedTableaux(3); SSVT
        Semistandard set-valued tableaux of size 3

        sage: SSVT = SemistandardSetValuedTableaux(3, max_entry=2); SSVT
        Semistandard set-valued tableaux of size 3 and max entry 2
        sage: list(SSVT)
        [[[[1], [1], [1]]],
         [[[1], [1], [1, 2]]],
         [[[1], [1], [2]]],
         [[[1], [1, 2], [2]]],
         [[[1], [2], [2]]],
         [[[1, 2], [2], [2]]],
         [[[2], [2], [2]]],
         [[[1], [1]], [[2]]],
         [[[1], [1, 2]], [[2]]],
         [[[1], [2]], [[2]]]]

        sage: SSVT = SemistandardSetValuedTableaux([2,2], max_entry=3); SSVT
        Semistandard set-valued tableaux of shape [2, 2] and max entry 3
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        """
        Normalize and process input to return the correct parent and
        ensure a unique representation.

        TESTS::

            sage: SemistandardSetValuedTableaux()
            Semistandard set-valued tableaux
            sage: SemistandardSetValuedTableaux(max_entry=4)
            Semistandard set-valued tableaux of max entry 4
            sage: SemistandardSetValuedTableaux(2)
            Semistandard set-valued tableaux of size 2
            sage: SemistandardSetValuedTableaux(2, max_entry=4)
            Semistandard set-valued tableaux of size 2 and max entry 4
            sage: SemistandardSetValuedTableaux([4,2,1])
            Semistandard set-valued tableaux of shape [4, 2, 1]
            sage: SemistandardSetValuedTableaux([4,2,1], max_entry=3)
            Semistandard set-valued tableaux of shape [4, 2, 1] and max entry 3
            sage: SemistandardSetValuedTableaux([1], max_entry=None)
            Semistandard set-valued tableaux of shape [1]
            sage: from sage.rings.infinity import PlusInfinity
            sage: SemistandardSetValuedTableaux([1], max_entry=PlusInfinity())
            Semistandard set-valued tableaux of shape [1]
            sage: SemistandardSetValuedTableaux([])
            Semistandard set-valued tableaux of shape []

            sage: SemistandardSetValuedTableaux(7, max_entry=0)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must be a positive integer, None, or PlusInfinity
            sage: SemistandardSetValuedTableaux([2], max_entry=0)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must be a positive integer, None, or PlusInfinity
            sage: SemistandardSetValuedTableaux({1,2,3})
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a non-negative integer or a partition
        """
        max_entry = kwargs.get('max_entry')
        if args:
            p = args[0]

            # if p is size
            if isinstance(p,(int,Integer)) and p >= 0:
                if max_entry is None or max_entry == PlusInfinity():
                    #return SemistandardSetValuedTableaux_size_inf(Integer(p), max_entry=None)
                    return SemistandardSetValuedTableaux_size(Integer(p), max_entry=None)
                elif isinstance(max_entry, (int,Integer)) == False or max_entry <= 0:
                    raise ValueError("the maximum entry must be a positive integer, None, or PlusInfinity")
                #elif isinstance(max_entry, (int,Integer)) and max_entry <= 0:
                #    raise ValueError("the maximum entry must be a positive integer, None, or PlusInfinity")
                return SemistandardSetValuedTableaux_size(Integer(p), max_entry)

            # if p is shape
            elif p in _Partitions:
                if max_entry is None or max_entry == PlusInfinity():
                    #return SemistandardSetValuedTableaux_shape_inf(_Partitions(p))
                    return SemistandardSetValuedTableaux_shape(_Partitions(p), max_entry=None)
                elif isinstance(max_entry, (int,Integer)) == False or max_entry <= 0:
                    raise ValueError("the maximum entry must be a positive integer, None, or PlusInfinity")
                return SemistandardSetValuedTableaux_shape(_Partitions(p), max_entry)
            else:
                raise ValueError("the argument must be a non-negative integer or a partition")            
        else:
            if max_entry is None or max_entry == PlusInfinity():
                return SemistandardSetValuedTableaux_all(max_entry=None)
            elif isinstance(max_entry, (int,Integer)) == False or max_entry <= 0:
                raise ValueError("the maximum entry must be a positive integer, None, or PlusInfinity")
            return SemistandardSetValuedTableaux_all(max_entry=max_entry)

    Element = SemistandardSetValuedTableau

    def __init__(self, *args, **kwargs):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: S = SemistandardSetValuedTableaux(2, max_entry=3)
            sage: TestSuite(S).run()
        """
        if 'max_entry' in kwargs:
            if kwargs['max_entry']==PlusInfinity():
                self.max_entry = None
            else:
                self.max_entry = kwargs['max_entry']
            kwargs.pop('max_entry')
        else:
            self.max_entry = None

        Tableaux.__init__(self, *args, **kwargs)

    def __contains__(self, t):
        """
        Determine if ``t`` is an element of ``self``.

        TESTS::

            sage: T1 = [[[1], [1,2]], [[2]]]
            sage: T2 = Tableau([[[1], [1,3]], [[2]]])
            sage: T3 = Tableau([[[1,2], [2]], [[2]]])
            sage: T1 in SemistandardSetValuedTableaux(3)
            True
            sage: T2 in SemistandardSetValuedTableaux(3)
            True
            sage: T3 in SemistandardSetValuedTableaux(3)
            False
            sage: T2 in SemistandardSetValuedTableaux([2,1])
            True
            sage: T1 in SemistandardSetValuedTableaux([2,1,1])
            False
            sage: T1 in SemistandardSetValuedTableaux(3, max_entry=2)
            True
            sage: T2 in SemistandardSetValuedTableaux(3, max_entry=2)
            False

            sage: SSVT = []
            sage: for st in StandardTableaux([3,1,1]):
            ....:     SSVT.append([[[_] for _ in row] for row in st])
            ....:        
            sage: all(i in SemistandardSetValuedTableaux([3,1,1]) for i in SSVT)
            True
        """
        if isinstance(t,SemistandardSetValuedTableau):
            return True
        # t is assumed to be at least a list of lists with shape given by a partition
        elif Tableaux.__contains__(self,t):         
            for row in t:
                for cell in row:
                    # checks that cell contains a set
                    if not isinstance(cell, (list,tuple,set)) or len(cell) > len(set(cell)):
                        return False
                    # checks that the set consists of nonempty set of integers
                    if len(cell)==0 or not all(isinstance(elt,(int,Integer)) and elt>0 for elt in cell):
                        return False
                # checks cells are weakly increasing along rows
                for i in range(len(row)-1):
                    left,right = row[i],row[i+1]
                    if max(left) > min(right):
                        return False
            # checks cells are strictly increasing along columns
            for up, down in zip(t, t[1:]):
                if any(max(up[i]) >= min(down[i]) for i in range(len(down))):
                    return False 
            # if self has max_entry, checks if that all elements in each cell are within max_entry
            return self.max_entry is None or max(max(max(cell) for cell in row) for row in t) <= self.max_entry


class SemistandardSetValuedTableaux_all(SemistandardSetValuedTableaux):
    """
    Class of all semistandard set-valued tableaux.
    """
    def __init__(self, max_entry=None):
        """
        Initialize the class of all semistandard set-valued tableaux.

        TESTS::

            sage: SSVT = SemistandardSetValuedTableaux()
            sage: [[[1,2], [2]], [[3]]] in SSVT
            True
            sage: [[[1,2], [2]], [[2]]] in SSVT
            False
            sage: [[[1,2], [1,3]], [[3]]] in SSVT
            False
            sage: [[[1,2], []], [[3]]] in SSVT
            False
            sage: TestSuite(SSVT).run()  # long time
        """
        if max_entry == PlusInfinity() or max_entry is None:
            SemistandardSetValuedTableaux.__init__(self, category=InfiniteEnumeratedSets())
        else:
            SemistandardSetValuedTableaux.__init__(self, max_entry=max_entry, category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: SemistandardSetValuedTableaux()
            Semistandard set-valued tableaux
            sage: SemistandardSetValuedTableaux(max_entry=4)
            Semistandard set-valued tableaux of max entry 4
        """
        if self.max_entry is None:
            return "Semistandard set-valued tableaux"
        return "Semistandard set-valued tableaux of max entry {}".format(self.max_entry)

    def __contains__(self, x):
        """
        Determine if ``t`` is an element of ``self``.

        TESTS::     

            sage: SSVT1 = SemistandardSetValuedTableaux()
            sage: [[[1,2]],[[2]],[[4]]] in SSVT1
            False
            sage: [[[1,2],[1,3]],[[2]],[[4]]] in SSVT1
            False
            sage: [[[1,2]],[[3]],[[4]]] in SSVT1
            True
            sage: [[[1,2],[2]],[[3]],[[4]]] in SSVT1
            True
            sage: [[[1,2],[3],[40]]] in SSVT1
            True
            sage: [[[1]],[[3,40]],[[50]]] in SSVT1
            True

            sage: SSVT2 = SemistandardSetValuedTableaux(max_entry=4)
            sage: [[[1,2]],[[2]],[[4]]] in SSVT2
            False
            sage: [[[1,2],[1,3]],[[2]],[[4]]] in SSVT2
            False
            sage: [[[1,2]],[[3]],[[4]]] in SSVT2
            True
            sage: [[[1,2],[2]],[[3]],[[4]]] in SSVT2
            True
            sage: [[[1,2],[3],[4]]] in SSVT2
            True
            sage: [[[1]],[[3,4]],[[5]]] in SSVT2
            False
        """
        return SemistandardSetValuedTableaux.__contains__(self, x)


class SemistandardSetValuedTableaux_size(SemistandardSetValuedTableaux, DisjointUnionEnumeratedSets):
    """
    Class of all semistandard set-valued tableaux of a fixed size.
    """
    def __init__(self, n, max_entry=None):
        """
        Initializes a class of semistandard set-valued tableaux of size ``n``.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardSetValuedTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: TestSuite(SemistandardSetValuedTableaux(0,max_entry=4)).run()
            sage: TestSuite(SemistandardSetValuedTableaux(3,max_entry=2)).run()
        """
        #if max_entry is None:
        #    self.max_entry = n
        if max_entry is None:
            category = InfiniteEnumeratedSets()
            Parts = Partitions_n(n)
        else:
            category = FiniteEnumeratedSets()
            Parts = [P for P in Partitions_n(n) if len(P)<=max_entry]

        SSVT = lambda p: SemistandardSetValuedTableaux_shape(p, max_entry=max_entry)
        DisjointUnionEnumeratedSets.__init__(self, Family(Parts, SSVT),
                                            facade=True, keepkey=False)
        SemistandardSetValuedTableaux.__init__(self,max_entry=max_entry,category=category)
        self._size = Integer(n)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: SemistandardSetValuedTableaux(3,max_entry=4)
            Semistandard set-valued tableaux of size 3 and max entry 4
            sage: SemistandardSetValuedTableaux(5)
            Semistandard set-valued tableaux of size 5
        """
        if self.max_entry is None:
             return "Semistandard set-valued tableaux of size {}".format(self._size)
        return "Semistandard set-valued tableaux of size {} and max entry {}".format(self._size, self.max_entry)

    def __contains__(self, x):
        """
        Determine if ``x`` is an element of ``self``.

        TESTS::

            sage: SSVT = SemistandardSetValuedTableaux(3,max_entry=4)
            sage: [[[1,2]],[[2]],[[4]]] in SSVT
            False
            sage: [[[1,2],[1,3]],[[2]],[[4]]] in SSVT
            False
            sage: [[[1,2]],[[3]],[[4]]] in SSVT
            True
            sage: [[[1,2],[2]],[[3]],[[4]]] in SSVT
            False
            sage: [[[1,2],[3],[4]]] in SSVT
            True
            sage: [[[1]],[[3,4]],[[5]]] in SSVT
            False

            sage: SSVT = SemistandardSetValuedTableaux(3)
            sage: [[[1,2]],[[2]],[[3]]] in SSVT
            False
            sage: [[[1,2]],[[3]],[[4]]] in SSVT
            True
            sage: [[[1,2],[2]],[[3]],[[4]]] in SSVT
            False
            sage: [[[1,2],[2],[3]]] in SSVT
            True
            sage: [[[1]],[[3,4]],[[50]]] in SSVT
            True
        """
        return SemistandardSetValuedTableaux.__contains__(self, x) and \
        sum(map(len, x)) == self._size

    def __iter__(self):
        """
        Iterate over all semistandard set-valued tableaux associated to the
        size of ``self``.

        EXAMPLES::

            sage: SSVT = SemistandardSetValuedTableaux(3, max_entry=2)
            sage: SSVT.list()
            [[[[1], [1], [1]]],
            [[[1], [1], [1, 2]]],
            [[[1], [1], [2]]],
            [[[1], [1, 2], [2]]],
            [[[1], [2], [2]]],
            [[[1, 2], [2], [2]]],
            [[[2], [2], [2]]],
            [[[1], [1]], [[2]]],
            [[[1], [1, 2]], [[2]]],
            [[[1], [2]], [[2]]]]
            sage: len(SSVT)
            10
        """
        for part in Partitions(self._size):
            for ssvt in SemistandardSetValuedTableaux_shape(part, self.max_entry):
                yield self.element_class(self, ssvt)


class SemistandardSetValuedTableaux_shape(SemistandardSetValuedTableaux):
    """
    Class of all semistandard set-valued tableaux of a fixed shape.
    """
    @staticmethod
    def __classcall_private__(cls, p, max_entry=None):
        """
        Normalize the attributes for the class.

        TESTS::

            sage: SSVT = SemistandardSetValuedTableaux([2,1])
            sage: SSVT._shape.parent()
            Partitions

            sage: SSVT1 = SemistandardSetValuedTableaux(shape=(2,1), max_entry=3)
            sage: SSVT2 = SemistandardSetValuedTableaux(shape=[2,1], max_entry=3)
            sage: SSVT1 is SSVT2
            True
        """
        shape = _Partitions(p)
        return SemistandardSetValuedTableaux.__classcall__(cls, shape, \
               max_entry=max_entry)

    def __init__(self, p, max_entry):
        """
        Initialize a class of semistandard set-valued tableaux of given shape ``p``.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardSetValuedTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: TestSuite( SemistandardSetValuedTableaux([2,1,1],max_entry=3) ).run()
        """
        #if max_entry is None:
        #    self.max_entry = len(p)
        if max_entry is None:
            category = RegularCrystals().Infinite()
            self._cartan_type = CartanType(['A+oo'])
        else:
            category = ClassicalCrystals()
            self._cartan_type = CartanType(['A', max_entry-1])
        self._shape = p
        self.Element = CrystalElementSemistandardSetValuedTableau
        SemistandardSetValuedTableaux.__init__(self, max_entry=max_entry)
        Parent.__init__(self, category=category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: SemistandardSetValuedTableaux([3,1], max_entry=2)
            Semistandard set-valued tableaux of shape [3, 1] and max entry 2
            sage: SemistandardSetValuedTableaux([3,1])
            Semistandard set-valued tableaux of shape [3, 1]
        """
        if self.max_entry is None:
            return "Semistandard set-valued tableaux of shape {}".format(self._shape)
        return "Semistandard set-valued tableaux of shape {} and max entry {}".format(self._shape, self.max_entry)

    def __contains__(self, x):
        """
        Determine if ``x`` is an element of ``self``.

        EXAMPLES::
        
            sage: SSVT = SemistandardSetValuedTableaux([2,1,1],max_entry=6)
            sage: [[[1,2],[2,3]],[[2]],[[3]]] in SSVT
            False
            sage: [[[1,2],[2,3]],[[3,4]],[[6]]] in SSVT
            True
            sage: [[[1,2],[2,3]],[[3,6]],[[6,7]]] in SSVT
            False
            sage: [[[1,2],[1,3]],[[3,4]],[[6]]] in SSVT
            False
            sage: [[[1,2],[2,3]],[[3,4],[4,5,6]],[[6]]] in SSVT
            False

            sage: SSVT = SemistandardSetValuedTableaux([2,1,1])
            sage: [[[1,2],[2,3]],[[2]],[[3]]] in SSVT
            False
            sage: [[[1,2],[2,3]],[[3,4]],[[6]]] in SSVT
            True
            sage: [[[1],[2,3,4]],[[3,10]],[[100]]] in SSVT
            True
            sage: [[[1,2],[1,3]],[[3,4]],[[6]]] in SSVT
            False
            sage: [[[1,2],[2,3]],[[3,4],[4,5,6]],[[6]]] in SSVT
            False         
        """
        return SemistandardSetValuedTableaux.__contains__(self, x) and \
               [len(row) for row in x] == self._shape

    @lazy_attribute
    def module_generators(self):
        """
        Return the generators of ``self`` as a crystal.

        .. WARNING::

            If ``max_entry`` is None, this will raise a ValueError.

        TESTS::

            sage: SSVT = SemistandardSetValuedTableaux([2,2],max_entry=3)
            sage: SSVT.module_generators
            ([[[1], [1]], [[2], [2]]], [[[1], [1]], [[2], [2, 3]]], [[[1], [1, 2]], [[2], [3]]], [[[1], [1, 2]], [[2, 3], [3]]])
        """
        P = self._shape
        # if len(P) == 0:
        #     return tuple([SemistandardSetValuedTableau([])])
        m = self.max_entry
        if m is None:
            # L = _generate_pairs(P, m)
            # crowd = lambda i: L[i][0].crowding(L[i][1], mark=None)
            # return Family(L, crowd, lazy=True)
            raise ValueError("max_entry must be specified to produce module_generators")
        else:
            L = _generate_pairs(P, m)
            return tuple([S.crowding(F, mark=None, m=m) for S,F in L])

    def shape(self):
        """
        Return the shape of the semistandard set valued tableau ``self``.

        TESTS::

            sage: SemistandardSetValuedTableaux([3,1],max_entry=2).shape()
            [3, 1]
            sage: SemistandardSetValuedTableaux([2,2]).shape()
            [2, 2]
        """
        return self._shape

    def __iter__(self):
        """
        Iterate over all semistandard set-valued tableaux associated to the
        shape of ``self``.

        .. WARNING::

            If ``max_entry`` is None, this will raise a ValueError.

        EXAMPLES::

            sage: SSVT = SemistandardSetValuedTableaux([2,2],max_entry=3)
            sage: SSVT.list()
            [[[[1], [1]], [[2], [2]]],
            [[[1], [1]], [[2], [2, 3]]],
            [[[1], [1]], [[2], [3]]],
            [[[1], [1]], [[2, 3], [3]]],
            [[[1], [1, 2]], [[2], [3]]],
            [[[1], [1, 2]], [[2, 3], [3]]],
            [[[1], [1]], [[3], [3]]],
            [[[1], [1, 2]], [[3], [3]]],
            [[[1], [2]], [[2], [3]]],
            [[[1], [2]], [[2, 3], [3]]],
            [[[1], [2]], [[3], [3]]],
            [[[1, 2], [2]], [[3], [3]]],
            [[[2], [2]], [[3], [3]]]]
            sage: len(SSVT)
            13
        """
        if self.max_entry is None:
            raise ValueError("max_entry must be specified to iterate")

        def jumps(row):
            """
            Return list of indices i where row[i] < row[i+1] and also
            the last index of row.

            INPUT:
            
            - ``row`` -- a list of list of integers
            """
            return [i for i in range(len(row)-1) if max(row[i])<min(row[i+1])] + [len(row)-1]
        
        def addable(cell, right, below, max_entry):
            """
            Return a list of numbers that can added to the cell ``cell`` of a
            semistandard set-valued tableau.

            INPUT:

            - ``cell`` -- nonempty list of integers

            - ``right`` -- list of integers or None

            - ``below`` -- list of integers or None

            - ``max_entry`` -- nonnegative integer
            """
            cell_val = max(cell)+1
            right_val = min(right)+1 if right is not None else max_entry+1
            below_val = min(below) if below is not None else max_entry+1

            if right_val < below_val:
                return list(range(cell_val, right_val))
            else:
                return list(range(cell_val, below_val))

        from sage.misc.mrange import cartesian_product_iterator as CPI
        import copy
        for t in SemistandardTableaux(shape=self._shape, max_entry=self.max_entry):
            tab = [[[entry] for entry in row] for row in t]
            addable_indices = []
            addable_sets = []
            for i in range(len(tab)):
                J = jumps(tab[i])
                for j in J:
                    cell = tab[i][j]
                    right, below = None, None
                    if i+1 < len(tab):
                        if j < len(tab[i+1]):
                            below = tab[i+1][j]
                    if j+1 < len(tab[i]):
                        right = tab[i][j+1]

                    # if a cell can be added with extra letters
                    addable_indices += [(i,j)]
                    elts_to_add = [tuple(_) for _ in powerset(addable(cell,right,below,self.max_entry))]
                    addable_sets.append(elts_to_add)
            
            for cp in CPI(addable_sets):
                tab_copy = copy.deepcopy(tab)
                for k in range(len(addable_sets)):
                    i,j = addable_indices[k]
                    tab_copy[i][j] += cp[k]
                yield self.element_class(self, tab_copy)

######################
#  Helper functions  #
######################

def _insertion_sequence(T):
    """
    Returns a sequence of words to insert in the uncrowding map.

    The algorithm assumes that entries in cells of T are sorted in increasing order.
    
    EXAMPLES::

        sage: T = SemistandardSetValuedTableau([[[2,3,1], [6,3,5], [7], [11,8,12,10]], [[4,5,6,7], [7], [9,12,8]], [[8], [9,8], [13]], [[9,10]]])
        sage: from sage.combinat.semistandard_set_valued_tableau import _insertion_sequence
        sage: _insertion_sequence(T)
        [[10, 9], [8, 9, 13, 8], [7, 7, 12, 9, 8, 6, 5, 4], [3, 6, 7, 12, 11, 10, 8, 5, 3, 2, 1]]
    """
    if T not in SemistandardSetValuedTableaux():
        raise ValueError("Semistandard set-valued tableau not given")

    if len(T) == 0:
        return []
    
    seq = []
    for row in T[::-1]:
        S = []
        M = [cell[-1] for cell in row[:-1]]
        S += M
        H = row[-1][::-1]
        for cell in row[:-1][::-1]:
            H += cell[:-1][::-1]
        S += H
        seq += [S]
    return seq

def _uncrowding_insertion(seq, P=None, Q=None, mark='X'):
    """
    Return the image under the uncrowding insertion for ``seq`` and pair of
    tableaux ``P`` and ``Q``.

    INPUT:
    - ``seq`` -- a sequence of integers to be inserted. This should insert to
                 a hook shape.
        
    - ``P`` -- a semistandard Young tableau (default is empty)
        
    - ``Q`` -- a flagged increasing tableau with same shape as ``P`` 
      (default is empty)

    - ``mark`` -- symbol used to create a flag on tableau ``Q`` (default: 'X'). 
      If skew tableau is specified, ``mark`` should be assigned as None.

    OUTPUT: a tuple of
    
    - semistandard Young tableau `P'`
    
    - flagged increasing tableau `Q'` with same shape as `P'`

    EXAMPLES::

        sage: P = SemistandardTableau([[4],[6],[7]])
        sage: Q = Tableau([['X'],['X'],[1]])
        sage: from sage.combinat.semistandard_set_valued_tableau import _uncrowding_insertion
        sage: _uncrowding_insertion([3,3,9,8],P,Q)
        ([[3, 3, 8], [4, 9], [6], [7]], [['X', 'X', 'X'], ['X', 1], ['X'], [1]])

        sage: from sage.combinat.semistandard_set_valued_tableau import _uncrowding_insertion
        sage: _uncrowding_insertion([1,1,4,5,4,3,2])
        ([[1, 1, 2, 4], [3], [4], [5]], [['X', 'X', 'X', 'X'], [1], [2], [3]])
    """
    Pp = P
    Qq = Q
    if P is None:
        Pp = SemistandardTableau([])
    if Q is None:
        Qq = Tableau([])
    if not isinstance(Pp, SemistandardTableau):
        raise ValueError("P should be instance of SemistandardTableau")
    if not isinstance(Qq, Tableau):
        raise ValueError("Q should be instance of Tableau")
    if Pp.shape() != Qq.shape():
        raise ValueError("P and Q must be of same shape")        

    # Insertion is done on P first
    for x in seq:
        Pp = Pp.bump(x)
    shape = Pp.shape()

    # Recording is then done on Q
    temp = Tableau([[mark]*shape[0]] + Qq.to_list())
    cells = [pair for pair in Pp.cells() if pair not in temp.cells()]
    for cell in cells: 
        temp = temp.add_entry(cell, cell[0])
    Qq = temp
    return Pp, Qq

def _max_outer_shape(P, m):
    """
    Return the largest outer shape of flagged increasing tableau given 
    partition ``P`` and maximum entry ``m``.

    EXAMPLES::

        sage: P,m = [8,6,3,1], 6
        sage: from sage.combinat.semistandard_set_valued_tableau import _max_outer_shape
        sage: _max_outer_shape(P, m)
        [8, 7, 5, 4, 4, 4]

        sage: from sage.combinat.semistandard_set_valued_tableau import _max_outer_shape
        sage: P,m = [23,20,18,18,12,7,4,2,1],10
        sage: _max_outer_shape(P, m)
        [23, 21, 20, 20, 16, 12, 10, 9, 9, 9]
    """
    if P not in _Partitions or len(P) == 0:
        raise ValueError("P should be a nonempty integer partition")
    if len(P) > m:
        raise ValueError("m should be at the least length of partition P")
    Q = P + [0]*(m-len(P))
    L = [P[0]]
    for i in range(1, m):
        L += [min(L[-1], i+Q[i])]
    return Partition(L)

def _is_flagged_increasing(T):
    """
    Return True if `T` is a flagged increasing tableau and False otherwise.

    `T` needs to be a semistandard skew tableau.

    EXAMPLES::
        sage: T1 = SkewTableau([[None,None],[None,2],[1,1]])
        sage: from sage.combinat.semistandard_set_valued_tableau import _is_flagged_increasing
        sage: _is_flagged_increasing(T1)
        Traceback (most recent call last):
        ...
        ValueError: T needs to be a semistandard skew tableau

        sage: from sage.combinat.semistandard_set_valued_tableau import _is_flagged_increasing
        sage: T2 = SkewTableau([[None,None],[None,2],[1,3]])
        sage: _is_flagged_increasing(T2)
        False

        sage: from sage.combinat.semistandard_set_valued_tableau import _is_flagged_increasing
        sage: T3 = SkewTableau([[None,None],[None,1],[1,2]])
        sage: _is_flagged_increasing(T3)
        True

        sage: from sage.combinat.semistandard_set_valued_tableau import _is_flagged_increasing
        sage: T4 = SkewTableau([[None,None],[None,1],[2,2]])
        sage: _is_flagged_increasing(T4)
        False
    """
    if isinstance(T,SkewTableau) and T.is_semistandard():
        # Checks for flagged condition
        for i in range(1,len(T)):
            values = [T[x][y] for x,y in T.cells() if T[x][y] in T[i]]
            if len(values) > 0 and any([k > i for k in values]):
                return False
        return T.conjugate().is_semistandard()
    raise ValueError("T needs to be a semistandard skew tableau")

def _highest_weight_tableau(P):
    """
    Return the highest weight semistandard tableau of shape `P`.

        EXAMPLES::
            sage: P = [5,3,2,1]
            sage: from sage.combinat.semistandard_set_valued_tableau import _highest_weight_tableau
            sage: _highest_weight_tableau(P)
            [[1, 1, 1, 1, 1], [2, 2, 2], [3, 3], [4]]
    """
    return SemistandardTableau([[i+1]*P[i] for i in range(len(P))])

def _generate_pairs(P, m):
    r"""
    Return the set of all pairs of valid tableaux given partition and maximum entry.
    
    INPUT:

    - ``P`` -- partition
    
    - ``m`` -- maximum entry

    OUTPUT: List of all pairs `(S,F)` where
    
    - `S` is a highest weight semistandard tableau
    
    - `F` is a flagged increasing tableau
    """
    # Due to some quirk with SemistandardSkewTableaux, one needs separate 
    # initialization for a flagged increasing tableau with same inner and 
    # outer shape

    if m is None:
        # def pairs(k):
        #     if k ==len(P):
        #         S,F = _highest_weight_tableau(P),Tableau([[None]*P[i] for i in range(len(P))])
        #         return tuple([(S,F)])
        #     out = _max_outer_shape(P, k)
        #     maximum = lambda T: max([T[x][y] for x,y in T.cells()])
        #     for n in range(sum(P)+1, len(P)*k+1):
        #         for Q in Partitions(n,inner=P,outer=out):
        #             # Filters all flagged increasing tableau of shape Q/P with largest entry exactly k
        #             Sk = [F for F in SemistandardSkewTableaux([Q,P], max_entry=k) if _is_flagged_increasing(F) and maximum(F)==k]
        #             L += [(_highest_weight_tableau(Q),Tableau(F)) for F in Sk]
        #     return tuple(L)

        # NN = NonNegativeIntegers()
        # Int = Family(NN, lambda i:i+len(P), lazy=True)
        # F0 = Family(Int, pairs, lazy=True)
        # def straighten(F0, k):
        #     floor, i = 0, 0
        #     while k-floor >= len(F0[i]):
        #         floor += len(F0[i])
        #         i += 1
        #     return F0[i][k-floor]

        # L = Family(F0, straighten, lazy=True)
        # return L
        raise ValueError("m has to be a positive integer") 
        # currently, having no maximum entry is not supported

    else:
        if len(P) == 0:
            return tuple([(SemistandardTableau([]), Tableau([]))])

        S = _highest_weight_tableau(P)
        F = Tableau([[None]*P[i] for i in range(len(P))])
        L = [(S, F)]
        out = _max_outer_shape(P, m)
        for n in range(sum(P)+1, len(P)*m+1):
            for Q in Partitions(n, inner=P, outer=out):
                Sk = [F for F in SemistandardSkewTableaux([Q,P], max_entry=m) \
                            if _is_flagged_increasing(F)]
                L += [(_highest_weight_tableau(Q),Tableau(F)) for F in Sk]
        return tuple(L)