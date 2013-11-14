r"""
Immaculate Tableaux

AUTHORS:

- Travis Scrimshaw (11-13-2013): Initial implementation.

This file consists of the following major classes:

Element classes:

* :class:`ImmaculateTableau`
* :class:`StandardImmaculateTableau`

Factory class:

* :class:`ImmaculateTableaux`
* :class:`StandardImmaculateTableaux`

Parent classes:

* :class:`ImmaculateTableaux_all`
* :class:`ImmaculateTableaux_size`
* :class:`ImmaculateTableaux_size_weight`
* :class:`ImmaculateTableaux_shape`
* :class:`ImmaculateTableaux_shape_weight`
* :class:`StandardImmaculateTableaux_all`
* :class:`StandardImmaculateTableaux_size`
* :class:`StandardImmaculateTableaux_shape`
"""
#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
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
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.element import Element
from sage.structure.global_options import GlobalOptions
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.decorators import rename_keyword
from sage.rings.infinity import PlusInfinity
from sage.rings.arith import factorial
from sage.rings.integer import Integer
from sage.combinat.combinat import CombinatorialObject
from sage.combinat.composition import Composition, Compositions
from integer_vector import IntegerVectors
import sage.libs.symmetrica.all as symmetrica
import sage.misc.prandom as random
import copy
import permutation
from sage.misc.flatten import flatten
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.misc.misc import uniq, prod
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_cat import Sets
from sage.combinat.combinatorial_map import combinatorial_map
from sage.misc.superseded import deprecated_function_alias

class ImmaculateTableau(CombinatorialObject, Element):
    """
    An immaculate tableau.

    INPUT:

    - ``t`` -- an immaculate tableau, a list of lists, or an empty list

    OUTPUT:

    An immaculate tableau object constructed from ``t``.
    """
    def __init__(self, parent, t):
        r"""
        Initializes a tableau.

        TESTS::

            sage: t = Tableaux()([[1,1],[1]])
            sage: s = Tableaux(3)([[1,1],[1]])
            sage: s==t
            True
            sage: t.parent()
            Tableaux
            sage: s.parent()
            Tableaux of size 3
            sage: r = Tableaux()(s); r.parent()
            Tableaux
            sage: s is t # identical tableaux are distinct objects
            False
        """
        if isinstance(t, Tableau):
            Element.__init__(self, parent)
            # Since we are (suppose to be) immutable, we can share the underlying data
            CombinatorialObject.__init__(self, t._list)
            return

        # CombinatorialObject verifies that t is a list
        # We must verify t is a list of lists
        if not all(isinstance(row, list) for row in t):
            raise ValueError("A tableau must be a list of lists.")

        Element.__init__(self, parent)
        CombinatorialObject.__init__(self, t)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: Tableaux.global_options(display="list")
            sage: t
            [[1, 2, 3], [4, 5]]
            sage: Tableaux.global_options(display="array")
            sage: t
              1  2  3
              4  5
            sage: Tableaux.global_options(display="compact"); t
            1,2,3/4,5
            sage: Tableaux.global_options.reset()
        """
        return self._repr_list()
        #return self.parent().global_options.dispatch(self,'_repr_','display')

    def _repr_list(self):
        """
        Return a string representation of ``self`` as a list.

        EXAMPLES::

            sage: T = Tableau([[1,2,3],[4,5]])
            sage: T._repr_list()
            '[[1, 2, 3], [4, 5]]'
        """
        return repr(self._list)

    def _repr_diagram(self):
        """
        Return a string representation of ``self`` as an array.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: print t._repr_diagram()
              1  2  3
              4  5
            sage: Tableaux.global_options(convention="french")
            sage: print t._repr_diagram()
              4  5
              1  2  3
            sage: Tableaux.global_options.reset()
        """
        return '\n'.join(["".join(map(lambda x: "%3s"%str(x) , row)) for row in self])
        #if self.parent().global_options('convention') == "English":
        #    return '\n'.join(["".join(map(lambda x: "%3s"%str(x) , row)) for row in self])
        #else:
        #    return '\n'.join(["".join(map(lambda x: "%3s"%str(x) , row)) for row in reversed(self)])

    def _ascii_art_(self):
        """
        TESTS::

            sage: ascii_art(list(StandardTableaux(3)))
            [                              1 ]
            [              1  3    1  2    2 ]
            [   1  2  3,   2   ,   3   ,   3 ]
            sage: Tableaux.global_options(ascii_art="compact")
            sage: ascii_art(list(StandardTableaux(3)))
            [                        |3| ]
            [          |2|    |3|    |2| ]
            [ |1|2|3|, |1|3|, |1|2|, |1| ]
            sage: Tableaux.global_options(ascii_art="table")
            sage: ascii_art(list(StandardTableaux(3)))
            [                                      +---+ ]
            [                                      | 3 | ]
            [                +---+      +---+      +---+ ]
            [                | 2 |      | 3 |      | 2 | ]
            [ +---+---+---+  +---+---+  +---+---+  +---+ ]
            [ | 1 | 2 | 3 |  | 1 | 3 |  | 1 | 2 |  | 1 | ]
            [ +---+---+---+, +---+---+, +---+---+, +---+ ]
            sage: Tableaux.global_options(convention="french")
            sage: Tableaux.global_options(ascii_art="repr")
            sage: ascii_art(list(StandardTableaux(3)))
            [                              3 ]
            [              2       3       2 ]
            [   1  2  3,   1  3,   1  2,   1 ]
            sage: Tableaux.global_options.reset()
        """
        #ascii = self.parent().global_options.dispatch(self,'_ascii_art_','ascii_art')
        ascii = self._ascii_art_compact()
        from sage.misc.ascii_art import AsciiArt
        return AsciiArt(ascii.splitlines())

    _ascii_art_repr = _repr_diagram

    def _ascii_art_table(self):
        """
        TESTS::

            sage: t = Tableau([[1,2,3],[4,5]]); print t._ascii_art_table()
            +---+---+
            | 4 | 5 |
            +---+---+---+
            | 1 | 2 | 3 |
            +---+---+---+
            sage: t = Tableau([]); print t._ascii_art_table()
            ++
            ++
        """
        if len(self) == 0: return "++\n++"
        matr = ""
        for row in self:
            l1 = ""; l2 =  ""
            for e in row:
                l1 += "+---"
                l2 += "| " + str(e) + " "
            l1 += "+"; l2 += "|"
            matr = l1 + "\n" + l2 + "\n" + matr
        matr += "+---"** Integer(len(self[0])) + "+"
        return matr

    def _ascii_art_compact(self):
        """
        TESTS::

            sage: t = Tableau([[1,2,3],[4,5]]); print t._ascii_art_compact()
            |4|5|
            |1|2|3|
        """
        if len(self) == 0:
            return "."
        matr = ''
        for row in self:
            l1 = ""
            for e in row:
                l1 += "|" + str(e)
            l1 += "|"
            matr = l1 + "\n" + matr
        return matr

    def _latex_(self):
        r"""
        Return a LaTeX version of ``self``.

        EXAMPLES::

            sage: t = Tableau([[1,1,2],[2,3],[3]])
            sage: latex(t)    # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\cline{1-3}
            \lr{1}&\lr{1}&\lr{2}\\\cline{1-3}
            \lr{2}&\lr{3}\\\cline{1-2}
            \lr{3}\\\cline{1-1}
            \end{array}$}
            }
            sage: Tableaux.global_options(convention="french")
            sage: latex(t)    # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[t]{*{3}c}\cline{1-1}
            \lr{3}\\\cline{1-2}
            \lr{2}&\lr{3}\\\cline{1-3}
            \lr{1}&\lr{1}&\lr{2}\\\cline{1-3}
            \end{array}$}
            }
            sage: Tableaux.global_options.reset()
        """
        return self._latex_diagram()
        #return self.parent().global_options.dispatch(self,'_latex_', 'latex')

    _latex_list = _repr_list

    def _latex_diagram(self):
        r"""
        Return a LaTeX representation of ``self`` as a Young diagram.

        EXAMPLES::

            sage: t = Tableau([[1,1,2],[2,3],[3]])
            sage: print t._latex_diagram()
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\cline{1-3}
            \lr{1}&\lr{1}&\lr{2}\\\cline{1-3}
            \lr{2}&\lr{3}\\\cline{1-2}
            \lr{3}\\\cline{1-1}
            \end{array}$}
            }
        """
        if len(self) == 0:
            return "{\\emptyset}"
        from output import tex_from_array
        return tex_from_array(self)

    def __call__(self, *cell):
        r"""

        INPUT:

        - ``cell`` -- a pair of integers, tuple, or list specifying a cell in
          the tableau

        OUTPUT:

        - The value in the corresponding cell.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: t(1,0)
            4
            sage: t((1,0))
            4
            sage: t(3,3)
            Traceback (most recent call last):
            ...
            IndexError: The cell (3,3) is not contained in [[1, 2, 3], [4, 5]]
        """
        try:
            i,j = cell
        except ValueError:
            i,j = cell[0]

        try:
            return self[i][j]
        except IndexError:
            raise IndexError, "The cell (%d,%d) is not contained in %s"%(i,j,self)

    @combinatorial_map(name='shape')
    def shape(self):
        r"""
        Return the shape of ``self``.

        EXAMPLES::

            sage: Tableau([[1,2,3],[4,5],[6]]).shape()
            [3, 2, 1]
        """
        from sage.combinat.composition import Composition
        return Composition([len(row) for row in self])

    def size(self):
        """
        Return the size of the shape of ``self``.

        EXAMPLES::

            sage: Tableau([[1, 4, 6], [2, 5], [3]]).size()
            6
            sage: Tableau([[1, 3], [2, 4]]).size()
            4
        """
        return sum([len(row) for row in self])

    def pp(self):
        """
        Returns a pretty print string of the tableau.

        EXAMPLES::

            sage: T = Tableau([[1,2,3],[3,4],[5]])
            sage: T.pp()
              1  2  3
              3  4
              5
            sage: Tableaux.global_options(convention="french")
            sage: T.pp()
              5
              3  4
              1  2  3
            sage: Tableaux.global_options.reset()
        """
        print self._repr_diagram()

    def to_word_by_row(self):
        """
        Return the word obtained from a row reading of the tableau ``self``
        (starting with the lowermost row, reading every row from left
        to right).

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).to_word_by_row()
            word: 3412
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word_by_row()
            word: 325146
        """
        from sage.combinat.words.word import Word
        w = []
        for row in reversed(self):
            w += row
        return Word(w)

    to_word = to_word_by_row

    def standardization(self, check=True):
        r"""
        Return the standardization of ``self``, assuming ``self`` is a
        semistandard immaculate tableau.
        """
        if check and self not in SemistandardTableaux():
            raise ValueError("the tableau must be semistandard")
        T = from_shape_and_word(self.shape(), self.to_word_by_row().standard_permutation())
        return StandardTableaux()(T)

    @combinatorial_map(name ='reading word permutation')
    def reading_word_permutation(self):
        """
        Returns a permutation with the entries of ``self`` obtained by reading
        ``self`` in the given reading order.

        EXAMPLES::

            sage: StandardTableau([[1,2],[3,4]]).reading_word_permutation()
            [3, 4, 1, 2]

        Check that :trac:`14724` is fixed::

            sage: SemistandardTableau([[1,1]]).reading_word_permutation()
            [1, 2]
        """
        return permutation.Permutation(self.standardization().to_word())

    def entries(self):
        """
        Returns a list of all entries of self, in the order obtained
        by reading across the rows.

        EXAMPLES::

            sage: t = Tableau([[1,3], [2]])
            sage: t.entries()
            [1, 3, 2]
        """
        return sum(self, [])

    def entry(self, cell):
        """
        Returns the entry of cell in self. Cell is a tuple (i,j) of
        coordinates.

        EXAMPLES::

            sage: t = Tableau([[1,2],[3,4]])
            sage: t.entry( (0,0) )
            1
            sage: t.entry( (1,1) )
            4
        """
        i,j = cell
        return self[i][j]

    def weight(self):
        r"""
        Return the weight of the tableau ``self``. Trailing zeroes are
        omitted when returning the weight.

        The weight of a tableau `T` is the sequence `(a_1, a_2, a_3, \ldots )`,
        where `a_k` is the number of entries of `T` equal to `k`. This
        sequence contains only finitely many nonzero entries.

        The weight of a tableau `T` is the same as the weight of the
        reading word of `T`, for any reading order.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).weight()
            [1, 1, 1, 1]

            sage: Tableau([]).weight()
            []

            sage: Tableau([[1,3,3,7],[4,2],[2,3]]).weight()
            [1, 2, 3, 1, 0, 0, 1]

        TESTS:

        We check that this agrees with going to the word::

            sage: t = Tableau([[1,3,4,7],[6,2],[2,3]])
            sage: def by_word(T):
            ....:     ed = T.to_word().evaluation_dict()
            ....:     m = max(ed.keys()) + 1
            ....:     return [ed.get(k,0) for k in range(1,m)]
            sage: by_word(t) == t.weight()
            True
            sage: SST = SemistandardTableaux(shape=[3,1,1])
            sage: all(by_word(t) == t.weight() for t in SST)
            True
        """
        if len(self) == 0:
            return []
        m = max(max(row) for row in self)
        res = [0] * m
        for row in self:
            for i in row:
                res[i - 1] += 1
        return res

    evaluation = weight

    def is_row_strict(self):
        """
        Return ``True`` if ``self`` is a row strict immaculate tableau
        and ``False`` otherwise.

        A tableau is row strict if the entries in each row are in increasing
        order.

        EXAMPLES::

            sage: Tableau([[1, 3], [2, 4]]).is_row_strict()
            True
            sage: Tableau([[1, 2], [2, 4]]).is_row_strict()
            True
            sage: Tableau([[2, 3], [2, 4]]).is_row_strict()
            True
            sage: Tableau([[5, 3], [2, 4]]).is_row_strict()
            False
        """
        return all(row[i]<row[i+1] for row in self for i in range(len(row)-1))

    def is_column_strict(self):
        """
        Return ``True`` if ``self`` is a column strict tableau and ``False``
        otherwise.

        A tableau is column strict if the entries in each column are in
        increasing order.

        EXAMPLES::

            sage: Tableau([[1, 3], [2, 4]]).is_column_strict()
            True
            sage: Tableau([[1, 2], [2, 4]]).is_column_strict()
            True
            sage: Tableau([[2, 3], [2, 4]]).is_column_strict()
            False
            sage: Tableau([[5, 3], [2, 4]]).is_column_strict()
            False
        """
        return all(self[r-1][c]<self[r][c] for (r,c) in self.cells() if r>0)

    def is_standard(self):
        """
        Return ``True`` if ``self`` is a standard immaculate tableau
        and ``False`` otherwise.

        EXAMPLES::

            sage: Tableau([[1, 3], [2, 4]]).is_standard()
            True
            sage: Tableau([[1, 2], [2, 4]]).is_standard()
            False
            sage: Tableau([[2, 3], [2, 4]]).is_standard()
            False
            sage: Tableau([[5, 3], [2, 4]]).is_standard()
            False
        """
        entries=self.entries()
        entries.sort()
        return entries == range(1,self.size()+1) and self.is_row_strict() and self.is_column_strict()

    def is_increasing(self):
        """
        Return ``True`` if ``self`` is an increasing tableau and
        ``False`` otherwise.

        A tableau is increasing if it is both row strict and column strict.

        EXAMPLES::

            sage: Tableau([[1, 3], [2, 4]]).is_increasing()
            True
            sage: Tableau([[1, 2], [2, 4]]).is_increasing()
            True
            sage: Tableau([[2, 3], [2, 4]]).is_increasing()
            False
            sage: Tableau([[5, 3], [2, 4]]).is_increasing()
            False
            sage: Tableau([[1, 2, 3], [2, 3], [3]]).is_increasing()
            True
        """
        return self.is_row_strict() and self.is_column_strict()

    def is_rectangular(self):
        """
        Return ``True`` if the tableau ``self`` is rectangular and
        ``False`` otherwise.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).is_rectangular()
            True
            sage: Tableau([[1,2,3],[4,5],[6]]).is_rectangular()
            False
            sage: Tableau([]).is_rectangular()
            True
        """
        if len(self) == 0:
            return True
        width = len(self[0])
        for row in self:
            if len(row) != width:
                return False
        return True

    def cells(self):
        """
        Return a list of the coordinates of the cells of ``self``.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).cells()
            [(0, 0), (0, 1), (1, 0), (1, 1)]
        """
        s = []
        for i in range(len(self)):
            s += [ (i,j) for j in range(len(self[i])) ]
        return s

    def cells_containing(self, i):
        r"""
        Return the list of cells in which the letter `i` appears in the
        tableau ``self``. The list is ordered with cells appearing from
        left to right.

        Cells are given as pairs of coordinates `(a, b)`, where both
        rows and columns are counted from `0` (so `a = 0` means the cell
        lies in the leftmost column of the tableau, etc.).

        EXAMPLES::

            sage: t = Tableau([[1,1,3],[2,3,5],[4,5]])
            sage: t.cells_containing(5)
            [(2, 1), (1, 2)]
            sage: t.cells_containing(4)
            [(2, 0)]
            sage: t.cells_containing(6)
            []

            sage: t = Tableau([[1,1,2,4],[2,4,4],[4]])
            sage: t.cells_containing(4)
            [(2, 0), (1, 1), (1, 2), (0, 3)]

            sage: t = Tableau([[1,1,2,8,9],[2,5,6,11],[3,7,7,13],[4,8,9],[5],[13],[14]])
            sage: t.cells_containing(8)
            [(3, 1), (0, 3)]

            sage: Tableau([]).cells_containing(3)
            []
        """
        cell_list = []
        for r in range(len(self)-1, -1, -1):
            rth_row = self[r]
            for c,val in enumerate(rth_row):
                if val == i:
                    cell_list.append((r,c))
        return cell_list

    def restrict(self, n):
        """
        Return the restriction of the (standard) tableau to `n`. If possible,
        the restricted tableau will have the same parent as this tableau.

        EXAMPLES::

            sage: Tableau([[1,2],[3],[4]]).restrict(3)
            [[1, 2], [3]]
            sage: StandardTableau([[1,2],[3],[4]]).restrict(2)
            [[1, 2]]

        If possible the restricted tableau will belong to the same category as
        the original tableau::

            sage: S=StandardTableau([[1,2,4,7],[3,5],[6]]); S.category()
            Category of elements of Standard tableaux
            sage: S.restrict(4).category()
            Category of elements of Standard tableaux
            sage: SS=StandardTableaux([4,2,1])([[1,2,4,7],[3,5],[6]]); SS.category()
            Category of elements of Standard tableaux of shape [4, 2, 1]
            sage: SS.restrict(4).category()
            Category of elements of Standard tableaux

            sage: Tableau([[1,2],[3],[4]]).restrict(3)
            [[1, 2], [3]]
            sage: Tableau([[1,2],[3],[4]]).restrict(2)
            [[1, 2]]
            sage: SemistandardTableau([[1,1],[2]]).restrict(1)
            [[1, 1]]
            sage: _.category()
            Category of elements of Semistandard tableaux
        """
        res = [ [y for y in row if y <= n] for row in self]
        res = [row for row in res if row != []]
        # attempt to return a tableau of the same type
        try:
            return self.parent()( res )
        except StandardError:
            try:
                return self.parent().Element( res )
            except StandardError:
                return Tableau(res)

    def to_chain(self):
        """
        Return the chain of compositions corresponding to the (semi)standard
        immaculate tableau.

        EXAMPLES::

            sage: Tableau([[1,2],[3],[4]]).to_chain()
            [[], [1], [2], [2, 1], [2, 1, 1]]
            sage: Tableau([[1,1],[2]]).to_chain()
            [[], [2], [2, 1]]
            sage: Tableau([[1,1],[3]]).to_chain()
            [[], [2], [2], [2, 1]]
            sage: Tableau([]).to_chain()
            [[]]
        """
        if self == []:
            return [self.shape()]
        m = max(self.to_word())
        return [self.restrict(k).shape() for k in range(m+1)]

    def to_list(self):
        """
        EXAMPLES::

            sage: t = Tableau([[1,2],[3,4]])
            sage: l = t.to_list(); l
            [[1, 2], [3, 4]]
            sage: l[0][0] = 2
            sage: t
            [[1, 2], [3, 4]]
        """
        return [row[:] for row in self]

    def slide_multiply(left, right):
        """
        Multiply two tableaux using jeu de taquin.

        This product makes the set of tableaux into an associative monoid.
        The empty tableau is the unit in this monoid.

        See pp. 15 of [Ful1997]_.

        EXAMPLES::

            sage: t = Tableau([[1,2,2,3],[2,3,5,5],[4,4,6],[5,6]])
            sage: t2 = Tableau([[1,2],[3]])
            sage: t.slide_multiply(t2)
            [[1, 1, 2, 2, 3], [2, 2, 3, 5], [3, 4, 5], [4, 6, 6], [5]]
        """
        st = []
        if len(left) == 0:
            return right
        else:
            l = len(left[0])

        for row in range(len(right)):
            st.append([None]*l + right[row])
        for row in range(len(left)):
            st.append(left[row])

        from sage.combinat.skew_tableau import SkewTableau
        return SkewTableau(st).rectify()

    def _slide_up(self, c):
        r"""
        Auxiliary method used for promotion, which removes cell `c` from ``self``,
        slides the letters of ``self`` up using jeu de taquin slides, and
        then fills the empty cell at `(0,0)` with the value `0`.

        TESTS::

            sage: t = Tableau([[1,1,2],[2,3,5],[4,5]])
            sage: t._slide_up((2,1))
            [[0, 1, 2], [1, 3, 5], [2, 4]]

            sage: t._slide_up((1,2))
            [[0, 1, 2], [1, 2, 3], [4, 5]]

            sage: t = Tableau([[1,1,3],[2,3,5],[4,5]])
            sage: t._slide_up((1,2))
            [[0, 1, 1], [2, 3, 3], [4, 5]]
        """
        new_st = [x[:] for x in self]
        spotl, spotc = c
        while [spotl, spotc] != [0,0]:
            #once moving box is in first column, just move letters up
            #(French notation!)
            if spotc == 0:
                new_st[spotl][spotc] = new_st[spotl-1][spotc]
                spotl -= 1
                continue
            #once moving box is in first row, just move letters up
            elif spotl == 0:
                new_st[spotl][spotc] = new_st[spotl][spotc-1]
                spotc -= 1
                continue
            else:
                #If we get to this stage, we need to compare
                below = new_st[spotl-1][spotc]
                left = new_st[spotl][spotc-1]
                if below >= left:
                    #Swap with the cell below
                    new_st[spotl][spotc] = new_st[spotl-1][spotc]
                    spotl -= 1
                    continue
                else:
                    #Swap with the cell to the left
                    new_st[spotl][spotc] = new_st[spotl][spotc-1]
                    spotc -= 1
                    continue
        #set box in position (0,0) to 0
        new_st[0][0] = 0
        return Tableau(new_st)

    def _slide_down(self, c, n):
        r"""
        Auxiliary method used for promotion, which removes cell `c` from ``self``,
        slides the letters of ``self`` down using jeu de taquin slides, and
        then fills the empty cell with the value `n + 2`.

        When the entries of ``self`` are positive integers, and cell `c` is
        filled with `1`, then the position of `c` is irrelevant.

        TESTS::

            sage: t = Tableau([[1,1,2],[2,3,5],[4,5]])
            sage: t._slide_down((0, 0), 8)
            [[1, 2, 5], [2, 3, 10], [4, 5]]

            sage: t._slide_down((0, 1), 8)
            [[1, 2, 5], [2, 3, 10], [4, 5]]

            sage: t = Tableau([[1,1,2,2,2,3],[2,2,4,6,6],[4,4,5,7],[5,8]])
            sage: t._slide_down((0, 1), 9)
            [[1, 2, 2, 2, 2, 3], [2, 4, 4, 6, 6], [4, 5, 7, 11], [5, 8]]
        """
        new_st = [x[:] for x in self]
        #new_st is a deep copy of self, so as not to mess around with self.
        new_st_shape = [len(x) for x in self]
        spotl, spotc = c
        #spotl and spotc are the coordinates of the wandering hole.
        #All comments and variable names below refer to French notation.
        while True:
            #"right_neighbor" and "upper_neighbor" refer to neighbors of the
            #hole.
            go_right = None
            if len(new_st_shape) > spotl + 1 and new_st_shape[spotl + 1] >= spotc + 1:
                upper_neighbor = new_st[spotl + 1][spotc]
                go_right = False
            if new_st_shape[spotl] != spotc + 1:
                right_neighbor = new_st[spotl][spotc + 1]
                if go_right is None or upper_neighbor > right_neighbor:
                    go_right = True
            if go_right == True:
                new_st[spotl][spotc] = right_neighbor
                spotc += 1
            elif go_right == False:
                new_st[spotl][spotc] = upper_neighbor
                spotl += 1
            else:
                break
        new_st[spotl][spotc] = n + 2
        return Tableau(new_st)

    def promotion_inverse(self, n):
        """
        Return the image of ``self`` under the inverse promotion operator.

        The inverse promotion operator, applied to a tableau `t`, does the
        following:

        Iterate over all letters `1` in the tableau `t`, from right to left.
        For each of these letters, do the following:

        - Remove the letter from `t`, thus leaving a hole where it used to be.

        - Apply jeu de taquin to move this hole northeast (in French notation)
          until it reaches the outer boundary of `t`.

        - Fill `n+2` into the hole once jeu de taquin has completed.

        Once this all is done, subtract `1` from each letter in the tableau.
        This is not always well-defined. Restricted to the class of
        semistandard tableaux whose entries are all `\leq n + 1`, this is the
        usual inverse promotion operator defined on this class.

        When ``self`` is a standard tableau of size ``n + 1``, this definition of
        inverse promotion is the map called "promotion" in [Sg2011]_ (p. 23) and
        in [St2009]_, and is the inverse of the map called "promotion" in
        [Hai1992]_ (p. 90).

        EXAMPLES::

            sage: t = Tableau([[1,2],[3,3]])
            sage: t.promotion_inverse(2)
            [[1, 2], [2, 3]]

            sage: t = Tableau([[1,2],[2,3]])
            sage: t.promotion_inverse(2)
            [[1, 1], [2, 3]]

            sage: t = Tableau([[1,2,5],[3,3,6],[4,7]])
            sage: t.promotion_inverse(8)
            [[1, 2, 4], [2, 5, 9], [3, 6]]

            sage: t = Tableau([])
            sage: t.promotion_inverse(2)
            []

        TESTS:

        We check the equivalence of two definitions of inverse promotion
        on semistandard tableaux::

            sage: ST = SemistandardTableaux(shape=[4,2,1], max_entry=7)
            sage: def bk_promotion_inverse7(st):
            ....:     st2 = st
            ....:     for i in range(1, 7):
            ....:         st2 = st2.bender_knuth_involution(i, check=False)
            ....:     return st2
            sage: all( bk_promotion_inverse7(st) == st.promotion_inverse(6) for st in ST ) # long time
            True
            sage: ST = SemistandardTableaux(shape=[2,2,2], max_entry=7)
            sage: all( bk_promotion_inverse7(st) == st.promotion_inverse(6) for st in ST ) # long time
            True

        A test for :trac:`13203`::

            sage: T = Tableau([[1]])
            sage: type(T.promotion_inverse(2)[0][0])
            <type 'sage.rings.integer.Integer'>
        """
        if self.is_rectangular():
            n = Integer(n)
            if self.size() == 0:
                return self
            s = self.shape()[0]
            l = self.weight()[0]
            word = [i for i in self.to_word() if i>1]
            word = [i-1 for i in word]
            t = Tableau([])
            t = t.insert_word(word)
            t = t.to_list()
            if l < s:
                for i in range(l):
                    t[len(t)-1].append(n+1)
            else:
                t.append([n+1 for i in range(s)])
            return Tableau(t)
        # Now, the non-rectangular case.
        p = self
        for c in reversed(self.cells_containing(1)):
            p = p._slide_down(c, n)
        return Tableau([[i-1 for i in row] for row in p])

    def promotion(self, n):
        r"""
        Return the image of ``self`` under the promotion operator.

        The promotion operator, applied to a tableau `t`, does the following:

        Iterate over all letters `n+1` in the tableau `t`, from left to right.
        For each of these letters, do the following:

        - Remove the letter from `t`, thus leaving a hole where it used to be.

        - Apply jeu de taquin to move this hole southwest (in French notation)
          until it reaches the inner boundary of `t`.

        - Fill `0` into the hole once jeu de taquin has completed.

        Once this all is done, add `1` to each letter in the tableau.
        This is not always well-defined. Restricted to the class of
        semistandard tableaux whose entries are all `\leq n + 1`, this is the
        usual promotion operator defined on this class.

        When ``self`` is a standard tableau of size ``n + 1``, this definition of
        promotion is precisely the one given in [Hai1992]_ (p. 90). It is the
        inverse of the maps called "promotion" in [Sg2011]_ (p. 23) and in [St2009]_.

        REFERENCES:

        .. [Hai1992] Mark D. Haiman,
           *Dual equivalence with applications, including a conjecture of Proctor*,
           Discrete Mathematics 99 (1992), 79-113,
           http://www.sciencedirect.com/science/article/pii/0012365X9290368P

        .. [Sg2011] Bruce E. Sagan,
           *The cyclic sieving phenomenon: a survey*,
           :arXiv:`1008.0790v3`

        EXAMPLES::

            sage: t = Tableau([[1,2],[3,3]])
            sage: t.promotion(2)
            [[1, 1], [2, 3]]

            sage: t = Tableau([[1,1,1],[2,2,3],[3,4,4]])
            sage: t.promotion(3)
            [[1, 1, 2], [2, 2, 3], [3, 4, 4]]

            sage: t = Tableau([[1,2],[2]])
            sage: t.promotion(3)
            [[2, 3], [3]]

            sage: t = Tableau([[1,1,3],[2,2]])
            sage: t.promotion(2)
            [[1, 2, 2], [3, 3]]

            sage: t = Tableau([[1,1,3],[2,3]])
            sage: t.promotion(2)
            [[1, 1, 2], [2, 3]]

            sage: t = Tableau([])
            sage: t.promotion(2)
            []

        TESTS:

        We check the equivalence of two definitions of promotion on
        semistandard tableaux::

            sage: ST = SemistandardTableaux(shape=[3,2,2,1], max_entry=6)
            sage: def bk_promotion6(st):
            ....:     st2 = st
            ....:     for i in range(5, 0, -1):
            ....:         st2 = st2.bender_knuth_involution(i, check=False)
            ....:     return st2
            sage: all( bk_promotion6(st) == st.promotion(5) for st in ST ) # long time
            True
            sage: ST = SemistandardTableaux(shape=[4,4], max_entry=6)
            sage: all( bk_promotion6(st) == st.promotion(5) for st in ST ) # long time
            True

        We also check :meth:`promotion_inverse()` is the inverse
        of :meth:`promotion()`::

            sage: ST = SemistandardTableaux(shape=[3,2,1], max_entry=7)
            sage: all( st.promotion(6).promotion_inverse(6) == st for st in ST ) # long time
            True
        """
        if self.is_rectangular():
            t = self.rotate_180()
            t = [[n+2-i for i in row] for row in t.to_list()]
            t = Tableau(t).promotion_inverse(n)
            t = [[n+2-i for i in row] for row in t.to_list()]
            return Tableau(t).rotate_180()
        p = self
        for c in self.cells_containing(n+1):
            p = p._slide_up(c)
        return Tableau([[i+1 for i in row] for row in p])

    def row_stabilizer(self):
        """
        Return the PermutationGroup corresponding to the row stabilizer of
        ``self``.

        This assumes that every integer from `1` to the size of ``self``
        appears exactly once in ``self``.

        EXAMPLES::

            sage: rs = Tableau([[1,2,3],[4,5]]).row_stabilizer()
            sage: rs.order() == factorial(3)*factorial(2)
            True
            sage: PermutationGroupElement([(1,3,2),(4,5)]) in rs
            True
            sage: PermutationGroupElement([(1,4)]) in rs
            False
            sage: rs = Tableau([[1, 2],[3]]).row_stabilizer()
            sage: PermutationGroupElement([(1,2),(3,)]) in rs
            True
            sage: rs.one().domain()
            [1, 2, 3]
            sage: rs = Tableau([[1],[2],[3]]).row_stabilizer()
            sage: rs.order()
            1
            sage: rs = Tableau([[2,4,5],[1,3]]).row_stabilizer()
            sage: rs.order()
            12
            sage: rs = Tableau([]).row_stabilizer()
            sage: rs.order()
            1
        """

        # Ensure that the permutations involve all elements of the
        # tableau, by including the identity permutation on the set [1..k].
        k = self.size()
        gens = [range(1,k+1)]
        for i in range(len(self)):
            for j in range(0, len(self[i])-1):
                gens.append( (self[i][j], self[i][j+1]) )
        return PermutationGroup( gens )


    def height(self):
        """
        Return the height of ``self``.

        EXAMPLES::

            sage: Tableau([[1,2,3],[4,5]]).height()
            2
            sage: Tableau([[1,2,3]]).height()
            1
            sage: Tableau([]).height()
            0
        """
        return len(self)

    def _heights(self):
        """
        EXAMPLES::

            sage: Tableau([[1,2,3,4],[5,6],[7],[8]])._heights()
            [1, 3, 4, 4]
            sage: Tableau([])._heights()
            []
            sage: Tableau([[1]])._heights()
            [1]
            sage: Tableau([[1,2]])._heights()
            [1, 1]
            sage: Tableau([[1,2],[3],[4]])._heights()
            [1, 3]
        """
        cor = self.corners()
        ncor = len(cor)
        if ncor == 0:
            return []
        k = len(self)
        cor = [ [k-i,j+1]  for i,j in reversed(cor)]

        heights = [1]*(cor[0][1])
        for i in range(1, ncor):
            heights += [ cor[i][0] ]*(cor[i][1]-cor[i-1][1])

        return heights

    def last_letter_lequal(self, tab2):
        """
        Return ``True`` if ``self`` is less than or equal to ``tab2`` in the last
        letter ordering.

        EXAMPLES::

            sage: st = StandardTableaux([3,2])
            sage: f = lambda b: 1 if b else 0
            sage: matrix( [ [ f(t1.last_letter_lequal(t2)) for t2 in st] for t1 in st] )
            [1 1 1 1 1]
            [0 1 1 1 1]
            [0 0 1 1 1]
            [0 0 0 1 1]
            [0 0 0 0 1]
        """
        n = self.size()
        if not isinstance(tab2, Tableau):
            try:
                tab2 = Tableau(tab2)
            except StandardError:
                raise TypeError("tab2 must be a standard tableau")

        if tab2.size() != n:
            raise ValueError("tab2 must be the same size as self")

        if self == tab2:
            return True

        for j in range(n, 1, -1):
            self_j_pos = None
            for i in range(len(self)):
                if j in self[i]:
                    self_j_pos = i
                    break

            tab2_j_pos = None
            for i in range(len(tab2)):
                if j in tab2[i]:
                    tab2_j_pos = i
                    break

            if self_j_pos < tab2_j_pos:
                return True
            if tab2_j_pos < self_j_pos:
                return False

    def charge(self):
        r"""
        Return the charge of the reading word of ``self``.  See
        :meth:`sage.combinat.words.finite_word.FiniteWord_class.charge`
        for more information.

        EXAMPLES::

            sage: Tableau([[1,1],[2,2],[3]]).charge()
            0
            sage: Tableau([[1,1,3],[2,2]]).charge()
            1
            sage: Tableau([[1,1,2],[2],[3]]).charge()
            1
            sage: Tableau([[1,1,2],[2,3]]).charge()
            2
            sage: Tableau([[1,1,2,3],[2]]).charge()
            2
            sage: Tableau([[1,1,2,2],[3]]).charge()
            3
            sage: Tableau([[1,1,2,2,3]]).charge()
            4
        """
        return self.to_word().charge()

    def cocharge(self):
        r"""
        Return the cocharge of the reading word of ``self``.  See
        :meth:`sage.combinat.words.finite_word.FiniteWord_class.cocharge`
        for more information.

        EXAMPLES::

            sage: Tableau([[1,1],[2,2],[3]]).cocharge()
            4
            sage: Tableau([[1,1,3],[2,2]]).cocharge()
            3
            sage: Tableau([[1,1,2],[2],[3]]).cocharge()
            3
            sage: Tableau([[1,1,2],[2,3]]).cocharge()
            2
            sage: Tableau([[1,1,2,3],[2]]).cocharge()
            2
            sage: Tableau([[1,1,2,2],[3]]).cocharge()
            1
            sage: Tableau([[1,1,2,2,3]]).cocharge()
            0
        """
        return self.to_word().cocharge()


    def add_entry(self,cell,m):
        """
        Set the entry in ``cell`` equal to ``m``. If the cell does not exist then
        extend the tableau, otherwise just replace the entry.

        EXAMPLES::

            sage: s=StandardTableau([[1,2,5],[3,4]]); s.pp()
              1  2  5
              3  4
            sage: t=s.add_entry( (1,2), 6); t.pp()
              1  2  5
              3  4  6
            sage: t.category()
            Category of elements of Standard tableaux
            sage: s.add_entry( (2,0), 6).pp()
              1  2  5
              3  4
              6
            sage: u=s.add_entry( (1,2), 3); u.pp()
              1  2  5
              3  4  3
            sage: u.category()
            Category of elements of Tableaux
            sage: s.add_entry( (2,2),3)
            Traceback (most recent call last):
            ...
            IndexError: (2, 2) is not an addable cell of the tableau

        """
        tab=self.to_list()
        (r,c)=cell
        try:
            tab[r][c]=m   # will work if we are replacing an entry
        except IndexError:
            # Only add a new row if (r,c) is an addable cell (previous code
            # added added m to the end of row r independently of the value of c)
            if (r,c) in self.shape().outside_corners():  # an addable node
                if r==len(tab):
                    tab.append([])

                tab[r].append(m)

            else:
                raise IndexError, '%s is not an addable cell of the tableau' % ((r,c),)

        # attempt to return a tableau of the same type as self
        if tab in self.parent():
            return self.parent()(tab)
        else:
            try:
                return self.parent().Element(tab)
            except StandardError:
                return Tableau(tab)


    def promotion_operator(self, i):
        """
        EXAMPLES::

            sage: t = Tableau([[1,2],[3]])
            sage: t.promotion_operator(1)
            [[[1, 2], [3], [4]], [[1, 2], [3, 4]], [[1, 2, 4], [3]]]
            sage: t.promotion_operator(2)
            [[[1, 1], [2, 3], [4]],
             [[1, 1, 2], [3], [4]],
             [[1, 1, 4], [2, 3]],
             [[1, 1, 2, 4], [3]]]
            sage: Tableau([[1]]).promotion_operator(2)
            [[[1, 1], [2]], [[1, 1, 2]]]
            sage: Tableau([[1,1],[2]]).promotion_operator(3)
            [[[1, 1, 1], [2, 2], [3]],
             [[1, 1, 1, 2], [2], [3]],
             [[1, 1, 1, 3], [2, 2]],
             [[1, 1, 1, 2, 3], [2]]]

        TESTS::

            sage: Tableau([]).promotion_operator(2)
            [[[1, 1]]]
            sage: Tableau([]).promotion_operator(1)
            [[[1]]]
        """
        chain = self.to_chain()
        part = self.shape()
        weight = self.weight()
        perm = permutation.from_reduced_word(range(1, len(weight)+1))
        l = part.add_horizontal_border_strip(i)
        ltab = [ from_chain( chain + [next] ) for next in l ]
        return [ x.symmetric_group_action_on_values(perm) for x in ltab]


    ##################################
    # actions on tableaux from words #
    ##################################
    def raise_action_from_words(self, f, *args):
        """
        EXAMPLES::

            sage: from sage.combinat.tableau import symmetric_group_action_on_values
            sage: import functools
            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: f = functools.partial(t.raise_action_from_words, symmetric_group_action_on_values)
            sage: f([1,2,3])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: f([3,2,1])
            [[1, 1, 1, 1], [2, 3], [3]]
            sage: f([1,3,2])
            [[1, 1, 2, 2], [2, 2], [3]]
        """
        w = self.to_word()
        w = f(w, *args)
        return from_shape_and_word(self.shape(), w)

    def symmetric_group_action_on_values(self, perm):
        """
        EXAMPLES::

            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: t.symmetric_group_action_on_values([1,2,3])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t.symmetric_group_action_on_values([3,2,1])
            [[1, 1, 1, 1], [2, 3], [3]]
            sage: t.symmetric_group_action_on_values([1,3,2])
            [[1, 1, 2, 2], [2, 2], [3]]
        """
        return self.raise_action_from_words(symmetric_group_action_on_values, perm)


class SemistandardTableau(Tableau):
    """
    A class to model a semistandard tableau.

    INPUT:

    - ``t`` -- a tableau, a list of lists, or an empty list

    OUTPUT:

    - A SemistandardTableau object constructed from ``t``.

    A semistandard tableau is a tableau whose entries are positive integers,
    which are weakly increasing in rows and strictly increasing down columns.

    EXAMPLES::

        sage: t = SemistandardTableau([[1,2,3],[2,3]]); t
        [[1, 2, 3], [2, 3]]
        sage: t.shape()
        [3, 2]
        sage: t.pp() # pretty print
        1 2 3
        2 3
        sage: t = Tableau([[1,2],[2]])
        sage: s = SemistandardTableau(t); s
        [[1, 2], [2]]
        sage: SemistandardTableau([]) # The empty tableau
        []

    When using code that will generate a lot of tableaux, it is slightly more
    efficient to construct a SemistandardTableau from the appropriate
    :class:`Parent` object::

        sage: SST = SemistandardTableaux()
        sage: SST([[1, 2, 3], [4, 5]])
        [[1, 2, 3], [4, 5]]

    .. SEEALSO:

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`

    TESTS::

        sage: SemistandardTableau([[1,2,3],[1]])
        Traceback (most recent call last):
        ...
        ValueError: [[1, 2, 3], [1]] is not a column strict tableau

        sage: SemistandardTableau([[1,2,1]])
        Traceback (most recent call last):
        ...
        ValueError: The rows of [[1, 2, 1]] are not weakly increasing

        sage: SemistandardTableau([[0,1]])
        Traceback (most recent call last):
        ...
        ValueError: entries must be positive integers
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a SemistandardTableau is only ever constructed as an
        element_class call of an appropriate parent.

        TESTS::

            sage: t = SemistandardTableau([[1,1],[2]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Semistandard tableaux
            sage: t.category()
            Category of elements of Semistandard tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.SemistandardTableaux_all_with_category.element_class'>
        """
        if isinstance(t, SemistandardTableau):
            return t
        elif t in SemistandardTableaux():
            return SemistandardTableaux_all().element_class(SemistandardTableaux_all(), t)

        # t is not a semistandard tableau so we give an appropriate error message
        if t not in Tableaux():
            raise ValueError('%s is not a tableau' % t)

        if not all(isinstance(c,(int,Integer)) and c>0 for row in t for c in row):
            raise ValueError("entries must be positive integers"%t)

        valid_entries=range(1,1+max(sum((list(row) for row in t),[])))
        if not all(c in valid_entries for row in t for c in row):
            raise ValueError("the entries must be in %s"%(t,valid_entries))

        if any(row[c]>row[c+1] for row in t for c in range(len(row)-1)):
            raise ValueError("The rows of %s are not weakly increasing"%t)

        # If we're still here ``t`` cannot be column strict
        raise ValueError('%s is not a column strict tableau' % t)


    def __init__(self, parent, t):
        r"""
        Initializes a semistandard tableau.

        TESTS::

            sage: t = Tableaux()([[1,1],[2]])
            sage: s = SemistandardTableaux(3)([[1,1],[2]])
            sage: s==t
            True
            sage: s.parent()
            Semistandard tableaux of size 3 and maximum entry 3
            sage: r = SemistandardTableaux(3)(t); r.parent()
            Semistandard tableaux of size 3 and maximum entry 3
            sage: isinstance(r, Tableau)
            True
        """
        super(SemistandardTableau, self).__init__(parent, t)

        # Tableau() has checked that t is tableau, so it remains to check that
        # the entries of t are positive integers
        from sage.sets.positive_integers import PositiveIntegers
        if any(c not in PositiveIntegers() for row in t for c in row):
            raise ValueError("the entries of a semistandard tableau must be non-negative integers")

        # which are weakly increasing along rows
        if any(row[c]>row[c+1] for row in t for c in xrange(len(row)-1)):
            raise ValueError("the entries in each row of a semistandard tableau must be weakly increasing")

        # and strictly increasing down columns
        if len(t)>0 and any(t[r][c] >= t[r+1][c] for c in xrange(len(t[0])) for r in xrange(len(t)-1) if len(t[r+1])>c):
            raise ValueError("the entries of each column of a semistandard tableau must be strictly increasing")

class StandardTableau(SemistandardTableau):
    """
    A class to model a standard tableau.

    INPUT:

    - ``t`` -- a Tableau, a list of lists, or an empty list

    OUTPUT:

    - A StandardTableau object constructed from ``t``.

    A standard tableau is a semistandard tableau whose entries are exactly the
    positive integers from 1 to `n`, where `n` is the size of the tableau.

    EXAMPLES::

        sage: t = StandardTableau([[1,2,3],[4,5]]); t
        [[1, 2, 3], [4, 5]]
        sage: t.shape()
        [3, 2]
        sage: t.pp() # pretty print
        1 2 3
        4 5
        sage: t.is_standard()
        True
        sage: StandardTableau([]) # The empty tableau
        []

    When using code that will generate a lot of tableaux, it is slightly more
    efficient to construct a StandardTableau from the appropriate
    :class:`Parent` object::

        sage: ST = StandardTableaux()
        sage: ST([[1, 2, 3], [4, 5]])
        [[1, 2, 3], [4, 5]]

    .. SEEALSO:

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`

        sage: StandardTableau([[1,2,3],[4,4]])
        Traceback (most recent call last):
        ...
        ValueError: the entries in each row of a standard tableau must be strictly increasing
        sage: StandardTableau([[1,3,2]])
        Traceback (most recent call last):
        ...
        ValueError: the entries in each row of a semistandard tableau must be weakly increasing
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a :class:`StandardTableau` is only ever constructed
        as an ``element_class`` call of an appropriate parent.

        TESTS::

            sage: t = StandardTableau([[1,2],[3]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Standard tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.StandardTableaux_all_with_category.element_class'>
        """
        if isinstance(t, StandardTableau):
            return t

        return StandardTableaux_all().element_class(StandardTableaux_all(), t)

    def __init__(self, parent, t):
        r"""
        Initializes a standard tableau.

        TESTS::

            sage: t = Tableaux()([[1,2],[3]])
            sage: s = StandardTableaux(3)([[1,2],[3]])
            sage: s==t
            True
            sage: s.parent()
            Standard tableaux of size 3
            sage: r = StandardTableaux(3)(t); r.parent()
            Standard tableaux of size 3
            sage: isinstance(r, Tableau)
            True
        """
        super(StandardTableau, self).__init__(parent, t)

        # t is semistandard so we only need to check that it is standard
        if any(row[c]==row[c+1] for row in self for c in xrange(len(row)-1)):
            raise ValueError("the entries in each row of a standard tableau must be strictly increasing")

        # and that the entries are in bijection with {1,2,...,n}
        if sorted(flatten(self._list))!=range(1,self.size()+1):
            raise ValueError("the entries in a standard tableau must be in bijection with 1,2,...,n")



    def content(self, k, multicharge=[0]):
        """
        Returns the content of ``k`` in a standard tableau. That is, if
        ``k`` appears in row `r` and column `c` of the tableau then we
        return `c-r`.

        The ``multicharge`` is a list of length 1 which gives an offset for
        all of the contents. It is included mainly for compatibility with
        :class:`TableauTuple`.

        EXAMPLES::

            sage: StandardTableau([[1,2],[3,4]]).content(3)
            -1

            sage: StandardTableau([[1,2],[3,4]]).content(6)
            Traceback (most recent call last):
            ...
            ValueError: 6 does not appear in tableau
        """
        for r in range(len(self)):
          try:
            return self[r].index(k) - r + multicharge[0]
          except ValueError:
            pass
        raise ValueError("%d does not appear in tableau"%k)

    def dominates(self, t):
        r"""
        Return ``True`` if ``self`` dominates the tableau ``t``. That is,
        if the shape of the tableau restricted to `k` dominates the shape of
        ``t`` restricted to `k`, for `k = 1, 2, \ldots, n`.

        When the two tableaux have the same shape, then this ordering
        coincides with the Bruhat ordering for the corresponding permutations.

        INPUT:

        - ``t`` -- A tableau

        EXAMPLES::

            sage: s=StandardTableau([[1,2,3],[4,5]])
            sage: t=StandardTableau([[1,2],[3,5],[4]])
            sage: s.dominates(t)
            True
            sage: t.dominates(s)
            False
            sage: all(StandardTableau(s).dominates(t) for t in StandardTableaux([3,2]))
            True
            sage: s.dominates([[1,2,3,4,5]])
            False

        """
        t=StandardTableau(t)
        return all(self.restrict(m).shape().dominates(t.restrict(m).shape())
                        for m in xrange(1,1+self.size()))

    def is_standard(self):
        """
        Return ``True`` since ``self`` is a standard tableau.

        EXAMPLES::

            sage: StandardTableau([[1, 3], [2, 4]]).is_standard()
            True
        """
        return True

    def up(self):
        """
        An iterator for all the standard tableaux that can be
        obtained from ``self`` by adding a cell.

        EXAMPLES::

            sage: t = StandardTableau([[1,2]])
            sage: [x for x in t.up()]
            [[[1, 2, 3]], [[1, 2], [3]]]
        """
        #Get a list of all places where we can add a cell
        #to the shape of self

        outside_corners = self.shape().outside_corners()

        n = self.size()

        #Go through and add n+1 to the end of each
        #of the rows
        for row, _ in outside_corners:
            new_t = map(list, self)
            if row != len(self):
                new_t[row] += [n+1]
            else:
                new_t.append([n+1])
            yield StandardTableau(new_t)

    def up_list(self):
        """
        Return a list of all the standard tableaux that can be obtained
        from ``self`` by adding a cell.

        EXAMPLES::

            sage: t = StandardTableau([[1,2]])
            sage: t.up_list()
            [[[1, 2, 3]], [[1, 2], [3]]]
        """
        return list(self.up())

    def down(self):
        """
        An iterator for all the standard tableaux that can be obtained
        from ``self`` by removing a cell. Note that this iterates just
        over a single tableau (or nothing if ``self`` is empty).

        EXAMPLES::

            sage: t = StandardTableau([[1,2],[3]])
            sage: [x for x in t.down()]
            [[[1, 2]]]
            sage: t = StandardTableau([])
            sage: [x for x in t.down()]
            []
        """
        if len(self) > 0:
            yield self.restrict( self.size() - 1 )

    def down_list(self):
        """
        Return a list of all the standard tableaux that can be obtained
        from ``self`` by removing a cell. Note that this is just a singleton
        list if ``self`` is nonempty, and an empty list otherwise.

        EXAMPLES::

            sage: t = StandardTableau([[1,2],[3]])
            sage: t.down_list()
            [[[1, 2]]]
            sage: t = StandardTableau([])
            sage: t.down_list()
            []
        """
        return list(self.down())

    def standard_descents(self):
        """
        Return a list of the integers `i` such that `i` appears
        strictly further north than `i + 1` in ``self`` (this is not
        to say that `i` and `i + 1` must be in the same column). The
        list is sorted in increasing order.

        EXAMPLES::

            sage: StandardTableau( [[1,3,4],[2,5]] ).standard_descents()
            [1, 4]
            sage: StandardTableau( [[1,2],[3,4]] ).standard_descents()
            [2]
            sage: StandardTableau( [[1,2,5],[3,4],[6,7],[8],[9]] ).standard_descents()
            [2, 5, 7, 8]
            sage: StandardTableau( [] ).standard_descents()
            []
        """
        descents = []
        #whatpart gives the number for which self is a partition
        whatpart = sum(i for i in self.shape())
        #now find the descents
        for i in range(1, whatpart):
            #find out what row i and i+1 are in (we're using the
            #standardness of self here)
            for j in range(len(self)):
                if self[j].count(i+1) > 0:
                    break
                if self[j].count(i) > 0:
                    descents.append(i)
                    break
        return descents

    def standard_number_of_descents(self):
        """
        Return the number of all integers `i` such that `i` appears
        strictly further north than `i + 1` in ``self`` (this is not
        to say that `i` and `i + 1` must be in the same column). A
        list of these integers can be obtained using the
        :meth:`standard_descents` method.

        EXAMPLES::

            sage: StandardTableau( [[1,2],[3,4],[5]] ).standard_number_of_descents()
            2
            sage: StandardTableau( [] ).standard_number_of_descents()
            0
            sage: tabs = StandardTableaux(5)
            sage: all( t.standard_number_of_descents() == t.schuetzenberger_involution().standard_number_of_descents() for t in tabs )
            True
        """
        return len(self.standard_descents())

    def standard_major_index(self):
        """
        Return the major index of the standard tableau ``self`` in the
        standard meaning of the word. The major index is defined to be
        the sum of the descents of ``self`` (see :meth:`standard_descents`
        for their definition).

        EXAMPLES::

            sage: StandardTableau( [[1,4,5],[2,6],[3]] ).standard_major_index()
            8
            sage: StandardTableau( [[1,2],[3,4]] ).standard_major_index()
            2
            sage: StandardTableau( [[1,2,3],[4,5]] ).standard_major_index()
            3
        """
        return sum(self.standard_descents())

    def promotion_inverse(self, m=None):
        """
        Return the image of ``self`` under the inverse promotion operator.
        The optional variable `m` should be set to the size of ``self`` minus
        `1` for a minimal speedup; otherwise, it defaults to this number.

        The inverse promotion operator, applied to a standard tableau `t`,
        does the following:

        Remove the letter `1` from `t`, thus leaving a hole where it used to be.
        Apply jeu de taquin to move this hole northeast (in French notation)
        until it reaches the outer boundary of `t`. Fill `n + 1` into this hole,
        where `n` is the size of `t`. Finally, subtract `1` from each letter in
        the tableau. This yields a new standard tableau.

        This definition of inverse promotion is the map called "promotion" in
        [Sg2011]_ (p. 23) and in [St2009]_, and is the inverse of the map
        called "promotion" in [Hai1992]_ (p. 90).

        See the :meth:`sage.combinat.tableau.promotion_inverse` method for a
        more general operator.

        EXAMPLES::

            sage: t = StandardTableau([[1,3],[2,4]])
            sage: t.promotion_inverse()
            [[1, 2], [3, 4]]

        We check the equivalence of two definitions of inverse promotion on
        standard tableaux::

            sage: ST = StandardTableaux(7)
            sage: def bk_promotion_inverse7(st):
            ....:     st2 = st
            ....:     for i in range(1, 7):
            ....:         st2 = st2.bender_knuth_involution(i, check=False)
            ....:     return st2
            sage: all( bk_promotion_inverse7(st) == st.promotion_inverse() for st in ST ) # long time
            True
        """
        if m is None:
            m = self.size() - 1
        return StandardTableau(Tableau(self.to_list()).promotion_inverse(m))

    def promotion(self, m=None):
        r"""
        Return the image of ``self`` under the promotion operator.

        The promotion operator, applied to a standard tableau `t`, does the
        following:

        Remove the letter `n` from `t`, thus leaving a hole where it used to be.
        Apply jeu de taquin to move this hole southwest (in French notation)
        until it reaches the inner boundary of `t`. Fill `0` into the hole once
        jeu de taquin has completed. Finally, add `1` to each letter in the
        tableau. The resulting standard tableau is the image of `t` under the
        promotion operator.

        This definition of promotion is precisely the one given in [Hai1992]_
        (p. 90). It is the inverse of the maps called "promotion" in [Sg2011]_
        (p. 23) and in [St2009]_.

        See the :meth:`sage.combinat.tableau.promotion` method for a
        more general operator.

        EXAMPLES::

            sage: ST = StandardTableaux(7)
            sage: all( st.promotion().promotion_inverse() == st for st in ST ) # long time
            True
            sage: all( st.promotion_inverse().promotion() == st for st in ST ) # long time
            True
            sage: st = StandardTableau([[1,2,5],[3,4]])
            sage: parent(st.promotion())
            Standard tableaux
        """
        if m is None:
            m = self.size() - 1
        return StandardTableau(Tableau(self.to_list()).promotion(m))

def from_chain(chain):
    """
    Returns a semistandard tableau from a chain of partitions.

    EXAMPLES::

        sage: from sage.combinat.tableau import from_chain
        sage: from_chain([[], [2], [2, 1], [3, 2, 1]])
        [[1, 1, 3], [2, 3], [3]]
    """
    res = [[0]*chain[-1][i] for i in range(len(chain[-1]))]
    for i in reversed(range(2, len(chain)+1)):
        for j in range(len(chain[i-1])):
            for k in range(chain[i-1][j]):
                res[j][k] = i -1
    return Tableau(res)

@rename_keyword(deprecation=13605, order='convention')
def from_shape_and_word(shape, w, convention="French"):
    r"""
    Returns a tableau from a shape and word.

    INPUT:

    - ``shape`` -- a partition

    - ``w`` -- a word whose length equals that of the partition

    - ``convention`` -- a string which can take values ``"French"`` or
      ``"English"``; the default is ``"French"``

    OUTPUT:

    A tableau, whose shape is ``shape`` and whose reading word is ``w``.
    If the ``convention`` is specified as ``"French"``, the reading word is to be read
    starting from the top row in French convention (= the bottom row in English
    convention). If the ``convention`` is specified as ``"English"``, the reading word
    is to be read starting with the top row in English convention.

    EXAMPLES::

        sage: from sage.combinat.tableau import from_shape_and_word
        sage: t = Tableau([[1, 3], [2], [4]])
        sage: shape = t.shape(); shape
        [2, 1, 1]
        sage: word = t.to_word(); word
        word: 4213
        sage: from_shape_and_word(shape, word)
        [[1, 3], [2], [4]]
        sage: word = Word(flatten(t))
        sage: from_shape_and_word(shape, word, convention = "English")
        [[1, 3], [2], [4]]
    """
    res = []
    j = 0
    if convention == "French":
        shape = reversed(shape)
    for l in shape:
        res.append( list(w[j:j+l]) )
        j += l
    if convention == "French":
        res.reverse()
    return Tableau(res)

class ImmaculateTableaux(UniqueRepresentation, Parent):
    """
    The class of immaculate tableaux.

    INPUT:

    - ``n`` (optional) -- a non-negative integer

    OUTPUT:

    - If ``n`` is specified, the class of tableaux of size ``n``. Otherwise,
      the class of all tableaux.

    A tableau in Sage is a finite list of lists, whose lengths are weakly
    decreasing, or an empty list, representing the empty tableau.  The entries
    of a tableau can be any sage object. Because of this, no enumeration
    through the set of Tableaux is possible.

    EXAMPLES::

        sage: T = Tableaux(); T
        Tableaux
        sage: T3 = Tableaux(3); T3
        Tableaux of size 3
        sage: [['a','b']] in T
        True
        sage: [['a','b']] in T3
        False
        sage: t = T3([[1,1,1]]); t
        [[1, 1, 1]]
        sage: t in T
        True
        sage: t.parent()
        Tableaux of size 3
        sage: T([]) # the empty tableau
        []
        sage: T.category()
        Category of sets

    .. SEEALSO:

        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`

    TESTS::

        sage: TestSuite( Tableaux() ).run()
        sage: TestSuite( Tableaux(5) ).run()
        sage: t = Tableaux(3)([[1,2],[3]])
        sage: t.parent()
        Tableaux of size 3
        sage: Tableaux(t)
        Traceback (most recent call last):
        ...
        ValueError: The argument to Tableaux() must be a non-negative integer.
        sage: Tableaux(3)([[1, 1]])
        Traceback (most recent call last):
        ...
        ValueError: [[1, 1]] is not an element of Tableaux of size 3.

        sage: t0 = Tableau([[1]])
        sage: t1 = Tableaux()([[1]])
        sage: t2 = Tableaux()(t1)
        sage: t0 == t1 == t2
        True
        sage: t1 in Tableaux()
        True
        sage: t1 in Tableaux(1)
        True
        sage: t1 in Tableaux(2)
        False

        sage: [[1]] in Tableaux()
        True
        sage: [] in Tableaux(0)
        True

    Check that trac:`14145` has been fixed::

        sage: 1 in Tableaux()
        False
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`Tableaux` for more
        information.

        TESTS::

            sage: Tableaux()
            Tableaux
            sage: Tableaux(3)
            Tableaux of size 3
        """
        if args:
            n = args[0]
        elif 'n' in kwargs:
            n = kwargs[n]
        else:
            n = None

        if n == None:
            return Tableaux_all()
        else:
            if not isinstance(n,(int, Integer)) or n < 0:
                raise ValueError( "The argument to Tableaux() must be a non-negative integer." )
            return Tableaux_size(n)

    Element = Tableau
    global_options = TableauOptions

    def _element_constructor_(self, t):
        r"""
        Constructs an object from ``t`` as an element of ``self``, if
        possible. This is inherited by all Tableaux, SemistandardTableaux, and
        StandardTableaux classes.

        INPUT:

        - ``t`` -- Data which can be interpreted as a tableau

        OUTPUT:

        - The corresponding tableau object

        TESTS::

            sage: T = Tableaux(3)
            sage: T([[1,2,1]]).parent() is T     # indirect doctest
            True
            sage: T( StandardTableaux(3)([[1, 2, 3]])).parent() is T
            True
            sage: T([[1,2]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 2]] is not an element of Tableaux of size 3.
        """
        if not t in self:
            raise ValueError("%s is not an element of %s."%(t, self))

        return self.element_class(self, t)

    def __contains__(self, x):
        """
        TESTS::

            sage: T = sage.combinat.tableau.Tableaux()
            sage: [[1,2],[3,4]] in T
            True
            sage: [[1,2],[3]] in T
            True
            sage: [] in T
            True
            sage: [['a','b']] in T
            True
            sage: Tableau([['a']]) in T
            True

            sage: [1,2,3] in T
            False
            sage: [[1],[1,2]] in T
            False

        Check that :trac:`14145` is fixed::

            sage: 1 in sage.combinat.tableau.Tableaux()
            False
        """
        from sage.combinat.partition import _Partitions
        if isinstance(x, Tableau):
            return True
        elif isinstance(x, list) and all(isinstance(y, list) for y in x):
            # any list of lists of partition shape is a tableau
            return map(len,x) in _Partitions
        else:
            return False

class ImmaculateTableaux_all(ImmaculateTableaux):
    """
    All immaculate tableaux.
    """
    def __init__(self):
        r"""
        Initializes the class of all Immaculate tableaux.

        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_all()
            sage: TestSuite(T).run()

        """
        super(ImmaculateTableaux_all, self).__init__(category=Sets())

    def _repr_(self):
        """
        TESTS::

            sage: ImmaculateTableaux()
            Immaculate tableaux
        """
        return "Immaculate tableaux"

    def an_element(self):
        r"""
        Return a particular element of the class.

        TESTS::

            sage: T = Tableaux()
            sage: T.an_element()
            [[1, 1], [1]]
        """
        return self.element_class(self, [[1, 1], [1]])


class ImmaculateTableaux_size(ImmaculateTableaux):
    """
    Immaculate tableaux of a fixed size `n`.
    """
    def __init__(self, n):
        r"""
        Initialize ``self``.

        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_size(3)
            sage: TestSuite(T).run()

            sage: T = sage.combinat.tableau.Tableaux_size(0)
            sage: TestSuite(T).run()
        """
        super(ImmaculateTableaux_size, self).__init__(category=Sets())
        self.size = n

    def __contains__(self,x):
        """
        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_size(3)
            sage: [[2,4], [1]] in T
            True

            sage: [[2,4],[1,3]] in T
            False

        Check that :trac:`14145` is fixed::

            sage: 1 in sage.combinat.tableau.Tableaux_size(3)
            False
        """
        return ImmaculateTableaux.__contains__(self, x) and sum(map(len,x)) == self.size

    def _repr_(self):
        """
        TESTS::

            sage: ImmaculateTableaux(4)
            Immaculate tableaux of size 4
        """
        return "Immaculate tableaux of size {}".format(self.size)

    def an_element(self):
        r"""
        Return a particular element of the class.

        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_size(3)
            sage: T.an_element()
            [[1, 1], [1]]
            sage: T = sage.combinat.tableau.Tableaux_size(0)
            sage: T.an_element()
            []
        """
        if self.size==0:
            return self.element_class(self, [])

        if self.size==1:
            return self.element_class(self, [[1]])

        return self.element_class(self, [[1]*(self.size-1),[1]])


##########################
# Semi-standard tableaux #
##########################
class SemistandardTableaux(Tableaux):
    """
    A factory class for the various classes of semistandard tableaux.

    INPUT:

    Keyword arguments:

    - ``size`` -- The size of the tableaux
    - ``shape`` -- The shape of the tableaux
    - ``eval`` -- The weight (also called content or weight) of the tableaux
    - ``max_entry`` -- A maximum entry for the tableaux.  This can be a
      positive integer or infinity (``oo``). If ``size`` or ``shape`` are
      specified, ``max_entry`` defaults to be ``size`` or the size of
      ``shape``.

    Positional arguments:

    - The first argument is interpreted as either ``size`` or ``shape``
      according to  whether it is an integer or a partition
    - The second keyword argument will always be interpreted as ``eval``

    OUTPUT:

    - The appropriate class, after checking basic consistency tests. (For
      example, specifying ``eval`` implies a value for `max_entry`).

    A semistandard tableau is a tableau whose entries are positive integers,
    which are weakly increasing in rows and strictly increasing down columns.
    Note that Sage uses the English convention for partitions and tableaux;
    the longer rows are displayed on top.

    Classes of semistandard tableaux can be iterated over if and only if there
    is some restriction.

    EXAMPLES::

        sage: SST = SemistandardTableaux([2,1]); SST
        Semistandard tableaux of shape [2, 1] and maximum entry 3
        sage: SST.list()
        [[[1, 1], [2]],
         [[1, 1], [3]],
         [[1, 2], [2]],
         [[1, 2], [3]],
         [[1, 3], [2]],
         [[1, 3], [3]],
         [[2, 2], [3]],
         [[2, 3], [3]]]

        sage: SST = SemistandardTableaux(3); SST
        Semistandard tableaux of size 3 and maximum entry 3
        sage: SST.list()
        [[[1, 1, 1]],
         [[1, 1, 2]],
         [[1, 1, 3]],
         [[1, 2, 2]],
         [[1, 2, 3]],
         [[1, 3, 3]],
         [[2, 2, 2]],
         [[2, 2, 3]],
         [[2, 3, 3]],
         [[3, 3, 3]],
         [[1, 1], [2]],
         [[1, 1], [3]],
         [[1, 2], [2]],
         [[1, 2], [3]],
         [[1, 3], [2]],
         [[1, 3], [3]],
         [[2, 2], [3]],
         [[2, 3], [3]],
         [[1], [2], [3]]]

        sage: SST = SemistandardTableaux(3, max_entry=2); SST
        Semistandard tableaux of size 3 and maximum entry 2
        sage: SST.list()
        [[[1, 1, 1]],
         [[1, 1, 2]],
         [[1, 2, 2]],
         [[2, 2, 2]],
         [[1, 1], [2]],
         [[1, 2], [2]]]

        sage: SST = SemistandardTableaux(3, max_entry=oo); SST
        Semistandard tableaux of size 3
        sage: SST[123]
        [[3, 4], [6]]

        sage: SemistandardTableaux(max_entry=2)[11]
        [[1, 1], [2]]

        sage: SemistandardTableaux()[0]
        []

    .. SEEALSO:

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`SemistandardTableaux`
        for more information.

        TESTS::

            sage: SemistandardTableaux()
            Semistandard tableaux
            sage: SemistandardTableaux(3)
            Semistandard tableaux of size 3 and maximum entry 3
            sage: SemistandardTableaux(size=3)
            Semistandard tableaux of size 3 and maximum entry 3
            sage: SemistandardTableaux(0)
            Semistandard tableaux of size 0 and maximum entry 0
            sage: SemistandardTableaux([2,1])
            Semistandard tableaux of shape [2, 1] and maximum entry 3
            sage: SemistandardTableaux(shape=[2,1])
            Semistandard tableaux of shape [2, 1] and maximum entry 3
            sage: SemistandardTableaux([])
            Semistandard tableaux of shape [] and maximum entry 0
            sage: SemistandardTableaux(eval=[2,1])
            Semistandard tableaux of size 3 and weight [2, 1]
            sage: SemistandardTableaux(max_entry=3)
            Semistandard tableaux with maximum entry 3
            sage: SemistandardTableaux(3, [2,1])
            Semistandard tableaux of size 3 and weight [2, 1]
            sage: SemistandardTableaux(3, shape=[2,1])
            Semistandard tableaux of shape [2, 1] and maximum entry 3
            sage: SemistandardTableaux(3, [2,1], shape=[2,1])
            Semistandard tableaux of shape [2, 1] and weight [2, 1]
            sage: SemistandardTableaux(3, max_entry=4)
            Semistandard tableaux of size 3 and maximum entry 4
            sage: SemistandardTableaux(3, max_entry=oo)
            Semistandard tableaux of size 3
            sage: SemistandardTableaux([2, 1], max_entry=oo)
            Semistandard tableaux of shape [2, 1]
            sage: SemistandardTableaux([2, 1], [2, 1])
            Semistandard tableaux of shape [2, 1] and weight [2, 1]
            sage: mu = Partition([2,1]); SemistandardTableaux(mu, mu)
            Semistandard tableaux of shape [2, 1] and weight [2, 1]
            sage: SemistandardTableaux(3, [2, 1], max_entry=2)
            Semistandard tableaux of size 3 and weight [2, 1]

            sage: SemistandardTableaux(3, shape=[2])
            Traceback (most recent call last):
            ...
            ValueError: size and shape are different sizes

            sage: SemistandardTableaux(3, [2])
            Traceback (most recent call last):
            ...
            ValueError: size and eval are different sizes

            sage: SemistandardTableaux([2],[3])
            Traceback (most recent call last):
            ...
            ValueError: shape and eval are different sizes

            sage: SemistandardTableaux(2,[2], max_entry=4)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must match the weight

            sage: SemistandardTableaux(eval=[2], max_entry=oo)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must match the weight

            sage: SemistandardTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: shape must be a (skew) partition
        """
        from sage.combinat.partition import Partition, _Partitions
        # Process the keyword arguments -- allow for original syntax where
        #   n == size,  p== shape and mu == eval
        n = kwargs.get('n', None)
        size = kwargs.get('size', n)

        p = kwargs.get('p', None)
        shape = kwargs.get('shape', p)

        mu = kwargs.get('eval', None)
        mu = kwargs.get("mu", mu)

        max_entry = kwargs.get('max_entry', None)

        # Process the positional arguments
        if args:
            # The first arg could be either a size or a shape
            if isinstance(args[0], (int, Integer)):
                if size is not None:
                    raise ValueError( "size was specified more than once" )
                else:
                    size = args[0]
            else:
                if shape is not None:
                    raise ValueError( "the shape was specified more than once" )
                shape = args[0] # we check it's a partition later

        if len(args) == 2:
            # The second non-keyword argument is the weight
            if mu is not None:
                raise ValueError( "the weight was specified more than once" )
            else:
                mu = args[1]

        # Consistency checks
        if size is not None:
            if not isinstance(size, (int, Integer)):
                raise ValueError( "size must be an integer" )
            elif size < 0:
                raise ValueError( "size must be non-negative" )

        if shape is not None:
            from sage.combinat.skew_partition import SkewPartitions
            # use in (and not isinstance) below so that lists can be used as
            # shorthand
            if shape in _Partitions:
                shape = Partition(shape)
            elif shape in SkewPartitions():
                from sage.combinat.skew_tableau import SemistandardSkewTableaux
                return SemistandardSkewTableaux(shape, mu)
            else:
                raise ValueError( "shape must be a (skew) partition" )

        if mu is not None:
            if (not mu in Compositions()) and\
                    (not mu in _Partitions):
                raise ValueError( "mu must be a composition" )
            mu = Composition(mu)

        is_inf = max_entry is PlusInfinity()

        if max_entry is not None:
            if not is_inf and not isinstance(max_entry, (int, Integer)):
                raise ValueError( "max_entry must be an integer or PlusInfinity" )
            elif max_entry <= 0:
                raise ValueError( "max_entry must be positive" )

        if (mu is not None) and (max_entry is not None):
            if max_entry != len(mu):
                raise ValueError( "the maximum entry must match the weight" )

        if (size is not None) and (shape is not None):
            if sum(shape) != size:
                # This could return an empty class instead of an error
                raise ValueError( "size and shape are different sizes" )

        if (size is not None) and (mu is not None):
            if sum(mu) != size:
                # This could return an empty class instead of an error
                raise ValueError( "size and eval are different sizes" )

        # Dispatch appropriately
        if (shape is not None) and (mu is not None):
            if sum(shape) != sum(mu):
                # This could return an empty class instead of an error
                raise ValueError( "shape and eval are different sizes" )
            else:
                return SemistandardTableaux_shape_weight(shape, mu)

        if (shape is not None):
            if is_inf:
                return SemistandardTableaux_shape_inf(shape)
            return SemistandardTableaux_shape(shape, max_entry)

        if (mu is not None):
            return SemistandardTableaux_size_weight(sum(mu), mu)

        if (size is not None):
            if is_inf:
                return SemistandardTableaux_size_inf(size)
            return SemistandardTableaux_size(size, max_entry)

        return SemistandardTableaux_all(max_entry)

    Element = SemistandardTableau

    def __init__(self, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SemistandardTableaux()
            sage: TestSuite(S).run()
        """
        if kwds.has_key('max_entry'):
            self.max_entry = kwds['max_entry']
            kwds.pop('max_entry')
        else:
            self.max_entry = None
        Tableaux.__init__(self, **kwds)

    def __getitem__(self, r):
        r"""
        The default implementation of ``__getitem``__ for enumerated sets
        does not allow slices so we override it.

        EXAMPLES::

            sage: StandardTableaux([4,3,3,2])[10:20]     # indirect doctest
            [[[1, 3, 9, 12], [2, 5, 10], [4, 6, 11], [7, 8]],
             [[1, 2, 9, 12], [3, 5, 10], [4, 6, 11], [7, 8]],
             [[1, 3, 9, 12], [2, 4, 10], [5, 6, 11], [7, 8]],
             [[1, 2, 9, 12], [3, 4, 10], [5, 6, 11], [7, 8]],
             [[1, 5, 8, 12], [2, 6, 10], [3, 7, 11], [4, 9]],
             [[1, 4, 8, 12], [2, 6, 10], [3, 7, 11], [5, 9]],
             [[1, 3, 8, 12], [2, 6, 10], [4, 7, 11], [5, 9]],
             [[1, 2, 8, 12], [3, 6, 10], [4, 7, 11], [5, 9]],
             [[1, 4, 8, 12], [2, 5, 10], [3, 7, 11], [6, 9]],
             [[1, 3, 8, 12], [2, 5, 10], [4, 7, 11], [6, 9]]]

            sage: SemistandardTableaux(size=2, max_entry=oo)[5]
            [[2, 3]]

            sage: SemistandardTableaux([2,1], max_entry=oo)[3]
            [[1, 2], [3]]

            sage: SemistandardTableaux(3, max_entry=2)[0:5]    # indirect doctest
            [[[1, 1, 1]],
            [[1, 1, 2]],
            [[1, 2, 2]],
            [[2, 2, 2]],
            [[1, 1], [2]]]

            sage: SemistandardTableaux([2,2], [2, 1, 1])[0]    # indirect doctest
            [[1, 1], [2, 3]]

            sage: SemistandardTableaux([1,1,1], max_entry=4)[0:4]
            [[[1], [2], [3]],
             [[1], [2], [4]],
             [[1], [3], [4]],
             [[2], [3], [4]]]

            sage: SemistandardTableaux(3, [2,1])[1]    # indirect doctest
            [[1, 1], [2]]

            sage: StandardTableaux(3)[:]  # indirect doctest
            [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]]]

            sage: StandardTableaux([2,2])[1]   # indirect doctest
            [[1, 2], [3, 4]]

        TESTS::

            sage: SemistandardTableaux()[5]
            [[1], [2]]

            sage: SemistandardTableaux(max_entry=2)[5]
            [[2, 2]]

            sage: SemistandardTableaux()[:]
            Traceback (most recent call last):
            ...
            ValueError: infinite set

            sage: SemistandardTableaux(size=2, max_entry=oo)[:]
            Traceback (most recent call last):
            ...
            ValueError: infinite set
        """
        if isinstance(r,(int,Integer)):
            return self.unrank(r)
        elif isinstance(r,slice):
            start=0 if r.start is None else r.start
            stop=r.stop
            if stop is None and not self.is_finite():
                raise ValueError( 'infinite set' )
        else:
            raise ValueError( 'r must be an integer or a slice' )
        count=0
        tabs=[]
        for t in self:
            if count==stop:
                break
            if count>=start:
                tabs.append(t)
            count+=1

        # this is to cope with empty slices endpoints like [:6] or [:}
        if count==stop or stop is None:
            return tabs
        raise IndexError('value out of range')

    def __contains__(self, t):
        """
        Return ``True`` if ``t`` can be interpreted as a
        :class:`SemistandardTableau`.

        TESTS::

            sage: T = sage.combinat.tableau.SemistandardTableaux_all()
            sage: [[1,2],[2]] in T
            True
            sage: [] in T
            True
            sage: Tableau([[1]]) in T
            True
            sage: StandardTableau([[1]]) in T
            True

            sage: [[1,2],[1]] in T
            False
            sage: [[1,1],[5]] in T
            True

        Check that :trac:`14145` is fixed::

            sage: 1 in sage.combinat.tableau.SemistandardTableaux()
            False
        """
        if isinstance(t, SemistandardTableau):
            return self.max_entry is None or \
                    len(t) == 0 or \
                    max(flatten(t)) <= self.max_entry
        elif t == []:
            return True
        elif Tableaux.__contains__(self, t) and all(c>0 for row in t for c in row) \
                and all(row[i] <= row[i+1] for row in t for i in range(len(row)-1)) \
                and all(t[r][c] < t[r+1][c] for c in range(len(t[0]))
                        for r in range(len(t)-1) if len(t[r+1]) > c):
            return self.max_entry is None or max(flatten(t)) <= self.max_entry
        else:
            return False

class SemistandardTableaux_all(SemistandardTableaux, DisjointUnionEnumeratedSets):
    """
    All semistandard tableaux.

    .. WARNING::

        Input is not checked; please use :class:`SemistandardTableaux` to
        ensure the options are properly parsed.
    """
    def __init__(self, max_entry=None):
        r"""
        Initializes the class of all semistandard tableaux.

        TESTS::

            sage: T = sage.combinat.tableau.SemistandardTableaux_all()
            sage: TestSuite(T).run()

            sage: T=sage.combinat.tableau.SemistandardTableaux_all(max_entry=3)
            sage: TestSuite(T).run()
        """
        if max_entry is not PlusInfinity():
            self.max_entry = max_entry
            SST_n = lambda n: SemistandardTableaux_size(n, max_entry)
            DisjointUnionEnumeratedSets.__init__( self,
                    Family(NonNegativeIntegers(), SST_n),
                    facade=True, keepkey = False)

        else:
            self.max_entry = None

    def _repr_(self):
        """
        TESTS::

            sage: SemistandardTableaux()    # indirect doctest
            Semistandard tableaux

            sage: SemistandardTableaux(max_entry=3)
            Semistandard tableaux with maximum entry 3
        """
        if self.max_entry is not None:
            return "Semistandard tableaux with maximum entry %s"%str(self.max_entry)
        return "Semistandard tableaux"


    def list(self):
        """
        TESTS::

            sage: SemistandardTableaux().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


class SemistandardTableaux_size_inf(SemistandardTableaux):
    """
    Semistandard tableaux of fixed size `n` with no maximum entry.

    .. WARNING::

        Input is not checked; please use :class:`SemistandardTableaux` to
        ensure the options are properly parsed.
    """
    def __init__(self, n):
        r"""
        Initializes the class of semistandard tableaux of size ``n`` with no
        maximum entry.

        TESTS::

            sage: T = sage.combinat.tableau.SemistandardTableaux_size_inf(3)
            sage: TestSuite(T).run()
        """
        super(SemistandardTableaux_size_inf, self).__init__(
              category = InfiniteEnumeratedSets())
        self.size = n


    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux(3, max_entry=oo))    # indirect doctest
            'Semistandard tableaux of size 3'
        """
        return "Semistandard tableaux of size %s"%str(self.size)

    def __contains__(self, t):
        """
        Return ``True`` if ``t`` can be interpreted as an element of this
        class.

        TESTS::

            sage: T = SemistandardTableaux(3, max_entry=oo)
            sage: [[1,2],[5]] in T
            True
            sage: StandardTableau([[1, 2], [3]]) in T
            True

            sage: [] in T
            False
            sage: Tableau([[1]]) in T
            False

        Check that :trac:`14145` is fixed::

            sage: 1 in SemistandardTableaux(3, max_entry=oo)
            False
        """
        return SemistandardTableaux.__contains__(self, t) and sum(map(len, t)) == self.size

    def __iter__(self):
        """
        EXAMPLES::

            sage: sst = SemistandardTableaux(3, max_entry=oo)
            sage: [sst[t] for t in range(0,5)]
            [[[1, 1, 1]],
             [[1, 1, 2]],
             [[1, 2, 2]],
             [[2, 2, 2]],
             [[1, 1], [2]]]
            sage: sst[1000]
            [[2, 12], [7]]
            sage: sst[0].parent() is sst
            True
        """
        from sage.combinat.partition import Partitions
        # Iterates through with maximum entry as order
        i = 1
        while(True):
            for part in Partitions(self.size):
                if i != 1:
                    for k in range(1, self.size+1):
                        for c in IntegerVectors(self.size - k, i-1):
                            c.append(k)
                            for sst in SemistandardTableaux_shape_weight(part, Composition(c)):
                                yield self.element_class(self, sst)
                else:
                    for sst in SemistandardTableaux_shape_weight(part, Composition([self.size])):
                        yield self.element_class(self, sst)
            i += 1


    def list(self):
        """
        TESTS::

            sage: SemistandardTableaux(3, max_entry=oo).list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


class SemistandardTableaux_shape_inf(SemistandardTableaux):
    """
    Semistandard tableaux of fixed shape `p` and no maximum entry.

    .. WARNING::

        Input is not checked; please use :class:`SemistandardTableaux` to
        ensure the options are properly parsed.
    """
    def __init__(self, p):
        r"""
        Initializes the class of semistandard tableaux of shape ``p`` and no
        maximum entry.

        TESTS::

            sage: SST = SemistandardTableaux([2,1], max_entry=oo)
            sage: type(SST)
            <class 'sage.combinat.tableau.SemistandardTableaux_shape_inf_with_category'>
            sage: TestSuite(SST).run()
        """
        super(SemistandardTableaux_shape_inf, self).__init__(
              category = InfiniteEnumeratedSets())
        self.shape = p


    def __contains__(self, x):
        """
        EXAMPLES::

            sage: SST = SemistandardTableaux([2,1], max_entry=oo)
            sage: [[13, 67], [1467]] in SST
            True
            sage: SST = SemistandardTableaux([3,1], max_entry=oo)
            sage: [[13, 67], [1467]] in SST
            False

        Check that :trac:`14145` is fixed::

            sage: SST = SemistandardTableaux([3,1], max_entry=oo)
            sage: 1 in SST
            False
        """
        return SemistandardTableaux.__contains__(self, x) and map(len,x)==self.shape

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux([2,1], max_entry=oo))    # indirect doctest
            'Semistandard tableaux of shape [2, 1]'
        """
        return "Semistandard tableaux of shape %s" %str(self.shape)


    def __iter__(self):
        """
        An iterator for the semistandard partitions of shape ``p`` and no
        maximum entry. Iterates through with maximum entry as order.

        EXAMPLES::

            sage: SST = SemistandardTableaux([3, 1], max_entry=oo)
            sage: SST[1000]
            [[1, 1, 10], [6]]
            sage: [ SST[t] for t in range(0, 5) ]
            [[[1, 1, 1], [2]],
             [[1, 1, 2], [2]],
             [[1, 2, 2], [2]],
             [[1, 1, 1], [3]],
             [[1, 1, 2], [3]]]
            sage: SST[0].parent() is SST
            True
        """
        # Iterates through with maximum entry as order
        i = 1
        n = sum(self.shape)
        while(True):
            if i != 1:
                for k in range(1, n+1):
                    for c in IntegerVectors(n - k, i-1):
                        c.append(k)
                        for sst in SemistandardTableaux_shape_weight(self.shape, Composition(c)):
                            yield self.element_class(self, sst)
            else:
                for sst in SemistandardTableaux_shape_weight(self.shape, Composition([n])):
                    yield self.element_class(self, sst)
            i += 1


class SemistandardTableaux_size(SemistandardTableaux):
    """
    Semistandard tableaux of fixed size `n`.

    .. WARNING::

        Input is not checked; please use :class:`SemistandardTableaux`
        to ensure the options are properly parsed.
    """
    def __init__(self, n, max_entry=None):
        r"""
        Initializes the class of semistandard tableaux of size ``n``.

        TESTS::

            sage: SST = SemistandardTableaux(3); SST
            Semistandard tableaux of size 3 and maximum entry 3
            sage: type(SST)
            <class 'sage.combinat.tableau.SemistandardTableaux_size_with_category'>
            sage: TestSuite(SST).run()

            sage: SST = SemistandardTableaux(3, max_entry=6)
            sage: type(SST)
            <class 'sage.combinat.tableau.SemistandardTableaux_size_with_category'>
            sage: TestSuite(SST).run()
        """

        if max_entry is None:
            max_entry = n
        super(SemistandardTableaux_size, self).__init__(max_entry = max_entry,
                  category = FiniteEnumeratedSets())
        self.size = n

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux(3))    # indirect doctest
            'Semistandard tableaux of size 3 and maximum entry 3'

            sage: repr(SemistandardTableaux(3, max_entry=6))
            'Semistandard tableaux of size 3 and maximum entry 6'
        """
        return "Semistandard tableaux of size %s and maximum entry %s"%(str(self.size), str(self.max_entry))

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[1,2],[3,3]] in SemistandardTableaux(3)
            False
            sage: [[1,2],[3,3]] in SemistandardTableaux(4)
            True
            sage: [[1,2],[3,3]] in SemistandardTableaux(4, max_entry=2)
            False
            sage: SST = SemistandardTableaux(4)
            sage: all([sst in SST for sst in SST])
            True

        Check that :trac:`14145` is fixed::

            sage: SST = SemistandardTableaux(4)
            sage: 1 in SST
            False
        """
        if self.size==0:
            return x == []

        return SemistandardTableaux.__contains__(self, x) \
            and sum(map(len,x)) == self.size and max(flatten(x)) <= self.max_entry

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: SemistandardTableaux(3).cardinality()
            19
            sage: SemistandardTableaux(4).cardinality()
            116
            sage: SemistandardTableaux(4, max_entry=2).cardinality()
            9
            sage: SemistandardTableaux(4, max_entry=10).cardinality()
            4225
            sage: ns = range(1, 6)
            sage: ssts = [ SemistandardTableaux(n) for n in ns ]
            sage: all([sst.cardinality() == len(sst.list()) for sst in ssts])
            True
        """
        from sage.combinat.partition import Partitions
        c = 0
        for part in Partitions(self.size):
            c += SemistandardTableaux_shape(part, self.max_entry).cardinality()
        return c


    def __iter__(self):
        """
        EXAMPLES::

            sage: [ t for t in SemistandardTableaux(2) ]
            [[[1, 1]], [[1, 2]], [[2, 2]], [[1], [2]]]
            sage: [ t for t in SemistandardTableaux(3) ]
            [[[1, 1, 1]],
             [[1, 1, 2]],
             [[1, 1, 3]],
             [[1, 2, 2]],
             [[1, 2, 3]],
             [[1, 3, 3]],
             [[2, 2, 2]],
             [[2, 2, 3]],
             [[2, 3, 3]],
             [[3, 3, 3]],
             [[1, 1], [2]],
             [[1, 1], [3]],
             [[1, 2], [2]],
             [[1, 2], [3]],
             [[1, 3], [2]],
             [[1, 3], [3]],
             [[2, 2], [3]],
             [[2, 3], [3]],
             [[1], [2], [3]]]

            sage: [ t for t in SemistandardTableaux(3, max_entry=2) ]
            [[[1, 1, 1]],
             [[1, 1, 2]],
             [[1, 2, 2]],
             [[2, 2, 2]],
             [[1, 1], [2]],
             [[1, 2], [2]]]

            sage: sst = SemistandardTableaux(3)
            sage: sst[0].parent() is sst
            True
        """
        from sage.combinat.partition import Partitions
        for part in Partitions(self.size):
            for sst in SemistandardTableaux_shape(part, self.max_entry):
                yield self.element_class(self, sst)

class SemistandardTableaux_shape(SemistandardTableaux):
    """
    Semistandard tableaux of fixed shape `p` with a given max entry.

    INPUT:

    - ``p`` -- A partition

    - ``max_entry`` -- The max entry; defaults to the size of ``p``.

    .. WARNING::

        Input is not checked; please use :class:`SemistandardTableaux` to
        ensure the options are properly parsed.
    """
    def __init__(self, p, max_entry=None):
        r"""
        Initializes the class of semistandard tableaux of shape ``p``, with a
        given ``max_entry``.

        TESTS::

            sage: SST = SemistandardTableaux([2,1])
            sage: TestSuite(SST).run()

            sage: SST = SemistandardTableaux([2,1], max_entry=5)
            sage: TestSuite(SST).run()
        """
        if max_entry is None:
            max_entry = sum(p)
        super(SemistandardTableaux_shape, self).__init__(max_entry = max_entry,
              category = FiniteEnumeratedSets())
        self.shape = p

    def __iter__(self):
        """
        An iterator for the semistandard partitions of the specified shape.

        EXAMPLES::

            sage: [ t for t in SemistandardTableaux([3]) ]
            [[[1, 1, 1]],
             [[1, 1, 2]],
             [[1, 1, 3]],
             [[1, 2, 2]],
             [[1, 2, 3]],
             [[1, 3, 3]],
             [[2, 2, 2]],
             [[2, 2, 3]],
             [[2, 3, 3]],
             [[3, 3, 3]]]
            sage: [ t for t in SemistandardTableaux([2,1]) ]
            [[[1, 1], [2]],
             [[1, 1], [3]],
             [[1, 2], [2]],
             [[1, 2], [3]],
             [[1, 3], [2]],
             [[1, 3], [3]],
             [[2, 2], [3]],
             [[2, 3], [3]]]
            sage: [ t for t in SemistandardTableaux([1,1,1]) ]
            [[[1], [2], [3]]]

            sage: [ t for t in SemistandardTableaux([1,1,1], max_entry=4) ]
            [[[1], [2], [3]],
             [[1], [2], [4]],
             [[1], [3], [4]],
             [[2], [3], [4]]]

            sage: sst = SemistandardTableaux([3])
            sage: sst[0].parent() is sst
            True
        """
        for c in IntegerVectors(sum(self.shape), self.max_entry):
            for sst in SemistandardTableaux_shape_weight(self.shape, Composition(c)):
                yield self.element_class(self, sst)


    def __contains__(self, x):
        """
        EXAMPLES::

            sage: SST = SemistandardTableaux([2,1])
            sage: all([sst in SST for sst in SST])
            True
            sage: len(filter(lambda x: x in SST, SemistandardTableaux(3)))
            8
            sage: SST.cardinality()
            8

            sage: SST = SemistandardTableaux([2,1], max_entry=4)
            sage: all([sst in SST for sst in SST])
            True
            sage: SST.cardinality()
            20
        """
        return SemistandardTableaux.__contains__(self, x) and map(len, x) == self.shape

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux([2,1]))    # indirect doctest
            'Semistandard tableaux of shape [2, 1] and maximum entry 3'

            sage: repr(SemistandardTableaux([2,1], max_entry=5))
            'Semistandard tableaux of shape [2, 1] and maximum entry 5'
        """
        return "Semistandard tableaux of shape %s and maximum entry %s" %(str(self.shape), str(self.max_entry))

    def cardinality(self):
        """
        Returns the cardinality of ``self``.

        EXAMPLES::

            sage: SemistandardTableaux([2,1]).cardinality()
            8
            sage: SemistandardTableaux([2,2,1]).cardinality()
            75
            sage: SymmetricFunctions(QQ).schur()([2,2,1]).expand(5)(1,1,1,1,1) # cross check
            75
            sage: SemistandardTableaux([5]).cardinality()
            126
            sage: SemistandardTableaux([3,2,1]).cardinality()
            896

            sage: SemistandardTableaux([3,2,1], max_entry=7).cardinality()
            2352
        """
        c = 0
        for comp in IntegerVectors(sum(self.shape), self.max_entry):
            c += SemistandardTableaux_shape_weight(self.shape, Composition(comp)).cardinality()
        return c

class SemistandardTableaux_shape_weight(SemistandardTableaux_shape):
    r"""
    Semistandard tableaux of fixed shape `p` and weight `\mu`.

    .. WARNING::

        Input is not checked; please use :class:`SemistandardTableaux` to
        ensure the options are properly parsed.
    """
    def __init__(self, p, mu):
        r"""
        Initializes the class of all semistandard tableaux of shape ``p`` and
        weight ``mu``.

        TESTS::

            sage: SST = SemistandardTableaux([2,1], [2,1])
            sage: TestSuite(SST).run()
        """
        super(SemistandardTableaux_shape_weight, self).__init__(p, len(mu))
        self.weight = mu

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux([2,1],[2,1]))    # indirect doctest
            'Semistandard tableaux of shape [2, 1] and weight [2, 1]'
        """
        return "Semistandard tableaux of shape %s and weight %s"%(self.shape, self.weight)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: SST = SemistandardTableaux([2,1], [2,1])
            sage: all([sst in SST for sst in SST])
            True
            sage: len(filter(lambda x: x in SST, SemistandardTableaux(3)))
            1
            sage: SST.cardinality()
            1
        """
        if x not in SemistandardTableaux_shape(self.shape, self.max_entry):
            return False
        n = sum(self.shape)

        if n == 0 and len(x) == 0:
            return True

        content = {}
        for row in x:
            for i in row:
                content[i] = content.get(i, 0) + 1
        content_list = [0]*int(max(content))

        for key in content:
            content_list[key-1] = content[key]

        if content_list != self.weight:
            return False

        return True


    def cardinality(self):
        """
        Returns the number of semistandard tableaux of the given shape and
        weight, as computed by ``kostka_number`` function of symmetrica.

        EXAMPLES::

            sage: SemistandardTableaux([2,2], [2, 1, 1]).cardinality()
            1
            sage: SemistandardTableaux([2,2,2], [2, 2, 1,1]).cardinality()
            1
            sage: SemistandardTableaux([2,2,2], [2, 2, 2]).cardinality()
            1
            sage: SemistandardTableaux([3,2,1], [2, 2, 2]).cardinality()
            2
        """
        return symmetrica.kostka_number(self.shape,self.weight)

    def __iter__(self):
        """
        TESTS::

            sage: sst = SemistandardTableaux([3,1],[2,1,1])
            sage: [sst[i] for i in range(2)]
            [[[1, 1, 2], [3]], [[1, 1, 3], [2]]]
            sage: sst[0].parent() is sst
            True
        """
        for t in symmetrica.kostka_tab(self.shape, self.weight):
            yield self.element_class(self, t)


    def list(self):
        """
        Return a list of semistandard tableau in ``self`` generated by
        semmetrica.

        EXAMPLES::

            sage: SemistandardTableaux([2,2], [2, 1, 1]).list()
            [[[1, 1], [2, 3]]]
            sage: SemistandardTableaux([2,2,2], [2, 2, 1,1]).list()
            [[[1, 1], [2, 2], [3, 4]]]
            sage: SemistandardTableaux([2,2,2], [2, 2, 2]).list()
            [[[1, 1], [2, 2], [3, 3]]]
            sage: SemistandardTableaux([3,2,1], [2, 2, 2]).list()
            [[[1, 1, 2], [2, 3], [3]], [[1, 1, 3], [2, 2], [3]]]
        """
        return symmetrica.kostka_tab(self.shape, self.weight)


class SemistandardTableaux_size_weight(SemistandardTableaux):
    r"""
    Semistandard tableaux of fixed size `n` and weight `\mu`.

    .. WARNING::

        Input is not checked; please use :class:`SemistandardTableaux` to
        ensure the options are properly parsed.
    """
    def __init__(self, n, mu):
        r"""
        Initializes the class of semistandard tableaux of size ``n`` and
        weight ``mu``.

        TESTS::

            sage: SST = SemistandardTableaux(3, [2,1])
            sage: TestSuite(SST).run()
        """
        super(SemistandardTableaux_size_weight, self).__init__(max_entry=len(mu),
              category = FiniteEnumeratedSets())
        self.size = n
        self.weight = mu

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux(3, [2,1]))    # indirect doctest
            'Semistandard tableaux of size 3 and weight [2, 1]'
        """
        return "Semistandard tableaux of size %s and weight %s"%(self.size, self.weight)

    def __iter__(self):
        """
        EXAMPLES::

            sage: [ t for t in SemistandardTableaux(3, [2,1]) ]
            [[[1, 1, 2]], [[1, 1], [2]]]
            sage: [ t for t in SemistandardTableaux(4, [2,2]) ]
            [[[1, 1, 2, 2]], [[1, 1, 2], [2]], [[1, 1], [2, 2]]]
            sage: sst = SemistandardTableaux(4, [2,2])
            sage: sst[0].parent() is sst
            True
        """
        from sage.combinat.partition import Partitions
        for p in Partitions(self.size):
            for sst in SemistandardTableaux_shape_weight(p, self.weight):
                yield self.element_class(self, sst)


    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: SemistandardTableaux(3, [2,1]).cardinality()
            2
            sage: SemistandardTableaux(4, [2,2]).cardinality()
            3
        """
        from sage.combinat.partition import Partitions
        c = 0
        for p in Partitions(self.size):
            c += SemistandardTableaux_shape_weight(p, self.weight).cardinality()
        return c

    def __contains__(self, x):
        """
        TESTS::

            sage: SST = SemistandardTableaux(6, [2,2,2])
            sage: all([sst in SST for sst in SST])
            True
            sage: all([sst in SST for sst in SemistandardTableaux([3,2,1],[2,2,2])])
            True
        """
        from sage.combinat.partition import Partition
        return x in SemistandardTableaux_shape_weight(Partition(map(len,
            x)), self.weight)

########################
# Standard Tableaux    #
########################

class StandardTableaux(SemistandardTableaux):
    """
    A factory for the various classes of standard tableaux.

    INPUT:

    - Either a non-negative integer (possibly specified with the keyword ``n``)
      or a partition.

    OUTPUT:

    - With no argument, the class of all standard tableaux

    - With a non-negative integer argument, ``n``, the class of all standard
      tableaux of size ``n``

    - With a partition argument, the class of all standard tableaux of that
      shape.

    A standard tableau is a semistandard tableaux which contains each of the
    entries from 1 to ``n`` exactly once.

    All classes of standard tableaux are iterable.

    EXAMPLES::

        sage: ST = StandardTableaux(3); ST
        Standard tableaux of size 3
        sage: ST.first()
        [[1, 2, 3]]
        sage: ST.last()
        [[1], [2], [3]]
        sage: ST.cardinality()
        4
        sage: ST.list()
        [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]]]

    .. SEEALSO:

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableau`
        - :class:`StandardSkewTableaux`

    TESTS::

        sage: StandardTableaux()([])
        []
        sage: ST = StandardTableaux([2,2]); ST
        Standard tableaux of shape [2, 2]
        sage: ST.first()
        [[1, 3], [2, 4]]
        sage: ST.last()
        [[1, 2], [3, 4]]
        sage: ST.cardinality()
        2
        sage: ST.list()
        [[[1, 3], [2, 4]], [[1, 2], [3, 4]]]
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`StandardTableaux` for
        more information.

        TESTS::

            sage: StandardTableaux()
            Standard tableaux
            sage: StandardTableaux(3)
            Standard tableaux of size 3
            sage: StandardTableaux([2,1])
            Standard tableaux of shape [2, 1]
            sage: StandardTableaux(0)
            Standard tableaux of size 0

            sage: StandardTableaux(-1)
            Traceback (most recent call last):
            ...
            ValueError: The argument must be a non-negative integer or a partition.
            sage: StandardTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: The argument must be a non-negative integer or a partition.
        """
        from sage.combinat.partition import _Partitions, Partition
        from sage.combinat.skew_partition import SkewPartitions

        if args:
            n = args[0]
        elif 'n' in kwargs:
            n = kwargs[n]
        else:
            n = None

        if n is None:
            return StandardTableaux_all()

        elif n in _Partitions:
            return StandardTableaux_shape(Partition(n))

        elif n in SkewPartitions():
            from sage.combinat.skew_tableau import StandardSkewTableaux
            return StandardSkewTableaux(n)

        if not isinstance(n,(int, Integer)) or n < 0:
            raise ValueError( "The argument must be a non-negative integer or a partition." )

        return StandardTableaux_size(n)

    Element = StandardTableau

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[1,1],[2,3]] in StandardTableaux()
            False
            sage: [[1,2],[3,4]] in StandardTableaux()
            True
            sage: [[1,3],[2,4]] in StandardTableaux()
            True
            sage: [[1,3],[2,5]] in StandardTableaux()
            False
            sage: [] in StandardTableaux()
            True

        Check that :trac:`14145` is fixed::

            sage: 1 in StandardTableaux()
            False
        """
        if isinstance(x, StandardTableau):
            return True
        elif Tableaux.__contains__(self, x):
            flatx = sorted(sum((list(row) for row in x),[]))
            return flatx == range(1,len(flatx)+1) and (len(x)==0 or
                     (all(row[i]<row[i+1] for row in x for i in range(len(row)-1)) and
                       all(x[r][c]<x[r+1][c] for c in range(len(x[0]))
                                             for r in range(len(x)-1) if len(x[r+1])>c )
                     ))
        return False

class StandardTableaux_all(StandardTableaux, DisjointUnionEnumeratedSets):
    """
    All standard tableaux.
    """
    def __init__(self):
        r"""
        Initializes the class of all standard tableaux.

        TESTS::

            sage: ST = StandardTableaux()
            sage: TestSuite(ST).run()
        """
        DisjointUnionEnumeratedSets.__init__( self,
                Family(NonNegativeIntegers(), StandardTableaux_size),
                facade=True, keepkey = False)

    def _repr_(self):
        """
        TESTS::

            sage: repr(StandardTableaux())    # indirect doctest
            'Standard tableaux'
        """
        return "Standard tableaux"


class StandardTableaux_size(StandardTableaux):
    """
    Semistandard tableaux of fixed size `n`.

    .. WARNING::

        Input is not checked; please use :class:`StandardTableaux` to ensure
        the options are properly parsed.
    """
    def __init__(self, n):
        r"""
        Initializes the class of all standard tableaux of size ``n``.

        TESTS::

            sage: TestSuite( StandardTableaux(0) ).run()
            sage: TestSuite( StandardTableaux(3) ).run()
        """
        super(StandardTableaux_size, self).__init__(
              category = FiniteEnumeratedSets())
        self.size = Integer(n)


    def _repr_(self):
        """
        TESTS::

            sage: repr(StandardTableaux(3))    # indirect doctest
            'Standard tableaux of size 3'
        """
        return "Standard tableaux of size %s"%self.size

    def __contains__(self, x):
        """
        TESTS::

            sage: ST3 = StandardTableaux(3)
            sage: all([st in ST3 for st in ST3])
            True
            sage: ST4 = StandardTableaux(4)
            sage: filter(lambda x: x in ST3, ST4)
            []

        Check that :trac:`14145` is fixed::

            sage: 1 in StandardTableaux(4)
            False
        """
        return StandardTableaux.__contains__(self, x) and sum(map(len, x)) == self.size

    def __iter__(self):
        """
        EXAMPLES::

            sage: [ t for t in StandardTableaux(1) ]
            [[[1]]]
            sage: [ t for t in StandardTableaux(2) ]
            [[[1, 2]], [[1], [2]]]
            sage: [ t for t in StandardTableaux(3) ]
            [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]]]
            sage: [ t for t in StandardTableaux(4) ]
            [[[1, 2, 3, 4]],
             [[1, 3, 4], [2]],
             [[1, 2, 4], [3]],
             [[1, 2, 3], [4]],
             [[1, 3], [2, 4]],
             [[1, 2], [3, 4]],
             [[1, 4], [2], [3]],
             [[1, 3], [2], [4]],
             [[1, 2], [3], [4]],
             [[1], [2], [3], [4]]]
            sage: ST4 = StandardTableaux(4)
            sage: ST4[0].parent() is ST4
            True
        """
        from sage.combinat.partition import Partitions
        for p in Partitions(self.size):
            for st in StandardTableaux(p):
                yield self.element_class(self, st)

    def cardinality(self):
        r"""
        Return the cardinality of number of standard tableaux of size
        ``n``.

        The number of standard tableaux of `n` is equal to the number of
        involutions of size `n`. This is a consequence of the symmetry of
        the RSK correspondence, that if `\sigma \mapsto (P, Q)`, then
        `\sigma^{-1} \mapsto (Q, P)`. For more information, see
        :wikipedia:`Robinson-Schensted-Knuth_correspondence#Symmetry`.

        ALGORITHM:

        The algorithm uses the fact that standard tableaux of size
        ``n`` are in bijection with the involutions of size ``n``,
        (see page 41 in section 4.1 of [Ful1997]_).  For each number of
        fixed points, you count the number of ways to choose those
        fixed points multiplied by the number of perfect matchings on
        the remaining values.

        REFERENCES:

        .. [Ful1997] Fulton, William.  Young Tableaux.
           Cambridge University Press, 1997

        EXAMPLES::

            sage: StandardTableaux(3).cardinality()
            4
            sage: ns = [1,2,3,4,5,6]
            sage: sts = [StandardTableaux(n) for n in ns]
            sage: all([st.cardinality() == len(st.list()) for st in sts])
            True
            sage: StandardTableaux(50).cardinality()
            27886995605342342839104615869259776

        TESTS::

            sage: def cardinality_using_hook_formula(n):
            ....:     c = 0
            ....:     for p in Partitions(n):
            ....:         c += StandardTableaux(p).cardinality()
            ....:     return c
            sage: all([cardinality_using_hook_formula(i) == StandardTableaux(i).cardinality() for i in range(10)])
            True
        """
        tableaux_number = self.size % 2  # identity involution
        fixed_point_numbers = xrange(tableaux_number, self.size + 1 - tableaux_number, 2)

        # number of involution of size "size" (number of way to put
        # "fixed_point_number" in "size" box * number of involutions
        # without fixed point of size "size" - "fixed_point_number"
        for fixed_point_number in fixed_point_numbers:
            tableaux_number += (self.size.binomial(fixed_point_number) *
                                prod(range(1, self.size - fixed_point_number, 2)))

        return tableaux_number


class StandardTableaux_shape(StandardTableaux):
    """
    Semistandard tableaux of a fixed shape `p`.

    .. WARNING::

        Input is not checked; please use :class:`SemistandardTableaux` to
        ensure the options are properly parsed.
    """
    def __init__(self, p):
        r"""
        Initializes the class of all semistandard tableaux of a given shape.

        TESTS::

            sage: TestSuite( StandardTableaux([2,1,1]) ).run()
        """
        from sage.combinat.partition import Partition
        super(StandardTableaux_shape, self).__init__(category = FiniteEnumeratedSets())
        self.shape = Partition(p)


    def __contains__(self, x):
        """
        EXAMPLES::

            sage: ST = StandardTableaux([2,1,1])
            sage: all([st in ST for st in ST])
            True
            sage: len(filter(lambda x: x in ST, StandardTableaux(4)))
            3
            sage: ST.cardinality()
            3

        Check that :trac:`14145` is fixed::

            sage: 1 in StandardTableaux([2,1,1])
            False
        """
        return StandardTableaux.__contains__(self, x) and map(len,x) == self.shape

    def _repr_(self):
        """
        TESTS::

            sage: repr(StandardTableaux([2,1,1]))    # indirect doctest
            'Standard tableaux of shape [2, 1, 1]'
        """
        return "Standard tableaux of shape %s"%str(self.shape)

    def cardinality(self):
        r"""
        Returns the number of standard Young tableaux of this shape.

        A formula for the number of Young tableaux associated with a given
        partition. In each cell, write the sum of one plus the number of
        cells horizontally to the right and vertically below the cell (the
        hook length). The number of tableaux is then n! divided by the
        product of all hook lengths.

        For example, consider the partition [3,2,1] of 6 with Ferrers
        Diagram::

            # # #
            # #
            #

        When we fill in the cells with the hook
        lengths, we obtain::

            5 3 1
            3 1
            1

        The hook length formula returns

        .. MATH::

            \frac{6!}{(5 \cdot 3 \cdot 1 \cdot 3 \cdot 1 \cdot 1} = 16.

        EXAMPLES::

            sage: StandardTableaux([3,2,1]).cardinality()
            16
            sage: StandardTableaux([2,2]).cardinality()
            2
            sage: StandardTableaux([5]).cardinality()
            1
            sage: StandardTableaux([6,5,5,3]).cardinality()
            6651216

        REFERENCES:

        - http://mathworld.wolfram.com/HookLengthFormula.html
        """
        pi = self.shape

        number = factorial(sum(pi))
        hook = pi.hook_lengths()

        for row in range(len(pi)):
            for col in range(pi[row]):
                #Divide the hook length by the entry
                number /= hook[row][col]

        return Integer(number)

    def __iter__(self):
        r"""
        An iterator for the standard Young tableaux associated to the
        shape `p` of ``self``.

        EXAMPLES::

            sage: [t for t in StandardTableaux([2,2])]
            [[[1, 3], [2, 4]], [[1, 2], [3, 4]]]
            sage: [t for t in StandardTableaux([3,2])]
            [[[1, 3, 5], [2, 4]],
             [[1, 2, 5], [3, 4]],
             [[1, 3, 4], [2, 5]],
             [[1, 2, 4], [3, 5]],
             [[1, 2, 3], [4, 5]]]
            sage: st = StandardTableaux([2,1])
            sage: st[0].parent() is st
            True
        """

        pi = self.shape
        #Set the initial tableaux by filling it in going down the columns
        tableau = [[None]*n for n in pi]
        size = sum(pi)
        row = 0
        col = 0
        for i in range(size):
            tableau[row][col] = i+1

            #If we can move down, then do it;
            #otherwise, move to the next column over
            if ( row + 1 < len(pi) and col < pi[row+1]):
                row += 1
            else:
                row = 0
                col += 1

        yield self.element_class(self, tableau)

        if self.cardinality() == 1:
            last_tableau = True
        else:
            last_tableau = False

        while not last_tableau:
            #Convert the tableau to "vector format"
            #tableau_vector[i] is the row that number i
            #is in
            tableau_vector = [None]*size
            for row in range(len(pi)):
                for col in range(pi[row]):
                    tableau_vector[tableau[row][col]-1] = row

            #Locate the smallest integer j such that j is not
            #in the lowest corner of the subtableau T_j formed by
            #1,...,j.  This happens to be first j such that
            #tableau_vector[j]<tableau_vector[j-1].
            #l will correspond to the shape of T_j
            l = [0]*size
            l[0] = 1
            j = 0
            for i in range(1,size):
                l[tableau_vector[i]] += 1
                if ( tableau_vector[i] < tableau_vector[i-1] ):
                    j = i
                    break

            #Find the last nonzero row of l and store it in k
            i = size - 1
            while ( l[i] == 0 ):
                i -= 1
            k = i

            #Find a new row for the letter j (next lowest corner)
            t = l[ 1 + tableau_vector[j] ]
            i = k
            while ( l[i] != t ):
                i -= 1

            #Move the letter j to row i
            tableau_vector[j] = i
            l[i] -= 1

            #Fill in the columns of T_j using 1,...,j-1 in increasing order
            m = 0
            while ( m < j ):
                r = 0
                while ( l[r] != 0 ):
                    tableau_vector[m] = r
                    l[r] -= 1
                    m += 1
                    r += 1

            #Convert the tableau vector back to the regular tableau
            #format
            row_count= [0]*len(pi)
            tableau = [[None]*n for n in pi]

            for i in range(size):
                tableau[tableau_vector[i]][row_count[tableau_vector[i]]] = i+1
                row_count[tableau_vector[i]] += 1

            yield self.element_class(self, tableau)

            #Check to see if we are at the last tableau
            #The last tableau if given by filling in the
            #partition along the rows.  For example, the
            #last partition corresponding to [3,2] is
            #[[1,2,3],
            # [4,5]]
            last_tableau = True
            i = 1
            for row in range(len(pi)):
                for col in range(pi[row]):
                    if tableau[row][col] != i:
                        last_tableau = False
                    i += 1

        return


    def list(self):
        r"""
        Returns a list of the standard Young tableaux of the specified shape.

        EXAMPLES::

            sage: StandardTableaux([2,2]).list()
            [[[1, 3], [2, 4]], [[1, 2], [3, 4]]]
            sage: StandardTableaux([5]).list()
            [[[1, 2, 3, 4, 5]]]
            sage: StandardTableaux([3,2,1]).list()
            [[[1, 4, 6], [2, 5], [3]],
             [[1, 3, 6], [2, 5], [4]],
             [[1, 2, 6], [3, 5], [4]],
             [[1, 3, 6], [2, 4], [5]],
             [[1, 2, 6], [3, 4], [5]],
             [[1, 4, 5], [2, 6], [3]],
             [[1, 3, 5], [2, 6], [4]],
             [[1, 2, 5], [3, 6], [4]],
             [[1, 3, 4], [2, 6], [5]],
             [[1, 2, 4], [3, 6], [5]],
             [[1, 2, 3], [4, 6], [5]],
             [[1, 3, 5], [2, 4], [6]],
             [[1, 2, 5], [3, 4], [6]],
             [[1, 3, 4], [2, 5], [6]],
             [[1, 2, 4], [3, 5], [6]],
             [[1, 2, 3], [4, 5], [6]]]
        """
        return [y for y in self]


    def random_element(self):
        """
        Returns a random standard tableau of the given shape using the
        Green-Nijenhuis-Wilf Algorithm.

        EXAMPLES::

            sage: StandardTableaux([2,2]).random_element()
            [[1, 2], [3, 4]]
        """

        p = self.shape

        t = [[None]*n for n in p]


        #Get the cells in the
        cells = []
        for i in range(len(p)):
            for j in range(p[i]):
                cells.append((i,j))

        m = sum(p)
        while m > 0:

            #Choose a cell at random
            cell = random.choice(cells)


            #Find a corner
            inner_corners = p.corners()
            while cell not in inner_corners:
                hooks = []
                for k in range(cell[1], p[cell[0]]):
                    hooks.append((cell[0], k))
                for k in range(cell[0], len(p)):
                    if p[k] > cell[1]:
                        hooks.append((k, cell[1]))

                cell = random.choice(hooks)


            #Assign m to cell
            t[cell[0]][cell[1]] = m

            p = p.remove_cell(cell[0])

            cells.remove(cell)

            m -= 1

        return self.element_class(self, t)


##########################
# Symmetric group action #
##########################
def unmatched_places(w, open, close):
    """
    EXAMPLES::

        sage: from sage.combinat.tableau import unmatched_places
        sage: unmatched_places([2,2,2,1,1,1],2,1)
        ([], [])
        sage: unmatched_places([1,1,1,2,2,2],2,1)
        ([0, 1, 2], [3, 4, 5])
        sage: unmatched_places([], 2, 1)
        ([], [])
        sage: unmatched_places([1,2,4,6,2,1,5,3],2,1)
        ([0], [1])
        sage: unmatched_places([2,2,1,2,4,6,2,1,5,3], 2, 1)
        ([], [0, 3])
        sage: unmatched_places([3,1,1,1,2,1,2], 2, 1)
        ([1, 2, 3], [6])
    """
    lw = len(w)
    places_open = []
    places_close = []
    for i in range(lw):
        letter = w[i]
        if letter == open:
            places_open.append(i)
        elif letter == close:
            if places_open == []:
                places_close.append(i)
            else:
                places_open.pop()
    return places_close, places_open


def symmetric_group_action_on_values(word, perm):
    """
    EXAMPLES::

        sage: from sage.combinat.tableau import symmetric_group_action_on_values
        sage: symmetric_group_action_on_values([1,1,1],[1,3,2])
        [1, 1, 1]
        sage: symmetric_group_action_on_values([1,1,1],[2,1,3])
        [2, 2, 2]
        sage: symmetric_group_action_on_values([1,2,1],[2,1,3])
        [2, 2, 1]
        sage: symmetric_group_action_on_values([2,2,2],[2,1,3])
        [1, 1, 1]
        sage: symmetric_group_action_on_values([2,1,2],[2,1,3])
        [2, 1, 1]
        sage: symmetric_group_action_on_values([2,2,3,1,1,2,2,3],[1,3,2])
        [2, 3, 3, 1, 1, 2, 3, 3]
        sage: symmetric_group_action_on_values([2,1,1],[2,1])
        [2, 1, 2]
        sage: symmetric_group_action_on_values([2,2,1],[2,1])
        [1, 2, 1]
        sage: symmetric_group_action_on_values([1,2,1],[2,1])
        [2, 2, 1]
    """
    w = list(word)
    ts = permutation.Permutation(perm).reduced_word()
    for j in reversed(range(len(ts))):
        r = ts[j]
        l = r + 1
        places_r, places_l = unmatched_places(w, l, r)

        #Now change the number of l's and r's in the new word
        nbl = len(places_l)
        nbr = len(places_r)
        ma = max(nbl, nbr)
        dif = ma - min(nbl, nbr)
        if ma == nbl:
            for i in range(dif):
                w[places_l[i]] = r
        else:
            for i in range(nbr-dif,ma):
                w[places_r[i]] = l
    return w

