r"""
AbstractTableau Element class.

See `class:AbstractTableaux` for the corresponding Parent class.

AUTHORS:

- Josh Swanson (and others) (2015): initial version
"""
#*****************************************************************************
#       Copyright (C) 2015 Jan Keitel
#                     2015 Josh Swanson
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

import six
from collections               import defaultdict
from sage.structure.element    import Element
from sage.misc.abstract_method import abstract_method

class AbstractTableau(Element):
    r"""
    Abstract class for the various Element classes of tableaux.

    A tableau is thought of as a mapping which sends some pairs
    ``(x, y)`` (commonly, but not necessarily, pairs of nonnegative
    integers) to some arbitrary objects. Any two x-coordinates are
    assumed to be comparable via `>` and `<`, and so are any two
    y-coordinates.

    Subclasses are welcome to add further data (e.g., a skew shape).
    Tableaux are assumed to be immutable; see :trac:`15862`.

    The pairs ``(x, y)`` in the domain of a tableau are called its
    *cells*, and their images under the tableaux are referred to as
    the *entries* in those cells.

    Subclasses must implement :meth:`_dict_unsafe`.
    """
    def __init__(self, parent, *args, **kwds):
        r"""
        Initialize the `class:AbstractTableau`.

        Input normalization and validation should be done in parent classes,
        likely in their ``_element_constructor_`` method or coercions.
        Element class initialization should be quite minimal.

        We need to either call Element's ``__init__`` method or
        set `_parent` by hand.
        
        TESTS::

            sage: SkewTableau([[None, None, 1, 3], [None, 4, 4], [None]]).parent()
            Skew tableaux 
        """
        self._parent = parent

    def __hash__(self):
        r""""
        Return the hash of ``self`` based on the underlying mapping of
        cells to values.
        
        TESTS::

            sage: b1 = BadShapeTableau({(2, 2): 1, (0, 0): 0})
            sage: b2 = BadShapeTableau({(2, 2): 2, (0, 0): 0})
            sage: b3 = BadShapeTableau({(2, 1): 1, (0, 0): 0})
            sage: b1.__hash__() == b1.__hash__()
            True
            sage: b1.__hash__() != b2.__hash__()
            True
            sage: b2.__hash__() != b3.__hash__()
            True
            sage: b3.__hash__() != b1.__hash__()
            True
        """
        return hash(tuple(six.iteritems(self._dict_unsafe())))

    def _richcmp_(self, other, op):
        r"""
        Provide rich comparison.

        Equality and inequality are tested by comparing the underlying `_dict_unsafe`.

        TESTS::

            sage: b1 = BadShapeTableau({(1, 1): 1, (3, 4): 5})
            sage: b2 = BadShapeTableau({(1, 1): 2, (3, 4): 5})
            sage: b1 == b2 # indirect doctest
            False
            sage: b1 != b2
            True
        """
        # TODO: find a way to get these magic constants from outside of Cython
        if op == 2:
            return self._dict_unsafe() == other._dict_unsafe()
        elif op == 3:
            return self._dict_unsafe() != other._dict_unsafe()
        raise TypeError("unorderable types")

    def _repr_(self):
        r"""
        Return the representation string.

        TESTS::

            sage: BadShapeTableau({(-1, -2): 'a'})
            {(-1, -2): 'a'}
        """
        return repr(self._dict_unsafe())

    def dict(self):
        r"""
        Return a dictionary which is a copy of the underlying data.

        Values are not deeply copied. For instance, if entries are lists,
        one could change the underlying tableau by modifying those lists.
        Tableaux are assumed to be immutable, so this is okay.

        OUTPUT:

        A dictionary. The value of this dictionary at the key
        ``(x, y)`` is the entry of the tableau ``self`` in the cell
        ``(x, y)``.

        TEST::

            sage: st = BadShapeTableau({(4, 3): 'a', (1, 2): 3})
            sage: st.dict() == {(4, 3): 'a', (1, 2): 3}
            True
        """
        return dict(self._dict_unsafe())

    @abstract_method
    def _dict_unsafe(self):
        r"""
        Return a dictionary representing the underlying data.

        It is unsafe to alter this dictionary or pass it along, since
        this might or might not mutate ``self`` and tableaux are assumed
        to be immutable.
        """
        pass

    def cells(self):
        r"""
        Return an iterable over the cells in this tableau in no particular order.

        TESTS::

            sage: b = BadShapeTableau({(4, -1): -2, (3, 3): 'cow'})
            sage: set(b.cells()) == set(((4, -1), (3, 3)))
            True
        """
        return six.iterkeys(self._dict_unsafe())

    def _group_by(self, a, b):
        r"""
        Return a list of lists given by grouping the values of ``self``.

        Inner lists are the fibers of projection onto index `a` of keys,
        sorted by index `b`. The outer list is sorted by index `a`.

        TESTS::

            sage: st = SkewTableau([[None, 1, 2], [None, 3], [None], [1]])
            sage: st._group_by(0, 1)
            [[1, 2], [3], [1]]
            sage: st._group_by(1, 0)
            [[1], [1, 3], [2]]
        """
        # We imagine a=0 and b=1, meaning we group into rows.
        rows = defaultdict(list)
        _dict = self._dict_unsafe()
        for cell, label in six.iteritems(_dict):
            rows[cell[a]].append((cell[b], label))
        # Sort inside each row by column index
        for c, row in six.iteritems(rows):
            row_sorted = sorted(row, key=lambda v: v[0])
            rows[c] = [k[1] for k in row_sorted]
        # Sort rows themselves
        rows = sorted(six.iteritems(rows), key=lambda v: v[0])
        return [row[1] for row in rows]

    def rows(self):
        r"""
        Return the values of ``self`` as a list of rows in order of
        increasing row index, where each row is a list and has been sorted
        by increasing column index.

        TESTS::
        
            sage: s = SkewTableau([[None, 1, 1, 2, 4], [None, 2, 3, 3, 5], [1, 3, 5, 7]])
            sage: s.rows()
            [[1, 1, 2, 4], [2, 3, 3, 5], [1, 3, 5, 7]]
        """
        return self._group_by(0, 1)

    def columns(self):
        r"""
        Return the values of ``self`` as a list of columns in order of
        increasing columns index, where each column is a list and has been
        sorted by increasing row index.

        TESTS::

            sage: s = SkewTableau([[None, 1, 1, 2, 4], [None, 2, 3, 3, 5], [1, 3, 5, 7]])
            sage: s.columns()
            [[1], [1, 2, 3], [1, 3, 5], [2, 3, 7], [4, 5]]
        """
        return self._group_by(1, 0)

    def iter_cells(self):
        r"""
        Iterate over the cells of ``self`` in no particular order.

        TESTS::

            sage: s = SkewTableau([[None, 5, 5, 6], [None, 1, 2, 3]])
            sage: set(s.iter_cells()) == {(0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (1, 3)}
            True
        """
        return six.iterkeys(self._dict_unsafe())

    def iter_entries(self):
        r"""
        Iterate over the entries of ``self`` in no particular order.

        TESTS::

            set(SkewTableau([[None, 7, 8, 9], [None, 1, 2, 4]]).iter_entries()) == {1, 2, 4, 7, 8, 9}
            True
        """
        return six.itervalues(self._dict_unsafe())

    def cells_by_content(self, c):
        r"""
        Return a list of the coordinates of the cells in ``self`` with content ``c``.

        The content of a cell `(a, b)` is defined as `b-a`.

        EXAMPLES::

            sage: s = SkewTableau([[None,1,2],[3,4,5],[6]])
            sage: s.cells_by_content(0)
            [(1, 1)]
            sage: s.cells_by_content(1)
            [(0, 1), (1, 2)]
            sage: s.cells_by_content(2)
            [(0, 2)]
            sage: s.cells_by_content(-1)
            [(1, 0)]
            sage: s.cells_by_content(-2)
            [(2, 0)]
        """
        return [(a, b) for (a, b) in self.iter_cells()
                                  if b-a == c]

    def entries_by_content(self, c):
        r"""
        Return a list of the entries in ``self`` with content ``c``.

        The content of a cell `(a, b)` is defined as `b-a`.

        EXAMPLES::

            sage: s = SkewTableau([[None,1,2],[3,4,5],[6]])
            sage: s.entries_by_content(0)
            [4]
            sage: s.entries_by_content(1)
            [1, 5]
            sage: s.entries_by_content(2)
            [2]
            sage: s.entries_by_content(-1)
            [3]
            sage: s.entries_by_content(-2)
            [6]
        """
        return [v for (a,b), v in six.iteritems(self._dict_unsafe())
                               if b-a==c]

    def cells(self):
        r"""
        Return a list of cells in ``self`` in no particular order.

        EXAMPLES::

            sage: s = SkewTableau([[None,1,2],[3],[6]])
            sage: set(s.cells()) == set([(0, 1), (0, 2), (1, 0), (2, 0)])
            True
        """
        return list(self.iter_cells())

    def cells_containing(self, i):
        r"""
        Return the list of cells which have entry ``i``.

        The list is ordered by column index, with no guarantee on row
        index ordering.

        EXAMPLES::

            sage: t = SkewTableau([[None,None,3],[None,3,5],[4,5]])
            sage: t.cells_containing(5)
            [(2, 1), (1, 2)]
            sage: t.cells_containing(4)
            [(2, 0)]
            sage: t.cells_containing(2)
            []

            sage: t = SkewTableau([[None,None,None,None],[None,4,5],[None,5,6],[None,9],[None]])
            sage: t.cells_containing(2)
            []
            sage: t.cells_containing(4)
            [(1, 1)]
            sage: t.cells_containing(5)
            [(2, 1), (1, 2)]

            sage: SkewTableau([]).cells_containing(3)
            []

            sage: SkewTableau([[None,None],[None]]).cells_containing(3)
            []
        """

        return sorted([c for c, v
                       in six.iteritems(self._dict_unsafe())
                       if v==i], key=lambda c: c[1])

    def to_word_by_row(self):
        r"""
        Return the row reading word.

        More precisely, return a :class:`~sage.combinat.words.word.FiniteWord_list`
        obtained from reading entries row-by-row, from the topmost row to
        the bottommost, in increasing order of column index within each row.

        TESTS::

            sage: s = SkewTableau([[None, 'fruit', 1], ['fly', 0, 0]])
            sage: s.to_word_by_row()
            word: fly,0,0,fruit,1

            sage: Tableau([[1,2],[3,4]]).to_word_by_row()
            word: 3412
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word_by_row()
            word: 325146

            sage: Tableau([[1,2],[3,4]]).to_word()
            word: 3412
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word()
            word: 325146
        """
        from sage.combinat.words.word import Word
        return Word(v for x in reversed(self.rows()) for v in x)

    def to_word_by_column(self):
        r"""
        Return the column reading word.

        More precisely, return a :class:`~sage.combinat.words.word.FiniteWord_list`
        obtained from reading entries row-by-row, from the topmost row to
        the bottommost, in increasing order of column index within each row.

        TESTS::

            sage: s = SkewTableau([[None, 'fruit', 1], ['fly', 0, 0]])
            sage: s.to_word_by_column()
            word: 1,0,fruit,0,fly

            sage: Tableau([[1,2],[3,4]]).to_word_by_column()
            word: 3142
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word_by_column()
            word: 321546
        """
        from sage.combinat.words.word import Word
        return Word(v for x in reversed(self.columns()) for v in x)

    # Alias
    to_word = to_word_by_row

    def weight_counter(self):
        r"""
        Return a `class:Counter` mapping values of ``self`` to multiplicities.

        TESTS::

            sage: from collections import Counter
            sage: s = SkewTableau([[None, None, 2, 2, 3, 4], [None, 1, 2, 3]])
            sage: s.weight_counter() == Counter({2: 3, 3: 2, 1: 1, 4: 1})
            True
        """
        from collections import Counter
        return Counter(self.iter_entries())

    def size(self):
        r"""
        Return the number of cells in ``self``.

        TESTS::

            sage: b = BadShapeTableau({(3, 3): 2, (2, 2): 1, (1, 1): -1})
            sage: b.size()
            3
            sage: SkewTableau([[None, 2, 4], [None, 3], [1]]).size()
            4
            sage: SkewTableau([[None, 2], [1, 3]]).size()
            3
        """
        return len(self._dict_unsafe())

    # Alias
    __len__ = size
