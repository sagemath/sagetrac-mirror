r"""
AbstractTableau Element class.

See AbstractTableaux for the corresponding Parent class.

AUTHORS:

- Josh Swanson (and others) (2015): initial version
"""
#*****************************************************************************
#       Copyright (C) 2015 Josh Swanson,
#                     2015 Jan Keitel
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
from collections import defaultdict
from copy import deepcopy
from sage.structure.element import Element
from sage.misc.abstract_method import abstract_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.lazy_import import LazyImport

class AbstractTableau(Element):
    r"""
    Abstract class for the various element classes of tableaux.

    A tableau is thought of as a mapping which sends some pairs
    ``(x, y)`` (commonly, but not necessarily, pairs of nonnegative
    integers) to some arbitrary objects. Any two x-coordinates are
    assumed to be comparable via `>` and `<`, and so are any two
    y-coordinates.
    Subclasses are welcome to add further data (e.g., a skew shape).
    Tableaux are assumed to be immutable; see :trac:`15862`.

    The pairs ``(x, y)`` in the domain of a tableau are called its
    *cells*, and their images under the tableaux are referred to as
    its *entries* (in these cells).

    Subclasses must implement :meth:`_dict_unsafe`.
    """
    # See __classcall__.
    __metaclass__ = ClasscallMetaclass

    # "Default" parent class. See __classcall__.
    _generic_parent = LazyImport(
                      'sage.combinat.tableaux.abstract_tableaux',
                      'AbstractTableaux')

    def __init__(self, parent, *args, **kwds):
        r"""
        Initialize the AbstractTableau.

        Input validation should be done in parent classes, likely in
        their ``_element_constructor_`` method or coercions.
        Element class initialization should be quite minimal.

        We need to either call Element's __init__ method or
        set _parent by hand.
        """
        self._parent = parent

    @staticmethod
    def _gp(m, c):
        r"""
        Lazily import a parent class by name.

        Allows separating element and parent class files while also
        using __classcall__ shortcuts on element classes
        without circular references. One could simply call LazyImport
        directly in subclasses, but this way allows subclasses to
        ignore more details.
        """
        return LazyImport(m, c)

    @staticmethod
    def __classcall__(cls, *args, **kwds):
        r"""
        Provide shortcut syntax like ``Tableau([[1, 2], [3]])`` for
        ``Tableaux()([[1, 2], [3]])``.

        Input validation should be done in parent classes, likely in
        their ``_element_constructor_`` method or coercions.
        Element class initialization should be quite minimal.

        Inherited by child classes, unlike ``__classcall_private__``.
        """
        return cls._generic_parent()(*args, **kwds)

    def __hash__(self):
        r""""
        Return the hash of ``self`` based on the underlying mapping of
        cells to values.
        """
        return hash(tuple(six.iteritems(self._dict_unsafe())))

    def _repr_(self):
        r"""
        Return the representation string.

        OUTPUT:

        A string.
        """
        return repr(self._dict_unsafe())

    def dict(self):
        r"""
        Return a dictionary which is a copy of the underlying data.

        OUTPUT:

        A dictionary. The value of this dictionary at the key
        ``(x, y)`` is the entry of the tableau ``self`` in the cell
        ``(x, y)``.
        """
        return deepcopy(self._dict_unsafe())

    @abstract_method
    def _dict_unsafe(self):
        r"""
        Return a dictionary representing the underlying data.

        It is unsafe to alter this dictionary or pass it along, since
        this might (and might not) mutate ``self``. That said, it is
        also unsafe to try and mutate ``self`` by altering this
        dictionary; this will often fail or result in broken objects.
        The main use of ``_dict_unsafe`` should be getting entries
        of ``self``.
        """
        pass

    def cells(self):
        r"""
        Return an iterable over the cells in this tableau.
        """
        return six.iterkeys(self._dict_unsafe())

    def _group_by(self, a, b):
        r"""
        Return a list of lists given by grouping the values of ``self``.

        Inner lists come from fibers of projection onto index `a` of keys,
        sorted by index `b`. The outer list is sorted by index `a`.
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
        by column index.
        """
        return self._group_by(0, 1)

    def cols(self):
        r"""
        Return the values of ``self`` as a list of columns in order of
        increasing columns index, where each column is a list and has been
        sorted by row index.
        """
        return self._group_by(1, 0)

    def iter_by_row_reversed(self):
        r"""
        Iterate over values of ``self``, row-by-row from higher row indices to lower,
        and within each row from lower column indexes to higher.
        """
        # Row index is more important than column index; higher row indices
        #   are earlier; lower column indices are earlier
        def mixed_lex((x1, y1), (x2, y2)):
            if x1 < x2:
                return 1
            if x1 > x2:
                return -1
            # Assume x1 == x2
            if y1 < y2:
                return -1
            if y1 > y2:
                return 1
            # Assume y1 == y2
            return 0

        _dict = self._dict_unsafe()
        # TODO: replace cmp with key for Python 3 compatibility
        sorted_cells = sorted(six.iterkeys(_dict), cmp=mixed_lex)
        for cell in sorted_cells:
            yield _dict[cell]

    def iter_by_col_reversed(self):
        r"""
        Iterate over values of ``self``, column-by-column from higher column indices to
        lower, and within each column from lower row indexes to higher.
        """
        # Column index is more important than row index; higher column indices
        #   are earlier; lower row indices are earlier
        def mixed_lex((x1, y1), (x2, y2)):
            if y1 < y2:
                return 1
            if y1 > y2:
                return -1
            # Assume y1 == y2
            if x1 < x2:
                return -1
            if x1 > x2:
                return 1
            # Assume x1 == x2
            return 0

        _dict = self._dict_unsafe()
        # TODO: replace cmp with key for Python 3 compatibility
        sorted_cells = sorted(six.iterkeys(_dict), cmp=mixed_lex)
        for cell in sorted_cells:
            yield _dict[cell]

    def word_by_row(self):
        r"""
        Return a flattened version of ``self.iter_by_row'' as a
        :class:`sage.combinat.words.word.FiniteWord_list`.
        """
        from sage.combinat.words.word import Word
        return Word(sum(self.iter_by_row_reversed(), []))

    def word_by_col(self):
        r"""
        Return a flattened version of ``self.iter_by_col'' as a
        :class:`sage.combinat.words.word.FiniteWord_list`.
        """
        from sage.combinat.words.word import Word
        return Word(sum(self.iter_by_col_reversed(), []))

    # Alias
    word = word_by_row

    def weight_dict(self):
        r"""
        Return a Counter mapping values of ``self`` to multiplicities.
        """
        from collections import Counter
        return Counter(six.itervalues(self._dict_unsafe()))

    def size(self):
        r"""
        Return the number of cells in ``self``.
        """
        return len(self._dict_unsafe())

    # Alias
    __len__ = size
