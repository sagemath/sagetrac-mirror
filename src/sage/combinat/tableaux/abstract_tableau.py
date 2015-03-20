r"""
Abstract Tableau(x) classes.
"""

from sage.structure.element import Element
from sage.misc.abstract_method import abstract_method
from sage.misc.classcall_metaclass import ClasscallMetaclass

from collections import defaultdict

class AbstractTableau(Element):
    r"""
    Abstract class for the various element classes of Tableau.
    
    Tableau are thought of as a mapping from pairs (x, y) to some
    arbitrary labels. Two x-coordinates or two y-coordinates are assumed
    to be comparable via `>` and `<`. Subclasses are welcome to add further
    data (ex. shape).
    
    Tableau are assumed to be immutable; see `trac`:15862.
    """
    __metaclass__ = ClasscallMetaclass

    def __hash__(self):
        r""""
        Return the hash of ``self``.
        """
        return hash(tuple(self._dict_unsafe().iteritems()))

    def dict(self):
        r"""
        Return a dictionary which is a copy of the underlying data.
        
        Keys are assumed to be pairs. 
        """
        return {k:v for k, v in self._dict_unsafe().iteritems()}

    @abstract_method
    def _dict_unsafe(self):
        r"""
        Return a dictionary representing the underlying data.
        
        May or may not be a copy. Child classes will likely override.
        """
        pass

    def cells(self):
        r"""
        Return a list of the cells in this tableau.
        """
        return self._dict_unsafe().keys()

    def rows(self):
        r"""
        Return a list of rows of ``self`` in order of increasing row
        index, where each row is a list.
        """
        cols = defaultdict(list)
        _dict = self._dict_unsafe()
        for cell, label in _dict.iteritems():
            cols[cell[0]].append((cell[1], label))
        # Sort the columns
        for x, col in cols.iteritems():
            col_sorted = sorted(col, key=lambda v:v[0])
            cols[x] = [k[1] for k in col_sorted]
        # Sort the rows
        cols = sorted(((x, col) for x, col in cols.iteritems()),
                     key=lambda v: v[0])
        return [col[1] for col in cols]
            
    def cols(self):
        r"""
        Return a list of columns of ``self`` in order of increasing column
        index, where each column is a list.
        """
        rows = defaultdict(list)
        _dict = self._dict_unsafe()
        for cell, label in _dict.iteritems():
            rows[cell[1]].append((cell[0], label))
        # Sort the columns
        for x, row in rows.iteritems():
            row_sorted = sorted(row, key=lambda v:v[0])
            rows[x] = [k[1] for k in row_sorted]
        # Sort the rows
        rows = sorted(((x, row) for x, row in rows.iteritems()),
                      key=lambda v: v[0])
        return [row[1] for row in rows]

    def iter_by_row_reversed(self):
        r"""
        Iterate over entries of ``self``, row-by-row from higher row indices to lower,
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
        sorted_cells = sorted(_dict.keys(), cmp=mixed_lex)
        for cell in sorted_cells:
            yield _dict[cell]

    def iter_by_col_reversed(self):
        r"""
        Iterate over entries of ``self``, column-by-column from higher column indices to
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
        sorted_cells = sorted(_dict.keys(), cmp=mixed_lex)
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
        Return a Counter mapping values of ``self'' to multiplicities.
        """
        from collections import Counter
        return Counter(self._dict_unsafe().values())

    def size(self):
        r"""
        Return the number of cells in ``self``.
        """
        return len(self._dict_unsafe())

#    filter, restrict, append, reflect, transpose

class ExampleTableau(AbstractTableau):
    def __init__(self, d=None):
        if d is None:
            d = {(0,1): 3, (-3, 5): 17, (1,1): "a"}
        self._dict = d
        
    def _dict_unsafe(self):
        return self._dict

    def _repr_(self):
        return repr(self._dict)