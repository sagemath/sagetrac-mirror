r"""
Abstract Tableau(x) classes.
"""
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
        _dict = self._dict_unsafe()
        _cells = _dict.keys()
        row_indexes = sorted(set(x for x, _ in cells))
        row_trans = {x:i for i, x in enumerate(row_indexes)}
        rows = [[] for i in range(len(row_indexes))]
        for x, y in _cells:
            rows[row_trans[x]].append(_dict[(x, y)])
        return rows

    def cols(self):
        r"""
        Return a list of columns of ``self`` in order of increasing column
        index, where each column is a list.
        """
        _dict = self._dict_unsafe()
        _cells = _dict.keys()
        col_indexes = sorted(set(y for _, y in cells))
        col_trans = {y:i for i, x in enumerate(col_indexes)}
        cols = [[] for i in range(len(col_indexes))]
        for x, y in _cells:
            cols[col_trans[y]].append(_dict[(x, y)])
        return cols

    def iter_by_row(self):
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

    def iter_by_col(self):
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
        Return a flattened version of ``self.iter_by_row'' as a ``Word``.
        """
        from sage.combinat.words.word import Word
        return Word(sum(self.iter_by_row(), []))

    def word_by_col(self):
        r"""
        Return a flattened version of ``self.iter_by_col'' as a ``Word``.
        """
        from sage.combinat.words.word import Word
        return Word(sum(self.iter_by_col(), []))

    # Alias
    word = word_by_row

#    filter, restrict, append, size, weight, reflect, transpose

class ExampleTableau(AbstractTableau):
    def __init__(self, d=None):
        if d is None:
            d = {(0,1): 3, (-3, 5): 17, (1,1): "a"}
        self._dict = d
        
    def _dict_unsafe(self):
        return self._dict