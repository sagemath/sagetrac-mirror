r"""
Basic abstract tableau element classes.
"""

import six
from collections import defaultdict
from copy import deepcopy
from sage.structure.element import Element
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.lazy_import import LazyImport


def parent_class(s):
    r"""
    Lazily import a parent class by name.
    
    Allows separating element and parent class files without circular
    references.
    """
    return LazyImport('sage.combinat.tableaux.abstract_tableaux', s)

class AbstractTableau(six.with_metaclass(ClasscallMetaclass, Element)):
    r"""
    Abstract class for the various element classes of Tableau.

    A Tableau is thought of as a mapping from pairs (x, y) to some
    arbitrary labels. Two x-coordinates or two y-coordinates are assumed
    to be comparable via `>` and `<`. Subclasses are welcome to add further
    data (ex. skew shape). Tableau are assumed to be immutable; see
    :trac:`15862`.

    Child classes must implement :meth:`_dict_unsafe`.
    """
    # "Default" parent class. See __classcall__.
    _generic_parent = parent_class('AbstractTableaux')

    @staticmethod
    def __classcall__(cls, *args, **kwds):
        r"""
        Provide shortcut syntax like ``Tableau([[1, 2], [3]])`` for
        ``Tableaux()([[1, 2], [3]])``.
        
        Input validation should be done in parent classes, likely in
        their own __classcall__ methods. Element class initialization
        should be quite minimal.
        
        Inherited by child classes, unlike __classcall_private__.
        """
        return cls._generic_parent(*args, **kwds)

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
        
        Keys are assumed to be pairs. 
        """
        return deepcopy(self._dict_unsafe())

    @abstract_method
    def _dict_unsafe(self):
        r"""
        Return a dictionary representing the underlying data.
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
        Return a Counter mapping values of ``self'' to multiplicities.
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


class BadShapeTableau(AbstractTableau):
    r"""
    A tableau of bad shape.
    
    A BadShapeTableau is specified by a dictionary
    whose keys are pairs of integers and whose values are
    arbitrary.
    """
    _generic_parent = parent_class('BadShapeTableaux')

    def __init__(self, parent, dct, check=True):
        r"""
        Initialize the BadShapeTableau.
        
        INPUT:
        
        - ``dct`` -- a dictionary (or more generally something
          passable to ``dict``) whose keys are pairs of integers
        - ``check`` -- (default: ``True``) if ``True``, then check that
          the keys of ``dct`` are in fact pairs of integers
        """
        self._dict = deepcopy(dict(dct))
        
        if check:
            if not all(x in ZZ and y in ZZ for x, y
                       in six.iterkeys(self._dict)):
                raise ValueError('keys must be pairs of integers')

        AbstractTableau.__init__(self, parent)

    def _dict_unsafe(self):
        r"""
        Return the dictionary containing the defining data of ``self``.
        
        OUTPUT:
        
        A dictionary.
        """
        return self._dict

    def filter_by_cells(self, predicate):
        r"""
        Return a tableau whose cells satisfy the given
        ``predicate``.
            
        INPUT:
            
        - ``predicate`` -- a function accepting two parameters and returning
        ``True`` or ``False``
            
        OUTPUT:
            
        A BadShapeTableau
        """
        data = {k:v for k, v in six.iteritems(self._dict_unsafe()) 
                if predicate(*k)}
        return self.__class__(dictionary=data, check=True)

    def filter_by_values(self, predicate):
        r"""
        Return a tableau whose values satisfy the given
        ``predicate``.
        
        INPUT:
        
        - ``predicate`` -- a function accepting one parameter and returning
          ``True`` or ``False``
          
        OUTPUT:
        
        A BadShapeTableau
        """
        data = {k:v for k, v in six.iteritems(self._dict_unsafe()) 
                if predicate(v)}
        return self.__class__(dictionary=data, check=True)

    def transpose(self):
        r"""
        Return the transposition of ``self``.
        
        OUTPUT:
        
        A BadShapeTableau.
        """
        data = {(k[1], k[0]):v for k, v in self._dict_unsafe()}
        return self.__class__(dictionary=data, check=True)


class SkewTableau(BadShapeTableau):
    r"""
    A tableau of skew shape.

    Cell locations are pairs of positive integers, though values are
    unrestricted. Cells must form a skew shape.
    """
    def __init__(self, parent, inner, rows):
        r"""
        Initialize the SkewTableau.

        INPUT:
        
        - ``inner`` -- an iterable of weakly decreasing non-negative
          integers
        - ``check`` -- (default: ``True``) if ``True``, then ensure
          keys are positive integers and the shape is a skew shape.

        Input can be specified either in terms of a
        dictionary or a list. If both are specified,
        then the dictionary data is skipped.
        """
        if not l and dictionary:
            rows = defaultdict(dict)
            min_row_index = None
            max_row_index = None
            min_col_index = None
            for i, j in six.iteritems(dictionary):
                rows[i[0]][i[1]] = j
                if (min_row_index is None) or (min_row_index > i[0]):
                    min_row_index = i[0]
                if (max_row_index is None) or (min_row_index < i[0]):
                    max_row_index = i[0]
                if (min_col_index is None) or (min_col_index > i[1]):
                    min_col_index = i[1]
            # Create list
            l = []
            for row_index in range(min_row_index, max_row_index + 1):
                row = rows[row_index]
                if not row:
                    l.append([])
                    continue
                tmp = [None]*(max(row.keys()) - min_col_index + 1)
                for col_index, value in six.iteritems(row):
                    tmp[col_index - min_col_index] = value
                l.append(tmp)

        # We don't have to do the checking of BadShapeTableau, that
        # happens above anyway, since the cells are used as list indices
        BadShapeTableau.__init__(self, parent, dictionary={}, check=False)
        # Check
        if check:
            raise NotImplementedError
        self._list = tuple(map(tuple, l))

    @cached_method
    def _dict_unsafe(self):
        d = {(i, j): k for i, row in enumerate(self._tuple())
            for j, k in enumerate(row)}
        return d

    def _tuple(self):
        r"""
        Return the tuple of tuples that defines ``self``.
        """
        return self._list

    def _repr_(self):
        r"""
        Return the representation string.
        
        OUTPUT:
        
        A string.
        """
        return repr(self.to_list())

    def to_list(self):
        r"""
        Return the list of lists that defines ``self``.
        
        This is a copy of the defining data.
        
        OUTPUT:
        
        A list of lists.
        """
        return [list(i) for i in self._tuple()]


class StraightTableau(SkewTableau):
    def __init__(self, parent, l=[], dictionary=[], check=False):
        SkewTableau.__init__(self, parent, l, dictionary, check)









class ExampleTableau(BadShapeTableau):
    def __init__(self):
        BadShapeTableau.__init__(self, d = {(0,1): 3, (-3, 5): 17, (1,1): "a"})



