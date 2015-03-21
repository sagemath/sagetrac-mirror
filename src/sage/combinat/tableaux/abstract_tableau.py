r"""
Abstract Tableau(x) classes.
"""

from sage.structure.element import Element
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.integer_ring import ZZ
from sage.categories.sets_cat import Sets
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from collections import defaultdict

class AbstractTableau(Element):
    r"""
    Abstract class for the various element classes of Tableau.
    
    A Tableau is thought of as a mapping from pairs (x, y) to some
    arbitrary labels. Two x-coordinates or two y-coordinates are assumed
    to be comparable via `>` and `<`. Subclasses are welcome to add further
    data (ex. shape).
    
    Tableau are assumed to be immutable; see :trac:`15862`.
    
    This is an abstract class. Child classes must implement
    :meth:`_dict_unsafe`.
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

class BadShapeTableau(AbstractTableau):
    r"""
    Tableaux of bad shape
    
    A BadShapeTableau is specified by a dictionary
    whose keys are pairs of integers.
    """
    @staticmethod
    def __classcall_private__(cls, dictionary, check=True):
        r"""
        Return the BadShapeTableau defined by ``dictionary``.
        """
        if isinstance(dictionary, cls):
            return dictionary
        
        return BadShapeTableaux()(dictionary, check)
    
    def __init__(self, parent, dictionary, check=True):
        r"""
        Initialize the BadShapeTableau.
        
        INPUT:
        
        - ``dictionary`` -- a dictionary whose keys are assumed to be pairs
          of integers
        - ``check`` -- (default: ``True``) if ``True``, then check that
          the keys of ``dict`` are in fact pairs of integers
        """
        if check and not all(i in ZZ and j in ZZ for i, j in
                             dictionary.keys()):
            raise ValueError('keys must be pairs of integers')
        self._dict = dictionary
        Element.__init__(self, parent)
        
    def _dict_unsafe(self):
        r"""
        Return the dictionary containing the defining data of ``self``.
        
        OUTPUT:
        
        A dictionary.
        """
        return self._dict

    def _repr_(self):
        r"""
        Return the representation string.
        
        OUTPUT:
        
        A string.
        """
        return repr(self._dict_unsafe())
    
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
        data = {i:j for i, j in self._dict_unsafe() if predicate(*i)}
        return self.__class__(dictionary=data, check=False)
    
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
        data = {i:j for i, j in self._dict_unsafe() if predicate(j)}
        return self.__class__(dictionary=data, check=False)

    def transpose(self):
        r"""
        Return the transposition of ``self``.
        
        OUTPUT:
        
        A BadShapeTableau.
        """
        data = {(i[1], i[0]):j for i, j in self._dict_unsafe()}
        return self.__class__(dictionary=data, check=False)



class SkewTableau(BadShapeTableau):
    
    @staticmethod
    def __classcall_private__(cls, l=[], dictionary=[], check=False):
        r"""
        Return the SkewTableau corresponding to the defining data.
        """
        if isinstance(l, cls):
            return l
        return SkewTableaux()(l, dictionary, check)

    def __init__(self, parent, l=[], dictionary={}, check=False):
        r"""
        Initialize the SkewTableau.
        
        Input can be specified either in terms of a
        dictionary or a list. If both are specified,
        then the dictionary data is skipped.
        """
        if not l and dictionary:
            rows = defaultdict(dict)
            min_row_index = None
            max_row_index = None
            min_col_index = None
            for i, j in dictionary.iteritems():
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
                for col_index, value in row.iteritems():
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









class AbstractTableaux(UniqueRepresentation, Parent):
    Element = AbstractTableau

    def __init__(self):
        r"""
        Initialize the parent.
        """
        Parent.__init__(self, category=Sets())
                        
    def _repr_(self):
        r"""
        Return the representation string.
        
        OUTPUT:
        
        A string.
        """
        return "Abstract Tableaux"
                        

                        
class BadShapeTableaux(AbstractTableaux):
    Element = BadShapeTableau
    
    def _repr_(self):
        r"""
        Return the representation string.
        
        OUTPUT:
        
        A string.
        """
        return "Bad Shape Tableaux"

class SkewTableaux(BadShapeTableaux):
    Element = SkewTableau
    
    def _repr_(self):
        r"""
            Return the representation string.
            
            OUTPUT:
            
            A string.
            """
        return "Skew Tableaux"


class StraightTableaux(SkewTableaux):
    Element = SkewTableau
    
    def _repr_(self):
        r"""
            Return the representation string.
            
            OUTPUT:
            
            A string.
            """
        return "Straight Tableaux"