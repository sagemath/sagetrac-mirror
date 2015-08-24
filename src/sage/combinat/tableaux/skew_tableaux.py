r"""
SkewTableax Parent class and subclasses.

See SkewTableau and its subclasses for the corresponding Element classes.

AUTHORS:

- Mike Hansen: Initial version
- Travis Scrimshaw, Arthur Lubovsky (2013-02-11): Factored out
  ``CombinatorialClass``
- Josh Swanson (and others) (2015): tableaux refactoring/cleanup

CLASSES:

.. TODO:: List classes.

Orphaned tests::

    sage: T = SkewTableau([[None, None, 1], [3], [4]])
    sage: T in SkewTableaux()
    True
    sage: [[None,1],[2,3]] in SkewTableaux()
    True

    sage: [[None, 2], [1, 3]] in StandardSkewTableaux()
    True
    sage: [[None, 2], [2, 4]] in StandardSkewTableaux()
    False
    sage: [[None, 3], [2, 4]] in StandardSkewTableaux()
    False
    sage: [[None, 2], [1, 4]] in StandardSkewTableaux()
    False

    sage: [[None, 2], [1, 3]] in SemistandardSkewTableaux()
    True
    sage: [[None, 2], [2, 4]] in SemistandardSkewTableaux()
    True
    sage: [[None, 3], [2, 4]] in SemistandardSkewTableaux()
    True
    sage: [[None, 2], [2, 4]] in SemistandardSkewTableaux()
    True
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#                     2013 Arthur Lubovsky
#                     2015 Josh Swanson
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
from itertools import ifilter
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.skew_partition import SkewPartition, SkewPartitions
from sage.matrix.all import zero_matrix
from sage.rings.all import Integer, ZZ, QQ
from sage.rings.arith import factorial
from sage.rings.infinity import PlusInfinity
from sage.structure.global_options import GlobalOptions

# TODO: rethink imports here vs. in-function, and likewise in skew_tableau.py

from sage.combinat.tableaux.bad_shape_tableaux import BadShapeTableaux
from sage.combinat.tableaux.skew_tableau import (
     SkewTableau, SemistandardSkewTableau, StandardSkewTableau)

TableauOptions=GlobalOptions(name='skew tableaux',
    doc=r"""
    Sets the global options for elements of the SkewTableau class.
    By default, they are displayed as a list, LaTeXed as a Young
    diagram, and English convention is used.
    """,
    end_doc=r"""

    .. NOTE::

        Changing the ``convention`` for tableaux also changes the
        ``convention`` for partitions.

    TODO: GlobalOptions docs are not currently doctested. Manually make sure
    these tests (and possibly the Partitions tests?) work after refactoring
    has been completed.

    If no parameters are set, then the function returns a copy of the
    options dictionary.

    EXAMPLES::

        sage: T = Tableau([[1,2,3],[4,5]])
        sage: T
        [[1, 2, 3], [4, 5]]
        sage: Tableaux.global_options(display="array")
        sage: T
          1  2  3
          4  5
        sage: Tableaux.global_options(convention="french")
        sage: T
          4  5
          1  2  3

    Changing the ``convention`` for tableaux also changes the ``convention``
    for partitions and vice versa::

        sage: P = Partition([3,3,1])
        sage: print P.ferrers_diagram()
        *
        ***
        ***
        sage: Partitions.global_options(convention="english")
        sage: print P.ferrers_diagram()
        ***
        ***
        *
        sage: T
          1  2  3
          4  5

    The ASCII art can also be changed::

        sage: t = Tableau([[1,2,3],[4,5]])
        sage: ascii_art(t)
          1  2  3
          4  5
        sage: Tableaux.global_options(ascii_art="table")
        sage: ascii_art(t)
        +---+---+
        | 4 | 5 |
        +---+---+---+
        | 1 | 2 | 3 |
        +---+---+---+
        sage: Tableaux.global_options(ascii_art="compact")
        sage: ascii_art(t)
        |4|5|
        |1|2|3|
        sage: Tableaux.global_options.reset()
    """,
    display=dict(default="list",
                 description='Controls the way in which tableaux are printed',
                 values=dict(list='print tableaux as lists',
                             diagram='display as Young diagram (similar to :meth:`~sage.combinat.tableau.Tableau.pp()`',
                             compact='minimal length string representation'),
                 alias=dict(array="diagram", ferrers_diagram="diagram", young_diagram="diagram"),
                 case_sensitive=False),
    ascii_art=dict(default="repr",
                 description='Controls the ascii art output for tableaux',
                 values=dict(repr='display using the diagram string representation',
                             table='display as a table',
                             compact='minimal length ascii art'),
                 case_sensitive=False),
    latex=dict(default="diagram",
               description='Controls the way in which tableaux are latexed',
               values=dict(list='as a list', diagram='as a Young diagram'),
               alias=dict(array="diagram", ferrers_diagram="diagram", young_diagram="diagram"),
               case_sensitive=False),
    convention=dict(default="English",
                    description='Sets the convention used for displaying tableaux and partitions',
                    values=dict(English='use the English convention',French='use the French convention'),
                    case_sensitive=False),
    notation = dict(alt_name="convention")
)

class SkewTableaux(BadShapeTableaux):
    r"""
    Parent class of all skew tableaux. See `meth:_element_constructor_` for
    element class construction options and examples.
    
    See SkewTableau for the Element class.
    """
    Element = SkewTableau

    # Convenient shortcut to global options
    global_options = TableauOptions

    def _element_constructor_(self, st=0, expr=None, shape_word=None,
                                    dct=None, check=True):
        r"""
        Construct a new SkewTableau by converting from one of several input
        formats, optionally validating and normalizing input.

        These are conversions. If multiple formats are specified, the
        left-most is used. If no format is specified, the "trivial" skew
        tableau with no entries is returned.

        INPUT:
        - ``x`` -- used in :meth:`Parent.__call__` and can be ignored here.
        - ``st`` -- an iterable of rows from top to bottom in English
          notation, where each row is an iterable of entries from left
          to right but where ``None``'s do not give cells of the skew
          tableau
        - ``expr`` -- a pair (``inner``, ``rows``) where ``inner`` is
          the inner partition and ``rows`` is an iterable of rows from
          bottom to top in English notation, where each row is an iterable
          of the entries in that row from left to right. Provided for
          compatibility with MuPAD-Combinat.
        - ``shape_word'' -- a pair (``shape``, ``word``) where ``shape``
          is a skew partition and the word ``word`` is obtained from the
          row reading
        - ``dct`` -- a dictionary whose keys are pairs of non-negative
          integers
        - ``check`` -- (default: ``True``) if ``True``, then validate 
          and normalize input: ensure ``st`` or the cells of ``dct``
          actually form a skew shape, remove empty rows from ``st``,
          etc.
        """
        # TODO: normalize st input by removing empty rows
        if st is not 0:
            return self.from_st(st, check)

        if expr is not None:
            return self.from_expr(expr, check)

        if shape_word is not None:
            shape, word = shape_word
            return self.from_shape_and_word(shape, word, check)

        if dct is not None:
            return self.from_dict(dct, check)

        return self._new_element(())

    def from_st(self, st, check=True):
        try:
            # Remove empty rows
            st = tuple( ifilter(lambda t: t,
                        (tuple(row) for row in st)) )
        except TypeError:
            raise TypeError("each element of the skew tableau must be an iterable")

        if check:
            # Make sure rows begin with blocks of nones
            inner = [0]*len(st)
            outer = [0]*len(st)
            for i in range(len(st)):
                row = st[i]

                finishedNones = False
                noneCount = 0
                for j, val in enumerate(row):
                    if val is None:
                        if finishedNones:
                            raise ValueError("Row %s does not start with a block of Nones"%str(row))
                        noneCount += 1
                    else:
                        finishedNones = True

                inner[i] = noneCount
                outer[i] = len(row)

            # Make sure inner and outer are partitions
            if (any(inner[i+1] > inner[i] for i in range(len(st) - 1)) or
                any(outer[i+1] > outer[i] for i in range(len(st) - 1))):
                raise ValueError("Input must be of skew shape")

        return self._new_element(st)

    def from_expr(self, expr, check=True):
        r"""
        Return a :class:`SkewTableau` from a MuPAD-Combinat expr for a skew
        tableau. The first list in ``expr`` is the inner shape of the skew
        tableau. The second list are the entries in the rows of the skew
        tableau from bottom to top.

        Provided primarily for compatibility with MuPAD-Combinat.

        EXAMPLES::

            sage: SkewTableaux().from_expr([[1,1],[[5],[3,4],[1,2]]])
            [[None, 1, 2], [None, 3, 4], [5]]
        """
        inner, outer = expr
        inner = inner+[0]*(len(outer)-len(inner))

        if check:
            # TODO: make sure inner is weakly decreasing; make sure shape
            # is skew
            pass

        st = []
        for i in range(len(outer)):
            st.append( [None]*(inner[i]) + outer[-(i+1)] )

        return self._new_element(st)

    def from_shape_and_word(self, shape, word, check=True):
        r"""
        Return the skew tableau corresponding to the skew partition ``shape``
        and the word ``word`` obtained from the row reading.

        EXAMPLES::

            sage: t = SkewTableau([[None, 1, 3], [None, 2], [4]])
            sage: shape = t.shape()
            sage: word  = t.to_word()
            sage: SkewTableaux().from_shape_and_word(shape, word)
            [[None, 1, 3], [None, 2], [4]]
        """
        # TODO: add validation

        inner, outer = shape

        st = [ [None]*row_length for row_length in inner ]
        w_count = 0
        for i in reversed(range(len(inner))):
            for j in range(inner[i]):
                if i >= len(outer) or j >= outer[i]:
                    st[i][j] = word[w_count]
                    w_count += 1

        return self._new_element(st)
    
    def from_dict(self, dct, check=True):
        r"""
        Return the skew tableau obtained from the given dictionary.
        """
        # TODO: add validation
        
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

        st = []
        for row_index in range(min_row_index, max_row_index + 1):
            row = rows[row_index]
            if not row:
                st.append([])
                continue
            tmp = [None]*(max(row.keys()) - min_col_index + 1)
            for col_index, value in six.iteritems(row):
                tmp[col_index - min_col_index] = value
            st.append(tmp)
        
        return self._new_element(st)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SkewTableaux()
            Skew tableaux
        """
        return "Skew tableaux"

    def __contains__(self, other):
        # The default contains method doesn't trust conversions
        #   enough. TODO: make ticket to add coersions from
        #   SkewTableaux to List.
        try:
            self(0, other)
        except Exception:
            return False
        return True

# TODO: make this into a factory function rather than a class
class SemistandardSkewTableaux(SkewTableaux):
    r"""
    Semistandard skew tableaux.

    TODO: standardize Parent docstrings documenting input types.

    This class can be initialized with several optional variables:
    the size of the skew tableaux (as a nameless integer variable),
    their shape (as a nameless skew partition variable), their
    weight (:meth:`~sage.combinat.skew_tableau.SkewTableau.weight`,
    as a nameless second variable after either the size or the
    shape) and their maximum entry (as an optional keyword variable
    called ``max_entry``, unless the weight has been specified). If
    neither the weight nor the maximum entry is specified, the
    maximum entry defaults to the size of the tableau.

    Note that "maximum entry" does not literally mean the highest
    entry; instead it is just an upper bound that no entry is
    allowed to surpass.

    EXAMPLES:

    The (infinite) class of all semistandard skew tableaux::

        sage: SemistandardSkewTableaux()
        Semistandard skew tableaux

    The (still infinite) class of all semistandard skew tableaux
    with maximum entry `2`::

        sage: SemistandardSkewTableaux(max_entry=2)
        Semistandard skew tableaux with maximum entry 2

    The class of all semistandard skew tableaux of given size `3`
    and maximum entry `3`::

        sage: SemistandardSkewTableaux(3)
        Semistandard skew tableaux of size 3 and maximum entry 3

    To set a different maximum entry::

        sage: SemistandardSkewTableaux(3, max_entry = 7)
        Semistandard skew tableaux of size 3 and maximum entry 7

    Specifying a shape::

        sage: SemistandardSkewTableaux([[2,1],[]])
        Semistandard skew tableaux of shape [2, 1] / [] and maximum entry 3

    Specifying both a shape and a maximum entry::

        sage: S = SemistandardSkewTableaux([[2,1],[1]], max_entry = 3); S
        Semistandard skew tableaux of shape [2, 1] / [1] and maximum entry 3
        sage: S.list()
        [[[None, 1], [1]],
         [[None, 2], [1]],
         [[None, 1], [2]],
         [[None, 3], [1]],
         [[None, 1], [3]],
         [[None, 2], [2]],
         [[None, 3], [2]],
         [[None, 2], [3]],
         [[None, 3], [3]]]

        sage: for n in range(5):
        ....:     print n, len(SemistandardSkewTableaux([[2,2,1],[1]], max_entry = n))
        0 0
        1 0
        2 1
        3 9
        4 35

    Specifying a shape and a weight::

        sage: SemistandardSkewTableaux([[2,1],[]],[2,1])
        Semistandard skew tableaux of shape [2, 1] / [] and weight [2, 1]

    (the maximum entry is redundant in this case and thus is ignored).

    Specifying a size and a weight::

        sage: SemistandardSkewTableaux(3, [2,1])
        Semistandard skew tableaux of size 3 and weight [2, 1]

    .. WARNING::

        If the shape is not specified, the iterator of this class
        yields only skew tableaux whose shape is reduced, in the
        sense that there are no empty rows before the last nonempty
        row, and there are no empty columns before the last
        nonempty column. (Otherwise it would go on indefinitely.)

    .. WARNING::

        This class acts as a factory. The resulting classes are mainly
        useful for iteration. Do not rely on their containment tests,
        as they are not correct, e. g.::

            sage: SkewTableau([[None]]) in SemistandardSkewTableaux(2)
            True
    """
    Element = SemistandardSkewTableau

    # TODO: move to ..._all class when this is made into a factory
    # TODO: add coersion to SkewTableaux
    def _element_constructor_(self, *args, **kwds):
        ret = super(SemistandardSkewTableaux, self)._element_constructor_(*args, **kwds)
        if "check" not in kwds or ("check" in kwds and kwds["check"]):
            if not ret.is_semistandard():
                raise ValueError("Input is not semistandard")
        return ret

    @staticmethod
    def __classcall_private__(cls, p=None, mu=None, max_entry=None):
        """
        Return the correct parent based upon the input.

        EXAMPLES::

            sage: SSST1 = SemistandardSkewTableaux([[3, 2, 1], [1, 1]])
            sage: SSST2 = SemistandardSkewTableaux(SkewPartition([[3, 2, 1], [1, 1]]))
            sage: SSST1 is SSST2
            True
        """
        if p is None:
            if mu is None:
                return SemistandardSkewTableaux_all(max_entry)
            raise ValueError("You must specify either a size or a shape")

        if isinstance(p, (int, Integer)):
            if mu is None:
                return SemistandardSkewTableaux_size(p, max_entry)
            else:
                return SemistandardSkewTableaux_size_weight(p, mu)

        if p in SkewPartitions():
            if mu is None:
                return SemistandardSkewTableaux_shape(p, max_entry)
            else:
                return SemistandardSkewTableaux_shape_weight(p, mu)

        raise ValueError("Invalid input")

class SemistandardSkewTableaux_all(SemistandardSkewTableaux):
    """
    Class of all semistandard skew tableaux, possibly with a given
    maximum entry.
    """
    def __init__(self, max_entry):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SemistandardSkewTableaux()
            sage: TestSuite(S).run()

            sage: S = SemistandardSkewTableaux(3)
            sage: TestSuite(S).run()
        """
        SemistandardSkewTableaux.__init__(self, category=InfiniteEnumeratedSets())
        if max_entry is None:
            self.max_entry = PlusInfinity()
        else:
            self.max_entry = max_entry

    def _repr_(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux()
            Semistandard skew tableaux
        """
        if self.max_entry == PlusInfinity():
            return "Semistandard skew tableaux"
        return "Semistandard skew tableaux with maximum entry {}".format(self.max_entry)

    def __iter__(self):
        """
        Iterate over the elements of ``self``.

        EXAMPLES::

            sage: it = SemistandardSkewTableaux(max_entry = 5).__iter__()
            sage: [next(it) for x in range(12)]
            [[],
             [[1]],
             [[2]],
             [[3]],
             [[4]],
             [[5]],
             [[1, 1]],
             [[1, 2]],
             [[1, 3]],
             [[1, 4]],
             [[1, 5]],
             [[2, 2]]]

        If no max entry is specified, the iteration goes over all
        semistandard skew tableaux of size `n` with max entry `n`,
        for all `n`::

            sage: it = SemistandardSkewTableaux().__iter__()
            sage: [next(it) for x in range(10)]
            [[],
             [[1]],
             [[1, 1]],
             [[1, 2]],
             [[2, 2]],
             [[1], [2]],
             [[None, 1], [1]],
             [[None, 2], [1]],
             [[None, 1], [2]],
             [[None, 2], [2]]]
        """
        if self.max_entry == PlusInfinity():
            # Old behavior, kept here for backwards compatibility.
            # The usefulness of this iterator is questionable.
            n = 0
            while True:
                for ssst in SemistandardSkewTableaux_size(n, n):
                    yield self._new_element(ssst)
                n += 1
        else:
            n = 0
            while True:
                for ssst in SemistandardSkewTableaux_size(n, self.max_entry):
                    yield self._new_element(ssst)
                n += 1

class SemistandardSkewTableaux_size(SemistandardSkewTableaux):
    """
    Class of all semistandard skew tableaux of a fixed size `n`,
    possibly with a given maximum entry which defaults to `n`.
    """
    def __init__(self, n, max_entry):
        """
        EXAMPLES::

            sage: S = SemistandardSkewTableaux(3)
            sage: TestSuite(S).run()
        """
        self.n = n
        if max_entry is None:
            self.max_entry = n
        else:
            self.max_entry = max_entry
        SemistandardSkewTableaux.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux(3)
            Semistandard skew tableaux of size 3 and maximum entry 3
            sage: SemistandardSkewTableaux(3, max_entry=8)
            Semistandard skew tableaux of size 3 and maximum entry 8
        """
        return ("Semistandard skew tableaux of size %s and maximum entry %s"
                %(repr(self.n), repr(self.max_entry)))

    def cardinality(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux(2).cardinality()
            8
        """
        return sum(SemistandardSkewTableaux_shape(p, self.max_entry).cardinality()
                   for p in SkewPartitions(self.n))

    def __iter__(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux(2).list()
            [[[1, 1]],
             [[1, 2]],
             [[2, 2]],
             [[1], [2]],
             [[None, 1], [1]],
             [[None, 2], [1]],
             [[None, 1], [2]],
             [[None, 2], [2]]]
        """
        for p in SkewPartitions(self.n):
            for ssst in SemistandardSkewTableaux_shape(p, self.max_entry):
                # TODO: why isn't this just yielding ssst?
                yield self._new_element(ssst)

class SemistandardSkewTableaux_size_weight(SemistandardSkewTableaux):
    r"""
    Class of semistandard tableaux of a fixed size `n` and weight `\mu`.
    """
    @staticmethod
    def __classcall_private__(cls, n, mu):
        """
        Normalize our input to ensure we have a unique representation.

        EXAMPLES::

            sage: S = SemistandardSkewTableaux(3, [2,1])
            sage: S2 = SemistandardSkewTableaux(int(3), (2,1))
            sage: S is S2
            True
        """
        # TODO: is the whole point of this method for use in the cardinality and
        # __iter__ methods? If so, remove in favor of a more direct call
        return super(SemistandardSkewTableaux_size_weight, cls).__classcall__(cls, n, tuple(mu))

    def __init__(self, n, mu):
        """
        EXAMPLES::

            sage: S = SemistandardSkewTableaux(3,[2,1])
            sage: TestSuite(S).run()
        """
        self.n = n
        self.mu = mu
        SemistandardSkewTableaux.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux(3,[2,1])
            Semistandard skew tableaux of size 3 and weight [2, 1]
        """
        return ("Semistandard skew tableaux of size %s and weight %s"
                %(repr(self.n),list(self.mu)))

    def cardinality(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux(2,[1,1]).cardinality()
            4
        """
        return sum(SemistandardSkewTableaux_shape_weight(p, self.mu).cardinality()
                   for p in SkewPartitions(self.n))

    def __iter__(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux(2,[1,1]).list()
            [[[1, 2]], [[1], [2]], [[None, 2], [1]], [[None, 1], [2]]]
        """
        for p in SkewPartitions(self.n):
            for ssst in SemistandardSkewTableaux_shape_weight(p, self.mu):
                # TODO: should this just yield ssst?
                yield self._new_element(ssst)

class SemistandardSkewTableaux_shape(SemistandardSkewTableaux):
    r"""
    Class of semistandard skew tableaux of a fixed skew shape
    `\lambda / \mu` with a given max entry.

    A semistandard skew tableau with max entry `i` is required to have all
    its entries less or equal to `i`. It is not required to actually
    contain an entry `i`.

    INPUT:

    - ``p`` -- A skew partition

    - ``max_entry`` -- The max entry; defaults to the size of ``p``.

    .. WARNING::

        Input is not checked; please use :class:`SemistandardSkewTableaux` to
        ensure the options are properly parsed.
    """
    @staticmethod
    def __classcall_private__(cls, p, max_entry=None):
        """
        Normalize our input to ensure we have a unique representation.

        EXAMPLES::

            sage: S = SemistandardSkewTableaux([[2,1],[]])
            sage: S2 = SemistandardSkewTableaux(SkewPartition([[2,1],[]]))
            sage: S is S2
            True
        """
        # TODO: is the whole point of this method for use in the cardinality and
        # __iter__ methods? If so, remove in favor of a more direct call
        if max_entry is None:
            max_entry = sum(p[0])-sum(p[1])
        return super(SemistandardSkewTableaux_shape, cls).__classcall__(cls, SkewPartition(p), max_entry)

    def __init__(self, p, max_entry):
        """
        EXAMPLES::

            sage: S = SemistandardSkewTableaux([[2,1],[]])
            sage: S == loads(dumps(S))
            True
            sage: TestSuite(S).run()
        """
        self.p = p
        self.max_entry = max_entry
        SemistandardSkewTableaux.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]])
            Semistandard skew tableaux of shape [2, 1] / [] and maximum entry 3
        """
        return ("Semistandard skew tableaux of shape %s and maximum entry %s"
                %(repr(self.p), repr(self.max_entry)))

    def cardinality(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]]).cardinality()
            8
            sage: SemistandardSkewTableaux([[2,1],[]], max_entry=2).cardinality()
            2
        """
        return sum(SemistandardSkewTableaux_shape_weight(self.p, mu).cardinality()
                   for mu in IntegerVectors(self.p.size(), self.max_entry))

    def __iter__(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]]).list()
            [[[1, 1], [2]],
             [[1, 1], [3]],
             [[1, 2], [2]],
             [[1, 3], [2]],
             [[1, 2], [3]],
             [[1, 3], [3]],
             [[2, 2], [3]],
             [[2, 3], [3]]]
            sage: from sage.combinat.skew_tableau import SemistandardSkewTableaux_shape
            sage: SemistandardSkewTableaux_shape([[2,1],[]], max_entry=2).list()
            [[[1, 1], [2]], [[1, 2], [2]]]
        """
        for mu in IntegerVectors(self.p.size(), self.max_entry):
            for ssst in SemistandardSkewTableaux_shape_weight(self.p, mu):
                # TODO: should this just yield ssst?
                yield self._new_element(ssst)

class SemistandardSkewTableaux_shape_weight(SemistandardSkewTableaux):
    r"""
    Class of semistandard skew tableaux of a fixed skew shape `\lambda / \nu`
    and weight `\mu`.
    """
    @staticmethod
    def __classcall_private__(cls, p, mu):
        """
        Normalize our input to ensure we have a unique representation.

        EXAMPLES::

            sage: S = SemistandardSkewTableaux([[2,1],[]], [2,1])
            sage: S2 = SemistandardSkewTableaux(SkewPartition([[2,1],[]]), (2,1))
            sage: S is S2
            True
        """
        # TODO: is the whole point of this method for use in the cardinality and
        # __iter__ methods? If so, remove in favor of a more direct call
        p = SkewPartition(p)
        mu = tuple(mu)
        return super(SemistandardSkewTableaux_shape_weight, cls).__classcall__(cls, p, mu)

    def __init__(self, p, mu):
        """
        EXAMPLES::

            sage: S = SemistandardSkewTableaux([[2,1],[]],[2,1])
            sage: S == loads(dumps(S))
            True
            sage: TestSuite(S).run()
        """
        self.p = p
        self.mu = mu
        SemistandardSkewTableaux.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]],[2,1])
            Semistandard skew tableaux of shape [2, 1] / [] and weight [2, 1]
        """
        return "Semistandard skew tableaux of shape %s and weight %s"%(repr(self.p), list(self.mu))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]],[2,1]).list()
            [[[1, 1], [2]]]
        """
        from sage.combinat.ribbon_tableau import RibbonTableaux_shape_weight_length
        for x in RibbonTableaux_shape_weight_length(self.p, self.mu, 1):
            # TODO: should this just yield x?
            yield self._new_element(x)

class StandardSkewTableaux(SemistandardSkewTableaux):
    """
    Standard skew tableaux.

    EXAMPLES::

        sage: S = StandardSkewTableaux(); S
        Standard skew tableaux
        sage: S.cardinality()
        +Infinity

    ::

        sage: S = StandardSkewTableaux(2); S
        Standard skew tableaux of size 2
        sage: S.cardinality()
        4

    ::

        sage: StandardSkewTableaux([[3, 2, 1], [1, 1]]).list()
        [[[None, 1, 2], [None, 3], [4]],
         [[None, 1, 2], [None, 4], [3]],
         [[None, 1, 3], [None, 2], [4]],
         [[None, 1, 4], [None, 2], [3]],
         [[None, 1, 3], [None, 4], [2]],
         [[None, 1, 4], [None, 3], [2]],
         [[None, 2, 3], [None, 4], [1]],
         [[None, 2, 4], [None, 3], [1]]]
    """
    Element = StandardSkewTableau

    # TODO: move to ..._all class when this is made into a factory
    # TODO: add coersion to SemistandardSkewTableaux
    def _element_constructor_(self, *args, **kwds):
        ret = super(StandardSkewTableaux, self)._element_constructor_(*args, **kwds)
        if "check" not in kwds or ("check" in kwds and kwds["check"]):
            # Semistandardness has already been checked
            if not ret.has_standard_entries():
                raise ValueError("Input is not standard")
        return ret

    @staticmethod
    def __classcall_private__(cls, skp=None):
        """
        Return the class of standard skew tableaux of skew shape ``skp``.

        EXAMPLES::

            sage: SST1 = StandardSkewTableaux([[3, 2, 1], [1, 1]])
            sage: SST2 = StandardSkewTableaux(SkewPartition([[3, 2, 1], [1, 1]]))
            sage: SST1 is SST2
            True
        """
        if skp is None:
            return type.__call__(StandardSkewTableaux)
        elif isinstance(skp, (int, Integer)):
            return StandardSkewTableaux_size(skp)
        elif skp in SkewPartitions():
            return StandardSkewTableaux_shape(skp)
        else:
            raise TypeError("Invalid argument")

    def __init__(self):
        """
        EXAMPLES::

            sage: s = StandardSkewTableaux()
            sage: TestSuite(s).run()
        """
        super(StandardSkewTableaux, self).__init__(category=InfiniteEnumeratedSets())

    def _coerce_map_from_(self, S):
        return isinstance(S, StandardSkewTableaux)

    def _repr_(self):
        """
        EXAMPLES::

            sage: StandardSkewTableaux()
            Standard skew tableaux
        """
        return "Standard skew tableaux"

    def __iter__(self):
        """
        Iterate through all standard skew tableaux having
        no empty rows (before nonempty rows) and no empty columns
        (before nonempty columns).

        EXAMPLES::

            sage: it = StandardSkewTableaux().__iter__()
            sage: [next(it) for x in range(10)]
            [[],
             [[1]],
             [[1, 2]], [[1], [2]], [[None, 1], [2]], [[None, 2], [1]],
             [[1, 2, 3]], [[1, 2], [3]], [[1, 3], [2]],
             [[None, 1, 2], [3]]]
        """
        n = 0
        while True:
            for st in StandardSkewTableaux_size(n):
                # TODO: should this just yield st?
                yield self._new_element(st)
            n += 1

class StandardSkewTableaux_size(StandardSkewTableaux):
    """
    Standard skew tableaux of a fixed size `n`.
    """
    def __init__(self, n):
        """
        EXAMPLES::

            sage: S = StandardSkewTableaux(3)
            sage: TestSuite(S).run()
        """
        self.n = n
        super(StandardSkewTableaux, self).__init__(category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: StandardSkewTableaux(3)
            Standard skew tableaux of size 3
        """
        return "Standard skew tableaux of size %s"%self.n

    def cardinality(self):
        """
        EXAMPLES::

            sage: StandardSkewTableaux(1).cardinality()
            1
            sage: StandardSkewTableaux(2).cardinality()
            4
            sage: StandardSkewTableaux(3).cardinality()
            24
            sage: StandardSkewTableaux(4).cardinality()
            194
        """
        return sum(StandardSkewTableaux_shape(skp).cardinality()
                   for skp in SkewPartitions(self.n))

    def __iter__(self):
        """
        Iterate through all standard skew tableaux of size `n` having
        no empty rows (before nonempty rows) and no empty columns
        (before nonempty columns). (The last two requirements
        ensure that the iterator terminates after finitely many steps.)

        EXAMPLES::

            sage: StandardSkewTableaux(2).list()
            [[[1, 2]], [[1], [2]], [[None, 1], [2]], [[None, 2], [1]]]

            sage: StandardSkewTableaux(3).list()
            [[[1, 2, 3]],
             [[1, 2], [3]], [[1, 3], [2]],
             [[None, 1, 2], [3]], [[None, 1, 3], [2]],
             [[None, 2, 3], [1]],
             [[None, 1], [2, 3]], [[None, 2], [1, 3]],
             [[None, None, 1], [2, 3]], [[None, None, 2], [1, 3]], [[None, None, 3], [1, 2]],
             [[1], [2], [3]],
             [[None, 1], [None, 2], [3]], [[None, 1], [None, 3], [2]], [[None, 2], [None, 3], [1]],
             [[None, 1], [2], [3]], [[None, 2], [1], [3]], [[None, 3], [1], [2]],
             [[None, None, 1], [None, 2], [3]], [[None, None, 1], [None, 3], [2]],
             [[None, None, 2], [None, 1], [3]], [[None, None, 3], [None, 1], [2]],
             [[None, None, 2], [None, 3], [1]], [[None, None, 3], [None, 2], [1]]]
        """
        for skp in SkewPartitions(self.n):
            for sst in StandardSkewTableaux_shape(skp):
                # TODO: should this just yield sst?
                yield self._new_element(sst)

class StandardSkewTableaux_shape(StandardSkewTableaux):
    r"""
    Standard skew tableaux of a fixed skew shape `\lambda / \mu`.
    """
    @staticmethod
    def __classcall_private__(cls, skp):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S = StandardSkewTableaux([[3, 2, 1], [1, 1]])
            sage: S2 = StandardSkewTableaux(SkewPartition([[3, 2, 1], [1, 1]]))
            sage: S is S2
            True
        """
        # TODO: is the whole point of this method for use in the cardinality and
        # __iter__ methods? If so, remove in favor of a more direct call
        return super(StandardSkewTableaux_shape, cls).__classcall__(cls, SkewPartition(skp))

    def __init__(self, skp):
        """
        TESTS::

            sage: S = StandardSkewTableaux([[3, 2, 1], [1, 1]])
            sage: TestSuite(S).run()
        """
        self.skp = skp
        super(StandardSkewTableaux, self).__init__(category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: StandardSkewTableaux([[3, 2, 1], [1, 1]])
            Standard skew tableaux of shape [3, 2, 1] / [1, 1]
        """
        return "Standard skew tableaux of shape %s"%repr(self.skp)

    def cardinality(self):
        """
        Return the number of standard skew tableaux with shape of the skew
        partition ``skp``. This uses a formula due to Aitken
        (see Cor. 7.16.3 of [Sta-EC2]_).

        EXAMPLES::

            sage: StandardSkewTableaux([[3, 2, 1], [1, 1]]).cardinality()
            8
        """
        outer, inner = self.skp
        m = len(outer)
        n = sum(outer) - sum(inner)
        outer = list(outer)
        inner = list(inner) + [0]*(m-len(inner))
        a = zero_matrix(QQ, m)
        for i in range(m):
            for j in range(m):
                v = outer[i] - inner[j] - i + j
                if v < 0:
                    a[i,j] = 0
                else:
                    a[i,j] = 1/factorial(v)
        return ZZ(factorial(n) * a.det())

    def __iter__(self):
        """
        An iterator for all the standard skew tableaux with shape of the
        skew partition ``skp``. The standard skew tableaux are ordered
        lexicographically by the word obtained from their row reading.

        EXAMPLES::

            sage: [st for st in StandardSkewTableaux([[3, 2, 1], [1, 1]])]
            [[[None, 1, 2], [None, 3], [4]],
             [[None, 1, 2], [None, 4], [3]],
             [[None, 1, 3], [None, 2], [4]],
             [[None, 1, 4], [None, 2], [3]],
             [[None, 1, 3], [None, 4], [2]],
             [[None, 1, 4], [None, 3], [2]],
             [[None, 2, 3], [None, 4], [1]],
             [[None, 2, 4], [None, 3], [1]]]
        """
        dag = self.skp.to_dag(format="tuple")
        le_list = list(dag.topological_sort_generator())

        empty = [[None]*row_length for row_length in self.skp.outer()]

        for le in le_list:
            # TODO: switch to _new_element to avoid a circular reference
            # and work in the _label_skew code.
            yield self._new_element(_label_skew(le, empty))

def _label_skew(list_of_cells, sk):
    """
    Return a filled-in standard standard skew tableau given an
    ordered list ``list_of_cells`` of the coordinates to fill in
    (as pairs) and an empty shape ``sk``.

    EXAMPLES::

        sage: import sage.combinat.skew_tableau as skew_tableau
        sage: l = [(0, 0), (1, 1), (1, 0), (0, 1)]
        sage: empty = [[None,None],[None,None]]
        sage: skew_tableau._label_skew(l, empty)
        [[1, 4], [3, 2]]
    """
    i = 1
    skew = [list(row) for row in sk]
    for row, column in list_of_cells:
            skew[row][column] = i
            i += 1
    return skew

# TODO: deal with pickling/unpickling. See the end of the old skew_tableau.py.