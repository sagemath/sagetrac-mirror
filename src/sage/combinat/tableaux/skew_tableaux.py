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

from sage.combinat.tableaux.bad_shape_tableaux import BadShapeTableaux
from sage.combinat.tableaux.skew_tableau import SkewTableau

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
    Parent class of all skew tableaux.
    
    See SkewTableau for the Element class.
    """
    Element = SkewTableau

    # Convenient shortcut to global options
    global_options = TableauOptions

    def _element_constructor_(self, st=None, expr=None, shape_word=None,
                                    dct=None, check=True):
        r"""
        Construct a new SkewTableau using one of several input formats,
        optionally validating input.
        
        If multiple formats are specified, the left-most is used. If
        no format is specified, the "trivial" skew tableau is returned.

        INPUT:
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
        - ``check`` -- (default: ``True``) if ``True``, then validate input:
          ensure ``st`` or the cells of ``dct`` actually form a skew shape,
          etc.
        """
        if st is not None:
            # TODO: normalize st input by removing rows with all None's and
            # empty rows.
            if check:
                try:
                    st = tuple(tuple(row) for row in st)
                except TypeError:
                    raise TypeError("each element of the skew tableau must be an iterable")

                # TODO: make sure None's are contiguous and left-justified;
                # make sure shape is a skew tableau

            return self._new_element(st)

        if expr is not None:
            return self.from_expr(expr, check)

        if shape_word is not None:
            shape, word = shape_word
            return self.from_shape_and_word(shape, word, check)

        if dct is not None:
            return self.from_dict(dct, check)

        return self._new_element([[]])

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

