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

class SkewTableaux(BadShapeTableaux):
    r"""
    Parent class of all skew tableaux.
    
    See SkewTableau for the Element class.
    """
    Element = SkewTableau

    def _element_constructor_(self, skp=None, expr=None, shape_word=None,
                                    dct=None, check=True):
        r"""
        Construct a new SkewTableau using one of several input formats,
        optionally validating input.
        
        If multiple formats are specified, the left-most is used. If
        no format is specified, the "trivial" skew tableau is returned.

        INPUT:
        - ``skp`` -- an iterable of rows from top to bottom in English
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
        if skp is not None:
            if check:
                try:
                    skp = tuple(tuple(row) for row in skp)
                except TypeError:
                    raise TypeError("each element of the skew tableau must be an iterable")

                # TODO: make sure None's are contiguous and left-justified;
                # make sure shape is a skew tableau

            return self._new_element(skp)
        
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

        skp = []
        for i in range(len(outer)):
            skp.append( [None]*(inner[i]) + outer[-(i+1)] )

        return self._new_element(skp)

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
        
        
        skp = [ [None]*row_length for row_length in inner ]
        w_count = 0
        for i in reversed(range(len(inner))):
            for j in range(inner[i]):
                if i >= len(outer) or j >= outer[i]:
                    skp[i][j] = word[w_count]
                    w_count += 1

        return self._new_element(skp)
    
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

        skp = []
        for row_index in range(min_row_index, max_row_index + 1):
            row = rows[row_index]
            if not row:
                skp.append([])
                continue
            tmp = [None]*(max(row.keys()) - min_col_index + 1)
            for col_index, value in six.iteritems(row):
                tmp[col_index - min_col_index] = value
            skp.append(tmp)
        
        return self._new_element(skp)

    def _repr_(self):
        r"""
        Return the representation string.

        OUTPUT:

        A string.
        """
        return "Skew Tableaux"
