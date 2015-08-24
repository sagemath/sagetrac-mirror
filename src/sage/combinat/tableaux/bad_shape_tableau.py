r"""
BadShapeTableau Element class.

See BadShapeTableaux for the corresponding Parent class.

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

from sage.combinat.tableaux.abstract_tableau import AbstractTableau

class BadShapeTableau(AbstractTableau):
    r"""
    A tableau of bad shape.

    See Parent `class:BadShapeTableaux` for construction options.

    A tableau of bad shape is a tableau whose cells are pairs of
    integers and whose entries are arbitrary.
    """
    def __init__(self, parent, dct):
        r"""
        Initialize the BadShapeTableau.

        INPUT:

        - ``dct`` -- a dictionary whose keys are pairs of integers
        """
        self._set_parent(parent)
        self._dict = dict(dct)

    def _dict_unsafe(self):
        r"""
        Return the dictionary containing the defining data of ``self``.

        OUTPUT:

        A dictionary.
        """
        return self._dict

    def filter_by_cells(self, predicate):
        r"""
        Return the subtableau of ``self`` which consists only of the
        entries in the cells satisfying the given predicate
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
        Return the subtableau of ``self`` which consists only of the
        entries satisfying the given ``predicate``.

        INPUT:

        - ``predicate`` -- a function accepting one parameter and returning
          ``True`` or ``False``

        OUTPUT:

        A BadShapeTableau
        """
        data = {k:v for k, v in six.iteritems(self._dict_unsafe())
                if predicate(v)}
        return self.__class__(dictionary=data, check=True)

    def conjugate(self):
        r"""
        Return the conjugate of ``self``.

        If `T` is a tableau of bad shape, then the conjugate of `T`
        is the tableau of bad shapes whose cells are the pairs
        `(x, y)` for `(y, x)` being cells of `T`, and which sends
        every `(x, y)` to `T(y, x)`.

        OUTPUT:

        A BadShapeTableau.
        """
        data = {(k[1], k[0]): v for k, v in self._dict_unsafe()}
        return self.__class__(dictionary=data, check=True)

    # Alias
    transpose = conjugate

# Use a factory method to create `class:BadShapeTableau`.
def BadShapeTableauFactory(*args, **kwds):
    r"""
    (See `class:BadShapeTableau` for docstring.)
    """
    from bad_shape_tableaux import BadShapeTableaux
    return BadShapeTableaux()(*args, **kwds)
BadShapeTableauFactory.__doc__ = BadShapeTableau.__doc__
