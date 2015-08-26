r"""
BadShapeTableau Element class.

See :class:`BadShapeTableaux` for the corresponding Parent class.

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

from sage.combinat.tableaux.abstract_tableau import AbstractTableau

class BadShapeTableau(AbstractTableau):
    r"""
    A tableau of bad shape.

    See Parent class:`BadShapeTableaux` for construction options.

    A tableau of bad shape is a tableau whose cells are pairs of
    integers and whose entries are arbitrary.

    EXAMPLES::

        sage: b = BadShapeTableau({(7, 7): 7, (9, 6): 3});
        sage: b.to_word()
        word: 37
    """
    def __init__(self, parent, dct):
        r"""
        Initialize the BadShapeTableau.

        INPUT:

        - ``dct`` -- a dictionary whose keys are pairs of integers

        TESTS::
        
            sage: b = BadShapeTableau({(4, 4): 4, (3, 2): -1})
            sage: TestSuite(b).run()
        """
        self._set_parent(parent)
        self._dict = dict(dct)

    def _dict_unsafe(self):
        r"""
        Return the dictionary containing the defining data of ``self``.

        TESTS::

            sage: b = BadShapeTableau({(2, 3): 4, (5, 6): 7})
            sage: b._dict_unsafe() == {(2, 3): 4, (5, 6): 7}
            True
        """
        return self._dict

    def filter_by_cells(self, predicate):
        r"""
        Return the subtableau of ``self`` obtained by removing cells
        which do not satisfy the given predicate ``predicate``.

        INPUT:

        - ``predicate`` -- a function accepting two parameters, namely the
          row and column indexes of a cell, and returning ``True`` or ``False``

        OUTPUT:

        A class:`BadShapeTableau`

        TEST::

            sage: b = BadShapeTableau({(6, 7): 'cat', (2, 2): 5, (3, 4): 0})
            sage: b2 = b.filter_by_cells(lambda r, c: r <= 3)
            sage: b2.dict() == {(2, 2): 5, (3, 4): 0}
            True
        """
        data = {k:v for k, v in six.iteritems(self._dict_unsafe())
                if predicate(*k)}
        return self.parent()(data, check=True)

    def filter_by_entries(self, predicate):
        r"""
        Return the subtableau of ``self`` obtained by removing cells
        whose entries do not satisfy the given predicate ``predicate``.

        INPUT:

        - ``predicate`` -- a function accepting one parameter, namely the
          entry of a cell, and returning ``True`` or ``False``

        OUTPUT:

        A class:`BadShapeTableau`

        TEST::

            sage: b = BadShapeTableau({(6, 7): 'cat', (2, 2): 5, (3, 4): 0})
            sage: b.filter_by_entries(lambda v: type(v) is str)
            {(6, 7): 'cat'}
        """
        data = {k:v for k, v in six.iteritems(self._dict_unsafe())
                if predicate(v)}
        return self.parent()(data, check=True)

    def conjugate(self):
        r"""
        Return the conjugate of ``self``.

        The conjugate of a tableaux of bad shape is the tableau of
        bad shape obtained by swapping row and column indexes in
        each cell.

        OUTPUT:

        A class:`BadShapeTableau`

        TEST::

            sage: b = BadShapeTableau({(6, 7): 'cat', (2, 2): 5, (3, 4): 0})
            sage: b.conjugate().dict() == {(7, 6): 'cat', (2, 2): 5, (4, 3): 0}
            True
        """
        data = {(k[1], k[0]): v for k, v in six.iteritems(self._dict_unsafe())}
        return self.parent()(data, check=True)

    # Alias
    transpose = conjugate

# Use a factory method to create class:`BadShapeTableau`.
def BadShapeTableauFactory(*args, **kwds):
    r"""
    Construct a new BadShapeTableau, optionally validating input.

    INPUT:

    - ``dct`` -- a dictionary (or more generally something
      passable to ``dict``) whose keys are pairs of integers
    - ``check`` -- (default: ``True``) if ``True``, then check that
      the keys of ``dct`` are in fact pairs of integers

    EXAMPLES::

        sage: BadShapeTableau({(9, 8): 7, (6, 5): 3}).dict()[(9, 8)]
        7

    TESTS::

        sage: BadShapeTableaux()({('fly', 1): 5}, check=False) # indirect doctest
        {('fly', 1): 5}
        sage: BadShapeTableaux()({('fly', 1): 5})
        Traceback (most recent call last):
        ...
        ValueError: keys must be pairs of integers
        sage: BadShapeTableaux()('horse')
        Traceback (most recent call last):
        ...
        ValueError: dct must be a dictionary
    """
    from bad_shape_tableaux import BadShapeTableaux
    return BadShapeTableaux()(*args, **kwds)
