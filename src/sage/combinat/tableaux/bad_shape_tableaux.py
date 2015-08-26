r"""
BadShapeTableaux Parent class.

See BadShapeTableau for the corresponding Element class.

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
from sage.rings.integer_ring                  import ZZ

from sage.combinat.tableaux.abstract_tableaux import AbstractTableaux
from sage.combinat.tableaux.bad_shape_tableau import BadShapeTableau
from sage.combinat.tableaux.bad_shape_tableau import BadShapeTableauFactory

class BadShapeTableaux(AbstractTableaux):
    r"""
    Parent class of all bad shape tableaux.
    
    See BadShapeTableau for the Element class.

    EXAMPLES::

        sage: B = BadShapeTableaux(); B
        Bad Shape Tableaux
        sage: B.an_element()

    TESTS::

        sage: B = BadShapeTableaux()
        sage: TestSuite(B).run()
    """
    Element = BadShapeTableau

    def _element_constructor_(self, dct, check=True):
        try:
            dct = dict(dct)
        except:
            raise ValueError('dct must be a dictionary')

        if check:
            if not all(x in ZZ and y in ZZ for x, y in six.iterkeys(dct)):
                raise ValueError('keys must be pairs of integers')

        return self._new_element(dct)
    # For user convenience, documentation for this method has
    #    been placed on the bad shape tableau factory function.
    _element_constructor_.__doc__ = BadShapeTableauFactory.__doc__

    def _an_element_(self):
        r"""
        Return a typical element of ``self``.

        TESTS::

            sage: b = BadShapeTableaux().an_element() # indirect doctest
            sage: sorted(b.dict().items())
            [((-1, -2), (1, 2)), ((1, 1), 0), ((1, 2), 4), ((2, -2), 'cow')]
        """
        return self({(1, 1): 0, (2, -2): 'cow',
                    (-1, -2): (1, 2), (1, 2): 4})

    def _repr_(self):
        r"""
        Return the representation string.

        TESTS::

            sage: BadShapeTableaux()
            Bad Shape Tableaux
        """
        return "Bad Shape Tableaux"
