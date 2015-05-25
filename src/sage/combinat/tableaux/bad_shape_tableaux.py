r"""
BadShapeTableaux Parent class.

See BadShapeTableau for the corresponding Element class.

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
from sage.rings.integer_ring import ZZ

from sage.combinat.tableaux.abstract_tableaux import AbstractTableaux
from sage.combinat.tableaux.bad_shape_tableau import BadShapeTableau

class BadShapeTableaux(AbstractTableaux):
    r"""
    Parent class of all bad shape tableaux.
    
    See BadShapeTableau for the Element class.
    """
    Element = BadShapeTableau

    def _element_constructor_(self, dct, check=True):
        r"""
        Construct a new BadShapeTableau, optionally validating input.

        INPUT:

        - ``dct`` -- a dictionary (or more generally something
          passable to ``dict``) whose keys are pairs of integers
        - ``check`` -- (default: ``True``) if ``True``, then check that
          the keys of ``dct`` are in fact pairs of integers
        """
        try:
            dct = dict(dct)
        except:
            raise ValueError('dct must be a dictionary')

        if check:
            try:
                all(x in ZZ and y in ZZ for x, y in six.iterkeys(dct))
            except:
                raise ValueError('keys must be pairs of integers')

        return self.element_class(self, dct)

    def _repr_(self):
        r"""
        Return the representation string.

        OUTPUT:

        A string.
        """
        return "Bad Shape Tableaux"
