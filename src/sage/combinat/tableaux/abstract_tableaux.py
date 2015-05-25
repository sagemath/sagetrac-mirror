r"""
AbstractTableaux Parent class.

This is an abstract base class for numerous tableaux-like
Parent subclasses. See AbstractTableau for the corresponding
Element class.

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
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import Sets

from sage.combinat.tableaux.abstract_tableau import AbstractTableau

class AbstractTableaux(UniqueRepresentation, Parent):
    r"""
    Parent class of all abstract tableaux.
    
    See AbstractTableau for the Element class.
    """
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

    def _element_constructor_(self, t):
        r"""
        Constructs an object from ``t`` as an element of ``self``,
        if possible. Input validation and normalization should be done here.

        INPUT:

        - ``t`` -- Data which can be interpreted as a tableau

        OUTPUT:

        - The corresponding tableau object
        """
        if not t in self:
            raise ValueError("%s is not an element of %s."%(t, self))

        return self.element_class(self, t)
