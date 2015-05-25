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

    def _element_constructor_(self, *args, **kwds):
        r"""
        Constructs an Element of ``self``.
        
        Input validation and normalization should be done here.

        OUTPUT:

        - The corresponding tableau object
        """
        return self._new_element(*args, **kwds)

    def _new_element(self, *args, **kwds):
        r"""
        Constructs an Element of ``self``.
        
        We assume the Element class will implement __classcall__, which
        will pick an appropriate Parent class to construct a new Element.
        Hence we must use __call__ directly to bypass using __classcall__,
        which would be circular.
        """
        return type.__call__(self.element_class, self, *args, **kwds)