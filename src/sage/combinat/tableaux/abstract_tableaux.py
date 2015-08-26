r"""
AbstractTableaux Parent class.

This is an abstract base class for numerous tableaux-like
Parent classes. See class:`AbstractTableau` for the
corresponding Element class.

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
from sage.structure.parent                   import Parent
from sage.structure.unique_representation    import UniqueRepresentation
from sage.categories.sets_cat                import Sets

from sage.combinat.tableaux.abstract_tableau import AbstractTableau

class AbstractTableaux(UniqueRepresentation, Parent):
    r"""
    Parent class of all abstract tableaux.

    See AbstractTableau for the Element class.

    EXAMPLES::

        sage: BadShapeTableaux()
        Bad Shape Tableaux
        sage: b = BadShapeTableaux()({(1, 2): -1, (2, 3): 'kitty'})
        sage: set(b.iter_entries()) == set(-1, 'kitty')
        True
    """
    Element = AbstractTableau

    def __init__(self, category=None):
        r"""
        Initialize the parent.

        TESTS::

            sage: AbstractTableaux().category() # indirect doctest
            Category of sets
        """
        if category is None:
            Parent.__init__(self, category=Sets())
        else:
            Parent.__init__(self, category=category)

    def _element_constructor_(self, *args, **kwds):
        r"""
        Constructs an Element of ``self``.
        
        Input validation and normalization should be done here.

        TESTS::

            sage: type(AbstractTableaux()()) # indirect doctest
            <class 'sage.combinat.tableaux.abstract_tableau.AbstractTableaux_with_category.element_class'>
        """
        return self._new_element(*args, **kwds)

    def _new_element(self, *args, **kwds):
        r"""
        Constructs an Element of ``self``.
        
        Does not do input validation or normalization.

        TESTS::

            sage: type(AbstractTableaux()()) # indirect doctest
            <class 'sage.combinat.tableaux.abstract_tableau.AbstractTableaux_with_category.element_class'>
        """
        return self.element_class(self, *args, **kwds)

    def _repr_(self):
        r"""
        Return the representation string.

        TESTS::

            sage: AbstractTableaux()
            Abstract tableaux
        """
        return "Abstract tableaux"
