"""
This file define the categorie of Objects with options.
"""
#*****************************************************************************
#       Copyright (C) 2014 Adrien Boussicault <bousica@labri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.global_options import GlobalOptions
from sage.categories.category import Category

class ObjectsWithOptions( Category ):
    def super_categories( self ):
        return []

    class ParentMethods:
        def __init_extra__(self):
            self.reset_otions()
            self.reset_otions_for_element()

        def get_options_for_element( self ):
            if self._options_for_element is None:
                return self._default_options_for_element()
            return self._options_for_element

        def set_options_for_element( self, *get_value, **set_value ):
            if self._options_for_element is None:
                self._options_for_element = deepcopy(
                    self.get_options_for_element()
                )
            self._options_for_element( *get_value, **set_value )

        def reset_otions_for_element( self ):
            self._options_for_element = None

        def get_options( self ):
            if self._options is None:
                return self._default_options()
            return self._options

        def set_options_for_element( self, *get_value, **set_value ):
            if self._options is None:
                self._options = deepcopy( self.get_options() )
            self._options( *get_value, **set_value )

        def reset_otions( self ):
            self._options = None


    class ElementMethods:
        def __init_extra__(self):
            self._options = None

        def get_options( self ):
            r"""
            Return all the options of the object.
            """
            if self._options is None:
                return self.parent().get_options_for_element()
            return self._options

        def set_options( self, *get_value, **set_value ):
            r"""
            Set new options to the object.
            """
            if self._options is None:
                self._options = deepcopy( self.get_options() )
            self._options( *get_value, **set_value )

        def reset_options():
            self._options = None
