"""
This file define some tools to obtain ascii and latex outputs for any tableaux.
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
from sage.sets.set import Set
from sage.categories.object_with_options import ObjectsWithOptions
from sage.misc.ascii_art import AsciiArt
from sage.misc.ascii_art import ascii_art
from sage.misc.graphic_tableau import GraphicTableau

class GraphicOutputsFor2DTableaux( Category ):

    def _repr_( self ):
        return "Category of graphic outputs for 2d tableaux"

    graphic_tableaux_options = GlobalOptions(
        name = 'Options for set of tableaux',
        doc=r"""
        """,
        end_doc=r"""
        """
    )

    graphic_tableaux_options_for_element = GlobalOptions(
        name = 'Options for displaying 2D tableaux',
        doc=r"""
        """,
        end_doc=r"""
        """,
        drawing_options=dict(
            default = GraphicTableau.default_graphic_options(),
            description='the drawing options',
            checker=lambda x: isinstance(x, GlobalOtions)
        ),
        ascii_options=dict(
            default = GraphicTableau.default_ascii_options(),
            description='the ascii art options',
            checker=lambda x: isinstance(x, GlobalOtions)
        ),
        repr=dict(
            default="list",
            values= dict(
                list='displayed as list',
                ascii='as a ascii art'
            )
        ),
        latex=dict(
            default="drawing",
            values= dict(
                list='displayed as list',
                drawing='as a drawing',
            )
        )
    )

    def super_categories( self ):
        return [ ObjectsWithOptions() ]

    class ParentMethods:
        def _default_options( self ):
            return GraphicOutputsFor2DTableaux.graphic_tableaux_options 

        def _default_options_for_element( self ):
            return GraphicOutputsFor2DTableaux.graphic_tableaux_options_for_element 

    class ElementMethods:
        def _repr_tableau(self):
            return self.get_options().dispatch(self, '_repr_', 'repr')

        def _repr_list(self):
            return list( self ).__repr__()

        def _repr_ascii(self):
            return ascii_art ( self )._repr_()

        def _to_graphic_tableau( self ):
            if 'separators' in dir( self ):
                separators = self.separators()
            else:
                separators = None
            if '_filter_for_ascii' in dir( self ):
                filter_for_ascii = self._filter_for_ascii()
            else:
                filter_for_ascii = None
            return GraphicTableau( 
                self.array(), separators, filter_for_ascii,
                self.get_options()['drawing_options'],
                self.get_options()['ascii_options']
            )

        def _ascii_art( self ):
            return ascii_art( self._to_graphic_tableau() )

        def _ascii_art_( self ):
            return self._ascii_art()

        def _latex_(self):
            return self.get_options().dispatch(self, '_latex_', 'latex')

        def _latex_drawing(self):
            return latex( self._to_graphic_tableau() )

        def _latex_list(self):
            return latex( list(self) )

