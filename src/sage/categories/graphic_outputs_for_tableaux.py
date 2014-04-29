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

class TextAlignement:
    H_CENTER = 1
    H_LEFT = 2
    H_RIGHT = 3
    V_TOP = 4
    V_BOTTOM = 5
    V_CENTER = 6

class Ascii_array:
    def __init__( self, ascii_art ):
        if isinstance( ascii_art, list ):
            lines = ascii_art
        else:
            art = ascii_art._repr_()
            lines = art.split('\n')
        height = len( lines )
        width = len( lines[0] )
        for l in lines:
            if width < len( l ):
                width = len( l )
        self._lines = [ [ ' ' for w in range( width ) ] for h in range(height) ]
        for h in range( len(lines) ):
            for w in range( len( lines[h] ) ):
                self._lines[h][w] = lines[h][w]

    class Line:
        def __init__( self, lines, h ):
            self._lines = lines
            self._h = h
        def height(self):
            return len( self._lines )
        def width(self):
            return len( self._lines[0] )
        def __getitem__( self, w ):
            if (0 <= self._h < self.height()) and (0 <= w < self.width()) :
                return self._lines[self._h][w]
            return None
        def __setitem__( self, w, value ):
            if (0 <= self._h < self.height())  and (0 <= w < self.width()) :
                self._lines[self._h][w] = value

    def __getitem__( self, h ):
        return self.Line( self._lines, h )

    def width( self ):
        return len( self._lines[0] )

    def height( self ):
        return len( self._lines )

    def clip(
        self, array, pos_no, pos_se, 
        h_align = TextAlignement.H_CENTER, v_align = TextAlignement.V_CENTER 
    ):
        clip_width = pos_se[1] - pos_no[1]
        clip_height = pos_se[0] - pos_no[0]

        if not isinstance( array, Ascii_array ):
            array = Ascii_array( array  )
        array_width = array.width()
        array_height = array.height()

        if (v_align == TextAlignement.V_CENTER) :
            def H( h ):
                return int(array_height/2) - int(clip_height/2) + (h - pos_no[0])
        elif v_align == TextAlignement.V_TOP :
            def H( h ):
                return h - pos_no[0]
        elif v_align == TextAlignement.V_BOTTOM :
            def H( h ):
                return array_height - clip_height + (h - pos_no[0])

        if h_align == TextAlignement.H_CENTER :
            def W( w ):
                return int(array_width/2) - int(clip_width/2) + (w - pos_no[1])
        elif h_align == TextAlignement.H_LEFT :
            def W( w ):
                return w - pos_no[1]
        elif h_align == TextAlignement.H_RIGHT :
            def W( w ):
                return array_width - clip_width + (w - pos_no[1])

        for h in range( max(0, pos_no[0]), min(self.height(), pos_se[0]) ):
            for w in range( max(0, pos_no[1]), min(self.width(), pos_se[1]) ):
                array_h = H(h)
                array_w = W(w)
                if ( 0 <= array_h < array_height ) and ( 0 <= array_w < array_width ) :
                    self[h][w] = array[ array_h ][ array_w ]

    def get_new_geometry(
        self, height, width,
        h_align = TextAlignement.H_CENTER, v_align = TextAlignement.V_CENTER 
    ):
        res = Ascii_array( [ [ ' ' for w in range(width) ] for h in range(height)] )
        res.clip(
            self, [0,0], [width, height],
            h_align = h_align, v_align = v_align
        )
        return res

    def change_geometry(
        self, height, width, 
        h_align = TextAlignement.H_CENTER, v_align = TextAlignement.V_CENTER 
    ):
        res = self.get_new_geometry( 
            height, width, h_align = h_align, v_align = v_align
        )
        self._lines = res._lines

    def _ascii_art_( self ):
        return ascii_art( self._repr_() )

    def _repr_(self):
        res = ""
        for i in range( len(self._lines)-1 ):
            for c in self._lines[i]:
                res += c
            res += '\n'
        for c in self._lines[-1]:
            res += c
        return res

    def __repr__(self):
        return self._repr_()

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
            default= dict(
                scale=1, line_size=1, point_size=3.5
                , line_color='black', point_color='black'
                , translation=[0,0], rotation=0
            ),
            description='the drawing options',
            checker=lambda x: Set(x.keys()).issubset(
                Set( [
                    'scale', 'line_size', 'point_size'
                    , 'line_color', 'point_color', 'translation', 
                    'rotation'
                ] )
            )
        ),
        ascii_options=dict(
            default= dict(
                row_separator=None, column_separator=None,
                horizontal_padding=None, vertical_padding=None,
                packed=False, compact=False, use_separator=True
            ),
            description='the ascii art options',
            checker=lambda x: Set(x.keys()).issubset(
                Set( [
                    'row_separator', 'column_separator',
                    'horizontal_padding', 'vertical_padding',
                    'packed', 'compact', 'use_separator'
                ] )
            )
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
            return "list repr"

        def _repr_ascii(self):
            return "ascii repr"

        def _ascii_art_( self ):
            ascii_options = self.get_options()['ascii_options']
            array = self.array()
            if 'separator' in dir( self ):
                separator = self.separator()
            else:
                separator = None
            height_array = len( array )
            def to_ascii_array( h, w ):
                if array[h][w] is None:
                    return None
                else:
                    return Ascii_array( ascii_art( array[h][w] ) )
            lines_ascii = [
                [ to_ascii_array(h,w) for w in range( len(array[h]) ) ]
                for h in range(height_array)
            ]
            max_width = 0
            max_height = 0
            width_array = 0
            for line in lines_ascii :
                if len(line) > width_array:
                    width_array = len(line)
                for e in line:
                    if not e is None :
                        if max_width < e.width():
                            max_width = e.width()
                        if max_height < e.height():
                            max_height = e.height()
            if not ascii_options['compact']:
                max_width += 2
                max_height += 0
            print( max_width )
            print( max_height )
            res_height = height_array * ( max_height + 1 ) + 1
            res_width = width_array * ( max_width + 1 ) + 1
            res = Ascii_array( 
                [ [ ' ' for w in range(res_width) ] for h in range(res_height) ]
            )
            for h in range(height_array):
                for w in range( len(array[h]) ):
                    if not separator is None:
                        sep = separator[h][w]
                    if not lines_ascii[h][w] is None :
                        pos_no = [1+h*(1+max_height), 1+w*(1+max_width)]
                        pos_se = [(h+1)*(1+max_height), (w+1)*(1+max_width)]
                        res.clip( lines_ascii[h][w], pos_no, pos_se )
                        if separator is None:
                            sep = [ True, True, True, True ]
                    if sep[0]:
                        for w1 in range( pos_no[1]-1, pos_se[1]+1 ):
                            res[ pos_no[0]-1 ][ w1 ] = '-'
                        res[ pos_no[0]-1 ][ pos_no[1]-1 ] = '+'
                        res[ pos_no[0]-1 ][ pos_se[1] ] = '+'
                    if sep[1]:
                        for h1 in range( pos_no[0]-1, pos_se[0]+1 ):
                            res[ h1 ][ pos_se[1] ] = '|'
                        res[ pos_se[0] ][ pos_se[1] ] = '+'
                        res[ pos_no[0]-1 ][ pos_se[1] ] = '+'
                    if sep[2]:
                        for w1 in range( pos_no[1]-1, pos_se[1]+1 ):
                            res[ pos_se[0] ][ w1 ] = '-'
                        res[ pos_se[0] ][ pos_no[1]-1 ] = '+'
                        res[ pos_se[0] ][ pos_se[1] ] = '+'
                    if sep[3]:
                        for h1 in range( pos_no[0]-1, pos_se[0]+1 ):
                            res[ h1 ][ pos_no[1]-1 ] = '|'
                        res[ pos_no[0]-1 ][ pos_no[1]-1 ] = '+'
                        res[ pos_se[0] ][ pos_no[1]-1 ] = '+'

            return res._ascii_art_()

        def _latex_(self):
            return self.get_options().dispatch(self, '_latex_', 'latex')

        def _latex_drawing(self):
            return "drawing latex"

        def _latex_list(self):
            return "list latex"

