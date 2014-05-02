"""
This file define some tools to obtain ascii and latex outputs for any tableaux.
"""
#*****************************************************************************
#       Copyright (C) 2014 Adrien Boussicault <boussica@labri.fr>
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
from sage.all import Reals
from sage.all import Integers

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
        h_align = 'center', v_align = 'center'
    ):
        clip_width = pos_se[1] - pos_no[1]
        clip_height = pos_se[0] - pos_no[0]

        if not isinstance( array, Ascii_array ):
            array = Ascii_array( array  )
        array_width = array.width()
        array_height = array.height()

        if (v_align == 'center') :
            def H( h ):
                return int(array_height/2) - int(clip_height/2) + (h - pos_no[0])
        elif v_align == 'top' :
            def H( h ):
                return h - pos_no[0]
        elif v_align == 'bottom' :
            def H( h ):
                return array_height - clip_height + (h - pos_no[0])

        if h_align == 'center' :
            def W( w ):
                return int(array_width/2) - int(clip_width/2) + (w - pos_no[1])
        elif h_align == 'left' :
            def W( w ):
                return w - pos_no[1]
        elif h_align == 'right' :
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
        h_align = 'center', v_align = 'center' 
    ):
        res = Ascii_array( [ [ ' ' for w in range(width) ] for h in range(height)] )
        res.clip(
            self, [0,0], [width, height],
            h_align = h_align, v_align = v_align
        )
        return res

    def change_geometry(
        self, height, width, 
        h_align = 'center', v_align = 'center' 
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


def _default_graphic_options( ):
    return GlobalOptions(
        'grahic representation of tableau',
        doc='', end_doc='',
        scale=dict(default=10, description='The scale of the figure', checker=lambda sc: sc in Reals() and sc>0 ),
        line_size=dict(default=1, description='The size of the lines of th tableau', checker=lambda sc: sc in Reals() and sc>0),
        point_size=dict(default=3.5, description='The size of the points inside boxes', checker=lambda sc: sc in Reals() and sc>0),
        line_color=dict(default='black', description='The color of the lines of th tableau', checker=lambda v: isinstance( v, str )),
        point_color=dict(default='black', description='The color of the points inside boxes', checker=lambda v: isinstance( v, str )),
        translation=dict(default=[0,0], description='The translation', checker=lambda v: isinstance( v, list ) and len(v)==2 ),
        rotation=dict(default=0, description='The rotation with center is O.', checker=lambda sc: sc in Reals() )
    )

def _default_ascii_options( ):
    return GlobalOptions(
        'ascii representation of tableau',
        doc='', end_doc='',
        row_separator=dict(default='-', description='A character for the row separators', checker=lambda v: isinstance( v, str ) and len(v)==1 ),
        column_separator=dict(default='|', description='A character for the column separators', checker=lambda v: isinstance( v, str ) and len(v)==1 ),
        corner_separator=dict(default='+', description='A character for the corner separators', checker=lambda v: isinstance( v, str ) and len(v)==1 ),
        horizontal_padding=dict(default=1, description='The size of space beetween the separator and the box contain.', checker=lambda v: v in Integers()),
        vertical_padding=dict(default=0, description='The size of space beetween the separator and the box contain.', checker=lambda v: v in Integers()),
        packed=dict(default=False, description='If set to True, the tableau is packed.', checker=lambda v: isinstance(v, bool )),
        compact=dict(default=False, description='If set to True, the tableau is compacted.', checker=lambda v: isinstance(v, bool )),
        draw_separator=dict(default=True, description='Seprarators is draw.', checker=lambda v: isinstance(v, bool )),
        horizontal_alignement=dict(default='center', values=dict( center='aligned to center', left='aligned to left', right='aligned to right'  )),
        vertical_alignement=dict(default='center', values=dict( center='aligned to center', top='aligned to top', bottom='aligned to bottom'  )),
    )

class GraphicTableau:

    @staticmethod
    def default_graphic_options( ):
        return _default_graphic_options( )

    @staticmethod
    def default_ascii_options( ):
        return _default_ascii_options( )

    def __init__(
        self, array, separators = None,
        entry_map = lambda x : x,
        graphic_options =  _default_graphic_options(),
        ascii_options = _default_ascii_options()
    ):
        r"""
            TODO

            EXAMPLES::

                sage: from sage.misc.graphic_tableau import GraphicTableau
                sage: g = GraphicTableau( [
                ....:     [ None, None, 'S', 'A', 'G',  'E' ],
                ....:     [ 'C', 'O', 'M', 'B', 'I', 'N', 'A', 'T'],
                ....:     [ ], 
                ....:     [ None, Partitions(4)[0], Partitions(4)[1], Partitions(4)[2],  Partitions(4)[3], Partitions(4)[4]  ] 
                ....: ] )
                sage: ascii_art( g )
                              +------+------+------+------+              
                              |      |      |      |      |              
                              |      |      |      |      |              
                              |   S  |   A  |   G  |   E  |              
                              |      |      |      |      |              
                +------+------+------+------+------+------+------+------+
                |      |      |      |      |      |      |      |      |
                |      |      |      |      |      |      |      |      |
                |   C  |   O  |   M  |   B  |   I  |   N  |   A  |   T  |
                |      |      |      |      |      |      |      |      |
                +------+------+------+------+------+------+------+------+
                <BLANKLINE>
                <BLANKLINE>
                <BLANKLINE>
                <BLANKLINE>
                       +------+------+------+------+------+              
                       |      |      |      |      |   *  |              
                       |      |  *** |  **  |  **  |   *  |              
                       | **** |  *   |  **  |  *   |   *  |              
                       |      |      |      |  *   |   *  |              
                       +------+------+------+------+------+              

                sage: g = GraphicTableau( [
                ....:     [], 
                ....:     [None, 'B'],
                ....:     [None], 
                ....: ] )
                sage: ascii_art( g )
                <BLANKLINE>
                <BLANKLINE>
                   +---+
                   | B |
                   +---+
                <BLANKLINE>
                <BLANKLINE>
                sage: ascii_art( g )._repr_()
                '        \n        \n   +---+\n   | B |\n   +---+\n        \n        '

                sage: options = GraphicTableau.default_ascii_options()
                sage: options( compact=True )
                sage: g = GraphicTableau(
                ....:     [[tree for tree in BinaryTrees(i)] for i in range(4)],
                ....:     ascii_options = options
                ....: )
                sage: ascii_art( g )
                +-----+                        
                |     |                        
                |     |                        
                |     |                        
                |     |                        
                |     |                        
                +-----+                        
                |     |                        
                |     |                        
                |  o  |                        
                |     |                        
                |     |                        
                +-----+-----+                  
                |     |     |                  
                | o   |   o |                  
                |  \  |  /  |                  
                |   o | o   |                  
                |     |     |                  
                +-----+-----+-----+-----+-----+
                |o    | o   |     |   o |    o|
                | \   |  \  |  o  |  /  |   / |
                |  o  |   o | / \ | o   |  o  |
                |   \ |  /  |o   o|  \  | /   |
                |    o| o   |     |   o |o    |
                +-----+-----+-----+-----+-----+
        """

        def maximal_dimension( array ):
            if array is None:
                return [0,0]
            width = 0
            for line in array :
                if len(line) > width:
                    width = len(line)
            return [ len(array), width ]

        def normalized_array( array, height, width, fct=None ):
            def get_entry( array, h , w, fct ):
                if h < len( array ) and w < len( array[h] ):
                    if fct is None:
                        return array[h][w]
                    else:
                        return fct( array[h][w] )
                else:
                    return None
            return [
                [ get_entry( array, h, w, entry_map ) for w in range( width) ]
                for h in range( height )
            ]

        height_array, width_array = maximal_dimension( array )
        height_separators, width_separators = maximal_dimension( separators )
        height = max( height_array, height_separators )
        width = max( width_array, width_separators )

        self._array = normalized_array( array, height, width, entry_map )

        def normalized_separators( separators, array, height, width ):
            def get_entry( h , w ):
                if separators is None:
                    if not array[h][w] is None:
                        return [True, True, True, True]
                else:
                    if h < len( separators ) and w < len( separators[h] ):
                        return separators[h][w]
                return [False, False, False, False]
            return [
                [ get_entry( h, w ) for w in range( width) ]
                for h in range( height )
            ]
        self._separators = normalized_separators(
            separators, self._array, height, width
        )
        self.graphic_options = graphic_options
        self.ascii_options = ascii_options

    def __getitem__( self, h ):
        class Line :
            def __init__( self, array, h ):
                self._array = array
                self._h = h
            def __getitem__( self, w ):
                return self._array[self._h][w]
            def __repr__( self ):
                return str( self._array[h] )
        return Line( self._array, h )

    def _repr_(self):
        return str( self._ascii_art_() )

    def __repr__(self):
        return self._repr_()

    def array( self ):
        return self._array

    def separators( self ):
        return self._separators

    def height( self ):
        return len( self._array )

    def width( self ):
        return len( self._array[0] )

    def get_ascii_entry( self, h, w ):
        return Ascii_array( ascii_art( self[h][w] ) )

    def dimensions( self ):
        heights = [ 0 for h in range( self.height() ) ]
        widths = [ 0 for w in range( self.width() ) ]

        w_offset = 0
        h_offset = 0
        if not self.ascii_options['compact']:
            w_offset += 2
            h_offset += 0

        for h in range( self.height() ):
            for w in range( self.width() ):
                if not self[h][w] is None:
                    figure = self.get_ascii_entry( h, w )
                    heights[h] = max( figure.height() + h_offset, heights[h] )
                    widths[w] = max( figure.width() + w_offset, widths[w] )
        return [heights, widths]

    def max_of_dimensions( self ):
        heights, widths = self.dimensions()
        return [ max(heights), max(widths) ]

    def picture_dimension( self ):
        max_height, max_width = self.max_of_dimensions()
        height_picture = self.height() * ( max_height + 1 ) + 1
        width_picture = self.width() * ( max_width + 1 ) + 1
        return [ height_picture, width_picture ]

    def horizontal_alignement( self ):
        return self.ascii_options['horizontal_alignement']

    def vertical_alignement( self ):
        return self.ascii_options['vertical_alignement']

    def _ascii_art_( self ):
        height_picture, width_picture = self.picture_dimension()

        picture = Ascii_array( [
            [ ' ' for w in range(width_picture) ] for h in range(height_picture)
        ] )

        max_height, max_width = self.max_of_dimensions()

        def draw_a_box( h, w ):
            # we draw inside the box
            if not self[h][w] is None :
                pos_no = [1+h*(1+max_height), 1+w*(1+max_width)]
                pos_se = [(h+1)*(1+max_height), (w+1)*(1+max_width)]
                picture.clip(
                    self.get_ascii_entry( h, w ), pos_no, pos_se,
                    self.horizontal_alignement(), self.vertical_alignement()
                )

            #we draw the borders

            row_sep = self.ascii_options['row_separator']
            col_sep = self.ascii_options['column_separator']
            corner_sep = self.ascii_options['corner_separator']
            if self.separators()[h][w][0]:
                for w1 in range( pos_no[1]-1, pos_se[1]+1 ):
                    picture[ pos_no[0]-1 ][ w1 ] = row_sep
                picture[ pos_no[0]-1 ][ pos_no[1]-1 ] = corner_sep
                picture[ pos_no[0]-1 ][ pos_se[1] ] = corner_sep
            if self.separators()[h][w][1]:
                for h1 in range( pos_no[0]-1, pos_se[0]+1 ):
                    picture[ h1 ][ pos_se[1] ] = col_sep
                picture[ pos_se[0] ][ pos_se[1] ] = corner_sep
                picture[ pos_no[0]-1 ][ pos_se[1] ] = corner_sep
            if self.separators()[h][w][2]:
                for w1 in range( pos_no[1]-1, pos_se[1]+1 ):
                    picture[ pos_se[0] ][ w1 ] = row_sep
                picture[ pos_se[0] ][ pos_no[1]-1 ] = corner_sep
                picture[ pos_se[0] ][ pos_se[1] ] = corner_sep
            if self.separators()[h][w][3]:
                for h1 in range( pos_no[0]-1, pos_se[0]+1 ):
                    picture[ h1 ][ pos_no[1]-1 ] = col_sep
                picture[ pos_no[0]-1 ][ pos_no[1]-1 ] = corner_sep
                picture[ pos_se[0] ][ pos_no[1]-1 ] = corner_sep

        # We draw all the boxes
        for h in range( self.height() ):
            for w in range( self.width() ):
                draw_a_box( h, w )


        return picture._ascii_art_()

    def _latex_(self):
        pass
