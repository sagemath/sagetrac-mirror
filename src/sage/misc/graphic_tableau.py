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
from sage.all import Reals
from sage.all import Integers

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


def _default_graphic_options( ):
    return GlobalOptions(
        'Graphic Option for tableau',
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
        'Ascii Option for tableau',
        doc='', end_doc='',
        row_separator=dict(default='-', description='A character for the row separators', checker=lambda v: isinstance( v, str ) and len(v)==1 ),
        column_separator=dict(default='|', description='A character for the column separators', checker=lambda v: isinstance( v, str ) and len(v)==1 ),
        corner_separator=dict(default='+', description='A character for the corner separators', checker=lambda v: isinstance( v, str ) and len(v)==1 ),
        horizontal_padding=dict(default=1, description='The size of space beetween the separator and the box contain.', checker=lambda v: v in Integers()),
        vertical_padding=dict(default=0, description='The size of space beetween the separator and the box contain.', checker=lambda v: v in Integers()),
        packed=dict(default=False, description='If set to True, the tableau is packed.', checker=lambda v: isinstance(v, bool )),
        compact=dict(default=False, description='If set to True, the tableau is compacted.', checker=lambda v: isinstance(v, bool )),
        draw_separator=dict(default=True, description='Seprarators is draw.', checker=lambda v: isinstance(v, bool )),
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
            ....:     [], 
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
        """
        self._array = array
        self._separators = separators
        self._entry_map = entry_map
        self.graphic_options = graphic_options
        self.ascii_options = ascii_options

    def _repr_(self):
        return str( self._ascii_art_() )

    def __repr__(self):
        return self._repr_()

    def _ascii_art_( self ):
        height_array = len( self._array )
        def to_ascii_array( h, w ):
            if self._entry_map( self._array[h][w] ) is None:
                return None
            else:
                return Ascii_array( 
                    ascii_art( self._entry_map( self._array[h][w] ) ) 
                )
        lines_ascii = [
            [ to_ascii_array(h,w) for w in range( len(self._array[h]) ) ]
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
        if not self.ascii_options['compact']:
            max_width += 2
            max_height += 0
        res_height = height_array * ( max_height + 1 ) + 1
        res_width = width_array * ( max_width + 1 ) + 1
        res = Ascii_array( 
            [ [ ' ' for w in range(res_width) ] for h in range(res_height) ]
        )
        for h in range(height_array):
            for w in range( len(self._array[h]) ):
                if not self._separators is None:
                    sep = self._separators[h][w]
                else:
                    sep = [ False, False, False, False ]
                if not lines_ascii[h][w] is None :
                    pos_no = [1+h*(1+max_height), 1+w*(1+max_width)]
                    pos_se = [(h+1)*(1+max_height), (w+1)*(1+max_width)]
                    res.clip( lines_ascii[h][w], pos_no, pos_se )
                    if self._separators is None:
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
            pass
