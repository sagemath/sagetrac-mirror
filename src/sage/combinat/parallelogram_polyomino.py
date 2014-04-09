r"""
The Parallelogram Polyominoes 
=====================

The goal of this module is to give some tools to manipulate the 
parallelogram polyominoes.
"""
#*****************************************************************************
#  Copyright (C) 2014 Adrien Boussicault (boussica@labri.fr), 
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element_wrapper import ElementWrapper
from sage.structure.list_clone import ClonableList
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.set_factories import (
    SetFactory, ParentWithSetFactory, TopMostParentPolicy
)
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.global_options import GlobalOptions
from sage.sets.set import Set
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.rings.integer import Integer
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers 
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from copy import deepcopy
from sage.matrix.constructor import matrix

default_tikz_options = dict(
    scale=1, line_size=1, point_size=3
    , color_line='black', color_point='black'
    , color_bounce_0='red', color_bounce_1='blue'
    , translation=[0,0], rotation=0
)

ParallelogramPolyominoesOptions=GlobalOptions(
    name = 'Parallelogram Polyominoes',
    doc=r"""
    """,
    end_doc=r"""
    """,
    tikz_options=dict(
        default= default_tikz_options,
        description='the tikz options',
        checker=lambda x: Set(x.keys()).issubset(
            Set( [
                'scale', 'line_size', 'point_size'
                , 'color_line', 'color_point', 'translation', 
                'rotation', 'color_bounce_0', 'color_bounce_1'
            ] )
        )
    ),
    drawing_components=dict(
        default= dict( diagram=True ),
        description='Different tree-like tableaux components to draw',
        checker=lambda x: Set(x.keys()).issubset(
            Set( [
                'diagram', 'tree', 'bounce_0', 'bounce_1'
            ] )
        )
    ),
    display=dict(
        default="list",
        values= dict(
            list='displayed as list',
            drawing='as a drawing'
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

class ParallelogramPolyomino( ClonableList ):
    r"""
    The class of Parallelogram Polyomino
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        r"""
        """
        return cls._auto_parent._element_constructor_( *args, **opts )

    @lazy_class_attribute
    def _auto_parent(cls):
        r"""
        """
        return ParallelogramPolyominoes()

    def check( self ):
        r"""
        This method raise an error if the internal data of the class doesn't
        represent a parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,0,0,1,0,1,0,1,1], [1,0,1,1,0,0,1,0,0] ]
            ....: )
            sage: pp = ParallelogramPolyomino( [[0,1],[1,0]] )
            sage: pp = ParallelogramPolyomino( [[1],[1]] )
            sage: pp = ParallelogramPolyomino( [[1,0],[0,1]] )
            Traceback (most recent call last):
            ...
            ValueError: Lower and upper path are crossing.
            sage: pp = ParallelogramPolyomino( [[1],[0,1]] )
            Traceback (most recent call last):
            ...
            ValueError: Lower upper paht have different size ( 2 != 1 ).
            sage: pp = ParallelogramPolyomino( [[1],[0]] )
            Traceback (most recent call last):
            ...
            ValueError: The two paths don't join together at the end.
            sage: pp = ParallelogramPolyomino( [[0],[1]] )
            Traceback (most recent call last):
            ...
            ValueError: The two paths don't join together at the end.
            sage: pp = ParallelogramPolyomino( [[0],[0]] )
            Traceback (most recent call last):
            ...
            ValueError: A Parallelogam Polyomino can have the path [[0],[0]].
            sage: pp = ParallelogramPolyomino( [[],[0]] )
            Traceback (most recent call last):
            ...
            ValueError: A Parallelogam Polyomino can have lower or upper path equals to [].
            sage: pp = ParallelogramPolyomino( [[0],[]] )
            Traceback (most recent call last):
            ...
            ValueError: A Parallelogam Polyomino can have lower or upper path equals to [].
            sage: pp = ParallelogramPolyomino( [[],[]] )
            Traceback (most recent call last):
            ...
            ValueError: A Parallelogam Polyomino can have lower or upper path equals to [].
        """
        lower_path = self.lower_path()
        upper_path = self.upper_path()
        if( lower_path==[0] and upper_path==[0] ):
            raise ValueError, "A Parallelogam Polyomino can have the path [[0],[0]]."
        if( lower_path==[] or upper_path==[] ):
            raise ValueError, "A Parallelogam Polyomino can have lower or upper path equals to []."
        if len(upper_path) != len(lower_path):
            raise ValueError, "Lower upper paht have different size ( %s != %s )."%(
                len(upper_path), len(lower_path)
            )
        p_up = [0,0]
        p_down = [0,0]
        for i in range( len(upper_path)-1 ):
            p_up[ 1-upper_path[i] ] += 1
            p_down[ 1-lower_path[i] ] += 1
            if( p_up[0] <= p_down[0] or p_down[1] <= p_up[1] ):
                raise ValueError, "Lower and upper path are crossing."
        p_up[ 1-upper_path[-1] ] += 1
        p_down[ 1-lower_path[-1] ] += 1
        if( p_up[0] != p_down[0] or p_up[1] != p_down[1] ):
            raise ValueError, "The two paths don't join together at the end."

    def __init__(self, parent, value, check=True):
        r"""
        Construct a parallelogram polyomino.
        
        The input is a pair of lower path and upper path.

        The lower and upper path of the empty parallelogram polyomino is
        [1] and [1].

        EXAMPLES::

            sage: lower_path = [0, 0, 1, 0, 1, 1]
            sage: upper_path = [1, 1, 0, 1, 0, 0]
            sage: pp = ParallelogramPolyomino( [ lower_path, upper_path] )
            sage: pp
            [[0, 0, 1, 0, 1, 1], [1, 1, 0, 1, 0, 0]]
            sage: pp = ParallelogramPolyomino( [ [0, 1], [1, 0]] )
            sage: pp
            [[0, 1], [1, 0]]
            sage: pp = ParallelogramPolyomino( [ [1], [1]] )
            sage: pp
            [[1], [1]]
        """
        ClonableList.__init__(self, parent, value )
        if check:
            if not isinstance( value, (list, tuple) ) :
                raise ValueError, "Value %s must be a list or a tuple."%(value)
            self.check()
        self._options = None

    def get_options( self ):
        if self._options is None:
            return self.parent().get_options()
        return self._options

    def set_options( self, *get_value, **set_value ):
        if self._options is None:
            self._options = deepcopy( self.get_options() )
        self._options( *get_value, **set_value )

    def upper_path( self ):
        r"""
        Get the upper path of the parallelogram polyomino
        
        EXAMPLES::
        
            sage: lower_path = [0, 0, 1, 0, 1, 1]
            sage: upper_path = [1, 1, 0, 1, 0, 0]
            sage: pp = ParallelogramPolyomino( [ lower_path, upper_path ] )
            sage: pp.upper_path()
            [1, 1, 0, 1, 0, 0]
        """
        return list( ClonableList.__getitem__( self, 1 ) )

    def lower_path( self ):
        r"""
        Get the lower path of the parallelogram polyomino
        
        EXAMPLES::
        
            sage: lower_path = [0, 0, 1, 0, 1, 1]
            sage: upper_path = [1, 1, 0, 1, 0, 0]
            sage: pp = ParallelogramPolyomino( [ lower_path, upper_path ] )
            sage: pp.lower_path()
            [0, 0, 1, 0, 1, 1]
        """
        return list( ClonableList.__getitem__( self, 0 ) )

    def _path_heights( self, word, up ):
        res = []
        h = 0
        for e in word:
            if e == up :
                h += 1
            else:
                res.append( h )
        return res

    def upper_heights( self ):
        return self._path_heights( self.upper_path(), 0 )

    def lower_heights( self ):
        return self._path_heights( self.lower_path(), 0 )

    def upper_widths( self ):
        return self._path_heights( self.upper_path(), 1 )

    def lower_widths( self ):
        return self._path_heights( self.lower_path(), 1 )

    def widths( self ):
        r"""
        This method return a list of the withds of the parallelogram 
        polyomino : the parallelogram polyomino is splited row by row and the
        method return a the list containing the sizes of the rows.

        Examples::
        
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,0,0,1,0,1,0,1,1], [1,0,1,1,0,0,1,0,0] ]
            ....: )
            sage: pp.widths()
            [1, 3, 3, 3, 2]
            sage: pp = ParallelogramPolyomino( [ [0, 1], [1, 0] ] )
            sage: pp.widths()
            [1]
            sage: pp = ParallelogramPolyomino( [ [1], [1] ] )
            sage: pp.widths()
            []
        """
        widths = []
        uw = self.upper_widths()
        lw = self.lower_widths()
        for i in range( len(lw) ):
            widths.append( uw[i] - lw[i] )
        return widths

    def heights( self ):
        r"""
        This method return a list of heights of the parallelogram 
        polyomino : the parallelogram polyomino is splited column by column and
        the method return a the list containing the sizes of the rows.

        Examples::
        
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,0,0,1,0,1,0,1,1], [1,0,1,1,0,0,1,0,0] ]
            ....: )
            sage: pp.heights()
            [3, 3, 4, 2]
            sage: pp = ParallelogramPolyomino( [ [0, 1], [1, 0] ] )
            sage: pp.heights()
            [1]
            sage: pp = ParallelogramPolyomino( [ [1], [1] ] )
            sage: pp.heights()
            [0]
        """
        heights = []
        uh = self.upper_heights()
        lh = self.lower_heights()
        for i in range( len(uh) ):
            heights.append( lh[i] - uh[i] )
        return heights

    def width( self ):
        r"""
        Return the width of the parallelogram polyomino.

        EXAMPLES::
        
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,1,0,0,1,1,0,1,1,1], [1,1,1,0,1,0,0,1,1,0] ]
            ....: )
            sage: pp.width()
            6
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,1], [1,0] ]
            ....: )
            sage: pp.width()
            1
            sage: pp = ParallelogramPolyomino(
            ....:     [ [1], [1] ]
            ....: )
            sage: pp.width()
            1
        """
        if( self.size()==0 ):
            return 1
        return self.upper_widths()[-1]

    def height( self ):
        r"""
        Return the height of the parallelogram polyomino.

        EXAMPLES::
        
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,1,0,0,1,1,0,1,1,1], [1,1,1,0,1,0,0,1,1,0] ]
            ....: )
            sage: pp.height()
            4
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,1], [1,0] ]
            ....: )
            sage: pp.height()
            1
            sage: pp = ParallelogramPolyomino(
            ....:     [ [1], [1] ]
            ....: )
            sage: pp.height()
            0
        """
        if( self.size()==0 ):
            return 0
        return self.lower_heights()[-1]

    @cached_method
    def get_array( self ):
        r"""
            Return an array of 0 and 1 such that the represent the boxes of
            the parallelogram polyomino
            
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,0,0,0,1,0,1,0,1], [1,0,0,0,1,1,0,0,0] ]
            ....: )
            sage: matrix( pp.get_array() )
            [1 0 0]
            [1 0 0]
            [1 0 0]
            [1 1 1]
            [0 1 1]
            [0 0 1]
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,1], [1,0] ]
            ....: )
            sage: pp.get_array()
            [[1]]
            sage: pp = ParallelogramPolyomino(
            ....:     [ [1], [1] ]
            ....: )
            sage: pp.get_array()
            []
        """
        width = self.width()
        height = self.height()
        lower_widths = self.lower_widths()
        widths = self.widths()
        def val( w, h ):
            if w >= len( widths ) or w < 0:
                return 0
            if lower_widths[w] <= h and h < lower_widths[w]+widths[w]:
                return 1
            return 0
        return [ [ val(h,w) for w in range(width) ] for h in range(height) ]

    class _polyomino_row:
        def __init__( self, polyomino, row ):
            self.polyomino = polyomino
            self.row = row
        def __getitem__( self, column ):
            if( 
                self.is_inside()
                and 0<= column and column < self.polyomino.width()
            ):
                return self.polyomino.get_array()[ self.row ][ column ]
            return 0
        def is_inside( self ):
            return 0 <= self.row and self.row < self.polyomino.height()
        def is_outside( self ):
            return not( self.is_inside() )
        def __repr__( self ):
            if self.is_outside():
                return "The (outside) row %s of the parallelogram"%(self.row)
            else:
                return "The row %s of the parallelogram polyomino"%(self.row)

    def __getitem__(self, row):
        r"""
            Return the row of the parallelogram polyomino.
            The index of the tow can be out of range of the height of 
            the parallelogram polyomino.
            In that case, the row returned is outside the parallelogram
            polyomino.

            EXAMPLES::
                sage: pp = ParallelogramPolyomino(
                ....:     [ [0,0,1,0,1,1], [1,1,0,0,1,0] ]
                ....: )
                sage: pp[0].is_inside()
                True
                sage: pp[3].is_outside()
                True
                sage: pp[-1].is_outside()
                True
                sage: pp[0][0]
                1
                sage: pp[ [0,1] ]
                1
                sage: pp[2][0]
                0
                sage: pp[-1][0]
                0
                sage: pp[ [4,4] ]
                0
        """
        if isinstance( row, list ):
            return self._polyomino_row( self, row[0] )[ row[1] ]
        return self._polyomino_row( self, row )

    def bounce_path(self, direction = 1 ):
        r"""
            Return the bounce path of the parallelogram polyomino.

            The bounce is a path with two steps (1,0) and (0,1).

            If 'direction' is 1 (resp. 0), the bounce path is the path
            starting at position position (h=1,w=0) (resp. (h=0,w=1)) with 
            initial direction, the vector (0,1) (resp. (1,0)), and turning 
            each time the path crosses the permiter of the parallelogram
            polyomino.

            The path is coded by a list of integer. Each integer represent
            the size of the path beetween two turnings.

            You can visualize the two bounce paths by using the following 
            commands :
            
            EXAMPLES::
                sage: PP = ParallelogramPolyomino(
                ....:     [ [0,0,1,0,1,1], [1,1,0,0,1,0] ]
                ....: )
                sage: PP.bounce_path( direction=1)
                [2, 2, 1]
                sage: PP.bounce_path( direction=0)
                [2, 1, 1, 1]

                sage: PP = ParallelogramPolyomino(
                ....:     [ [0,0,1,1,1,0,0,1,1], [1,1,1,0,1,1,0,0,0] ]
                ....: )
                sage: PP.bounce_path( direction=1)
                [3, 1, 2, 2]
                sage: PP.bounce_path( direction=0)
                [2, 4, 2]

                sage: #PP = ParallelogramPolyomino(
                ....: #    [ [0,0,1,0,1,1], [1,1,0,0,1,0] ]
                ....: #)
                sage: #PP.set_options(
                ....: #    drawing_components=dict(
                ....: #        diagram = True
                ....: #        ,bounce_0 = True
                ....: #        , bounce_1 = True
                ....: #    )
                ....: #)
                sage: #view( PP )

                sage: PP = ParallelogramPolyomino(
                ....:     [ [0,1], [1,0] ]
                ....: )
                sage: PP.bounce_path( direction=1)
                [1]
                sage: PP.bounce_path( direction=0)
                [1]

                sage: PP = ParallelogramPolyomino(
                ....:     [ [1], [1] ]
                ....: )
                sage: PP.bounce_path( direction=1)
                []
                sage: PP.bounce_path( direction=0)
                []

        """
        result = []
        pos = [0,0]
        pos[direction] -= 1
        old = list( pos ) 
        ne = list( pos )
        ne[direction] += 1
        while( self[ ne ] == 1 ):
            pos[direction]+=1
            while( self[pos] == 1  ):
                pos[direction]+=1
            pos[direction] -= 1
            result.append( pos[direction]-old[direction] )
            direction = 1 - direction
            old[0], old[1] = pos
            ne[0], ne[1] = pos
            ne[ direction ] += 1
        return result

    def bounce( self, direction = 1 ):
        r"""
            Return the bounce of the parallelogram polyomino.

            Les p be the bounce path of the parallelogram 
            polyomino. ( p=self.bounce_path() )
            The bounce is defined by :
            sum( [ (1+ floor(i/2))*p[i] for i in range(len(p)) ] )

            EXAMPLES::

                sage: PP = ParallelogramPolyomino(
                ....:     [ [0,0,1,0,1,1], [1,1,0,0,1,0] ]
                ....: )
                sage: PP.bounce( direction=1 )
                6
                sage: PP.bounce( direction=0)
                7

                sage: PP = ParallelogramPolyomino(
                ....:     [ [0,0,1,1,1,0,0,1,1], [1,1,1,0,1,1,0,0,0] ]
                ....: )
                sage: PP.bounce( direction=1)
                12
                sage: PP.bounce( direction=0)
                10

                sage: PP = ParallelogramPolyomino(
                ....:     [ [0,1], [1,0] ]
                ....: )
                sage: PP.bounce( direction=1)
                1
                sage: PP.bounce( direction=0)
                1

                sage: PP = ParallelogramPolyomino(
                ....:     [ [1], [1] ]
                ....: )
                sage: PP.bounce( direction=1)
                0
                sage: PP.bounce( direction=0)
                0
        """
        result = 0
        path = self.bounce_path( direction )
        for i in range( len(path) ):
            result += (1+(int(i/2)))*path[i]
        return result

    def area( self ):
        r"""
            Returns the are of the parallelogram polyomino.
            
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,0,1,0,1,1,1,0,0,1,1], [1,1,0,1,0,1,1,0,1,0,0] ]
            ....: )
            sage: pp.area()
            13
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,1], [1,0] ]
            ....: )
            sage: pp.area()
            1
            sage: pp = ParallelogramPolyomino(
            ....:     [ [1], [1] ]
            ....: )
            sage: pp.area()
            0
        """
        res = 0
        for h in self.heights():
            res += h
        return res

    def _repr_(self):
        r"""
        Return a string representation of the parallelogram polyomino.

        EXAMPLES::
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,0,1,0,1,1], [1,1,0,0,1,0] ]
            ....: )
            sage: pp
            [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            sage: pp.set_options( display='drawing' )
            sage: pp
            [1 1 0]
            [1 1 0]
            [0 1 1]
        """
        return self.get_options().dispatch(self, '_repr_', 'display')

    def _repr_list( self ):
        r"""
        Return a string representation with list style.

        EXAMPLES::
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,0,1,0,1,1], [1,1,0,0,1,0] ]
            ....: )
            sage: pp._repr_list()
            '[[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]'
        """
        return ClonableList._repr_(self)

    def _repr_drawing( self ):
        r"""
        Return a string representing a drawing od the parallelogram polyomino.

        EXAMPLES::
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,0,1,0,1,1], [1,1,0,0,1,0] ]
            ....: )
            sage: pp._repr_drawing()
            '[1 1 0]\n[1 1 0]\n[0 1 1]'
        """
        return str( matrix( self.get_array() ) )

    def get_tikz_options( self ):
        return self.get_options()['tikz_options']

    def _to_tikz_diagram( self ):
        tikz_options = self.get_tikz_options()
        grid_width = self.width() + 1
        grid_height = self.height() + 1
        def X( x ):
            return x
        def Y( y ):
            return grid_height-1-y
        res = ""
        if self.size() == 0:
            res += "\n  \\draw[color_line=%s, line width=%s] (%s, %s) -- (%s,%s);"%(
                tikz_options['color'], tikz_options['line_size'],
                X(0),Y(0),
                X(1),Y(0)
            )
            return res 
        res += "\n  \\draw[color=%s, line width=%s] (%s, %s) -- (%s,%s);"%(
            tikz_options['color_line'], tikz_options['line_size'],
            X(0),Y(0),
            X(0),Y(self.lower_heights()[0])
        )
        res += "\n  \\draw[color=%s, line width=%s] (%s, %s) -- (%s,%s);"%(
            tikz_options['color_line'], tikz_options['line_size'],
            X(grid_width-1), Y(self.upper_heights()[ grid_width-2 ]),
            X(grid_width-1), Y(self.lower_heights()[ grid_width-2 ])
        )
        res += "\n  \\draw[color=%s, line width=%s] (%s, %s) -- (%s,%s);"%(
            tikz_options['color_line'], tikz_options['line_size'],
            X(0), Y(0), X(self.upper_widths()[0]), Y(0)
        )
        res += "\n  \\draw[color=%s, line width=%s] (%s, %s) -- (%s,%s);"%(
            tikz_options['color_line'], tikz_options['line_size'],
            X(self.lower_widths()[ grid_height-2 ]),
            Y(grid_height -1),
            X(self.upper_widths()[ grid_height-2 ]),
            Y(grid_height -1)
        )
        for w in xrange( 1, grid_width-1 ):
            h1 = self.upper_heights()[w-1]
            h2 = self.lower_heights()[w]
            res += "\n  \\draw[color=%s, line width=%s] (%s, %s) -- (%s,%s);"%(
                tikz_options['color_line'], tikz_options['line_size'],
                X(w),Y(h1), X(w),Y(h2)
            )
        for h in xrange( 1, grid_height-1 ):
            w1 = self.lower_widths()[h-1]
            w2 = self.upper_widths()[h]
            res += "\n  \\draw[color=%s, line width=%s] (%s, %s) -- (%s,%s);"%(
                tikz_options['color_line'], tikz_options['line_size'],
                X(w1),Y(h), X(w2),Y(h)
            )
        return res

    def to_tikz( self ):
        r"""
        Return the tikz code of the parallelogram polyomino.

        This code is the code present inside a tikz latex environemet.
        """
        res = ""
        drawing_components = self.get_options()['drawing_components']
        print drawing_components
        if 'diagram' in  drawing_components :
            res += self._to_tikz_diagram()
        return res

    def geometry(self):
        r"""
            Returns a pair [h,w] containing the height and the width of the 
            parallelogram polyomino

            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,1,1,1,1], [1,1,1,1,0] ] 
            ....: )
            sage: pp.geometry()
            [1, 4]
            sage: pp = ParallelogramPolyomino(
            ....:     [ [0,1], [1,0] ] 
            ....: )
            sage: pp.geometry()
            [1, 1]
            sage: pp = ParallelogramPolyomino(
            ....:     [ [1], [1] ] 
            ....: )
            sage: pp.geometry()
            [0, 1]
        """
        return [ self.height(), self.width() ]

    def size( self ):
        r"""
        Returns the size of the parallelogram polyomino.
        
        sage: pp = ParallelogramPolyomino(
        ....:     [ [0,0,0,0,1,0,1,1], [1,0,0,0,1,1,0,0] ] 
        ....: )
        sage: pp.size()
        7
        sage: pp = ParallelogramPolyomino(
        ....:     [ [0,1], [1,0] ] 
        ....: )
        sage: pp.size()
        1
        sage: pp = ParallelogramPolyomino(
        ....:     [ [1], [1] ] 
        ....: )
        sage: pp.size()
        0
        """
        return len( self.upper_path() )-1

    def _latex_(self):
        r"""
        Return a LaTeX version of ``self``.

        For more on the latex options, see 
        :meth:`ParallelogramPolyominoes.global_options`.
        """
        return self.get_options().dispatch(self, '_latex_', 'latex')

    def _latex_drawing( self ):
        r"""
        Return a LaTeX version of ``self`` in a drawing style.
        """
        latex.add_package_to_preamble_if_available("tikz")
        tikz_options = self.get_tikz_options()
        res = "\n\\begin{tikzpicture}[scale=%s]"%(tikz_options['scale'])
        res += self.to_tikz()
        res += "\n\\end{tikzpicture}"
        return res

    def _latex_list( self ):
        r"""
        Return a LaTeX version of ``self`` in a list style.
        """
        return "\\[%s\\]"%(self._repr_list())
        NotImplemented


class ParallelogramPolyominoesFactory(SetFactory):
    r"""
    The parallelogram polyominoes factory.
    """
    def __call__(self, size=None, policy=None):
        r"""
        """
        if policy is None:
            policy = self._default_policy

        if isinstance(size, (Integer, int)):
            return ParallelogramPolyominoes_size(size, policy)
        if size is None:
            return ParallelogramPolyominoes_all(policy)
        raise ValueError, "Invalide argument for Tre-like tableaux Factory."

    def add_constraints(self, cons, (args, opts)):
        r"""
        """
        return cons+args

    @lazy_attribute
    def _default_policy(self):
        return TopMostParentPolicy(self, (), ParallelogramPolyomino)

    def _repr_(self):
        """
        """
        return "Factory for parallelogram polyominoes"

ParallelogramPolyominoes = ParallelogramPolyominoesFactory()
ParallelogramPolyominoes.__doc__ = ParallelogramPolyominoesFactory.__call__.__doc__


class ParallelogramPolyominoes_size(ParentWithSetFactory, UniqueRepresentation):
    r"""
    The parallelogram polyominoes of size `n`.
    """
    def __init__(self, size, policy):
        r"""
        """
        self._size = size
        ParentWithSetFactory.__init__(
            self, (size,), policy, category = FiniteEnumeratedSets()
        )

    def _repr_(self):
        r"""
        """
        return "parallelogram polyominoes of size %s"%(self._size)

    def an_element(self):
        r"""
        """
        return self.first()

    def first( self ):
        return next( self.__iter__() )

    def check_element(self, el, check):
        r"""
        """
        if el.size() != self._size:
            raise ValueError, "The parallelogram polyomino have a Wrong size : %s"%(el.size())

    def __iter__(self):
        r"""
        """
        NotImplemented

    def get_options( self ):
        return self.global_options

    global_options = ParallelogramPolyominoesOptions



class ParallelogramPolyominoes_all( ParentWithSetFactory, DisjointUnionEnumeratedSets ):
    r"""
    This class enumerate all the parallelogram polyominoes.
    """
    def __init__(self, policy):
        r"""
        """
        ParentWithSetFactory.__init__(
            self, (), policy, category = FiniteEnumeratedSets()
        )
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), ParallelogramPolyominoes_size),
            facade=True, keepkey = False,
            category = self.category()
        )

    def _repr_(self):
        r"""
        """
        return "ParallelogramPolyominoes_all"

    def check_element(self, el, check):
        r"""
        """
        pass

    def get_options( self ):
        return self.global_options

    global_options = ParallelogramPolyominoesOptions
