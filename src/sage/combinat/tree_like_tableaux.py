r"""
The Tree-like tableaux
=====================

The goal of this module is to give some tools to manipulate the
tree-like tableaux.

TO SAY:
    The code is based on tableau.py,
"""
#*****************************************************************************
#  Copyright (C) 2014  Patxi Laborde Zubieta (patxi.laborde.zubieta@gmail.com),
#
#          co-auteur : Adrien Boussicault (boussica@labri.fr)
#                        - Tikz options
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.all import factorial
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.family import Family
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.structure.list_clone import ClonableList
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.set import Set
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.set_factories import (
    SetFactory, ParentWithSetFactory, TopMostParentPolicy
)
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.global_options import GlobalOptions
from sage.combinat.partition import Partition
from sage.rings.integer import Integer
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.misc.cachefunc import cached_method
from sage.combinat.permutation import (Permutation, Permutations)
from sage.combinat.binary_tree import (BinaryTree, BinaryTrees, LabelledBinaryTree, LabelledBinaryTrees)

#from sage.combinat.non_ambiguous_tree import *

from copy import deepcopy


default_tikz_options = dict(
    scale=1, line_size=1, ribbon_size=1, point_size=3,
    color_array='black', color_ribbon='red', color_point='black'
    , color_special_point='black', color_tabular='black'
    , translation=[0,0], rotation=0
)

TreeLikeTableauxOptions=GlobalOptions(
    name = 'tree-like tableaux',
    doc=r"""
    """,
    end_doc=r"""
    """,
    tikz_options=dict(
        default= default_tikz_options,
        description='the tikz options',
        checker=lambda x: Set(x.keys()).issubset(
            Set( [
                'scale', 'line_size', 'ribbon_size', 'point_size', 'color_array'
                , 'color_ribbon', 'color_point', 'color_special_point'
                , 'color_tabular', 'translation', 'rotation'
            ] )
        )
    ),
    drawing_components=dict(
        default= dict( diagram=True ),
        description='Different tree-like tableaux components to draw',
        checker=lambda x: Set(x.keys()).issubset(
            Set( [
                'diagram', 'tree', 'insertion_history', 'ribbons', '2+12'
                , '1+21'
            ] )
        )
    ),
    display=dict(
        default="list",
        values= dict(
            list='displayed as list',
            drawing='as a drawing',
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

class _drawing_tool:
    r"""
    TODO
    Technical class to produce TIKZ drawing.

    This class contains some 2D geometric tools to produce some TIKZ 
    drawings.
    With that classes you can use options to set up drawing informations.
    The the class will produce a drawin by using those informations.

    EXAMPLES::
        sage: from sage.combinat.parallelogram_polyomino import (
        ....:     _drawing_tool, default_tikz_options,
        ....:     ParallelogramPolyominoesOptions
        ....: )
        sage: opt = ParallelogramPolyominoesOptions['tikz_options']
        sage: dt = _drawing_tool(opt)
        sage: dt.draw_line([1, 1], [-1, -1])
        '\n  \\draw[color=black, line width=1] (1.000000, 1.000000) -- (-1.000000, -1.000000);'

        sage: fct = lambda vec: [2*vec[0], vec[1]]
        sage: dt = _drawing_tool(opt, fct)
        sage: dt.draw_line([1, 1], [-1, -1])
        '\n  \\draw[color=black, line width=1] (2.000000, 1.000000) -- (-2.000000, -1.000000);'

        sage: import copy
        sage: opt = copy.deepcopy(opt)
        sage: opt['mirror'] = [0,1]
        sage: dt = _drawing_tool(opt)
        sage: dt.draw_line([1, 1], [-1, -1])
        '\n  \\draw[color=black, line width=1] (-1.000000, 1.000000) -- (1.000000, -1.000000);'

    """
    def __init__(self, options, XY=lambda v: v):
        self._XY = lambda v: XY([float(v[0]), float(v[1])])
        self._translation = options['translation']
        self._mirror = options['mirror']
        self._rotation = options['rotation']
        self._color_line = options['color_line']
        self._line_size = options['line_size']
        self._point_size = options['point_size']
        self._color_point = options['color_point']

    def XY(self, v):
        """
        This function give the image of v by some transformation given by the 
        drawingoption of ``_drawing_tool``.

        The transformation is the composition of rotation, mirror, translation 
        and XY user function.
        First we apply XY function, then the translation, then the mirror and 
        finaly the rotation.

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     _drawing_tool, ParallelogramPolyominoesOptions
            ....: )
            sage: opt = ParallelogramPolyominoesOptions['tikz_options']
            sage: dt = _drawing_tool(opt)
            sage: dt.XY( [1, 1] )
            [1.0, 1.0]

            sage: fct = lambda vec: [2*vec[0], vec[1]]
            sage: dt = _drawing_tool(opt, fct)
            sage: dt.XY([1, 1])
            [2.0, 1.0]

            sage: import copy
            sage: opt = copy.deepcopy(opt)
            sage: opt['mirror'] = [0, 1]
            sage: dt = _drawing_tool(opt)
            sage: dt.XY([1, 1])
            [-1.0, 1.0]
        """
        def translate(pos, v):
            return [pos[0]+v[0], pos[1]+v[1]]

        def rotate(pos, angle):
            [x, y] = pos
            return [x*cos(angle) - y*sin(angle), x*sin(angle) + y*cos(angle)]

        def mirror(pos,  axe):
            if axe is None:
                return pos
            if not isinstance(axe, (list, tuple)):
                raise ValueError(
                    "mirror option should be None or a list of two real" +
                    " encoding a 2D vector."
                )
            n = float(sqrt(axe[0]**2 + axe[1]**2))
            axe[0] = float(axe[0]/n)
            axe[1] = float(axe[1]/n)
            sp = (pos[0]*axe[0] + pos[1]*axe[1])
            sn = (- pos[0]*axe[1] + pos[1]*axe[0])
            return [
                sp*axe[0] + sn*axe[1],
                sp*axe[1] - sn*axe[0]
            ]
        return rotate(
            mirror(
                translate(self._XY(v), self._translation),
                self._mirror
            ), self._rotation
        )

    def draw_line(self, v1, v2, color=None, size=None):
        """
        Return the TIKZ code for a line according the drawing option given 
        to ``_drawing_tool``.

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     _drawing_tool, ParallelogramPolyominoesOptions
            ....: )
            sage: opt = ParallelogramPolyominoesOptions['tikz_options']
            sage: dt = _drawing_tool(opt)
            sage: dt.draw_line([1, 1], [-1, -1])
            '\n  \\draw[color=black, line width=1] (1.000000, 1.000000) -- (-1.000000, -1.000000);'

        """
        if color is None:
            color = self._color_line
        if size is None:
            size = self._line_size
        [x1, y1] = self.XY(v1)
        [x2, y2] = self.XY(v2)
        return "\n  \\draw[color=%s, line width=%s] (%f, %f) -- (%f, %f);" % (
            color, size, float(x1), float(y1), float(x2), float(y2)
        )

    def draw_point(self, p1, color=None, size=None):
        """
        Return the TIKZ code for a point according the drawing option given 
        to ``_drawing_tool``.

        EXAMPLES::

            sage: from sage.combinat.parallelogram_polyomino import (
            ....:     _drawing_tool, ParallelogramPolyominoesOptions
            ....: )
            sage: opt = ParallelogramPolyominoesOptions['tikz_options']
            sage: dt = _drawing_tool(opt)
            sage: dt.draw_point([1, 1])
            '\n  \\filldraw[color=black] (1.000000, 1.000000) circle (3.5pt);'

        """
        if color is None:
            color = self._color_point
        if size is None:
            size = self._point_size
        [x1, y1] = self.XY(p1)
        return "\n  \\draw[color=black, fill=%s] (%f, %f) circle (%spt);" % (
            color, float(x1), float(y1), size
        )

class TreeLikeTableau( ClonableList ):
    r"""
    The class of Tree-Like tableaux.

    Ref: J.C. Aval, A. Boussicault, P. Nadeau, Tree-like tableau,
    arXiv:1109.0371v2
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, value, check=True):
        r"""
        """
        if isinstance(value, cls):
            return value

        # We must verify ``t`` is a list of iterables, and also
        # normalize it to be a list of tuples.
        try:
            value = [list(_) for _ in value]
        except TypeError:
            raise ValueError("A tableau must be a list of iterables.")

        return TreeLikeTableaux_all().element_class( TreeLikeTableaux_all(), value, check )

    def __init__(self, parent, value, check=True):
        r"""
        The class for a tree-like tableau element.

        EXAMPLES::

            sage: tlt = TreeLikeTableau([[1,1,0,1],[1,0,1,0],[0,1],[0,1],[1,0]])
            sage: tlt
            [[1, 1, 0, 1], [1, 0, 1, 0], [0, 1], [0, 1], [1, 0]]
            sage: tlt = TreeLikeTableau( [[]] )
            sage: tlt
            [[]]
        """
        if isinstance(value, TreeLikeTableau):
            # Since we are (supposed to be) immutable, we can share the underlying data
            ClonableList.__init__(self, parent, value, check)
            return

        # Normalize t to be a list of lists.
        value = [list(_) for _ in value]

        ClonableList.__init__(self, parent, value, check)
        # This dispatches the input verification to the :meth:`check`
        # method.
        self._options = None

    def __eq__(self, other):
        r"""
        Check whether ``self`` is equal to ``other``.

        For now, two elements are equal if their underlying
        defining lists compare equal.

        INPUT:

        ``other`` -- the element that ``self`` is compared to

        OUTPUT:

        A Boolean.

        """
        if isinstance(other, TreeLikeTableau):
            return list(self) == list(other)
        else:
            return list(self) == other

    def __copy__(self):
        r"""
        Return a copy of the object
        """
        res = TreeLikeTableau( deepcopy( list(self) ) )
        res._set_mutable()
        return res

    def __hash__(self):
        return hash(tuple(map(tuple,list(self))))

    def _check_it_is_a_list_of_lists_of_0s_1s( self ):
        r"""
        Check if the the elements of the list that was given, contains lists
        with only only O's and 1's inside them.
        TESTS::
            sage: tlt = TreeLikeTableau( [[1,0,1,1],[8,1],[0,1]] )
            Traceback (most recent call last):
            ...
            ValueError: The element at the intersection of the column 0 and the row 1, i.e 8, is not a 0 or a 1.
            sage: tlt = TreeLikeTableau( [[1,0,1,1],[1,1],[0,1]] )
            sage: tlt
            [[1, 0, 1, 1], [1, 1], [0, 1]]
        """

        for irow in range(len(self)):
            element = ClonableList.__getitem__( self, irow)
            for icol in range(len(element)):
                if element[icol] != 0 and element[icol] !=1:
                    raise ValueError,\
                    "The element at the intersection of the column %s "%(icol)+\
                    "and the row %s, i.e %s, is not a 0 or a 1."%(irow,\
                    element[icol])

    def _check_there_is_a_root( self ):
        r"""
        Check if there is a 1 in the cell self[0][0].
        TESTS::
            sage: tlt = TreeLikeTableau( [[0,0,1,1],[1,1],[0,1]] )
            Traceback (most recent call last):
            ...
            ValueError: There is no root.
            sage: tlt = TreeLikeTableau( [[1,0,1,1],[1,1],[0,1]] )
            sage: tlt
            [[1, 0, 1, 1], [1, 1], [0, 1]]
        """
        if not self[0][0]==1:
            raise ValueError, "There is no root."

    def _check_column_have_a_point( self, icol ):
        r"""
        Check if each column have a 1.
        TESTS::
            sage: tlt = TreeLikeTableau( [[1,0,1,1],[1,0],[1,0]] )
            Traceback (most recent call last):
            ...
            ValueError: The column 1 is empty.
            sage: tlt = TreeLikeTableau( [[1,0,1,1],[1,1],[0,1]] )
            sage: tlt
            [[1, 0, 1, 1], [1, 1], [0, 1]]
        """
        irow = 0
        while( self[irow][icol] != -1 ):
            if self[irow][icol] == 1:
                return
            irow += 1
        raise ValueError, "The column %s is empty."%(icol)

    def _check_row_have_a_point( self, irow ):
        r"""
        Check if each row have a 1.
        TESTS::
            sage: tlt = TreeLikeTableau( [[1,0,1,1],[1,0],[1,0]] )
            Traceback (most recent call last):
            ...
            ValueError: The column 1 is empty.
            sage: tlt = TreeLikeTableau( [[1,0,1,1],[1,1],[0,1]] )
            sage: tlt
            [[1, 0, 1, 1], [1, 1], [0, 1]]
        """
        if not ( 1 in ClonableList.__getitem__( self, irow )):
            raise ValueError, "The row %s is empty."%(irow)

    def _check_cell_is_not_a_forbidden_pattern( self, irow, icol ):
        r"""
        Check if the cell self[irow][icol] is not a forbideen pattern.
        A forbidden pattern is a pointed cell fulfilling one of the
        two following conditions :
            - there is at least a point above him in the same column and
            at least a point to its left in the same row.
            - there is no point above him in the same column and
            no point to its left in the same row.

        TESTS::
            sage: tlt=TreeLikeTableau([[1,1,0,1],[1,0,1],[0,1,1]])
            Traceback (most recent call last):
            ...
            ValueError: The cell at the intersection of the row 2 and the column 2 has a point above him and at his left.
        """
        if irow == 0 and icol == 0:
            return
        point_above = False
        point_to_the_left = False
        c=0
        while(point_to_the_left == False and c<icol):
            point_to_the_left = (self[irow][c] == 1)
            c += 1
        r=0
        while(point_above == False and r<irow):
            point_above = (self[r][icol] == 1)
            r += 1
        if not point_above and not point_to_the_left:
            raise ValueError, "The cell at the intersection of the row %s and the column %s has no point above him and at his left."%(irow,icol)
        if point_above and point_to_the_left:
            raise ValueError, "The cell at the intersection of the row %s and the column %s has a point above him and at his left."%(irow,icol)

    def check( self ):
        r"""
        Check if the internal data contain valid information for a tree-like
        tableaux.

        EXAMPLES::

        TESTS::
        """
        if self._get_list() == [[]]:
            return
        self._check_it_is_a_list_of_lists_of_0s_1s()
        self._check_there_is_a_root()
        for irow in range( self.height() ):
            self._check_row_have_a_point( irow )
        for icol in range( self.width() ):
            self._check_column_have_a_point( icol )
        for irow in range( self.height() ):
            icol = 0
            while self[irow][icol] !=-1:
                if self[irow][icol] == 1:
                    self._check_cell_is_not_a_forbidden_pattern(irow, icol)
                icol += 1

    def options(self,  *get_value, **set_values):
        r"""
        Return all the options of the object.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.get_options()
            options for Parallelogram Polyominoes
        """
        if self._options is None:
            return self.parent()._options( *get_value, **set_values)
        return self._options( *get_value, **set_values)

    def _get_options(self):
        if self._options is None:
            return self.parent()._get_options()

    def __repr__(self):
        r"""
        TODO
        Return a string representation of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: pp
            [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            sage: pp.set_options(display='drawing')
            sage: pp
            [1 1 0]
            [1 1 0]
            [0 1 1]
        """
        return self._get_options().dispatch(self, '_repr_', 'display')

    def _repr_list(self):
        r"""
        TODO
        Return a string representation with list style.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: pp._repr_list()
            '[[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]'
        """
        return ClonableList._repr_(self)

    def _repr_drawing(self):
        res = ""
        for row in self:
            for cell in row:
                res+=str(cell)+" "
            res=res[-1]
        return res[:-1]


    def height( self ):
        r"""
        Give the number of rows of the tree-like tableau.

        EXAMPLES::

            sage: tlt = TreeLikeTableau( [[1,0,1,1],[1,1],[0,1]] )
            sage: tlt.height()
            3
            sage: tlt = TreeLikeTableau( [[1]] )
            sage: tlt.height()
            1
            sage: tlt = TreeLikeTableau( [[]] )
            sage: tlt.height()
            1
        """
        return ClonableList.__len__( self )

    def width( self ):
        r"""
        Give the number of columns of the tree-like tableau.

        EXAMPLES::
            sage: tlt = TreeLikeTableau( [[1,0,1,1],[1,1],[0,1]] )
            sage: tlt.width()
            4
            sage: tlt = TreeLikeTableau( [[1]] )
            sage: tlt.width()
            1
            sage: tlt = TreeLikeTableau( [[]] )
            sage: tlt.width()
            0
        """
        return len( ClonableList.__getitem__( self, 0 ) )

    def __getitem__(self, irow):
        r"""
        Get the entry of a tree-like-tableau.

        EXAMPLES::
            sage: tlt = TreeLikeTableau( [[1,1,0,0],[1,0,1,1],[1,0,0]] )
            sage: tlt[0][0]
            1
            sage: tlt[0][2]
            0
            sage: tlt[2][3]
            -1
            sage: tlt[-1][0]
            -1
            sage: tlt[0]
            [1, 1, 0, 0]
            sage: tlt[-1]
            Row outside the tree-like tableau.
            sage: tlt[3]
            Row outside the tree-like tableau.

        TESTS::
            sage: tlt[-1][-1]
            -1
            sage: tlt[0][-1]
            -1
            sage: tlt[2][4]
            -1
            sage: tlt[2][2]
            0
            sage: tlt[ [2, 2] ]
            0
        """
        if( isinstance( irow, list )):
            return self._row( self, irow[0] )[ irow[1] ]
        return self._row( self, irow )

    def __setitem__(self , irow, lst):
        r"""
        Modify an entry of the tree-like-tableau

        EXAMPLES::

            sage: tlt = TreeLikeTableau( [[1,0,1,0],[1,1,0,1],[1,0]] )
            sage: with tlt.clone() as tlt1:
            ....:     tlt1[2][0] = 0
            ....:     tlt1[2][1] = 1
            sage: tlt1
            [[1, 0, 1, 0], [1, 1, 0, 1], [0, 1]]
            sage: l =  [0, 1, 0]
            sage: with tlt.clone() as tlt2:
            ....:     tlt2[2] = l
            sage: tlt2
            [[1, 0, 1, 0], [1, 1, 0, 1], [0, 1, 0]]
            sage: tlt2[2] == [0, 1, 0]
            True
            sage: tlt[2] is l
            False
        """
        self._require_mutable()
        ClonableList.__setitem__(self , irow, list(lst) )

    class _row:
        r"""
        Class used for the functioning of __getitem__ and __setitem__
        """

        def __init__( self, tlt, irow ):
            self._tlt = tlt
            self._irow = irow
            if( irow < 0 or irow >= self._tlt.height() ):
                self._lst = None
            else:
                self._lst = ClonableList.__getitem__( self._tlt, irow )

        def __getitem__( self, icol ):
            if not self._lst:
                return -1
            elif icol < 0 or icol >= len( self._lst ):
                return -1
            else:
                return self._lst[icol]

        def __setitem__( self, icol , value):
            self._tlt._require_mutable()
            if self[icol]==-1:
                raise ValueError, "The cell at the intersection of the row %s and the column %s is outside the tree-like tableau"%(irow,icol)
            else:
                self._lst[icol] = value

        def __repr__(self):
            if self._lst:
                return str(self._lst)
            else:
                return "Row outside the tree-like tableau."

        def __len__(self):
            r"""
            Return the length of the row.
            """
            list=self._lst
            if list==None:
                return 0
            return len(list)

        def __eq__(self,other):
            r"""
            When == is used on tlt[i], it compares the list describing the row.

            TESTS::

                sage: tlt = TreeLikeTableau( [[1,0,1,0],[1,1,0,1],[1,0]] )
                sage: tlt[0]==[1,0,1,0]
                True
                sage: tlt = TreeLikeTableau( [[1,0,1,0],[1,1,0,1],[1,0]] )
                sage: tlt[1]==[1,0,1,0]
                False
                sage: tlt = TreeLikeTableau( [[]] )
                sage: tlt[0]==[]
                True
            """
            return self._lst==other


    def border_cells(self):
        r"""
        Give the list of the border cells starting from south-ouest to north-east.

        EXAMPLES::
            sage: tlt=TreeLikeTableau([[1,0,0,1],[1,0,1],[1,1]])
            sage: tlt
            [[1, 0, 0, 1], [1, 0, 1], [1, 1]]
            sage: tlt.border_cells()
            [(2, 0), (2, 1), (1, 1), (1, 2), (0, 2), (0, 3)]
            sage: tlt=TreeLikeTableau([[1,1,1,1]])
            sage: tlt
            [[1, 1, 1, 1]]
            sage: tlt.border_cells()
            [(0, 0), (0, 1), (0, 2), (0, 3)]
            sage: tlt=TreeLikeTableau([[]])
            sage: tlt.border_cells()
            []
        """
        lst = []
        if self == TreeLikeTableau([[]]):
            return lst
        irow = self.height()-1
        icol = 0

        while irow != -1:
            if lst != []:
                lst.append( (irow, icol) )
                icol += 1
            while self[irow][icol] != -1:
                lst.append( (irow ,icol) )
                icol += 1
            icol -= 1
            irow -=1
            while self[irow][icol+1] == -1 and irow != -1:
                lst.append( (irow, icol) )
                irow -= 1
        return lst


    def border_path(self):
        r"""
        Give the border path code by a list of 0's and 1's. The integer 0
        represents a east step whereas 1 represents a north step.

        EXAMPLES::

            sage: tlt = TreeLikeTableau([[1,1,0,1,1,1],[0,0,0,0,0,1],[1,0,1,0]
            ....: ,[1,0,0],[0,0,1]])
            sage: tlt.height()
            5
            sage: tlt.width()
            6
            sage: tlt.border_path()
            [0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1]

            sage: tlt = TreeLikeTableau([[]])
            sage: tlt.border_path()
            [1]
        """

        if self == TreeLikeTableau([[]]):
            return [1]
        part = self.partition()
        lst = []
        for irow in range(self.height()-1):
            lst.insert(0,1)
            for icol in range(part[irow]-part[irow+1]):
                lst.insert(0,0)
        lst.insert(0,1)
        for icol in range(part[self.height()-1]):
            lst.insert(0,0)
        return lst

    def _insert_ribbon(self, pos1, pos2):
        r"""
        Insert a ribbon of empty cells from the position 'pos1' to the position
        'pos2'.

        TESTS:: TODO clone check=false or enlever la mutabilite il y a une fonction qui le permet

            sage: tlt = TreeLikeTableau( [[1,1,1,0],[0,1,0,1],[1,0],[1]] )
            sage: tlt._set_mutable()
            sage: tlt._insert_ribbon(1, 5)
            sage: tlt
            [[1, 1, 1, 0], [0, 1, 0, 1], [1, 0, 0, 0], [1, 0, 0]]


        """
        self._require_mutable()
        [irow, icol, direction] = self._border_position_to_coordinates_and_direction(pos1)

        k = pos1
        border_path = self.border_path()
        while( k <= pos2):
            while( k <= pos2 and border_path[k] == 1):
                self[irow]._lst.append(0)
                irow -= 1
                k += 1
            if k <= pos2:
                irow += 1
                self[irow]._lst.pop(-1)
            while( k <= pos2 and border_path[k] == 0 ):
                self[irow]._lst.append(0)
                k += 1
            if k <= pos2:
                self[irow]._lst.append(0)
                irow -= 1

    def _insert_column(self, irow, icol):
        r"""
        Insert a column of 'irow+1' cells, with only the lowest one dotted, between
        the column 'icol-1' and 'icol'.

        TESTS::

            sage: tlt = TreeLikeTableau( [[1,1,0,1],[1,0,1,0],[0,1],[0,1],[1,0]] )
            sage: tlt._set_mutable()
            sage: tlt._insert_column(4,2)
            sage: tlt
            [[1, 1, 0, 0, 1], [1, 0, 0, 1, 0], [0, 1, 0], [0, 1, 0], [1, 0, 1]]

            sage: tlt = TreeLikeTableau( [[]] )
            sage: tlt._set_mutable()
            sage: tlt._insert_column(0,0)
            sage: tlt
            [[1]]
        """
        self._require_mutable()
        for k in range(irow):
            self[k]._lst.insert(icol,0)
        self[irow]._lst.insert(icol,1)

    def _insert_row(self, irow, icol):
        r"""
        Insert a row of 'icol+1' cells, with only the right most cell dotted,
        between the rows 'irow-1' and 'irow'.

        TESTS::

            sage: tlt = TreeLikeTableau( [[1,0,1],[1,0],[1,1]] )
            sage: tlt._set_mutable()
            sage: tlt._insert_row(1,2)
            sage: tlt
            [[1, 0, 1], [0, 0, 1], [1, 0], [1, 1]]
        """
        self._require_mutable()
        lst=[]
        for k in range(icol):
            lst.append(0)
        lst.append(1)
        return self._get_list().insert(irow,lst)

    def _border_position_to_coordinates_and_direction(self, position):
        r"""
        Give the coodinates of the cell outside the tree-like tableau
        having the border labelled by 'position' as an edge, and
        the direction of this border, i.e 'horizontal' or 'vertical'.
        """
        irow = self.height()
        icol = -1
        k = 0
        border_path = self.border_path()
        while( k <= position):
            while( k <= position and border_path[k] == 0 ):
                icol += 1
                k += 1
            if k <= position:
                icol += 1
            while( k <= position and border_path[k] == 1):
                irow -= 1
                k += 1
            if k <= position:
                icol -= 1
        if border_path[position] == 0:
            direction = 'horizontal'
        else:
            direction = 'vertical'
        return [irow, icol, direction]

    def special_point(self):
        r"""
        Give the position of the special point in the border path.

        EXAMPLES::

            sage:
        """
        border_cells = self.border_cells()
        for k in range(len(border_cells)-1,0,-1):
            (irow, icol) = border_cells[k]
            if self[irow][icol] == 1 and border_cells[k-1][0] == irow:
                return k
        return 0
    def insert_point(self, position ):
        r"""
        This is the insertion algortihm.
        This method insert a new point inside the tree-like tableau on the
        border number 'position' (the border are numbered from left to right)

        ( see definition 2.2 of the article "tree-like tableaux", J.C. Aval,
          A. Boussicault, P. Nadeau, arXiv:1109.0371v2 )

        EXAMPLES::
            sage: tlt = TreeLikeTableau( [[1,1,0,1],[1,0,1,0],[0,1],[0,1],[1,0]] )
            sage: with tlt.clone() as tlt1:
            ....:     tlt1.insert_point( 3 )
            sage: tlt1
            [[1, 1, 0, 0, 1], [1, 0, 0, 1, 0], [0, 1, 0, 0], [0, 1, 1, 0], [1, 0]]
        """
        self._require_mutable()
        [irow, icol, direction] = self._border_position_to_coordinates_and_direction( position)
        position_special_point=self.special_point()
        if direction == 'horizontal':
            self._insert_row(irow, icol)
        else:
            self._insert_column(irow, icol)
        if position < position_special_point:
            self._insert_ribbon(position + 1, position_special_point + 1)

    def _remove_ribbon(self):
        r"""
        This method correspond to the first step of remove_point.

        TESTS::
            sage: tlt = TreeLikeTableau([[1,0,1,1,1],[0,0,0,0,1],[1,1,0,0,0],[1,0,0,0,0],[0,1,0],[1]])
            sage: tlt._set_mutable()
            sage: tlt._remove_ribbon()
            sage: tlt
            [[1, 0, 1, 1, 1], [0, 0, 0, 0, 1], [1, 1, 0, 0], [1, 0], [0, 1], [1]]

            sage: tlt = TreeLikeTableau([[1,1,1,0,1],[0,0,1,0,0],[0,1,0,1],[1,0,0],[1]])
            sage: tlt._set_mutable()
            sage: tlt._remove_ribbon()
            sage: tlt
            [[1, 1, 1, 0, 1], [0, 0, 1, 0, 0], [0, 1, 0, 1], [1, 0, 0], [1]]
        """
        if self == TreeLikeTableau([[1]]):
            return
        border_cells = self.border_cells()
        k = self.special_point()
        if k+1 == len(border_cells):
            return
        if border_cells[k][0] != border_cells[k+1][0]:
            return
        k=k+1
        (irow, icol) = border_cells[k]
        old_self = self.__copy__()
        while old_self[irow][icol] == 0:
            self._get_list()[irow].pop()
            k = k+1
            (irow, icol) = border_cells[k]

    def _remove_column(self, icol):
        r"""
        Remove the column number ''icol''.

        TESTS::

            sage: tlt = TreeLikeTableau([[1,1,1,0,1],[0,0,1,0,0],[0,1,0,1],[1,0,0],[1]])
            sage: tlt._set_mutable()
            sage: tlt._remove_column(3)
            sage: tlt
            [[1, 1, 1, 1], [0, 0, 1, 0], [0, 1, 0], [1, 0, 0], [1]]
        """
        irow = 0
        while self[irow][icol] != -1:
            self[irow]._lst.pop(icol)
            irow +=1

    def _remove_row(self, irow):
        r"""
        Remove the row number ''irow''.

        TESTS::

            sage: tlt = TreeLikeTableau([[1,1,1,1],[0,0,0,1],[0,1],[1,0]])
            sage: tlt._set_mutable()
            sage: tlt._remove_row(1)
            sage: tlt
            [[1, 1, 1, 1], [0, 1], [1, 0]]
        """
        self.pop(irow)

    def remove_point(self):
        r"""
        This is the inverse of the insertion algortihm.
        This method remove the special point from the tree-like tableau.
        This method return the position border of the removed point.

        ( see article "tree-like tableaux", J.C. Aval, A. Boussicault,
          P. Nadeau, arXiv:1109.0371v2 )

        EXAMPLES::
            sage: tlt = TreeLikeTableau( [[1,1,0,1,1,1],[0,1,0,0,0,0],[1,0,0,0],[0,1,1,0],[1,0]] )
            sage: with tlt.clone() as tlt1:
            ....:     tlt1.remove_point()
            3
            sage: tlt1
            [[1, 1, 1, 1, 1], [0, 1], [1, 0], [0, 1], [1, 0]]
        """
        self._require_mutable()
        border_sp = self.special_point()
        self._remove_ribbon()
        [irow_sp, icol_sp, direction] = self._border_position_to_coordinates_and_direction(border_sp)
        for k in range(irow_sp-1):
            if self[k][icol_sp] == 1:
                self._remove_row(irow_sp-1) #irow_sp correspond to the line under ''border_sp'' horizontal border.
                return border_sp
        self._remove_column(icol_sp)
        return border_sp

    def partition(self):
        r"""
        Give the partition given by the rows.

        TESTS::
            sage: tlt = TreeLikeTableau( [[1,0,1],[1,1],[1,0],[1]] )
            sage: tlt.partition()
            [3, 2, 2, 1]
            sage: tlt = TreeLikeTableau( [[]] )
            sage: tlt.partition()
            []
        """
        return Partition([len(self[irow]) for irow in range(self.height())])

    def transpose(self):
        r"""
        Transpose the tree-like tableau 'self'.

        'self' needs to be mutable.

        EXAMPLES::

            sage: tlt = TreeLikeTableau( [[1,1,0,1],[1,0,1],[0,1,0],[0,1],[1]])
            sage: tlt
            [[1, 1, 0, 1], [1, 0, 1], [0, 1, 0], [0, 1], [1]]
            sage: with tlt.clone() as tlt1:
            ....:   tlt1.transpose()
            ....:   tlt1
            [[1, 1, 0, 0, 1], [1, 0, 1, 1], [0, 1, 0], [1]]
            sage: tlt.partition() == tlt1.partition().conjugate()
            True
            sage: tlt = TreeLikeTableau( [[]] )
            sage: tlt
            [[]]
            sage: with tlt.clone() as tlt1:
            ....:   tlt1.transpose()
            ....:   tlt1
            [[]]
        """
        self._require_mutable()
        if self==TreeLikeTableau([[]]):
            return
        lst = []
        for icol in range(self.width()):
            lst.append([])
            for irow in range(self.height()):
                tmp = self[irow][icol]
                if tmp == -1:
                    break
                lst[icol].append(self[irow][icol])
        self._set_list(lst)

    def _repr_(self):
        r"""
        The text representation of a tree-like tableau

        EXAMPLES::
            sage: tlt = TreeLikeTableau( [[1,0,1,0],[1,1,0,1],[0,1],[0,1],[1,0]] )
            sage: tlt
            [[1, 0, 1, 0], [1, 1, 0, 1], [0, 1], [0, 1], [1, 0]]
        """
        return self._repr_list()
        return self.parent().global_options.dispatch(self, '_repr_', 'display')

    def _repr_list( self ):
        return ClonableList._repr_(self)

    def _repr_drawing( self ):
        res=""
        for line in self._get_list():
            res+=str(line)+'\n'
        return res[:-1] #remove the last '\n'

    def get_tikz_options( self ):
        res = default_tikz_options
        user = self.options()['tikz']
        for opt in user:
            res[opt] = user[opt]
        return res

    def to_tikz( self ):
        tikz_options = get_tikz_options()
        drawing_components = self.options()['drawing_components']
        NotImplemented

    def size( self ):
        r"""
        Return the size of the tree-like tableaux : it is the number of point
        (the 1's inside the tableau) inside the tableau.

        EXAMPLES::
            sage: tlt = TreeLikeTableau( [[1,0,1,0],[1,1,0,1],[0,1],[0,1],[1,0]] )
            sage: tlt.size()
            8
            sage: tlt = TreeLikeTableau( [[1]] )
            sage: tlt.size()
            1
            sage: tlt = TreeLikeTableau( [[]] )
            sage: tlt.size()
            0
        """
        return self.height().__add__(self.width()) -1

    def _latex_(self):
        r"""
        Return a LaTeX version of ``self``.

        For more on the latex options, see :meth:`Partitions.global_options`.

        """
        return self.parent().global_options.dispatch(self, '_latex_', 'latex')

    def _latex_drawing( self ):
        latex.add_package_to_preamble_if_available("tikz")
        tikz_options = self.get_tikz_options()
        res = "\n\\begin{tikzpicture}[scale=%s]"%(tikz_options['scale'])
        tikz_code = self.to_tikz()
        res += "\n\\end{tikzpicture}"
        return res

    def _latex_list( self ):
        NotImplemented

    def top_points( self ):
        r"""
        Top points are the points of the first row except the root point.

        EXAMPLES::

            sage: tlt = TreeLikeTableau( [[1,0,1,0],[1,1,0,1],[0,1],[0,1],[1,0]] )
            sage: tlt.top_points()
            [(0, 2)]
            sage: tlt = TreeLikeTableau( [[1,1,0,1],[1,0,1],[0,1,0],[0,1],[1]])
            sage: tlt.top_points()
            [(0, 1), (0, 3)]
            sage: tlt = TreeLikeTableau( [[1,0,0],[1,0,1],[1,1]] )
            sage: tlt.top_points()
            []

        TESTS::

            sage: tlt = TreeLikeTableau([[]])
            sage: tlt.top_points()
            []
        """
        L=[]
        row0=self[0]
        for i in range(1,len(row0)):
            if row0[i]==1:
                L.append((0,i))
        return L

    def left_points( self ):
        r"""
        Left points are the points of the first column except the root point.

        EXAMPLES::

            sage: tlt = TreeLikeTableau( [[1,0,1,0],[1,1,0,1],[0,1],[0,1],[1,0]] )
            sage: tlt.left_points()
            [(1, 0), (4, 0)]
            sage: tlt = TreeLikeTableau( [[1,1,0,1],[1,0,1],[0,1,0],[0,1]])
            sage: tlt.left_points()
            [(1, 0)]
            sage: tlt = TreeLikeTableau( [[1,1,0],[0,1,1],[0,1]] )
            sage: tlt.left_points()
            []

        TESTS::

            sage: tlt = TreeLikeTableau([[]])
            sage: tlt.left_points()
            []
        """
        L=[]
        i=1
        while self[i][0]!=-1:
            if self[i][0]==1:
                L.append((i,0))
            i=i+1
        return L

    def crossings( self ):
        r"""
        A crossing is a cell which have a point above him in the same column
        and at its left in the same row.

        EXAMPLES::

            sage: tlt = TreeLikeTableau( [[1,0,1,0],[1,1,0,1],[0,1],[0,1],[1,0]] )
            sage: tlt.crossings()
            [(1, 2), (4, 1)]
            sage: tlt = TreeLikeTableau( [[1,1,1,1],[1,0,0,0],[1,0,0,0],[1,0,0,0]] )
            sage: tlt.crossings()
            [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]
            sage: tlt = TreeLikeTableau( [[1,0,1,1,1],[0,0,0,0,1],[0,0,1],[1,0],[1,1]])
            sage: tlt.crossings()
            []

        TESTS::

            sage: tlt = TreeLikeTableau([[]])
            sage: tlt.crossings()
            []


        """
        L=[]
        irow=0
        icol=0
        while self[irow][icol]!=-1:
            while self[irow][icol]!=-1:
                test_row=False
                test_col=False
                j=irow-1
                while not test_row and self[j][icol]!=-1:
                    test_row = self[j][icol]==1
                    j=j-1
                j=icol-1
                while not test_col and self[irow][j]!=-1:
                    test_col = self[irow][j]==1
                    j=j-1
                if test_row and test_col:
                    L.append((irow,icol))
                icol=icol+1
            icol=0
            irow=irow+1
        return L

    def nb_ribbons( self ):
        r"""
        Give the number of ribbons that are inserted when we construct ``self``
        from the tree-like tableau of size 1 by doing successive insertions.

        EXAMPLES::
            sage: tlt = TreeLikeTableau( [[1,1,1,1],[1,0,0,0],[1,0,0,0],[1,0,0,0]] )
            sage: tlt.nb_ribbons()
            3
            sage: tlt=TreeLikeTableau([[1,1,0,0,0,0],[0,1,1,1,1,0],[0,0,0,0,1,1],[0,0,0,0,0,1],[0,0,0,1]])
            sage: tlt.nb_ribbons()
            0

        TESTS::
            sage: tlt=TreeLikeTableau([[1]])
            sage: tlt.nb_ribbons()
            0
            sage: tlt=TreeLikeTableau([[]])
            sage: tlt.nb_ribbons()
            0
        """
        code=self.code()
        res=0
        for i in range(len(code)-1):
            if code[i+1]<code[i]:
                res+=1
        return res

    def corners( self ):
        r"""
        A corner is a cell for which the bottom and the right edges are border edges.
        A corner correspond de the pattern 01 in the border path.

        EXAMPLES::

            sage: tlt = TreeLikeTableau( [[1,1,0,1],[1,0,1],[0,1,0],[0,1],[1]])
            sage: tlt.corners()
            [(0, 3), (2, 2), (3, 1), (4, 0)]
            sage: tlt = TreeLikeTableau([[1]])
            sage: tlt.corners()
            [(0, 0)]
            sage: tlt=TreeLikeTableau([[]])
            sage: tlt.corners()
            []
        """
        return self.partition().corners()

    def occupied_corners( self ):
        r"""
        An occupied corner is a corner which has a point.

        EXAMPLES::

            sage: tlt = TreeLikeTableau( [[1,1,0,1],[1,0,1],[0,1,0],[0,1],[1]])
            sage: tlt.occupied_corners()
            [(0, 3), (3, 1), (4, 0)]
            sage: tlt = TreeLikeTableau([[1]])
            sage: tlt.occupied_corners()
            [(0, 0)]
            sage: tlt=TreeLikeTableau([[]])
            sage: tlt.occupied_corners()
            []
        """
        corners= self.corners()
        j=0
        for i in range(len(corners)):
            (irow,icol)=corners[i-j]
            if self[irow][icol]==0:
                corners.pop(i-j)
                j=j+1
        return  corners

    def empty_corners( self ):
        r"""
        An empty corner is a corner which doesn't have a point.

        EXAMPLES::

            sage: tlt = TreeLikeTableau( [[1,1,0,1],[1,0,1],[0,1,0],[0,1],[1]])
            sage: tlt.empty_corners()
            [(2, 2)]
            sage: tlt = TreeLikeTableau([[1]])
            sage: tlt.empty_corners()
            []
            sage: tlt=TreeLikeTableau([[]])
            sage: tlt.empty_corners()
            []
        """
        corners= self.corners()
        j=0
        for i in range(len(corners)):
            (irow,icol)=corners[i-j]
            if self[irow][icol]==1:
                corners.pop(i-j)
                j=j+1
        return  corners

    def is_symmetric( self ):
        r"""
        A symetric tree-like tableau is a tree-like tableau invariant by
        the axial symmetry with respect to its main diagonal.

        EXAMPLES::

            sage: tlt = TreeLikeTableau([[1,0,1,1],[0,0,0,1],[1,0,0],[1,1]])
            sage: tlt.is_symmetric()
            True
            sage: tlt = TreeLikeTableau([[1,1],[1]])
            sage: tlt.is_symmetric()
            True
            sage: tlt = TreeLikeTableau([[1,1,1],[1,0,0],[0,1]])
            sage: tlt.is_symmetric()
            False
            sage: tlt = TreeLikeTableau([[1]])
            sage: tlt.is_symmetric()
            True
            sage: tlt = TreeLikeTableau([[]])
            sage: tlt.is_symmetric()
            True
        """
        if self == TreeLikeTableau([[]]):
            return True
        i=0
        while i+1<=len(self[i]):
            for j in range(i+1,len(self[i])+1):
                if self[i][j]!=self[j][i]:
                    return False
            i=i+1
        return True

    def code( self ):
        r"""
        The code of a tree-like tableau T is the sequence of the indexes of
        the border edges at which the consecutive insertions were done in order to
        obtain T starting with the tree-like tableau of size 1.

        EXAMPLES::
            sage:tlt=TreeLikeTableau([[]]])
            sage:with tlt.clone() as tlt1:
            ....:   tlt1.insert_point(0)
            ....:   tlt1.insert_point(1)
            ....:   tlt1.insert_point(0)
            ....:   tlt1.insert_point(3)
            ....:   tlt1.insert_point(2)
            ....:   tlt1.insert_point(4)
            ....:   tlt1.insert_point(3)
            sage:tlt1.code()
            [0, 1, 0, 3, 2, 4, 3]

        TESTS::
            sage:tlt=TreeLikeTableau([[1]])
            sage:tlt.code()
            [0]
            sage:tlt=TreeLikeTableau([[]])
            sage:tlt.code()
            []

        """
        L=[]
        tlt=self.__copy__()
        for i in range(tlt.size()):
            L.insert(0,tlt.remove_point())
        return L

    def to_permutations_by_CN( self ):
        r"""
        Give the permutation associated to self using the bijection
        defined by S. Corteel and P. Nadeau Sylvie in :
            Bijections for permutation tableaux.
            European J. Combin., 30(1):295-310. 2009.

        EXAMPLES::

            sage: tlt = TreeLikeTableau([
            ....: [1,0,0,0,1],
            ....: [0,0,0,0,1],
            ....: [1,0,1,0],
            ....: [1,1,0,1],
            ....: [0,1,0]])
            sage: tlt.to_permutations_by_CN()
            [2, 3, 1, 8, 4, 7, 9, 6, 5]

            sage: tlt = TreeLikeTableau([
            ....: [1,1,0,1,0,0,0,0,0,1],
            ....: [0,0,0,1,0,0,0,0,0],
            ....: [1,0,1,0,0,1,1,0,0],
            ....: [0,0,0,0,0,1,0,1,0],
            ....: [0,0,1,0,0,0,0,0,0],
            ....: [0,1,0,0,0,0,0,0,1],
            ....: [0,0,0,0,0,1,0,0],
            ....: [1,0,0,0,1,0,0,0],
            ....: [0,0,1,0,0,0,0,0],
            ....: [0,0,0,0,0,0,0,1],
            ....: [0,1]])
            sage: tlt.to_permutations_by_CN()
            [8, 7, 19, 20, 3, 17, 2, 1, 6, 11, 18, 12, 13, 5, 9, 15, 14, 4, 16, 10]

            sage: for tlt in TreeLikeTableaux(6):
            ....:     p = tlt.to_permutations_by_CN()
            ....:     tlt1 = TreeLikeTableau.from_permutation_by_CN(p)
            ....:     if tlt != tlt1:
            ....:         print tlt
            ....:         print tlt1

        """
        if self == TreeLikeTableau([[]]):
            return Permutation(None)
        row_numbering=[]
        column_numbering=[]
        path = self.border_path()
        for i in range(1, self.size()+1):
            if path[-i] == 1:
                row_numbering.append(i)
            else:
                column_numbering.append(i)
        P=[]
        for i in range(self.height()): #initialisation with left points
            if self[i][0]==1:
                P.append(row_numbering[i])
        for i in range(1,self.width()):
            j=0
            while self[j][i] == 0: #looking for the highest point in column i
                j+=1
            index = P.index(row_numbering[j])
            j+=1
            while self[j][i] != -1:
                if self[j][i] == 1:
                    P.insert(index,row_numbering[j])
                    index+=1
                j+=1
            P.insert(index,column_numbering[-i])
        return Permutation(P)

    def to_permutations_by_the_code( self ):
        r"""
        Give the permutation associated to self by the bijection using the code.
        The code of self is the non invertion table of the returned permutation.
        ( see Theorem 4.2 of the article "tree-like tableaux", J.C. Aval,
          A. Boussicault, P. Nadeau, arXiv:1109.0371v2 )

        EXAMPLES::
            sage: tlt=TreeLikeTableau([[1,0,0,1],[1,1,0,0],[0,1,1],[1,0]])
            sage: tlt.code()
            [0, 0, 2, 1, 1, 0, 3]
            sage: tlt.to_permutations_by_the_code()
            [6, 2, 7, 5, 3, 1, 4]
            sage: tlt=TreeLikeTableau([[1,1,1,1,1,1]])
            sage: tlt.code()
            [0, 1, 2, 3, 4, 5]
            sage: tlt.to_permutations_by_the_code()
            [1, 2, 3, 4, 5, 6]
            sage: tlt=TreeLikeTableau([[1]])
            sage: tlt.code()
            [0]
            sage: tlt.to_permutations_by_the_code()
            [1]
            sage: tlt=TreeLikeTableau([[]])
            sage: tlt.code()
            []
            sage: tlt.to_permutations_by_the_code()
            []

        """
        L=range(1,self.size()+1)
        code=self.code()
        P=[]
        for i in code[::-1]:
            P.insert(0,L.pop(i))
        return Permutation(P)

    def to_increasing_trees( self ):
        r"""
        Return the increasing tree associated with self.
        A tree-like tableau is uniquely obtained by performing successive insertion
        starting from the tree-like tableau of size 0. If we label a point
        inserted during the i-th insertion, by i, then the labelled tree that appear
        inside the tree-like tableau, is the increasing tree returned by this function.
        ( see Theorem 4.3 of the article "tree-like tableaux", J.C. Aval,
          A. Boussicault, P. Nadeau, arXiv:1109.0371v2 )

        EXAMPLES::
            sage: tlt=TreeLikeTableau([[1,0,0,1],[1,1,0,0],[0,1,1],[1,0]])
            sage: tlt.to_increasing_trees()
            1[2[6[., .], 4[5[., 7[., .]], .]], 3[., .]]
            sage: tlt = TreeLikeTableau( [[1,1,0,1],[1,0,1],[0,1,0],[0,1],[1]] )
            sage: tlt.to_increasing_trees()
            1[3[4[., .], 5[., .]], 2[6[7[., .], .], 8[., .]]]
            sage: tlt=TreeLikeTableau([[1]])
            sage: tlt.to_increasing_trees()
            1[., .]
            sage: tlt=TreeLikeTableau([[]])
            sage: tlt.to_increasing_trees()
            .
        """
        if self == TreeLikeTableau([[]]):
            return LabelledBinaryTree(None)
        else:
            columns=range(self.width())
            rows=range(self.height())
            L=deepcopy( list(self) )
            tlt=self.clone()
            tlt._set_mutable()
            while tlt.size()>1:
                (sp_row,sp_col)=tlt.border_cells()[tlt.remove_point()]
                L[rows[sp_row]][columns[sp_col]]=tlt.size()+1
                if tlt.height()<len(rows):
                    rows.pop(sp_row)
                else:
                    columns.pop(sp_col)

            def rec(self,cell):
                lson=None
                rson=None
                (row,col)=cell
                irow=row+1
                icol=col+1
                while self[irow][col]!=-1:
                    if self[irow][col]==1:
                        lson = rec(self,(irow,col))
                        break
                    irow+=1
                while self[row][icol]!=-1:
                    if self[row][icol]==1:
                        rson = rec(self,(row,icol))
                        break
                    icol+=1
                return LabelledBinaryTree([lson,rson],label=L[row][col])
            return rec(self,(0,0))

    def to_binary_trees( self ):
        r"""
        Give the binary tree defined by the points inside self.
        It is a surjection from tree-like tableaux of size n
        to binary trees of size n.
        """
        if self == TreeLikeTableau([[]]):
            return BinaryTree()
        else:
            def rec(tlt,cell):
                lson=None
                rson=None
                (row,col)=cell
                irow=row+1
                icol=col+1
                while tlt[irow][col]!=-1:
                    if tlt[irow][col]==1:
                        lson = rec(tlt,(irow,col))
                        break
                    irow+=1
                while tlt[row][icol]!=-1:
                    if tlt[row][icol]==1:
                        rson = rec(tlt,(row,icol))
                        break
                    icol+=1
                return BinaryTree([lson,rson])
            return rec(self,(0,0))

    def to_permutations_by_increasing_tree( self ):
        r"""
        Give the permutation associated to self by doing an in order traversal
        of self.to_increasing_trees().
        """
        l=[]
        inct=self.to_increasing_trees()
        inct.in_order_traversal(lambda node: l.append( node.label() ))
        return Permutation(l)

#    def to_non_ambiguous_trees( self ):
#        r"""
#        Give the non-ambiguous tree obtained by completing the Ferrer diagram.
#        """
#        w = self.width()
#        l = deepcopy( self._get_list() )
#        for i in range(len(l)):
#            length = len(l[i])
#            l[i] = l[i] + [0]*(w-length)
#        return NonAmbiguousTree( l )

    @classmethod
    def from_permutation(cls, perm):
        r"""
        Give the tree-like tableau T associated to `perm` using the bijection
        given by `bijection` :
            code : see from_permutation_by_code
            inc_tree : see from_permutation_by_increasing_tree
            CN : see from_permutation_by_CN
        """
        bijection = 'code'
        if not (perm in Permutations()):
            raise ValueError, str(perm)+" is not a permutation"

        if bijection=='code':
            return TreeLikeTableau.from_permutation_by_code(perm)
        if bijection=='inc':
            return TreeLikeTableau.from_permutation_by_increasing_tree(perm)
        if bijection=='CN':
            return TreeLikeTableau.from_permutation_by_CN(perm)

    @classmethod
    def from_permutation_by_increasing_tree(cls, perm):
        r"""
        Give the tree-like tableau T such that perm=T.to_permutations_by_increasing_tree().
        """
        inct = perm.increasing_tree()
        return from_increasing_trees(inct)

    @classmethod
    def from_increasing_trees(cls, inct):
        r"""
        Give the tree-like tableau T such that `inct`=T.to_increasing_trees().
        """
        n = inct.node_number()
        L = [0]*(n+1) # The father of the node i will be in L[i]. L[0] and L[1] are not used.
        def tree_rec( t ):
            if not t[0].is_empty():
                L[t[0].label()] = [t.label(),0] #0 <-> left
                tree_rec(t[0])
            if not t[1].is_empty():
                L[t[1].label()] = [t.label(),1] #1 <-> right
                tree_rec(t[1])
        tree_rec(inct)

        #Construct the tree-like tableau with the inc tree inside
        tlt = TreeLikeTableau([[1]])
        tlt._set_mutable()
        for n in len(L)[2:]:
            father, son_type = L[n]
            #find father in tlt
            i=0
            while not father in tlt[i]:
                i+=1
            j=list(tlt[i]._lst).index(father)

            #find border cell
            border_cells = tlt.border_cells()
            if son_type==0:
                while tlt[i+1][j]!=-1:
                    i+=1
                insertion_index = border_cells.index(tlt[i][j])
                tlt.insert_point(insertion_index)
                tlt[i+1][j]=n
            if son_type==1:
                while tlt[i][j+1]!=-1:
                    j+=1
                insertion_index = border_cells.index(tlt[i][j])
                tlt.insert_point(insertion_index)
                tlt[i][j+1]=n

        #Replace all the positive integers by 1
        for row in range(tlt.height()):
            for col in range(len(tlt[row])):
                if tlt[row][col] != 0:
                    tlt[row][col]=1
        tlt.set_immutable()
        return tlt

    @classmethod
    def from_permutation_by_CN(cls, p):
        r"""
        Give the tree-like tableau T such that `p`=T.to_permutations_by_CN().
        """
        n = p.size()
        #we start by constructing the shape of the TLT using descents of `p`.
        descents=[ p(i) for i in range(1,n) if p(i)>p(i+1) ]
        descents.sort()
        descents.append(n+1)
        l = len(descents)
        i=1
        shape = []
        for j in descents:
            while i < j:
                shape.append( [0]*l )
                i+=1
            i+=1
            l=l-1
        descents.pop(-1)
        #we add points to the empty shape
        l = list(p)
        column_numbering = descents
        width = len(column_numbering)+1
        row_numbering = [i for i in range(1,n+1) if i not in column_numbering]
        for j in range(width-1):
            col = column_numbering[j]
            index_col = l.index(col)
            r = l[index_col+1]
            while index_col>0 and l[index_col-1]>r and not l[index_col-1] in column_numbering:
                row = l.pop(index_col-1)
                shape[row_numbering.index(row)][(width-1)-j] = 1
                index_col -= 1
            row = r
            shape[row_numbering.index(row)][width-1-j] = 1
            l.pop(index_col)
        for row in l:#the points of the first column are left
            shape[row_numbering.index(row)][0] = 1
        return TreeLikeTableau(shape)

    @classmethod
    def from_permutation_by_code(cls, perm):
        r"""
        Give the tree-like tableau T such that ``perm``=T.to_permutations_by_the_code().

        TESTS::
            sage:for tlt in TreeLikeTableaux(6):
            ....:   if not tlt==TreeLikeTableau.from_permutation_by_code(tlt.to_permutations_by_the_code()):
            ....:       print tlt
        """
        p=list(perm)
        tlt=TreeLikeTableau([[]])
        tlt._set_mutable()
        for i in range(len(perm)):
            cmp=0
            for j in range(i+1):
                pi=p[i]
                if p[j]<pi:
                    cmp+=1
                j+=1
            tlt.insert_point(cmp)
        tlt.set_immutable()
        return tlt

    @classmethod
    def from_non_ambiguous_trees( nat ):
        r"""
        Return the list of tree-like tableaux T such that ``nat``=T.to_non_ambiguous_trees().
        """
        if nat == NonAmbiguousTree([[]]):
            return [TreeLikeTableau([[]])]

        inner_partition = [0]
        array = nat.get_array()
        height = len(array)
        width = len(array[0])

        for i in range(1,height+1):
            row = array[-i]
            last = 0
            for i in range(len(row)):
                if row[i]==1:
                    last = i
            last += 1
            last = max(last,inner_partition[0])
            inner_partition.insert(0,last)
        inner_partition.pop()
        size_inner = sum(inner_partition)
        size_outer = width*height
        outer_partitions = []
        for size in range(size_inner,size_outer+1):
            outer_partitions = outer_partitions + Partitions(
                size, inner=inner_partition, outer=[width]*height
            ).list()

        tlt_list = []
        for p in outer_partitions:
            tlt = [array[i][:p[i]] for i in range(len(p))]
            tlt_list.append(TreeLikeTableau(tlt))
        return tlt_list

class TreeLikeTableaux(UniqueRepresentation, Parent):
    r"""
    The tree-like tableaux factory.
    """
    @staticmethod
    def __classcall_private__(self, size=None, symmetric=False):
        """
        TESTS::

            sage: from sage.combinat.tree_like_tableaux import (TreeLikeTableaux_all,
            ....:  TreeLikeTableaux_size, SymmetricTreeLikeTableaux_all,
            ....:  SymmetricTreeLikeTableaux_size)
            sage: isinstance(TreeLikeTableaux(2), TreeLikeTableaux)
            True
            sage: isinstance(TreeLikeTableaux(), TreeLikeTableaux)
            True
            sage: isinstance(TreeLikeTableaux(3, symmetric=True), TreeLikeTableaux)
            True
            sage: isinstance(TreeLikeTableaux(symmetric=True), TreeLikeTableaux)
            True
            sage: TreeLikeTableaux(2) is TreeLikeTableaux_size(2)
            True
            sage: TreeLikeTableaux(5).cardinality()
            120
            sage: TreeLikeTableaux() is TreeLikeTableaux_all()
            True
            sage: TreeLikeTableaux(3, symmetric=True) is SymmetricTreeLikeTableaux_size(3)
            True
            sage: TreeLikeTableaux(5, symmetric=True).cardinality()
            8
            sage: TreeLikeTableaux(symmetric=True) is SymmetricTreeLikeTableaux_all()
            True
        """
        if size is None:
            if symmetric:
                return SymmetricTreeLikeTableaux_all()
            else:
                return TreeLikeTableaux_all()
        else:
            if not (isinstance(size, (Integer, int)) and size >= 0):
                raise ValueError("size must be a nonnegative integer")
            if not symmetric:
                return TreeLikeTableaux_size(Integer(size))
            else:
                if size % 2 == 1 or size == 0:
                    return SymmetricTreeLikeTableaux_size(Integer(size))
                else:
                    raise ValueError("size must be 0 or odd")

    def _repr_(self):
        """
        """
        return "Factory for tree-like tableaux"

#################################################################
# Enumerated set of all full tree-like tableaux of a given size
#################################################################

class TreeLikeTableaux_size(TreeLikeTableaux):
    r"""
    The tree-like tableaux of size `n`.
    """
    def __init__(self, size):
        r"""
        """
        self._size = size
        super(TreeLikeTableaux_size, self).__init__(facade=TreeLikeTableaux_all(),
                                               category=FiniteEnumeratedSets())

    def _repr_(self):
        r"""
        """
        return "Tree-like tableaux of size %s"%(self._size)

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
            raise ValueError, "The tree-like tableau have a wrong size : %s"%(el.size())

    def __iter__(self):
        r"""
        """
        n = self._size

        if n==0:
            yield TreeLikeTableau([[]])
        if n==1:
            yield TreeLikeTableau([[1]])
        else:
            for tlt in TreeLikeTableaux(n-1):
                for i in range(n):
                    with tlt.clone() as tlt1:
                        tlt1.insert_point(i)
                    yield tlt1

    def options(self, *get_value, **set_values):
        return self._options( *get_value, **set_values)

    def _get_options(self):
        return self._options

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: TLT1 = TreeLikeTableaux(1)
            sage: TLT1([[]])   # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Wrong size : a tree-like tableau of size 1 is expected and one of size 0 was given
            sage: TLT1([[1]])   # indirect doctest
            [[1]]
            sage: TLT0 = TreeLikeTableaux(0)   # indirect doctest
            sage: TLT0([[]])
            [[]]
        """
        res = TreeLikeTableau(*args, **keywords)
        if res.size() != self._size:
            raise ValueError("Wrong size : a tree-like tableau of size %d is expected and one of size %d was given"%(self._size,res.size()))
        return res

    def random_element(self):
        r"""
        Return a random tree-like tableau  with uniform probability.

        This method generates a random permutation of the same size of self
        and convert it to a tree-like tableau.

        EXAMPLES::

            sage: TreeLikeTableaux(5).random_element() # random
            [[1, 1, 0, 0], [0, 1, 0, 1], [1, 0]]
            sage: TreeLikeTableaux(0).random_element()
            [[]]
            sage: TreeLikeTableaux(1).random_element()
            [[1]]

        TESTS::

            sage: TLT = TreeLikeTableaux(19)
            sage: all([TLT.random_element() in TLT for i in range(20)])
            True
        """
        from sage.combinat.permutation import Permutations
        if self._size == 0:
            return TreeLikeTableau([[]])
        return TreeLikeTableau.from_permutation(Permutations(self._size).random_element())

    def cardinality(self):
        r"""
        The cardinality of self.

        It is equal to size!

        TESTS::

            sage: TreeLikeTableaux(0).cardinality()
            1
            sage: TreeLikeTableaux(5).cardinality()
            120
            sage: TreeLikeTableaux(6).cardinality()
            720
        """
        return factorial(self._size)
    _options = TreeLikeTableauxOptions



#################################################################
# Enumerated set of all full tree-like tableaux
#################################################################

class TreeLikeTableaux_all( TreeLikeTableaux, DisjointUnionEnumeratedSets ):
    r"""
    This class enumerate all the tree-like tableaux.
    """
    def __init__(self):
        r"""
        """
        DisjointUnionEnumeratedSets.__init__(
            self,
            Family(NonNegativeIntegers(), TreeLikeTableaux_size),
            facade=True, keepkey = False
        )

    def _repr_(self):
        r"""
        """
        return "Tree Like Tableaux"

    def __contains__(self, x):
        """
        TESTS::

            sage: TLTS = TreeLikeTableaux()
            sage: 1 in TLTS
            False
            sage: TLTS([[]]) in TLTS
            True
        """
        if isinstance(x, self.element_class):
            return True

    def __call__(self, x=[[]], *args, **keywords):
        """
        Ensure that ``[[]]`` instead of ``0`` is passed by default.

        TESTS::

            sage: TLTS = TreeLikeTableaux()
            sage: TLTS()
            [[]]
        """
        return super(TreeLikeTableaux, self).__call__(x, *args, **keywords)

    def check_element(self, el, check):
        r"""
        TESTS::
            sage: TLTS = TreeLikeTableaux()
            sage: TLTS.check_element(TLTS.an_element(), True)
        """
        pass

    def options(self, *get_value, **set_values):
        return self._options( *get_value, **set_values)

    def _get_options(self):
        return self._options

    _options = TreeLikeTableauxOptions
    Element = TreeLikeTableau

class SymmetricTreeLikeTableau(TreeLikeTableau):
    """
    A symmmetric tree-like tableau.
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, value, check=True):
        r"""
        """
        if isinstance(value, SymmetricTreeLikeTableau):
            return value
        return SymmetricTreeLikeTableaux_all().element_class(SymmetricTreeLikeTableaux_all(), value)

        # value is not a a symmetric tree-like tableau so we give an appropriate error message

    def __init__(self, parent, value, check=True):
        """
        """
        super(SymmetricTreeLikeTableau, self).__init__(parent, value, check)
        if check:
            if not self.is_symmetric():
                raise ValueError("the tree-like tableau is not symmetric")

    def special_point_sym( self ):
        r"""
        EXAMPLES::

            sage: tlt = SymmetricTreeLikeTableau([[1,0,1,1],[0,0,0,1],[1,0,0],[1,1]])
            sage: tlt.special_point_sym()
            1
        """
        if not self.is_symmetric:
            raise ValueError, "The tree-like tableau is not symmetric."

        border_cells = self.border_cells()
        for k in range((len(border_cells)-1)/2-1,0,-1):
            (irow, icol) = border_cells[k]
            if self[irow][icol] == 1 and border_cells[k-1][0] == irow:
                return k
        return 0

    def insert_point_sym(self, position, epsilon ):
        r"""
        EXAMPLES::

            sage: tlt = TreeLikeTableau([[1,0,1,1],[0,0,0,1],[1,0,0],[1,1]])
        """
        self._require_mutable()

        [irow, icol, direction] = self._border_position_to_coordinates_and_direction(position)
        position_special_point_sym=self.special_point_sym()

        if direction == 'horizontal':
            self._insert_column(icol, irow)
            self._insert_row(irow, icol)
        else:
            self._insert_row(icol, irow)
            self._insert_column(irow+1, icol)
        if epsilon == 1:
            if position < position_special_point_sym:
                self._insert_ribbon(position + 1, position_special_point_sym + 1)
                self._insert_ribbon(self.size()-(position_special_point_sym+1),self.size()-(position+1))
        else:
                self._insert_ribbon(position+1,self.size()-(position+1))

    def __copy__(self):
        r"""
        Return a copy of the object
        """
        res = SymmetricTreeLikeTableau( deepcopy( list(self) ) )
        res._set_mutable()
        return res



#################################################################
# Enumerated set of all symmetric tree-like taleaux
#################################################################

class SymmetricTreeLikeTableaux_all(DisjointUnionEnumeratedSets, TreeLikeTableaux):
    """
    All symmetric tree-like taleaux.
    """
    def __init__(self):
        """
        TESTS::

            sage: from sage.combinat.tree_like_tableaux import SymmetricTreeLikeTableaux_all
            sage: STLTS = SymmetricTreeLikeTableaux_all()
            sage: STLTS.cardinality()
            +Infinity

            sage: it = iter(STLTS)
            sage: STLTS([[]])
            [[]]

            sage: STLTS is TreeLikeTableaux(symmetric=True)
            True
            sage: TestSuite(FB).run() # long time
        """
#            sage: (next(it), next(it), next(it), next(it), next(it))
#            (., [., .], [[., .], [., .]], [[., .], [[., .], [., .]]], [[[., .], [., .]],  [., .]])
#            sage: next(it).parent()
#            Binary trees
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), SymmetricTreeLikeTableaux_size),
            facade=True, keepkey=False)


    def _repr_(self):
        """
        TEST::

            sage: TreeLikeTableaux(symmetric=True)   # indirect doctest
            Symmetric tree-like tableaux
        """
        return "Symmetric tree-like tableaux"

    def __contains__(self, x):
        """
        TESTS::

            sage: STLTS = TreeLikeTableaux(symmetric=True)
            sage: 1 in STLTS
            False
            sage: STLTS([[]]) in STLTS
            True
        """
        return isinstance(x, TreeLikeTableau) and x.is_symmetric()

    def __call__(self, x=[[]], *args, **keywords):
        """
        Ensure that ``[[]]`` instead of ``0`` is passed by default.

        TESTS::

            sage: STLTS = TreeLikeTableaux(symmetric=True)
            sage: STLTS()
            [[]]
        """
        return super(TreeLikeTableaux, self).__call__(x, *args, **keywords)

    Element = SymmetricTreeLikeTableau

#################################################################
# Enumerated set of full binary trees of a given size
#################################################################

class SymmetricTreeLikeTableaux_size(TreeLikeTableaux):
    """
    Symmetric tree-like tableaux with a fixed size.
    """
    def __init__(self, size):
        r"""
        TESTS::

            sage: from sage.combinat.tree_like_tableaux import SymmetricTreeLikeTableaux_size
            sage: for i in range(1,6):
            ....:     TestSuite(TreeLikeTableaux(2*i-1, symmetric=True)).run()
        """
        super(SymmetricTreeLikeTableaux_size, self).__init__(category=FiniteEnumeratedSets())
        self._size = Integer(size)

    def _repr_(self):
        r"""
        TESTS::

            sage: TreeLikeTableaux(3, symmetric=True)
            Symmetric tree-like tableaux of size 3
        """
        return "Symmetric tree-like tableaux of size %s" % self._size

    def __contains__(self, x):
        r"""
        TESTS::

            sage: STLT3 = TreeLikeTableaux(3, symmetric=True)
            sage: 1 in STLT3
            False
            sage: STLT3([[1,1], [1]]) in STLT3
            True
            sage: TreeLikeTableau([[1,1], [1]]) in STLT3
            True
            sage: TreeLikeTableau([[1,1]]) in STLT3
            False
        """
        return (isinstance(x, TreeLikeTableau)
                and x.size() == self._size
                and x.is_symmetric())

    def _an_element_(self):
        """
        TESTS::

            sage: TreeLikeTableaux(0, symmetric=True).an_element()  # indirect doctest
            [[]]
        """
        return self.first()

    def cardinality(self):
        r"""
        The cardinality of self.

        TESTS::

            sage: TreeLikeTableaux(0, symmetric=True).cardinality()
            1
            sage: TreeLikeTableaux(5, symmetric=True).cardinality()
            8
            sage: TreeLikeTableaux(11, symmetric=True).cardinality()
            3840
        """
        if self._size == 0:
            return Integer(1)
        n = (self._size-1)/2
        return Integer(2**n*factorial(n))

    def random_element(self):
        r"""
        Return a random symmetric tree-like tableau  with uniform probability.

        This method generates a random permutation of the same size of self
        and convert it to a tree-like tableau.

        EXAMPLES::

            sage: TreeLikeTableaux(5, symmetric=True).random_element() # random
            [[1, 0, 1], [0, 0, 1], [1, 1]]]
            sage: TreeLikeTableaux(0, symmetric=True).random_element()
            [[]]
            sage: TreeLikeTableaux(1, symmetric=True).random_element()
            [[1]]

        TESTS::

            sage: STLT = TreeLikeTableaux(19, symmetric=True)
            sage: all([STLT.random_element() in STLT for i in range(20)])
            True
        """
        from sage.combinat.permutation import Permutations
        if self._size == 0:
            return SymmetricTreeLikeTableau([[]])

        if self._size == 1:
            return SymmetricTreeLikeTableau([[1]])

        n = Integer((self._size-1)/2)
        p = Permutations(n).random_element()
        tlt=SymmetricTreeLikeTableau([[1]])
        tlt._set_mutable()
        for i in range(0,len(p)):
            cmp=0
            pi=p[i]
            for j in range(i+1):
                if p[j]<pi:
                    cmp+=1
                j+=1
            epsilon = Set([-1,1]).random_element()
            tlt.insert_point_sym(cmp, epsilon)
        tlt.set_immutable()
        return tlt

    def __iter__(self):
        """
        """
        n = self._size

        if n==0:
            yield SymmetricTreeLikeTableau([[]])
        if n==1:
            yield SymmetricTreeLikeTableau([[1]])
        else:
            for tlt in TreeLikeTableaux(n-2,symmetric=True):
                for i in range((n-1)/2):
                    with tlt.clone() as tlt1:
                        tlt1.insert_point_sym(i,1)
                    yield tlt1

                    with tlt.clone() as tlt1:
                        tlt1.insert_point_sym(i,-1)
                    yield tlt1

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: STLT0 = TreeLikeTableaux(0,symmetric=True)
            sage: STLT0([[]])
            [[]]
            sage: STLT3 = TreeLikeTableaux(3,symmetric=True)
            sage: STLT3._element_constructor_([[1,1],[1,0]])
            [[1, 1], [1, 0]]
            sage: STLT3([[1,1],[1]])
            [[1, 1], [1]]
            sage: STLT3([[1, 1]])
            Traceback (most recent call last):
            ...
            ValueError: not symmetric
            sage: STLT3([[1]])
            Traceback (most recent call last):
            ...
            ValueError: Wrong size : a symmetric tree-like tableau of size 3 is expected and one of size 1 was given
        """
        res = TreeLikeTableau(*args, **keywords)
        if not res.is_symmetric():
            raise ValueError("not symmetric")
        if res.size() != self._size:
            raise ValueError("Wrong size : a symmetric tree-like tableau of size %d is expected and one of size %d was given"%(self._size,res.size()))
        return res
