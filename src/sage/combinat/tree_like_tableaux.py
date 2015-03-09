r"""
The Tree-like tableaux
=====================

The goal of this module is to give some tools to manipulate the
tree-like tableaux.
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

from sage.structure.element_wrapper import ElementWrapper
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.family import Family

from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.structure.list_clone import ClonableList
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

class TreeLikeTableau( ClonableList ):
    r"""
    The class of Tree-Like tableaux.

    Ref: J.C. Aval, A. Boussicault, P. Nadeau, Tree-like tabelau,
    arXiv:1109.0371v2
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        r"""
        """
        return cls._auto_parent._element_constructor_( *args, **opts )

    @lazy_class_attribute
    def _auto_parent(cls):
        """
        """
        return TreeLikeTableaux()

    def __copy__(self):
#TODO : choisir une solution pour la fonction copy, parce que si on herite de cette classe, et qu'on ne surcharge pas __copy__ l'objet obtenue sera de la classe TreeLikeTableau. Donc il faudrait surcharger la fonction copy a chaque heritage. Une solution alternative serait de rajouter une option dans le constructeur pour savoir si on cree un objet mutable ou pas a utiliser en duo avec x.__class__ mais il faut que tous les contructeurs suivent une syntaxe particuliere.
        r"""
        Return a copy of the object
        """
        from copy import deepcopy
        res = TreeLikeTableau( deepcopy( list(self) ) )
        res._set_mutable()
        return res

    def __hash__(self):
        return hash(tuple(map(tuple,list(self))))

    def _check_it_is_a_list_of_lists_of_0s_1s( self ):
        r"""
        Check if the the elements of the list that was given, contains list
        with only only O's and 1's inside them.
        TESTS::
            sage: tlt = TreeLikeTableau( [[1,0,1,1],5,[0,1]] )
            Traceback (most recent call last):
            ...
            ValueError: The element of the list indexed by 1, i.e 5 , is not a list.
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
            if type(element) != list:
                raise ValueError,\
                "The element of the list indexed by %s, i.e %s , is not a list."\
                 %(irow,element)
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
        if not ( 1 in ClonableList.__getitem__( self, irow )):
            raise ValueError, "The row %s is empty."%(irow)

    def _check_cell_is_not_a_forbidden_pattern( self, irow, icol ):
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
        ClonableList.__init__( self, parent, value )
        if check:
            if not isinstance( value, (list, tuple) ) :
                raise ValueError, "Value %s must be a list or a tuple.i"%(value)
            self.check()

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
        Give the coodinates of the cell outside the tree-like tableau having the border labelled by 'position' as an edge and and the direction of this border, i.e 'horizontal' or 'vertical'.
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
        border number 'position' ((the border are numbered from left to right)

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
        self.set_list(lst)

    def _repr_(self):
        r"""
        The text representation of a tree-like tableau

        EXAMPLES::
            sage: tlt = TreeLikeTableau( [[1,0,1,0],[1,1,0,1],[0,1],[0,1],[1,0]] )
            sage: tlt
            [[1, 0, 1, 0], [1, 1, 0, 1], [0, 1], [0, 1], [1, 0]]
        """
        return self.parent().global_options.dispatch(self, '_repr_', 'display')

    def _repr_list( self ):
        return ClonableList._repr_(self)

    def _repr_drawing( self ):
        NotImplemented

    def get_tikz_options():
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
            False


        """
        if self == TreeLikeTableau([[]]):
            return False
        i=0
        while i+1<=len(self[i]):
            for j in range(i+1,len(self[i])+1):
                if self[i][j]!=self[j][i]:
                    return False
            i=i+1
        return True

    def special_point_sym( self ):
        r"""
        EXAMPLES::

            sage: tlt=TreeLikeTableau([[1,0,1,1],[0,0,0,1],[1,0,0],[1,1]])
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

    def insert_point_sym( self, position, epsilon ):
        r"""
        EXAMPLES::

            sage: tlt=TreeLikeTableau([[1,0,1,1],[0,0,0,1],[1,0,0],[1,1]])
        """
        if not self.is_symmetric:
            raise ValueError, "The tree-like tableau is not symmetric."

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

    def code( self ):
        r"""
        The code of a tree-like tableau T is the sequence of the indexes of
        the border edges at which the consecutive insertions were done in order to
        obtain T starting with the tree-like tableau of size 1.

        EXAMPLES::
        """
        L=[]
        tlt=self.__copy__()
        for i in range(tlt.size()):
            L.insert(0,tlt.remove_point())
        return L

    def to_permutations_by_the_code( self ):
        NotImplemented
class TreeLikeTableauxFactory(SetFactory):
    r"""
    The tree-like tableaux factory.
    """
    def __call__(self, size=None, policy=None):
        r"""
        """
        if policy is None:
            policy = self._default_policy

        if isinstance(size, (Integer, int)):
            return TreeLikeTableaux_size(size, policy)
        if size is None:
            return TreeLikeTableaux_all(policy)
        raise ValueError, "Invalide argument for Tree-like tableaux Factory."

    def add_constraints(self, cons, (args, opts)):
        r"""
        """
        return cons+args

    @lazy_attribute
    def _default_policy(self):
        return TopMostParentPolicy(self, (), TreeLikeTableau)

    def _repr_(self):
        """
        """
        return "Factory for tree-like tableaux"

TreeLikeTableaux = TreeLikeTableauxFactory()
TreeLikeTableaux.__doc__ = TreeLikeTableauxFactory.__call__.__doc__


class TreeLikeTableaux_size(ParentWithSetFactory, UniqueRepresentation):
    r"""
    The tree-like tableaux of size `n`.
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
        return "tree-like tableaux of size %s"%(self._size)

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
        if self._size == 0:
            yield TreeLikeTableau([[]])
        elif self._size == 1:
            yield TreeLikeTableau([[1]])
        else:
            for tlt in TreeLikeTableaux(self._size - 1):#meilleure notation a trouver
                for k in range(self._size):
                    tlt._set_mutable()
                    tlt.insert_point(k)
                    tlt.set_immutable()
                    yield tlt
                    tlt.remove_point()


    global_options = TreeLikeTableauxOptions




class TreeLikeTableaux_all( ParentWithSetFactory, DisjointUnionEnumeratedSets ):
    r"""
    This class enumerate all the tree-like tableaux.
    """
    def __init__(self, policy):
        r"""
        """
        ParentWithSetFactory.__init__(
            self, (), policy, category = FiniteEnumeratedSets()
        )
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), self._tree_like_tableaux_size),
            facade=True, keepkey = False,
            category = self.category()
        )

    def _tree_like_tableaux_size( self, n ):
        return TreeLikeTableaux_size( n, policy=self.facade_policy() )

    def _repr_(self):
        r"""
        """
        return "TreeLikeTableaux_all"

    def check_element(self, el, check):
        r"""
        TESTS::
            sage: TLTS = TreeLikeTableaux()
            sage: TLTS.check_element(TLTS.an_element(), True)
        """
        pass

    global_options = TreeLikeTableauxOptions
