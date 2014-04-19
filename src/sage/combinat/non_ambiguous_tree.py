r"""
The Non-ambiguous trees
=====================

The goal of this module is to give some tools to manipulate the 
non-ambiguous trees.
"""
#*****************************************************************************
#  Copyright (C) 2014 Adrien Boussicault (boussica@labri.fr), 
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element_wrapper import ElementWrapper
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
from sage.combinat.combinat import catalan_number
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.binary_tree import BinaryTree
from sage.combinat.binary_tree import LabelledBinaryTree, LabelledBinaryTrees

default_non_ambiguous_tikz_options = dict(
    scale=1, line_size=1, point_size=3.5
    , color_line='black', color_point='black'
    , translation=[0,0], rotation=0
)

NonAmbiguousTreesOptions=GlobalOptions(
    name = 'Non-ambiguous Trees',
    doc=r"""
    """,
    end_doc=r"""
    """,
    tikz_options=dict(
        default= default_non_ambiguous_tikz_options,
        description='the tikz options',
        checker=lambda x: Set(x.keys()).issubset(
            Set( [
                'scale', 'line_size', 'point_size'
                , 'color_line', 'color_point', 'translation', 
                'rotation'
            ] )
        )
    ),
    drawing_components=dict(
        default= dict( tree=True ),
        description='Different tree-like tableaux components to draw',
        checker=lambda x: Set(x.keys()).issubset(
            Set( [
                'diagram', 'tree'
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

class _poset_s:
    def __init__( self ):
        self.vertices = []
        self.edges = []
    def to_poset( self ):
        return Poset( [ self.vertices, self.edges ], cover_relation=True )

def _posets_of_nodes(
    tree, left_poset, right_poset, path=(), lfather=None, rfather=None
):
    r"""
    This function caclulate the left poset and the right poset of the node
    of a tree.

    a right node r1 is lesser than  another right node r2 if r1 is a parent of
    r2.

    EXAMPLES::
        sage: from sage.combinat.non_ambiguous_tree import _poset_s 
        sage: from sage.combinat.non_ambiguous_tree import _posets_of_nodes 
        sage: left_poset = _poset_s()
        sage: right_poset = _poset_s()
        sage: T = BinaryTree( [ [None,[[],None]],[[None,[]],None] ] )
        sage: _posets_of_nodes( T, left_poset, right_poset )
        sage: left_poset.vertices
        [(0,), (0, 1, 0), (1, 0)]
        sage: left_poset.edges
        [[(0,), (0, 1, 0)]]
        sage: right_poset.vertices
        [(0, 1), (1,), (1, 0, 1)]
        sage: right_poset.edges
        [[(1,), (1, 0, 1)]]

    TESTS::

        sage: from sage.combinat.non_ambiguous_tree import _poset_s 
        sage: left_poset = _poset_s()
        sage: right_poset = _poset_s()
        sage: T = BinaryTree( [ [None,[[],[]]],[[[],[]],None] ] )
        sage: _posets_of_nodes( T, left_poset, right_poset )
        sage: left_poset.vertices
        [(0,), (0, 1, 0), (1, 0), (1, 0, 0)]
        sage: left_poset.edges
        [[(0,), (0, 1, 0)], [(1, 0), (1, 0, 0)]]
        sage: right_poset.vertices
        [(0, 1), (0, 1, 1), (1,), (1, 0, 1)]
        sage: right_poset.edges
        [[(0, 1), (0, 1, 1)], [(1,), (1, 0, 1)]]
    """
    if tree == BinaryTree():
        return
    if tree[0] != BinaryTree():
        son = tree[0]
        son_path = path+(0,)
        left_poset.vertices.append( son_path )
        if not lfather is None:
            left_poset.edges.append( [lfather, son_path] )
        _posets_of_nodes(
            son, left_poset, right_poset, son_path, son_path, rfather
        )
    if tree[1] != BinaryTree():
        son = tree[1]
        son_path = path+(1,)
        right_poset.vertices.append( son_path )
        if not rfather is None:
            right_poset.edges.append( [rfather, son_path] )
        _posets_of_nodes(
            son, left_poset, right_poset, son_path, lfather, son_path
        )

class NonAmbiguousTree( LabelledBinaryTree ):
    r"""
    The class of Non-ambiguous Trees.
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
        return NonAmbiguousTrees()

    def posets_of_nodes( self ):
        r"""
        Return a pair [pl, pr] of poset where pl is the poset of left nodes of 
        the tree and pr is the poset of right nodes.

        A left node r1 is lesser than another left node r2 if r1 is a parent of
        r2 in the tree.

        EXAMPLES::
        
            sage:  NotImplemented
        """
        left_poset = _poset_s()
        right_poset = _poset_s()
        _posets_of_nodes( self, left_poset, right_poset )
        return [ left_poset.to_poset(), right_poset.to_poset() ]

    def _recursive_check( tree, llabels, rlabels, lfather=None, rfather=None ):
        print "######"
        print tree
        if tree == LabelledBinaryTree( None ):
            return
        if tree[0] != LabelledBinaryTree( None ):
            if not lfather is None:
                if tree[0].label() >= lfather.label():
                    raise ValueError, "This is not a valid Non-ambiguous Tree." 
                llabels.append( tree[0].label() )
            self._recursive_check( tree[0], llabels, rlabels, tree[0], rfather )
        if tree[1] != LabelledBinaryTree( None ):
            if not rfather is None:
                if tree[1].label() >= rfather.label():
                    raise ValueError, "This is not a valid Non-ambiguous Tree." 
                rlabels.append( tree[1].label() )
            self._recursive_check( tree[1], llabels, rlabels, lfather, tree[1] )

    def check( self ):
        r"""
        This method raise an error if the internal data of the class doesn't
        represent a Non-ambiguous trees.

        EXAMPLES::

            sage: LB = LabelledBinaryTree
            sage: nat = NonAmbiguousTree(
            ....:     LB( [ 
            ....:         LB( [None, LB([None, LB([], 4)],2)] ,1),
            ....:         LB( [LB([],2), LB([],3) ], 1)
            ....:     ] )
            ....: )
            sage: nat = NonAmbiguousTree(
            ....:     LB( None )
            ....: )
            sage: nat = NonAmbiguousTree(
            ....:     LB( [], None )
            ....: )

        """
        llabels = []
        rlabels = []
        self._recursive_check( self, llabels, rlabels )
        if len(llabels) != len( Set(llabels) ):
            raise ValueError, "This is not a valid Non-ambiguous Tree." 
        if len(rlabels) != len( Set(rlabels) ):
            raise ValueError, "This is not a valid Non-ambiguous Tree." 

    def __init__(self, parent, value, check=True):
        r"""
        Construct a non-ambiguous tree.
        
        Two types of input is allowed :
        1) A Labelled binary trees.
            The label of left nodes contains the height position of the node.
            THe label of right nodes contains the width potition of the node.
            The height of a right node is the same as its father in the tree.
            The width of a left nodes is equals to the height of its father.
        2) An array containing of 1 and 0. The 1s represent the nodes of the the
            non-ambiguous tree and 0s nothing.
            Each row have to contain at least one 1.
            Each column have to contain at least one 1.
            There is a 1 situated at position (0,0) called the root.
            To each node, which are not the root, there is 
                a node on it left in the same row or 
                a node on the top in the same column 
                but not the both.

        For example :

        1 1 0 1 0
        1 0 1 0 1
        0 1 0 0 0

        is a non-ambiguous tree and can be represented by the following 
        non-ambguous tree:

                x
              /   \
             1      1
              \    / \
               2  2   3
                \
                 4

        If we replace each 1 of the previous matrix by the label of its node,
        in the non-ambiguous tree, we obtain :

                               Label of
                             right nodes
                               1 2 3 4
                               | | | |
                               V V V V

                             x 1 0 3 0
        Label of       1 ->  1 0 2 0 4         x is the root
        left nodes     2 ->  0 2 0 0 0         0 is nothing

        Ref : Combinatorics of non-ambiguous trees, J.C. Aval, A. Boussicault, 
              M. Bouvel, M. Silimbani, arXiv:1305.3716

        EXAMPLES::

            sage: nat = NonAmbiguousTree(
            ....:     [ [1,1,0,1,0], [1,0,1,0,1], [0,1,0,0,0] ]
            ....: )
            sage: nat
            None[1[., 2[., 4[., .]]], 1[ 2[., .], 3[., .] ]]

            sage: LB = LabelledBinaryTree
            sage: nat = NonAmbiguousTree(
            ....:     LB( [ 
            ....:         LB( [None, LB([None, LB([], 4)],2)] ,1),
            ....:         LB( [LB([],2), LB([],3) ], 1)
            ....:     ] )
            ....: )
            sage: nat
            None[1[., 2[., 4[., .]]], 1[ 2[., .], 3[., .] ]]
        """
        def _recursive_binary_tree( array, position=[0,0], node_type=None ):
            height = len( array )
            width = len( array[0] )
            if width == 0:
                return LabelledBinaryTree( None )
            lposition = [ position[0]+1, position[1] ]
            while(
                lposition[0] < height
                and array[ lposition[0], lposition[1] ] != 0
            ):
                lposition[0] += 1
            if lposition[0] < height:
                ltree = _recursive_binary_tree( array, lposition,  )
            else:
                ltree = None
            rposition = [ position[0], position[1]+1 ]
            while(
                rposition[1] < width
                and array[ rposition[0], rposition[1] ] != 0
            ):
                rposition[1] += 1
            if rposition[0] < width:
                rtree = _recursive_binary_tree( array, rposition )
            else:
                rtree = None
            if node_type is None:
                return LabelledBinaryTree( [ltree, rtree], label=None )
            return LabelledBinaryTree( [ltree, rtree], label=position[node_type] )

        if check:
            if value in LabelledBinaryTrees():
                print "@@@@111"
                print value
                print "@@@@222"
                super(NonAmbiguousTree, self).__init__(self, parent, value, check=check )
                #LabelledBinaryTree.__init__(self, parent, value, check=check)
                print "@@@@333"
                print( self )
                print "@@@@444"
            elif isinstance( value, (list, tuple) ) :
                LabelledBinaryTree.__init__(
                    self, parent, self._recursive_binary_tree( value )
                ) 
            else:
                raise ValueError, "Value %s must be a list or a tuple."%(value)
            self.check()
        self._options = None


    @cached_method
    def get_array( self ):
        r"""
        Return the array associated with the non ambiguous tree.
        
        EXAMPLES::

            sage: nat = NonAmbiguousTree(
            ....:     [ [1,1,0,1,0], [1,0,1,0,1], [0,1,0,0,0] ]
            ....: )
            sage: nat.get_array()
            [[1, 1, 0, 1, 0], [1, 0, 1, 0, 1], [0, 1, 0, 0, 0]]

            sage: LB = LabelledBinaryTree
            sage: nat = NonAmbiguousTree(
            ....:     LB( [ 
            ....:         LB( [None, LB([None, LB([], 4)],2)] ,1),
            ....:         LB( [LB([],2), LB([],3) ], 1)
            ....:     ] )
            ....: )
            sage: nat.get_array()
            [[1, 1, 0, 1, 0], [1, 0, 1, 0, 1], [0, 1, 0, 0, 0]]
        """
        NotImplemented

    def _repr_(self):
        r"""
        Return a string representation of the non-ambious tree.

        EXAMPLES::
            sage: nat = NonAmbiguousTree(
            ....:     [ [1,1,0,1,0], [1,0,1,0,1], [0,1,0,0,0] ]
            ....: )
            sage: nat
            None[1[., 2[., 4[., .]]], 1[ 2[., .], 3[., .] ]]
            sage: nat.set_options( display='drawing' )
            sage: nat
            [1 1 0 1 0]
            [1 0 1 0 1]
            [0 1 0 0 0]
        """
        return self.get_options().dispatch(self, '_repr_', 'display')

    def _repr_list( self ):
        r"""
        Return a string representation with list style.

        EXAMPLES::
            sage: nat = NonAmbiguousTree(
            ....:     [ [1,1,0,1,0], [1,0,1,0,1], [0,1,0,0,0] ]
            ....: )
            sage: nat._repr_list()
            'None[1[., 2[., 4[., .]]], 1[ 2[., .], 3[., .] ]]'
        """
        return LabelledBinaryTree._repr_(self)

    def _repr_drawing( self ):
        r"""
        Return a string representing a drawing of the non-ambiguous tree.

        EXAMPLES::
            sage: nat = NonAmbiguousTree(
            ....:     [ [1,1,0,1,0], [1,0,1,0,1], [0,1,0,0,0] ]
            ....: )
            sage: nat._repr_drawing()
            '[0 0 1 0 1 1]\n[1 1 0 0 1 0]'
        """
        return str( matrix( self.get_array() ) )

    def get_tikz_options( self ):
        return self.get_options()['tikz_options']

    def _to_tikz_diagram( self ):
        tikz_options = self.get_tikz_options()
        NotImplemented

    def _to_tikz_tree( self ):
        res = ""
        tikz_options = self.get_tikz_options()
        NotImplemented

    def to_tikz( self ):
        r"""
        Return the tikz code of the non-ambiguous tree.

        This code is the code present inside a tikz latex environemet.
        """
        res = ""
        drawing_components = self.get_options()['drawing_components']
        if 'diagram' in  drawing_components :
            res += self._to_tikz_diagram()
        if 'tree' in  drawing_components :
            res += self._to_tikz_tree()
        return res

    def _latex_(self):
        r"""
        Return a LaTeX version of ``self``.

        For more on the latex options, see 
        :meth:`NonAmbiguousTrees.global_options`.
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


class NonAmbiguousTreesFactory(SetFactory):
    r"""
    The non-ambiguous trees factory.
    """
    def __call__(self, size=None, policy=None):
        r"""
        """
        if policy is None:
            policy = self._default_policy

        if isinstance(size, (Integer, int)):
            return NonAmbiguousTrees_size(size, policy)
        if size is None:
            return NonAmbiguousTrees_all(policy)
        raise ValueError, "Invalide argument for non-ambiguous tee Factory."

    def add_constraints(self, cons, (args, opts)):
        r"""
        """
        return cons+args

    @lazy_attribute
    def _default_policy(self):
        return TopMostParentPolicy(self, (), NonAmbiguousTree)

    def _repr_(self):
        """
        """
        return "Factory for non-ambiguous trees."

NonAmbiguousTrees = NonAmbiguousTreesFactory()
NonAmbiguousTrees.__doc__ = NonAmbiguousTreesFactory.__call__.__doc__


class NonAmbiguousTrees_size(ParentWithSetFactory, UniqueRepresentation):
    r"""
    The non-ambiguous tree of size `n`.
    """
    def __init__(self, size, policy):
        r"""
        Construct a set of non-ambiguous trees of a given size.
        """
        self._size = size
        ParentWithSetFactory.__init__(
            self, (size,), policy, category = FiniteEnumeratedSets()
        )

    def _repr_(self):
        r"""
        Return the string representation of the set of non-ambiguous trees.

        EXAMPLES::
        
            sage: NonAmbiguousTrees( 3 )
            non-ambiguous trees of size 3
        """
        return "non-ambiguous trees of size %s"%(self._size)

    def an_element(self):
        r"""
        """
        return next( self.__iter__() )

    def check_element(self, el, check):
        r"""
        """
        if el.size() != self.size():
            raise ValueError, "The non-ambiguous trees have a Wrong size : %s"%(el.size())

    def cardinality( self ):
        r"""
        Return the number of non-ambiguous trees.

        EXAMPLES::

            sage: NonAmbiguousTrees(1).cardinality()
            1
            sage: NonAmbiguousTrees(2).cardinality()
            2
            sage: NonAmbiguousTrees(3).cardinality()
            5

            sage: all( [
            ....:     NonAmbiguousTrees(i).cardinality() 
            ....:     == catalan_number(i)
            ....:     for i in range( 6 )
            ....: ] )
            True

            sage: all( [
            ....:     NonAmbiguousTrees(i).cardinality() 
            ....:     == len( list( NonAmbiguousTrees(i) ) )
            ....:     for i in range( 6 )
            ....: ] )
            True
        """
        NotImplemented

    def __iter__(self):
        r"""
        Rerturn a non-ambiguous trees generator.
        
        EXAMPLES::
            sage: len( list( NonAmbiguousTrees(3) ) ) == 5
            True
            sage: all( [ 
            ....:     nat in NonAmbiguousTrees()
            ....:     for nat in NonAmbiguousTrees(3)
            ....: ] )
            True
        """
        NotImplemented

    def get_options( self ):
        return self.global_options

    def size( self ):
        r"""
        Return the size of the non-ambiguous trees generated by this parent.
        
        EXAMPLES::

            sage: NonAmbiguousTrees(0).size()
            0
            sage: NonAmbiguousTrees(1).size()
            1
            sage: NonAmbiguousTrees(5).size()
            5
        """
        return self._size

    def set_options( self, *get_value, **set_value ):
        self.global_options( *get_value, **set_value )

    global_options = NonAmbiguousTreesOptions

class NonAmbiguousTrees_all( ParentWithSetFactory, DisjointUnionEnumeratedSets ):
    r"""
    This class enumerate all the non-ambiguous trees.
    """
    def __init__(self, policy):
        r"""
        Construct the set af all the non-ambiguous trees.

        EXAMPLES::
        
            sage: NATS = NonAmbiguousTrees()
            sage: NATS
            non-ambiguous trees
        
            sage: NonAmbiguousTree( [[0,1,1],[1,1,0]] )  in NATSS
            True

            sage: NATS = NonAmbiguousTrees()
            sage: next( NATS.__iter__() ) in NATS
            True
        """
        ParentWithSetFactory.__init__(
            self, (), policy, category = FiniteEnumeratedSets()
        )
        DisjointUnionEnumeratedSets.__init__(
            self, Family(
                NonNegativeIntegers(), self._non_ambiguous_trees_size
            ),
            facade=True, keepkey = False,
            category = self.category()
        )

    def _non_ambiguous_trees_size( self, n ):
        return NonAmbiguousTrees_size( n, policy=self.facade_policy() )

    def _repr_(self):
        r"""
        Returns a string representation of the set of non-ambiguous trees.

        EXAMPLES::
        
            sage: NATS = NonAmbiguousTrees()
            sage: NATS
            non-ambiguous trees
        """
        return "non-ambiguous trees"

    def check_element(self, el, check):
        r"""
        Check is a given element `el` is in the set of non-ambiguous trees.

        EXAMPLES::

            sage: NATS = NonAmbiguousTrees()
            sage: NonAmbiguousTree( [[0,1,1],[1,1,0]] ) in NATS
            True
        """
        pass

    def get_options( self ):
        r"""
        Returns all the aptions associated with the set of non-ambiguous trees.
        
        EXAMPLES::

            sage: NATS = NonAmbiguousTrees()
            sage: options = NATS.get_options()
            sage: options
            options for Non-ambiguous Trees
            sage: options()
            Current options for Non-ambiguous Trees
              - display:            list
            ...
        """
        return self.global_options

    def set_options( self, *get_value, **set_value ):
        self.global_options( *get_value, **set_value )

    global_options = NonAmbiguousTreesOptions
