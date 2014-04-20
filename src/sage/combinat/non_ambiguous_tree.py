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
from sage.combinat.combinat import catalan_number
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.binary_tree import BinaryTrees, BinaryTree
from sage.combinat.binary_tree import LabelledBinaryTree, LabelledBinaryTrees
from sage.functions.other import factorial
from sage.combinat.posets.posets import Poset

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
        return Poset( [ self.vertices, self.edges ], cover_relations=True )

def _rec_posets_of_nodes(
    tree, left_poset, right_poset, path=(), lfather=None, rfather=None
):
    r"""
    This function caclulate the left poset and the right poset of the node
    of a tree.

    a right node r1 is lesser than  another right node r2 if r1 is a parent of
    r2.

    EXAMPLES::
        sage: from sage.combinat.non_ambiguous_tree import _poset_s
        sage: from sage.combinat.non_ambiguous_tree import _rec_posets_of_nodes
        sage: left_poset = _poset_s()
        sage: right_poset = _poset_s()
        sage: T = BinaryTree( [ [None,[[],None]],[[None,[]],None] ] )
        sage: _rec_posets_of_nodes( T, left_poset, right_poset )
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
        sage: _rec_posets_of_nodes( T, left_poset, right_poset )
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
        _rec_posets_of_nodes(
            son, left_poset, right_poset, son_path, son_path, rfather
        )
    if tree[1] != BinaryTree():
        son = tree[1]
        son_path = path+(1,)
        right_poset.vertices.append( son_path )
        if not rfather is None:
            right_poset.edges.append( [rfather, son_path] )
        _rec_posets_of_nodes(
            son, left_poset, right_poset, son_path, lfather, son_path
        )


def _posets_of_nodes( tree ):
    left_poset = _poset_s()
    right_poset = _poset_s()
    _rec_posets_of_nodes( tree, left_poset, right_poset )
    return [ left_poset.to_poset(), right_poset.to_poset() ]

class NonAmbiguousTree( ClonableList ):
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

        """
        return _posets_of_nodes( self.get_tree() )

    def get_tree( self ):
        return self[0]

    def _recursive_check(
        self, tree, llabels, rlabels, lfather=None, rfather=None
    ):
        if tree == LabelledBinaryTree( None ):
            return
        if tree[0] != LabelledBinaryTree( None ):
            if not lfather is None:
                if tree[0].label() <= lfather.label():
                    raise ValueError, "This is not a valid Non-ambiguous Tree."
                llabels.append( tree[0].label() )
            self._recursive_check( tree[0], llabels, rlabels, tree[0], rfather )
        if tree[1] != LabelledBinaryTree( None ):
            if not rfather is None:
                if tree[1].label() <= rfather.label():
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
            ....:     ], 0 )
            ....: )
            sage: nat = NonAmbiguousTree(
            ....:     LB( None )
            ....: )
            sage: nat = NonAmbiguousTree(
            ....:     LB( [], 0 )
            ....: )

        """
        llabels = []
        rlabels = []
        self._recursive_check( self.get_tree(), llabels, rlabels )
        if len(llabels) != len( Set(llabels) ):
            raise ValueError, "This is not a valid Non-ambiguous Tree."
        if len(rlabels) != len( Set(rlabels) ):
            raise ValueError, "This is not a valid Non-ambiguous Tree."
        if not self.get_tree().is_empty() and self.get_tree().label() != 0:
            raise ValueError, "The label of the root have to 0 and not %s."%(
                self.get_tree().label()
            )

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

        By convention, we will label the root by 0.

        Ref : Combinatorics of non-ambiguous trees, J.C. Aval, A. Boussicault,
              M. Bouvel, M. Silimbani, arXiv:1305.3716

        EXAMPLES::

            sage: nat = NonAmbiguousTree(
            ....:     [ [1,1,0,1,0], [1,0,1,0,1], [0,1,0,0,0] ]
            ....: )
            sage: nat
            0[1[., 2[., 4[., .]]], 1[2[., .], 3[., .]]]

            sage: LB = LabelledBinaryTree
            sage: nat = NonAmbiguousTree(
            ....:     LB( [
            ....:         LB( [None, LB([None, LB([], 4)],2)] ,1),
            ....:         LB( [LB([],2), LB([],3) ], 1)
            ....:     ], 0)
            ....: )
            sage: nat
            0[1[., 2[., 4[., .]]], 1[2[., .], 3[., .]]]
        """
        def _recursive_binary_tree( array, position=[0,0], node_type=None ):
            height = len( array )
            width = len( array[0] )
            if width == 0:
                return LabelledBinaryTree( None )
            lposition = [ position[0]+1, position[1] ]
            while(
                lposition[0] < height
                and array[ lposition[0] ][ lposition[1] ] == 0
            ):
                lposition[0] += 1
            if lposition[0] < height:
                ltree = _recursive_binary_tree( array, lposition, 0 )
            else:
                ltree = None
            rposition = [ position[0], position[1]+1 ]
            while(
                rposition[1] < width
                and array[ rposition[0] ][ rposition[1] ] == 0
            ):
                rposition[1] += 1
            if rposition[1] < width:
                rtree = _recursive_binary_tree( array, rposition, 1 )
            else:
                rtree = None
            if node_type is None:
                return LabelledBinaryTree( [ltree, rtree], label=0 )
            return LabelledBinaryTree(
                [ltree, rtree], label=position[node_type]
            )

        if check:
            if value in LabelledBinaryTrees():
                tree = value
            elif isinstance( value, (list, tuple) ) :
                tree = _recursive_binary_tree( value )
            else:
                raise ValueError, "Value %s must be a list or a tuple."%(value)
            ClonableList.__init__(self, parent, [tree] )
            self.check()
        self._options = None

    def left_node_number( self ):
        return self.get_tree().left_node_number()

    def right_node_number( self ):
        return self.get_tree().right_node_number()

    def width( self ):
        r"""
        Return the width of the non-ambiguous tree.

        EXAMPLES::

            sage: nat = NonAmbiguousTree(
            ....:     [ [1,1,0,1,0], [1,0,1,0,1], [0,1,0,0,0] ]
            ....: )
            sage: nat.width()
            5
            sage: nat = NonAmbiguousTree( [[1]] )
            sage: nat.width()
            1
            sage: nat = NonAmbiguousTree( [[]] )
            sage: nat.width()
            0
        """
        if self.get_tree().is_empty():
            return 0
        return 1 + self.right_node_number()

    def height( self ):
        r"""
        Return the height of the non-ambiguous tree.

        EXAMPLES::

            sage: nat = NonAmbiguousTree(
            ....:     [ [1,1,0,1,0], [1,0,1,0,1], [0,1,0,0,0] ]
            ....: )
            sage: nat.height()
            3
            sage: nat = NonAmbiguousTree( [[1]] )
            sage: nat.height()
            1
            sage: nat = NonAmbiguousTree( [[]] )
            sage: nat.height()
            0
        """
        if self.get_tree().is_empty():
            return 0
        return 1 + self.left_node_number()

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
            ....:     ], 0 )
            ....: )
            sage: nat.get_array()
            [[1, 1, 0, 1, 0], [1, 0, 1, 0, 1], [0, 1, 0, 0, 0]]

            sage: nat = NonAmbiguousTree( [[1]] )
            sage: nat.get_array()
            [[1]]

            sage: nat = NonAmbiguousTree( [[]] )
            sage: nat.get_array()
            [[]]
        """
        if self.get_tree().is_empty():
            return [[]]
        res = [
            [ 0 for w in range(self.width()) ]
            for h in range(self.height())
        ]
        res[0][0] = 1
        def init_array( array, tree, lfather=0, rfather=0 ):
            if tree.is_empty():
                return
            if not tree[0].is_empty():
                array[ tree[0].label() ][ lfather ] = 1
                init_array( array, tree[0], lfather, tree[0].label() )
            if not tree[1].is_empty():
                array[ rfather ][ tree[1].label() ] = 1
                init_array( array, tree[1], tree[1].label(), rfather )
        init_array( res, self.get_tree() )
        return res

    def _repr_(self):
        r"""
        Return a string representation of the non-ambious tree.

        EXAMPLES::
            sage: nat = NonAmbiguousTree(
            ....:     [ [1,1,0,1,0], [1,0,1,0,1], [0,1,0,0,0] ]
            ....: )
            sage: nat
            0[1[., 2[., 4[., .]]], 1[2[., .], 3[., .]]]
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
            '0[1[., 2[., 4[., .]]], 1[2[., .], 3[., .]]]'
        """
        return self.get_tree()._repr_()

    def _repr_drawing( self ):
        r"""
        Return a string representing a drawing of the non-ambiguous tree.

        EXAMPLES::
            sage: nat = NonAmbiguousTree(
            ....:     [ [1,1,0,1,0], [1,0,1,0,1], [0,1,0,0,0] ]
            ....: )
            sage: nat._repr_drawing()
            '[1 1 0 1 0]\n[1 0 1 0 1]\n[0 1 0 0 0]'
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

    def get_options( self ):
        r"""
        Return all the opitons of the object.
        """
        if self._options is None:
            return self.parent().get_options()
        return self._options

    def set_options( self, *get_value, **set_value ):
        r"""
        Set new options to the object.
        """
        if self._options is None:
            self._options = deepcopy( self.get_options() )
        self._options( *get_value, **set_value )


class NonAmbiguousTreesFactory(SetFactory):
    r"""
    The non-ambiguous trees factory.
    """
    def __call__(self, size=None, tree=None, policy=None):
        r"""
        """
        if policy is None:
            policy = self._default_policy

        if not tree is None:
            return NonAmbiguousTrees_binarytree(tree, policy)
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


class NonAmbiguousTrees_binarytree(ParentWithSetFactory, UniqueRepresentation):
    r"""
    The non-ambiguous tree of a given binary tree.
    """
    def __init__(self, btree, policy):
        r"""
        Construct a set of non-ambiguous trees of a given binary tree.
        """
        self._btree = btree
        ParentWithSetFactory.__init__(
            self, (btree,), policy, category = FiniteEnumeratedSets()
        )

    def _repr_(self):
        r"""
        Return the string representation of the set of non-ambiguous trees.

        EXAMPLES::

            sage: NonAmbiguousTrees( tree = BinaryTree([[],[]]) )
            non-ambiguous trees of binary tree [[., .], [., .]]
        """
        return "non-ambiguous trees of binary tree %s"%(self.tree())

    def an_element(self):
        r"""
        """
        return next( self.__iter__() )

    def check_element(self, el, check):
        r"""
        """
        if BinaryTree( el.get_tree() ) == self.tree():
            raise ValueError, "The non-ambiguous trees have a Wrong binary tree : %s"%(
                BinaryTree(el.get_tree())
            )

    def cardinality( self ):
        r"""
        Return the number of non-ambiguous trees.

        EXAMPLES::

            sage: NonAmbiguousTrees(
            ....:     tree = BinaryTree( [[[], []], [[None, []], None]] )
            ....: ).cardinality()
            9
            sage: NonAmbiguousTrees( tree = BinaryTree( [] ) ).cardinality()
            1
            sage: NonAmbiguousTrees( tree = BinaryTree() ).cardinality()
            1

        TESTS::
            sage: NonAmbiguousTrees(0).cardinality()
            1
            sage: NonAmbiguousTrees(1).cardinality()
            1
            sage: NonAmbiguousTrees(2).cardinality()
            2
            sage: NonAmbiguousTrees(3).cardinality()
            5
            sage: NonAmbiguousTrees(4).cardinality()
            16
            sage: NonAmbiguousTrees(5).cardinality()
            63
            sage: NonAmbiguousTrees(6).cardinality()
            294
            sage: NonAmbiguousTrees(7).cardinality()
            1585
            sage: NonAmbiguousTrees(8).cardinality()
            9692
        """
        def rec_cardinality( tree ):
            res = [ 0, 0, 1, 1 ]
            if tree.is_empty():
                return res
            if not tree[0].is_empty() :
                lres = rec_cardinality( tree[0] )
                res[0] += ( 1 + lres[0] )
                res[1] += lres[1]
                res[2] *= ( lres[2] * ( 1 + lres[0] ) )
                res[3] *= lres[3]
            if not tree[1].is_empty() :
                rres = rec_cardinality( tree[1] )
                res[0] += rres[0]
                res[1] += ( 1 + rres[1] )
                res[2] *= rres[2]
                res[3] *= ( rres[3] * ( 1 + rres[1] ) )
            return res
        lnode_nb, rnode_nb, prod_lnode_nb, prod_rnode_nb = rec_cardinality(
            self.tree()
        )
        return (
            ( factorial( lnode_nb )/prod_lnode_nb )
            * ( factorial( rnode_nb )/prod_rnode_nb )
        )

    def __iter__(self):
        r"""
        Rerturn a non-ambiguous trees generator.

        EXAMPLES::
            sage: btree = BinaryTree( [[[], []], [[None, []], None]] )
            sage: NAS = NonAmbiguousTrees( tree = btree )
            sage: generator = iter( NAS )
            sage: Set( generator) == Set( [
            ....:     NonAmbiguousTree( [[1,0,1,0], [0,0,1,1], [1,1,0,0], [1,0,0,0]] ),
            ....:     NonAmbiguousTree( [[1,1,0,0], [0,1,0,1], [1,0,1,0], [1,0,0,0]] ),
            ....:     NonAmbiguousTree( [[1,1,0,0], [0,1,1,0], [1,0,0,1], [1,0,0,0]] ),
            ....:     NonAmbiguousTree( [[1,0,1,0], [1,1,0,0], [0,0,1,1], [1,0,0,0]] ),
            ....:     NonAmbiguousTree( [[1,1,0,0], [1,0,1,0], [0,1,0,1], [1,0,0,0]] ),
            ....:     NonAmbiguousTree( [[1,1,0,0], [1,0,0,1], [0,1,1,0], [1,0,0,0]] ),
            ....:     NonAmbiguousTree( [[1,0,1,0], [1,1,0,0], [1,0,0,0], [0,0,1,1]] ),
            ....:     NonAmbiguousTree( [[1,1,0,0], [1,0,1,0], [1,0,0,0], [0,1,0,1]] ),
            ....:     NonAmbiguousTree( [[1,1,0,0], [1,0,0,1], [1,0,0,0], [0,1,1,0]] ),
            ....: ] )
            True

            sage: btree = BinaryTree( [[], None] )
            sage: NAS = NonAmbiguousTrees( tree = btree )
            sage: generator = iter( NAS )
            sage: generator.next()
            0[1[., .], .]

            sage: btree = BinaryTree( [None, []] )
            sage: NAS = NonAmbiguousTrees( tree = btree )
            sage: generator = iter( NAS )
            sage: generator.next()
            0[., 1[., .]]

            sage: btree = BinaryTree( [] )
            sage: NAS = NonAmbiguousTrees( tree = btree )
            sage: generator = iter( NAS )
            sage: generator.next()
            0[., .]

            sage: btree = BinaryTree( )
            sage: NAS = NonAmbiguousTrees( tree = btree )
            sage: generator = iter( NAS )
            sage: generator.next()
            .
        """
        if self.tree().is_empty() :
            yield NonAmbiguousTree( LabelledBinaryTree( None ) )
        else:
            lposet, rposet =  _posets_of_nodes( self.tree() )
            for lw in lposet.linear_extensions():
                for rw in rposet.linear_extensions():
                    lbtree = LabelledBinaryTree( self.tree() )
                    lbtree._set_mutable()
                    for h in range( len(lw) ):
                        lbtree.set_label( lw[h], h+1 )
                    for w in range( len(rw) ):
                        lbtree.set_label( rw[w], w+1 )
                    lbtree.set_label( (), 0 )
                    lbtree.set_immutable()
                    yield NonAmbiguousTree( lbtree )

    def get_options( self ):
        return self.global_options

    def set_options( self, *get_value, **set_value ):
        self.global_options( *get_value, **set_value )

    def tree( self ):
        return self._btree

    global_options = NonAmbiguousTreesOptions


class NonAmbiguousTrees_size(ParentWithSetFactory, DisjointUnionEnumeratedSets):
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
        DisjointUnionEnumeratedSets.__init__(
            self, Family(
                BinaryTrees(size), self._non_ambiguous_trees_tree
            ),
            facade=True, keepkey = False,
            category = self.category()
        )

    def _non_ambiguous_trees_tree( self, btree ):
        return NonAmbiguousTrees_binarytree(
            btree, policy = self.facade_policy()
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
            raise ValueError, "The non-ambiguous trees have a Wrong size : %s"%(
                el.size()
            )

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

class NonAmbiguousTrees_all(ParentWithSetFactory, DisjointUnionEnumeratedSets):
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

            sage: NonAmbiguousTree( [[0,1,1],[1,1,0]] )  in NATS
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
