r"""
The Tree-like tableaux
=====================

The goal of this module is to give some tools to manipulate the 
tree-like tableaux.
"""
#*****************************************************************************
#  Copyright (C) 2014 Adrien Boussicault (boussica@labri.fr), 
#                     Patxi Laborde-Zubieta (patxi.laborde.zubieta@gmail.com)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element_wrapper import ElementWrapper
from sage.structure.list_clone import ClonableList
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.set_factories import SetFactory, SetFactoryParent, TopMostParentPolicy

class TreeLikeTableau( ClonableList ):
    r"""
    The class of Tree-Like tableaux.
    
    Ref: J.C. Aval, A. Boussicault, P. Nadeau, Tree-like tabelau, arXiv:1109.0371v2
    """

    def check():
        r"""
        Check if the internal data contain valid information for a tree-like 
        tableaux.

        EXAMPLES::
            sage: NotImplemented

        TESTS::
            sage: NotImplemented
        """
        NotImplemented

    def __init__(self, parent, value, check=True):
        r"""
        The class for a tree-like tableau element.
        
        EXAMPLES::
            sage: from sage.combinat.tree_like_tableaux import TreeLikeTableau
            sage: tlt = TreeLikeTableau( [[1,1,0,1],[1,0,1,0],[0,1],[0,1],[1,0]] )
            sage: tlt
            [[1, 1, 0, 1], [1, 0, 1, 0], [0, 1], [0, 1], [1, 0]]
            sage: tlt = TreeLikeTableau( [] )
            sage: tlt
            []
        """
        if check:
            if not isinstance( value, (list, tuple) ) :
                raise ValueError, "Value %s must be a list or a tuple.i"%(value)
            self.check()
        ClonableList.__init__(self, parent, value )

    def __getitem__(self):
        r"""
        Get the entry of a tree-like-tableau.
        
        EXAMPLES::
            sage: from sage.combinat.tree_like_tableaux import TreeLikeTableau
            sage: tlt = TreeLikeTableau( [[1,1,0,0],[1,0,1,1],[1,0,0]] )
            sage: tlt[0][0]
            1
            sage: tlt[0][2]
            0
            sage: tlt[2][3]
            -1
            sage: tlt[-1][0]
            -1

        TESTS::
            sage: tlt[-1][-1]
            -1
            sage: tlt[0][-1]
            -1
            sage: tlt[2][4]
            -1
            sage: tlt[2][2]
            -1
        """
        NotImplemented

    def __setitem__(self):
        r"""
        Modidy an entry of the tree-like-tableau

        EXAMPLES::
            sage: from sage.combinat.tree_like_tableaux import TreeLikeTableau
            sage: tlt = TreeLikeTableau( [[1,0,1,0],[1,1,0,1],[1,0]] )
            sage: with tlt.clone() as tlt1:
            sage:     tlt1[2][0] = 0
            sage:     tlt1[2][1] = 1
            sage: tlt1
            [[1, 0, 1, 0], [1, 1, 0, 1], [0, 1]]
            sage: with tlt.clone() as tlt2:
            sage:     tlt2[2] = [0, 1, 0]
            sage: tlt2
            [[1, 0, 1, 0], [1, 1, 0, 1], [0, 1, 0]]
        """
        self.check_mutable()
        NotImplemented

    def insert_point( self, position ):
        r"""
        This is the insertion algortihm.
        This method insert a new point inside the tree-like tableau on the 
        bordr number 'position' ((the border are numbered from left to right)

        ( see definition 2.2 of the article "tree-like tableaux", J.C. Aval, 
          A. Boussicault, P. Nadeau, arXiv:1109.0371v2 )

        EXAMPLES::
            sage: from sage.combinat.tree_like_tableaux import TreeLikeTableau
            sage: tlt = TreeLikeTableau( [[1,1,0,1],[1,0,1,0],[0,1],[0,1],[1,0]] )
            sage: with tlt.clone() as tlt1:
            sage:     tlt1.insert_point( 3 )
            sage: tlt1
            [[1, 1, 0, 0, 1], [1, 0, 0, 1, 0], [0, 1, 0, 0], [0, 1, 1, 0], [1, 0]]
        """
        self.check_mutable()
        NotImplemented

    def remove_point( self  ):
        r"""
        This is the inverse of the insertion algortihm.
        This method remove the special point from the tree-like tableau.
        This method return the position border of the removed point.

        ( see article "tree-like tableaux", J.C. Aval, A. Boussicault, 
          P. Nadeau, arXiv:1109.0371v2 )

        EXAMPLES::
            sage: from sage.combinat.tree_like_tableaux import TreeLikeTableau
            sage: tlt = TreeLikeTableau( [[1,1,0,1],[1,0,1,0],[0,1],[0,1],[1,0]] )
            sage: with tlt.clone() as tlt1:
            sage:     tlt1.insert_point( 3 )
            sage: tlt1
            [[1, 1, 0, 0, 1], [1, 0, 0, 1, 0], [0, 1, 0, 0], [0, 1, 1, 0], [1, 0]]
        """
        self.check_mutable()
        NotImplemented

    def __repr__(self):
        r"""
        The text representation of a tree-lke tableau
        
        EXAMPLES::
            sage: from sage.combinat.tree_like_tableaux import TreeLikeTableau
            sage: tlt = TreeLikeTableau( [[1,0,1,0],[1,1,0,1],[0,1],[0,1],[1,0]] )
            sage: tlt
        """
        NotImplemented

    def size( self ):
        r"""
        Return the size of the tree-like tableaux : it is the number of point 
        (the 1 inside the tableau) inside the tableau.

        EXAMPLES::
            sage: from sage.combinat.tree_like_tableaux import TreeLikeTableau
            sage: tlt = TreeLikeTableau( [[1,0,1,0],[1,1,0,1],[0,1],[0,1],[1,0]] )
            sage: tlt.size()
            8
        """
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
        raise ValueError, "Invalide argument for Tre-like tableaux Factory."

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


class TreeLikeTableaux_size(SetFactoryParent, UniqueRepresentation):
    r"""
    The tree-like tableaux of size `n`.
    """
    def __init__(self, size, policy):
        r"""
        """
        self._size = size
        SetFactoryParent.__init__(
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
            raise ValueError, "The tree-like tableau have a Wrong size : %s"%(el.size())

    def __iter__(self):
        r"""
        """
        NotImplemented



class TreeLikeTableaux_all( SetFactoryParent, DisjointUnionEnumeratedSets ):
    r"""
    This class enumerate all the tree-like tableaux.
    """
    def __init__(self, policy):
        r"""
        """
        SetFactoryParent.__init__(
            self, (), policy, category = FiniteEnumeratedSets()
        )
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), TreeLikeTableaux_size),
            facade=True, keepkey = False,
            category = self.category()
        )

    def _repr_(self):
        r"""
        """
        return "TreeLikeTableaux_all"

    def check_element(self, el, check):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: TLTS = TreeLikeTableaux()
            sage: TLTS.check_element(TLTS.an_element(), True)
        """
        pass

