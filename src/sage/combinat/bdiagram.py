r"""
The BDiagram Hopf Algebra and morphisms
=======================================

This module deals with a new combinatorial objects called BDiagrams.

AUTHORS:

- Imad Eddine Bousbaa 
- Adrien Boussicault
- Zakaria Chemli

REFERENCES:

..  Imad Eddine Bousbaa, Ali Chouria, Jean-Gabriel Luque.
   *A combinatorial Hopf algebra for the boson normal ordering problem *,
   :arxiv:`1512.05937`.
"""
#******************************************************************************
#  Copyright (C) 2016    
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.combinat.subset import Subsets
from sage.graphs.digraph import DiGraph
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.categories.sets_with_grading import SetsWithGrading
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.list_clone import ClonableList
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.set_factories import (
    SetFactory, ParentWithSetFactory, TopMostParentPolicy
)
from sage.combinat.composition import Composition
from sage.combinat.free_module import CombinatorialFreeModule
from sage.structure.parent import Parent
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.global_options import GlobalOptions
from sage.sets.set import Set
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.disjoint_union_enumerated_sets \
    import DisjointUnionEnumeratedSets
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
from sage.functions.trig import cos, sin
from sage.functions.other import sqrt
from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis


class BDiagram(ClonableList):
    r"""

    TESTS::
        sage: BDiagram( [[], [], [], []] )
        [[], {}, {}, ()]
    """
    __metaclass__ = InheritComparisonClasscallMetaclass
    @staticmethod
    def __classcall_private__(cls, value):
        r"""
        """
        return BDiagrams()(value) 
    def check(self):
        pass
    def __init__(self, parent, value):
        r"""
        """
        if not isinstance(value, (list, tuple)):
            raise ValueError(
                "Value %s must be a list or a tuple." % (value)
            )
        if not len(value) == 4:
            raise ValueError(
                "Wrong number of arguments"
            )
        ClonableList.__init__(
            self, parent, (
                Composition( value[0] ),
                Set( value[1] ),
                Set( value[2] ),
                tuple( map( tuple, value[3] ) )
            )
        )
        self.check()
            
    def size(self):
        r"""
        Return the number of half edges.

        EXAMPLES::

            sage: exem2 = BDiagram(
            ....:     ( [3,3,4,2,3,2],
            ....:     {1,2,3,4,6,5,11}, {0,1,3,5,6,7,9,10,13,14,16},
            ....:     [(1,7),(2,13),(3,10),(11,16)])
            ....: )
            sage: exem2.size()
            17
            sage: exem2 = BDiagram([[], [], [], []])
            sage: exem2.size()
            0
        """
        return self.composition().size()
        ###########################################################
        ###############  Geters of BDiagram   #####################
        ###########################################################

# set of non-cut outer  half-edges
    def outer_set(self):
        r"""
        EXAMPLES::

            sage: exem2 = BDiagram(
            ....:     ( [3,3,4,2,3,2],
            ....:     {1,2,3,4,6,5,11}, {0,1,3,5,6,7,9,10,13,14,16},
            ....:     [(1,7),(2,13),(3,10),(11,16)])
            ....: )
            sage: exem2.outer_set()
            {1, 2, 3, 4, 5, 6, 11}
            sage: exem2 = BDiagram([[], [], [], []])
            sage: exem2.outer_set()
            {}
        """
        return self[1]
   
# set of non-cut inner half-edges
    def inner_set(self):
        r"""
        EXAMPLES::

            sage: exem2 = BDiagram(
            ....:     ( [3,3,4,2,3,2],
            ....:     {1,2,3,4,6,5,11}, {0,1,3,5,6,7,9,10,13,14,16},
            ....:     [(1,7),(2,13),(3,10),(11,16)])
            ....: )
            sage: exem2.inner_set()
            {0, 1, 3, 5, 6, 7, 9, 10, 13, 14, 16}
            sage: exem2 = BDiagram([[], [], [], []])
            sage: exem2.inner_set()
            {}
        """
        return self[2]

# set of edges
    def edges_set(self):
        r"""
            sage: exem2 = BDiagram(
            ....:     ( [3,3,4,2,3,2],
            ....:     {1,2,3,4,6,5,11}, {0,1,3,5,6,7,9,10,13,14,16},
            ....:     [(1,7),(2,13),(3,10),(11,16)])
            ....: )
            sage: exem2.edges_set()
            ((1, 7), (2, 13), (3, 10), (11, 16))
            sage: exem2 = BDiagram([[], [], [], []])
            sage: exem2.edges_set()
            ()
        """
        return self[3]

# set of outer cut half-edges
    def outer_cut_set(self):
        return [e for e in range (n) if e not in self[1]]
# set of inner cut half-edges
    def inner_cut_set(self):
        return [e for e in range (n) if e not in self[2]]

    def free_outer_set(self):
        res= []
        for edge in self[1]:
            if edge not in map(lambda x: x[0],self.edges_set()):
                res.append( edge )
        return res

    def free_inner_set(self):
        res = []
        for edge in self.inner_set():
            if edge not in map(lambda x: x[1],self.edges_set()):
                res.append( edge )
        return res
    def bug_number(self):
        r"""
            sage: exem2 = BDiagram(
            ....:     ( [3,3,4,2,3,2],
            ....:     {1,2,3,4,6,5,11}, {0,1,3,5,6,7,9,10,13,14,16},
            ....:     [(1,7),(2,13),(3,10),(11,16)])
            ....: )
            sage: exem2.bug_number()
            6
            sage: exem2 = BDiagram([[], [], [], []])
            sage: exem2.bug_number()
            0
        """
        return len (self.composition())
    
    def composition(self):
        return self[0]
    
    def _edge_to_bug_number( self, half_edge):
        r"""
            sage: exem2 = BDiagram(
            ....:     ( [3,3,4,2,3,2],
            ....:     {1,2,3,4,6,5,11}, {0,1,3,5,6,7,9,10,13,14,16},
            ....:     [(1,7),(2,13),(3,10),(11,16)])
            ....: )
            sage: exem2._edge_to_bug_number(10)
            3
            sage: exem2._edge_to_bug_number(16)
            5
            sage: exem2._edge_to_bug_number(0)
            0
            sage: exem2._edge_to_bug_number(1)
            0
            sage: exem2._edge_to_bug_number(2)
            0
            sage: exem2._edge_to_bug_number(3)
            1
            sage: exem2._edge_to_bug_number(4)
            1
            sage: exem2._edge_to_bug_number(5)
            1
            sage: exem2._edge_to_bug_number(6)
            2
            sage: exem2 = BDiagram([[],[],[],[]])
            sage: exem2._edge_to_bug_number(0)
            Traceback (most recent call last):
            ...
            AssertionError
        """
        assert( self.bug_number() > 0 ) 
        i=0; compo=self.composition()[i] 
        while half_edge>=compo:
            i+=1            
            compo=compo+self.composition()[i]
        return i


        ###########################################################
        ###############  Sub Bdiagram  ############################
        ###########################################################

    def _bug_digraph(self):
        r"""
            sage: exem2 = BDiagram(
            ....:     ( [3,3,4,2,3,2],
            ....:     {1,2,3,4,6,5,11}, {0,1,3,5,6,7,9,10,13,14,16},
            ....:     [(1,7),(2,13),(3,10),(11,16)])
            ....: )
            sage: d=exem2._bug_digraph()
            sage: d.vertices()
            [0, 1, 2, 3, 4, 5]
            sage: d.edges()
            [(0, 2, None), (0, 4, None), (1, 3, None), (3, 5, None)]
        """
        liste_adj = {}
        for vertex in range( self.bug_number() ):
            liste_adj[vertex] = []
        for edge in self.edges_set():
            liste_adj[
                self._edge_to_bug_number(edge[0])
            ].append( self._edge_to_bug_number(edge[1]) )
        return DiGraph( liste_adj )

    def _sub_intervals(self, list_of_bug):
        r"""
        Return a list of interval of half edges id
        associated with each bug of ``list_of_bug``.

        The ``list_of_bug`` have to be sorted !

        EXAMPLES::

            sage: bd = BDiagram([[3,7,5,11,13,2], [], [], []])
            sage: bd._sub_intervals([2,4,5])
            [(10, 15), (26, 39), (39, 41)]
            sage: bd._sub_intervals([])
            []
            sage: bd = BDiagram([[], [], [], []])
            sage: bd._sub_intervals([])
            []
        """
        res = []
        prev=0
        mi=0
        for i in list_of_bug:
            for j in range(prev, i):
                mi+= self.composition()[j]
            prev = i
            res.append( (mi,mi+self.composition()[i]) )
        return res

    def sub_bdiagram( self, list_of_bug ):
        r"""
        EXAMPLES::
            sage: exem2 = BDiagram(
            ....:     ( [3,3,4,2,3,2],
            ....:     {1,2,3,4,5,6,11}, {0,1,3,5,6,7,9,10,13,14,16},
            ....:     [(1,7),(2,13),(3,10),(11,16)])
            ....: )
            sage: exem2.sub_bdiagram([1, 3, 5])
            [[3, 2, 2], {0, 1, 2, 4}, {0, 2, 3, 6}, ((0, 3), (4, 6))]
            sage: exem2.sub_bdiagram([0, 2, 4])
            [[3, 4, 3], {1, 2, 3}, {0, 1, 3, 4, 6, 8, 9}, ((1, 4), (2, 8))]
            sage: exem2.sub_bdiagram([0, 1, 2, 3, 4, 5])
            [[3, 3, 4, 2, 3, 2], {1, 2, 3, 4, 5, 6, 11}, {0, 1, 3, 5, 6, 7, 9, 10, 13, 14, 16}, ((1, 7), (2, 13), (3, 10), (11, 16))]
            sage: exem2.sub_bdiagram([])
            [[], {}, {}, ()]
            sage: exem2 = BDiagram([[], [], [], []])
            sage: exem2.sub_bdiagram([])
            [[], {}, {}, ()]
            
        """
        sub_edges = []
        sub_compo = []
        sub_outers= []
        sub_inners= []
        dictio_of_outers = dict()
        dictio_of_inners = dict()
        # get intervals
        list_of_bug.sort()
        intervals = self._sub_intervals(list_of_bug)
        # get sub composition
        for i in list_of_bug:
            sub_compo.append( self.composition()[i] )
        # get sub outers set and inner outer set
        index = 0
        lenght = 0

        for interval in intervals:
            for i in self.outer_set():
                if interval[0] <= i and i< interval[1]:
                    dictio_of_outers[i] = i-interval[0]+ lenght
                    sub_outers.append( dictio_of_outers[i] )
            for i in self.inner_set():
                if interval[0] <= i and i< interval[1]:
                    dictio_of_inners [i] = i-interval[0] + lenght 
                    sub_inners.append( dictio_of_inners[i] )
            lenght+=sub_compo[index]
            index+=1
        # get sub edges
        for edge in self.edges_set():
            if self._edge_to_bug_number(edge[0]) in list_of_bug:
                sub_edges.append( ( dictio_of_outers[edge[0]], dictio_of_inners[edge[1]] ) )
        # return the sub diagram 
        return BDiagram ( [sub_compo,sub_outers,sub_inners,sub_edges]  )


###############################################################################################################
################################ The Factory of BDiagrams  ####################################################
###############################################################################################################

class BDiagramsFactory():
    def __call__(self, size=None):
        if size is None:
            return BDiagrams_all()
        else:
            assert isinstance(size,(Integer, int)) and size >= 0, '%s is not a non-negative integer.' %size
            return BDiagrams_size(size)

class BDiagrams_all(DisjointUnionEnumeratedSets):
    Element = BDiagram
    def __init__(self):
        DisjointUnionEnumeratedSets.__init__(
            self,            
            Family(NonNegativeIntegers(), self.graded_component),
            facade=True, keepkey=False
        )

    def __repr__(self):
        r"""
        retourne le nom de la classe parent.
        EXAMPLES::

        """
        return "B Diagrams"

    def graded_component(self, taille):
        return BDiagrams_size(taille)

class BDiagrams_size(UniqueRepresentation, Parent):
    Element = BDiagram
    def __init__(self, size):
        assert isinstance(size, (int, Integer)) and size >= 0
        self._size = size
        Parent.__init__(self, category = FiniteEnumeratedSets())

    def _repr_(self):
        r"""
        EXAMPLES::
        """
        return 'BDiagrams of size %s' %(self._size)

    def __contains__(self, obj):
        pass

    def __iter__(self):
        r"""
        EXAMPLES::
        """
        size=self._size 
        for comp in Compositions(size):
            for c in cartesian_product( [Subsets(size), Subsets(size)] ) :
                for edges in self._allowed_matchings( c[0], c[1], comp ):
                    yield BDiagram( (comp, c[0], c[1], tuple(edges)) )


    def _allowed_matchings2_rec( self, pos_inner_edge, inner_set, outer_set, comp ):
        if len(inner_set) <= pos_inner_edge or len(outer_set) == 0:
            yield []
            return
        for edges in self._allowed_matchings2_rec( pos_inner_edge+1, inner_set, outer_set, comp ):
            yield edges
        for i in range( len(outer_set) ):
            edge = ( inner_set[pos_inner_edge], outer_set[i] )
            if self._edge_is_allowed( comp, edge ):
                for edges in self._allowed_matchings2_rec( pos_inner_edge+1, inner_set, outer_set[:i] + outer_set[i+1:], comp ):
                    yield edges + [edge]


    def _allowed_matchings( self, inner_set, outer_set, comp ):
        if len(outer_set) == 0:
            yield []
            return    
        for m in self._allowed_matchings2_rec( 0, inner_set, outer_set, comp ):
            yield m

    def _get_bug_position( self, comp, half_edge):
        i=0; compo=comp[i] 
        while half_edge>compo:
            i+=1            
            compo=compo+comp[i]
        return i

    def _edge_is_allowed(self, comp, edge):
        if edge[0]>=edge[1]:
            return False
        if self._get_bug_position( comp,edge[0]) < self._get_bug_position( comp,edge[1]):
            return True 


BDiagrams = BDiagramsFactory()

#
#
#def diagram_size( d ):
#    return d[0].size()
#
#def diagram_edges( d ):
#    return d[3]
#
#def diagram_outer_half_edges( d ):
#    return d[1]
# 
#def diagram_inner_half_edges( d ):
#    return d[2]
#
#def nb_outer_free_half_edges( d ):
#    return len( diagram_outer_half_edges( d ) )-len( diagram_edges( d ) ) 
#
#var('x,y')
#
#def poids( d ):
#    return x**( diagram_size( d ) ) * y**( nb_outer_free_half_edges( d ) )
#
#


###############################################################################################################
################################ BDiagrams Hopf Algebra  ######################################################
###############################################################################################################


class BDiagramHopfAlgebra(CombinatorialFreeModule):
    def __init__(self, base_ring):
        CombinatorialFreeModule.__init__(self, base_ring, BDiagrams(),category=GradedHopfAlgebrasWithBasis(base_ring).Connected())

    @cached_method
    def one_basis(self):
        r"""
            sage: hbd = BDiagramHopfAlgebra(QQ)
            sage: hbd.one_basis()
            [[], {}, {}, ()]
            sage: hbd.one()
            B[[], {}, {}, ()]
            sage: isinstance(hbd.one_basis(), BDiagram)
            True
        """
        return self.basis().keys()( ([], [], [], []) )

    def degree_on_basis(self, bd):
        return bd.size()

    def _repr_(self):
        return "BDiagram Hopf algebrawith basis over %s" % self.base_ring()

    def _repr_term(self, bd):
        return 'B' + repr(bd)

        ###########################################################
        ###############  Product of Bdiagrams #####################
        ###########################################################


    def _shift_list(self,l,sh):
        return map(lambda x: x+sh, l )
    
    def _matchings_rec( self, pos_inner_edge, inner_set, outer_set):
        if len(inner_set) <= pos_inner_edge or len(outer_set) == 0:
            yield []
            return
        for edges in self._matchings_rec( pos_inner_edge+1, inner_set, outer_set ):
            yield edges
        for i in range( len(outer_set) ):
            edge = ( inner_set[pos_inner_edge], outer_set[i] )
            for edges in self._matchings_rec( pos_inner_edge+1, inner_set, outer_set[:i] + outer_set[i+1:] ):
                yield edges + [edge]


    def _matchings( self, inner_set, outer_set):
        if len(outer_set) == 0 or len(inner_set) == 0:
            yield []
            return    
        for m in self._matchings_rec( 0, inner_set, outer_set):
            yield m

    def product_on_basis(self, bdiag_1, bdiag_2):
        res = 0
        comp = bdiag_1[0]+bdiag_2[0]
        outer_set =list(bdiag_1.outer_set())+list(self._shift_list(bdiag_2.outer_set(),bdiag_1.size()))
        inner_set = list(bdiag_1.inner_set())+list(self._shift_list(bdiag_2.inner_set(),bdiag_1.size()))
        edges = list(bdiag_1.edges_set())+list(self._shift_list(bdiag_2.edges_set(),bdiag_1.size()))

        for match in self._matchings(bdiag_1.free_outer_set(),self._shift_list(bdiag_2.free_inner_set(), bdiag_1.size() )):
            res += self( BDiagram( (comp,outer_set,inner_set,edges+match) ) )
        return res

        ###########################################################
        ###############  Coproduct of Bdiagram  ###################
        ###########################################################

    def coproduct_on_basis(self, bdiagram):
        """
        The coproduct of bdiagram.

        .. MATH::

            \Delta(P_i) = \sum_{j=0}^i P_{i-j} \otimes P_j

        INPUT:

        OUTPUT:

        TESTS::

            sage: hbd = BDiagramHopfAlgebra(QQ)
            sage: b = BDiagram(
            ....:     ( [3,3,4,2,3,2],
            ....:     {1,2,3,4,5,6,11}, {0,1,3,5,6,7,9,10,13,14,16},
            ....:     [(1,7),(2,13),(3,10),(11,16)])
            ....: )
            sage: hbd.coproduct_on_basis( b )
            B[[], {}, {}, ()] # B[[3, 3, 4, 2, 3, 2], {1, 2, 3, 4, 5, 6, 11}, {0, 1, 3, 5, 6, 7, 9, 10, 13, 14, 16}, ((1, 7), (2, 13), (3, 10), (11, 16))] + B[[3, 2, 2], {0, 1, 2, 4}, {0, 2, 3, 6}, ((0, 3), (4, 6))] # B[[3, 4, 3], {1, 2, 3}, {0, 1, 3, 4, 6, 8, 9}, ((1, 4), (2, 8))] + B[[3, 3, 4, 2, 3, 2], {1, 2, 3, 4, 5, 6, 11}, {0, 1, 3, 5, 6, 7, 9, 10, 13, 14, 16}, ((1, 7), (2, 13), (3, 10), (11, 16))] # B[[], {}, {}, ()] + B[[3, 4, 3], {1, 2, 3}, {0, 1, 3, 4, 6, 8, 9}, ((1, 4), (2, 8))] # B[[3, 2, 2], {0, 1, 2, 4}, {0, 2, 3, 6}, ((0, 3), (4, 6))] 
        """
        res = 0
        graph = bdiagram._bug_digraph()
        connected_components = graph.connected_components()
        nb_compo = len(connected_components)
        id_compo = range(nb_compo)
        for left_compo in Subsets( id_compo ):
            right_compo = Set(id_compo) - left_compo
            left = []
            for i in left_compo:
                left += connected_components[i]
            right = []
            for i in right_compo:
                right += connected_components[i]
            left = bdiagram.sub_bdiagram( left )
            right = bdiagram.sub_bdiagram( right ) 
            res += self( left ).tensor( self(right) )
        return res
