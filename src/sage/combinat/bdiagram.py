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


class BDiagram(ClonableList):
    r"""
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
        """
        return len(self[0]) 
################################################################################
######################  Product of Elements  ###################################
################################################################################
    def outer_set(self):
        return self[1]
    
    def inner_set(self):
        return self[2]
    def outer_cut_set(self):
        return [e for e in range (n) if e not in self[1]]
    def inner_cut_set(self):
        return [e for e in range (n) if e not in self[2]]
    def free_outer_set(self):
        res= []
        for edge in self[1]:
            if edge not in map(lambda x: x[0],self.edges_set()):
                res.append( edge )
        return res
    def edges_set(self):
        return self[3]
    def free_inner_set(self):
        res = []
        for edge in self.inner_set():
            if edge not in map(lambda x: x[1],self.edges_set()):
                res.append( edge )
        return res
    def _bugs_number(self):
        return len (self[0])

    def _edge_bug_number( self, half_edge):
        i=0; compo=self[0][i] 
        while half_edge>compo:
            i+=1            
            compo=compo+self[0][i]
        return i

    def _matrix_connection(self):
        Mat = Matrix( self._bugs_number() )
        for edge in self.edges_set():
            print ( (self._edge_bug_number(edge[0]),self._edge_bug_number(edge[1])) )
            Mat[self._edge_bug_number(edge[0]),self._edge_bug_number(edge[1])]=1
        return Mat


#
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
#
#
#
#
#
#
#
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
#
#def get_position_first_allowed_inner(set_inner,comp,edge):
#    i=0; compo=comp[i]; liste_valid_inner=[]
#    while edge>compo:
#        i+=1
#        compo=compo+comp[i]
#    for pos in range(len(set_inner)):
#        if set_inner[pos]>=compo:
#            return pos
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
##class BDiagramsFactory():
##    def __call__(self):
##        pass
##class BDiagrams_all():
##    pass 
#
##class BDiagrams_size():
##    def __init__(self.size):
##        pass
##    def __iter__(self):
##        pass
#
##class BDiagram():
##    pass
#
##BDiagrams=BdiagramsFactory
#
#
#
#
##def get_free_outer_set(outer_set,connection_set):
##        set_free_outer=set()
##    for edge in outer_set:
##        if edge not in map(lambda x: x[0],connection_set):
##            set_free_outer.add(edge)
##    return set_free_outer
#
##def get_free_inner_set(inner_set,connection_set):
##        set_free_inner=set()
##    for edge in inner_set:
##        if edge not in map(lambda x: x[1],connection_set):
##            set_free_inner.add(edge)
##    return set_free_inner
#
#


class BDiagramHopfAlgebra(CombinatorialFreeModule):
    def __init__(self, base_ring):
        CombinatorialFreeModule.__init__(self, base_ring, BDiagrams(),category=GradedHopfAlgebrasWithBasis(base_ring).Connected())

    @cached_method
    def one_basis(self):
        return self.basis().keys()( ([], [], [], []) )

    def degree_on_basis(self, bd):
        return bd.size()

    def _repr_(self):
        return "BDiagram Hopf algebrawith basis over %s" % self.base_ring()

    def _repr_term(self, bd):
        return 'B' + repr(bd)

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

    def coproduct_on_basis(self, i):
        """
        The coproduct of bdiagram.

        .. MATH::

            \Delta(P_i) = \sum_{j=0}^i P_{i-j} \otimes P_j

        INPUT:

        - ``i`` -- a non-negative integer

        OUTPUT:

        - an element of the tensor square of ``self``

        TESTS::

            sage: H = GradedHopfAlgebrasWithBasis(QQ).Connected().example()
            sage: H.monomial(3).coproduct()
            P0 # P3 + 3*P1 # P2 + 3*P2 # P1 + P3 # P0

        """
        pass
#        return self.sum_of_terms(
#            ((i-j, j), binomial(i, j))
#            for j in range(i+1)
#        )






