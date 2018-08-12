# distutils: language = c++
"""
Weisfeiler Lehman Algorithm

This method implements the Weisfeiler Lehman algorithm itself and a correctness checker for its results.
Furthermore it offers the internal method used to compute the orbitals for the correctness checker.
"""


from __future__ import print_function
from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string
from libcpp.utility cimport pair
from libcpp cimport bool
from sys import version_info
from copy import copy 
import random
from sys import version_info    
from sage.graphs.graph import GenericGraph as SageGraph
from sage.misc.latex import latex, dict_function

def _getiterator(l):
    if version_info >= (3):
        return l.items()
    else:
        return l.iteritems()

def _edge_multiplicity_to_label(G, edge_labels=False):
    from collections import Counter
    edge_multiplicity_dict = Counter(G.multiple_edges(labels=False))
   
    G.remove_multiple_edges()
    G.allow_multiple_edges(False)
    for u, v, m in G.edges():
        if (u,v) not in edge_multiplicity_dict:
            l = 1
        else:
            l = edge_multiplicity_dict[(u,v)]
        if edge_labels:
            new_label = (m, l)
        else:
            new_label = l
        G.set_edge_label(u, v, new_label)


class _Graph(object):
    '''
    For each multiple edge e whose elements have labels l_0,l_1,...,l_(m-1), 
    delete the multiple edge and substitute it with a single edge, with label 
    (l_i, m) where i can be any integer between 0 and m-1
    '''
            
    
    
    '''
    Graph instantiates an adjacency dictionary based graph object.
    It can represent vertex colored, directed or undirected graphs.
    '''
    def __init__(self, graph,
                 edge_labels=False):
        '''
        *graph*
            Sage Graph object.  Mandatory argument.

        *edge_labels*
            Boolean indicating if the method must keep edge labels into account
        '''
        G = copy(graph)
        multiedged = G.has_multiple_edges()
        self._first_level_vertices = {}
        self._vertex_coloring = []
        self.graph = G
        if multiedged:
            _edge_multiplicity_to_label(self.graph, edge_labels)
            edge_labels = True
        self.number_of_vertices = self.graph.order()
        self.directed = self.graph.is_directed()
        self._relabel_map = {}
        self._normalized = False

    def relabel(self):
        if self._normalized:
            return
        self._normalized = True
        m = self.graph.relabel(complete_partial_function=False, return_map=True)
        self._inverted_relabel_map = {v: k for k, v in _getiterator(m)}
        return m

    def _get_vertex_coloring(self):
        return self._vertex_coloring

    vertex_coloring = property(_get_vertex_coloring)

    def set_vertex_coloring(self, vertex_coloring, relabel_map={}):
        '''
        Define a vertex coloring of the Graph.

        *vertex_coloring*
            A list of disjoint sets of vertices representing a
            partition of the vertex set; vertices not listed are
            placed into a single additional part.
        '''
        self._vertex_coloring = []
        if vertex_coloring == None:
            vertex_coloring = []
        vertex_coloring = [set(v) for v in vertex_coloring]
        if vertex_coloring or self._first_level_vertices:
            
            if self._first_level_vertices:
                vc = [[int(self._first_level_vertices[v]) for v in k] for k in vertex_coloring]
            else:
                vc = vertex_coloring
            
            if relabel_map:
                vc = [set(int(relabel_map[v]) for v in k) for k in vc]
            
            if self._first_level_vertices:
                vs = set([int(relabel_map[v]) for k,v in _getiterator(self._first_level_vertices)])
                self._vertex_coloring.append(set([int(k) for k in range(self.number_of_vertices) if not k in vs]))
            else:
                vs = set(range(self.number_of_vertices))
                
            for p in vc:
                if p <= vs:
                    self._vertex_coloring.append([int(v) for v in p])
                    vs -= p
                else:
                    raise ValueError('Invalid partition: %s' % vertex_coloring)
            if vs:
                self._vertex_coloring.append(vs)
                
            if len(self._vertex_coloring) == 1 and not self._first_level_vertices:
                self._vertex_coloring = []

    def _repr_(self):
        
        s = ['Graph(sage_graph=%s,' %
             (self.graph._repr_)]
        s.append(' vertex_coloring = [')
        for x in self._vertex_coloring:
            s.append('  set(%s),' % list(x))
        s.append(' ],')
        s.append(')')
        return '\n'.join(s)
    
    def _latex_(self):
    
        s = latex(self.graph)
        s.append("\n");
        s.append("Vertex coloring: %s" % (dict_function(s.vertex_coloring)))
        return '\n'.joins(s)
    
    def undo_relabel(self):
        '''
        Beware, undoing the label doesn't fix the vertex coloring accordingly, at this stage
        '''        
        if not self._normalized:
            return
        self.graph.relabel(complete_partial_function=False, perm=self.inverted_relabel_map)
        del self._inverted_relabel_map
        self._normalized = False


cdef vector[GraphNode] _sageGraphToLists(G, partition = [], has_edge_labels=False):
    partition_dict = {v:i for i,x in enumerate(partition,1) for v in x}
    
    cdef vector[GraphNode] nodeArray;
    nodeArray.resize(G.order())
    
    for v in G: #Graph is relabeled already
        nodeArray[v].idx = v
        nodeArray[v].color = partition_dict.get(v, 0)
        nodeArray[v].adj_list.reserve(len(G[v]))
    directed = G.is_directed()
    
    edge_labels = G.edge_labels()
    edge_labels_dict = {}
    counter = 1
    for el in edge_labels:
        if(has_edge_labels):
            if(edge_labels_dict.setdefault(el, counter) == counter):
                counter += 1
        else:
            edge_labels_dict[el] = counter
    
    for v,u,l in G.edges(labels=True):
        nodeArray[v].adj_list.push_back([<int>u, <int>edge_labels_dict[l]])
        if not directed:
            nodeArray[u].adj_list.push_back([<int>v, <int>edge_labels_dict[l]])
    return nodeArray

def _sameOrbital(a, b, p_to_o, v_to_so):
    u,v = a
    x,y = b
    for perm, orbit_dict in p_to_o.items():
        if u == orbit_dict[x]:
            g = perm
            break
    else:
        return False
    p = g(y)
    return v_to_so[v][u] == v_to_so[p][u]

def compute_orbitals(G, partition, edge_labels):
    """
    Compute the orbitals of the graph G, as a list of lists

    INPUT:

    -  ``G`` - Graph of which to compute the orbitals

    - ``partition`` - Default is the unit partition, otherwise computes the orbitals of the subgroup of the full automorphism group respecting the partition.
    
    - ``edge_labels`` - Default False, otherwise computes the orbitals of the subgroup of the full automorphism groups whose permutations respect edge labels

    OUTPUT:
    
    A list of lists representing the partition induced by the orbitals on `V^2`, where `V` is the vertex set of ``G``
    """
    import itertools
    aut = G.automorphism_group(edge_labels=edge_labels, partition=partition)
    orbits = aut.orbits()
    perms_to_orbits = {}
    vertex_to_stabilizer_orbit = {}
    for perm in aut:
        d = perms_to_orbits.setdefault(perm, {})
        for v in G:
            d[v] = perm(v)
    orbitals = {}
    for v,u in itertools.product(G,G):
        for rep in orbitals:
            if _sameOrbital(rep, (v,u), perms_to_orbits, vertex_to_stabilizer_orbit):
                orbitals[rep].append((v,u))
                break
        else:
            orbitals[(v,u)] = [(v,u)]
            _addStabilizer(aut, v, vertex_to_stabilizer_orbit)
    return orbitals.values()

def _addStabilizer(aut, v, vertex_to_stabilizer_orbit):
    stabO = aut.stabilizer(v).orbits()
    for i, o in enumerate(stabO):
        for u in o:
            stabilizer_orbit = vertex_to_stabilizer_orbit.setdefault(u, {})
            stabilizer_orbit[v] = i
    
def check_orbit_correctness(color_classes, g, partition = None, edge_labels=False, cardinality=2):
    """
    Check the results of the ``WeisfeilerLehman`` method output against the color_classes of the graph ``g``

    The algorithm computes either the orbits the graph ``g``, if the ``color_classes`` are on vertices,
    or its orbitals otherwise (through the ``compute_orbitals`` method), and compares them against the ``color_classes`` returned by the
    ``WeisfeilerLehman`` method.

    INPUT:

    -  ``color_classes`` - List of lists representing a partition of either the vertices or the edges of g, output by the ``WeisfeilerLehman`` method

    -  ``g`` - Graph on which the ``WeisfeilerLehman`` method was run

    - ``partition`` - Default is the unit partition, otherwise checks ``color_classes`` against the orbits (or orbitals) of the subgroup of the full automorphism group respecting the partition.
    
    - ``edge_labels`` - Default False, otherwise checks ``color_classes`` only against the orbits (or orbitals) of automorphism subgroups whose permutations respect edge labels

    -  ``cardinality`` - 1 if ``color_classes`` is a partition on vertices, 2 if it's on edges, default is 2

    OUTPUT:
    
    It outputs one of three strings:
    
    - "Correct" - if ``color_classes`` represents the orbits of  ``g``
    - "Refinable" - if ``color_classes`` represents a (partial) union of the orbits of  ``g``, and thus a higher value of k passed to the ``WeisfeilerLehman`` method might output the correct result
    - "Wrong" - otherwise
    """
    if(cardinality != 1 and cardinality != 2):
        raise ValueError("Cardinality must be either 1 or 2")
    G = g
    if not partition:
        partition = [G.vertices()]
    n = G.order()
    checker = sorted(map(sorted,color_classes))
    set_checker = map(set, checker)
    #Check the color_classes are well formatted
    union_checker = set.union(*set_checker)
    if len(union_checker) != sum(map(len,color_classes)):
        raise ValueError("The color_classes are not disjoint")
    if cardinality == 2 and len(union_checker) != n*n:
        if len(union_checker) == (n*(n+1))/2:
            raise ValueError("The color_classes don't contain both directions of every edge")
        else:
            raise ValueError("The color_classes are not complete")
    elif cardinality == 1 and len(union_checker) != n:
        raise ValueError("The color_classes are not complete")
    #Check if they are correct
    if cardinality == 2:
        O = compute_orbitals(G, partition, edge_labels)
        O.sort()
    elif cardinality == 1:
        O = G.automorphism_group(partition=partition, edge_labels=edge_labels, orbits=True, return_group=False)
        O = map(sorted, O)
        O.sort()
    if O == checker: return "Correct"
    #Check if they only need to be refined
    O = map(set,O)
    elements_left = set(range(len(O)))
    for s in set_checker:
        temp = set()
        for el in elements_left:
            if O[el].issubset(s):
                temp.add(el)
        elements_left = elements_left.difference(temp)
        if len(elements_left) == 0:
            break
    if len(elements_left) == 0:
        return "Refinable"
    else:
        return "Wrong"



def _check_orbit_correctness2(color_classes, g, cardinality=2):
    """
    Check the results of the ``WeisfeilerLehman`` method output against the color_classes of the graph ``g``

    The algorithm computes the orbits of either the graph ``g``, if the ``color_classes`` are on vertices,
    or of a suitable line graph otherwise, and compares them against the ``color_classes`` returned by the
    ``WeisfeilerLehman`` method.
    
    .. WARNING::
    
        This implementation, while at times faster than ``check_orbit_correctness`` method's one, consumes a lot more memory
    
    INPUT:

    -  ``color_classes`` - List of lists representing a partition of either the vertices or the edges of g, output by the ``WeisfeilerLehman`` method

    -  ``g`` - Graph on which the ``WeisfeilerLehman`` method was run

    -  ``cardinality`` - 1 if ``color_classes`` is a partition on vertices, 2 if it's on edges

    OUTPUT:
    
    It outputs one of three strings:
    
    - "Correct" - if ``color_classes`` represents the orbits of  ``g``
    - "Refinable" - if ``color_classes`` represents a (partial) union of the orbits of  ``g``, and thus a higher value of k passed to the ``WeisfeilerLehman`` method might output the correct result
    - "Wrong" - otherwise
    """
    if(cardinality != 1 and cardinality != 2):
        raise ValueError("Cardinality must be either 1 or 2")
    h = g.copy()
    if h.has_multiple_edges():
        _edge_multiplicity_to_label(h, True)
    G = h.copy()
    n = G.order()
    checker = sorted(map(sorted,color_classes))
    set_checker = map(set, checker)
    #Check the color_classes are well formatted
    union_checker = set.union(*set_checker)
    if len(union_checker) != sum(map(len,color_classes)):
        raise ValueError("The color_classes are not disjoint")
    if cardinality == 2 and len(union_checker) != n*n:
        if len(union_checker) == (n*(n+1))/2:
            raise ValueError("The color_classes don't contain both directions of every edge")
        else:
            raise ValueError("The color_classes are not complete")
    elif cardinality == 1 and len(union_checker) != n:
        raise ValueError("The color_classes are not complete")
    #Check if they are correct
    if cardinality == 2:
        import sage.graphs.line_graph
        from collections import Counter
        notG = G.complement()
        newlabel = 0
        labels = G.edge_labels()
        while newlabel in labels:
            newlabel += 1
        for (v,u) in notG.edges(labels=False):
            G.add_edge(v,u,newlabel)
        lineG = G.to_directed().line_graph()
        oG = h.automorphism_group(edge_labels=True, orbits=True, return_group=False)
        oG = map(sorted, oG)
        O = []
        for o in oG:
            O.append([])
            for el in o:
                O[-1].append((el,el))
        part_dict = {}
        for v in lineG:
            part_dict.setdefault(v[2], []).append(v)
        o_lineG = lineG.automorphism_group(partition=part_dict.values(), orbits=True, return_group=False)
        o_lineG = map(lambda l: [(v[0],v[1]) for v in l], o_lineG)
        O = O + map(sorted, o_lineG)
        O.sort()
    elif cardinality == 1:
        O = g.automorphism_group(orbits=True, return_group=False)
        O.sort()
    if O == checker: return "Correct"
    #Check if they only need to be refined
    O = map(set,O)
    elements_left = set(range(len(O)))
    subsets = []
    for s in set_checker:
        subsets.append([])
        temp = set()
        for el in elements_left:
            if O[el].issubset(s):
                temp.add(el)
                subsets[-1].append(el)
        elements_left = elements_left.difference(temp)
        if len(elements_left) == 0:
            break
    if len(elements_left) == 0:
        return "Refinable"
    else:
        return "Wrong"

def WeisfeilerLehman(G, k, partition=[], edge_labels=False, result='edge_classes'):
    """
    Return the coloring found by applying the k-th order of the Weisfeiler Lehman algorithm.

    The Weisfeiler Lehman method is a way of generating a partition of the nodes (or edges) of a graph ``G``
    such that the orbits of the automorphism group on G's vertices (respectively edges) refine said partition.
    This method was first described in [LW1968]_ (in russian, link to translation in the reference)
    and is also better and succintly described in [GKMS2017]_.
    The idea can be summarised as doing multiple rounds of vertex (resp. edge) coloring, where each round's color
    is determined by the colors of the (hyper)edges in the previous round, until the coloring is stable,
    that is the color classes stay the same between two subsequent rounds of coloring.
    Higher values of ``k`` (or higher orders of the WL method) base the choice of color for each round
    on interactions between a higher number of vertices, and thus correctly return the orbits of strictly more
    graphs than lower values of ``k``.
    
    This method can be used as a negative oracle for isomorphism between two graphs ``G`` and ``H`` by running it
    on their disjoint union: if the two components have different coloring, the graphs are surely not isomorphic; nothing
    can be said if their coloring is the same instead, since WL doesn't give any guarantee of returning the correct orbits
    instead of just a set of sets of vertices that is by them refined.
    
    A particular case is that of planar graphs, which are always distinguished by k-WL for `k \geq 3`, as proved in [KPS2017]_

    INPUT:

    -  ``G`` - Graph to be colored

    -  ``k`` - Order of the Weisfeiler Lehman algorithm to be used

    -  ``partition`` - A list of lists representing the partition of the vertices induced by the initial coloring of the vertices of ``G``

    -  ``edge_labels`` - If True, take into account the labeling of the edges of `G``
    
    -  ``result`` - Value that allows to format the output differently. Currently accepted values are 'vertex_classes' and 'edge_classes', which return
                    the color classes on, respectively, vertices and edges as a list of lists, and 'graph', which instead returns a coloured looped weighted complete graph ``C``
                    with the same vertex set as ``G``, where each edge is coloured according to its computed color class.
                    The color classes on vertices are represented by colouring the self loops on each vertex of ``C``

    OUTPUT:
    
    A list of lists of vertices (resp. edges) representing how the vertices (resp. edges) have been divided in color classes, 
    or a coloured looped weighted complete digraph, depending on the value of the ``result`` parameter
    
    .. WARNING::
        
        The current version of the method, while supporting multiple edges and labels on both vertices and edges, does NOT support self-loops.
        If the graph to be colored contains loops, either remove the loops and color the corresponding vertices differently, or create a support
        graph G' where each vertex with a self loop is transformed into two vertices that are connected to the same vertices as the original one, 
        and that have an edge (or two opposing edges if it's a DiGraph) between them

    EXAMPLES:
    
    The method can't distinguish the CycleGraph with 7 vertices and the disjoint union of two CycleGraphs with 3 and 4 vertices for ``k`` = 1, but can for ``k`` = 2.  ::
    
        sage: import sage.graphs.weisfeiler_lehman
        sage: G = graphs.CycleGraph(7)
        sage: H = graphs.CycleGraph(4)
        sage: H = H.disjoint_union(graphs.CycleGraph(3))
        sage: H.relabel()
        sage: res1 = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(G, 1, result='vertex_classes')
        sage: res2 = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(H, 1, result='vertex_classes')
        sage: res1 == res2
        True
        sage: res1 = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(G, 2, result='vertex_classes')
        sage: res2 = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(H, 2, result='vertex_classes')
        sage: res1 == res2
        False
        sage: len(res1)
        1
        sage: len(res2)
        2
    
    The orbits on vertices for the Shrikhande are immediately found for ``k`` = 1 and stay stable for increasing values of ``k``.  ::
        
        sage: import sage.graphs.weisfeiler_lehman
        sage: g = graphs.ShrikhandeGraph()
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(g, 1, result='vertex_classes')
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, g, cardinality=1)
        'Correct'
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(g, 2, result='vertex_classes')
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, g, cardinality=1)
        'Correct'
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(g, 3, result='vertex_classes')
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, g, cardinality=1)
        'Correct'
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(g, 4, result='vertex_classes')
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, g, cardinality=1)
        'Correct'

    The color classes on edges for the Shrikhande aren't equal to the orbits until ``k`` = 3, but stay stable for increasing values of ``k``.  ::

        sage: import sage.graphs.weisfeiler_lehman
        sage: g = graphs.ShrikhandeGraph()
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(g, 1)
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, g)
        'Refinable'
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(g, 2)
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, g)
        'Refinable'
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(g, 3)
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, g)
        'Correct'
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(g, 4)
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, g)
        'Correct'
        
    Check that the counter-example in :trac:`10482` is correctly recognised. ::
        
        sage: import sage.graphs.weisfeiler_lehman
        sage: G = Graph(':~?@?_??C????_W???@?__S?????OGCCA@?_WKF???OGGCA@_oWOH??_OOKECA@OgU?@?__WKGCA`OoY?????GCA@@?_OKEBA@?_OIDA`_oWMFBo?GGCB@`?_OIDA`_oWKFB`p?a???_OGGCB@_o_OGDA`OgWKEB`ow_OHCaW?CCA@_o_OGDA`OoWKFB`ow_OHCa`S????_OGGCA@_oWKGCA`OgWKFB`p?_QHDA`OkUJ`?_OKEBA@?gSKEB`ow_QIDA`WkUJEBGGCAB@_o_SIDB@ow[MGCAPGgSJDap_oYL`_oWOIDB@ow[MGCa`OkUJDb@_oYLEb`oy????OGCCA@?oWKGDB@p?_QHDA`Wm?@@?o_SKFCA@GcSIDap_oYLFB`w{_OGCG???@?_OOGCB@_o_SKFCA@GcSIDaq@?_PGcQK?AA@`?gWMGCAPGgSJDb@_sYMFBpx?aPHO???A@?__OGEB@`?gWMGCAPGgSJDcA@?aPGcaPKf??__WOIEBa@?cQIDApWoWLEb`o{]OGCQHGcQHcqXOh????GCAA@?_WKECA`_w_OHCa`OkUOGCAHCaQHCqXOgTIdW?CCBA@Oo[OGCaPOgUJEB@gs[MFbq@CaQHcq`OiTJDg??A@@?_WKGCA`OoWMFCAPOkWLFBq@?_OGcQHCcQHCqXKgSIDAhSiUJDax[nCA@OgWKFBa@GgUKEB@gsYMFB`w{]OGCAHCcQHCqXOgSIdQhWkUJDqx[oW_??GCCA@_o_OIDB@_w[OHDAp_s[NGCAHCaQHCqXKgSIDQhSiUJDax[mWKERHeGCA`OoWMFCAPOkWKEBPgs[MFBpw{aRIDA`SiTJDap[mVJeB@cqYLO??CAA@?oWOGDA`_o[MGCa`WoYMFcAHCaQHCqXKeSIdQhWkVJdr@_qXLEbPkvCA@OgWKFBa@GgUKEB@gsYMFB`w{]PGcaPGeRIdQpWkVJdr@cqYLEbXkx???OGGCB@`?_SIEB@ow_QIDb@gw]OGcQHGcRHcq`SiTJDax[mVKEBHcqYLErXkw[MfRhuGCA`OoWMFCAPOkWKEBPgs[MFBpw{aRIdQpWkVJdqx_qXKebPkuZMFRhs{]NFg??A@@?_WKGDB@p?cSJGCA@CaQHCaXKgSIDQhWkUJdr@gw]_?OOKOGCaPOgUJEBPo{_OGcaPKgSIdap[oWKebPkw[Mfbp|?_OGG?CCBA@?gSKEB`p_s[NGCQHGeRIDQhWmVKERHguZMFRhw}^OGCA@A`OgW_SKFCAPOkWKEbPow]NGcqh[oWKeRHgsZLer`oy\\Mfbp{}^OGCADAaPO??CAA@?oWOIEBa@GgUOGCAHCcQHCqXOgSIdQpWkVJeBPo{_OGCIDAaPGcQLEb_?OOKOGCaPOgUJEBPo{_OGcaPKgSIdap[oWKebPkw[Mfbp|?_OgSQHCbPgsaQ?@@?o_OIDB@_w[WLFBq@CaQHcq`SiUJdr@cqYLer`sy]NfsADA`PGcYLGdQhW_SKFCAPOkWKEbPow]NGcqh[oWKeRHgsZLer`oy\\Mfbp{}^OgSQLEdQhg??A@@?_WKGDB@p?cSJGCA@CaQHCaXKgSIDQhWkUJdr@gw]OGCADA`OgcQHEbPhCaTIdRHcy\\N??__W_OHCa`OkUKEb`x?_PHCaXOgTJDax_oXLEbXow\\NFby@A`PGcQLEcQHSiTKeRHcy\\P??__WOGDA`_o[MKEb`x?aPHCqXOiTJDqx_qXLErXoy\\NFry@?`OgSQHEbPhCaPIdQhSqXKfRhtA`QhSp?gWMGCa`WoWLEb`o{]PHdQx_oXKeRPguZLfB`sy\\NFbx{}`PGcYLGcQHSiTKeRhtIhT???OGGCB@`?gWMGCa`X?_OGcQPGcRHdA`OiTJDap[mWLFBq@?_OgSIHCaPgsYLGcQhSiXKfRhtA`QhTIdQlUj_?OOKOGCaPOgUJEBPo{_OGcaPKgSIdap[oWKebPkw[Mfbp|?`OgcQHEbPhCaTIeRHcy\\MgSiTQhSjTitYpX??__WOGDA`_o[MKEb`x?aPHCqXOiTJDqx_qXLErXoy\\NFry@A`OgcQLEbQHSiTKeRhsy\\OgSiTIhSjTitapYlUjW_SKFCAPOkWKEbPow]NGcqh[oWKeRHgsZLer`oy\\Mfbp{}^OgcYLEdQhcy\\MhTIdYlUkUJDitYlVJdqy????OGCCA@?oWKGCA`OoWMFCAPOk_OGcaPKgSIdap[oYMFcA@?`PGcaPGdRHdA`OhTIdapWlVJg???@?_OOGCB@_o_OIDB@_w[OHDAq@?aQHCq`OiUJDr@gw]OGCADCaQHCaTKeSIDAdSiUJDat[mWKEC?AA@`?_SIEB@ow_QIDb@_sYMFBpx?aPHCqXOiTJDqxcu\\NgCADAaPhCaTIeRiDAdQiTjDatYmVkEB@_pWkUK?AA@`?_SIEB@ow_QIDb@_sYMFBpx?aPHCqXOiTJDqxcu\\NgCADAaPhCaTIeRiDAdQiTjDatYmVkEJHd???OGGCB@`?gWMGCAPGgSJDb@gw]OHDAp_oXLEbXow\\NFby@AaPGsYPIeRHsy`QiTItYpYmVJtz@_oWKUJDcqXKeZLf???OGGCB@`?gWMGCAPGgSJDb@gw]OHDAp_oXLEbXow\\NFby@AaPGsYPIeRHsy`QiTItYpYmVJtz@_pWkeRHerXlEcOIEBa@?cQIDApWoWKEbPgw[MFbpxCeTJeBHcsZLfBhs{^NgSIHEbPhSiXMfRiTIhUjTjTix]nVkEB@apWkeRLesYLUjUGDB@p?_QHDA`WkWKEBPgs[MFBpw{aRIdr@cqYLer`sy]NfsIDCbPgsiTKfRhtIdSjTitit[nVjuRLgsYlUrY????OGCCA@?oWKGCA`OoWMFCAPOk_OGcaPKgSIdap[oYMFcA@?`PGcaPGdRHdA`OhTIdapWlVJeB@_pWkeRHerYLEbPitYlerXmvZo???A@?__OGEB@`?_SIEB@ow_QIDcA@CcQHdA`SkUJeBPo{_OGCIHCcQHCiXKgSIDIhSkUJDix[oWKEJDcqXKuZPgsYlUjXkuZLuz\\ow_?OOKGCA`OoWMFCAPOkWKEbPow]NGCQHGeRIDQhWmVKerh|?_OgSQLGcQhSq\\OgSiTQlWkUjTq|_oWkUJHcrXkubPgtYlUjXkuZluz`ox[mW?CCBA@?gSKEB`p?cSJEB@gs[MFbq@CaQHcq`SiUJdrHky^OGCIDCbQHCiTKfSIDIdSjUJDit[nWKURHerYLEbTitZLer\\mvZmFBdqy\\O??CAA@?oWOIEBa@?cQIDApWoYMFcAPOkWKERPgu[MFRpw}_OgcQLEcQhcq\\MgSidQlUkUjdq|]oWkURHcrXkuZPgtYlUrXmv[MVJhsy\\mv[??@?__OKECA`_w_OHCa`OkUKEb`x?cSJEB@csYLfB`s{]NgCIHCbPhCiXKfRiDIhSjTjDix[nVkEJDcqXKuZLgsYlUrXkvZlvBdqy\\MfZluz]NG_SKFCA@GcSIDap_oWLEbPow[NFbqHKiVKERHguZMFRhw}^OgSQLEbQhSq\\MfSiTQlUjUjTq|]nWkUJHcrXlEjTiuZLuzdqx\\MfZlu{]nVjuGDB@p?_QHDA`WkWKEBPgs[MFBpw{aRIdr@cqYLer`sy]NfsIDCbPgsiTKfRhtIdSjTitit[nVjuRLgtYlerXmvZmVJdsy\\mvZlw|]nVjx{}', loops=False)
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(G, 1, result='vertex_classes')
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, G, cardinality=1)
        'Refinable'
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(G, 2, result='vertex_classes')
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, G, cardinality=1)
        'Correct'
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(G, 3, result='vertex_classes') # long time
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, G, cardinality=1) # long time
        'Correct'
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(G, 1)
        sage: len(res) == 6
        True
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, G) # long time
        'Refinable'
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(G, 2)
        sage: len(res) == 272
        True
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, G) # long time
        'Correct'
        sage: res = sage.graphs.weisfeiler_lehman.WeisfeilerLehman(G, 3) # long time
        sage: sage.graphs.weisfeiler_lehman.check_orbit_correctness(res, G) # long time
        'Correct'

    """
    if not isinstance(G, SageGraph):
        raise TypeError
    if result not in ['edge_classes', 'vertex_classes', 'graph']:
        raise ValueError("Unknown output format")
    g_temp = _Graph(G, edge_labels)
    if(G.has_multiple_edges()):
        edge_labels = True
    g = g_temp
    #print(g.graph.edges())
    g._relabel_map = g.relabel()
    g.set_vertex_coloring(partition, g._relabel_map)
    cdef vector[GraphNode] res = _sageGraphToLists(g.graph, g.vertex_coloring, edge_labels)
    cdef unordered_map[int, vector[pair[int,int]]] coloring = k_WL(res, k, g.vertex_coloring)
    inverted_relabel_map = {v:k for k,v in g._relabel_map.iteritems()}
    if result == 'graph':
        formatted_result = G.union(G.complement())
        formatted_result.allow_loops(True)
        formatted_result.weighted(True)
        for v in formatted_result:
            formatted_result.add_edge(v,v)
        formatted_result = formatted_result.to_directed()
        c = 0
    else:
        formatted_result = []
    for p in coloring:
        if result == 'vertex_classes':
            l = [inverted_relabel_map[el.first] for el in p.second if el.first == el.second]
            if not l:
                continue
            formatted_result.append(l)
        elif result == 'edge_classes':
            l = [(inverted_relabel_map[el.first],inverted_relabel_map[el.second]) for el in p.second]
            formatted_result.append(l)
        elif result == 'graph':
            for el in p.second:
                formatted_result.set_edge_label(inverted_relabel_map[el.first], inverted_relabel_map[el.second], c)
            c += 1
    return formatted_result
