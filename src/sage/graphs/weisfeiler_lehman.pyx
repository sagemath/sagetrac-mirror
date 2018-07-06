# distutils: language = c++
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
'''
    graph.py

Module graph contains the definition of the Graph class
and utilities dealing with graph objects.
'''

def getiterator(l):
    if version_info >= (3):
        return l.items()
    else:
        return l.iteritems()
        
class Graph(object):
    '''
    For each multiple edge e whose elements have labels l_0,l_1,...,l_(m-1), 
    delete the multiple edge and substitute it with a single edge, with label 
    (l_i, m) where i can be any integer between 0 and m-1
    '''
    def _edge_multiplicity_to_label(self, edge_labels=False):
        from collections import Counter
        G = self.graph
        edge_multiplicity_dict = Counter(G.multiple_edges(labels=False))
       
        G.remove_multiple_edges()
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
            self._edge_multiplicity_to_label(edge_labels)
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
        self._inverted_relabel_map = {v: k for k, v in getiterator(m)}
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
                vc = [[int(relabel_map[v]) for v in k] for k in vc]
            
            if self._first_level_vertices:
                vs = set([int(relabel_map[v]) for k,v in getiterator(self._first_level_vertices)])
                self._vertex_coloring.append(set([int(k) for k in range(self.number_of_vertices) if not k in vs]))
            else:
                vs = set(range(self.number_of_vertices))
                
            for p in vertex_coloring:
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

def _array_to_disjoint_representation(perm):
    if not isinstance(perm,dict):
        perm = {k:v for k,v in enumerate(perm)}
    result = []
    remaining = set(perm)
    while remaining:
        el = remaining.pop()
        cycle = [el]
        el = perm[el]
        while el not in cycle:
            remaining.remove(el)
            cycle.append(el)
            el = perm[el]
        if len(cycle) == 1:
            continue
        result.append(cycle)
    return result
            

#def autgrp(g, partition=[], edge_labels=False):
#    '''
#    Compute the automorphism group of a graph.

#    *g*
#        A PyNauty Graph object.

#    return -> (generators, grpsize1, grpsize2, orbits, numorbits)
#        For the detailed description of the returned components, see
#        Nauty's documentation.
#    '''
#    if not isinstance(g, SageGraph):
#        raise TypeError
        
    

#    g_temp = Graph(g, edge_labels)
#    if(g.has_multiple_edges()):
#        edge_labels = True
#    g = g_temp
#    if edge_labels:
#        g.graph, g._first_level_vertices = g.graph.remove_labels()
#        g.number_of_vertices = g.graph.order()
#    g._relabel_map = g.relabel()
#    g.set_vertex_coloring(partition, g._relabel_map)

#    #gens, grpsize1, grpsize2, orbits, numorbits = nautywrap.graph_autgrp(g)
#    gens = [_array_to_disjoint_representation(gen) for gen in gens]
#    if g._normalized:
#        gens = [[tuple([g._inverted_relabel_map[v] for v in cycle]) for cycle in gen] for gen in gens]
#    if g._first_level_vertices:
#        inverted_first_level_vertihttps://vittorioromeo.info/index/blog/passing_functions_to_functions.htmlces = {v:k for k,v in getiterator(g._first_level_vertices)}
#        gens = [[tuple([inverted_first_level_vertices[v] for v in cycle]) for cycle in gen if cycle[0] in inverted_first_level_vertices] for gen in gens]
    
#    result = {}
#    for k,v in enumerate(orbits):
#        if g._first_level_vertices:
#            if g._inverted_relabel_map[k] not in inverted_first_level_vertices:
#                continue
#        if(g._normalized):
#            k = g._inverted_relabel_map[k]
#            if g._first_level_vertices:
#                k = inverted_first_level_vertices[k]
#        result.setdefault(v, []).append(k)
#    orbits = result.values()
#    numorbits = len(orbits)
#    return gens,grpsize1,grpsize2,orbits,numorbits


def certificate(g):
    '''
    Compute a certificate based on the canonical labeling of vertices.

    *g*
        A Graph object.

    return ->
        The certificate as a byte string.
    '''
    if not isinstance(g, Graph):
        raise TypeError
    #return nautywrap.graph_cert(g)


#def isomorphic(a, b, a_partition=[], b_partition=[], edge_labels=False):
#    '''
#    Determine if two graphs are isomorphic.

#    *a,b*
#        Two Graph objects.

#    return ->
#        True if *a* and *b* are isomorphic graphs, False otherwise,
#    '''
#    if not isinstance(a, SageGraph):
#        raise TypeError
#    if not isinstance(b, SageGraph):
#        raise TypeError
#    a_temp = Graph(a, edge_labels)
#    b_temp = Graph(b, edge_labels)
#    if a.has_multiple_edges() or b.has_multiple_edges():
#        edge_labels = True
#    a = a_temp
#    b = b_temp
#    if edge_labels:
#        edge_labels_list = a.graph.edge_labels()
#        edge_labels_list.extend(b.graph.edge_labels())
#        a.graph, a._first_level_vertices = a.graph.remove_labels(edge_labels_list)
#        b.graph, b._first_level_vertices = b.graph.remove_labels(edge_labels_list)
#        a.number_of_vertices = a.graph.order()
#        b.number_of_vertices = b.graph.order()
#    a._relabel_map = a.relabel()
#    b._relabel_map = b.relabel()
#    a.set_vertex_coloring(a_partition, a._relabel_map)
#    b.set_vertex_coloring(b_partition, b._relabel_map)
#    if a.number_of_vertices != b.number_of_vertices:
#        return False
#    elif list(map(len, a.vertex_coloring)) != list(map(len, b.vertex_coloring)):
#        return False
#    else:
#        return certificate(a) == certificate(b)
    
cdef vector[GraphNode] sageGraphToLists(G, partition = [], relabel_map={}, has_edge_labels=False):
    partition_dict = {v:i for i,x in enumerate(partition,1) for v in x}
    
    cdef vector[GraphNode] nodeArray;
    nodeArray.resize(G.order())
    
    for v in G:
        if v in relabel_map:
            v = int(relabel_map[v])
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
        v = int(relabel_map[v])
        u = int(relabel_map[u])
        nodeArray[v].adj_list.push_back([<int>u, <int>edge_labels_dict[l]])
        if not directed:
            nodeArray[u].adj_list.push_back([<int>v, <int>edge_labels_dict[l]])
    return nodeArray

def prova(G, k, partition=[], edge_labels=False):
    if not isinstance(G, SageGraph):
        raise TypeError

    g_temp = Graph(G, edge_labels)
    if(G.has_multiple_edges()):
        edge_labels = True
    g = g_temp
    #Use the part below if using edge labels is not possible
    #if edge_labels:
    #    g.graph, g._first_level_vertices = g.graph.remove_labels()
    #    g.number_of_vertices = g.graph.order()
    g._relabel_map = g.relabel()
    g.set_vertex_coloring(partition, g._relabel_map)
    
    #I should add support for labels, and convert vertex labels and edge labels to integers. Call the functions you developed
    cdef vector[GraphNode] res = sageGraphToLists(g.graph, g.vertex_coloring, g._relabel_map, edge_labels)
    print(g.graph.vertices())
    print(g.graph.edges())
    print(g.vertex_coloring)
    print(res)
    print(g._relabel_map)
    cdef unordered_map[Tuple[int], int] coloring = k_WL(res, k, g.vertex_coloring)
    
    resultDict = {}
    
    #Normally, one would check for equality of the results. If using initial partitions though, one should check for an automorphism between the colors.
    #That is, c1(tuple) = g(c2(tuple)) with g bijective, since it could happen that the initial coloring produces a different initial order that doesn't cause any issue with k-WL,
    #but produces a different permutation of the final coloring
    for p in coloring:
        l = [el for el in p.first]
        c = p.second
        resultDict[tuple(l)] = c
    return resultDict
