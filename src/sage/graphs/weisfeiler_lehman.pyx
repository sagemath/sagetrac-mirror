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
                vc = [set(int(relabel_map[v]) for v in k) for k in vc]
            
            if self._first_level_vertices:
                vs = set([int(relabel_map[v]) for k,v in getiterator(self._first_level_vertices)])
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


cdef vector[GraphNode] sageGraphToLists(G, partition = [], has_edge_labels=False):
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

def WeisfeilerLehman(G, k, partition=[], edge_labels=False, result_cardinality=1):
    if not isinstance(G, SageGraph):
        raise TypeError
    g_temp = Graph(G, edge_labels)
    if(G.has_multiple_edges()):
        edge_labels = True
    g = g_temp
    
    g._relabel_map = g.relabel()
    print(g._relabel_map)
    g.set_vertex_coloring(partition, g._relabel_map)
    
    #I should add support for labels, and convert vertex labels and edge labels to integers. Call the functions you developed
    cdef vector[GraphNode] res = sageGraphToLists(g.graph, g.vertex_coloring, edge_labels)
#    print(g.graph.vertices())
#    print(g.graph.edges())
#    print(g.vertex_coloring)
#    print(res)
#    print(g._relabel_map)
    cdef unordered_map[int, vector[pair[int,int]]] coloring = k_WL(res, k, g.vertex_coloring)
    inverted_relabel_map = {v:k for k,v in g._relabel_map.iteritems()}
    resultDict = {}
    
    #Normally, one would check for equality of the results. If using initial partitions though, one should check for an automorphism between the colors.
    #That is, c1(tuple) = g(c2(tuple)) with g bijective, since it could happen that the initial coloring produces a different initial order that doesn't cause any issue with k-WL,
    #but produces a different permutation of the final coloring
    for p in coloring:
        if result_cardinality == 1:
            l = [inverted_relabel_map[el.first] for el in p.second if el.first == el.second]
        elif result_cardinality == 2:
            edge = (inverted_relabel_map[el.first],inverted_relabel_map[el.second])
            l = [tuple(el) for el in p.second]
        else:
            raise ValueError("Cardinality above 2 not yet implemented")
        c = p.first
        resultDict[c] = l
    return resultDict
