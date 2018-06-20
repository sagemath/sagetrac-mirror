# distutils: language = c++
from __future__ import print_function
from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free
from sys import version_info

def getiterator(el):
    if version_info >= (3):
        return el.items()
    else:
        return el.iteritems()
    
cdef vector[GraphNode] sageGraphToLists(G, partition = [], relabel_map={}):
    cdef int maxIdx = 0
    partition_dict = {v:i for i,x in enumerate(partition,1) for v in x}
    
    cdef vector[GraphNode] nodeArray;
    nodeArray.resize(G.order())
    
    for v in G:
        l = relabel_map.setdefault(v, maxIdx)
        if l == maxIdx:
            maxIdx += 1
        nodeArray[l].idx = l
        nodeArray[l].color = partition_dict.get(v, 0)
        nodeArray[l].adj_list.reserve(len(G[v]))
    directed = G.is_directed()
    
    for v,u in G.edges(labels=False):
        v = relabel_map[v]
        u = relabel_map[u]
        nodeArray[v].adj_list.push_back(<int>u)
        if not directed:
            nodeArray[u].adj_list.push_back(<int>v)
    return nodeArray

def prova(G, partition=[]):
    relabel_map = {}
    cdef vector[GraphNode] res = sageGraphToLists(G, partition, relabel_map)
    print("%d" % prova2(res))
    
