# distutils: language = c++
from __future__ import print_function
from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string
from libcpp.utility cimport pair
from libcpp cimport bool
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

def prova(G, k, partition=[]):
    relabel_map = {}
    cdef vector[GraphNode] res = sageGraphToLists(G, partition, relabel_map)
    cdef unordered_map[Tuple[int], int] coloring = k_WL(res, k)
    print("Finished")
    resultDict = {}
    
    #Normally, one would check for equality of the results. If using initial partitions though, one should check for an automorphism between the colors.
    #That is, c1(tuple) = g(c2(tuple)) with g bijective, since it could happen that the initial coloring produces a different initial order that doesn't cause any issue with k-WL,
    #but produces a different permutation of the final coloring
    for p in coloring:
        l = [el for el in p.first]
        c = p.second
        resultDict[tuple(l)] = c
    return resultDict
