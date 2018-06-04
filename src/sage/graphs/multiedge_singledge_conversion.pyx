# cython: binding=True
"""
Multiedge to singledge

This module contains the following methods:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`multiedge_to_singledge` | Transforms a graph with multiple edges in an equivalent graph with single edges

AUTHORS:

- Dario Asprone - original implementation

REFERENCE:

.. 


Methods
-------
"""

from sage.graphs.graph import Graph, DiGraph
from math import floor, log
from sys import version_info

def getiterator(l):
    if version_info >= (3):
        return l.items()
    else:
        return l.iteritems()

def get_bits_list(x):
    return [b=='1' for b in bin(x)[2:][::-1]]
    
def multiedge_to_singledge(G):
    if G.is_directed():
        newG = DiGraph(loops=True, multiedges=False)
    else:
        newG = Graph(loops=True, multiedges=False)
    G_dict = G.to_dictionary(edge_labels=True, multiple_edges=True)
    edge_list = {}
    max_multiplicity = 0
    for v, neighbor in getiterator(G_dict):
        for u, multiedge in getiterator(neighbor):
            if (u,v) in edge_list:
                continue
            m = len(multiedge)
            edge_list[(v,u)] = m
            max_multiplicity = max(m, max_multiplicity)
    del G_dict
    k = max_multiplicity
    d = int(log(k, 2)) + 1
    first_level_vertices = {}
    for v in G:
        first_level_vertices[v] = (v,0)
        for l in range(1,d):
            if(G.is_directed()): 
                newG.add_edges([((v,l-1),(v,l)),((v,l),(v,l-1))])
            else:
                newG.add_edge((v,l-1),(v,l))
    for (v,u),label in getiterator(edge_list):
        for idx, bit in enumerate(get_bits_list(label)):
            if(bit):
                newG.add_edge((v,idx),(u,idx))
                if G.is_directed():
                    newG.add_edge((u,idx),(v,idx))
    return newG, first_level_vertices

        
        
    
    
