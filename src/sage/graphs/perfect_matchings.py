"""
Perfect Matchings

This module contains the following methods:

AUTHORS:

- Arvind Ayyer - original implementation

REFERENCE:



Methods
-------
"""

#*****************************************************************************
#                       Copyright (C) 2014 Arvind Ayyer
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

def perfect_matchings(G):
    r"""
    Computes all perfect matchings of the graph G



    """
    n = G.num_verts()

    if n == 0:
        return [[]]

    if n % 2 == 1:
        raise ValueError("there is no perfect matching for a graph with an odd number of vertices")

    g = G.vertices()[0]
    N = G.neighbors(g)
    PP = []
    for h in N:
        H = copy(G)
        H.delete_vertices([g,h])
        P = perfect_matchings(H)
        for p in P:
            p.append((g,h))
        PP += P
    return PP


    
