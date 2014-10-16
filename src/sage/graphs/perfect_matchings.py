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

from sage.misc.decorators import options
from sage.misc.cachefunc import cached_method


def perfect_matchings(G):
    r"""
    Computes all perfect matchings or dimer coverings or 1-factors of the graph G

    INPUT:

    * ``G`` - a graph with an even number of vertices

    EXAMPLES::

    sage: G = graphs.GridGraph([3,2])
    sage: perfect_matchings(G)
    [[((1, 1), (2, 1)), ((1, 0), (2, 0)), ((0, 0), (0, 1))],
    [((2, 0), (2, 1)), ((1, 0), (1, 1)), ((0, 0), (0, 1))],
    [((2, 0), (2, 1)), ((0, 1), (1, 1)), ((0, 0), (1, 0))]]
    sage: G = graphs.GridGraph([3,3])
    sage: perfect_matchings(G)
    Traceback (most recent call last)
    ...
    ValueError: there is no perfect matching for a graph with an odd number of vertices

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
        H = G.copy()
        H.delete_vertices([g,h])
        P = perfect_matchings(H)
        for p in P:
            p.append((g,h))
        PP += P
    return PP


    
