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
    n = G.vertices()

    if n == 0:
        return [[]]

    if n % 2 == 1:
        raise ValueError("there is no perfect matching for a graph with an odd number of vertices")

    g = G.vertices()[0]
    
