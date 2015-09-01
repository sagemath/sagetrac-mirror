r"""
Documentation for the network coding instances in the catalog

This module contains implementations for many of the functions accessible
through :mod:`matroids. <sage.matroids.networks_catalog>` 
(type those lines in Sage and hit ``tab`` for a list).

The docstrings include educational information about each named network with
the hopes that this class can be used as a reference. However, for a more
comprehensive list of properties we refer to [DFZ]_.


AUTHORS:

- Jayant Apte (2015-08-12): initial version

REFERENCES
==========

..  [DFZ] Dougherty, R.; Freiling, C.; Zeger, K., Networks, Matroids, and Non-Shannon Information Inequalities, Information Theory, IEEE Transactions on , vol.53, no.6, pp.1949,1969, June 2007

Functions
=========
"""
#*****************************************************************************
#       Copyright (C) 2015 Jayant Apte <jayant91089@gmail.com >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matroids.network_coding import NCinstance
from sage.matroids.network_coding_helpers import*

def FanoNet():
    r"""
    Return Fano network associated with the Fano matroid. 

    The Fano network is a size 7 network. Scalar solvable only 
    over finite fields of even characteristic.
    
    EXAMPLES::

        sage: F=FanoNet()
        sage: F.constraints
        {0: [{1, 2}, {1, 2, 4}],
         1: [{2, 3}, {2, 3, 5}],
         2: [{4, 5}, {4, 5, 6}],
         3: [{3, 4}, {3, 4, 7}],
         4: [{1, 6}, {1, 3, 6}],
         5: [{6, 7}, {2, 6, 7}],
         6: [{5, 7}, {1, 5, 7}]}
        sage: F.size
        7
        sage: F.nsrc
        3
    """
    size=7
    nsrc=3
    cons=[[[1,2],[1,2,4]],[[2,3],[2,3,5]],[[4,5],[4,5,6]],[[3,4],[3,4,7]],[[1,6],[3,1,6]],[[6,7],[2,6,7]],[[5,7],[1,5,7]]]
    return NCinstance(size,nsrc,conslist2dict(cons))


def ButterflyNet():
    r"""
    Return Butterfly network associated with the $U^2_3$ matroid.

    The Butterfly network is a size 3 network with 2 sources. Smallest example
    where network coding outperforms routing.
    
    EXAMPLES::
    
        sage: from sage.matroids.networks_catalog import *
        sage: B=ButterflyNet()
        sage: B.constraints
        {0: [{1, 2}, {1, 2, 3}], 1: [{1, 3}, {1, 2, 3}], 2: [{2, 3}, {1, 2, 3}]}
        sage: B.size
        3
        sage: B.nsrc
        2  
    """
    size=3
    nsrc=2
    cons=[[[1,2],[1,2,3]],[[1,3],[1,2,3]],[[2,3],[1,2,3]]]
    #cons={0:[set([1]),set([1,4,6])],1:[set([2]),set([2,5,7])],2:[set([4,5]),set([3,4,5])],3:[set([3]),set([3,8,9])],4:[set([6,8]),set([2,6,8])],5:[set([7,9]),set([1,7,9])]}
    return NCinstance(size,nsrc,conslist2dict(cons))

def NonFanoNet():
    r"""
    Return Non-fano network associated with the Non-fano matroid.

    The Non-fano network is a size 7 network with 3 sources. Scalar solvable 
    only over finite fields of odd characteristic.
    
    EXAMPLES::
    
        sage: NF=NonFanoNet()
        sage: NF.constraints
        {0: [{1, 2, 3}, {1, 2, 3, 4}],
         1: [{1, 2}, {1, 2, 5}],
         2: [{1, 3}, {1, 3, 6}],
         3: [{2, 3}, {2, 3, 7}],
         4: [{4, 5}, {3, 4, 5}],
         5: [{4, 6}, {2, 4, 6}],
         6: [{4, 7}, {1, 4, 7}]}
        sage: NF.size
        7
        sage: NF.nsrc
        3
    """
    size=7
    nsrc=3
    cons=[[[1,2,3],[1,2,3,4]],[[1,2],[1,2,5]],[[1,3],[1,3,6]],[[2,3],[2,3,7]],[[4,5],[3,4,5]],[[4,6],[2,4,6]],[[4,7],[1,4,7]]]#,[[5,6,7],[1,2,3,5,6,7]]]
    return NCinstance(size,nsrc,conslist2dict(cons))

def VamosNet():
    r"""
    Return Vamos network associated with the Vamos matroid.

    The Vamos network is a size 8 network with 4 sources. Not scalar 
    (or vector) linear solvable over any finite field.
    
    EXAMPLES::

        sage: from sage.matroids.networks_catalog import *
        sage: V=VamosNet()
        sage: V.constraints
        {0: [{1, 2, 3, 4}, {1, 2, 3, 4, 5}],
         1: [{1, 2, 5}, {1, 2, 5, 6}],
         2: [{2, 3, 6}, {2, 3, 6, 7}],
         3: [{3, 4, 7}, {3, 4, 7, 8}],
         4: [{4, 8}, {2, 4, 8}],
         5: [{2, 3, 4, 8}, {1, 2, 3, 4, 8}],
         6: [{1, 4, 5, 8}, {1, 2, 3, 4, 5, 8}],
         7: [{1, 2, 3, 7}, {1, 2, 3, 4, 7}],
         8: [{1, 5, 7}, {1, 3, 5, 7}]}
        sage: V.size
        8
        sage: V.nsrc
        4
        """
    size=8
    nsrc=4
    cons=[[[1,2,3,4],[1,2,3,4,5]],[[1,2,5],[1,2,5,6]],[[2,3,6],[2,3,6,7]],[[3,4,7],[3,4,7,8]],[[4,8],[2,4,8]],[[2,3,4,8],[1,2,3,4,8]],[[1,4,5,8],[1,2,3,4,5,8]],[[1,2,3,7],[1,2,3,4,7]],[[1,5,7],[1,3,5,7]]]
    return NCinstance(size,nsrc,conslist2dict(cons))

def MNet():
    r"""
    Return M-network.
    
    The M-network is a size 12 network with 4 sources. Not scalar linear 
    solvable over any finite field but is vector linear solvable over any
    field (as it is routing solvable).
    
    EXAMPLES::

        sage: from sage.matroids.networks_catalog import *
        sage: M=MNet()
        sage: M.constraints
        {0: [{1, 2}, {1, 2, 5, 6}],
         1: [{3, 4}, {3, 4, 7, 8}],
         2: [{6, 7}, {6, 7, 9, 10, 11, 12}],
         3: [{5, 8, 9}, {1, 3, 5, 8, 9}],
         4: [{5, 8, 10}, {1, 4, 5, 8, 10}],
         5: [{5, 8, 11}, {2, 3, 5, 8, 11}],
         6: [{5, 8, 12}, {2, 4, 5, 8, 12}]}
        sage: M.size
        12
        sage: M.nsrc
        4
    """
    size=12
    nsrc=4
    cons=[[[1,2],[1,2,5,6]],[[3,4],[3,4,7,8]],[[6,7],[6,7,9,10,11,12]],[[5,8,9],[1,3,5,8,9]],[[5,8,10],[1,4,5,8,10]],[[5,8,11],[2,3,5,8,11]],[[5,8,12],[2,4,5,8,12]]]
    return NCinstance(size,nsrc,conslist2dict(cons))

def U2kNet(k):
    r"""
    Return the $U^2_k$ network associated with the uniform matroid $U^2_k$.
    
    The $U^2_k$ network is scalar linear solvable over $GF(q)$ iff $q>=k-1$

    INPUT:
        
    - ``k`` -- An integer
    
    EXAMPLES::
    
        sage: U23=U2kNet(3)
        sage: B=ButterflyNet()
        sage: U23.constraints==B.constraints and U23.nsrc==B.nsrc and U23.size==B.size
        True
        sage: U24=U2kNet(4)
        sage: U24.constraints
        {0: [{1, 2}, {1, 2, 3}],
         1: [{1, 2}, {1, 2, 4}],
         2: [{1, 3}, {1, 2, 3}],
         3: [{1, 3}, {1, 3, 4}],
         4: [{1, 4}, {1, 2, 4}],
         5: [{1, 4}, {1, 3, 4}],
         6: [{2, 3}, {1, 2, 3}],
         7: [{2, 3}, {2, 3, 4}],
         8: [{2, 4}, {1, 2, 4}],
         9: [{2, 4}, {2, 3, 4}],
         10: [{3, 4}, {1, 3, 4}],
         11: [{3, 4}, {2, 3, 4}]}
        sage: U24.size
        4
        sage: U24.nsrc
        2
    """
    size=k
    nsrc=2
    cons=[]
    for s2 in Combinations(range(1,k+1),2):
        for s3 in Combinations(range(1,k+1),3):
            if set(s2).issubset(set(s3)):
                cons.append([s2,s3])
    return NCinstance(size,nsrc,conslist2dict(cons))


