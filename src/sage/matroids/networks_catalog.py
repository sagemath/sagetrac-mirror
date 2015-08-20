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
    Return Fano network associated with the Fano matroid

    The Fano network is a size 7 network.

    EXAMPLES::
    """
    size=7
    nsrc=3
    cons=[[[1,2],[1,2,4]],[[2,3],[2,3,5]],[[4,5],[4,5,6]],[[3,4],[3,4,7]],[[1,6],[3,1,6]],[[6,7],[2,6,7]],[[5,7],[1,5,7]]]
    return NCinstance(size,nsrc,conslist2dict(cons))


def ButterflyNet():
    r"""
    Return Butterfly network associated with the $U^2_3$ matroid

    The Butterfly network is a size 3 network.

    EXAMPLES::
    """
    size=3
    nsrc=2
    cons=[[[1,2],[1,2,3]],[[1,3],[1,2,3]],[[2,3],[1,2,3]]]
    #cons={0:[set([1]),set([1,4,6])],1:[set([2]),set([2,5,7])],2:[set([4,5]),set([3,4,5])],3:[set([3]),set([3,8,9])],4:[set([6,8]),set([2,6,8])],5:[set([7,9]),set([1,7,9])]}
    return NCinstance(size,nsrc,conslist2dict(cons))

def NonFanoNet():
    size=7
    nsrc=3
    cons=[[[1,2,3],[1,2,3,4]],[[1,2],[1,2,5]],[[1,3],[1,3,6]],[[2,3],[2,3,7]],[[4,5],[3,4,5]],[[4,6],[2,4,6]],[[4,7],[1,4,7]]]#,[[5,6,7],[1,2,3,5,6,7]]]
    return NCinstance(size,nsrc,conslist2dict(cons))

def VamosNet():
    size=8
    nsrc=4
    cons=[[[1,2,3,4],[1,2,3,4,5]],[[1,2,5],[1,2,5,6]],[[2,3,6],[2,3,6,7]],[[3,4,7],[3,4,7,8]],[[4,8],[2,4,8]],[[2,3,4,8],[1,2,3,4,8]],[[1,4,5,8],[1,2,3,4,5,8]],[[1,2,3,7],[1,2,3,4,7]],[[1,5,7],[1,3,5,7]]]
    return NCinstance(size,nsrc,conslist2dict(cons))

def MNet():
    size=12
    nsrc=4
    cons=[[[1,2],[1,2,5,6]],[[3,4],[3,4,7,8]],[[6,7],[6,7,9,10,11,12]],[[5,8,9],[1,3,5,8,9]],[[5,8,10],[1,4,5,8,10]],[[5,8,11],[2,3,5,8,11]],[[5,8,12],[2,4,5,8,12]]]
    return NCinstance(size,nsrc,conslist2dict(cons))

def U23Net():
    return ButterflyNet()

def U24Net():
    size=4
    nsrc=2
    return [[[1,2],[1,2,3]],[[1,3],[1,3,4]],[[2,3],[1,2,3]],[[3,4],[1,3,4]],[[1,4],[1,2,4]],[[2,4],[1,2,4]]]

def U2kNet(k):
    size=k
    nsrc=2
    cons=[]
    for s2 in Combinations(range(1,k+1),2):
        for s3 in Combinations(range(1,k+1),3):
            if set(s2).issubset(set(s3)):
                cons.append([s2,s3])
    return NCinstance(size,nsrc,conslist2dict(cons))


