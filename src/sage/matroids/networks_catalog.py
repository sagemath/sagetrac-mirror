r"""
Documentation for the network coding instances in the catalog

This module contains implementations for many of the functions accessible
through :mod:`matroids. <sage.matroids.network_coding.networks_catalog>` and
:mod:`matroids.network_coding.named_networks. <sage.matroids.networks_catalog>`
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
#       Copyright (C) 2015 Jayant Apte <michael@welsh.co.nz >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matroids.network_coding import NCinstance


def FanoNet():
    """
    Return Fano network associated with the Fano matroid

    The Fano network is a size 7 network.

    EXAMPLES::
    """
    size=7
    nsrc=3
    cons=[[[1],[1,2,4]],[[2],[2,4,5]],[[3],[3,5,7]],[[4,5],[4,5,6]],[[3,4],[3,4,7]],[[1,6],[3,1,6]],[[6,7],[2,6,7]],[[5,7],[1,5,7]]]
    return 

def ButterflyNet():
    """
    Return Fano network associated with the Fano matroid

    The Butterfly network is a size 9 network.

    EXAMPLES::
    """
    size=9
    nsrc=2
    cons=[[[1],[1,4,6]],[[2],[2,5,7]],[[4,5],[3,4,5]],[[3],[3,8,9]],[[6,8],[2,6,8]],[[7,9],[1,7,9]]]
    return NCinstance(9,2,cons)

def NonFanoNet():
    return 


