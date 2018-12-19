r"""
Several algorithms working implicitly with the hasse_diagram (of a polytope), including calculating the f_vector.

This is a wrapper for the functions in hasse_diagram.cc.

This computes implicitely a finite atomic and coatomic lattices, where every interval of length two has at least 4 elements.
(Exactly 4 is known as the diamond property).
In particular this module calculates quickly the f_vector of polytopes. The input must be a tuple of coatoms given each by a tuple of atoms. The atoms must be labeled 0,...,n.

AUTHOR:

- Jonathan Kliem (2018-12)
"""


#*****************************************************************************
#       Copyright (C) 2018 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import

from .hasse_diagram cimport   CombinatorialType_ptr, init_CombinatorialType, dimension, edges, f_vector, ridges, delete_CombinatorialType


cdef class CombinatorialType:
    cdef CombinatorialType_ptr _C
    r"""
    The class of a Combinatorial type of an atomic and coatiomic lattice, where every interval of length 2 has at least 4 elements.
    
    One must give
    - an incidence_matrix (with rows corresponding to the facets)
    or
    - facets as list of vertices, where the vertices are labeld 0,...n
    - nr_vertices or vertices
    
    EXAMPLE::
    
        sage: P = polytopes.permutahedron(7)
        sage: C = sage.geometry.combinatorial_type.base.CombinatorialType(incidence_matrix=P.incidence_matrix())
        sage: C.f_vector()
        (1L, 5040L, 15120L, 16800L, 8400L, 1806L, 126L, 1L)
    """
    def __init__(self,facets=None,vertices=None,nr_vertices=None,incidence_matrix=None):
        
        if incidence_matrix:
            rg = range(incidence_matrix.nrows())
            tup =  tuple(tuple(incidence_matrix[i,j] for i in rg) for j in range(incidence_matrix.ncols()) if not all(incidence_matrix[i,j] for i in rg))#transpose and get rid of trivial inequalites (which all vertices satisfie)
            self._C = init_CombinatorialType(tup)
        else:
            if vertices:
                nr_vertices = len(vertices)
            if facets and nr_vertices:
                try:
                    facets = tuple(tuple(int(i) for i in j) for j in facets)
                except:
                    raise ValueError("facets must be given as tuple of tuples of vertices")
                self._C = init_CombinatorialType(facets,nr_vertices)
            else:
                raise ValueError("Not sufficient information provided to obtain a CombinatorialType")
    def __del__(self):
        delete_CombinatorialType(self._C)
    def edges(self):
        r"""
        Calculates the edges of the CombinatorialType, i.e. the rank 2 faces.
        
        NOTE: If you want to compute edges and f_vector it is recommended to compute edges first.
        """
        return edges(self._C)
    def dimension(self):
        return dimension(self._C)
    def ridges(self):
        r"""
        Calculates the ridges of the CombinatorialType, i.e. the rank 2 faces. Those are given as tuples of facets.
        
        E.g. a ridge (1,2) corresponds to the meet of facet[1] and facet[2].
        
        NOTE: If you want to compute ridges and f_vector it is recommended to compute ridges first.
        """
        return ridges(self._C)
    def f_vector(self):
        r"""
        Calculates the f_vector of the CombinatorialType, i.e. the vector containing the nr of faces of each rank.
        
        NOTE: If you also want to compute edges or ridges, it is recommended to do that first.
        """
        return f_vector(self._C)
