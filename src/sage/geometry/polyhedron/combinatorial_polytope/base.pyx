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
from sage.rings.integer import Integer
from sage.graphs.graph                      import Graph

from .hasse_diagram cimport   CombinatorialPolytope_ptr, init_CombinatorialPolytope, dimension, edges, f_vector, ridges, delete_CombinatorialPolytope

#TODO take care of the empty polyhedron, which does not have vertices
cdef class CombinatorialPolytope:
    cdef CombinatorialPolytope_ptr _C
    cdef tuple _V
    cdef dict _Vinv
    r"""
    A class of atomic and coatiomic Eulerian lattices.

    One must give
    - an incidence_matrix (with rows corresponding to the facets)
    or
    - facets as list of vertices, where the vertices are labeld 0,...n
    - nr_vertices or vertices
    
    EXAMPLE::
    
        sage: P = polytopes.permutahedron(7)
        sage: C = sage.geometry.combinatorial_polytope.base.CombinatorialPolytope(incidence_matrix=P.incidence_matrix())
        sage: C.f_vector()
        (1, 5040, 15120, 16800, 8400, 1806, 126, 1)
    """
    def __init__(self, data, vertices=None):
        if vertices:
            self._V    = tuple(vertices)
            self._Vinv = { v : i for i,v in enumerate(self._V) }
        else:
            self._V = None
            self._Vinv = None

        if hasattr(data,"incidence_matrix"):#TODO: Better check for incidence_matrix
            data = data.incidence_matrix()

        if hasattr(data,"nrows"):#TODO: Better check for matrix
            rg = range(data.nrows())
            tup =  tuple(tuple(data[i,j] for i in rg) for j in range(data.ncols()) if not all(data[i,j] for i in rg))#transpose and get rid of trivial inequalites (which all vertices satisfie)
            self._C = init_CombinatorialPolytope(tup)
        else:
            if self._V is None:
                vertices = sorted(set.union(*map(set,data)))
                nr_vertices = len(vertices)
                if vertices != range(len(vertices)):
                    self._V    = tuple(vertices)
                    self._Vinv = { v : i for i,v in enumerate(self._V) }
            else:
                nr_vertices = len(self._V)

            if not self._V is None:
                f = lambda v: self._Vinv[v]
            else:
                f = lambda v: int(v)
            facets = tuple(tuple(f(i) for i in j) for j in data)
            self._C = init_CombinatorialPolytope(facets,nr_vertices)

    def __dealloc__(self):
        r"""
        This function deallocates all the memomory used by the underlying C++-class
        """
        delete_CombinatorialPolytope(self._C)

    def edges(self):
        r"""
        Calculates the edges of the CombinatorialPolytope, i.e. the rank 2 faces.
        
        NOTE: If you want to compute edges and f_vector it is recommended to compute edges first.
        """
        if self._V is not None:
            f = lambda i : self._V[i]
        else:
            f = lambda i : Integer(i)

        return tuple((f(i),f(j)) for i,j in edges(self._C))
    def edge_graph(self):
        return Graph(self.edges(),format="list_of_edges")
    def dimension(self):
        return dimension(self._C)

    def ridges(self):
        r"""
        Calculates the ridges of the CombinatorialPolytope, i.e. the rank 2 faces. Those are given as tuples of facets.
        
        E.g. a ridge (1,2) corresponds to the meet of facet[1] and facet[2].
        
        NOTE: If you want to compute ridges and f_vector it is recommended to compute ridges first.
        """
        return tuple((Integer(i),Integer(j)) for i,j in ridges(self._C))
    def ridge_graph(self):
        return Graph(self.ridges(),format="list_of_edges")
    def f_vector(self):
        r"""
        Calculates the f_vector of the CombinatorialPolytope, i.e. the vector containing the nr of faces of each rank.
        
        NOTE: If you also want to compute edges or ridges, it is recommended to do that first.
        """
        return tuple(Integer(i) for i in f_vector(self._C))
        
    def face_lattice(self):
        pass
        
