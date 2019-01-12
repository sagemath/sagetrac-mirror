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
from sage.graphs.digraph import DiGraph
from sage.combinat.posets.lattices import FiniteLatticePoset


from .hasse_diagram cimport   CombinatorialPolyhedron_ptr, init_CombinatorialPolyhedron, dimension, edges, f_vector, ridges, incidences, record_all_faces, get_faces, get_flag, delete_CombinatorialPolyhedron

#TODO take care of the empty polyhedron, which does not have vertices
cdef class CombinatorialPolyhedron:
    cdef CombinatorialPolyhedron_ptr _C
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
        sage: C = sage.geometry.combinatorial_polytope.base.CombinatorialPolyhedron(incidence_matrix=P.incidence_matrix())
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
            self._C = init_CombinatorialPolyhedron(tup)
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
            self._C = init_CombinatorialPolyhedron(facets,nr_vertices)

    def __dealloc__(self):
        r"""
        This function deallocates all the memomory used by the underlying C++-class
        """
        delete_CombinatorialPolyhedron(self._C)

    def edges(self):
        r"""
        Calculates the edges of the CombinatorialPolyhedron, i.e. the rank 2 faces.
        
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
        Calculates the ridges of the CombinatorialPolyhedron, i.e. the rank 2 faces. Those are given as tuples of facets.
        
        E.g. a ridge (1,2) corresponds to the meet of facet[1] and facet[2].
        
        NOTE: If you want to compute ridges and f_vector it is recommended to compute ridges first.
        """
        return tuple((Integer(i),Integer(j)) for i,j in ridges(self._C))
    def ridge_graph(self):
        return Graph(self.ridges(),format="list_of_edges")
    def f_vector(self):
        r"""
        Calculates the f_vector of the CombinatorialPolyhedron, i.e. the vector containing the nr of faces of each rank.
        
        NOTE: If you also want to compute edges or ridges, it is recommended to do that first.
        """
        return tuple(Integer(i) for i in f_vector(self._C))
    def _record_all_faces(self):
        record_all_faces(self._C)
    def faces(self,dimension, facet_repr=False):#TODO fix the cpp-function for dimension 0,1,self.dimension -1, self.dimension
        r"""
        Gets all k-faces for specified dimenion k.

        By default faces are given as tuple of vertices, but one may also choose tuple of facets instead.
        """
        if (facet_repr):
            facet_repr = 1
        else:
            facet_repr = 0
        return tuple(tuple(Integer(j) for j in i) for i in get_faces(self._C, dimension, facet_repr))
    def incidences(self,dimension_one,dimension_two):
        r"""
        Gets a tuple of all incidens between faces of dimension dimension_one and dimension_two.

        Incidences are given as tuple of integers, where the integer corresponds to the order according to self.faces(dimension)
        """
        return tuple((Integer(i),Integer(j)) for i,j in incidences(self._C,dimension_one,dimension_two))
    def face_lattice(self,vertices=True,facets=False):
        r"""
        Generates the face-lattice.

        The faces in the face-lattice will be labeled by default by tuples of vertices. They can also be labeled by tuple of facets or both or as 0, 1, 2, ... where the order corresponds to the order of self.faces[0],self.faces[1],self.faces[2],...
        """

        f_vector = self.f_vector()
        self._record_all_faces()
        dimension = self.dimension()
        dic = {}
        for k in range(-1,dimension+1):
            faces = (self.faces(k),self.faces(k,facet_repr=True))
            dic[k] = tuple((faces[0][i],faces[1][i]) for i in range(f_vector[k+1]))
        edges = tuple((i[0] + sum(f_vector[:k+1]),i[1] + sum(f_vector[:k+2])) for k in range(-1,dimension) for i in self.incidences(k,k+1))
        all_faces = tuple(i for k in range(-1,dimension+1) for i in dic[k])
        if vertices and facets:
            f = lambda i : all_faces[i]
        elif vertices:  
            f = lambda i : all_faces[i][0]
        elif facets:
            f = lambda i : all_faces[i][1]
        else:
            f = lambda i : i
        V = tuple(range(sum(f_vector)))
        D = DiGraph([V, edges], format='vertices_and_edges')
        D.relabel(f)
        return FiniteLatticePoset(D);
    def get_flag(self,flag):
        flag = sorted(set(flag))
        flag = tuple(i for i in flag if i in range(-1,self.dimension()+1))
        if flag == (-1,):
            return 1
        flag = tuple(i for i in flag if i in range(0,self.dimension()+1))
        return Integer(get_flag(self._C,flag))
