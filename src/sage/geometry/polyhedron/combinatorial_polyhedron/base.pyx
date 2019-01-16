r"""
Several algorithms working implicitly with the hasse_diagram (of a polytope), including calculating the f_vector, dimension, flags, level-sets and even the face lattice.

This is a wrapper for the functions in hasse_diagram.cc.

This computes implicitely a finite atomic and coatomic lattices, where every interval of length two has at least 4 elements.
(Exactly 4 is known as the diamond property).
In particular this module calculates quickly the f_vector of polytopes. The input must be a tuple of coatoms given each by a tuple of atoms. The atoms must be labeled 0,...,n.

EXAMPLES:

Initialization::

    sage: P = polytopes.cube()
    sage: CombinatorialPolyhedron(P)
    The Combinatorial Type of a Polyhedron of dimension 3 with 8 vertices
    
    sage: points = [(1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1)]
    sage: L = LatticePolytope(points)
    sage: CombinatorialPolyhedron(L)
    The Combinatorial Type of a Polyhedron of dimension 3 with 6 vertices

You can also give the facets explicitely::

    sage: CombinatorialPolyhedron(((1,2,3),(1,2,4),(1,3,4),(2,3,4)))
    The Combinatorial Type of a Polyhedron of dimension 3 with 4 vertices

Compute the dimension::

    sage: C = CombinatorialPolyhedron(((1,2,3),(1,2,4),(1,3,4),(2,3,4)))
    sage: C.dimension()
    3

Compute the f-vector::

    sage: P = polytopes.permutahedron(7)
    sage: C = CombinatorialPolyhedron(P)
    sage: C.f_vector()
    (1, 5040, 15120, 16800, 8400, 1806, 126, 1)
    

Retrieve the vertices::

    sage: P = polytopes.cube()
    sage: C = CombinatorialPolyhedron(P)
    sage: C.vertices()
    (A vertex at (-1, -1, -1),
     A vertex at (-1, -1, 1),
     A vertex at (-1, 1, -1),
     A vertex at (-1, 1, 1),
     A vertex at (1, -1, -1),
     A vertex at (1, -1, 1),
     A vertex at (1, 1, -1),
     A vertex at (1, 1, 1))
    
    sage: points = [(1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1)]
    sage: L = LatticePolytope(points)
    sage: C = CombinatorialPolyhedron(L)
    sage: C.vertices()
    (M(1, 0, 0), M(0, 1, 0), M(0, 0, 1), M(-1, 0, 0), M(0, -1, 0), M(0, 0, -1))

Get the edges::

    sage: P = polytopes.cyclic_polytope(3,5)
    sage: C = CombinatorialPolyhedron(P)
    sage: C.edges()
    ((A vertex at (0, 0, 0), A vertex at (4, 16, 64)),
     (A vertex at (1, 1, 1), A vertex at (4, 16, 64)),
     (A vertex at (2, 4, 8), A vertex at (4, 16, 64)),
     (A vertex at (3, 9, 27), A vertex at (4, 16, 64)),
     (A vertex at (0, 0, 0), A vertex at (3, 9, 27)),
     (A vertex at (1, 1, 1), A vertex at (3, 9, 27)),
     (A vertex at (0, 0, 0), A vertex at (2, 4, 8)),
     (A vertex at (1, 1, 1), A vertex at (2, 4, 8)),
     (A vertex at (0, 0, 0), A vertex at (1, 1, 1)))

Get the ridges::

    sage: P = polytopes.cyclic_polytope(3,5)
    sage: C = CombinatorialPolyhedron(P)
    sage: C.edges()
    sage: C.ridges()
    ((An inequality (12, -7, 1) x + 0 >= 0, An inequality (-4, 5, -1) x + 0 >= 0),
     (An inequality (-14, 7, -1) x + 8 >= 0, An inequality (-4, 5, -1) x + 0 >= 0),
     (An inequality (-26, 9, -1) x + 24 >= 0,
      An inequality (-14, 7, -1) x + 8 >= 0),
     (An inequality (12, -7, 1) x + 0 >= 0,
      An inequality (-26, 9, -1) x + 24 >= 0),
     (An inequality (6, -5, 1) x + 0 >= 0, An inequality (12, -7, 1) x + 0 >= 0),
     (An inequality (6, -5, 1) x + 0 >= 0, An inequality (-26, 9, -1) x + 24 >= 0),
     (An inequality (2, -3, 1) x + 0 >= 0, An inequality (6, -5, 1) x + 0 >= 0),
     (An inequality (2, -3, 1) x + 0 >= 0, An inequality (-14, 7, -1) x + 8 >= 0),
     (An inequality (2, -3, 1) x + 0 >= 0, An inequality (-4, 5, -1) x + 0 >= 0))

Compute the edge graph::

    sage: P = polytopes.buckyball()
    sage: C = CombinatorialPolyhedron(P)
    sage: G = C.edge_graph()
    sage: G
    Graph on 60 vertices
    sage: G.diameter()
    9

Compute the ridge graph::

    sage: P = polytopes.permutahedron(7)
    sage: C = CombinatorialPolyhedron(P)
    sage: C.ridge_graph()
    Graph on 126 vertices

Get the k-faces::

    sage: P = polytopes.simplex()
    sage: C = CombinatorialPolyhedron(P)
    sage: C.faces(1)
    ((A vertex at (0, 0, 1, 0), A vertex at (0, 1, 0, 0)),
     (A vertex at (0, 0, 1, 0), A vertex at (1, 0, 0, 0)),
     (A vertex at (0, 1, 0, 0), A vertex at (1, 0, 0, 0)),
     (A vertex at (0, 0, 0, 1), A vertex at (0, 1, 0, 0)),
     (A vertex at (0, 0, 0, 1), A vertex at (1, 0, 0, 0)),
     (A vertex at (0, 0, 0, 1), A vertex at (0, 0, 1, 0)))
    sage: C.faces(1,facet_repr=True)
    ((An equation (1, 1, 1, 1) x - 1 == 0, 0, 3),
     (An equation (1, 1, 1, 1) x - 1 == 0, 1, 3),
     (An equation (1, 1, 1, 1) x - 1 == 0, 2, 3),
     (An equation (1, 1, 1, 1) x - 1 == 0, 0, 2),
     (An equation (1, 1, 1, 1) x - 1 == 0, 1, 2),
     (An equation (1, 1, 1, 1) x - 1 == 0, 0, 1))

Get the incidences between j-faces and k-faces,
this is given according to the enumaration of faces(j) and faces(k)::

    sage: j = 2
    sage: k = 3
    sage: P = polytopes.simplex()
    sage: C.incidences(2,3)
    ((0, 0), (1, 0), (2, 0), (3, 0))
    sage: C.incidences(j,k)
    ((0, 0), (1, 0), (2, 0), (3, 0))

Get the face lattice as Finite Poset::

    sage: P = polytopes.permutahedron(7)
    sage: C = CombinatorialPolyhedron(P)
    sage: C.face_lattice()
    Finite lattice containing 47294 elements

One can simplify the labels of the vertices in the face lattice::

    sage: P = polytopes.permutahedron(6)
    sage: %time CombinatorialPolyhedron(P).face_lattice(vertices=False,facets=False)
    CPU times: user 400 ms, sys: 20 ms, total: 420 ms
    Wall time: 380 ms
    Finite lattice containing 4684 elements
    sage: %time CombinatorialPolyhedron(P).face_lattice(vertices=True,facets=True)
    CPU times: user 2.48 s, sys: 60 ms, total: 2.54 s
    Wall time: 2.4 s
    Finite lattice containing 4684 elements
    sage: %time CombinatorialPolyhedron(P).face_lattice(vertices=True,facets=False) #default
    CPU times: user 1.87 s, sys: 48 ms, total: 1.92 s
    Wall time: 1.76 s
    Finite lattice containing 4684 elements

Calculate the entries of the flag-vector::

    sage: P = polytopes.permutahedron(7)
    sage: C = CombinatorialPolyhedron(P)
    sage: C.flag(1,3)
    151200
    sage: C.flag(5,2)
    67200
    sage: C.flag(2,4,6)
    100800


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
from sage.geometry.polyhedron.base import is_Polyhedron
from sage.geometry.lattice_polytope import is_LatticePolytope

from .hasse_diagram cimport   CombinatorialPolyhedron_ptr, init_CombinatorialPolyhedron, dimension, edges, f_vector, ridges, incidences, record_all_faces, get_faces, get_flag, delete_CombinatorialPolyhedron, get_maxnumberedges, get_maxnumberincidences

from cpython cimport array
import array
from cpython.mem cimport PyMem_Malloc, PyMem_Free

#TODO take care of the empty polyhedron, which does not have vertices
cdef class CombinatorialPolyhedron:
    cdef CombinatorialPolyhedron_ptr _C
    cdef tuple _V
    cdef tuple _H
    cdef tuple _equalities
    cdef dict _Hinv
    cdef dict _Vinv
    cdef int is_empty
    cdef unsigned int _length_Hrep
    cdef unsigned int _length_Vrep
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
    def __init__(self, data, vertices=None, facets=None, is_unbounded=False, nr_lines=0):
        self.is_empty = 0#TODO full-dimensional polyhedron
        if nr_lines:
            is_unbounded = True
        self._equalities = ()
        if is_Polyhedron(data):
            if data.is_empty():
                self.is_empty = 1
                return
            vertices = data.Vrepresentation()
            facets = tuple(inequality for inequality in data.Hrepresentation() if inequality.is_inequality())
            self._equalities = tuple(inequality for inequality in data.Hrepresentation() if not inequality.is_inequality())
            is_unbounded = not data.is_compact()
            nr_lines = data.n_lines()
            data = data.incidence_matrix()
            self._length_Hrep = data.ncols()
            self._length_Vrep = data.nrows()
        if is_LatticePolytope(data):
            if data.npoints() == 1:
                self.is_empty = 1
                return
            vertices = data.vertices()
            self._length_Vrep = len(vertices)
            facets = data.facets()
            self._length_Hrep = len(facets)
            data = tuple(tuple(vert for vert in facet.vertices()) for facet in facets)
        if vertices:
            self._V    = tuple(vertices)
            self._Vinv = { v : i for i,v in enumerate(self._V) }
        else:
            self._V = None
            self._Vinv = None
        if facets:
            self._H    = tuple(facets)
            self._Hinv = { v : i for i,v in enumerate(self._H) }
        else:
            self._H = None
            self._Hinv = None
        if hasattr(data,"incidence_matrix"):#TODO: Better check for incidence_matrix
            data = data.incidence_matrix()
        if hasattr(data,"nrows"):#TODO: Better check for matrix
            self._length_Hrep = data.ncols()
            self._length_Vrep = data.nrows()
            rg = range(data.nrows())
            tup =  tuple(tuple(data[i,j] for i in rg) for j in range(data.ncols()) if not all(data[i,j] for i in rg))#transpose and get rid of trivial inequalites (which all vertices satisfie)
            if tup == ((),):
                self.is_empty = 1
                return
            self._C = init_CombinatorialPolyhedron(tup, int(is_unbounded), int (nr_lines))
        else:
            if self._V is None:
                if len(data) == 0:
                    self.is_empty = 1
                    return
                vertices = sorted(set.union(*map(set,data)))
                nr_vertices = len(vertices)
                if vertices != range(len(vertices)):
                    self._V    = tuple(vertices)
                    self._Vinv = { v : i for i,v in enumerate(self._V) }
            else:
                nr_vertices = len(self._V)
            self._length_Vrep = nr_vertices
            if not self._V is None:
                f = lambda v: self._Vinv[v]
            else:
                f = lambda v: int(v)
            facets = tuple(tuple(f(i) for i in j) for j in data)
            self._length_Hrep = len(facets)
            self._C = init_CombinatorialPolyhedron(facets, nr_vertices, int(is_unbounded), int (nr_lines))
            
    def __dealloc__(self):
        r"""
        This function deallocates all the memomory used by the underlying C++-class
        """
        if self.is_empty == 1:
            return
        delete_CombinatorialPolyhedron(self._C)
    
    def __repr__(self):
        if self.is_empty == 1:
            return "The Combinatorial Type of the empty Polyhedron"
        return "The Combinatorial Type of a Polyhedron of dimension %s with %s vertices"%(self.dimension(),len(self.vertices()))
    
    def vertices(self):
        if self.is_empty == 1:
            return ()
        if self._V is not None:
            return self._V
        else:
            return tuple(i[0] for i in self.faces(0))
    
    def edges(self):
        r"""
        Calculates the edges of the CombinatorialPolyhedron, i.e. the rank 2 faces.
        
        NOTE: If you want to compute edges and f_vector it is recommended to compute edges first.
        """
        if self.is_empty == 1:
            return ()
        cdef unsigned int ** edgepointer = edges(self._C) 
        cdef unsigned long maxnumberedges = get_maxnumberedges()
        #the edges are being saved in a list basically with the first entry the first vertex of the first edges, the second entry the second vertex of that edge
        #as there might be many edges they are saved in an array of arrays, with each array containing maxnumberedges of edges
        if self.dimension() <= 0:
            return ()
        cdef unsigned long nr_edges = self.f_vector()[2]
        if nr_edges > maxnumberedges*maxnumberedges:
            raise ValueError("Cannot calculate %s edges"%nr_edges)
        if self._V is not None:
            f = lambda i : self._V[i]
        else:
            f = lambda i : Integer(i)
        vertex_one = lambda i : f(edgepointer[i / maxnumberedges][2*(i % maxnumberedges)])
        vertex_two = lambda i : f(edgepointer[i / maxnumberedges][2*(i % maxnumberedges)+1])
        return tuple((vertex_one(i), vertex_two(i)) for i in range(nr_edges))
        
    def edge_graph(self):
        return Graph(self.edges(),format="list_of_edges")
        
    def dimension(self):
        if self.is_empty == 1:
            return -1
        return Integer(dimension(self._C))
        
    def ridges(self, add_equalities=False):
        r"""
        Calculates the ridges of the CombinatorialPolyhedron, i.e. the rank 2 faces. Those are given as tuples of facets.
        
        E.g. a ridge (1,2) corresponds to the meet of facet[1] and facet[2].
        
        NOTE: If you want to compute ridges and f_vector it is recommended to compute ridges first.
        """
        if self.is_empty == 1:
            return ()
        cdef unsigned int ** ridgepointer = ridges(self._C) 
        cdef unsigned long maxnumberedges = get_maxnumberedges()
        #the ridges are being saved in a list basically, with the first entry the first facet of the first ridge, the second entry the second facet of that ridges
        #as there might be many ridges they are saved in an array of arrays, with each array containing maxnumberedges of ridges
        if self.dimension() <= 0:
            return ()
        cdef unsigned long nr_ridges = self.f_vector()[-3]
        if nr_ridges > maxnumberedges*maxnumberedges:
            raise ValueError("Cannot calculate %s ridges"%nr_ridges)
        if self._H is not None:
            f = lambda i : self._H[i]
        else:
            f = lambda i : Integer(i)
        facet_one = lambda i : f(ridgepointer[i / maxnumberedges][2*(i % maxnumberedges)])
        facet_two = lambda i : f(ridgepointer[i / maxnumberedges][2*(i % maxnumberedges)+1])
        if add_equalities:
            return tuple(((self._equalities + (facet_one(i),)),(self._equalities + (facet_two(i),))) for i in range(nr_ridges))
        else:
            return tuple((facet_one(i),facet_two(i)) for i in range(nr_ridges))
        
    def ridge_graph(self):
        return Graph(self.ridges(),format="list_of_edges")
        
    def f_vector(self):
        if self.is_empty == 1:
            return (1,)
        r"""
        Calculates the f_vector of the CombinatorialPolyhedron, i.e. the vector containing the nr of faces of each rank.
        
        NOTE: If you also want to compute edges or ridges, it is recommended to do that first.
        """
        cdef array.array vector = array.array('L', [0 for _ in range(self.dimension()+2)])
        f_vector(self._C, vector.data.as_ulongs)
        return tuple(Integer(i) for i in vector)
        
    def _record_all_faces(self):
        if self.is_empty == 1:
            return 
        r"""
        Records all faces of the Polyhedron, such that you can quickly get them with self.faces(dimension)
        """
        record_all_faces(self._C)
        
    def faces(self, dimension, facet_repr=False):#TODO fix the cpp-function for dimension 0,1,self.dimension -1, self.dimension
        r"""
        Gets all k-faces for specified dimenion k.

        By default faces are given as tuple of vertices, but one may also choose tuple of facets instead.
        """
        
        if self.is_empty == 1:
            if dimension == -1:
                return ((),)
            else:
                return ()
        dim = self.dimension()
        if not dimension in range(-1,dim + 1):
            return ()
        cdef unsigned long number_of_faces = self.f_vector()[dimension + 1]
        cdef unsigned int length_of_face
        if (facet_repr == False):
            size_of_face = self._length_Vrep
        else:
            size_of_face = self._length_Hrep
        cdef unsigned int **faces_to_return = <unsigned int**> PyMem_Malloc(number_of_faces * sizeof(unsigned int *))
        for i in range(number_of_faces):
            faces_to_return[i] = <unsigned int*> PyMem_Malloc(size_of_face * sizeof(unsigned int))
        cdef unsigned int *length_of_faces = <unsigned int*> PyMem_Malloc(number_of_faces * sizeof(unsigned int))
        addtuple = ()
        
        if (facet_repr):
            facet_repr = 1
            if self._H is not None:
                f = lambda i : self._H[i]
            else:
                f = lambda i : Integer(i)
            addtuple = self._equalities
        else:
            facet_repr = 0
            if self._V is not None:
                f = lambda i : self._V[i]
            else:
                f = lambda i : Integer(i)
        get_faces(self._C, dimension, facet_repr, faces_to_return, length_of_faces)
        returntuple = tuple(addtuple + tuple(f(faces_to_return[facecounter][j]) for j in range(length_of_faces[facecounter])) for facecounter in range(number_of_faces))
        for i in range(number_of_faces):
            PyMem_Free(faces_to_return[i])
        PyMem_Free(faces_to_return)
        PyMem_Free(length_of_faces)
        return returntuple
        
    def incidences(self,dimension_one,dimension_two):
        r"""
        Gets a tuple of all incidens between faces of dimension dimension_one and dimension_two.

        Incidences are given as tuple of integers, where the integer corresponds to the order according to self.faces(dimension)
        """
        if self.is_empty == 1:
            if dimension_one == dimension_two == -1:
                return ((),())
            else:
                return ()
        cdef unsigned long maxnumberincidences = get_maxnumberincidences()
        cdef unsigned long nr_incidences[1]
        nr_incidences[:] = [0]
        cdef unsigned int twisted[1]
        twisted[:] = [0]
        cdef unsigned long **incidencepointer = incidences(self._C,dimension_one,dimension_two, nr_incidences, twisted)
        if nr_incidences[0] > maxnumberincidences*maxnumberincidences:
            raise ValueError("Cannot calculate %s incidences"%nr_incidences[0])
        if twisted[0] == 0:
            incidence_one = lambda i : Integer(incidencepointer[i / maxnumberincidences][2*(i % maxnumberincidences)])
            incidence_two = lambda i : Integer(incidencepointer[i / maxnumberincidences][2*(i % maxnumberincidences)+1])
        else:
            incidence_two = lambda i : Integer(incidencepointer[i / maxnumberincidences][2*(i % maxnumberincidences)])
            incidence_one = lambda i : Integer(incidencepointer[i / maxnumberincidences][2*(i % maxnumberincidences)+1])
        return tuple((incidence_one(i), incidence_two(i)) for i in range(nr_incidences[0]))
    def face_lattice(self,vertices=True,facets=False):
        r"""
        Generates the face-lattice.

        The faces in the face-lattice will be labeled by default by tuples of vertices. They can also be labeled by tuple of facets or both or as 0, 1, 2, ... where the order corresponds to the order of self.faces[0],self.faces[1],self.faces[2],...
        """

        if self.is_empty == 1:
            D = DiGraph([((),),()], format='vertices_and_edges')
            return FiniteLatticePoset(D);
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
        
    def flag(self,*flag):
        r"""
        Return counts the elements in that entry of the flag vector.
        
        flag(i) is equivalent to f_vector(i)
        
        flag(i_1,...,i_n) will count the number of tuples (face_1,...,face_n), where face_j is and i_j face and face_1 \subset face_2 ...
        
        The implementation sorts the i_1,...i_n.
        """
        for number in flag:
            if not isinstance(number,Integer):
                return TypeError("All arguments of combinatorialPolyhedron.flag() must be integers.")
        dim = self.dimension()
        flag = set(number for number in flag if number in range(-1,dim+1))
        if flag == set():
            return 0
        if flag  <= set([-1,dim]):
            return 1
        if self.is_empty == 1:
            return 0
        cdef array.array flagarray = array.array('I', sorted(number for number in flag if number in range(0,dim)))
        return Integer(get_flag(self._C, flagarray.data.as_uints, len(flagarray)))
