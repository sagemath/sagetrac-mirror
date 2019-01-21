r"""
Several algorithms working implicitly with the hasse_diagram (of a polytope), including calculating the f_vector, dimension, flags, level-sets and even the face lattice.

This is a wrapper for the functions in hasse_diagram.cc.

This computes implicitely a finite atomic and coatomic lattices, where every interval of length two has at least 4 elements.
(Exactly 4 is known as the diamond property).
In particular this module calculates quickly the f_vector of polytopes. The input must be a tuple of coatoms given each by a tuple of atoms. The atoms must be labeled 0,...,n.


AUTHOR:

- Jonathan Kliem (2019-01)


"""


#*****************************************************************************
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
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

from .hasse_diagram cimport   CombinatorialPolyhedron_ptr, init_CombinatorialPolyhedron, dimension, edges, f_vector, ridges, incidences, record_all_faces, get_faces, face_iterator, face_iterator_init, get_flag, delete_CombinatorialPolyhedron, get_maxnumberedges, get_maxnumberincidences

from cpython cimport array
import array
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from sage.structure.sage_object cimport SageObject

cdef class CombinatorialPolyhedron(SageObject):
    r"""
    The class of the Combinatorial Type of a Polyehdron, a Polytope.

    INPUT:

    - ``data`` -- a ``Polyhedron``, i.e. an instance of
      :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`.

    or

    - ``data`` -- a ``LatticePolytope``, i.e. an instance of
      :class:`~sage.geometry.lattice_polytope.LatticePolytopeClass`.

    or

    - ``data`` -- an ``incidence_matrix`` as in
      :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`
      of :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`.

      * ``vertices`` -- a list of ``[vertices, rays, lines]``, if
        the rows in the incidence_matrix should correspond to names.

      * ``facets`` -- a list of facets, if
        the columns in the incidence_matrix should correspond to names.

      * ``nr_lines`` -- for bounded Polyhedra, this should be the 
      default: ``None``. For unbounded Polyhedra, this needs to be set
      to the correct nr of lines, i.e. the the maximum nr of lines with
      linearly independent directions in the Polyehdron.

    or

    - ``data`` -- a list of facets,
      each facet given as a list of ``[vertices, rays, lines]``.
      If the Polyhedron is unbounded, then rays and lines are required.
      If the Polyehdron contains no lines the rays can be thought of as
      the vertices of the facets deleted from a bounded Polyhedron. See
      :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`
      on how to use rays and lines.

      * ``facets`` -- a list of names of the facets, if
        the facets given should correspond to names.

      * ``unbounded`` -- for bounded Polyhedra, this should be the 
      default ``None``. For unbounded Polyhedra, this needs to be set
      to the correct nr of lines, i.e. the the maximum nr of lines with
      linearly independent directions in the Polyehdron.


    EXAMPLES:

    Input is Polyhedron::

        sage: P = polytopes.cube()
        sage: CombinatorialPolyhedron(P)
        Combinatorial Type of a Polyhedron of dimension 3 with 8 vertices

    Input is a LatticePolytope::

        sage: points = [(1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1)]
        sage: L = LatticePolytope(points)
        sage: CombinatorialPolyhedron(L)
        Combinatorial Type of a Polyhedron of dimension 3 with 6 vertices

    Input is an incidence matrix::

        sage: data = Polyhedron(rays=[[0,1]]).incidence_matrix()
        sage: CombinatorialPolyhedron(data, nr_lines=0)
        Combinatorial Type of a half-space of dimension 1
        sage: C = CombinatorialPolyhedron(data, vertices=['myvertex'], facets=['myfacet'], nr_lines=0)
        sage: C.Vrepresentation()
        ('myvertex',)
        sage: C.Hrepresentation()
        ('myfacet',)

    You can also give the facets explicitely::

        sage: CombinatorialPolyhedron(((1,2,3),(1,2,4),(1,3,4),(2,3,4)))
        Combinatorial Type of a Polyhedron of dimension 3 with 4 vertices
        sage: facetnames = ['facet0', 'facet1', 'facet2', 'myfacet3']
        sage: facetinc = ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        sage: C = CombinatorialPolyhedron(facetinc, facets=facetnames)
        sage: C.Vrepresentation()
        (1, 2, 3, 4)
        sage: C.Hrepresentation()
        ('facet0', 'facet1', 'facet2', 'myfacet3')
        
    Specifying the nr of lines is important::
    
        sage: P = Polyhedron(ieqs=[[1,-1,0],[1,1,0]])
        sage: C = CombinatorialPolyhedron(P) #this works fine
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices
        sage: data = P.incidence_matrix()
        sage: C = CombinatorialPolyhedron(data) #wrong usage!
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 1 vertices
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C = CombinatorialPolyhedron(data, nr_lines=1) #correct usage!
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices
        sage: C.f_vector()
        (1, 0, 2, 1)
        

    
    """
    cdef CombinatorialPolyhedron_ptr _C
    cdef tuple _V
    cdef tuple _H
    cdef tuple _equalities
    cdef dict _Hinv
    cdef dict _Vinv
    cdef int is_trivial #in some instances the polyhedron might not have facets or otherwise produce errors in the C function
    cdef int _dimension #in the case of is_trivial we will manually tell the dimension
    cdef unsigned int _length_Hrep
    cdef unsigned int _length_Vrep
    cdef int _unbounded #this is 0, if the Polyhedron is unbounded, otherwise it is 1 + nr_lines
    def __init__(self, data, vertices=None, facets=None, nr_lines=None):
        r"""
        Initializes the combinatorial polyhedron.

        See :class:`CombinatorialPolyhedron` for a description of the input
        data.

        TESTS::

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],[0,2,3],[1,2,3]])    # indirect doctests
        """
        cdef unsigned int **incidence_matrix
        cdef unsigned int **facets_pointer
        cdef unsigned int *len_facets
        if nr_lines is None:
            self._unbounded = 0
        else:
            self._unbounded = 1 + int(nr_lines)
        self.is_trivial = 0
        self._equalities = ()
        if is_Polyhedron(data):
            if data.is_empty():
                self.is_trivial = 1
                self._dimension = -1
                return
            vertices = data.Vrepresentation()
            facets = tuple(inequality for inequality in data.Hrepresentation())
            if (len(vertices) == data.n_lines() + 1) and (data.n_lines > 0):#in this case the Polyhedron does not have facets
                self.is_trivial = 1
                self._dimension = data.n_lines()
                self._V = tuple(vertices)
                return
            self._unbounded
            if not data.is_compact():
                self._unbounded = 1 + data.n_lines()
            data = data.incidence_matrix()
            self._length_Hrep = data.ncols()
            self._length_Vrep = data.nrows()
        if is_LatticePolytope(data):
            if data.npoints() == 0:
                self.is_trivial = 1
                self._dimension = -1
                return
            if data.npoints() == 1:
                self.is_trivial = 1
                self._dimension = 0
                self._V = data.vertices()
                return
            self._unbounded = 0
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
            facets = tuple(facets)
            test = [0] * len(facets)
            for i in range(len(facets)):
                test[i] = 1
                if hasattr(facets[i], "is_inequality"):
                    if not facets[i].is_inequality():
                        test[i] = 0
            self._H    = tuple(facets[i] for i in range(len(facets)) if test[i]) #only keeping those that are actual inequalities
            self._equalities = tuple(facets[i] for i in range(len(facets)) if not test[i]) #the inequalities are saved here
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
            tup =  tuple(tuple(data[i,j] for i in rg) for j in range(data.ncols()) if not all(data[i,j] for i in rg))#transpose and get rid of equalities (which all vertices satisfie)
            if len(tup) == 0:#the case of the empty Polyhedron
                self.is_trivial = 1
                self._dimension = -1  + data.nrows()#the elements in the Vrep are assumed to be one vertex and otherwise lines
                return
            if len(tup) == 1:#the case of a half space
                self.is_trivial = 2
                self._dimension = -1 + data.nrows()#the elements in the Vrep are assumed to be one vertex, one ray and otherwise lines
                return
            incidence_matrix = <unsigned int**> PyMem_Malloc(len(tup) * sizeof(unsigned int *))
            for i in range(len(tup)):
                incidence_matrix[i] = <unsigned int*> PyMem_Malloc(self._length_Vrep * sizeof(unsigned int))
                for j in range(self._length_Vrep):
                    incidence_matrix[i][j] = tup[i][j]
            self._C = init_CombinatorialPolyhedron(incidence_matrix, len(tup), self._length_Vrep, self._unbounded)
            for i in range(len(tup)):
                PyMem_Free(incidence_matrix[i])
            PyMem_Free(incidence_matrix)
        elif isinstance(data,Integer):
            if data < -1:
                TypeError("A polyhedron must have dimension at least -1")
            self.is_trivial = 1
            self._dimension = data
        else:#assumes the facets are given as a list of vertices/rays/lines
            if len(data) == 0:
                    self.is_trivial = 1
                    self._dimension = -1
                    return
            if len(data) == 1:
                    self.is_trivial = 2
                    self._dimension = len(data[0]) - 1
                    if self._dimension <= 0:
                        self.is_trivial = 1 #we are treating a polyhedron equal to its affine hull
                    return
            if self._V is None:
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
            facets_pointer = <unsigned int**> PyMem_Malloc(len(facets) * sizeof(unsigned int *))
            len_facets = <unsigned int*> PyMem_Malloc(len(facets) * sizeof(unsigned int))
            for i in range(len(facets)):
                len_facets[i] = len(facets[i])
                facets_pointer[i] = <unsigned int*> PyMem_Malloc(len_facets[i]  * sizeof(unsigned int))
                for j in range(len_facets[i] ):
                    facets_pointer[i][j] = facets[i][j]
            self._C = init_CombinatorialPolyhedron(facets_pointer, len(facets), len_facets, nr_vertices, self._unbounded)
            for i in range(len(facets)):
                PyMem_Free(facets_pointer[i])
            PyMem_Free(facets_pointer)
            PyMem_Free(len_facets)

    def __dealloc__(self):
        r"""
        This function deallocates all the memomory used by the underlying C++-class
        """
        if self.is_trivial > 0:
            return
        delete_CombinatorialPolyhedron(self._C)

    def _repr_(self):
        if self.is_trivial == 1:
            if self._dimension == 0:
                return "Combinatorial Type of the Polyhedron with one vertex"
            if self._dimension > 0:
                return "Combinatorial Type of a trivial Polyhedron of dimension %s"%self._dimension
            return "Combinatorial Type of the empty Polyhedron"
        if self.is_trivial == 2:
            return "Combinatorial Type of a half-space of dimension %s"%self._dimension
        return "Combinatorial Type of a Polyhedron of dimension %s with %s vertices"%(self.dimension(),len(self.vertices()))

    def __reduce__(self):
        r"""
        Override __reduce__ to correctly pickle/unpickle.

        TESTS::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: tuple(i for i in C.face_iter(facet_repr=True)) == tuple(i for i in C1.face_iter(facet_repr=True))
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: tuple(i for i in C.face_iter(facet_repr=True)) == tuple(i for i in C1.face_iter(facet_repr=True))
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: C.f_vector() == C1.f_vector()
            True
        """
        unbounded = None
        if self._unbounded:
            unbounded = Integer(unbounded) - 1
        if self.is_trivial == 1:
            return (CombinatorialPolyhedron, (Integer(self._dimension), self.Vrepresentation(), self.Hrepresentation(), unbounded))
        if self.is_trivial == 2:
            pickletuple = tuple(0 for _ in range(self._dimension + 1))
            return (CombinatorialPolyhedron, ((pickletuple,), self.Vrepresentation(), self.Hrepresentation(), unbounded))#it is important that there is exactly one tuple of the correct length
        return (CombinatorialPolyhedron, (self.facets(), self.Vrepresentation(), self.Hrepresentation(), unbounded))

    def Vrepresentation(self):
        r"""
        Return a list of names of ``[vertices, rays, lines]``.

        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.Vrepresentation()
            (A line in the direction (0, 0, 1),
             A ray in the direction (1, 0, 0),
             A vertex at (0, 0, 0),
             A ray in the direction (0, 1, 0))

        """
        if self._V is not None:
            return self._V
        else:
            return tuple(Integer(i) for i in range(self._length_Vrep))

    def Hrepresentation(self):
        r"""
        Returns a list of names of facets and possibly some equalities.

        EXAMPLES::

            sage: P = polytopes.permutahedron(3)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.Hrepresentation()
            (An equation (1, 1, 1) x - 6 == 0,
             An inequality (0, -1, -1) x + 5 >= 0,
             An inequality (0, 0, -1) x + 3 >= 0,
             An inequality (0, -1, 0) x + 3 >= 0,
             An inequality (0, 1, 0) x - 1 >= 0,
             An inequality (0, 1, 1) x - 3 >= 0,
             An inequality (0, 0, 1) x - 1 >= 0)

        """
        if self._H is not None:
            return self._equalities + self._H
        else:
            return tuple(Integer(i) for i in range(self._length_Hrep))

    def vertices(self, names=True):
        r"""
        Returns the elements in the ``Vrepresentation`` that are vertices.

        In the case of an unbounded Polyhedron, there might be lines and
        rays in the Vrepresentation.

        If ``names`` is set to ``False``, then the vertices are given by
        their indices in the Vrepresentation.

        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[0,0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.vertices()
            (A vertex at (0, 0, 0),)
            sage: C.Vrepresentation()
            (A vertex at (0, 0, 0),
             A ray in the direction (0, 0, 1),
             A ray in the direction (0, 1, 0),
             A ray in the direction (1, 0, 0))
            sage: P = polytopes.cross_polytope(3)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.vertices()
            (A vertex at (-1, 0, 0),
             A vertex at (0, -1, 0),
             A vertex at (0, 0, -1),
             A vertex at (0, 0, 1),
             A vertex at (0, 1, 0),
             A vertex at (1, 0, 0))
            sage: C.vertices(names=False)
            (0, 1, 2, 3, 4, 5)

            sage: points = [(1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1)]
            sage: L = LatticePolytope(points)
            sage: C = CombinatorialPolyhedron(L)
            sage: C.vertices()
            (M(1, 0, 0), M(0, 1, 0), M(0, 0, 1), M(-1, 0, 0), M(0, -1, 0), M(0, 0, -1))
            sage: C.vertices(names=False)
            (0, 1, 2, 3, 4, 5)

        """

        return tuple(i[0] for i in self.faces(0, names=names))

    def facets(self, names=True):
        r"""
        Returns the facets as lists of ``[vertices, rays, lines]``.

        If ``names`` is ``False``, then the vertices in the facets
        are given by their indices in the Vrepresentation.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.facets()
            ((A vertex at (-1, -1, 1),
              A vertex at (-1, 1, 1),
              A vertex at (1, -1, 1),
              A vertex at (1, 1, 1)),
             (A vertex at (-1, 1, -1),
              A vertex at (-1, 1, 1),
              A vertex at (1, 1, -1),
              A vertex at (1, 1, 1)),
             (A vertex at (1, -1, -1),
              A vertex at (1, -1, 1),
              A vertex at (1, 1, -1),
              A vertex at (1, 1, 1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, -1, 1),
              A vertex at (-1, 1, -1),
              A vertex at (-1, 1, 1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, 1, -1),
              A vertex at (1, -1, -1),
              A vertex at (1, 1, -1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, -1, 1),
              A vertex at (1, -1, -1),
              A vertex at (1, -1, 1)))
            sage: C.facets(names=False)
            ((1, 3, 5, 7),
             (2, 3, 6, 7),
             (4, 5, 6, 7),
             (0, 1, 2, 3),
             (0, 2, 4, 6),
             (0, 1, 4, 5))


        """
        return self.faces(self.dimension()-1, names=names)

    def edges(self, names=True):
        r"""
        Returns the edges of the CombinatorialPolyhedron,
        i.e. the rank 1 faces, which contain 2 vertices.

        If ``names`` is set to ``False``, then the vertices in the edges
        are given by their indices in the Vrepresentation.

        If you want to compute all faces of dimension 1,
        use :meth:`CombinatorialPolyhedron.faces` instead.

        .. NOTE::

            To compute edges and f_vector, first compute edges.

        EXAMPLES::

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
            sage: C.edges(names=False)
            ((0, 4), (1, 4), (2, 4), (3, 4), (0, 3), (1, 3), (0, 2), (1, 2), (0, 1))

            sage: P = Polyhedron(rays=[[-1,0],[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ()
            sage: C.faces(1)
            ((A line in the direction (1, 0), A vertex at (0, 0)),)
            
            sage: P = Polyhedron(vertices=[[0,0],[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ((A vertex at (0, 0), A vertex at (1, 0)),)

        """
        if self.is_trivial > 0:
            return ()
        cdef unsigned int ** edgepointer = edges(self._C)
        cdef unsigned long maxnumberedges = get_maxnumberedges()
        #the edges are being saved in a list basically with the first entry the first vertex of the first edges, the second entry the second vertex of that edge
        #as there might be many edges they are saved in an array of arrays, with each array containing maxnumberedges of edges
        if self.dimension() <= 0:
            return ()
        cdef unsigned long nr_edges = self.f_vector()[2]
        if self.f_vector()[1] < 2:#in case of a cone or a Polyhedron there are no edges
            return ()
        if self.dimension() == 1:#the C-module assumes at least dimension 2 in order to produce any edges
            return self.faces(1, names=names)
        if nr_edges > maxnumberedges*maxnumberedges:
            raise ValueError("Cannot calculate %s edges"%nr_edges)
        if self._V is not None and names==True:
            f = lambda i : self._V[i]
        else:
            f = lambda i : Integer(i)
        vertex_one = lambda i : f(edgepointer[i / maxnumberedges][2*(i % maxnumberedges)])
        vertex_two = lambda i : f(edgepointer[i / maxnumberedges][2*(i % maxnumberedges)+1])
        return tuple((vertex_one(i), vertex_two(i)) for i in range(nr_edges))

    def edge_graph(self,names=True):
        r"""
        Returns the edge graph.

        If ``names`` is set to ``False``, the vertices will carry names
        according to the indexing of the Vrepresentation.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(3,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edge_graph()
            Graph on 5 vertices
            sage: G = C.edge_graph()
            sage: G.degree()
            [3, 4, 4, 3, 4]

        """

        return Graph(self.edges(names=names),format="list_of_edges")

    def dimension(self):
        r"""
        Returns the dimension of the ``CombinatorialPolyehdron``.

        EXAMPLES::

            sage: C = CombinatorialPolyhedron(((1,2,3),(1,2,4),(1,3,4),(2,3,4)))
            sage: C.dimension()
            3

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
            sage: CombinatorialPolyhedron(P).dimension()
            3


        """
        if self.is_trivial > 0:
            return Integer(self._dimension)
        return Integer(dimension(self._C))

    def ridges(self, add_equalities=False, names=True):
        r"""
        Returns the ridges.

        The ridges of the CombinatorialPolyhedron are the faces
        contained in exactly two facets.

        If you want to compute all faces of codimension 1,
        use :meth:`CombinatorialPolyhedron.faces` instead.

        The ridges will be given by the facets, they are contained in.

        - If ``add_equalities`` is ``True``, then equalities the entire
          Polyhedron satisfies, are added.

        - If ``names`` is ``False``, then the facets in the ridges are
          given by their indices in the Hrepresentation.

        .. NOTE::

            To compute ridges and f_vector, compute ridges first.

        EXAMPLES::

            sage: P = polytopes.permutahedron(2)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.ridges()
            ((An inequality (0, -1) x + 2 >= 0, An inequality (0, 1) x - 1 >= 0),)
            sage: C.ridges(add_equalities=True)
            (((An equation (1, 1) x - 3 == 0, An inequality (0, -1) x + 2 >= 0),
              (An equation (1, 1) x - 3 == 0, An inequality (0, 1) x - 1 >= 0)),)
            sage: P = polytopes.cyclic_polytope(4,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.ridges()
            ((An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-12, 19, -8, 1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (24, -26, 9, -1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (-12, 19, -8, 1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (-12, 19, -8, 1) x + 0 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (-12, 19, -8, 1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (8, -14, 7, -1) x + 0 >= 0))
            sage: C.ridges(names=False)
            ((0, 4),
             (1, 4),
             (2, 4),
             (3, 4),
             (0, 3),
             (1, 3),
             (2, 3),
             (0, 2),
             (1, 2),
             (0, 1))

            sage: P = Polyhedron(rays=[[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C
            Combinatorial Type of a half-space of dimension 1
            sage: C.ridges()
            ()
            sage: C.faces(0, facet_repr=True)
            ((An equation (0, 1) x + 0 == 0, An inequality (1, 0) x + 0 >= 0),)


        """

        if self.is_trivial > 1:
            return ()
        if self.f_vector()[-2] <= 1:
            return ()#in some instances there are no ridges as intersections of facets, i.e. if we treat a ray, or a halfspace
        cdef unsigned int ** ridgepointer = ridges(self._C)
        cdef unsigned long maxnumberedges = get_maxnumberedges()
        #the ridges are being saved in a list basically, with the first entry the first facet of the first ridge, the second entry the second facet of that ridges
        #as there might be many ridges they are saved in an array of arrays, with each array containing maxnumberedges of ridges
        if self.dimension() <= 0:
            return ()
        cdef unsigned long nr_ridges = self.f_vector()[-3]
        if nr_ridges > maxnumberedges*maxnumberedges:
            raise ValueError("Cannot calculate %s ridges"%nr_ridges)
        if self._H is not None and names==True:
            f = lambda i : self._H[i]
        else:
            f = lambda i : Integer(i)
        facet_one = lambda i : f(ridgepointer[i / maxnumberedges][2*(i % maxnumberedges)])
        facet_two = lambda i : f(ridgepointer[i / maxnumberedges][2*(i % maxnumberedges)+1])
        if add_equalities:
            return tuple(((self._equalities + (facet_one(i),)),(self._equalities + (facet_two(i),))) for i in range(nr_ridges))
        else:
            return tuple((facet_one(i),facet_two(i)) for i in range(nr_ridges))

    def ridge_graph(self, names=True):
        r"""
        Returns the ridge graph.

        The ridge graph of the CombinatorialPolyhedron consists of
        ridges as edges and facets as vertices.

        If ``names`` is ``False``, the ``vertices`` of the graph  will
        be the incidences of the facets in the Hrepresentation.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(4,6)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.ridge_graph()
            Graph on 9 vertices

        """
        return Graph(self.ridges(names=names),format="list_of_edges")

    def f_vector(self):
        r"""
        Calculates the ``f_vector`` of the CombinatorialPolyhedron.

        The ``f_vector`` contains the number of faces of dimension ``k``
        for each ``k`` in ``range(-1, self.dimension() + 1)``.

        .. NOTE::

            If you also want to compute edges and/or ridges, do so
            first.

        EXAMPLES::

            sage: P = polytopes.permutahedron(7)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.f_vector()
            (1, 5040, 15120, 16800, 8400, 1806, 126, 1)

            sage: P = polytopes.cyclic_polytope(6,40)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.f_vector()
            (1, 40, 780, 9880, 25940, 25200, 8400, 1)

        ALGORITHM:

        The number of facets is assumed to be at least two here.
        The algorithm to visit all proper faces exactly once is roughly
        equivalent to:

        facets = [set(facet) for facet in self.facets()]
        ComputeNextStep(facets, [])

        #this algorithm assumes at each step to receive all facets of
        #some face except those contained in a face of forbidden
        def ComputeNextStep(faces, forbidden):

            for face in faces:
                pass #do something here with that face

            while len(faces) > 1:
                one_face = faces.pop()
                newfaces = [one_face.intersection(face) for face in faces]
                #newfaces contains all intersection

                newfaces2 = []
                for face1 in newfaces:
                    #face1 is a facet of one_face iff
                    #it is not contained in another facet
                    if all(not face1 < face2 for face2 in newfaces):
                        newfaces2.append(face1)
                #newfaces2 contains all facets of one_face not contained
                #in any one of forbidden and maybe some that are
                #contained in one of forbidden

                newfaces3 = []
                for face1 in newfaces2:
                    if all(not face1 < face2 for face2 in forbidden):
                        newfaces3.append(face1)
                #newfaces3 contains exactly all facets of one_face but
                #those contained in one face of forbidden

                #visit all faces in one_face that are not contained in
                #one of forbidden
                ComputeNextStep(newfaces3, forbidden)

                #we have visited all faces in one_face, so we should not
                #visit one ever again

                forbidden.append(one_face)

            return


        """
        if self.is_trivial == 1:
            if self._dimension == -1:
                return (Integer(1),)
            return (Integer(1),) + tuple(Integer(0) for _ in range(self._dimension)) + (Integer(1),)
        if self.is_trivial == 2:
            return (Integer(1),) + tuple(Integer(0) for _ in range(self._dimension - 1)) + (Integer(1),Integer(1))
        cdef array.array vector = array.array('L', [0 for _ in range(self.dimension()+2)])
        f_vector(self._C, vector.data.as_ulongs)
        return tuple(Integer(i) for i in vector)

    def _record_all_faces(self):
        r"""
        Records all faces of the Polyhedron, such that they can quickly accessed later.
        """
        if self.is_trivial >= 1:
            return
        record_all_faces(self._C)

    def faces(self, dimension, facet_repr=False, names=True):
        r"""
        Gets all k-faces for specified dimenion k.

        By default faces are given as tuple of vertices.

        If ``facet_repr`` is set to ``True``, then vertices are given in
        as tuple of facets.

        If ``names`` is set to ``False``, then the vertices are given by
        their indices in the Vrepresentation.


        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.faces(2)
            ((A vertex at (0, 0, 0, 1),
              A vertex at (0, 0, 1, 0),
              A vertex at (0, 1, 0, 0)),
             (A vertex at (0, 0, 0, 1),
              A vertex at (0, 0, 1, 0),
              A vertex at (1, 0, 0, 0)),
             (A vertex at (0, 0, 0, 1),
              A vertex at (0, 1, 0, 0),
              A vertex at (1, 0, 0, 0)),
             (A vertex at (0, 0, 1, 0),
              A vertex at (0, 1, 0, 0),
              A vertex at (1, 0, 0, 0)))
            sage: C.faces(2, facet_repr=True)
            ((An equation (1, 1, 1, 1) x - 1 == 0,
              An inequality (0, -1, -1, -1) x + 1 >= 0),
             (An equation (1, 1, 1, 1) x - 1 == 0, An inequality (0, 1, 0, 0) x + 0 >= 0),
             (An equation (1, 1, 1, 1) x - 1 == 0, An inequality (0, 0, 1, 0) x + 0 >= 0),
             (An equation (1, 1, 1, 1) x - 1 == 0, An inequality (0, 0, 0, 1) x + 0 >= 0))

            sage: P = Polyhedron(rays=[[1,0],[0,1]], vertices=[[1,0],[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.faces(1)
            ((A ray in the direction (0, 1), A vertex at (0, 1)),
             (A ray in the direction (1, 0), A vertex at (1, 0)),
             (A vertex at (0, 1), A vertex at (1, 0)))
            sage: C.faces(1, names=False)
            ((0, 1), (2, 3), (1, 3))
            sage: C.faces(1, facet_repr=True)
            ((An inequality (1, 0) x + 0 >= 0,),
             (An inequality (0, 1) x + 0 >= 0,),
             (An inequality (1, 1) x - 1 >= 0,))

        ALGORITHM:

        See :meth:`f_vector` for a description
        on how all faces are visited.
        """
        if self.is_trivial == 1: #the case of a Polyhedron equal to its affine hull
            if dimension == -1:
                return ((),)
            if dimension == self._dimension:
                if facet_repr == True:
                    return ((),)
                if self._V is not None and names == True:
                    return (tuple(self._V),)
                else:
                    if dimension == 0:
                        return (Integer(0),)
                    return ((),)
            return ()

        if self.is_trivial == 2: #the case of a half-space
            if dimension == -1:
                if facet_repr == True:
                    if self._H is not None and names == True:
                        return (self._equalities + tuple(self._H),)
                    else:
                        return ((Integer(0),),)
                return ((),)
            if dimension == self._dimension - 1:
                if facet_repr == True:
                    if self._H is not None and names == True:
                        return (self._equalities + tuple(self._H),)
                    else:
                        return ((Integer(0),),)
                return ((),)
            if dimension == self._dimension:
                return ((),)
            return ()

        dim = self.dimension()
        if not dimension in range(-1,dim + 1):
            return ()

        #creating an array of arrays for the C function to store the faces in
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

        #filling the array
        get_faces(self._C, dimension, facet_repr, faces_to_return, length_of_faces)

        #translating the result to the desired representation
        addtuple = ()
        if (facet_repr):
            facet_repr = 1
            if self._H is not None and names == True:
                f = lambda i : self._H[i]
            else:
                f = lambda i : Integer(i)
            addtuple = self._equalities
        else:
            facet_repr = 0
            if self._V is not None and names == True:
                f = lambda i : self._V[i]
            else:
                f = lambda i : Integer(i)

        returntuple = tuple(addtuple + tuple(f(faces_to_return[facecounter][j]) for j in range(length_of_faces[facecounter])) for facecounter in range(number_of_faces))

        #cleaning up
        for i in range(number_of_faces):
            PyMem_Free(faces_to_return[i])
        PyMem_Free(faces_to_return)
        PyMem_Free(length_of_faces)

        return returntuple

    def face_iter(self, dimension=None, vertex_repr=True, facet_repr=False, names=True):
        r"""
        Iterator over all faces of specified dimension.

        If ``dimension`` is not specified then iterate over all faces.

        If ``vertex_repr`` is ``True``, then give the faces as lists of
        elements in ``Vrepresentation``.

        If ``facet_repr`` is ``True``, then give the faces as lists of
        elements in ``Hrepresentation``.

        - Both ``vertex_repr`` and ``facet_repr`` can be set to ``True,
          this will give each face as tuple of the form
          (``vertex_repr``, ``facet_rerpr``).

        If ``names`` is ``False`` then vertices and facets are labeled
        by their indexes.

        .. WARNING::

            There can only be one face iterator around. The second
            one will invalid the first one. Even worse, if you then call
            the first one again, the output will be incorrect and the
            second one will not yield all faces.

            This shares resources with other methods. Do not call other
            methods while iteration, s.t.

            - :meth:`faces`
            - :meth:`vertices`
            - :meth:`edges` (a second call is fine)
            - :meth:`ridges` (a second call is fine)
            - :meth:`f_vector` (a second call is fine)
            - :meth:`_record_all_faces`
            - :meth:`incidences`
            - :meth:`face_lattice`
            - :meth:`flag`
            - :meth:`k-simplicial`
            - :meth:`k-simple`

        EXAMPLES::

            sage: P = polytopes.permutahedron(7)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(dimension=2)
            sage: next(it)
            (A vertex at (5, 7, 1, 4, 6, 3, 2),
             A vertex at (5, 7, 2, 4, 6, 3, 1),
             A vertex at (6, 7, 1, 4, 5, 3, 2),
             A vertex at (6, 7, 2, 4, 5, 3, 1))
            sage: next(it)
            (A vertex at (5, 6, 1, 4, 7, 3, 2),
             A vertex at (5, 6, 2, 4, 7, 3, 1),
             A vertex at (5, 7, 1, 4, 6, 3, 2),
             A vertex at (5, 7, 2, 4, 6, 3, 1))
            sage: it = C.face_iter(dimension=2, names=False)
            sage: next(it)
            (3497, 3521, 4217, 4241)
            sage: next(it)
            (3377, 3401, 3497, 3521)
            sage: it = C.face_iter(dimension=2, facet_repr=True)
            sage: next(it)
            ((A vertex at (5, 7, 1, 4, 6, 3, 2),
              A vertex at (5, 7, 2, 4, 6, 3, 1),
              A vertex at (6, 7, 1, 4, 5, 3, 2),
              A vertex at (6, 7, 2, 4, 5, 3, 1)),
             (An equation (1, 1, 1, 1, 1, 1, 1) x - 28 == 0,
              An inequality (0, -1, 0, 0, 0, 0, 0) x + 7 >= 0,
              An inequality (0, 0, 1, 0, 0, 0, 1) x - 3 >= 0,
              An inequality (0, 0, 1, 0, 0, 1, 1) x - 6 >= 0,
              An inequality (0, 0, 1, 1, 0, 1, 1) x - 10 >= 0))
            sage: it = C.face_iter(dimension=2, vertex_repr=False, facet_repr=True, names=False)
            sage: next(it)
            (98, 106, 124, 125)
            sage: next(it)
            (80, 106, 124, 125)


        TESTS::

            sage: P = polytopes.permutahedron(7)
            sage: C = CombinatorialPolyhedron(P)
            sage: f = C.f_vector()
            sage: altf = tuple(len(tuple(C.face_iter(i))) for i in range(-1,C.dimension()+1))
            sage: altf == f
            True
            sage: allfaces = tuple(tuple(C.face_iter(i)) for i in range(-1,C.dimension()+1))
            sage: allfaces2 = tuple(tuple(C.face_iter(i)) for i in range(-1,C.dimension()+1))
            sage: all(sorted(sorted(j) for j in allfaces[i]) == sorted(sorted(j) for j in allfaces2[i]) for i in range(-1, C.dimension()+1))
            True

            sage: P = polytopes.cyclic_polytope(5,20)
            sage: C = CombinatorialPolyhedron(P)
            sage: f = C.f_vector()
            sage: altf = tuple(len(tuple(C.face_iter(i))) for i in range(-1,C.dimension()+1))
            sage: altf == f
            True
            sage: allfaces = tuple(tuple(C.face_iter(i)) for i in range(-1,C.dimension()+1))
            sage: C._record_all_faces()
            sage: allfaces2 = tuple(tuple(C.face_iter(i)) for i in range(-1,C.dimension()+1))
            sage: all(sorted(sorted(j) for j in allfaces[i]) == sorted(sorted(j) for j in allfaces2[i]) for i in range(-1, C.dimension()+1))
            True
            sage: C = CombinatorialPolyhedron(P)
            sage: allfaces = tuple(tuple(C.face_iter(i, vertex_repr=False, facet_repr=True) for i in range(-1, C.dimension()+1)))
            sage: C._record_all_faces()
            sage: allfaces2 = tuple(tuple(C.face_iter(i, vertex_repr=False, facet_repr=True) for i in range(-1, C.dimension()+1)))
            sage: all(sorted(sorted(j) for j in allfaces[i]) == sorted(sorted(j) for j in allfaces2[i]) for i in range(-1, C.dimension()+1))
            True



        ALGORITHM:

        See :meth:`f_vector` for a description
        on how all faces are visited.
        """

        if dimension is not None:
            dimension = int(dimension)
            dimensionrange = (dimension,)
        else:
            dimensionrange = range(-1,self.dimension()+1)
            dimension = -2

        if not facet_repr:
            vertex_repr == True

        if self.is_trivial > 0: #taking care of the trivial polynomial
            for dim in dimensionrange:
                if vertex_repr and facet_repr:
                    vert = self.faces(dim, names=names)
                    fac = self.faces(dim, facet_repr=True, names=names)
                    for i in range(len(vert)):
                        yield (vert[i],fac[i])
                elif vertex_repr:
                    vert = self.faces(dim, names=names)
                    for i in range(len(vert)):
                        yield vert[i]
                else:
                    fac = self.faces(dim, facet_repr=True, names=names)
                    for i in range(len(fac)):
                        yield fac[i]
            return

        if 0 == dimension:
            if vertex_repr and facet_repr:
                vert = self.faces(0, names=names)
                fac = self.faces(0, facet_repr=True, names=names)
                for i in range(len(vert)):
                    yield (vert[i],fac[i])
            elif vertex_repr:
                vert = self.faces(0, names=names)
                for i in range(len(vert)):
                    yield vert[i]
            else:
                fac = self.faces(0, facet_repr=True, names=names)
                for i in range(len(fac)):
                    yield fac[i]
            return

        if dimension == self.dimension() -1:
            if vertex_repr and facet_repr:
                vert = self.faces(dimension, names=names)
                fac = self.faces(dimension, facet_repr=True, names=names)
                for i in range(len(vert)):
                    yield (vert[i],fac[i])
            elif vertex_repr:
                vert = self.faces(dimension, names=names)
                for i in range(len(vert)):
                    yield vert[i]
            else:
                fac = self.faces(dimension, facet_repr=True, names=names)
                for i in range(len(fac)):
                    yield fac[i]
            return

        if -1 in dimensionrange:
            if vertex_repr and facet_repr:
                vert = self.faces(-1, names=names)
                fac = self.faces(-1, facet_repr=True, names=names)
                for i in range(len(vert)):
                    yield (vert[i],fac[i])
            elif vertex_repr:
                vert = self.faces(-1, names=names)
                for i in range(len(vert)):
                    yield vert[i]
            else:
                fac = self.faces(-1, facet_repr=True, names=names)
                for i in range(len(fac)):
                    yield fac[i]
            if -1 == dimension:
                return

        dim = self.dimension()
        if dim in dimensionrange:
            if vertex_repr and facet_repr:
                vert = self.faces(dim, names=names)
                fac = self.faces(dim, facet_repr=True, names=names)
                for i in range(len(vert)):
                    yield (vert[i],fac[i])
            elif vertex_repr:
                vert = self.faces(dim, names=names)
                for i in range(len(vert)):
                    yield vert[i]
            else:
                fac = self.faces(dim, facet_repr=True, names=names)
                for i in range(len(fac)):
                    yield fac[i]
            if dim == dimension:
                return

        dim = self.dimension()
        if not dimension in range(-2,dim + 1):
            return ()

        #creating two arrays for the C function to store the faces in
        cdef unsigned int successful
        size_of_faceV = self._length_Vrep
        size_of_faceH = self._length_Hrep
        cdef unsigned int VSize[1]
        VSize[0] = 0
        cdef unsigned int HSize[1]
        HSize[0] = 0
        cdef unsigned int *Vface_to_return = <unsigned int*> PyMem_Malloc(size_of_faceV * sizeof(unsigned int))
        cdef unsigned int *Hface_to_return = <unsigned int*> PyMem_Malloc(size_of_faceH * sizeof(unsigned int))

        #translating the result to the desired representation
        addtuple = self._equalities
        if self._H is not None and names == True:
            h = lambda i : self._H[i]
        else:
            h = lambda i : Integer(i)
        if self._V is not None and names == True:
            v = lambda i : self._V[i]
        else:
            v = lambda i : Integer(i)


        addtuple = ()
        if names:
            addtuple = self._equalities

        #init
        face_iterator_init(self._C, dimension, int(vertex_repr), int(facet_repr))

        #filling the array
        successful = face_iterator(self._C, Vface_to_return, VSize, Hface_to_return, HSize)
        while successful:
            if vertex_repr and facet_repr:
                yield (tuple(v(Vface_to_return[i]) for i in range(VSize[0])), addtuple +  tuple(h(Hface_to_return[i]) for i in range(HSize[0])))
            elif vertex_repr:
                yield tuple(v(Vface_to_return[i]) for i in range(VSize[0]))
            else:
                yield addtuple + tuple(h(Hface_to_return[i]) for i in range(HSize[0]))
            successful = face_iterator(self._C, Vface_to_return, VSize, Hface_to_return, HSize)

        #cleaning up
        PyMem_Free(Vface_to_return)
        PyMem_Free(Hface_to_return)

    def incidences(self,dimension_one, dimension_two):
        r"""
        Returns all incidences between ``dimension_one``-faces and
        ``dimension_two``-faces.

        Incidences are given as tuple of integers, where the integer
        corresponds to the enumeration according to :meth:`faces`

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.incidences(0,1)
            ((1, 0),
             (5, 0),
             (4, 1),
             (5, 1),
             (0, 2),
             (1, 2),
             (0, 3),
             (4, 3),
             (2, 4),
             (6, 4),
             (4, 5),
             (6, 5),
             (0, 6),
             (2, 6),
             (1, 7),
             (3, 7),
             (2, 8),
             (3, 8),
             (5, 9),
             (7, 9),
             (6, 10),
             (7, 10),
             (3, 11),
             (7, 11))
            sage: C.incidences(1,2)
            ((0, 0),
             (7, 0),
             (9, 0),
             (11, 0),
             (4, 1),
             (8, 1),
             (10, 1),
             (11, 1),
             (1, 2),
             (5, 2),
             (9, 2),
             (10, 2),
             (2, 3),
             (6, 3),
             (7, 3),
             (8, 3),
             (3, 4),
             (4, 4),
             (5, 4),
             (6, 4),
             (0, 5),
             (1, 5),
             (2, 5),
             (3, 5))
            sage: C.incidences(0,2)
            ((1, 0),
             (3, 0),
             (5, 0),
             (7, 0),
             (2, 1),
             (3, 1),
             (6, 1),
             (7, 1),
             (4, 2),
             (5, 2),
             (6, 2),
             (7, 2),
             (0, 3),
             (1, 3),
             (2, 3),
             (3, 3),
             (0, 4),
             (2, 4),
             (4, 4),
             (6, 4),
             (0, 5),
             (1, 5),
             (4, 5),
             (5, 5))
        """
        if self.is_trivial == 1:
            if dimension_one in (-1,self._dimension) and dimension_two in (-1,self._dimension):
                return ((0,0),)
            return ()
        if self.is_trivial == 2:
            if dimension_one in (-1,self._dimension,self._dimension -1) and dimension_two in (-1,self._dimension,self._dimension -1):
                return ((0,0),)
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

    def face_lattice(self, vertices=False, facets=False, names=False):
        r"""
        Generates the face-lattice.

        INPUT:

        - ``vertices`` -- if set to ``True`` the elements in the lattice
          will be named according to the vertices they contain.

        - ``facets`` -- if set to ``True`` the elements in the lattice
          will be named according to the facets they are contained in.

        - ``names`` -- if set to ``False``, facets and vertices will be
          named according to their index.

        Both vertices and facets can be set to ``True``.

        In the case of a trivial Polyhedron, which is equal to its own
        affine hull, ``facets`` will be set to ``False``, as the
        elements need distinct names.

        In the case of a half-space ``vertices`` and ``facets`` will be
        set to ``False``.

        OUTPUT:

        - a :class:'~sage.combinat.posets.lattices.FiniteLatticePoset'

        .. NOTE::

            As :class:'~sage.combinat.posets.lattices.FiniteLatticePoset'
            is awfully slow with elements having meaningful labels,
            the default of this function is to not do so.


        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0],[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.incidences(0,1)
            ((0, 0), (0, 1))
            sage: C.face_lattice()
            Finite lattice containing 5 elements
            sage: C.face_lattice(vertices=True).atoms()
            [(A vertex at (0, 0),)]

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: P1 = Polyhedron(rays=[[1,0], [-1,0]])
            sage: C1 = CombinatorialPolyhedron(P1)
            sage: C.face_lattice().is_isomorphic(C1.face_lattice())
            True

            sage: P = polytopes.permutahedron(7)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice()
            Finite lattice containing 47294 elements
        """

        #we must ignore part of the input to ensure an injective relabeling
        if self.is_trivial == 1:
            facets = False
        if self.is_trivial == 2:
            vertices = False
            facets = False

        f_vector = self.f_vector()
        self._record_all_faces()
        dimension = self.dimension()
        dic = {}
        range_f_vector = [k for k in range(-1, dimension + 1) if f_vector[k+1] > 0]
        range_f_vector1 = range_f_vector[1:-1]
        for k in range_f_vector:
            faces = (self.faces(k),self.faces(k,facet_repr=True))
            dic[k] = tuple((faces[0][i],faces[1][i]) for i in range(f_vector[k+1]))
        edges0 = ()
        if dimension >= 0:
            edges0 = tuple((i[0],i[1] + 1) for i in self.incidences(-1,range_f_vector[1])) #we must take care of the case, when there are no vertices
        edges = edges0 + tuple((i[0] + sum(f_vector[:k+1]),i[1] + sum(f_vector[:k+2])) for k in range_f_vector1 for i in self.incidences(k,k+1))
        all_faces = tuple(i for k in range_f_vector for i in dic[k])
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
        Returns the number of flags of given type.

        flag(i) is equivalent to f_vector(i).

        flag(i_1,...,i_n) will count the number of tuples
        (face_1,...,face_n), where each face_j is an i_j face and
        face_1 is contained in face_2 is contained in face_3 ...

        The implementation sorts the input arguments.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(8,30)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.flag(1,3,6)
            14490000

            sage: P = polytopes.permutahedron(7)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.flag(1,5)
            75600
            sage: C.flag(1,5,5)
            75600
            sage: C.flag(1,5,5,6)
            75600
            sage: C.flag(1,5,4)
            302400
            sage: C.flag(1,4,5)
            302400
            sage: C.flag(1,3,4,5)
            907200
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
        if self.is_trivial == 1:
            return 0
        if self.is_trivial == 2:
            if set(flag) <= set([-1,dim-1,dim]):
                return 1
            return 0
        cdef array.array flagarray = array.array('I', sorted(number for number in flag if number in range(0,dim)))
        return Integer(get_flag(self._C, flagarray.data.as_uints, len(flagarray)))


#ideas for later
#Error checking on intput!
#check for containments, shouldn't take long but is very nice to the user
#instead of is_unbounded, give nr_lines


#add k-simple, k-simplicial
#Example:
#sage: for i in Combinations(6,3): x.append(list(Integer(j in i) for j in range(6)))
# P = Polyhedron(vertices=x)
#this is 2-simplicial and 6-2 simple 6-1 dimensional polyhedron
#taken from lecture notes Guenter M. Ziegler
