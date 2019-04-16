# distutils: language = c++

r"""
CombinatorialPolyhedron gathers several algorithms of Polyhedra depending only
on the vertex-facet incidences.

Most importanly, this module offers a quick face iterator, which can be used
to calculate f-vector, edges, ridges and even the face lattice.

The FaceIterator would work on every atomic and coatomic lattice, where every
interval of length two has exactly 4 elements (known as the diamond property).

It also works on every lattice that is obtained by from such a lattice by
deleting all elements but the empty face contained in some of the coatoms.
Important examples are face lattices of unbounded polyhedra.

Representations in this module:

- Vertices              -- ``[vertices, rays, lines]`` of the polyhedron.
- Facets                -- facets of the polyhedron.
- Coatoms               -- the faces from which all others are constructed in
                           the algorithm. This will be facets or vertices.
                           In non-dual mode, faces are constructed as
                           intersections of the facets. In dual mode, the are
                           constructed theoretically as joins of vertices.
                           The coatoms are reprsented as incidences with the
                           atoms they contain.
- Atoms                 -- facets or vertices depending on application of algorithm.
                           Atoms are reprsented as incidences of coatoms they
                           are contained in.

- Vertex-Representation -- represents a face by a list of vertices it contains.
- Facet-Representation  -- represents a face by a list of facets it is contained in.
- Bit-Representation    -- represents incidences as ``uint64_t``-array, where
                           each Bit represents one incidences. There might
                           be trailing zeros, to fit alignment-requirements.
                           In most instances, faces are represented by the
                           Bit-representation, where each bit corresponds to
                           an atom.

AUTHOR:

- Jonathan Kliem (2019-03)
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

from __future__ import absolute_import, division, print_function
from sage.rings.integer import Integer
from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from sage.combinat.posets.lattices import FiniteLatticePoset
from sage.geometry.polyhedron.base import is_Polyhedron
from sage.geometry.lattice_polytope import is_LatticePolytope
from sage.structure.element import is_Matrix
from sage.misc.misc import is_iterator
from .list_of_faces \
        import incidence_matrix_to_bit_repr_of_facets, \
               incidence_matrix_to_bit_repr_of_vertices, \
               facets_tuple_to_bit_repr_of_facets, \
               facets_tuple_to_bit_repr_of_vertices

from sage.rings.integer cimport smallInteger
from cysignals.signals cimport sig_check, sig_block, sig_unblock

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef class CombinatorialPolyhedron(SageObject):
    r"""
    The class of the Combinatorial Type of a Polyehdron, a Polytope.

    INPUT:

    - ``data`` -- an instance of
      :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`

    or

    - ``data`` --  an instance of
      :class:`~sage.geometry.lattice_polytope.LatticePolytopeClass`

    or

    - ``data`` -- an ``incidence_matrix`` as in
      :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`

      * ``vertices`` -- a list of ``[vertices, rays, lines]``, if
        the rows in the incidence_matrix should correspond to names

      * ``facets`` -- a list of facets, if
        the columns in the incidence_matrix should correspond to names

      * ``nr_lines`` -- ``None`` for bounded Polyhedra,
        for unbounded Polyhedra, this needs to be set
        to the correct number of lines,
        i.e. the the maximum number of lines with
        linearly independent directions in the Polyehdron

    or

    - ``data`` -- a ``[list, tuple, iterator]`` of facets,
      each facet given as a list of ``[vertices, rays, lines]``
      if the Polyhedron is unbounded, then rays and lines are required
      if the Polyehdron contains no lines, the rays can be thought of as
      the vertices of the facets deleted from a bounded Polyhedron
      see :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`
      on how to use rays and lines

      * ``facets`` -- a list of names of the facets, if
        the facets given should correspond to names

      * ``nr_lines`` -- ``None`` for bounded Polyhedra,
        for unbounded Polyhedra, this needs to be set
        to the correct number of lines,
        i.e. the the maximum number of lines with
        linearly independent directions in the Polyehdron

    or

    - ``data`` -- an Integer, representing the dimension of a Polyhedron equal
      to its affine hull

    EXAMPLES:

    Input is Polyhedron::

        sage: P = polytopes.cube()
        sage: CombinatorialPolyhedron(P)
        Combinatorial Type of a Polyhedron of dimension 3 with 8 vertices

    Input is a LatticePolytope::

        sage: points = [(1,0,0), (0,1,0), (0,0,1),
        ....: (-1,0,0), (0,-1,0), (0,0,-1)]
        sage: L = LatticePolytope(points)
        sage: CombinatorialPolyhedron(L)
        Combinatorial Type of a Polyhedron of dimension 3 with 6 vertices

    Input is an incidence matrix::

        sage: data = Polyhedron(rays=[[0,1]]).incidence_matrix()
        sage: CombinatorialPolyhedron(data, nr_lines=0)
        Combinatorial Type of a Polyhedron of dimension 1 with 1 vertices
        sage: C = CombinatorialPolyhedron(data, vertices=['myvertex'],
        ....: facets=['myfacet'], nr_lines=0)
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

    Input is an integer::

        sage: CombinatorialPolyhedron(-1).f_vector()
        (1,)
        sage: CombinatorialPolyhedron(0).f_vector()
        (1, 1)
        sage: CombinatorialPolyhedron(5).f_vector()
        (1, 0, 0, 0, 0, 0, 1)

    Specifying the number of lines is important::

        sage: P = Polyhedron(ieqs=[[1,-1,0],[1,1,0]])
        sage: C = CombinatorialPolyhedron(P) #this works fine
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices
        sage: data = P.incidence_matrix()
        sage: vert = P.Vrepresentation()

    Incorrect due to missing number of lines::

        sage: C = CombinatorialPolyhedron(data, vertices=vert)
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 3 vertices
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C.vertices()
        (A line in the direction (0, 1),
         A line in the direction (0, 1),
         A line in the direction (0, 1))

    Correct usage with number of lines specified::

        sage: C = CombinatorialPolyhedron(data, vertices=vert, nr_lines=1)
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices
        sage: C.f_vector()
        (1, 0, 2, 1)
        sage: C.vertices()
        ()

    Initialization from Polyhedron will automatically specify number of lines::

        sage: P = Polyhedron(rays=[[1,0],[0,1]])
        sage: C = CombinatorialPolyhedron(P) # this works fine
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 1 vertices
        sage: data = P.incidence_matrix()
        sage: vert = P.Vrepresentation()

    Incorrect due to missing number of lines::

        sage: C = CombinatorialPolyhedron(data, vertices=vert)
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 3 vertices
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C.vertices()
        (A vertex at (0, 0), A vertex at (0, 0), A vertex at (0, 0))

    Correct usage with number of lines specified::

        sage: C = CombinatorialPolyhedron(data, vertices=vert, nr_lines=0)
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 1 vertices
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C.vertices()
        (A vertex at (0, 0),)
    """
    def __init__(self, data, vertices=None, facets=None, nr_lines=None):
        r"""
        Initialize :class:`CombinatorialPolyhedron`.

        See :class:`CombinatorialPolyhedron`.

        TESTS::

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],
            ....: [0,2,3],[1,2,3]])    # indirect doctests

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron).run()
        """
        self._dimension = -2  # a "NULL" value
        self._edges = NULL
        self._ridges = NULL
        self._face_lattice_incidences = NULL
        self._equalities = ()
        self._all_faces = None
        self._mem_tuple = ()

        # ``_length_edges_list`` should not be touched in an instance
        # of :class:`CombinatorialPolyhedron`. This number can be altered,
        # but should probably be a power of `2` (for memory usage).
        # ``_length_edges_list`` shouldn't be too small for speed and
        # shouldn't be too large, as ``ridges``, ``edges`` and ``incidences``
        # each have a memory overhead of
        # ``self._length_edges_list*2*sizeof(size_t *)``.
        self._length_edges_list = 16348

        if is_Polyhedron(data):
            # Input is ``Polyhedron``.
            vertices = data.Vrepresentation()
            facets = tuple(inequality for inequality in data.Hrepresentation())

            if not data.is_compact():
                self._unbounded = True
                self._nr_lines = int(data.n_lines())
            else:
                self._unbounded = False
                self._nr_lines = 0

            data = data.incidence_matrix()
        elif is_LatticePolytope(data):
            # Input is ``LatticePolytope``.
            self._unbounded = False
            self._nr_lines = 0
            vertices = data.vertices()
            self._length_Vrep = len(vertices)
            facets = data.facets()
            self._length_Hrep = len(facets)
            data = tuple(tuple(vert for vert in facet.vertices())
                         for facet in facets)
        else:
            # Input is different from ``Polyhedron`` and ``LatticePolytope``.
            if nr_lines is None:
                # According to input, the Polyhedron is bounded.
                self._unbounded = False
                self._nr_lines = 0
            else:
                # According to input, the Polyhedron is unbounded.
                self._unbounded = True
                self._nr_lines = int(nr_lines)

        if vertices:
            # The vertices have names, which the user might want later.
            self._V = tuple(vertices)
            self._Vinv = {v: i for i,v in enumerate(self._V)}
        else:
            self._V = None
            self._Vinv = None

        if facets:
            # The facets have names, which the user might want later.
            facets = tuple(facets)

            # Remove equalities from facets.
            test = [1] * len(facets)  # 0 if that facet is an equality
            for i in range(len(facets)):
                if hasattr(facets[i], "is_inequality"):
                    # At the moment this test only works for input being
                    if not facets[i].is_inequality():
                        test[i] = 0
            self._H = tuple(facets[i] for i in range(len(facets)) if test[i])

            # Saving the equalities.
            self._equalities = tuple(facets[i] for i in range(len(facets)) if not test[i])
        else:
            self._H = None

        if is_Matrix(data):
            # Input is incidence-matrix or was converted to it.
            self._length_Hrep = data.ncols()
            self._length_Vrep = data.nrows()

            # Initializing the facets in their Bit-representation.
            self.bitrep_facets = incidence_matrix_to_bit_repr_of_facets(data)

            # Initializing the vertices as their Bit-representation.
            self.bitrep_vertices = incidence_matrix_to_bit_repr_of_vertices(data)

            self._nr_facets = self.bitrep_facets.nr_faces

        elif isinstance(data, Integer):
            # To construct a trivial Polyhedron, equal to its affine hull,
            # one can give an Integer as Input.
            if data < -1:
                ValueError("any Polyhedron must have dimension at least -1")
            self._nr_facets = 0
            self._dimension = data

            # Initializing the facets in their Bit-representation.
            self.bitrep_facets = facets_tuple_to_bit_repr_of_facets((), 0)

            # Initializing the vertices as their Bit-representation.
            self.bitrep_vertices = facets_tuple_to_bit_repr_of_vertices((), 0)

        else:
            # Input is a "list" of facets.
            # The facets given by its ``[vertices, rays, lines]``.
            # Actually at least tuple, list, iterator will work.
            if is_iterator(data):
                data = tuple(data)

            if self._V is None:
                # Get the names of the vertices.
                vertices = sorted(set.union(*map(set, data)))
                nr_vertices = len(vertices)
                if vertices != range(len(vertices)):
                    self._V = tuple(vertices)
                    self._Vinv = {v: i for i,v in enumerate(self._V)}
            else:
                # Assuming the user gave as correct names for the vertices
                # and labeled them instead by `0,...,n`.
                nr_vertices = len(self._V)

            self._length_Vrep = nr_vertices

            # Relabel the vertices to be `0,...,n`.
            if self._V is not None:
                def f(v): return self._Vinv[v]
            else:
                def f(v): return int(v)
            facets = tuple(tuple(f(i) for i in j) for j in data)

            self._nr_facets = len(facets)
            self._length_Hrep = len(facets)

            # Initializing the facets in their Bit-representation.
            self.bitrep_facets = facets_tuple_to_bit_repr_of_facets(facets, nr_vertices)

            # Initializing the vertices as their Bit-representation.
            self.bitrep_vertices = facets_tuple_to_bit_repr_of_vertices(facets, nr_vertices)

    def _repr_(self):
        r"""
        Return a description of the Combinatorial Polyhedron.

        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a Polyhedron of dimension 3 with 4 vertices'

            sage: P = Polyhedron(vertices=[])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a Polyhedron of dimension -1 with 0 vertices'

            sage: P = Polyhedron(vertices=[[0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a Polyhedron of dimension 0 with 1 vertices'

            sage: P = Polyhedron(lines=[[0,0,1],[0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices'

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[-1,0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices'
        """
        return "Combinatorial Type of a Polyhedron of "\
               "dimension %s with %s vertices" \
               % (self.dimension(), self.nr_vertices())

    def __reduce__(self):
        r"""
        Override __reduce__ to correctly pickle/unpickle.

        TESTS::

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_iter()
            sage: it1 = C1.face_iter()
            sage: tup = tuple((it.vertex_repr(), it.facet_repr()) for _ in it)
            sage: tup1 = tuple((it1.vertex_repr(), it1.facet_repr()) for _ in it1)
            sage: tup == tup1
            True

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_iter()
            sage: it1 = C1.face_iter()
            sage: tup = tuple((it.vertex_repr(), it.facet_repr()) for _ in it)
            sage: tup1 = tuple((it1.vertex_repr(), it1.facet_repr()) for _ in it1)
            sage: tup == tup1
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_iter()
            sage: it1 = C1.face_iter()
            sage: tup = tuple((it.vertex_repr(), it.facet_repr()) for _ in it)
            sage: tup1 = tuple((it1.vertex_repr(), it1.facet_repr()) for _ in it1)
            sage: tup == tup1
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0],
            ....:                      [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_iter()
            sage: it1 = C1.face_iter()
            sage: tup = tuple((it.vertex_repr(), it.facet_repr()) for _ in it)
            sage: tup1 = tuple((it1.vertex_repr(), it1.facet_repr()) for _ in it1)
            sage: tup == tup1
            True
        """
        nr_lines = None
        if self._unbounded:
            nr_lines = smallInteger(self._nr_lines)
        # Give a constructor by list of facets.
        return (CombinatorialPolyhedron, (self.facets(),
                self.Vrepresentation(), self.Hrepresentation(), nr_lines))

    def Vrepresentation(self):
        r"""
        Return a list of names of ``[vertices, rays, lines]``.

        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0,0], [0,1,0], \
            ....:                      [0,0,1],[0,0,-1]])
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
            return tuple(smallInteger(i) for i in range(self._length_Vrep))

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
            return tuple(smallInteger(i) for i in range(self._length_Hrep))

    def dimension(self):
        r"""
        Return the dimension of the Polyhedron.

        EXAMPLES::

            sage: C = CombinatorialPolyhedron([(1,2,3), (1,2,4),
            ....:                              (1,3,4), (2,3,4)])
            sage: C.dimension()
            3

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
            sage: CombinatorialPolyhedron(P).dimension()
            3

        TESTS::

            sage: from itertools import combinations
            sage: N = combinations(range(1200), 1199)
            sage: C = CombinatorialPolyhedron(N)
            sage: try:
            ....:     alarm(0.1)
            ....:     C.dimension()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: try:
            ....:     alarm(0.1)
            ....:     C.f_vector()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
        """
        if self._dimension == -2:
            # Dimension not calculated yet.
            if self._nr_facets == 0:
                # The dimension of a trivial Polyhedron is assumed to contain
                # exactly one "vertex" and for each dimension one "line" as in
                # :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`
                self._dimension = self._length_Vrep - 1
            elif self._unbounded or self._nr_facets <= self._length_Vrep:
                self._dimension = self.bitrep_facets.calculate_dimension()
            else:
                # If the Polyhedron has many facets,
                # calculating the dimenion of the dual will be faster.
                # The dual exists, if the Polyhedron is bounded.
                self._dimension = self.bitrep_facets.calculate_dimension()
        return smallInteger(self._dimension)

    def nr_vertices(self):
        r"""
        Return the number of vertices.

        Is equivalent to ``len(self.vertices())``.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_vertices()
            8

            sage: P = polytopes.cyclic_polytope(4,20)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_vertices()
            20

            sage: P = Polyhedron(lines=[[0,1]], vertices=[[1,0], [-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_vertices()
            0

            sage: P = Polyhedron(rays=[[1,0,0], [0,1,0]], lines=[[0,0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_vertices()
            0

            sage: C = CombinatorialPolyhedron(4)
            sage: C.f_vector()
            (1, 0, 0, 0, 0, 1)
            sage: C.nr_vertices()
            0

            sage: C = CombinatorialPolyhedron(0)
            sage: C.f_vector()
            (1, 1)
            sage: C.nr_vertices()
            1
        """
        if self.dimension() == 0:
            # This specific trivial Polyhedron needs special attention.
            return Integer(1)
        elif not self._unbounded:
            # In the unbounded case, we need to actually calculated the vertices,
            # the the V-representation contains also ``lines`` and ``rays``.
            return Integer(self._length_Vrep)
        else:
            return len(self.vertices())

    def vertices(self, names=True):
        r"""
        Return the elements in the ``Vrepresentation`` that are vertices.

        In case of an unbounded Polyhedron, there might be lines and
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
            (A vertex at (1, 0, 0),
             A vertex at (0, 1, 0),
             A vertex at (0, 0, 1),
             A vertex at (0, 0, -1),
             A vertex at (0, -1, 0),
             A vertex at (-1, 0, 0))

            sage: C.vertices(names=False)
            (5, 4, 3, 2, 1, 0)

            sage: points = [(1,0,0), (0,1,0), (0,0,1),
            ....:           (-1,0,0), (0,-1,0), (0,0,-1)]
            sage: L = LatticePolytope(points)
            sage: C = CombinatorialPolyhedron(L)
            sage: C.vertices()
            (M(0, 0, -1), M(0, -1, 0), M(-1, 0, 0), M(0, 0, 1), M(0, 1, 0), M(1, 0, 0))
            sage: C.vertices(names=False)
            (5, 4, 3, 2, 1, 0)

            sage: P = Polyhedron(vertices=[[0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.vertices()
            (A vertex at (0, 0),)
        """
        if unlikely(self.dimension() == 0):
            # Handling the case of a trivial Polyhedron of dimension `0`.
            if names and self._V:
                return (self._V[0],)
            else:
                return (Integer(0),)
        dual = False
        if not self._unbounded:
            # In the bounded case, we already know all the vertices.
            # Whereas, in the unbounded case, some of those "vertices" might
            # be ``rays`` or ``lines``.
            dual = True

        # Get all `0`-dimensional faces from :meth:`face_iter`.
        face_iter = self.face_iter(0, dual=dual)
        return tuple(face_iter.vertex_repr(names=names)[0] for _ in face_iter)

    def nr_facets(self):
        r"""
        Return the number of facets.

        Is equivalent to ``len(self.facets())``.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_facets()
            6

            sage: P = polytopes.cyclic_polytope(4,20)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_facets()
            170

            sage: P = Polyhedron(lines=[[0,1]], vertices=[[1,0], [-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_facets()
            2

            sage: P = Polyhedron(rays=[[1,0], [-1,0], [0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_facets()
            1

            sage: C = CombinatorialPolyhedron(-1)
            sage: C.f_vector()
            (1,)
            sage: C.nr_facets()
            0

            sage: C = CombinatorialPolyhedron(0)
            sage: C.f_vector()
            (1, 1)
            sage: C.nr_facets()
            1
        """
        if unlikely(self.dimension() == 0):
            # This trivial Polyhedron needs special attention.
            return Integer(1)
        return Integer(self._nr_facets)

    def facets(self, names=True):
        r"""
        Return the facets as lists of ``[vertices, rays, lines]``.

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
        if unlikely(self.dimension() == 0):
            # Special attention for this trivial case.
            # There is actually one facet, but we have not initialized it.
            return ((),)

        # Get all facets from :meth:`face_iter`.
        face_iter = self.face_iter(self.dimension() - 1, dual=False)
        tup = tuple(face_iter.vertex_repr(names=names) for _ in face_iter)

        # It is important to have the facets in the exact same order as
        # on input, so that pickle/unpickle by :meth:`reduce` works.
        # Every facet knows its index by the facet-representation.
        face_iter = self.face_iter(self.dimension() - 1, dual=False)
        indices = tuple(face_iter.facet_repr(names=False)[0] for _ in face_iter)
        dic = {}
        for i in range(len(tup)):
            dic[indices[i]] = tup[i]
        return tuple(dic[i] for i in range(len(tup)))

    def edges(self, names=True):
        r"""
        Return the edges of the Polyhedron, i.e. the rank 1 faces.

        If ``names`` is set to ``False``, then the vertices in the edges
        are given by their indices in the Vrepresentation.

        .. NOTE::

            To compute edges and f_vector, first compute the edges.
            This might be faster.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(3,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ((A vertex at (3, 9, 27), A vertex at (4, 16, 64)),
             (A vertex at (2, 4, 8), A vertex at (4, 16, 64)),
             (A vertex at (1, 1, 1), A vertex at (4, 16, 64)),
             (A vertex at (0, 0, 0), A vertex at (4, 16, 64)),
             (A vertex at (2, 4, 8), A vertex at (3, 9, 27)),
             (A vertex at (0, 0, 0), A vertex at (3, 9, 27)),
             (A vertex at (1, 1, 1), A vertex at (2, 4, 8)),
             (A vertex at (0, 0, 0), A vertex at (2, 4, 8)),
             (A vertex at (0, 0, 0), A vertex at (1, 1, 1)))

            sage: C.edges(names=False)
            ((3, 4), (2, 4), (1, 4), (0, 4), (2, 3), (0, 3), (1, 2), (0, 2), (0, 1))

            sage: P = Polyhedron(rays=[[-1,0],[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ((A line in the direction (1, 0), A vertex at (0, 0)),)

            sage: P = Polyhedron(vertices=[[0,0],[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ((A vertex at (0, 0), A vertex at (1, 0)),)

            sage: from itertools import combinations
            sage: N = combinations(['a','b','c','d','e'], 4)
            sage: C = CombinatorialPolyhedron(N)
            sage: C.edges()
            (('d', 'e'),
             ('c', 'e'),
             ('b', 'e'),
             ('a', 'e'),
             ('c', 'd'),
             ('b', 'd'),
             ('a', 'd'),
             ('b', 'c'),
             ('a', 'c'),
             ('a', 'b'))

        TESTS::

            sage: from itertools import combinations
            sage: N = combinations(range(200),199)
            sage: C = CombinatorialPolyhedron(N)
            sage: try:
            ....:     alarm(0.1)
            ....:     C.edges()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: len(C.edges())
            19900
        """
        cdef size_t len_edge_list = self._length_edges_list
        if self._edges is NULL:
            # Calculate the edges.
            if self._unbounded:
                self._calculate_edges(dual=False)
            elif self._length_Vrep > self._nr_facets*self._nr_facets:
                # This is a wild estimate
                # that in this case it is better not to use the dual.
                self._calculate_edges(dual=False)
            else:
                # In most bounded cases, one should use the dual.
                self._calculate_ridges(dual=True)
        if self._edges is NULL:
            raise ValueError('could not determine edges')

        # The edges are being saved in a list basically.
        # The first entry represents the first vertex of the first edge.
        # The second entry represents the second vertex of that edge.
        # The third entry represents the first vertex of the second edge etc.

        # Possibly there are many edges.
        # Hence, edges are stored in an array of arrays,
        # with each array containing ``len_edge_list`` of edges.

        # Mapping the indices of the vertices to the names, if requested.
        if self._V is not None and names is True:
            def f(size_t i): return self._V[i]
        else:
            def f(size_t i): return smallInteger(i)

        # Getting the indices of the `i`-th edge.
        def vertex_one(size_t i):
            return f(self._edges[i // len_edge_list][2*(i % len_edge_list)])
        def vertex_two(size_t i):
            return f(self._edges[i // len_edge_list][2*(i % len_edge_list)+1])

        cdef size_t j
        return tuple((vertex_one(j), vertex_two(j)) for j in range(self._nr_edges))

    def edge_graph(self, names=True):
        r"""
        Return the edge graph.

        If ``names`` is set to ``False``, the vertices will carry names
        according to the indexing of the Vrepresentation.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(3,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edge_graph()
            Graph on 5 vertices
            sage: G = C.edge_graph()
            sage: G.degree()
            [4, 3, 4, 3, 4]
        """
        return Graph(self.edges(names=names), format="list_of_edges")

    def ridges(self, add_equalities=False, names=True):
        r"""
        Return the ridges.

        The ridges of a Polyhedron are the faces
        contained in exactly two facets.

        To obtain all faces of codimnesion 1 use
        :meth:`CombinatorialPolyhedron.face_iter` instead.

        The ridges will be given by the facets, they are contained in.

        INPUT:

        - ``add_equalities`` -- if ``True``, then equalities of the Polyhedron
          will be added

        - ``names`` -- if ``False``, then the facets are given by their indices

        .. NOTE::

            To compute ridges and f_vector, compute the ridges first.
            This might be faster.

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
            ((An inequality (24, -26, 9, -1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-12, 19, -8, 1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-12, 19, -8, 1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (-12, 19, -8, 1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (-12, 19, -8, 1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (8, -14, 7, -1) x + 0 >= 0))
            sage: C.ridges(names=False)
            ((3, 4),
             (2, 4),
             (1, 4),
             (0, 4),
             (2, 3),
             (1, 3),
             (0, 3),
             (1, 2),
             (0, 2),
             (0, 1))

            sage: P = Polyhedron(rays=[[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C
            Combinatorial Type of a Polyhedron of dimension 1 with 1 vertices
            sage: C.ridges()
            ()
            sage: it = C.face_iter(0)
            sage: for d in it: it.facet_repr()
            (An inequality (1, 0) x + 0 >= 0, An equation (0, 1) x + 0 == 0)

        TESTS::

            sage: from itertools import combinations
            sage: N = combinations(range(200),199)
            sage: C = CombinatorialPolyhedron(N)
            sage: try:
            ....:     alarm(0.1)
            ....:     C.ridges()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: len(C.ridges())
            19900
        """
        cdef size_t len_ridge_list = self._length_edges_list
        if self._ridges is NULL:
            # Calculate the ridges.
            if self._unbounded:
                self._calculate_ridges(dual=False)
            elif self._length_Vrep*self._length_Vrep < self._nr_facets:
                # This is a wild estimate
                # that in this case it is better to use the dual.
                self._calculate_edges(dual=True)
            else:
                # In most bounded cases, one should not use the dual.
                self._calculate_ridges(dual=False)
        if self._ridges is NULL:
            raise ValueError('could not determine ridges')
        nr_ridges = self._nr_ridges

        # The ridges are being saved in a list basically.
        # The first entry represents the first facet of the first ridge.
        # The second entry represents the second facet of that ridge.
        # The third entry represents the first facet of the second ridge etc.

        # Possibly there are many ridges.
        # Hence, ridges are stored in an array of arrays,
        # with each array containing ``len_ridge_list`` of ridges.

        # Mapping the indices of the vertices to the names, if requested.
        if self._H is not None and names is True:
            def f(size_t i): return self._H[i]
        else:
            def f(size_t i): return smallInteger(i)

        # Getting the indices of the `i`-th ridge.
        def facet_one(size_t i):
            return f(self._ridges[i // len_ridge_list][2*(i % len_ridge_list)])
        def facet_two(size_t i):
            return f(self._ridges[i // len_ridge_list][2*(i % len_ridge_list)+1])

        cdef size_t j
        if add_equalities:
            # Also getting the equalities for each facet.
            return tuple(
                ((self._equalities + (facet_one(i),)),
                 (self._equalities + (facet_two(i),)))
                for i in range(nr_ridges))
        else:
            return tuple((facet_one(i), facet_two(i))
                         for i in range(nr_ridges))

    def ridge_graph(self, names=True):
        r"""
        Return the ridge graph.

        The ridge graph of a Polyhedron consists of
        ridges as edges and facets as vertices.

        If ``names`` is ``False``, the ``vertices`` of the graph  will
        be the incidences of the facets in the Hrepresentation.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(4,6)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.ridge_graph()
            Graph on 9 vertices
        """
        return Graph(self.ridges(names=names), format="list_of_edges")

    def f_vector(self):
        r"""
        Calculate the ``f_vector`` of the Polyhedron.

        The ``f_vector`` contains the number of faces of dimension `k`
        for each `k` in ``range(-1, self.dimension() + 1)``.

        .. NOTE::

            To obtain edges and/or ridges as well, first do so. This might
            already calculate the ``f_vector``.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.f_vector()
            (1, 120, 240, 150, 30, 1)

            sage: P = polytopes.cyclic_polytope(6,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.f_vector()
            (1, 10, 45, 120, 185, 150, 50, 1)

        TESTS::

            sage: from itertools import combinations
            sage: N = combinations(range(25),24)
            sage: C = CombinatorialPolyhedron(N)
            sage: try:
            ....:     alarm(0.5)
            ....:     C.f_vector()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: C.f_vector()  # long time
            (1,
             25,
             300,
             2300,
             12650,
             53130,
             177100,
             480700,
             1081575,
             2042975,
             3268760,
             4457400,
             5200300,
             5200300,
             4457400,
             3268760,
             2042975,
             1081575,
             480700,
             177100,
             53130,
             12650,
             2300,
             300,
             25,
             1)
        """
        if not self._f_vector:
            self._calculate_f_vector()
        if not self._f_vector:
            raise ValueError("could not determine f_vector")
        return self._f_vector

    def face_iter(self, dimension=None, dual=None):
        r"""
        Iterator over all proper faces of specified dimension.

        INPUT:

        - ``dimension`` -- if specified, then iterate over only this dimension
        - ``dual`` -- if ``True``, iterate starting with the vertices,
          if ``False``, iterate starting with the facets,
          if ``None``, choose ``True`` or ``False`` for speed

        OUTPUT:

        - :class:`FaceIterator`

        .. NOTE::

            :class:`FaceIterator` is more than just a plain iterator.
            By default it will iterate over the dimensions of the faces, but
            more information can be received.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(dimension=2)
            sage: next(it)
            2
            sage: it.vertex_repr()
            (A vertex at (4, 1, 5, 2, 3),
             A vertex at (4, 2, 5, 1, 3),
             A vertex at (5, 1, 4, 2, 3),
             A vertex at (5, 2, 4, 1, 3))
            sage: next(it)
            2
            sage: it.vertex_repr()
            (A vertex at (4, 1, 5, 2, 3),
             A vertex at (4, 1, 5, 3, 2),
             A vertex at (5, 1, 4, 2, 3),
             A vertex at (5, 1, 4, 3, 2))
            sage: it.facet_repr()
            (An inequality (0, 1, 0, 0, 0) x - 1 >= 0,
             An inequality (0, 1, 0, 1, 1) x - 6 >= 0,
             An equation (1, 1, 1, 1, 1) x - 15 == 0)
            sage: it.facet_repr(names=False)
            (25, 29)
            sage: next(it)
            2
            sage: it.facet_repr(names=False)
            (12, 29)
            sage: it.vertex_repr(names=False)
            (76, 77, 82, 83, 88, 89)

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],[0,2,3],[1,2,3]])
            sage: it = C.face_iter()
            sage: for _ in it: it.vertex_repr()
            (1, 2, 3)
            (0, 2, 3)
            (0, 1, 3)
            (0, 1, 2)
            (2, 3)
            (1, 3)
            (1, 2)
            (3,)
            (2,)
            (1,)
            (0, 3)
            (0, 2)
            (0,)
            (0, 1)

            sage: P = Polyhedron(rays=[[1,0],[0,1]], vertices=[[1,0],[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(1)
            sage: for _ in it: it.vertex_repr()
            (A vertex at (0, 1), A vertex at (1, 0))
            (A ray in the direction (1, 0), A vertex at (1, 0))
            (A ray in the direction (0, 1), A vertex at (0, 1))

        .. SEEALSO::

            :class:`FaceIterator`.
        """
        cdef FaceIterator face_iter
        if dual is None:
            # Determine the faster way, to iterate through all faces.
            if self._unbounded or self._nr_facets <= self._length_Vrep:
                dual = False
            else:
                dual = True

        face_iter = self._face_iter(int(dual))
        if dimension is not None:
            # Setting the face iterator to return only requested dimension.
            face_iter.set_request_dimension(dimension)
        return face_iter

    cdef FaceIterator _face_iter(self, bint dual):
        r"""
        A method to obtain the FaceIterator as Cython object.

        See :meth:`CombinatorialPolyhedron.face_iter`
        """
        if dual and self._unbounded:
            raise ValueError("cannot iterate over dual of unbounded Polyhedron")
        return FaceIterator(self, dual)

    def face_lattice(self):
        r"""
        Generate the face-lattice.

        OUTPUT:

        - :class:'~sage.combinat.posets.lattices.FiniteLatticePoset'

        .. NOTE::

            Use :meth:`CombinatorialPolyhedron.face_lattice_dimension` to get
            the dimension for each element in the Face Lattice.
            Use :meth:`CombinatorialPolyhedron.face_lattice_vertex_repr` to get
            the vertex representation for each element in the Face Lattice.
            Use :meth:`CombinatorialPolyhedron.face_lattice_facet_repr` to get
            the facet_repr for each element in the Face Lattice.

        .. WARNING::

            The labeling of the face lattice might depend on archicture
            and implementation.
            Relabeling the face lattice with
            :meth:`CombinatorialPolyhedron.face_lattice_vertex_repr` and/or
            :meth:`CombinatorialPolyhedron.face_lattice_facet_repr` will be
            platform independent.

        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0],[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice()
            Finite lattice containing 5 elements

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: P1 = Polyhedron(rays=[[1,0], [-1,0]])
            sage: C1 = CombinatorialPolyhedron(P1)
            sage: C.face_lattice().is_isomorphic(C1.face_lattice())
            True

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice()
            Finite lattice containing 542 elements

        TESTS::

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice().is_isomorphic(P.face_lattice())
            True

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice().is_isomorphic(P.face_lattice())
            True
        """
        if not self._face_lattice_incidences:
            # Calculate all incidences.
            self._calculate_face_lattice_incidences()
        if self._face_lattice_incidences is NULL:
            raise TypeError("could not determine face lattice")

        cdef size_t **incidences = self._face_lattice_incidences
        cdef size_t nr_incidences = self._nr_face_lattice_incidences
        cdef size_t len_incidence_list = self._length_edges_list

        # The incidences are being saved in a list basically.
        # The first entry represents the first face of the first incidence.
        # The second entry represents the second face of that incidence.
        # The third entry represents the first face of the second incidence etc.

        # Possibly there are many incidences.
        # Hence, incidences are stored in an array of arrays,
        # with each array containing ``len_incidence_list`` of incidences.

        # Getting the indices of the `i`-th incidence.
        def face_one(size_t i):
            return Integer(incidences[i // len_incidence_list][2*(i % len_incidence_list)])
        def face_two(size_t i):
            return Integer(incidences[i // len_incidence_list][2*(i % len_incidence_list)+1])

        # Edges of the face-lattice/Hasse diagram.
        cdef size_t j
        edges = tuple((face_one(j), face_two(j))
                          for j in range(nr_incidences))

        V = tuple(range(sum(self._f_vector)))
        D = DiGraph([V, edges], format='vertices_and_edges')
        return FiniteLatticePoset(D)

    def face_lattice_dimension(self, index):
        r"""
        Return for each element in :meth:`CombinatorialPolyhedron.face_lattice`
        its dimension.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: def f(i):
            ....:     return (i, C.face_lattice_dimension(i))
            ....:
            sage: G = F.relabel(f)
            sage: set(G._elements)
            {(0, -1),
             (1, 0),
             (2, 0),
             (3, 0),
             (4, 0),
             (5, 0),
             (6, 0),
             (7, 0),
             (8, 0),
             (9, 1),
             (10, 1),
             (11, 1),
             (12, 1),
             (13, 1),
             (14, 1),
             (15, 1),
             (16, 1),
             (17, 1),
             (18, 1),
             (19, 1),
             (20, 1),
             (21, 2),
             (22, 2),
             (23, 2),
             (24, 2),
             (25, 2),
             (26, 2),
             (27, 3)}
        """
        f_vector = self.f_vector()
        dim = self.dimension()

        # Getting the dimension, by considering the following:
        # The level-set of dimension `d` will have indices `k, k+1, ..., k+n-1`,
        # where `n` is the number of faces of dimension `d` ( ``n = f_vector[d + 1]``)
        # and `k` is the number of face of dimension up to `d`, i.e.
        # ``k = sum(f_vector[:d])``.
        return max(d for d in range(dim+2) if sum(f_vector[:d]) <= index) - 1

    def face_lattice_vertex_repr(self, index, names=True):
        r"""
        Return for each element in :meth:`CombinatorialPolyhedron.face_lattice`
        its vertex-representation as in :meth:`ListOfAllFaces.vertex_repr` or
        :meth:`FaceIterator.vertex_repr`.

        If ``names`` is set to ``True``, then names of the
        ``[vertices, rays, lines]`` are used.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: F
            Finite lattice containing 28 elements
            sage: G = F.relabel(C.face_lattice_vertex_repr)
            sage: G._elements
            ((),
             (A vertex at (-1, -1, -1),),
             (A vertex at (-1, -1, 1),),
             (A vertex at (-1, -1, -1), A vertex at (-1, -1, 1)),
             (A vertex at (-1, 1, -1),),
             (A vertex at (-1, -1, -1), A vertex at (-1, 1, -1)),
             (A vertex at (-1, 1, 1),),
             (A vertex at (-1, -1, 1), A vertex at (-1, 1, 1)),
             (A vertex at (-1, 1, -1), A vertex at (-1, 1, 1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, -1, 1),
              A vertex at (-1, 1, -1),
              A vertex at (-1, 1, 1)),
             (A vertex at (1, -1, -1),),
             (A vertex at (-1, -1, -1), A vertex at (1, -1, -1)),
             (A vertex at (1, -1, 1),),
             (A vertex at (-1, -1, 1), A vertex at (1, -1, 1)),
             (A vertex at (1, -1, -1), A vertex at (1, -1, 1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, -1, 1),
              A vertex at (1, -1, -1),
              A vertex at (1, -1, 1)),
             (A vertex at (1, 1, -1),),
             (A vertex at (-1, 1, -1), A vertex at (1, 1, -1)),
             (A vertex at (1, -1, -1), A vertex at (1, 1, -1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, 1, -1),
              A vertex at (1, -1, -1),
              A vertex at (1, 1, -1)),
             (A vertex at (1, 1, 1),),
             (A vertex at (-1, 1, 1), A vertex at (1, 1, 1)),
             (A vertex at (1, -1, 1), A vertex at (1, 1, 1)),
             (A vertex at (-1, -1, 1),
              A vertex at (-1, 1, 1),
              A vertex at (1, -1, 1),
              A vertex at (1, 1, 1)),
             (A vertex at (1, 1, -1), A vertex at (1, 1, 1)),
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
              A vertex at (-1, 1, 1),
              A vertex at (1, -1, -1),
              A vertex at (1, -1, 1),
              A vertex at (1, 1, -1),
              A vertex at (1, 1, 1)))

            sage: P = Polyhedron(rays=[[0,1], [1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: G = F.relabel(C.face_lattice_vertex_repr)
            sage: G._elements
            ((),
             (A vertex at (0, 0),),
             (A vertex at (0, 0), A ray in the direction (0, 1)),
             (A vertex at (0, 0), A ray in the direction (1, 0)),
             (A vertex at (0, 0),
              A ray in the direction (0, 1),
              A ray in the direction (1, 0)))

            sage: def f(i): return C.face_lattice_vertex_repr(i, False)
            sage: G = F.relabel(f)
            sage: G._elements
            ((), (0,), (0, 1), (0, 2), (0, 1, 2))
        """
        self._record_all_faces()                            # Initalize ``_all_faces``, if not done yet.
        dim = self.face_lattice_dimension(index)            # Determine dimension to that index.
        newindex = index - sum(self._f_vector[:dim + 1])    # Index in that level-set.

        # Let ``_all_faces`` determine vertex-representation.
        return self._all_faces.vertex_repr(dim, newindex, names=names)

    def face_lattice_facet_repr(self, index, names=True):
        r"""
        Return for each element in :meth:`CombinatorialPolyhedron.face_lattice`
        its facet-representation as in :meth:`ListOfAllFaces.facet_repr` or
        :meth:`FaceIterator.facet_repr`.

        If ``names`` is set to ``True``, then names of the
        ``facets`` are used.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: F
            Finite lattice containing 28 elements
            sage: G = F.relabel(C.face_lattice_facet_repr)
            sage: G._elements
            ((An inequality (0, 0, -1) x + 1 >= 0,
              An inequality (0, -1, 0) x + 1 >= 0,
              An inequality (-1, 0, 0) x + 1 >= 0,
              An inequality (1, 0, 0) x + 1 >= 0,
              An inequality (0, 0, 1) x + 1 >= 0,
              An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0,
              An inequality (0, -1, 0) x + 1 >= 0,
              An inequality (1, 0, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0,
              An inequality (1, 0, 0) x + 1 >= 0,
              An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0, An inequality (1, 0, 0) x + 1 >= 0),
             (An inequality (1, 0, 0) x + 1 >= 0,
              An inequality (0, 0, 1) x + 1 >= 0,
              An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (1, 0, 0) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0,
              An inequality (0, -1, 0) x + 1 >= 0,
              An inequality (-1, 0, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0, An inequality (0, -1, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0,
              An inequality (-1, 0, 0) x + 1 >= 0,
              An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0, An inequality (-1, 0, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0,),
             (An inequality (0, -1, 0) x + 1 >= 0,
              An inequality (1, 0, 0) x + 1 >= 0,
              An inequality (0, 0, 1) x + 1 >= 0),
             (An inequality (1, 0, 0) x + 1 >= 0, An inequality (0, 0, 1) x + 1 >= 0),
             (An inequality (0, -1, 0) x + 1 >= 0, An inequality (1, 0, 0) x + 1 >= 0),
             (An inequality (1, 0, 0) x + 1 >= 0,),
             (An inequality (0, -1, 0) x + 1 >= 0,
              An inequality (-1, 0, 0) x + 1 >= 0,
              An inequality (0, 0, 1) x + 1 >= 0),
             (An inequality (0, -1, 0) x + 1 >= 0, An inequality (0, 0, 1) x + 1 >= 0),
             (An inequality (0, -1, 0) x + 1 >= 0, An inequality (-1, 0, 0) x + 1 >= 0),
             (An inequality (0, -1, 0) x + 1 >= 0,),
             (An inequality (-1, 0, 0) x + 1 >= 0,
              An inequality (0, 0, 1) x + 1 >= 0,
              An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (-1, 0, 0) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 0, 1) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 1, 0) x + 1 >= 0,),
             (An inequality (-1, 0, 0) x + 1 >= 0, An inequality (0, 0, 1) x + 1 >= 0),
             (An inequality (0, 0, 1) x + 1 >= 0,),
             (An inequality (-1, 0, 0) x + 1 >= 0,),
             ())

            sage: P = Polyhedron(rays=[[0,1], [1,0]], vertices=[[0,1], [1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: G = F.relabel(C.face_lattice_facet_repr)
            sage: G._elements
            ((An inequality (1, 0) x + 0 >= 0,
              An inequality (0, 1) x + 0 >= 0,
              An inequality (1, 1) x - 1 >= 0),
             (An inequality (1, 0) x + 0 >= 0, An inequality (1, 1) x - 1 >= 0),
             (An inequality (1, 0) x + 0 >= 0,),
             (An inequality (0, 1) x + 0 >= 0, An inequality (1, 1) x - 1 >= 0),
             (An inequality (0, 1) x + 0 >= 0,),
             (An inequality (1, 1) x - 1 >= 0,),
             ())

            sage: def f(i): return C.face_lattice_facet_repr(i, False)
            sage: G = F.relabel(f)
            sage: G._elements
            ((0, 1, 2), (0, 2), (1, 2), (1,), (2,), (0,), ())
        """
        self._record_all_faces()                            # Initalize ``_all_faces``, if not done yet.
        dim = self.face_lattice_dimension(index)            # Determine dimension to that index.
        newindex = index - sum(self._f_vector[:dim + 1])    # Index in that level-set.

        # Let ``_all_faces`` determine facet-representation.
        return self._all_faces.facet_repr(dim, newindex, names=names)

    cdef int _calculate_f_vector(self) except -1:
        r"""
        Calculate the ``f_vector`` of the Polyhedron.

        See :meth:`f_vector`.
        """
        if self._f_vector:
            return 0  # There is no need to recalculate the f_vector.

        cdef bint dual
        if self._unbounded or self._nr_facets <= self._length_Vrep:
            # In this case the non-dual approach is faster..
            dual = False
        else:
            # In this case the dual approach is faster.
            dual = True
        cdef FaceIterator face_iter = self._face_iter(dual)

        cdef int dim = self.dimension()
        cdef int d  # dimension of the current face of the iterator
        cdef MemoryAllocator mem = MemoryAllocator()

        # Initialize ``f_vector``.
        cdef size_t *f_vector = <size_t *> mem.calloc((dim + 2), sizeof(size_t))
        f_vector[0] = 1         # Face iterator will only visit proper faces.
        f_vector[dim + 1] = 1   # Face iterator will only visit proper faces.

        # For each face in the iterator, add `1` to the corresponding entry in
        # ``f_vector``.
        if self._nr_facets > 0 and dim > 0:
            d = face_iter.next_face()
            while (d < dim):
                sig_check()
                f_vector[d+1] += 1
                d = face_iter.next_face()

        # Copy ``f_vector``.
        if dual:
            # We have calculated the ``f_vector`` of the dual.
            # Reverse it:
            self._f_vector = \
                tuple(smallInteger(f_vector[dim+1-i]) for i in range(dim+2))
        else:
            self._f_vector = tuple(smallInteger(f_vector[i]) for i in range(dim+2))

    cdef int _calculate_edges(self, dual) except -1:
        r"""
        Calculate the edges of the Polyhedron.

        If ``dual`` is ``True``, calculate the edges of the dual. In this case,
        this will also calculate the ``f_vector``, if unknown.

        See :meth:`CombinatorialPolyhedron.edges` and :meth:`CombinatorialPolyhedron.ridges`.
        """
        if (self._edges is not NULL and not dual) or (self._ridges is not NULL and dual):
            return 0  # There is no need to recalculate.

        cdef MemoryAllocator mem = MemoryAllocator()
        cdef FaceIterator face_iter = self._face_iter(dual)
        cdef size_t len_edge_list = self._length_edges_list
        cdef int dim = self.dimension()
        cdef int d              # dimension of the current face of ``FaceIterator``
        cdef size_t *f_vector   # calculate f_vector, if not done already
        cdef bint is_f_vector   # True if f_vector was calculated previously

        # For each edge we determine its location in ``edges``
        # by ``edges[one][two]``
        cdef size_t **edges = <size_t**> mem.malloc(sizeof(size_t*))
        cdef size_t one, two

        cdef size_t counter = 0         # the number of edges so far
        cdef size_t current_length = 1  # dynamically enlarge **edges

        if self._f_vector:
            is_f_vector = True
        else:
            # in this case we will calculate the f_vector while we're at it
            is_f_vector = False

        if dim == 1:
            # In this case there is an edge, but its not a proper face.
            edges[0] = <size_t *> mem.allocarray(2, sizeof(size_t))
            edges[0][0] = 0
            edges[0][1] = 1
            counter = 1

            # Success, copy the data to ``CombinatorialPolyhedron``.
            if dual:
                # We have actually calculated the ridges.
                sig_block()
                self._nr_ridges = counter
                self._ridges = edges
                self._mem_tuple += (mem,)
                sig_unblock()
            else:
                sig_block()
                self._nr_edges = counter
                self._edges = edges
                self._mem_tuple += (mem,)
                sig_unblock()
            return 0

        if is_f_vector:
            # Only calculate the edges.

            if not dual:
                face_iter.set_request_dimension(1)
            else:
                # :meth:`FaceIterator.set_request_dimension`
                # requires the dimension of the original Polyhedron
                face_iter.set_request_dimension(dim - 2)

            if self._nr_facets > 0 and dim > 0:
                # If not, there won't even be any edges. Prevent error message.

                while (face_iter.next_face() == 1):

                    # Determine the position in ``edges``.
                    one = counter // len_edge_list
                    two = counter % len_edge_list

                    # Enlarge ``edges`` if needed.
                    if unlikely(two == 0):
                        if unlikely(one + 1 > current_length):
                            # enlarge **edges
                            current_length *= 2
                            edges = <size_t **> mem.reallocarray(edges, current_length, sizeof(size_t*))

                        edges[one] = <size_t *> mem.allocarray(2 * len_edge_list, sizeof(size_t))

                    # Set up face_iter.atom_repr
                    face_iter.set_atom_repr()

                    # Copy the information.
                    edges[one][2*two] = face_iter.atom_repr[0]
                    edges[one][2*two + 1] = face_iter.atom_repr[1]
                    counter += 1

            # Success, copy the data to ``CombinatorialPolyhedron``.
            if dual:
                sig_block()
                self._nr_ridges = counter
                self._ridges = edges
                self._mem_tuple += (mem,)
                sig_unblock()
            else:
                sig_block()
                self._nr_edges = counter
                self._edges = edges
                self._mem_tuple += (mem,)
                sig_unblock()
        else:
            # While doing the edges one might as well do the f-vector.
            f_vector = <size_t *> mem.calloc(dim + 2, sizeof(size_t))
            f_vector[0] = 1         # This is not a proper face.
            f_vector[dim + 1] = 1   # This is not a proper face.

            counter = 0
            if self._nr_facets > 0 and dim > 0:
                # If not, there won't even be any edges. Prevent error message.

                d = face_iter.next_face()
                while (d < dim):
                    f_vector[d+1] += 1

                    if d == 1:
                        # If it is an edge.

                         # Determine the position in ``edges``.
                        one = counter // len_edge_list
                        two = counter % len_edge_list

                        # Enlarge ``edges`` if needed.
                        if unlikely(two == 0):
                            if unlikely(one + 1 > current_length):
                                # enlarge **edges
                                current_length *= 2
                                edges = <size_t **> mem.reallocarray(edges, current_length, sizeof(size_t*))

                            edges[one] = <size_t *> mem.allocarray(2 * len_edge_list, sizeof(size_t))

                        # Set up face_iter.atom_repr
                        face_iter.set_atom_repr()

                        # Copy the information.
                        edges[one][2*two] = face_iter.atom_repr[0]
                        edges[one][2*two + 1] = face_iter.atom_repr[1]
                        counter += 1

                    d = face_iter.next_face()  # Go to next face.

            # Success, copy the data to ``CombinatorialPolyhedron``.
            if dual:
                sig_block()
                self._f_vector = \
                    tuple(smallInteger(f_vector[dim+1-i]) for i in range(dim+2))
                self._nr_ridges = counter
                self._ridges = edges
                self._mem_tuple += (mem,)
                sig_unblock()
            else:
                sig_block()
                self._f_vector = \
                    tuple(smallInteger(f_vector[i]) for i in range(dim+2))
                self._nr_edges = counter
                self._edges = edges
                self._mem_tuple += (mem,)
                sig_unblock()

    cdef int _calculate_ridges(self, dual) except -1:
        r"""
        Calculate the ridges of the Polyhedron.

        If ``dual`` is ``True``, calculate the ridges of the polar.

        See :meth:`edges` and :meth:`ridges`.
        """
        if (self._edges is not NULL and dual) or (self._ridges is not NULL and not dual):
            return 0  # There is no need to recalculate.

        cdef MemoryAllocator mem = MemoryAllocator()
        cdef FaceIterator face_iter = self._face_iter(dual)
        cdef size_t len_ridge_list = self._length_edges_list
        cdef int dim = self.dimension()

        # For each ridge we determine its location in ``ridges``
        # by ``ridges[one][two]``.
        cdef size_t **ridges = <size_t**> mem.malloc(sizeof(size_t*))
        cdef size_t one, two

        cdef size_t counter = 0         # the number of ridges so far
        cdef size_t current_length = 1  # dynamically enlarge **ridges

        if dim == 1 and self._nr_facets > 1:
            # In this case there is a ridge, but its not a proper face.
            ridges[0] = <size_t *> mem.allocarray(2, sizeof(size_t))
            ridges[0][0] = 0
            ridges[0][1] = 1
            counter = 1

            # Success, copy the data to ``CombinatorialPolyhedron``.
            if not dual:
                sig_block()
                self._nr_ridges = counter
                self._ridges = ridges
                self._mem_tuple += (mem,)
                sig_unblock()
            else:
                sig_block()
                self._nr_edges = counter
                self._edges = ridges
                self._mem_tuple += (mem,)
                sig_unblock()
            return 0

        if dual:
            # :meth:`FaceIterator.set_request_dimension`
            # requires the dimension of the original Polyhedron
            face_iter.set_request_dimension(1)
        else:
            face_iter.set_request_dimension(dim - 2)

        if self._nr_facets > 1 and dim > 0:
            # If not, there won't even be any ridges
            # as intersection of two distince facets.
            # Prevent error message.

            while (face_iter.next_face() == dim - 2):

                # Determine the position in ``ridges``.
                one = counter // len_ridge_list
                two = counter % len_ridge_list

                # Enlarge ``ridges`` if needed.
                if unlikely(two == 0):
                    if unlikely(one + 1 > current_length):
                        # enlarge **ridges
                        current_length *= 2
                        ridges = <size_t **> mem.reallocarray(ridges, current_length, sizeof(size_t*))

                    ridges[one] = <size_t *> mem.allocarray(2 * len_ridge_list, sizeof(size_t))

                # Set up face_iter.coatom_repr
                face_iter.set_coatom_repr()

                # Copy the information.
                ridges[one][2*two] = face_iter.coatom_repr[0]
                ridges[one][2*two + 1] = face_iter.coatom_repr[1]
                counter += 1

        # Success, copy the data to ``CombinatorialPolyhedron``.
        if not dual:
            sig_block()
            self._nr_ridges = counter
            self._ridges = ridges
            self._mem_tuple += (mem,)
            sig_unblock()
        else:
            sig_block()
            self._nr_edges = counter
            self._edges = ridges
            self._mem_tuple += (mem,)
            sig_unblock()

    cdef int _calculate_face_lattice_incidences(self) except -1:
        r"""
        Calculate all incidences for the face lattice.

        See :meth:`face_lattice`.
        """
        if self._face_lattice_incidences:
            return 1  # There is no need to recalculate the incidences.

        cdef size_t len_incidence_list = self._length_edges_list
        cdef int dim = self.dimension()
        f_vector = self.f_vector()
        cdef MemoryAllocator mem = MemoryAllocator()
        self._record_all_faces()  # set up ``self._all_faces``
        cdef ListOfAllFaces all_faces = self._all_faces

        # ``all_faces`` will store its incidences in ``first`` and ``second``.
        cdef size_t first = 0, second = 0

        # ``dimension_one`` and ``dimension_two`` will be the dimensions of the
        # incidences, we currently obtain from ``all_faces``.
        # Almost always ``dimension_two = dimension_one - 1``.
        cdef int dimension_one, dimension_two
        cdef int j  # an index for ``range(dimension_two + 1)``

        # The indices of the incidences in ``all_faces`` are levelwise.
        # Hence, we have to add to each index dependent on dimension:

        # For ``dimension_two`` we add:
        cdef size_t already_seen       # = sum(f_vector[j] for j in range(dimension_two + 1))

        # For ``dimension_one`` we add:
        cdef size_t already_seen_next  # = sum(f_vector[j] for j in range(dimension_two + 2))

        # For each incidence we determine its location in ``incidences``
        # by ``incidences[one][two]``.
        cdef size_t **incidences = <size_t**> mem.malloc(sizeof(size_t*))
        cdef size_t one, two

        cdef size_t counter = 0         # the number of incidences so far
        cdef size_t current_length = 1  # dynamically enlarge **incidences

        if all_faces is None:
            raise ValueError("could not determine a list of all faces")

        dimension_one = 0
        if dim > -1:
            while (f_vector[dimension_one + 1] == 0):
                # Taking care of cases, where there might be no faces
                # of dimension 0, 1, etc (``self._nr_lines > 0``).
                dimension_one += 1
            dimension_two = -1

        while (dimension_one < dim + 1):
            already_seen = sum(f_vector[j] for j in range(dimension_two + 1))
            already_seen_next = already_seen + f_vector[dimension_two + 1]

            if all_faces.dual:
                # If ``dual``, then ``all_faces`` has the dimensions reversed.
                all_faces.incidence_init(dim - 1 - dimension_two, dim - 1 - dimension_one)
            else:
                all_faces.incidence_init(dimension_one, dimension_two)

            # Get all incidences for fixed ``[dimension_one, dimension_two]``.
            while all_faces.next_incidence(&second, &first):

                # Determine the position in ``incidences``.
                one = counter // len_incidence_list
                two = counter % len_incidence_list

                # Enlarge ``incidences`` if needed.
                if unlikely(two == 0):
                    if unlikely(one + 1 > current_length):
                        # enlarge **incidences
                        current_length *= 2
                        incidences = <size_t **> mem.reallocarray(incidences, current_length, sizeof(size_t*))

                    incidences[one] = <size_t *> mem.allocarray(2 * len_incidence_list, sizeof(size_t))

                if all_faces.dual:
                    # If ``dual``, then ``second`` and ``first are flipped.
                    second += already_seen
                    first += already_seen_next
                    incidences[one][2*two] = second
                    incidences[one][2*two + 1] = first
                else:
                    second += already_seen_next
                    first += already_seen
                    incidences[one][2*two] = first
                    incidences[one][2*two + 1] = second

                counter += 1
                sig_check()

            # Increase dimensions.
            dimension_one += 1
            dimension_two = dimension_one - 1

        # Success, copy the data to ``CombinatorialPolyhedron``.
        self._nr_face_lattice_incidences = counter
        sig_block()
        self._mem_tuple += (mem,)
        self._face_lattice_incidences = incidences
        sig_unblock()

    def _record_all_faces(self):
        r"""
        Initialize :class:`ListOfAllFaces` for the Polyhedron.

        Record and sort all faces of the Polyhedron in that class.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C._record_all_faces()

        TESTS::

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: tup = tuple((it.vertex_repr(),it.facet_repr()) for _ in it)
            sage: rg = range(1,sum(C.f_vector()) - 1)
            sage: tup2 = tuple((C.face_lattice_vertex_repr(i),
            ....:               C.face_lattice_facet_repr(i)) for i in rg)
            sage: sorted(tup) == sorted(tup2)
            True

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: tup = tuple((it.vertex_repr(),it.facet_repr()) for _ in it)
            sage: rg = range(1,sum(C.f_vector()) - 1)
            sage: tup2 = tuple((C.face_lattice_vertex_repr(i),
            ....:               C.face_lattice_facet_repr(i)) for i in rg)
            sage: sorted(tup) == sorted(tup2)
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: tup = tuple((it.vertex_repr(),it.facet_repr()) for _ in it)
            sage: rg = range(1,sum(C.f_vector()) - 1)
            sage: tup2 = tuple((C.face_lattice_vertex_repr(i),
            ....:               C.face_lattice_facet_repr(i)) for i in rg)
            sage: sorted(tup) == sorted(tup2)
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0],
            ....:                      [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: tup = tuple((it.vertex_repr(),it.facet_repr()) for _ in it)
            sage: rg = range(1,sum(C.f_vector()) - 1)
            sage: tup2 = tuple((C.face_lattice_vertex_repr(i),
            ....:               C.face_lattice_facet_repr(i)) for i in rg)
            sage: sorted(tup) == sorted(tup2)
            True
        """
        if self._all_faces:
            return  # Have recorded all faces already.

        self._all_faces = ListOfAllFaces(self)
        if self._all_faces is None:
            raise ValueError("could not determine a list of all faces")
