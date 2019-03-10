"""
A class to keep information about faces of a polyhedron

This module gives you a tool to work with the faces of a polyhedron
and their relative position. First, you need to find the faces. To get
the faces in a particular dimension, use the
:meth:`~sage.geometry.polyhedron.base.face` method::

    sage: P = polytopes.cross_polytope(3)
    sage: P.faces(3)
    (<0,1,2,3,4,5>,)
    sage: P.faces(2)
    (<0,1,2>, <0,1,3>, <0,2,4>, <0,3,4>, <3,4,5>, <2,4,5>, <1,3,5>, <1,2,5>)
    sage: P.faces(1)
    (<0,1>, <0,2>, <1,2>, <0,3>, <1,3>, <0,4>, <2,4>, <3,4>, <2,5>, <3,5>, <4,5>, <1,5>)

or :meth:`~sage.geometry.polyhedron.base.face_lattice` to get the
whole face lattice as a poset::

    sage: P.face_lattice()
    Finite lattice containing 28 elements with distinguished linear extension

The faces are printed in shorthand notation where each integer is the
index of a vertex/ray/line in the same order as the containing
Polyhedron's :meth:`~sage.geometry.polyhedron.base.Vrepresentation` ::

    sage: face = P.faces(1)[3];  face
    <0,3>
    sage: P.Vrepresentation(0)
    A vertex at (-1, 0, 0)
    sage: P.Vrepresentation(3)
    A vertex at (0, 0, 1)
    sage: face.vertices()
    (A vertex at (-1, 0, 0), A vertex at (0, 0, 1))

The face itself is not represented by Sage's
:func:`sage.geometry.polyhedron.constructor.Polyhedron` class, but by
an auxiliary class to keep the information. You can get the face as a
polyhedron with the :meth:`PolyhedronFace.as_polyhedron` method::

    sage: face.as_polyhedron()
    A 1-dimensional polyhedron in ZZ^3 defined as the convex hull of 2 vertices
    sage: _.equations()
    (An equation (0, 1, 0) x + 0 == 0,
     An equation (1, 0, -1) x + 1 == 0)
"""

########################################################################
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################
from __future__ import print_function

from sage.structure.sage_object import SageObject
from sage.structure.richcmp import richcmp_method, richcmp
from sage.misc.all import cached_method
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix



#########################################################################
@richcmp_method
class PolyhedronFace(SageObject):
    r"""
    A face of a polyhedron.

    This class is for use in
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.face_lattice`.

    INPUT:

    No checking is performed whether the H/V-representation indices
    actually determine a face of the polyhedron. You should not
    manually create :class:`PolyhedronFace` objects unless you know
    what you are doing.

    OUTPUT:

    A :class:`PolyhedronFace`.

    EXAMPLES::

        sage: octahedron = polytopes.cross_polytope(3)
        sage: inequality = octahedron.Hrepresentation(2)
        sage: face_h = tuple([ inequality ])
        sage: face_v = tuple( inequality.incident() )
        sage: face_h_indices = [ h.index() for h in face_h ]
        sage: face_v_indices = [ v.index() for v in face_v ]
        sage: from sage.geometry.polyhedron.face import PolyhedronFace
        sage: face = PolyhedronFace(octahedron, face_v_indices, face_h_indices)
        sage: face
        <0,1,2>
        sage: face.dim()
        2
        sage: face.ambient_Hrepresentation()
        (An inequality (1, 1, 1) x + 1 >= 0,)
        sage: face.ambient_Vrepresentation()
        (A vertex at (-1, 0, 0), A vertex at (0, -1, 0), A vertex at (0, 0, -1))
    """

    def __init__(self, polyhedron, V_indices, H_indices):
        r"""
        The constructor.

        See :class:`PolyhedronFace` for more information.

        INPUT:

        - ``polyhedron`` -- a :class:`Polyhedron`. The ambient
          polyhedron.

        - ``V_indices`` -- list of sorted integers. The indices of the
          face-spanning V-representation objects in the ambient
          polyhedron.

        - ``H_indices`` -- list of sorted integers. The indices of the
          H-representation objects of the ambient polyhedron that are
          saturated on the face.

        TESTS::

            sage: from sage.geometry.polyhedron.face import PolyhedronFace
            sage: PolyhedronFace(Polyhedron(), [], [])   # indirect doctest
            <>
        """
        self._polyhedron = polyhedron
        self._ambient_Vrepresentation_indices = tuple(V_indices)
        self._ambient_Hrepresentation_indices = tuple(H_indices)
        self._ambient_Vrepresentation = tuple( polyhedron.Vrepresentation(i) for i in V_indices )
        self._ambient_Hrepresentation = tuple( polyhedron.Hrepresentation(i) for i in H_indices )

    def __hash__(self):
        r"""
        TESTS::

            sage: P = Polyhedron([[0,0],[0,1],[23,3],[9,12]])
            sage: list(map(hash, P.faces(1)))  # random
            [2377119663630407734,
             2377136578164722109,
             5966674064902575359,
             4795242501625591634]
        """
        return hash((self._polyhedron, self._ambient_Vrepresentation_indices))

    def vertex_generator(self):
        """
        Return a generator for the vertices of the face.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: face = triangle.faces(1)[0]
            sage: for v in face.vertex_generator(): print(v)
            A vertex at (0, 1)
            A vertex at (1, 0)
            sage: type(face.vertex_generator())
            <... 'generator'>
        """
        for V in self.ambient_Vrepresentation():
            if V.is_vertex():
                yield V

    @cached_method
    def vertices(self):
        """
        Return all vertices of the face.

        OUTPUT:

        A tuple of vertices.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: face = triangle.faces(1)[0]
            sage: face.vertices()
            (A vertex at (0, 1), A vertex at (1, 0))
        """
        return tuple(self.vertex_generator())

    def ray_generator(self):
        """
        Return a generator for the rays of the face.

        EXAMPLES::

            sage: pi = Polyhedron(ieqs = [[1,1,0],[1,0,1]])
            sage: face = pi.faces(1)[0]
            sage: next(face.ray_generator())
            A ray in the direction (1, 0)
        """
        for V in self.ambient_Vrepresentation():
            if V.is_ray():
                yield V

    @cached_method
    def rays(self):
        """
        Return the rays of the face.

        OUTPUT:

        A tuple of rays.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0],[1,1,0,0]])
            sage: face = p.faces(2)[0]
            sage: face.rays()
            (A ray in the direction (1, 0, 0), A ray in the direction (0, 1, 0))
        """
        return tuple(self.ray_generator())

    def line_generator(self):
        """
        Return a generator for the lines of the face.

        EXAMPLES::

            sage: pr = Polyhedron(rays = [[1,0],[-1,0],[0,1]], vertices = [[-1,-1]])
            sage: face = pr.faces(1)[0]
            sage: next(face.line_generator())
            A line in the direction (1, 0)
        """
        for V in self.ambient_Vrepresentation():
            if V.is_line():
                yield V

    @cached_method
    def lines(self):
        """
        Return all lines of the face.

        OUTPUT:

        A tuple of lines.

        EXAMPLES::

            sage: p = Polyhedron(rays = [[1,0],[-1,0],[0,1],[1,1]], vertices = [[-2,-2],[2,3]])
            sage: p.lines()
            (A line in the direction (1, 0),)
        """
        return tuple(self.line_generator())

    def __richcmp__(self, other, op):
        """
        Compare ``self`` and ``other``.

        INPUT:

        - ``other`` -- anything.

        OUTPUT:

        Two faces test equal if and only if they are faces of the same
        (not just isomorphic) polyhedron and their generators have the
        same indices.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: f = square.faces(1)
            sage: matrix(4,4, lambda i,j: ZZ(f[i] <= f[j]))
            [1 1 1 1]
            [0 1 1 1]
            [0 0 1 1]
            [0 0 0 1]
            sage: matrix(4,4, lambda i,j: ZZ(f[i] == f[j])) == 1
            True
        """
        if not isinstance(other, PolyhedronFace):
            return NotImplemented
        if self._polyhedron is not other._polyhedron:
            return NotImplemented
        return richcmp(self._ambient_Vrepresentation_indices,
                       other._ambient_Vrepresentation_indices, op)

    def ambient_Hrepresentation(self, index=None):
        r"""
        Return the H-representation objects of the ambient polytope
        defining the face.

        INPUT:

        - ``index`` -- optional. Either an integer or ``None``
          (default).

        OUTPUT:

        If the optional argument is not present, a tuple of
        H-representation objects. Each entry is either an inequality
        or an equation.

        If the optional integer ``index`` is specified, the
        ``index``-th element of the tuple is returned.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: for face in square.face_lattice():
            ....:     print(face.ambient_Hrepresentation())
            (An inequality (1, 0) x + 1 >= 0, An inequality (0, 1) x + 1 >= 0,
             An inequality (-1, 0) x + 1 >= 0, An inequality (0, -1) x + 1 >= 0)
            (An inequality (1, 0) x + 1 >= 0, An inequality (0, 1) x + 1 >= 0)
            (An inequality (1, 0) x + 1 >= 0, An inequality (0, -1) x + 1 >= 0)
            (An inequality (0, 1) x + 1 >= 0, An inequality (-1, 0) x + 1 >= 0)
            (An inequality (-1, 0) x + 1 >= 0, An inequality (0, -1) x + 1 >= 0)
            (An inequality (1, 0) x + 1 >= 0,)
            (An inequality (0, 1) x + 1 >= 0,)
            (An inequality (-1, 0) x + 1 >= 0,)
            (An inequality (0, -1) x + 1 >= 0,)
            ()
        """
        if index is None:
            return self._ambient_Hrepresentation
        else:
            return self._ambient_Hrepresentation[index]

    def ambient_Vrepresentation(self, index=None):
        r"""
        Return the V-representation objects of the ambient polytope
        defining the face.

        INPUT:

        - ``index`` -- optional. Either an integer or ``None``
          (default).

        OUTPUT:

        If the optional argument is not present, a tuple of
        V-representation objects. Each entry is either a vertex, a
        ray, or a line.

        If the optional integer ``index`` is specified, the
        ``index``-th element of the tuple is returned.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: for fl in square.face_lattice():
            ....:     print(fl.ambient_Vrepresentation())
            ()
            (A vertex at (-1, -1),)
            (A vertex at (-1, 1),)
            (A vertex at (1, -1),)
            (A vertex at (1, 1),)
            (A vertex at (-1, -1), A vertex at (-1, 1))
            (A vertex at (-1, -1), A vertex at (1, -1))
            (A vertex at (1, -1), A vertex at (1, 1))
            (A vertex at (-1, 1), A vertex at (1, 1))
            (A vertex at (-1, -1), A vertex at (-1, 1),
             A vertex at (1, -1), A vertex at (1, 1))
        """
        if index is None:
            return self._ambient_Vrepresentation
        else:
            return self._ambient_Vrepresentation[index]

    def n_ambient_Hrepresentation(self):
        """
        Return the number of objects that make up the ambient
        H-representation of the polyhedron.

        See also :meth:`ambient_Hrepresentation`.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(4)
            sage: face = p.face_lattice()[10]
            sage: face
            <0,2>
            sage: face.ambient_Hrepresentation()
            (An inequality (1, -1, 1, -1) x + 1 >= 0,
             An inequality (1, 1, 1, 1) x + 1 >= 0,
             An inequality (1, 1, 1, -1) x + 1 >= 0,
             An inequality (1, -1, 1, 1) x + 1 >= 0)
            sage: face.n_ambient_Hrepresentation()
            4
        """
        return len(self.ambient_Hrepresentation())

    def n_ambient_Vrepresentation(self):
        """
        Return the number of objects that make up the ambient
        V-representation of the polyhedron.

        See also :meth:`ambient_Vrepresentation`.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(4)
            sage: face = p.face_lattice()[10]
            sage: face
            <0,2>
            sage: face.ambient_Vrepresentation()
            (A vertex at (-1, 0, 0, 0), A vertex at (0, 0, -1, 0))
            sage: face.n_ambient_Vrepresentation()
            2
        """
        return len(self.ambient_Vrepresentation())

    def ambient_dim(self):
        r"""
        Return the dimension of the containing polyhedron.

        EXAMPLES::

            sage: P = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: face = P.faces(1)[0]
            sage: face.ambient_dim()
            4
        """
        return self._polyhedron.ambient_dim()

    @cached_method
    def dim(self):
        """
        Return the dimension of the face.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: fl = polytopes.dodecahedron().face_lattice()
            sage: [ x.dim() for x in fl ]
            [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3]
        """
        if self.n_ambient_Vrepresentation()==0:
            return -1
        else:
            origin = vector(self.ambient_Vrepresentation(0))
            v_list = [ vector(v)-origin for v in self.ambient_Vrepresentation() ]
            return matrix(v_list).rank()

    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT:

        A string listing the V-representation indices of the face.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: a_face = list( square.face_lattice() )[8]
            sage: a_face.__repr__()
            '<1,3>'
        """
        s = '<'
        s += ','.join([ str(v.index()) for v in self.ambient_Vrepresentation() ])
        s += '>'
        return s

    def polyhedron(self):
        """
        Return the containing polyhedron.

        EXAMPLES::

            sage: P = polytopes.cross_polytope(3); P
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
           sage: face = P.faces(2)[3]
            sage: face
            <0,3,4>
            sage: face.polyhedron()
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
        """
        return self._polyhedron

    @cached_method
    def as_polyhedron(self):
        """
        Return the face as an independent polyhedron.

        OUTPUT:

        A polyhedron.

        EXAMPLES::

            sage: P = polytopes.cross_polytope(3);  P
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
            sage: face = P.faces(2)[3]
            sage: face
            <0,3,4>
            sage: face.as_polyhedron()
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices

            sage: P.intersection(face.as_polyhedron()) == face.as_polyhedron()
            True
        """
        P = self._polyhedron
        parent = P.parent()
        Vrep = (self.vertices(), self.rays(), self.lines())
        return P.__class__(parent, Vrep, None)

    def normal_vector(self, algorithm='via_affine_hull', **kwds):
        r"""
        Compute the in-ward normal vector of this face of codimension 1.

        INPUT:

        - ``algorithm`` -- a string

          The following algorithms are available:

          * ``algorithm=via_affine_hull`` (default): Compute the normal vector by
            doing an orthonormal projection on the affine hull. The resulting vector
            will be normalized therefore.
          * ``algorithm=via_equation_system``: Compute the normal vector by
            solving a system of equations. This algorithm keeps rationality of the
            coordinates and the resulting vector is not normalized.

        OUTPUT:

        A vector.

        EXAMPLES::

            sage: [(face, face.normal_vector(algorithm='via_equation_system'))
            ....:  for face in polytopes.cube().faces(3-1)]
            [(<0,1,2,3>, (2, 0, 0)),
             (<0,1,4,5>, (0, 2, 0)),
             (<0,2,4,6>, (0, 0, 2)),
             (<1,3,5,7>, (0, 0, -2)),
             (<2,3,6,7>, (0, -2, 0)),
             (<4,5,6,7>, (-2, 0, 0))]
            sage: [(face, face.normal_vector(algorithm='via_equation_system'))
            ....:  for face in polytopes.simplex(2).faces(2-1)]
            [(<0,1>, (2, -1, -1)),
             (<0,2>, (-1, 2, -1)),
             (<1,2>, (-1, -1, 2))]
            sage: [(face, face.normal_vector(algorithm='via_equation_system'))
            ....:  for face in polytopes.simplex(3).faces(3-1)]
            [(<0,1,2>, (3, -1, -1, -1)),
             (<0,1,3>, (-1, 3, -1, -1)),
             (<0,2,3>, (-1, -1, 3, -1)),
             (<1,2,3>, (-1, -1, -1, 3))]

        ::

            sage: [(face, face.normal_vector(algorithm='via_affine_hull'))
            ....:  for face in polytopes.cube().faces(3-1)]
            [(<0,1,2,3>, (1, 0, 0)),
             (<0,1,4,5>, (0, 1, 0)),
             (<0,2,4,6>, (0, 0, 1)),
             (<1,3,5,7>, (0, 0, -1)),
             (<2,3,6,7>, (0, -1, 0)),
             (<4,5,6,7>, (-1, 0, 0))]
            sage: [(face, face.normal_vector(algorithm='via_affine_hull'))
            ....:  for face in polytopes.simplex(2).faces(2-1)]
            [(<0,1>, (0.6666666666666667?, -0.3333333333333334?, -0.3333333333333334?)),
             (<0,2>, (-0.3333333333333334?, 0.6666666666666667?, -0.3333333333333334?)),
             (<1,2>, (-0.3333333333333334?, -0.3333333333333334?, 0.6666666666666667?))]
            sage: [(face, face.normal_vector(algorithm='via_affine_hull'))
            ....:  for face in polytopes.simplex(3).faces(3-1)]
            [(<0,1,2>,
              (0.750000000000000?, -0.2500000000000000?, -0.2500000000000000?, -0.2500000000000000?)),
             (<0,1,3>,
              (-0.2500000000000000?, 0.750000000000000?, -0.2500000000000000?, -0.2500000000000000?)),
             (<0,2,3>,
              (-0.2500000000000000?, -0.2500000000000000?, 0.750000000000000?, -0.2500000000000000?)),
             (<1,2,3>,
              (-0.2500000000000000?, -0.2500000000000000?, -0.2500000000000000?, 0.750000000000000?))]

        TESTS::

            sage: polytopes.cube().faces(3)[0].normal_vector()
            Traceback (most recent call last):
            ...
            RuntimeError: only normal vectors of faces with codimension 1 can be computed
            sage: polytopes.cube().faces(1)[0].normal_vector()
            Traceback (most recent call last):
            ...
            RuntimeError: only normal vectors of faces with codimension 1 can be computed
            sage: polytopes.cube().faces(0)[0].normal_vector()
            Traceback (most recent call last):
            ...
            RuntimeError: only normal vectors of faces with codimension 1 can be computed
        """
        d = self.polyhedron().dim()
        if self.dim() != d-1:
            raise RuntimeError('only normal vectors of faces with codimension 1 '
                               'can be computed')

        if algorithm == 'via_equation_system':
            return self._normal_vector_via_equation_system_(**kwds)
        elif algorithm == 'via_affine_hull':
            return self._normal_vector_via_affine_hull_(**kwds)
        else:
            raise ValueError("unknown algorithm '{}'".format(algorithm))

    def _normal_vector_via_equation_system_(self):
        r"""
        Compute the in-ward normal vector of this face of codimension 1 by
        solving a system of equations.

        See :meth:`normal_vector` for details.

        OUTPUT:

        A vector (rational if the input was rational).

        TESTS::

            sage: [(face, face.normal_vector(algorithm='via_equation_system'))
            ....:  for face in polytopes.simplex(2).faces(2-1)]
            [(<0,1>, (2, -1, -1)),
             (<0,2>, (-1, 2, -1)),
             (<1,2>, (-1, -1, 2))]
        """
        from sage.functions.other import ceil
        from sage.matrix.constructor import Matrix
        from sage.rings.rational_field import QQ

        # preparation
        polyhedron = self.polyhedron()
        assert self.dim() == polyhedron.dim() - 1
        V_vectors = polyhedron.parametric_form()[1]
        V = Matrix(V_vectors).transpose()
        r = ceil(2*polyhedron.radius())

        # compute normal
        face_polyhedron = self.as_polyhedron()
        vectors = face_polyhedron.parametric_form()[1]
        W = Matrix(len(vectors), self.ambient_dim(), vectors)
        K = (W * V).right_kernel()
        assert K.dimension() == 1
        coefficients = K.gen()
        assert coefficients != K.zero()
        # we can now compute a normal vector
        n = sum(c*v for c, v in zip(coefficients, V_vectors))

        # determine the direction (in-/outward)
        m = r*n
        c = face_polyhedron.center()
        while True:
            if polyhedron.relative_interior_contains(c + m):
                return n
            if polyhedron.relative_interior_contains(c - m):
                return -n
            m = m / QQ(2)

    def _normal_vector_via_affine_hull_(self):
        r"""
        Compute the in-ward normal vector of this face of codimension 1 by
        doing an orthonormal projection on the affine hull.

        OUTPUT:

        A vector (normalized).

        TESTS::

            sage: [(face, face.normal_vector(algorithm='via_affine_hull'))
            ....:  for face in polytopes.simplex(2).faces(2-1)]
            [(<0,1>, (0.6666666666666667?, -0.3333333333333334?, -0.3333333333333334?)),
             (<0,2>, (-0.3333333333333334?, 0.6666666666666667?, -0.3333333333333334?)),
             (<1,2>, (-0.3333333333333334?, -0.3333333333333334?, 0.6666666666666667?))]
        """
        # preparation
        polyhedron = self.polyhedron()
        assert self.dim() == polyhedron.dim() - 1
        A = polyhedron.affine_hull(orthonormal=True, extend=True, as_affine_map=True)[0].matrix()

        def inequality(face):
            for H in face.ambient_Hrepresentation():
                if H.is_inequality():
                    return H

        # compute normal
        T = A * A.transpose()
        return T * inequality(self).vector()[1:]
