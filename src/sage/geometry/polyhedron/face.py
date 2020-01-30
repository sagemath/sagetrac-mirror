"""
A class to keep information about faces of a polyhedron

This module gives you a tool to work with the faces of a polyhedron
and their relative position. First, you need to find the faces. To get
the faces in a particular dimension, use the
:meth:`~sage.geometry.polyhedron.base.face` method::

    sage: P = polytopes.cross_polytope(3)
    sage: P.faces(3)
    (A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 6 vertices,)
    sage: [f.ambient_V_indices() for f in P.facets()]
    [(3, 4, 5),
     (2, 4, 5),
     (1, 3, 5),
     (1, 2, 5),
     (0, 3, 4),
     (0, 2, 4),
     (0, 1, 3),
     (0, 1, 2)]
    sage: [f.ambient_V_indices() for f in P.faces(1)]
    [(4, 5),
     (3, 5),
     (2, 5),
     (1, 5),
     (3, 4),
     (2, 4),
     (0, 4),
     (1, 3),
     (0, 3),
     (1, 2),
     (0, 2),
     (0, 1)]

or :meth:`~sage.geometry.polyhedron.base.face_lattice` to get the
whole face lattice as a poset::

    sage: P.face_lattice()
    Finite lattice containing 28 elements with distinguished linear extension

The faces are printed in shorthand notation where each integer is the
index of a vertex/ray/line in the same order as the containing
Polyhedron's :meth:`~sage.geometry.polyhedron.base.Vrepresentation` ::

    sage: face = P.faces(1)[8];  face
    A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices
    sage: face.ambient_V_indices()
    (0, 3)
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
        A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices
        sage: face.dim()
        2
        sage: face.ambient_V_indices()
        (0, 1, 2)
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
            A -1-dimensional face of a Polyhedron in ZZ^0
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
            sage: face = triangle.facets()[0]
            sage: for v in face.vertex_generator(): print(v)
            A vertex at (1, 0)
            A vertex at (1, 1)
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
            sage: face = triangle.faces(1)[2]
            sage: face.vertices()
            (A vertex at (0, 1), A vertex at (1, 0))
        """
        return tuple(self.vertex_generator())

    @cached_method
    def n_vertices(self):
        """
        Return the number of vertices of the face.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: Q = polytopes.cross_polytope(3)
            sage: face = Q.faces(2)[0]
            sage: face.n_vertices()
            3
        """
        return len(self.vertices())

    def ray_generator(self):
        """
        Return a generator for the rays of the face.

        EXAMPLES::

            sage: pi = Polyhedron(ieqs = [[1,1,0],[1,0,1]])
            sage: face = pi.faces(1)[1]
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
            sage: face = p.faces(2)[2]
            sage: face.rays()
            (A ray in the direction (1, 0, 0), A ray in the direction (0, 1, 0))
        """
        return tuple(self.ray_generator())

    @cached_method
    def n_rays(self):
        """
        Return the number of rays of the face.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0],[1,1,0,0]])
            sage: face = p.faces(2)[0]
            sage: face.n_rays()
            2
        """
        return len(self.rays())

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

    @cached_method
    def n_lines(self):
        """
        Return the number of lines of the face.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: p = Polyhedron(rays = [[1,0],[-1,0],[0,1],[1,1]], vertices = [[-2,-2],[2,3]])
            sage: p.n_lines()
            1
        """
        return len(self.lines())

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
            [1 1 0 0]
            [0 1 0 0]
            [1 1 1 0]
            [1 1 1 1]
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
            A 1-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 2 vertices
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
            A 1-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 2 vertices
            sage: face.ambient_Vrepresentation()
            (A vertex at (-1, 0, 0, 0), A vertex at (0, 0, -1, 0))
            sage: face.n_ambient_Vrepresentation()
            2
        """
        return len(self.ambient_Vrepresentation())

    def ambient_H_indices(self):
        """
        Return the indices of the H-representation objects of the
        ambient polyhedron that make up the H-representation of ``self``.

        See also :meth:`ambient_Hrepresentation`.

        OUTPUT:

        Tuple of indices

        EXAMPLES::

            sage: Q = polytopes.cross_polytope(3)
            sage: F = Q.faces(1)
            sage: [f.ambient_H_indices() for f in F]
            [(4, 5),
             (5, 6),
             (4, 7),
             (6, 7),
             (0, 5),
             (3, 4),
             (0, 3),
             (1, 6),
             (0, 1),
             (2, 7),
             (2, 3),
             (1, 2)]
        """
        return self._ambient_Hrepresentation_indices

    def ambient_V_indices(self):
        """
        Return the indices of the V-representation objects of the
        ambient polyhedron that make up the V-representation of ``self``.

        See also :meth:`ambient_Vrepresentation`.

        OUTPUT:

        Tuple of indices

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: F = P.faces(2)
            sage: [f.ambient_V_indices() for f in F]
            [(0, 1, 4, 5),
             (0, 2, 4, 6),
             (0, 1, 2, 3),
             (4, 5, 6, 7),
             (2, 3, 6, 7),
             (1, 3, 5, 7)]
        """
        return self._ambient_Vrepresentation_indices

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

        TESTS:

        Check that :trac:`28650` is fixed::

            sage: P = Polyhedron(vertices=[[1,0]], rays=[[1,0],[0,1]])
            sage: P.faces(2)
            (A 2-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,)
        """
        if self.n_ambient_Vrepresentation() == 0:
            return -1
        else:
            origin = self.vertices()[0].vector()
            v_list = [vector(v) - origin for v in
                     self.ambient_Vrepresentation() if v.is_vertex()]
            v_list += [vector(v) for v in self.ambient_Vrepresentation()
                      if v.is_ray() or v.is_line()]
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
            'A 1-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 2 vertices'
        """
        desc = ''
        desc += 'A ' + repr(self.dim()) + '-dimensional face'
        desc += ' of a Polyhedron in '
        desc += self.polyhedron().parent()._repr_ambient_module()

        if self.n_vertices() > 0:
            desc += ' defined as the convex hull of '
            desc += repr(self.n_vertices())
            if self.n_vertices() == 1: desc += ' vertex'
            else:                      desc += ' vertices'

            if self.n_rays() > 0:
                if self.n_lines() > 0: desc += ", "
                else:                  desc += " and "
                desc += repr(self.n_rays())
                if self.n_rays() == 1: desc += ' ray'
                else:                  desc += ' rays'

            if self.n_lines() > 0:
                if self.n_rays() > 0: desc += ", "
                else:                 desc += " and "
                desc += repr(self.n_lines())
                if self.n_lines() == 1: desc += ' line'
                else:                   desc += ' lines'

        return desc

    def polyhedron(self):
        """
        Return the containing polyhedron.

        EXAMPLES::

            sage: P = polytopes.cross_polytope(3); P
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
            sage: face = P.facets()[3]
            sage: face
            A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices
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
            A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: face.as_polyhedron()
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices

            sage: P.intersection(face.as_polyhedron()) == face.as_polyhedron()
            True
        """
        P = self._polyhedron
        parent = P.parent()
        Vrep = (self.vertices(), self.rays(), self.lines())
        return P.__class__(parent, Vrep, None)

def combinatorial_face_to_polyhedral_face(polyhedron, combinatorial_face):
    r"""
    Convert a combinatorial face to a face of a polyhedron.

    INPUT:

    - ``polyhedron`` -- a polyhedron containing ``combinatorial_face``
    - ``combinatorial_face`` -- a ``CombinatorialFace``

    OUTPUT: a ``PolyhedronFace``.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.face import combinatorial_face_to_polyhedral_face
        sage: P = polytopes.simplex()
        sage: C = P.combinatorial_polyhedron()
        sage: it = C.face_iter()
        sage: comb_face = next(it)
        sage: combinatorial_face_to_polyhedral_face(P, comb_face)
        A 2-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 3 vertices

    TESTS:

    Making sure that backends do not change their order of inequalites/equations
    without applying the changes to this method::

        sage: polytopes.simplex(backend='field').equations()[0].index()
        4
        sage: polytopes.simplex(backend='ppl').equations()[0].index()
        0
        sage: polytopes.simplex(backend='cdd').equations()[0].index()
        4
        sage: polytopes.simplex(backend='normaliz').equations()[0].index() # optional - pynormaliz
        4
        sage: polytopes.simplex(backend='polymake').equations()[0].index() # optional - polymake
        4
    """
    V_indices = combinatorial_face.ambient_V_indices()
    n_equations = polyhedron.n_equations()

    if polyhedron.backend() in ('ppl',):
        # Equations before inequalities in Hrep.
        H_indices = tuple(range(n_equations))
        H_indices += tuple(x+n_equations for x in combinatorial_face.ambient_H_indices())
    elif polyhedron.backend() in ('normaliz', 'cdd', 'field', 'polymake'):
        # Equations after the inequalities in Hrep.
        n_ieqs = polyhedron.n_inequalities()
        H_indices = tuple(range(n_ieqs, n_ieqs + n_equations))
        H_indices += tuple(x for x in combinatorial_face.ambient_H_indices())
    else:
        raise NotImplementedError("unknown backend")

    return PolyhedronFace(polyhedron, V_indices, H_indices)
