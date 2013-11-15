"""
Fast Lattice Polytopes using PPL.

The :func:`LatticePolytope_PPL` class is a thin wrapper around PPL
polyhedra. Its main purpose is to be fast to construct, at the cost of
being much less full-featured than the usual polyhedra. This makes it
possible to iterate with it over the list of all 473800776 reflexive
polytopes in 4 dimensions.

.. NOTE::

    For general lattice polyhedra you should use
    :func:`~sage.geometry.polyhedon.constructor.Polyhedron` with
    `base_ring=ZZ`.

The class derives from the PPL :class:`sage.libs.ppl.C_Polyhedron`
class, so you can work with the underlying generator and constraint
objects. However, integral points are generally represented by
`\ZZ`-vectors. In the following, we always use *generator* to refer
the PPL generator objects and *vertex* (or integral point) for the
corresponding `\ZZ`-vector.

EXAMPLES::

    sage: vertices = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (-9, -6, -1, -1)]
    sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
    sage: P = LatticePolytope_PPL(vertices);  P
    A 4-dimensional lattice polytope in ZZ^4 with 5 vertices
    sage: P.integral_points()
    ((-9, -6, -1, -1), (-3, -2, 0, 0), (-2, -1, 0, 0), (-1, -1, 0, 0),
     (-1, 0, 0, 0), (0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0))
    sage: P.integral_points_not_interior_to_facets()
    ((-9, -6, -1, -1), (-3, -2, 0, 0), (0, 0, 0, 0), (1, 0, 0, 0),
     (0, 1, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0))

Fibrations of the lattice polytopes are defined as lattice
sub-polytopes and give rise to fibrations of toric varieties for
suitable fan refinements. We can compute them using
:meth:`~LatticePolytope_PPL.fibration_generator` ::

    sage: F = P.fibration_generator(2).next()
    sage: F.vertices()
    ((1, 0, 0, 0), (0, 1, 0, 0), (-3, -2, 0, 0))

Finally, we can compute automorphisms and identify fibrations that
only differ by a lattice automorphism::

    sage: square = LatticePolytope_PPL((-1,-1),(-1,1),(1,-1),(1,1))
    sage: fibers = [ f.vertices() for f in square.fibration_generator(1) ];  fibers
    [((1, 0), (-1, 0)), ((0, 1), (0, -1)), ((-1, -1), (1, 1)), ((-1, 1), (1, -1))]
    sage: square.pointsets_mod_automorphism(fibers)
    (frozenset([(0, 1), (0, -1)]), frozenset([(1, 1), (-1, -1)]))

AUTHORS:

    - Volker Braun: initial version, 2012
"""

########################################################################
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

import copy
from sage.rings.integer import GCD_list
from sage.rings.all import ZZ, QQ
from sage.misc.all import union, cached_method, cached_function, prod, uniq
from sage.modules.all import (
    vector, zero_vector )
from sage.matrix.constructor import (
    matrix, column_matrix, diagonal_matrix )
from sage.libs.ppl import (
     C_Polyhedron, Linear_Expression, Variable,
    point, ray, line,
    Generator, Generator_System, Generator_System_iterator )
from sage.libs.ppl import (
    C_Polyhedron, Linear_Expression, Variable,
    point, ray, line, Generator, Generator_System,
    Constraint_System,
    Poly_Con_Relation )
from sage.geometry.polyhedron.lattice_euclidean_group_element import (
    LatticePolytopeNoEmbeddingError, LatticePolytopesNotIsomorphicError )




########################################################################
def _class_for_LatticePolytope(dim):
    """
    Return the appropriate class in the given dimension.

    Helper function for :func:`LatticePolytope_PPL`. You should not
    have to use this function manually.

    INPUT:

    - ``dim`` -- integer. The ambient space dimenson.

    OUTPUT:

    The appropriate class for the lattice polytope.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import _class_for_LatticePolytope
        sage: _class_for_LatticePolytope(2)
        <class 'sage.geometry.polyhedron.ppl_lattice_polygon.LatticePolygon_PPL_class'>
        sage: _class_for_LatticePolytope(3)
        <class 'sage.geometry.polyhedron.ppl_lattice_polytope.LatticePolytope_PPL_class'>
    """
    if dim <= 2:
        from sage.geometry.polyhedron.ppl_lattice_polygon import LatticePolygon_PPL_class
        return LatticePolygon_PPL_class
    return LatticePolytope_PPL_class


########################################################################
def LatticePolytope_PPL(*args):
    """
    Construct a new instance of the PPL-based lattice polytope class.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
        sage: LatticePolytope_PPL((0,0),(1,0),(0,1))
        A 2-dimensional lattice polytope in ZZ^2 with 3 vertices

        sage: from sage.libs.ppl import point, Generator_System, C_Polyhedron, Linear_Expression, Variable
        sage: p = point(Linear_Expression([2,3],0));  p
        point(2/1, 3/1)
        sage: LatticePolytope_PPL(p)
        A 0-dimensional lattice polytope in ZZ^2 with 1 vertex

        sage: P = C_Polyhedron(Generator_System(p));  P
        A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point
        sage: LatticePolytope_PPL(P)
        A 0-dimensional lattice polytope in ZZ^2 with 1 vertex

    A ``TypeError`` is raised if the arguments do not specify a lattice polytope::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
        sage: LatticePolytope_PPL((0,0),(1/2,1))
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

        sage: from sage.libs.ppl import point, Generator_System, C_Polyhedron, Linear_Expression, Variable
        sage: p = point(Linear_Expression([2,3],0), 5);  p
        point(2/5, 3/5)
        sage: LatticePolytope_PPL(p)
        Traceback (most recent call last):
         ...
        TypeError: generator is not a lattice polytope generator

        sage: P = C_Polyhedron(Generator_System(p));  P
        A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point
        sage: LatticePolytope_PPL(P)
        Traceback (most recent call last):
        ...
        TypeError: polyhedron has non-integral generators
    """
    polytope_class = LatticePolytope_PPL_class
    if len(args)==1 and isinstance(args[0], C_Polyhedron):
        polyhedron = args[0]
        polytope_class = _class_for_LatticePolytope(polyhedron.space_dimension())
        if not all(p.is_point() and p.divisor().is_one() for p in polyhedron.generators()):
            raise TypeError('polyhedron has non-integral generators')
        return polytope_class(polyhedron)
    if len(args)==1 \
            and isinstance(args[0], (list, tuple)) \
            and isinstance(args[0][0], (list,tuple)):
        vertices = args[0]
    else:
        vertices = args
    gs = Generator_System()
    for v in vertices:
        if isinstance(v, Generator):
            if (not v.is_point()) or (not v.divisor().is_one()):
                raise TypeError('generator is not a lattice polytope generator')
            gs.insert(v)
        else:
            gs.insert(point(Linear_Expression(v, 0)))
    if not gs.empty():
        dim = Generator_System_iterator(gs).next().space_dimension()
        polytope_class = _class_for_LatticePolytope(dim)
    return polytope_class(gs)



########################################################################
class LatticePolytope_PPL_class(C_Polyhedron):
    """
    The lattice polytope class.

    You should use :func:LatticePolytope_PPL` to construct instances.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
        sage: LatticePolytope_PPL((0,0),(1,0),(0,1))
        A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
    """

    def _repr_(self):
        """
        Return the string representation

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: P = LatticePolytope_PPL((0,0),(1,0),(0,1))
            sage: P
            A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
            sage: P._repr_()
            'A 2-dimensional lattice polytope in ZZ^2 with 3 vertices'

            sage: LatticePolytope_PPL()
            The empty lattice polytope in ZZ^0
        """
        desc = ''
        if self.n_vertices()==0:
            desc += 'The empty lattice polytope'
        else:
            desc += 'A ' + repr(self.affine_dimension()) + '-dimensional lattice polytope'
        desc += ' in ZZ^' + repr(self.space_dimension())

        if self.n_vertices()>0:
            desc += ' with '
            desc += repr(self.n_vertices())
            if self.n_vertices()==1: desc += ' vertex'
            else:                    desc += ' vertices'
        return desc



    def is_bounded(self):
        """
        Return whether the lattice polytope is compact.

        OUTPUT:

        Always ``True``, since polytopes are by definition compact.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((0,0),(1,0),(0,1)).is_bounded()
            True
        """
        return True

    @cached_method
    def n_vertices(self):
        """
        Return the number of vertices.

        OUTPUT:

        An integer, the number of vertices.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((0,0,0), (1,0,0), (0,1,0)).n_vertices()
            3
        """
        return len(self.minimized_generators())

    @cached_method
    def is_simplex(self):
        r"""
        Return whether the polyhedron is a simplex.

        OUTPUT:

        Boolean, whether the polyhedron is a simplex (possibly of
        strictly smaller dimension than the ambient space).

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((0,0,0), (1,0,0), (0,1,0)).is_simplex()
            True
        """
        return self.affine_dimension()+1==self.n_vertices()

    @cached_method
    def bounding_box(self):
        r"""
        Return the coordinates of a rectangular box containing the non-empty polytope.

        OUTPUT:

        A pair of tuples ``(box_min, box_max)`` where ``box_min`` are
        the coordinates of a point bounding the coordinates of the
        polytope from below and ``box_max`` bounds the coordinates
        from above.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((0,0),(1,0),(0,1)).bounding_box()
            ((0, 0), (1, 1))
        """
        box_min = []
        box_max = []
        if self.is_empty():
            raise ValueError('empty polytope is not allowed')
        for i in range(0, self.space_dimension()):
            x = Variable(i)
            coords = [ v.coefficient(x) for v in self.generators() ]
            max_coord = max(coords)
            min_coord = min(coords)
            box_max.append(max_coord)
            box_min.append(min_coord)
        return (tuple(box_min), tuple(box_max))

    @cached_method
    def n_integral_points(self):
        """
        Return the number of integral points.

        OUTPUT:

        Integer. The number of integral points contained in the
        lattice polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((0,0),(1,0),(0,1)).n_integral_points()
            3
        """
        if self.is_empty():
            return tuple()
        box_min, box_max = self.bounding_box()
        from sage.geometry.integral_points import rectangular_box_points
        return rectangular_box_points(box_min, box_max, self, count_only=True)

    @cached_method
    def integral_points(self):
        r"""
        Return the integral points in the polyhedron.

        Uses the naive algorithm (iterate over a rectangular bounding
        box).

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a ``ValueError`` is raised.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((-1,-1),(1,0),(1,1),(0,1)).integral_points()
            ((-1, -1), (0, 0), (0, 1), (1, 0), (1, 1))

            sage: simplex = LatticePolytope_PPL((1,2,3), (2,3,7), (-2,-3,-11))
            sage: simplex.integral_points()
            ((-2, -3, -11), (0, 0, -2), (1, 2, 3), (2, 3, 7))

        The polyhedron need not be full-dimensional::

            sage: simplex = LatticePolytope_PPL((1,2,3,5), (2,3,7,5), (-2,-3,-11,5))
            sage: simplex.integral_points()
            ((-2, -3, -11, 5), (0, 0, -2, 5), (1, 2, 3, 5), (2, 3, 7, 5))

            sage: point = LatticePolytope_PPL((2,3,7))
            sage: point.integral_points()
            ((2, 3, 7),)

            sage: empty = LatticePolytope_PPL()
            sage: empty.integral_points()
            ()

        Here is a simplex where the naive algorithm of running over
        all points in a rectangular bounding box no longer works fast
        enough::

            sage: v = [(1,0,7,-1), (-2,-2,4,-3), (-1,-1,-1,4), (2,9,0,-5), (-2,-1,5,1)]
            sage: simplex = LatticePolytope_PPL(v); simplex
            A 4-dimensional lattice polytope in ZZ^4 with 5 vertices
            sage: len(simplex.integral_points())
            49

        Finally, the 3-d reflexive polytope number 4078::

            sage: v = [(1,0,0), (0,1,0), (0,0,1), (0,0,-1), (0,-2,1),
            ...        (-1,2,-1), (-1,2,-2), (-1,1,-2), (-1,-1,2), (-1,-3,2)]
            sage: P = LatticePolytope_PPL(*v)
            sage: pts1 = P.integral_points()                     # Sage's own code
            sage: pts2 = LatticePolytope(v).points().columns()   # PALP
            sage: for p in pts1: p.set_immutable()
            sage: for p in pts2: p.set_immutable()
            sage: set(pts1) == set(pts2)
            True

            sage: timeit('Polyhedron(v).integral_points()')   # random output
            sage: timeit('LatticePolytope(v).points()')       # random output
            sage: timeit('LatticePolytope_PPL(*v).integral_points()')       # random output
        """
        if self.is_empty():
            return tuple()
        box_min, box_max = self.bounding_box()
        from sage.geometry.integral_points import rectangular_box_points
        points = rectangular_box_points(box_min, box_max, self)
        if not self.n_integral_points.is_in_cache():
            self.n_integral_points.set_cache(len(points))
        return points

    @cached_method
    def _integral_points_saturating(self):
        """
        Return the integral points together with information about
        which inequalities are saturated.

        See :func:`~sage.geometry.integral_points.rectangular_box_points`.

        OUTPUT:

        A tuple of pairs (one for each integral point) consisting of a
        pair ``(point, Hrep)``, where ``point`` is the coordinate
        vector of the intgeral point and ``Hrep`` is the set of
        indices of the :meth:`minimized_constraints` that are
        saturated at the point.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: quad = LatticePolytope_PPL((-1,-1),(0,1),(1,0),(1,1))
            sage: quad._integral_points_saturating()
            (((-1, -1), frozenset([0, 1])),
             ((0, 0), frozenset([])),
             ((0, 1), frozenset([0, 3])),
             ((1, 0), frozenset([1, 2])),
             ((1, 1), frozenset([2, 3])))
        """
        if self.is_empty():
            return tuple()
        box_min, box_max = self.bounding_box()
        from sage.geometry.integral_points import rectangular_box_points
        points= rectangular_box_points(box_min, box_max, self, return_saturated=True)
        if not self.n_integral_points.is_in_cache():
            self.n_integral_points.set_cache(len(points))
        if not self.integral_points.is_in_cache():
            self.integral_points.set_cache(tuple(p[0] for p in points))
        return points

    @cached_method
    def integral_points_not_interior_to_facets(self):
        """
        Return the integral points not interior to facets

        OUTPUT:

        A tuple whose entries are the coordinate vectors of integral
        points not interior to facets (codimension one faces) of the
        lattice polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: square = LatticePolytope_PPL((-1,-1),(-1,1),(1,-1),(1,1))
            sage: square.n_integral_points()
            9
            sage: square.integral_points_not_interior_to_facets()
            ((-1, -1), (-1, 1), (0, 0), (1, -1), (1, 1))
        """
        n = 1 + self.space_dimension() - self.affine_dimension()
        return tuple(p[0] for p in self._integral_points_saturating() if len(p[1])!=n)

    @cached_method
    def vertices(self):
        r"""
        Return the vertices as a tuple of `\ZZ`-vectors.

        OUTPUT:

        A tuple of `\ZZ`-vectors. Each entry is the coordinate vector
        of an integral points of the lattice polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((-9,-6,-1,-1),(0,0,0,1),(0,0,1,0),(0,1,0,0),(1,0,0,0))
            sage: p.vertices()
            ((-9, -6, -1, -1), (0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, 0), (1, 0, 0, 0))
            sage: p.minimized_generators()
            Generator_System {point(-9/1, -6/1, -1/1, -1/1), point(0/1, 0/1, 0/1, 1/1),
            point(0/1, 0/1, 1/1, 0/1), point(0/1, 1/1, 0/1, 0/1), point(1/1, 0/1, 0/1, 0/1)}
        """
        d = self.space_dimension()
        v = vector(ZZ, d)
        points = []
        for g in self.minimized_generators():
            for i in range(0,d):
                v[i] = g.coefficient(Variable(i))
            v_copy = copy.copy(v)
            v_copy.set_immutable()
            points.append(v_copy)
        return tuple(points)

    def vertices_saturating(self, constraint):
        """
        Return the vertices saturating the constraint

        INPUT:

        - ``constraint`` -- a constraint (inequality or equation) of
          the polytope.

        OUTPUT:

        The tuple of vertices saturating the constraint. The vertices
        are returned as `\ZZ`-vectors, as in :meth:`vertices`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((0,0),(0,1),(1,0))
            sage: ieq = iter(p.constraints()).next();  ieq
            x0>=0
            sage: p.vertices_saturating(ieq)
            ((0, 0), (0, 1))
        """
        from sage.libs.ppl import C_Polyhedron, Poly_Con_Relation
        result = []
        for i,v in enumerate(self.minimized_generators()):
            v = C_Polyhedron(v)
            if v.relation_with(constraint).implies(Poly_Con_Relation.saturates()):
                result.append(self.vertices()[i])
        return tuple(result)

    @cached_method
    def is_full_dimensional(self):
        """
        Return whether the lattice polytope is full dimensional.

        OUTPUT:

        Boolean. Whether the :meth:`affine_dimension` equals the
        ambient space dimension.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((0,0),(0,1))
            sage: p.is_full_dimensional()
            False
            sage: q = LatticePolytope_PPL((0,0),(0,1),(1,0))
            sage: q.is_full_dimensional()
            True
        """

        return self.affine_dimension() == self.space_dimension()

    def fibration_generator(self, dim):
        """
        Generate the lattice polytope fibrations.

        For the purposes of this function, a lattice polytope fiber is
        a sub-lattice polytope. Projecting the plane spanned by the
        subpolytope to a point yields another lattice polytope, the
        base of the fibration.

        INPUT:

        - ``dim`` -- integer. The dimension of the lattice polytope
          fiber.

        OUTPUT:

        A generator yielding the distinct lattice polytope fibers of
        given dimension.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((-9,-6,-1,-1),(0,0,0,1),(0,0,1,0),(0,1,0,0),(1,0,0,0))
            sage: list( p.fibration_generator(2) )
            [A 2-dimensional lattice polytope in ZZ^4 with 3 vertices]
        """
        assert self.is_full_dimensional()
        codim = self.space_dimension() - dim
        # "points" are the potential vertices of the fiber. They are
        # in the $codim$-skeleton of the polytope, which is contained
        # in the points that saturate at least $dim$ equations.
        points = [ p for p in self._integral_points_saturating() if len(p[1])>=dim ]
        points = sorted(points, key=lambda x:len(x[1]))

        # iterate over point combinations subject to all points being on one facet.
        def point_combinations_iterator(n, i0=0, saturated=None):
            for i in range(i0, len(points)):
                p, ieqs = points[i]
                if saturated is None:
                    saturated_ieqs = ieqs
                else:
                    saturated_ieqs = saturated.intersection(ieqs)
                if len(saturated_ieqs)==0:
                    continue
                if n == 1:
                    yield [i]
                else:
                    for c in point_combinations_iterator(n-1, i+1, saturated_ieqs):
                        yield [i] + c

        point_lines = [ line(Linear_Expression(p[0].list(),0)) for p in points ]
        origin = point()
        fibers = set()
        gs = Generator_System()
        for indices in point_combinations_iterator(dim):
            gs.clear()
            gs.insert(origin)
            for i in indices:
                gs.insert(point_lines[i])
            plane = C_Polyhedron(gs)
            if plane.affine_dimension() != dim:
                continue
            plane.intersection_assign(self)
            if (not self.is_full_dimensional()) and (plane.affine_dimension() != dim):
                continue
            try:
                fiber = LatticePolytope_PPL(plane)
            except TypeError:   # not a lattice polytope
                continue
            fiber_vertices = tuple(sorted(fiber.vertices()))
            if fiber_vertices not in fibers:
                yield fiber
                fibers.update([fiber_vertices])

    def pointsets_mod_automorphism(self, pointsets):
        """
        Return ``pointsets`` modulo the automorphisms of ``self``.

        INPUT:

        - ``polytopes`` a tuple/list/iterable of subsets of the
          integral points of ``self``.

        OUTPUT:

        Representatives of the point sets modulo the
        :meth:`lattice_automorphism_group` of ``self``.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: square = LatticePolytope_PPL((-1,-1),(-1,1),(1,-1),(1,1))
            sage: fibers = [ f.vertices() for f in square.fibration_generator(1) ]
            sage: square.pointsets_mod_automorphism(fibers)
            (frozenset([(0, 1), (0, -1)]), frozenset([(1, 1), (-1, -1)]))

            sage: cell24 = LatticePolytope_PPL(
            ...   (1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(1,-1,-1,1),(0,0,-1,1),
            ...   (0,-1,0,1),(-1,0,0,1),(1,0,0,-1),(0,1,0,-1),(0,0,1,-1),(-1,1,1,-1),
            ...   (1,-1,-1,0),(0,0,-1,0),(0,-1,0,0),(-1,0,0,0),(1,-1,0,0),(1,0,-1,0),
            ...   (0,1,1,-1),(-1,1,1,0),(-1,1,0,0),(-1,0,1,0),(0,-1,-1,1),(0,0,0,-1))
            sage: fibers = [ f.vertices() for f in cell24.fibration_generator(2) ]
            sage: cell24.pointsets_mod_automorphism(fibers)   # long time
            (frozenset([(1, 0, -1, 0), (-1, 0, 1, 0), (0, -1, -1, 1), (0, 1, 1, -1)]),
             frozenset([(-1, 0, 0, 0), (1, 0, 0, 0), (0, 0, 0, 1),
                        (1, 0, 0, -1), (0, 0, 0, -1), (-1, 0, 0, 1)]))
        """
        points = set()
        for ps in pointsets:
            points.update(ps)
        points = tuple(points)
        Aut = self.lattice_automorphism_group(points,
                                              point_labels=tuple(range(len(points))))
        indexsets = set([ frozenset([points.index(p) for p in ps]) for ps in pointsets ])
        orbits = []
        while len(indexsets)>0:
            idx = indexsets.pop()
            orbits.append(frozenset([points[i] for i in idx]))
            for g in Aut:
                g_idx = frozenset([g(i) for i in idx])
                indexsets.difference_update([g_idx])
        return tuple(orbits)

    @cached_method
    def ambient_space(self):
        r"""
        Return the ambient space.

        OUTPUT:

        The free module `\ZZ^d`, where `d` is the ambient space
        dimension.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: point = LatticePolytope_PPL((1,2,3))
            sage: point.ambient_space()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        from sage.modules.free_module import FreeModule
        return FreeModule(ZZ, self.space_dimension())

    def contains(self, point_coordinates):
        r"""
        Test whether point is contained in the polytope.

        INPUT:

        - ``point_coordinates`` -- a list/tuple/iterable of rational
          numbers. The coordinates of the point.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: line = LatticePolytope_PPL((1,2,3), (-1,-2,-3))
            sage: line.contains([0,0,0])
            True
            sage: line.contains([1,0,0])
            False
        """
        p = C_Polyhedron(point(Linear_Expression(list(point_coordinates), 1)))
        is_included = Poly_Con_Relation.is_included()
        for c in self.constraints():
            if not p.relation_with(c).implies(is_included):
                return False
        return True

    @cached_method
    def contains_origin(self):
        """
        Test whether the polytope contains the origin

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((1,2,3), (-1,-2,-3)).contains_origin()
            True
            sage: LatticePolytope_PPL((1,2,5), (-1,-2,-3)).contains_origin()
            False
        """
        return self.contains(self.ambient_space().zero())

    @cached_method
    def affine_space(self):
        r"""
        Return the affine space spanned by the polytope.

        OUTPUT:

        The free module `\ZZ^n`, where `n` is the dimension of the
        affine space spanned by the points of the polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: point = LatticePolytope_PPL((1,2,3))
            sage: point.affine_space()
            Free module of degree 3 and rank 0 over Integer Ring
            Echelon basis matrix:
            []
            sage: line = LatticePolytope_PPL((1,1,1), (1,2,3))
            sage: line.affine_space()
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [0 1 2]
        """
        vertices = self.vertices()
        if not self.contains_origin():
            v0 = vertices[0]
            vertices = [v-v0 for v in vertices]
        return self.ambient_space().span(vertices).saturation()

    def affine_lattice_polytope(self):
        """
        Return the lattice polytope restricted to
        :meth:`affine_space`.

        OUTPUT:

        A new, full-dimensional lattice polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: poly_4d = LatticePolytope_PPL((-9,-6,0,0),(0,1,0,0),(1,0,0,0));  poly_4d
            A 2-dimensional lattice polytope in ZZ^4 with 3 vertices
            sage: poly_4d.space_dimension()
            4
            sage: poly_2d = poly_4d.affine_lattice_polytope();  poly_2d
            A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
            sage: poly_2d.space_dimension()
            2
        """
        V = self.affine_space()
        if self.contains_origin():
            vertices = [ V.coordinates(v) for v in self.vertices() ]
        else:
            v0 = vertices[0]
            vertices = [ V.coordinates(v-v0) for v in self.vertices() ]
        return LatticePolytope_PPL(*vertices)

    @cached_method
    def base_projection(self, fiber):
        """
        The projection that maps the sub-polytope ``fiber`` to a
        single point.

        OUTPUT:

        The quotient module of the ambient space modulo the
        :meth:`affine_space` spanned by the fiber.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: poly = LatticePolytope_PPL((-9,-6,-1,-1),(0,0,0,1),(0,0,1,0),(0,1,0,0),(1,0,0,0))
            sage: fiber = poly.fibration_generator(2).next()
            sage: poly.base_projection(fiber)
            Finitely generated module V/W over Integer Ring with invariants (0, 0)
        """
        return self.ambient_space().quotient(fiber.affine_space())

    @cached_method
    def base_projection_matrix(self, fiber):
        """
        The projection that maps the sub-polytope ``fiber`` to a
        single point.

        OUTPUT:

        An integer matrix that represents the projection to the
        base.

        .. SEEALSO::

            The :meth:`base_projection` yields equivalent information,
            and is easier to use. However, just returning the matrix
            has lower overhead.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: poly = LatticePolytope_PPL((-9,-6,-1,-1),(0,0,0,1),(0,0,1,0),(0,1,0,0),(1,0,0,0))
            sage: fiber = poly.fibration_generator(2).next()
            sage: poly.base_projection_matrix(fiber)
            [0 0 1 0]
            [0 0 0 1]

        Note that the basis choice in :meth:`base_projection` for the
        quotient is usually different::

            sage: proj = poly.base_projection(fiber)
            sage: proj_matrix = poly.base_projection_matrix(fiber)
            sage: [ proj(p) for p in poly.integral_points() ]
            [(-1, -1), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (0, 1)]
            sage: [ proj_matrix*p for p in poly.integral_points() ]
            [(-1, -1), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 1), (1, 0)]
        """
        return matrix(ZZ, fiber.vertices()).right_kernel_matrix()

    def base_rays(self, fiber, points):
        """
        Return the primitive lattice vectors that generate the
        direction given by the base projection of points.

        INPUT:

        - ``fiber`` -- a sub-lattice polytope defining the
          :meth:`base_projection`.

        - ``points`` -- the points to project to the base.

        OUTPUT:

        A tuple of primitive `\ZZ`-vectors.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: poly = LatticePolytope_PPL((-9,-6,-1,-1),(0,0,0,1),(0,0,1,0),(0,1,0,0),(1,0,0,0))
            sage: fiber = poly.fibration_generator(2).next()
            sage: poly.base_rays(fiber, poly.integral_points_not_interior_to_facets())
            ((-1, -1), (0, 1), (1, 0))

            sage: p = LatticePolytope_PPL((1,0),(1,2),(-1,0))
            sage: f = LatticePolytope_PPL((1,0),(-1,0))
            sage: p.base_rays(f, p.integral_points())
            ((1),)
        """
        quo = self.base_projection(fiber)
        vertices = []
        for p in points:
            v = vector(ZZ,quo(p))
            if v.is_zero():
                continue
            d =  GCD_list(v.list())
            if d>1:
                for i in range(0,v.degree()):
                    v[i] /= d
            v.set_immutable()
            vertices.append(v)
        return tuple(uniq(vertices))

    @cached_method
    def has_IP_property(self):
        """
        Whether the lattice polytope has the IP property.

        That is, the polytope is full-dimensional and the origin is a
        interior point not on the boundary.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((-1,-1),(0,1),(1,0)).has_IP_property()
            True
            sage: LatticePolytope_PPL((-1,-1),(1,1)).has_IP_property()
            False
        """
        origin = C_Polyhedron(point(0*Variable(self.space_dimension())))
        is_included = Poly_Con_Relation.is_included()
        saturates = Poly_Con_Relation.saturates()
        for c in self.constraints():
            rel = origin.relation_with(c)
            if (not rel.implies(is_included)) or rel.implies(saturates):
                return False
        return True

    @cached_method
    def restricted_automorphism_group(self, vertex_labels=None):
        r"""
        Return the restricted automorphism group.

        First, let the linear automorphism group be the subgroup of
        the Euclidean group `E(d) = GL(d,\RR) \ltimes \RR^d`
        preserving the `d`-dimensional polyhedron. The Euclidean group
        acts in the usual way `\vec{x}\mapsto A\vec{x}+b` on the
        ambient space. The restricted automorphism group is the
        subgroup of the linear automorphism group generated by
        permutations of vertices. If the polytope is full-dimensional,
        it is equal to the full (unrestricted) automorphism group.

        INPUT:

        - ``vertex_labels`` -- a tuple or ``None`` (default). The
          labels of the vertices that will be used in the output
          permutation group. By default, the vertices are used
          themselves.

        OUTPUT:

        A
        :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic>`
        acting on the vertices (or the ``vertex_labels``, if
        specified).

        REFERENCES:

        ..  [BSS]
            David Bremner, Mathieu Dutour Sikiric, Achill Schuermann:
            Polyhedral representation conversion up to symmetries.
            http://arxiv.org/abs/math/0702239

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: Z3square = LatticePolytope_PPL((0,0), (1,2), (2,1), (3,3))
            sage: Z3square.restricted_automorphism_group(vertex_labels=(1,2,3,4))
            Permutation Group with generators [(2,3), (1,2)(3,4), (1,4)]
            sage: G = Z3square.restricted_automorphism_group(); G
            Permutation Group with generators [((1,2),(2,1)),
            ((0,0),(1,2))((2,1),(3,3)), ((0,0),(3,3))]
            sage: tuple(G.domain()) == Z3square.vertices()
            True
            sage: G.orbit(Z3square.vertices()[0])
            ((0, 0), (1, 2), (3, 3), (2, 1))

            sage: cell24 = LatticePolytope_PPL(
            ...   (1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(1,-1,-1,1),(0,0,-1,1),
            ...   (0,-1,0,1),(-1,0,0,1),(1,0,0,-1),(0,1,0,-1),(0,0,1,-1),(-1,1,1,-1),
            ...   (1,-1,-1,0),(0,0,-1,0),(0,-1,0,0),(-1,0,0,0),(1,-1,0,0),(1,0,-1,0),
            ...   (0,1,1,-1),(-1,1,1,0),(-1,1,0,0),(-1,0,1,0),(0,-1,-1,1),(0,0,0,-1))
            sage: cell24.restricted_automorphism_group().cardinality()
            1152
        """
        if not self.is_full_dimensional():
            return self.affine_lattice_polytope().\
                restricted_automorphism_group(vertex_labels=vertex_labels)
        if vertex_labels is None:
            vertex_labels = self.vertices()
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.graphs.graph import Graph
        # good coordinates for the vertices
        v_list = []
        for v in self.minimized_generators():
            assert v.divisor().is_one()
            v_coords = (1,) + v.coefficients()
            v_list.append(vector(v_coords))

        # Finally, construct the graph
        Qinv = sum( v.column() * v.row() for v in v_list ).inverse()
        G = Graph()
        for i in range(0,len(v_list)):
            for j in range(i+1,len(v_list)):
                v_i = v_list[i]
                v_j = v_list[j]
                G.add_edge(vertex_labels[i], vertex_labels[j], v_i * Qinv * v_j)
        return G.automorphism_group(edge_labels=True)

    @cached_method
    def lattice_automorphism_group(self, points=None, point_labels=None):
        """
        The integral subgroup of the restricted automorphism group.

        INPUT:

        - ``points`` -- A tuple of coordinate vectors or ``None``
          (default). If specified, the points must form complete
          orbits under the lattice automorphism group. If ``None`` all
          vertices are used.

        - ``point_labels`` -- A tuple of labels for the ``points`` or
          ``None`` (default). These will be used as labels for the do
          permutation group. If ``None`` the ``points`` will be used
          themselves.

        OUTPUT:

        The integral subgroup of the restricted automorphism group
        acting on the given ``points``, or all vertices if not
        specified.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: Z3square = LatticePolytope_PPL((0,0), (1,2), (2,1), (3,3))
            sage: Z3square.lattice_automorphism_group()
            Permutation Group with generators [(), ((1,2),(2,1)),
            ((0,0),(3,3)), ((0,0),(3,3))((1,2),(2,1))]

            sage: G1 = Z3square.lattice_automorphism_group(point_labels=(1,2,3,4));  G1
            Permutation Group with generators [(), (2,3), (1,4), (1,4)(2,3)]
            sage: G1.cardinality()
            4

            sage: G2 = Z3square.restricted_automorphism_group(vertex_labels=(1,2,3,4)); G2
            Permutation Group with generators [(2,3), (1,2)(3,4), (1,4)]
            sage: G2.cardinality()
            8

            sage: points = Z3square.integral_points();  points
            ((0, 0), (1, 1), (1, 2), (2, 1), (2, 2), (3, 3))
            sage: Z3square.lattice_automorphism_group(points, point_labels=(1,2,3,4,5,6))
            Permutation Group with generators [(), (3,4), (1,6)(2,5), (1,6)(2,5)(3,4)]
        """
        if not self.is_full_dimensional():
            return self.affine_lattice_polytope().lattice_automorphism_group()

        if points is None:
            points = self.vertices()
        if point_labels is None:
            point_labels = tuple(points)
        points = [ vector(ZZ, [1]+v.list()) for v in points ]
        map(lambda x:x.set_immutable(), points)

        vertices = [ vector(ZZ, [1]+v.list()) for v in self.vertices() ]
        pivots = matrix(ZZ, vertices).pivot_rows()
        basis = matrix(ZZ, [ vertices[i] for i in pivots ])
        Mat_ZZ = basis.parent()
        basis_inverse = basis.inverse()

        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
        lattice_gens = []
        G = self.restricted_automorphism_group(
            vertex_labels=tuple(range(len(vertices))))
        for g in G:
            image = matrix(ZZ, [ vertices[g(i)] for i in pivots ])
            m = basis_inverse*image
            if m not in Mat_ZZ:
                continue
            perm_list = [ point_labels[points.index(p*m)]
                          for p in points ]
            lattice_gens.append(perm_list)
        return PermutationGroup(lattice_gens, domain=point_labels)

    def affine_normal_form(self, output='matrix'):
        r"""
        Return the affine normal form of ``self``.

        Note that this method calls
        :meth:`~sage.geometry.lattice_polytope.LatticePolytopeClass.affine_normal_form`
        with parameter ``ignore_embedding=True`` in order to perform the
        computation, which is very slow.

        INPUT:

        - ``output`` -- string, (default: ``'matrix'``) 
            One of ``'matrix'``, ``'hash'``, or ``'tuple'``,
            determining how to represent the affine normal form
            of ``self``.
            If ``output='matrix'``, return a matrix whose columns correspond
            to the vertices of ``self`` in affine normal form.
            If ``output='tuple'``, a tuple of the same points is returned
            instead.
            If ``output='hash'``, the md5 hash of the string representing
            the affine normal form matrix is returned.

        OUTPUT:
 
        A matrix, a string, or a tuple. See ``output`` for more details.

        EXAMPLES::
            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((1, 2), (-1, 3), (0, 0))
            sage: p.affine_normal_form()
            [0 1 2]
            [0 0 5]
            sage: p.affine_normal_form(output='tuple')
            ((0, 0), (1, 0), (2, 5))
            sage: q = LatticePolytope_PPL((1, 0, 0), (0, 0, 1), (-1, 0, -1))
            sage: q.affine_normal_form()
            [0 1 2]
            [0 0 3]
            sage: r = LatticePolytope_PPL((1, 0, 0), (0, 1, 0),\
            ....: (0, 0, 1), (-1, -1, -1))
            sage: r.affine_normal_form()
            [0 1 0 3]
            [0 0 1 3]
            [0 0 0 4]
            sage: r.affine_normal_form(output='hash')
            '19768c02d2ebecaa737074dff0b91a5e'
        """
        # This is a terrible way of doing things
        # but for now it's the simplest thing to do
        if not hasattr(self, "_anf"):
            from sage.geometry.lattice_polytope import LatticePolytope
            lp = LatticePolytope(self.vertices())
            self._anf = lp.affine_normal_form(ignore_embedding = True)
        if output == 'matrix':
            return self._anf
        elif output == 'hash':
            import hashlib
            h = hashlib.md5(self._anf.str()).hexdigest()
            return h
        elif output == 'tuple':
            return tuple(self._anf.columns())
        else:
            raise ValueError(output + ' is not a valid parameter for output.')
    
    def normal_form(self, output='matrix'):
        r"""
        Return the normal form of ``self``.

        Note that this method calls
        :meth:`~sage.geometry.lattice_polytope.LatticePolytopeClass.normal_form`
        with parameter ``ignore_embedding=True`` in order to perform the
        computation, which is very slow.

        INPUT:

        - ``output`` -- string, (default: ``'matrix'``) 
            One of ``'matrix'``, ``'hash'``, or ``'tuple'``,
            determining how to represent the affine normal form
            of ``self``.
            If ``output='matrix'``, return a matrix whose columns correspond
            to the vertices of ``self`` in affine normal form.
            If ``output='tuple'``, a tuple of the same points is returned
            instead.
            If ``output='hash'``, the md5 hash of the string representing
            the affine normal form matrix is returned.

        OUTPUT:
 
        A matrix, a string, or a tuple. See ``output`` for more details.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((1, 2), (-1, 3), (0, 0))
            sage: p.normal_form()
            [0 1 4]
            [0 0 5]
            sage: q = LatticePolytope_PPL((1, 0, 0), (0, 0, 1), (-1, 0, -1))
            sage: q.normal_form()
            [ 1  0 -1]
            [ 0  1 -1]
        """
        # This is a terrible way of doing things
        # but for now it's the simplest thing to do
        if not hasattr(self, "_nf"):
            from sage.geometry.lattice_polytope import LatticePolytope
            lp = LatticePolytope(self.vertices())
            self._nf = lp.normal_form(ignore_embedding = True)
        if output == 'matrix':
            return self._nf
        elif output == 'hash':
            import hashlib
            h = hashlib.md5(self._nf.str()).hexdigest()
            return h
        elif output == 'tuple':
            return tuple(self._nf.columns())
        else:
            raise ValueError(output + ' is not a valid parameter for output.')

    def is_embeddable_in(self, other, affine=True):
        r"""
        Check whether ``self`` can be embedded in ``other``.

        INPUT:

        - ``other`` -- :class:`LatticePolytope_PPL_class`, the polytope
            ``self`` is to be embedded into.

        - ``affine`` -- (default: ``True``) If ``True``, check whether
            there exists an affine transformation embedding ``self``
            in ``other``. If ``False``, allow only linear transformations.

        OUTPUT:

        A Boolean.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: ReflexivePolytope_PPL
            sage: p, q, r = [ReflexivePolytope_PPL(2, i) for i in [0, 1, 15]]
            sage: p.is_embeddable_in(q)
            False
            sage: q.is_embeddable_in(p)
            False
            sage: p.is_embeddable_in(r)
            True
        """
        try:
            self.embed_in(other, output='sub_polytope')
            return True
        except LatticePolytopeNoEmbeddingError:
            return False

    def is_isomorphic(self, other, affine=True):
        r"""
        Test if ``self`` and ``polytope`` are isomorphic.

        INPUT:

        - ``other`` -- :class:`LatticePolytope_PPL_class`, the polytope
            ``self`` is compared to.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((1, 0, 0), (0, 0, 1), (-1, 0, -1))
            sage: q = LatticePolytope_PPL((1, 0), (0, 1), (-1, -1))
            sage: r = LatticePolytope_PPL((1, 0), (0, 1), (-2, -2))
            sage: p.is_isomorphic(q)
            True
            sage: p.is_isomorphic(r)
            False
            sage: s = q.embed_in(r, output='sub_polytope')
            sage: s
            A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
            sage: p.is_isomorphic(s)
            True            
        """
        if other.vertices() == self.vertices():
            return True
        if len(other.vertices()) != len(self.vertices()):
            return False
        if other.affine_dimension() != self.affine_dimension():
            return False
        if other.n_integral_points() != self.n_integral_points():
            return False
        # This is a terribly lazy way of doing things:
        if affine:
            return other.affine_normal_form() == self.affine_normal_form()
        else:
            return other.normal_form() == self.normal_form()

    def sub_polytope_generator(self):
        """
        Generate the maximal lattice sub-polytopes.

        OUTPUT:

        A generator yielding the maximal (with respect to inclusion)
        lattice sub polytopes. That is, each can be gotten as the
        convex hull of the integral points of ``self`` with one vertex
        removed.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: P = LatticePolytope_PPL((1,0,0), (0,1,0), (0,0,1), (-1,-1,-1))
            sage: for p in P.sub_polytope_generator():
            ....:     print p.vertices()
            ((0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0))
            ((-1, -1, -1), (0, 0, 0), (0, 1, 0), (1, 0, 0))
            ((-1, -1, -1), (0, 0, 0), (0, 0, 1), (1, 0, 0))
            ((-1, -1, -1), (0, 0, 0), (0, 0, 1), (0, 1, 0))
        """
        pointset = set(self.integral_points())
        for v in self.vertices():
            sub = list(pointset.difference([v]))
            yield LatticePolytope_PPL(*sub)

    def sub_polytopes(self, affine=True, return_hashes=False, max_batch=5000,
                      verbose=False, vertices_only=False):
        r"""
        Compute all subpolytopes of ``self``.

        INPUT:

        - ``affine`` -- (default: ``True``) If ``True``, identify sub polytopes
            that are related by affine transformation. If ``False``,
            then only those related by linear transformations are
            identified.

        - ``return_hashes`` -- (default: ``False``) If ``True``, then
            the return value is a tuple of tuples of dicts structured as
            follows::
            ``result[dim][integral_points - 1][md5_hash] = sub_polytope``
            If ``False``, just return a tuple of all subpolytopes.

        - ``max_batch`` -- (default: ``5000``) Integer, set the
            batch size submitted to PALP. Note that ``max_batch`` is
            not the number of polytopes submitted to PALP, but
            the number of polytopes whose subpolytopes are sent to PALP
            in one pass.
            Note that setting ``max_batch`` to large values may result
            in a single very large file being processed by PALP.

        - ``vertices_only`` -- (default: ``False``) If ``False``,
            return a matrix of vertices, if ``True`` return the corresponding
            ``LatticePolytope_PPL`` instead.

        - ``verbose`` -- (default: ``False``) If ``True``, print
            debugging information.

        OUTPUT:

        A tuple. See parameter ``return_hashes`` for more information.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: LatticePolytope_PPL, LatticePolytope_PPL_class
            sage: p = LatticePolytope_PPL((1, 0, 0), (0, 1, 0),\
            ....: (0, 0, 1), (-1, -1, -1))
            sage: p.sub_polytopes()
            (A 0-dimensional lattice polytope in ZZ^3 with 1 vertex,
             A 1-dimensional lattice polytope in ZZ^3 with 2 vertices,
             A 2-dimensional lattice polytope in ZZ^3 with 3 vertices,
             A 3-dimensional lattice polytope in ZZ^3 with 4 vertices,
             A 3-dimensional lattice polytope in ZZ^3 with 4 vertices)
            sage: p.sub_polytopes(affine=False)
            (A 0-dimensional lattice polytope in ZZ^3 with 1 vertex,
             A 0-dimensional lattice polytope in ZZ^3 with 1 vertex,
             A 1-dimensional lattice polytope in ZZ^3 with 2 vertices,
             A 1-dimensional lattice polytope in ZZ^3 with 2 vertices,
             A 2-dimensional lattice polytope in ZZ^3 with 3 vertices,
             A 2-dimensional lattice polytope in ZZ^3 with 3 vertices,
             A 3-dimensional lattice polytope in ZZ^3 with 4 vertices,
             A 3-dimensional lattice polytope in ZZ^3 with 4 vertices)
       
        Apart from the ordering, this method gives the same output
        as the method of LatticePolygon_PPL_class. However, it is
        significantly faster for large polytopes::

            sage: p = LatticePolytope_PPL((1,0),(0,1),(-1,-4),(2,5),(3,9))
            sage: s1 = p.sub_polytopes() # long time -- 5.1s on i7-2600
            sage: s2 = LatticePolytope_PPL_class.sub_polytopes(p) # long time -- 1s on i7-2600
            sage: len(s1) == len(s2) # long time
            True  
        """
        from sage.matrix.matrix import is_Matrix
        def process_poly(poly, force_process = False):
            if is_Matrix(poly) and not vertices_only:
                return LatticePolytope_PPL(*poly.columns())
            elif vertices_only or force_process:
                m = matrix(poly.vertices()).transpose()
                m.set_immutable()
                return m
            else:
                return poly

        if (affine and not hasattr(self, "_sub_polytopes_affine")) or \
           (not affine and not hasattr(self, "_sub_polytopes_normal")) :
            from sage.functions.other import ceil

            d = self.affine_dimension()
            n_pts = self.n_integral_points()
            # Sort hashes by dimension and by number of integral points
            hashes = tuple([[dict() for i in range(n_pts)] 
                            for j in range(d + 1)])
            # Add self
            if affine:
                h = self.affine_normal_form(output='hash')
            else:
                h = self.normal_form(output='hash')
            hashes[d][n_pts - 1][h] = process_poly(self)
            todo = [self]
            while n_pts > 1:
                n_pts -= 1
                n_set = len(todo)
                if verbose:
                    print 'The set with', n_pts, 'integral points has', \
                           n_set, 'polytopes.'
                # Ha! We can simply compute the subpolytopes of all polytopes in todo
                # But let's slice this up into smaller batches
                n_batches = ceil(n_set / float(max_batch))
                if verbose:
                    print 'The set is split up into', n_batches, 'batches.'
                batches = [todo[max_batch*k:max_batch*(k + 1)] for k in range(n_batches)]
                todo = []
                affine_hashes = [dict() for j in range(d + 1)]
                for k, batch in enumerate(batches):
                    if verbose:
                        print "  Batch", k + 1, ":"
                    sub_polys = [p for q in batch for p in q.sub_polytope_generator() \
                                 if not p.is_empty()]
                    if verbose:
                        print '    Assembled list of length', len(sub_polys), '.'
                    _compute_normal_forms(sub_polys)
                    if verbose:
                        print '    Computed normal forms.'
                    p_list = tuple([list() for j in range(d + 1)])
                    for p in sub_polys:
                        h = p.normal_form(output='hash')
                        pd = p.affine_dimension()
                        relevant_hashes = hashes[pd][n_pts - 1]
                        if relevant_hashes.has_key(h):
                            continue
                        if affine:
                            p_list[pd].append(p)
                            relevant_hashes[h] = None
                        else:
                            relevant_hashes[h] = process_poly(p)
                            todo.append(p)
                    if affine:
                        _compute_affine_normal_forms([i for bucket in p_list for i in bucket])
                        # Now check them again
                        # This may seem inefficient, but it actually
                        # leads to fewer calculations
                        for bd, bucket in enumerate(p_list):
                            for p in bucket:
                                h = p.affine_normal_form(output='hash')
                                if affine_hashes[bd].has_key(h):
                                    continue
                                affine_hashes[bd][h] = process_poly(p)
                                todo.append(p)
                    if verbose:
                        print '    Compared and hashed them.'
                # Overwrite the normal hashes
                if affine:
                    before = sum([len(hashes[bd][n_pts - 1]) for bd in range(d + 1)])
                    for bd in range(d + 1):
                        hashes[bd][n_pts - 1] = affine_hashes[bd]
                    after = sum([len(hashes[bd][n_pts - 1]) for bd in range(d + 1)])
                    if verbose: 
                        print 'Reduced', before, 'polytopes to', after, \
                              'polytopes using affine transformations.'
                if verbose:
                    print
            hashes = tuple([tuple(i) for i in hashes])
            if affine:
                self._sub_polytopes_affine = hashes
            else:
                self._sub_polytopes_normal = hashes
        if affine:
            sub_polys = self._sub_polytopes_affine
        else:
            sub_polys = self._sub_polytopes_normal
        # Convert between matrix and LatticePolytope_PPL
        poly_iterator = (k for i in sub_polys for j in i for k in j.values())
        first_poly = poly_iterator.next()
        if is_Matrix(first_poly) <> vertices_only:
            for i in sub_polys:
                for j in i:
                    for k in j:
                        j[k] = process_poly(j[k])
        # Choose output form
        if not return_hashes:
            return tuple([k for i in sub_polys for j in i for k in j.values()])
        else:
            return sub_polys

    def find_isomorphism(self, other):
        r"""
        Find isomorphism from ``self`` to ``other``.

        Note that this currently does not anything but call
        :meth:`~sage.geometry.lattice_polytope.LatticePolytopeClass.find_isomorphism`
        in order to perform the computations.

        INPUT:

        - ``other`` -- :class:`LatticePolytope_PPL_class`, the polytope
            ``self`` is to be mapped to.

        OUTPUT:

        The :class:`~sage.geometry.polyhedron.lattice_euclidean_group_element.LatticeEuclideanGroupElement`
        corresponding to the isomorphism from ``self`` to ``other``.
        If no isomorphism exists, a :class:`~sage.geometry.polyhedron.lattice_euclidean_group_element.LatticePolytopesNotIsomorphicError`
        is raised.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((1, 0, 0), (0, 0, 1), (-1, 0, -1))
            sage: q = LatticePolytope_PPL((1, 0), (0, 1), (-1, -1))
            sage: r = LatticePolytope_PPL((1, 0), (0, 1), (-2, -2))
            sage: iso = p.find_isomorphism(q)
            sage: iso
            The map A*x+b with A=
            [1 0 0]
            [0 0 1]
            b = 
            (0, 0)
            sage: iso(p) == q
            True
            sage: p.find_isomorphism(r)
            Traceback (most recent call last):
            ...
            LatticePolytopesNotIsomorphicError: The polytopes are not isomorphic.
        """
        # Very lazy way of doing things
        from sage.geometry.lattice_polytope import LatticePolytope
        lp1 = LatticePolytope(self.vertices())
        lp2 = LatticePolytope(other.vertices())
        iso = lp1.find_isomorphism(lp2)
        if not iso:
            raise LatticePolytopesNotIsomorphicError(
                'The polytopes are not isomorphic.')
        return iso

    @cached_method
    def _find_isomorphism_to_subreflexive_polygon(self):
        """
        Find an isomorphism to a sub-polytope of a maximal reflexive
        polygon.

        OUTPUT:

        A tuple consisting of the ambient reflexive polytope, the
        subpolytope, and the embedding of ``self`` into the ambient
        polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: polygon = LatticePolytope_PPL((0,0,2,1),(0,1,2,0),(2,3,0,0),(2,0,0,3))
            sage: polygon._find_isomorphism_to_subreflexive_polygon()
            (A 2-dimensional lattice polytope in ZZ^2 with 3 vertices,
             A 2-dimensional lattice polytope in ZZ^2 with 4 vertices,
             The map A*x+b with A=
             [ 1  1]
             [ 0  1]
             [-1 -1]
             [ 1  0]
             b =
             (-1, 0, 3, 0))
            sage: ambient, sub, embedding = _
            sage: ambient.vertices()
            ((0, 0), (0, 3), (3, 0))
            sage: sub.vertices()
            ((0, 1), (3, 0), (0, 3), (1, 0))
        """
        from ppl_lattice_polygon import sub_reflexive_polygons
        for p, ambient in sub_reflexive_polygons():
            try:
                return (ambient, p, p.find_isomorphism(self))
            except LatticePolytopesNotIsomorphicError:
                pass
        raise LatticePolytopeNoEmbeddingError('not a sub-polytope of a reflexive polygon')

    @cached_method
    def _find_isomorphism_to_subreflexive_polytope(self):
        """
        Find an isomorphism to a sub-polytope of a maximal reflexive
        polytope.

        OUTPUT:

        A tuple consisting of the ambient reflexive polytope, the
        subpolytope, and the embedding of ``self`` into the ambient
        polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((1, 0, 0), (0, 1, 0), (0, 0, 1),\
            ....: (-1, -1, -1))
            sage: iso = p._find_isomorphism_to_subreflexive_polytope() # long time
        """
        d = self.affine_dimension()
        npts = self.n_integral_points()
        # Check that the values are allowed.
        sub_polys = _reflexive_sub_polytopes(d, npts)
        h = self.affine_normal_form(output='hash')
        try:
            m, id = sub_polys[h]
            vertices = m.columns()
            sub_poly = LatticePolytope_PPL(*vertices)
            r = ReflexivePolytope_PPL(3, id)
            return (r, sub_poly, sub_poly.find_isomorphism(self))
        except KeyError:
            raise LatticePolytopeNoEmbeddingError('Not a sub-polytope ' +
                'of a 3-d reflexive polytope.')

    def embed_in(self, other, output='hom', inverse=False,
                   vertices_only=False):
        r"""
        Find an embedding as a sub-polytope of ``other``.

        INPUT:

        - ``other`` -- :class:`LatticePolytope_PPL_class`,
            the polytope in which ``self`` is to be embedded.

        - ``output`` -- string. One of ``'hom'`` (default),
            ``inverse_hom``, ``points``, or ``sub_polytope``.
            How the embedding is returned. See the output
            section for details.

        - ``inverse`` -- (default: ``False``). If ``True``,
            return the map from ``self`` to ``other`` instead
            of the other way round.

        OUTPUT:

        An embedding into ``other``. Depending on the
        ``output`` option slightly different data is returned.

        - If ``output='hom'``, a map from ``other`` onto
          ``self`` is returned.

        - If ``output='inverse_hom'``, a map from ``self`` onto
          ``other`` is returned.

        - If ``output='points'``, a dictionary containing the integral
          points of ``self`` as keys and the corresponding integral
          point of ``other`` as value.

        - If ``output='sub_polytope'``, the subpolytope of ``other``
          isomorphic to ``self`` is returned.

        If there is no such embedding, a
        :class:`~sage.geometry.polyhedron.lattice_euclidean_group_element.LatticePolytopeNoEmbeddingError`
        is raised.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
            ....: LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((1, 0), (0, 1), (-1, -1))
            sage: q = LatticePolytope_PPL((1, 2), (-1, 1), (-2, -1), (1, -1))
            sage: hom = p.embed_in(q)
            sage: sub_poly = p.embed_in(q, output='sub_polytope')
            sage: hom(sub_poly) == p
            True
        """
        d = self.affine_dimension()
        np = self.n_integral_points()
        if d > other.affine_dimension() or np > other.n_integral_points():
            raise LatticePolytopeNoEmbeddingError(
                'The polytope cannot be bigger than its embedding.')
        hash = self.affine_normal_form(output='hash')
        # Get the vertices only in order to speed up things for cases
        # with pre-computed subpolytopes.
        sub_polys = other.sub_polytopes(return_hashes=True, vertices_only=True)
        try:
            vertices = sub_polys[d][np - 1][hash]
            sub_poly = LatticePolytope_PPL(*vertices.columns())
        except KeyError:
            raise LatticePolytopeNoEmbeddingError('Cannot embed the polytope.')
        if output == 'sub_polytope':
            return sub_poly
        elif output == 'inverse_hom':
            iso = self.find_isomorphism(sub_poly)
            return iso
        else:
            iso = sub_poly.find_isomorphism(self)
            if output == 'hom':
                return iso
            elif output == 'points':
                points = dict()
                for p in sub_poly.integral_points():
                    points[ tuple(iso(p)) ] = p
                return points
            else:
                raise ValueError(output + ' is not a valid output parameter.')

    def embed_in_reflexive_polytope(self, output='hom'):
        r"""
        Find an embedding as a sub-polytope of a maximal reflexive
        polytope.

        INPUT:

        - ``output`` -- string. One of ``'hom'`` (default),
          ``'polytope'``, or ``points``. How the embedding is
          returned. See the output section for details.

        OUTPUT:

        An embedding into a reflexive polytope. Depending on the
        ``output`` option slightly different data is returned.

        - If ``output='hom'``, a map from a reflexive polytope onto
          ``self`` is returned.

        - If ``output='polytope'``, a reflexive polytope that contains
          ``self`` (up to a lattice linear transformation) is
          returned. That is, the domain of the ``output='hom'`` map is
          returned. If the affine span of ``self`` is less or equal
          2-dimensional, the output is one of the following three
          possibilities::
          :func:`~sage.geometry.polyhedron.ppl_lattice_polygon.polar_P2_polytope`,
          :func:`~sage.geometry.polyhedron.ppl_lattice_polygon.polar_P1xP1_polytope`,
          or
          :func:`~sage.geometry.polyhedron.ppl_lattice_polygon.polar_P2_112_polytope`.

        - If ``output='points'``, a dictionary containing the integral
          points of ``self`` as keys and the corresponding integral
          point of the reflexive polytope as value.

        If there is no such embedding, a
        :class:`~sage.geometry.polyhedron.lattice_euclidean_group_element.LatticePolytopeNoEmbeddingError`
        is raised. Even if it exists, the ambient reflexive polytope
        is usually not uniquely determined an a random but fixed
        choice will be returned.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: polygon = LatticePolytope_PPL((0,0,2,1),(0,1,2,0),(2,3,0,0),(2,0,0,3))
            sage: polygon.embed_in_reflexive_polytope()
            The map A*x+b with A=
            [ 1  1]
            [ 0  1]
            [-1 -1]
            [ 1  0]
            b =
            (-1, 0, 3, 0)
            sage: polygon.embed_in_reflexive_polytope('polytope')
            A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
            sage: polygon.embed_in_reflexive_polytope('points')
            {(0, 0, 2, 1): (1, 0), (2, 1, 0, 2): (2, 1),
             (0, 1, 2, 0): (0, 1), (1, 1, 1, 1): (1, 1),
             (1, 2, 1, 0): (0, 2), (2, 2, 0, 1): (1, 2),
             (2, 3, 0, 0): (0, 3), (1, 0, 1, 2): (2, 0),
             (2, 0, 0, 3): (3, 0)}

            sage: LatticePolytope_PPL((0,0), (4,0), (0,4)).embed_in_reflexive_polytope() # long time
            The map A*x+b with A=
            [0 1 0]
            [0 0 1]
            b = 
            (0, 0)
        """
        # First try embedding in a reflexive polygon
        try:
            ambient, subreflexive, hom = self._find_isomorphism_to_subreflexive_polygon()
        except LatticePolytopeNoEmbeddingError:
            ambient, subreflexive, hom = self._find_isomorphism_to_subreflexive_polytope()
        if output == 'hom':
            return hom
        elif output == 'polytope':
            return ambient
        elif output == 'points':
            points = dict()
            for p in subreflexive.integral_points():
                points[ tuple(hom(p)) ] = p
            return points
        else:
            raise ValueError('output='+str(output)+' is not valid.')


##############################################################################
# Methods to deal with large lists of PPLs at the same time                  #
##############################################################################
def _compute_normal_forms(polytopes):
    r"""
    Compute the normal forms for all ``polytopes``.

    Note that this is considerably faster than calling the method 
    ``normal_form`` of the individual polytopes.

    INPUT:

        - ``polytopes`` -- an iterable of :class:`LatticePolytope_PPL_class`.

    OUTPUT:

    Nothing, the function only fills the private attributes of the
    passed polytopes and raises a ``ValueError`` if something goes wrong.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
        ....: _compute_normal_forms, LatticePolytope_PPL
        sage: lp = LatticePolytope_PPL((0, 0), (1, 1), (-1, 1))
        sage: lp2 = LatticePolytope_PPL((0, 0, 0), (1, 0, 1), (-1, 0, 1))
        sage: _compute_normal_forms([lp, lp2])
        sage: lp._nf
        [0 1 1]
        [0 0 2]
        sage: lp2._nf == lp._nf
        True
    """
    from sage.geometry.lattice_polytope import _palp, read_palp_matrix

    def _compute_full_dim(polys):
        if not polys:
            return
        vertices = [matrix(p.vertices()).transpose() for p in polys]
        path = _palp('poly.x -fN', vertices=vertices)
        f = open(path)
        for p in polys:
            m = read_palp_matrix(f)
            m.set_immutable()
            if not m:
                raise ValueError("Empty normal form returned.")
            p._nf = m
        f.close()

    def _compute_codim(polys):
        add_origins = []
        vertices = []
        # Bring the vertices into desired form
        relevant_polys = []
        for p in polys:
            v = matrix(QQ, p.vertices())
            # A workaround for polytopes with vertices at origin
            if v.is_zero():
                m = matrix([[0]])
                m.set_immutable()
                p._nf = m
                continue
            embedding = v.image().basis_matrix()
            codim_embedding = embedding.ncols() - embedding.nrows()
            codim = p.space_dimension() - p.affine_dimension()
            add_origin = codim_embedding <> codim
            preimages = [embedding.solve_left(i) for i in v.rows()]
            if add_origin:
                preimages += [vector(embedding.nrows()*[0])]
            vertices.append(matrix(preimages).transpose())
            add_origins.append(add_origin)
            relevant_polys.append(p)
        # Process them
        if not relevant_polys:
            return
        path = _palp('poly.x -fN', vertices=vertices)
        del vertices
        f = open(path)
        # The afterwards treatment
        for p, add_origin in zip(polys, add_origins):
            m = read_palp_matrix(f)
            if not m:
                raise ValueError("Empty normal form returned.")
            cols = [i for i in m.columns() if not add_origin or not i.is_zero()]
            m = matrix(cols).transpose()
            m.set_immutable()
            p._nf = m
        f.close()

    polys_full_dim = []
    polys_codim = []
    for p in polytopes:
        if hasattr(p, "_nf"):
            continue
        if p.affine_dimension() <> p.space_dimension():
            polys_codim.append(p)
        else:
            polys_full_dim.append(p)

    _compute_full_dim(polys_full_dim)
    _compute_codim(polys_codim)

def _compute_affine_normal_forms(polytopes):
    r"""
    Compute the affine normal forms for all ``polytopes``.

    Note that this is considerably faster than calling the method 
    ``affine_normal_form`` of the individual polytopes.

    INPUT:

        - ``polytopes`` -- an iterable of :class:`LatticePolytope_PPL_class`.
    OUTPUT:

    Nothing, the function only fills the private attributes of the
    passed polytopes and raises a ``ValueError`` if something goes wrong.
    """
    from sage.geometry.lattice_polytope import _palp, read_palp_matrix
    from copy import copy

    def translated_vertices(p):
        v = matrix(QQ, p.vertices()).transpose()
        r = []
        for i in range(v.ncols()):
            tmp = copy(v)
            for j in range(v.ncols()):
                if i == j:
                    continue
                tmp.add_multiple_of_column(j, i, -1)
            tmp.add_multiple_of_column(i, i, -1)
            r.append(tmp)
        return r

    def _compute_full_dim(polys):
        if not polys:
            return
        # Prepare the vertices
        n_vertices = []
        vertices = []
        for p in polys:
            tmp = translated_vertices(p)
            n_vertices.append(len(tmp))
            vertices += tmp
        path = _palp('poly.x -fN', vertices=vertices)
        del vertices
        f = open(path)
        for p, n in zip(polys, n_vertices):
            tmp = []
            for i in range(n):
                m = read_palp_matrix(f)
                if not m:
                    raise ValueError("Empty normal form returned.")
                m.set_immutable()
                tmp.append(m)
            # Ours is the minimum            
            p._anf = min(tmp)
        f.close()

    def _compute_codim(polys):
        n_vertices = []
        vertices = []
        relevant_polys = []
        # Filter out one vertex polytopes and translate vertices
        for p in polys:
            if p.n_vertices() == 1:
                m = matrix([[0]])
                m.set_immutable()
                p._anf = m
                continue
            tmp = translated_vertices(p)
            n_vertices.append(len(tmp))
            vertices += tmp
            relevant_polys.append(p)

        embedded_vertices = []
        for v in vertices:
            v = v.transpose()
            embedding = v.image().basis_matrix()
            codim_embedding = embedding.ncols() - embedding.nrows()
            preimages = [embedding.solve_left(i) for i in v.rows()]
            embedded_vertices.append(matrix(preimages).transpose())
        vertices = embedded_vertices

        # Process them
        if not relevant_polys:
            return
        path = _palp('poly.x -fN', vertices=vertices)
        del vertices
        f = open(path)
        # The afterwards treatment
        for p, n in zip(relevant_polys, n_vertices):
            tmp = []
            for i in range(n):
                m = read_palp_matrix(f)
                if not m:
                    raise ValueError("Empty normal form returned.")
                m.set_immutable()
                tmp.append(m)
            # Ours is the maximum        
            p._anf = min(tmp)
        f.close()

    polys_full_dim = []
    polys_codim = []
    for p in polytopes:
        if hasattr(p, "_anf"):
            continue
        if p.affine_dimension() <> p.space_dimension():
            polys_codim.append(p)
        else:
            polys_full_dim.append(p)

    _compute_full_dim(polys_full_dim)
    _compute_codim(polys_codim)


########################################################################
#
#  Reflexive lattice polytopes and their subpolytopes
#
########################################################################

# Change the hard-coding at some point!
_subpolytope_path = '/home/pcl337b/jkeitel/Documents/Papers/BGK3/sub_polytopes/'

def maximal_polytopes_3d_indices():
    r"""
    Return the internal indices of the maximal three-dimensional reflexive
    polytopes.

    OUTPUT:
    
    A tuple containing the indices of the maximal three-dimensional reflexive
    polytopes in a somewhat random but fixed order. Note that the order
    of the indices determines which polytope is checked first when
    trying to embed polytopes in a reflexive polytope.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
        ....: maximal_polytopes_3d_indices
        sage: maximal_polytopes_3d_indices()
        (4311,
         4318,
         4317,
         4316,
         4314,
         4312,
         4309,
         4308,
         4302,
         4298,
         4286,
         4283,
         4281,
         4250,
         4237,
         3313)
    """
    return (4311, 4318, 4317, 4316, 4314, 4312,\
            4309, 4308, 4302, 4298, 4286, 4283,\
            4281, 4250, 4237, 3313)

def maximal_polytope(index):
    r"""
    Return the maximal reflexive polytope with index ``index``.

    INPUT:

        - ``index`` -- integer, the index of the polytope.

    OUTPUT:

    A reflexive polytope.
    """
    if index == 4311:
        return polar_P3_polytope()
    else:
        raise NotImplementedError(
            'At the moment we can only deal with polar P^3.')

def ReflexivePolytope_PPL(dim, index, load_sub_polytopes=False):
    r"""
    Return the reflexive polytope of dimension ``dim`` and index
    ``index``.

    INPUT:

    - ``dim`` -- integer, the dimension of the polytope. Allowed values
        are ``2`` and ``3``.

    - ``index`` -- integer, the index of the polytope in the PALP database.

    - ``load_sub_polytopes`` -- (default: ``False``) Whether to load the list
        of sub-polytopes from the hard drive. Only supported for
        3-dimensional polytopes.

    OUTPUT:

    A :class:`LatticePolytope_PPL_class`.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
        ....: ReflexivePolytope_PPL
        sage: p = ReflexivePolytope_PPL(3, 0, load_sub_polytopes=True) # long time
        sage: len(p._sub_polytopes_affine) # long time
        4
        sage: q = ReflexivePolytope_PPL(3, 1) # long time
        sage: hasattr(q, '_sub_polytopes_affine') # long time
        False
    """
    # For the time being read them from there
    from sage.geometry.lattice_polytope import ReflexivePolytope
    from sage.misc.all import load
    if dim == 2:
        vs = ReflexivePolytope(dim, index).vertices().columns()
        return LatticePolytope_PPL(*vs)
    elif dim == 3:
        vs = ReflexivePolytope(dim, index).vertices().columns()
        ppl = LatticePolytope_PPL(*vs)
        if load_sub_polytopes:
            # Load the subpolytopes
            ppl._sub_polytopes_affine = load(
                _subpolytope_path + 'a_' + str(index))
        return ppl
    else:
        raise NotImplementedError('Only reflexive polytopes of dimensions ' + \
                                  '2 and 3 are supported.')

def _reflexive_sub_polytopes(dim, npts):
    r"""
    Return the sub-polytopes of three-dimensional reflexive polytopes
    that are ``dim``-dimensional and have ``npts`` interior points.

    INPUT:

    - ``dim`` - integer, the dimension of the sub polytopes.

    - ``npts`` - integer, the number of integral points.

    OUTPUT:
    
    A dictionary of the structure
    ``{hash:(sub_polytope, index_reflexive_polytope)}``.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import \
        ....: _reflexive_sub_polytopes
        sage: d = _reflexive_sub_polytopes(3, 1)
        sage: len(d)
        1
    """
    from sage.misc.all import load
    if not dim in range(4) or not npts in range(1, 40):
        raise ValueError('The polytopes must have dimension <=3 and '
            + 'integral points <= 39.')
    fname = _subpolytope_path + 's_{0}_{1}'.format(dim, npts)
    try:
        data = load(fname)
    except IOError:
        data = {}
    return data

@cached_function
def polar_P3_polytope():
    r"""
    The polar of the `P^3` polytope

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import polar_P3_polytope
        sage: polar_P3_polytope()
        A 3-dimensional lattice polytope in ZZ^3 with 4 vertices
        sage: _.vertices()
        ((0, 0, 0), (0, 0, 4), (0, 4, 0), (4, 0, 0))
    """
    return LatticePolytope_PPL((0, 0, 0), (4, 0, 0), (0, 4, 0), (0, 0, 4))
