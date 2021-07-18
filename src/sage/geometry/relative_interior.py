r"""
Relative Interiors of Polyhedra and Cones
"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.geometry.convex_set import ConvexSet_relatively_open


class RelativeInterior(ConvexSet_relatively_open):
    r"""
    The relative interior of a polyhedron or cone

    This class should not be used directly. Use methods
    :meth:`~sage.geometry.polyhedron.Polyhedron_base.relative_interior`,
    :meth:`~sage.geometry.polyhedron.Polyhedron_base.interior`,
    :meth:`~sage.geometry.cone.ConvexRationalPolyhedralCone.relative_interior`,
    :meth:`~sage.geometry.cone.ConvexRationalPolyhedralCone.interior` instead.

    EXAMPLES::

        sage: segment = Polyhedron([[1, 2], [3, 4]])
        sage: segment.relative_interior()
        Relative interior of
         a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
        sage: octant = Cone([(1,0,0), (0,1,0), (0,0,1)])
        sage: octant.relative_interior()
        Relative interior of 3-d cone in 3-d lattice N
    """

    def __init__(self, polyhedron):
        r"""
        Initialize ``self``.

        INPUT:

        - ``polyhedron`` - an instance of :class:`Polyhedron_base` or
          :class:`ConvexRationalPolyhedralCone`.

        TESTS::

            sage: P = Polyhedron([[1, 2], [3, 4]])
            sage: from sage.geometry.relative_interior import RelativeInterior
            sage: TestSuite(RelativeInterior(P)).run()
        """
        self._polyhedron = polyhedron

    def __hash__(self):
        r"""
        TESTS::

            sage: P = Polyhedron([[1, 2], [3, 4]])
            sage: Q = Polyhedron([[3, 4], [1, 2]])
            sage: hash(P.relative_interior()) == hash(Q.relative_interior())
            True
        """
        return hash(self._polyhedron) ^ 1789

    def __contains__(self, point):
        r"""
        Return whether ``self`` contains ``point``.

        EXAMPLES::

            sage: octant = Cone([(1,0,0), (0,1,0), (0,0,1)])
            sage: ri_octant = octant.relative_interior(); ri_octant
            Relative interior of 3-d cone in 3-d lattice N
            sage: (1, 1, 1) in ri_octant
            True
            sage: (1, 0, 0) in ri_octant
            False
        """
        return self._polyhedron.relative_interior_contains(point)

    def ambient(self):
        r"""
        Return the ambient convex set or space.

        EXAMPLES::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: ri_segment = segment.relative_interior(); ri_segment
            Relative interior of
             a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: ri_segment.ambient()
            Vector space of dimension 2 over Rational Field
        """
        return self._polyhedron.ambient()

    def ambient_vector_space(self, base_field=None):
        r"""
        Return the ambient vector space.

        EXAMPLES::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: ri_segment = segment.relative_interior(); ri_segment
            Relative interior of
             a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: ri_segment.ambient_vector_space()
            Vector space of dimension 2 over Rational Field
        """
        return self._polyhedron.ambient_vector_space(base_field=base_field)

    def ambient_dim(self):
        r"""
        Return the dimension of the ambient space.

        EXAMPLES::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: segment.ambient_dim()
            2
            sage: ri_segment = segment.relative_interior(); ri_segment
            Relative interior of
             a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: ri_segment.ambient_dim()
            2
        """
        return self._polyhedron.ambient_dim()

    def dim(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: segment.dim()
            1
            sage: ri_segment = segment.relative_interior(); ri_segment
            Relative interior of
             a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: ri_segment.dim()
            1
        """
        return self._polyhedron.dim()

    def interior(self):
        r"""
        Return the interior of ``self``.

        EXAMPLES::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: ri_segment = segment.relative_interior(); ri_segment
            Relative interior of
             a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: ri_segment.interior()
            The empty polyhedron in ZZ^2

            sage: octant = Cone([(1,0,0), (0,1,0), (0,0,1)])
            sage: ri_octant = octant.relative_interior(); ri_octant
            Relative interior of 3-d cone in 3-d lattice N
            sage: ri_octant.interior() is ri_octant
            True
        """
        return self._polyhedron.interior()

    def relative_interior(self):
        r"""
        Return the relative interior of ``self``.

        As ``self`` is already relatively open, this method just returns ``self``.

        EXAMPLES::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: ri_segment = segment.relative_interior(); ri_segment
            Relative interior of
             a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: ri_segment.relative_interior() is ri_segment
            True
        """
        return self

    def closure(self):
        r"""
        Return the topological closure of ``self``.

        EXAMPLES::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: ri_segment = segment.relative_interior(); ri_segment
            Relative interior of
             a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: ri_segment.closure() is segment
            True
        """
        return self._polyhedron

    def is_universe(self):
        r"""
        Return whether ``self`` is the whole ambient space

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: ri_segment = segment.relative_interior(); ri_segment
            Relative interior of
             a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: ri_segment.is_universe()
            False
        """
        # Relies on ``self`` not set up for polyhedra that are already
        # relatively open themselves.
        assert not self._polyhedron.is_universe()
        return False

    def is_closed(self):
        r"""
        Return whether ``self`` is closed.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: ri_segment = segment.relative_interior(); ri_segment
            Relative interior of a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: ri_segment.is_closed()
            False
        """
        # Relies on ``self`` not set up for polyhedra that are already
        # relatively open themselves.
        assert not self._polyhedron.is_relatively_open()
        return False

    def _some_elements_(self):
        r"""
        Generate some points of ``self``.

        If ``self`` is empty, no points are generated; no exception will be raised.

        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: ri_P = P.relative_interior()
            sage: ri_P.an_element()              # indirect doctest
            (1/4, 1/4, 1/4, 1/4)
            sage: ri_P.some_elements()           # indirect doctest
            [(1/4, 1/4, 1/4, 1/4), (1/2, 1/4, 1/8, 1/8)]
        """
        for p in self._polyhedron.some_elements():
            if p in self:
                yield p

    def _repr_(self):
        r"""
        Return a description of ``self``.

        EXAMPLES::

            sage: P = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1]])
            sage: P.relative_interior()._repr_()
            'Relative interior of a 2-dimensional polyhedron in ZZ^4 defined as the convex hull of 3 vertices'
            sage: P.rename('A')
            sage: P.relative_interior()._repr_()
            'Relative interior of A'
        """
        repr_P = repr(self._polyhedron)
        if repr_P.startswith('A '):
            repr_P = 'a ' + repr_P[2:]
        return 'Relative interior of ' + repr_P

    def __eq__(self, other):
        r"""
        Compare ``self`` and ``other``.

        INPUT:

        - ``other`` -- any object

        EXAMPLES::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: ri_segment = segment.relative_interior(); ri_segment
            Relative interior of
             a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: segment2 = Polyhedron([[1, 2], [3, 4]], base_ring=AA)
            sage: ri_segment2 = segment2.relative_interior(); ri_segment2
            Relative interior of
             a 1-dimensional polyhedron in AA^2 defined as the convex hull of 2 vertices
            sage: ri_segment == ri_segment2
            True

        TESTS::

            sage: empty = Polyhedron(ambient_dim=2)
            sage: ri_segment == empty
            False
        """
        if type(self) != type(other):
            return False
        return self._polyhedron == other._polyhedron

    def __ne__(self, other):
        r"""
        Compare ``self`` and ``other``.

        INPUT:

        - ``other`` -- any object

        TESTS::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: ri_segment = segment.relative_interior(); ri_segment
            Relative interior of
             a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: segment2 = Polyhedron([[1, 2], [3, 4]], base_ring=AA)
            sage: ri_segment2 = segment2.relative_interior(); ri_segment2
            Relative interior of
             a 1-dimensional polyhedron in AA^2 defined as the convex hull of 2 vertices
            sage: ri_segment != ri_segment2
            False
        """
        return not (self == other)
