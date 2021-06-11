r"""
Convex Sets
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

from sage.structure.sage_object import SageObject
from sage.misc.abstract_method import abstract_method

class ConvexSet_base(SageObject):
    """
    Abstract base class for convex sets.
    """

    def is_empty(self):
        r"""
        Test whether ``self`` is the empty set.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = LatticePolytope([], lattice=ToricLattice(3).dual()); p
            -1-d lattice polytope in 3-d lattice M
            sage: p.is_empty()
            True
        """
        return self.dim() < 0

    def is_universe(self):
        r"""
        Test whether ``self`` is the whole ambient space.

        OUTPUT:

        Boolean.
        """
        if not self.is_full_dimensional():
            return False
        raise NotImplementedError

    @abstract_method
    def dim(self):
        r"""
        Return the dimension of ``self``.
        """

    def dimension(self):
        r"""
        Return the dimension of ``self``.

        This is the same as :meth:`dim`.

        EXAMPLES::

            sage: from sage.geometry.convex_set import ConvexSet_base
            sage: class ExampleSet(ConvexSet_base):
            ....:     def dim(self):
            ....:         return 42
            sage: ExampleSet().dimension()
            42
        """
        return self.dim()

    @abstract_method
    def ambient_dim(self):
        r"""
        Return the dimension of the ambient space.
        """

    def ambient_dimension(self):
        r"""
        Return the dimension of ``self``.

        This is the same as :meth:`ambient_dim`.

        EXAMPLES::

            sage: from sage.geometry.convex_set import ConvexSet_base
            sage: class ExampleSet(ConvexSet_base):
            ....:     def ambient_dim(self):
            ....:         return 91
            sage: ExampleSet().ambient_dimension()
            91
        """
        return self.ambient_dim()

    def codimension(self):
        r"""
        Return the codimension of ``self``.

        An alias is :meth:`codim`.

        EXAMPLES::

            sage: Polyhedron(vertices=[(1,2,3)], rays=[(1,0,0)]).codimension()
            2
        """
        return self.ambient_dim() - self.dim()

    codim = codimension

    def is_full_dimensional(self):
        r"""
        Return whether ``self`` is full dimensional.

        OUTPUT:

        Boolean. Whether the polyhedron is not contained in any strict
        affine subspace.

        EXAMPLES::

            sage: c = Cone([(1,0)])
            sage: c.is_full_dimensional()
            False

            sage: polytopes.hypercube(3).is_full_dimensional()
            True
            sage: Polyhedron(vertices=[(1,2,3)], rays=[(1,0,0)]).is_full_dimensional()
            False
        """
        return self.dim() == self.ambient_dim()

    def is_open(self):
        r"""
        Return whether ``self`` is open.

        The default implementation of this method only knows that the
        empty set and the ambient space are open.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.geometry.convex_set import ConvexSet_base
            sage: class ExampleSet(ConvexSet_base):
            ....:     def is_empty(self):
            ....:         return False
            ....:     def is_universe(self):
            ....:         return True
            sage: ExampleSet().is_open()
            True
        """
        if self.is_empty() or self.is_universe():
            return True
        raise NotImplementedError

    def is_relatively_open(self):
        r"""
        Return whether ``self`` is relatively open.

        The default implementation of this method only knows that open
        sets are also relatively open, and in addition singletons are
        relatively open.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.geometry.convex_set import ConvexSet_base
            sage: class ExampleSet(ConvexSet_base):
            ....:     def is_open(self):
            ....:         return True
            sage: ExampleSet().is_relatively_open()
            True
        """
        if self.is_open():
            return True
        if self.dim() == 0:
            return True
        raise NotImplementedError

    def is_closed(self):
        r"""
        Return whether ``self`` is closed.

        The default implementation of this method only knows that the
        empty set, a singleton set, and the ambient space are closed.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.geometry.convex_set import ConvexSet_base
            sage: class ExampleSet(ConvexSet_base):
            ....:     def dim(self):
            ....:         return 0
            sage: ExampleSet().is_closed()
            True
        """
        if self.is_empty() or self.dim() == 0 or self.is_universe():
            return True
        raise NotImplementedError

    def is_compact(self):
        r"""
        Return whether ``self`` is compact.

        The default implementation of this method only knows that a
        non-closed set cannot be compact, and that the empty set and
        a singleton set are compact.

        OUTPUT:

        Boolean.

            sage: from sage.geometry.convex_set import ConvexSet_base
            sage: class ExampleSet(ConvexSet_base):
            ....:     def dim(self):
            ....:         return 0
            sage: ExampleSet().is_compact()
            True
        """
        if not self.is_closed():
            return False
        if self.dim() < 1:
            return True
        raise NotImplementedError

    def closure(self):
        r"""
        Return the topological closure of ``self``.

        EXAMPLES::

            sage: from sage.geometry.convex_set import ConvexSet_closed
            sage: C = ConvexSet_closed()
            sage: C.closure() is C
            True
        """
        if self.is_closed():
            return self
        raise NotImplementedError

    def interior(self):
        r"""
        Return the topological interior of ``self``.

        EXAMPLES::

            sage: from sage.geometry.convex_set import ConvexSet_open
            sage: C = ConvexSet_open()
            sage: C.interior() is C
            True
        """
        if self.is_open():
            return self
        raise NotImplementedError

    def relative_interior(self):
        r"""
        Return the relative interior of ``self``.

        EXAMPLES::

            sage: from sage.geometry.convex_set import ConvexSet_relatively_open
            sage: C = ConvexSet_relatively_open()
            sage: C.relative_interior() is C
            True
        """
        if self.is_relatively_open():
            return self
        raise NotImplementedError

    def _test_convex_set(self, tester=None, **options):
        """
        Run some tests on the methods of :class:`ConvexSet_base`.

        TESTS::

            sage: from sage.geometry.convex_set import ConvexSet_open
            sage: class FaultyConvexSet(ConvexSet_open):
            ....:     def is_universe(self):
            ....:         return True
            ....:     def dim(self):
            ....:         return 42
            ....:     def ambient_dim(self):
            ....:         return 91
            sage: TestSuite(FaultyConvexSet()).run(skip='_test_pickling')
            Failure in _test_convex_set:
            ...
            The following tests failed: _test_convex_set

            sage: class BiggerOnTheInside(ConvexSet_open):
            ....:     def dim(self):
            ....:         return 100000
            ....:     def ambient_dim(self):
            ....:         return 3
            sage: TestSuite(BiggerOnTheInside()).run(skip='_test_pickling')
            Failure in _test_convex_set:
            ...
            The following tests failed: _test_convex_set

        """
        if tester is None:
            tester = self._tester(**options)
        dim = self.dim()
        codim = self.codim()
        tester.assertTrue(dim <= self.ambient_dim())
        if dim >= 0:
            tester.assertTrue(dim + codim == self.ambient_dim())
        if self.is_empty():
            tester.assertTrue(dim == -1)
        if self.is_universe():
            tester.assertTrue(self.is_full_dimensional())
        cl_self = self.closure()
        try:
            int_self = self.interior()
        except NotImplementedError:
            int_self = None
        try:
            relint_self = self.relative_interior()
        except NotImplementedError:
            relint_self = None
        if self.is_full_dimensional():
            tester.assertTrue(int_self == relint_self)
        if self.is_relatively_open():
            tester.assertTrue(self == relint_self)
        if self.is_open():
            tester.assertTrue(self == int_self)
        if self.is_closed():
            tester.assertTrue(self == cl_self)
        if self.is_compact():
            tester.assertTrue(self.is_closed())

    # Optional methods

    @abstract_method(optional=True)
    def affine_hull(self):
        r"""
        Return the affine hull of ``self``.
        """

    @abstract_method(optional=True)
    def cartesian_product(self, other):
        """
        Return the Cartesian product.

        INPUT:

        - ``other`` -- another convex set

        OUTPUT:

        The Cartesian product of ``self`` and ``other``.
        """

    @abstract_method(optional=True)
    def contains(self, point):
        """
        Test whether ``self`` contains the given ``point``.

        INPUT:

        - ``point`` -- a point or its coordinates
        """

    @abstract_method(optional=True)
    def intersection(self, other):
        r"""
        Return the intersection of ``self`` and ``other``.

        INPUT:

        - ``other`` -- another convex set

        OUTPUT:

        The intersection.
        """


    ### TO BE CONSIDERED


    @abstract_method(optional=True)
    def ambient_vector_space(self, base_ring=None):
        r"""
        Return the ambient vector space of ``self``.

        """

    @abstract_method(optional=True)
    def an_affine_basis(self):
        r"""
        Return points in ``self`` that form a basis for the affine span.
        """

    def _test_an_affine_basis(self, tester=None, **options):
        """
        Run tests on the method :meth:`.an_affine_basis`

        TESTS::

            sage: c = Cone([(1,0)])
            sage: c._test_an_affine_basis()
        """
        if tester is None:
            tester = self._tester(**options)
        try:
            if self.an_affine_basis == NotImplemented:
                raise NotImplementedError
            b = self.an_affine_basis()
        except NotImplementedError:
            pass
        else:
            m = matrix([1] + list(v) for v in b)
            tester.assertEqual(m.rank(), self.dim() + 1)
            for v in b:
                tester.assertIn(v, self)

    @abstract_method(optional=True)
    def centroid(self):
        """
        Return the centroid (center of mass) of ``self``.

        The mass is taken with respect to the induced Lebesgue measure,
        see :meth:`volume`.
        """

    @abstract_method(optional=True)
    def minkowski_sum(self, other):
        r"""
        Return the Minkowski sum.

        Minkowski addition of two subsets of a vector space is defined
        as

        .. MATH::

            X \oplus Y =
            \cup_{y\in Y} (X+y) =
            \cup_{x\in X, y\in Y} (x+y)

        See :meth:`minkowski_difference` for a partial inverse operation.
        """
        # default implementation could set this up as an ImageSet

    @abstract_method(optional=True)
    def minkowski_difference(self, other):
        r"""
        Return the Minkowski difference.

        Minkowski subtraction can equivalently be defined via
        Minkowski addition (see :meth:`minkowski_sum`) or as
        set-theoretic intersection via

        .. MATH::

            X \ominus Y =
            (X^c \oplus Y)^c =
            \cap_{y\in Y} (X-y)

        where superscript-"c" means the complement in the ambient
        vector space. The Minkowski difference of convex sets is
        convex, and the difference of polyhedra is again a
        polyhedron. We only consider the case of polyhedra in the
        following. Note that it is not quite the inverse of
        addition. In fact:

        * `(X+Y)-Y = X` for any polyhedra `X`, `Y`.

        * `(X-Y)+Y \subseteq X`

        * `(X-Y)+Y = X` if and only if Y is a Minkowski summand of X.
        """

    @abstract_method(optional=True)
    def is_minkowski_summand(self, Y):
        r"""
        Test whether ``Y`` is a Minkowski summand.

        See :meth:`minkowski_sum`.

        OUTPUT:

        Boolean. Whether there exists another convex set `Z` such that
        ``self`` can be written as `Y\oplus Z`.
        """

    @abstract_method(optional=True)
    def translation(self, displacement):
        # but AffineSpace.translation returns an automorphism
        """
        Return the translated convex set.

        INPUT:

        - ``displacement`` -- a displacement vector or a list/tuple of
          coordinates that determines a displacement vector

        OUTPUT:

        The translated convex set.
        """

    @abstract_method(optional=True)
    def dilation(self, scalar):
        """
        Return the dilated (uniformly stretched) convex set.
        """

    @abstract_method(optional=True)
    def volume(self, measure='ambient'):
        # but IntegerLattice.volume
        """
        Return the volume of ``self``.

        INPUT:

        - ``measure`` -- string. The measure to use. Allowed values are:

          * ``ambient`` (default): Lebesgue measure of ambient space (volume)
          * ``induced``: Lebesgue measure of the affine hull (relative volume)
        """

    @abstract_method(optional=True)
    def integrate(self, function, measure='ambient'):
        r"""
        Return the integral of ``function`` over ``self``.

        INPUT:

        - ``self`` -- Polyhedron

        - ``function`` -- a function

        - ``measure`` -- string, the measure to use

          Allowed values are:

          * ``ambient`` (default): Lebesgue measure of ambient space,
          * ``induced``: Lebesgue measure of the affine hull
        """

    # Other Polyhedron_base methods:

    #def ambient_space
    #def linear_transformation
    #def convex_hull
    #def affine_hull_projection
    #def as_manifold_subset
    #def integral_points
    #def integral_points_count

    # Minkowski-sum/dilation related operators: __add__, __sub__, _mul_, __truediv__, _acted_upon_
    # Subset-related operator: __and__, _richcmp_


    # Other ConvexRationalPolyhedralCone methods:

    #def ambient
    #def is_strictly_convex
    #def solid_restriction


    # Other LatticePolytope methods:

    #def affine_transform
    #def ambient
    #def points, interior_points, boundary_points ... all lattice points


    # Other FreeModule methods:

    #def scale
    #def ambient_vector_space
    #def vector_space

class ConvexSet_closed(ConvexSet_base):
    r"""
    Abstract base class for closed convex sets.
    """

    def is_closed(self):
        r"""
        Return whether ``self`` is closed.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: hcube = polytopes.hypercube(5)
            sage: hcube.is_closed()
            True
        """
        return True

    def is_open(self):
        r"""
        Return whether ``self`` is open.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: hcube = polytopes.hypercube(5)
            sage: hcube.is_open()
            False

            sage: zerocube = polytopes.hypercube(0)
            sage: zerocube.is_open()
            True
        """
        return self.is_empty() or self.is_universe()


class ConvexSet_compact(ConvexSet_closed):
    r"""
    Abstract base class for compact convex sets.
    """

    def is_universe(self):
        r"""
        Return whether ``self`` is the whole ambient space

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: cross3 = lattice_polytope.cross_polytope(3)
            sage: cross3.is_universe()
            False
            sage: point0 = LatticePolytope([[]]); point0
            0-d reflexive polytope in 0-d lattice M
            sage: point0.is_universe()
            True
        """
        return self.ambient_dim() == 0 and not self.is_empty()

    def is_compact(self):
        r"""
        Return whether ``self`` is compact.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: cross3 = lattice_polytope.cross_polytope(3)
            sage: cross3.is_compact()
            True
        """
        return True

    is_relatively_open = ConvexSet_closed.is_open


class ConvexSet_relatively_open(ConvexSet_base):
    r"""
    Abstract base class for relatively open convex sets.
    """

    def is_relatively_open(self):
        r"""
        Return whether ``self`` is relatively open.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: ri_segment = segment.relative_interior()
            sage: ri_segment.is_relatively_open()
            True
        """
        return True

    def is_open(self):
        r"""
        Return whether ``self`` is open.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: segment = Polyhedron([[1, 2], [3, 4]])
            sage: ri_segment = segment.relative_interior()
            sage: ri_segment.is_open()
            False
        """
        return self.is_empty() or self.is_full_dimensional()


class ConvexSet_open(ConvexSet_relatively_open):
    r"""
    Abstract base class for open convex sets.
    """

    def is_open(self):
        r"""
        Return whether ``self`` is open.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.geometry.convex_set import ConvexSet_open
            sage: b = ConvexSet_open()
            sage: b.is_open()
            True
        """
        return True

    def is_closed(self):
        r"""
        Return whether ``self`` is closed.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.geometry.convex_set import ConvexSet_open
            sage: class OpenBall(ConvexSet_open):
            ....:     def dim(self):
            ....:         return 3
            ....:     def is_universe(self):
            ....:         return False
            sage: OpenBall().is_closed()
            False
        """
        return self.is_empty() or self.is_universe()
