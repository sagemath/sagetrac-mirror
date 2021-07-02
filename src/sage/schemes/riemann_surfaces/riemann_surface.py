r"""
Riemann matrices and endomorphism rings of algebraic Riemann surfaces

This module provides a class, :class:`RiemannSurface`, to model the
Riemann surface determined by a plane algebraic curve over a subfield
of the complex numbers.

A homology basis is derived from the edges of a Voronoi cell decomposition based
on the branch locus. The pull-back of these edges to the Riemann surface
provides a graph on it that contains a homology basis.

The class provides methods for computing the Riemann period matrix of the
surface numerically, using a certified homotopy continuation method due to
[Kr2016].

The class also provides facilities for computing the endomorphism ring of the
period lattice numerically, by determining integer (near) solutions to the
relevant approximate linear equations.

AUTHORS:

- Alexandre Zotine, Nils Bruin (2017-06-10): initial version
- Nils Bruin, Jeroen Sijsling (2018-01-05): algebraization, isomorphisms
- Linden Disney-Hogg, Nils Bruin (2021-06-23): efficient integration

EXAMPLES:

We compute the Riemann matrix of a genus 3 curve::

    sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
    sage: R.<x,y> = QQ[]
    sage: f = x^4-x^3*y+2*x^3+2*x^2*y+2*x^2-2*x*y^2+4*x*y-y^3+3*y^2+2*y+1
    sage: S = RiemannSurface(f,prec=100)
    sage: M = S.riemann_matrix()

We test the usual properties, i.e., that the period matrix is symmetric and that
the imaginary part is positive definite::

    sage: all(abs(a) < 1e-20 for a in (M-M.T).list())
    True
    sage: iM = Matrix(RDF,3,3,[a.imag_part() for a in M.list()])
    sage: iM.is_positive_definite()
    True

We compute the endomorphism ring and check it has `\ZZ`-rank 6::

    sage: A = S.endomorphism_basis(80,8)
    sage: len(A) == 6
    True

In fact it is an order in a number field::

    sage: T.<t> = QQ[]
    sage: K.<a> = NumberField(t^6 - t^5 + 2*t^4 + 8*t^3 - t^2 - 5*t + 7)
    sage: all(len(a.minpoly().roots(K)) == a.minpoly().degree() for a in A)
    True
"""
# ****************************************************************************
#       Copyright (C) 2017 Alexandre Zotine, Nils Bruin
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from scipy.spatial import Voronoi
from sage.arith.misc import GCD, algdep
from sage.ext.fast_callable import fast_callable
from sage.graphs.graph import Graph
from sage.groups.matrix_gps.finitely_generated import MatrixGroup
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.matrix.constructor import Matrix
from sage.matrix.special import block_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.misc.misc_c import prod
from sage.modules.free_module import VectorSpace
from sage.numerical.gauss_legendre import integrate_vector, integrate_vector_N
from sage.rings.complex_mpfr import ComplexField, CDF
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.qqbar import number_field_elements_from_algebraics
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RealField
import sage.libs.mpmath.all as mpall


def voronoi_ghost(cpoints, n=6, CC=CDF):
    r"""
    Convert a set of complex points to a list of real tuples `(x,y)`, and
    appends n points in a big circle around them.

    The effect is that, with n >= 3, a Voronoi decomposition will have only
    finite cells around the original points. Furthermore, because the extra
    points are placed on a circle centered on the average of the given points,
    with a radius 3/2 times the largest distance between the center and the
    given points, these finite cells form a simply connected region.

    INPUT:

    - ``cpoints`` -- a list of complex numbers

    OUTPUT:

    A list of real tuples `(x,y)` consisting of the original points and a set of
    points which surround them.

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import voronoi_ghost
        sage: L = [1 + 1*I, 1 - 1*I, -1 + 1*I, -1 - 1*I]
        sage: voronoi_ghost(L) # abs tol 1e-6
        [(1.0, 1.0),
         (1.0, -1.0),
         (-1.0, 1.0),
         (-1.0, -1.0),
         (2.121320343559643, 0.0),
         (1.0606601717798216, 1.8371173070873836),
         (-1.060660171779821, 1.8371173070873839),
         (-2.121320343559643, 2.59786816870648e-16),
         (-1.0606601717798223, -1.8371173070873832),
         (1.06066017177982, -1.8371173070873845)]
    """
    cpoints = [CC(c) for c in cpoints]
    average = sum(cpoints) / len(cpoints)
    if len(cpoints) == 1:
        radius = 1
    else:
        radius = 3 * max(abs(c - average) for c in cpoints) / 2
    z = CC.zeta(n)
    extra_points = [average + radius * z**i for i in range(n)]
    return [tuple(c) for c in cpoints + extra_points]


def bisect(L, t):
    r"""
    Find position in a sorted list using bisection.

    Given a list `L = [(t_0,...),(t_1,...),...(t_n,...)]` with increasing `t_i`,
    find the index i such that `t_i <= t < t_{i+1}` using bisection. The rest of
    the tuple is available for whatever use required.

    INPUT:

    - ``L`` -- A list of tuples such that the first term of each tuple is a real
      number between 0 and 1. These real numbers must be increasing.

    - ``t`` -- A real number between `t_0` and `t_n`.

    OUTPUT:

    An integer i, giving the position in L where t would be in

    EXAMPLES:

    Form a list of the desired form, and pick a real number between 0 and 1::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import bisect
        sage: L = [(0.0, 'a'), (0.3, 'b'), (0.7, 'c'), (0.8, 'd'), (0.9, 'e'), (1.0, 'f')]
        sage: t = 0.5
        sage: bisect(L,t)
        1

    Another example which demonstrates that if t is equal to one of the t_i, it
    returns that index::

        sage: L = [(0.0, 'a'), (0.1, 'b'), (0.45, 'c'), (0.5, 'd'), (0.65, 'e'), (1.0, 'f')]
        sage: t = 0.5
        sage: bisect(L,t)
        3
    """
    # Defining starting indices for the loop.
    min = 0
    max = len(L) - 1
    # If the input t is not between 0 and 1, raise an error.
    if t < L[min][0] or t > L[max][0]:
        raise ValueError("value for t out of range")
    # Main loop.
    while (min < max-1):
        # Bisect.
        mid = (max+min)//2
        # If it's equal, return the index we bisected to.
        if t == L[mid][0]:
            return mid
        # If it's smaller, then we're on the left side.
        elif t < L[mid][0]:
            max = mid
        # Otherwise we're on the right side.
        else:
            min = mid
    # Once the loop terminates, we return what the indices converged to.
    return min


def numerical_inverse(C):
    """
    Compute numerical inverse of a matrix via LU decomposition

    INPUT:

    - ``C`` -- A real or complex invertible square matrix

    EXAMPLES::

        sage: C = matrix(CC,3,3,[-4.5606e-31 + 1.2326e-31*I,
        ....: -0.21313 + 0.24166*I,
        ....: -3.4513e-31 + 0.16111*I,
        ....: -1.0175 + 9.8608e-32*I,
        ....: 0.30912 + 0.19962*I,
        ....: -4.9304e-32 + 0.39923*I,
        ....: 0.96793 - 3.4513e-31*I,
        ....: -0.091587 + 0.19276*I,
        ....: 3.9443e-31 + 0.38552*I])
        sage: from sage.schemes.riemann_surfaces.riemann_surface import numerical_inverse
        sage: max(abs(c) for c in (C^(-1)*C-C^0).list()) < 1e-10
        False
        sage: max(abs(c) for c in (numerical_inverse(C)*C-C^0).list()) < 1e-10
        True
    """
    R = C.parent()
    prec = R.base_ring().prec()
    mpall.mp.prec = prec
    with mpall.workprec(prec):
        Cmp = mpall.matrix([mpall.sage_to_mpmath(list(c), prec) for c in C])
        PLU = mpall.lu(Cmp)
    P, L, U = [R([mpall.mpmath_to_sage(c, prec) for c in M]) for M in PLU]
    return U.inverse() * L.inverse() * P


class ConvergenceError(ValueError):
    r"""
    Error object suitable for raising and catching when Newton iteration fails.

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import ConvergenceError
        sage: raise ConvergenceError("test")
        Traceback (most recent call last):
        ...
        ConvergenceError: test
        sage: isinstance(ConvergenceError(),ValueError)
        True
    """
    pass


def differential_basis_baker(f):
    r"""
    Compute a differential bases for a curve that is nonsingular outside (1:0:0),(0:1:0),(0:0:1)

    Baker's theorem tells us that if a curve has its singularities at the coordinate vertices and meets
    some further easily tested genericity criteria,
    then we can read off a basis for the regular differentials from the interior of the
    Newton polygon spanned by the monomials. While this theorem only applies to special plane curves
    it is worth implementing because the analysis is relatively cheap and it applies to a lot of
    commonly encountered curves (e.g., curves given by a hyperelliptic model). Other advantages include
    that we can do the computation over any exact base ring (the alternative Singular based method for
    computing the adjoint ideal requires the rationals), and that we can avoid being affected by subtle bugs
    in the Singular code.

    ``None`` is returned when ``f`` does not describe a curve of the relevant type. If ``f`` is of the relevant
    type, but is of genus `0` then ``[]`` is returned (which are both False values, but they are not equal).

    INPUT:

    - `f` -- a bivariate polynomial

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import differential_basis_baker
        sage: R.<x,y> = QQ[]
        sage: f = x^3+y^3+x^5*y^5
        sage: differential_basis_baker(f)
        [y^2, x*y, x*y^2, x^2, x^2*y, x^2*y^2, x^2*y^3, x^3*y^2, x^3*y^3]
        sage: f = y^2-(x-3)^2*x
        sage: differential_basis_baker(f) is None
        True
        sage: differential_basis_baker(x^2+y^2-1)
        []

    TESTS::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import differential_basis_baker
        sage: R.<x,y> = QQ[]
        sage: f = y^12 - x*(x - 1)^7
        sage: differential_basis_baker(f) is None
        True

    """
    k = f.base_ring()
    R = PolynomialRing(k, 3, "x,y,z")
    x, y, z = R.gens()
    F = f(x / z, y / z).numerator()
    W = [F] + [F.derivative(v) for v in R.gens()]
    # we check that the singularities lie at (1:0:0),(0:1:0),(0:0:1)
    # by checking that the eliminations of x, y, z result in
    # (principal) ideals generated by a monomial. This is a sufficient
    # condition, but not completely necessary.
    # It's cheap to check, though.
    for c in R.gens():
        B = GCD([W[i].resultant(W[j], c) for i in range(4) for j in range(i)])
        if len(B.monomials()) > 1:
            return None
    from sage.geometry.polyhedron.constructor import Polyhedron
    D = {(k[0], k[1]): v for k, v in f.dict().items()}
    P = Polyhedron(D)
    kT = k['t']
    # here we check the additional genericity conditions: that the polynomials
    # along the edges of the Newton polygon are square-free.
    for e in P.bounded_edges():
        h = kT([D.get(tuple(c), 0) for c in Polyhedron(e).integral_points()])
        if not h.is_squarefree():
            return None
    x, y = f.parent().gens()
    return [x**(a[0] - 1) * y**(a[1] - 1) for a in P.integral_points()
            if P.interior_contains(a)]


class RiemannSurface(object):
    r"""
    Construct a Riemann Surface. This is specified by the zeroes of a bivariate
    polynomial with rational coefficients `f(z,w) = 0`.

    INPUT:

    - ``f`` -- a bivariate polynomial with rational coefficients. The surface is
      interpreted as the covering space of the coordinate plane in the first
      variable.

    - ``prec`` -- the desired precision of computations on the surface in bits
      (default: 53)

    - ``certification`` -- a boolean (default: True) value indicating whether
      homotopy continuation is certified or not. Uncertified homotopy
      continuation can be faster.

    - ``differentials`` -- (default: None). If specified, provides a list of
      polynomials `h` such that `h/(df/dw) dz` is a regular differential on the
      Riemann surface. This is taken as a basis of the regular differentials, so
      the genus is assumed to be equal to the length of this list. The results
      from the homology basis computation are checked against this value.
      Providing this parameter makes the computation independent from Singular.
      For a nonsingular plane curve of degree `d`, an appropriate set is given
      by the monomials of degree up to `d-3`.

    - ``integration_method`` -- (default: ``'heuristic'``). String specifing the 
      integration method to use when calculating the integrals of differentials. 
      The options are ``'heuristic'`` and ``'rigorous'``, the latter of
      which is often the most efficient. 

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
        sage: R.<z,w> = QQ[]
        sage: f = w^2 - z^3 + 1
        sage: RiemannSurface(f)
        Riemann surface defined by polynomial f = -z^3 + w^2 + 1 = 0, with 53 bits of precision

    Another Riemann surface with 100 bits of precision::

        sage: S = RiemannSurface(f, prec=100); S
        Riemann surface defined by polynomial f = -z^3 + w^2 + 1 = 0, with 100 bits of precision
        sage: S.riemann_matrix()^6 #abs tol 0.00000001
        [1.0000000000000000000000000000 - 1.1832913578315177081175928479e-30*I]

    We can also work with Riemann surfaces that are defined over fields with a
    complex embedding, but since the current interface for computing genus and
    regular differentials in Singular presently does not support extensions of
    QQ, we need to specify a description of the differentials ourselves. We give
    an example of a CM elliptic curve::

        sage: Qt.<t> = QQ[]
        sage: K.<a> = NumberField(t^2-t+3,embedding=CC(0.5+1.6*I))
        sage: R.<x,y> = K[]
        sage: f = y^2+y-(x^3+(1-a)*x^2-(2+a)*x-2)
        sage: S = RiemannSurface(f,prec=100,differentials=[1])
        sage: A = S.endomorphism_basis()
        sage: len(A)
        2
        sage: all( len(T.minpoly().roots(K)) > 0 for T in A)
        True

    TESTS:

    This elliptic curve has a relatively poorly conditioned set of branch
    points, so it challenges the path choice a bit. The code just verifies that
    the period is quadratic, because the curve has CM, but really the test is
    that the computation completes at all.::

        sage: prec = 50
        sage: Qx.<t> = QQ[]
        sage: CC = ComplexField(prec)
        sage: g = t^2-t-1
        sage: phiCC = g.roots(CC)[1][0]
        sage: K.<phi> = NumberField(g, embedding=phiCC)
        sage: R.<X,Y> = K[]
        sage: f = Y^2+X*Y+phi*Y-(X^3-X^2-2*phi*X+phi)
        sage: S = RiemannSurface(f,prec=prec, differentials=[1])
        sage: tau = S.riemann_matrix()[0, 0]
        sage: tau.algdep(6).degree() == 2
        True
    """
    def __init__(self, f, prec=53, certification=True, differentials=None, integration_method="heuristic"):
        r"""
        TESTS::

            sage: R.<z,w> = QQ[]
            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: S = RiemannSurface(w^2 - z^3 + 1)
            sage: TestSuite(S).run() #not tested; Unclear what pickling strategy is best.
        """
        # Initializations.
        self._prec = prec
        self._certification = certification
        self._integration_method = integration_method
        self._R = f.parent()
        if len(self._R.gens()) != 2:
            raise ValueError('only bivariate polynomials supported.')
        if f.degree() <= 1:
            raise ValueError('equation must be of degree at least 2.')
        z, w = self._R.gen(0), self._R.gen(1)
        self._CC = ComplexField(self._prec)
        self._RR = RealField(self._prec)
        self._CCz = PolynomialRing(self._CC, [self._R.gen(0)])
        self._CCw = PolynomialRing(self._CC, [self._R.gen(1)])
        self._RRz = PolynomialRing(self._RR, [self._R.gen(0)])
        self.f = f
        if differentials is not None:
            self._differentials = [self._R(a) for a in differentials]
            self.genus = len(self._differentials)
        else:
            B = differential_basis_baker(f)
            if B is not None:
                self._differentials = B
                self.genus = len(B)
            else:
                self._differentials = None
                self.genus = self._R.ideal(self.f).genus()
                if self.genus < 0:
                    raise ValueError("Singular reports negative genus. Specify differentials manually.")
        self.degree = self.f.degree(w)
        self._dfdw = self.f.derivative(w)
        self._dfdz = self.f.derivative(z)
        self._discriminant = self.f.resultant(self._dfdw, w)
        # Coefficients of the polynomial for use in homotopy continuation.
        self._a0 = self._CCz(self.f.coefficient({w: self.degree})(self._CCz.gen(), 0))
        self._a0roots = self._a0.roots(multiplicities=False)
        self._aks = [self._CCz(self.f.coefficient({w: self.degree - k - 1})
                               (self._CCz.gen(), 0)) for k in range(self.degree)]
        # Compute the branch locus. Takes the square-free part of the discriminant
        # because of numerical issues.
        self.branch_locus = []
        for x in self._discriminant.factor():
            self.branch_locus += self._CCz(x[0](self._CCz.gen(), 0)).roots(multiplicities=False)
        # Voronoi diagram and the important points associated with it
        self.voronoi_diagram = Voronoi(voronoi_ghost(self.branch_locus,
                                                     CC=self._CC))
        self._vertices = [self._CC(x0, y0)
                          for x0, y0 in self.voronoi_diagram.vertices]
        self._wvalues = [self.w_values(z0) for z0 in self._vertices]
        self._Sn = SymmetricGroup(range(self.degree))
        self._L = {}
        self._fastcall_f = fast_callable(f, domain=self._CC)
        self._fastcall_dfdw = fast_callable(self._dfdw, domain=self._CC)
        self._fastcall_dfdz = fast_callable(self._dfdz, domain=self._CC)

    def __repr__(self):
        r"""
        Return a string representation of the Riemann surface class.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: RiemannSurface(f)
            Riemann surface defined by polynomial f = -z^4 + w^2 + 1 = 0, with 53 bits of precision
        """
        s = 'Riemann surface defined by polynomial f = %s = 0, with %s bits of precision' % (self.f, self._prec)
        return s

    def w_values(self, z0):
        r"""
        Return the points lying on the surface above ``z0``.

        INPUT:

        - ``z0`` -- (complex) a point in the complex z-plane.

        OUTPUT:

        A set of complex numbers corresponding to solutions of `f(z_0,w) = 0`.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)

        Find the w-values above the origin, i.e. the solutions of `w^2 + 1 = 0`::

            sage: S.w_values(0) # abs tol 1e-14
            [-1.00000000000000*I, 1.00000000000000*I]
        """
        return self.f(z0,self._CCw.gen(0)).roots(multiplicities=False)

    @cached_method
    def downstairs_edges(self):
        r"""
        Compute the edgeset of the Voronoi diagram.

        OUTPUT:

        A list of integer tuples corresponding to edges between vertices in the
        Voronoi diagram.

        EXAMPLES:

        Form a Riemann surface, one with a particularly simple branch locus::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 + z^3 - z^2
            sage: S = RiemannSurface(f)

        Compute the edges::

            sage: S.downstairs_edges()
            [(0, 1), (0, 5), (1, 4), (2, 3), (2, 4), (3, 5), (4, 5)]

        This now gives an edgeset which one could use to form a graph.

        .. NOTE::

            The numbering of the vertices is given by the Voronoi package.
        """
        # Because of how we constructed the Voronoi diagram, the first n points
        # correspond to the branch locus points.
        # The regions of these points are all of the edges which don't go off
        # to infinity, which are exactly the ones we want.
        n = len(self.branch_locus)
        desired_edges = [self.voronoi_diagram.regions[self.voronoi_diagram.point_region[i]] for i in range(n)]
        # First construct the edges as a set because the regions will overlap
        # and we don't want to have two of the same edge.
        edges1 = set()
        for c in desired_edges:
            for j in range(len(c)-1):
                edges1.add(frozenset((c[j],c[j+1])))
            edges1.add(frozenset((c[0],c[-1])))
        # Then make it into a list and sort it.
        # The sorting is important - it will make computing the monodromy group
        # MUCH easier.
        # We orient all the edges so that we go from lower to higher
        # numbered vertex for the continuation.
        edges = [(i0,i1) if (i0 < i1) else (i1,i0) for (i0,i1) in edges1]
        edges.sort()
        return edges

    def downstairs_graph(self):
        r"""
        Return the Voronoi decomposition as a planar graph.

        The result of this routine can be useful to interpret the labelling of
        the vertices.

        OUTPUT:

        The Voronoi decomposition as a graph, with appropriate planar embedding.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: S.downstairs_graph()
            Graph on 11 vertices

        Similarly one can form the graph of the upstairs edges, which is
        visually rather less attractive but can be instructive to verify that a
        homology basis is likely correctly computed.::

            sage: G = Graph(S.upstairs_edges()); G
            Graph on 22 vertices
            sage: G.is_planar()
            False
            sage: G.genus()
            1
            sage: G.is_connected()
            True
        """
        G = Graph(self.downstairs_edges())
        G.set_pos(dict(enumerate(list(v) for v in self._vertices)))
        return G

    def _compute_delta(self, z1, epsilon, wvalues=None):
        r"""
        Compute a delta for homotopy continuation when moving along a path.

        INPUT:

        - ``z1`` -- a complex number in the z-plane

        - ``epsilon`` -- a real number, which is the minimum distance between
          the w-values above ``z1``

        - ``wvalues`` -- a list (default: ``None``). If specified, saves
          recomputation.

        OUTPUT:

        A real number, which is a step size for moving along a path.

        EXAMPLES:

        Form a Riemann Surface::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)

        Pick a point which lies on the Voronoi diagram, and compute an
        appropriate epsilon::

            sage: z1 = S._vertices[0]
            sage: currw = S.w_values(z1)
            sage: n = len(currw)
            sage: epsilon = min([abs(currw[i] - currw[n-j-1]) for i in range(n) for j in range(n-i-1)])/3
            sage: S._compute_delta(z1, epsilon) # abs tol 1e-8
            0.152628501142363

        If the Riemann surface does not have certified homotopy continuation,
        then the delta will just be the minimum distance away from a branch
        point::

            sage: T = RiemannSurface(f, certification=False)
            sage: z1 = T._vertices[0]
            sage: currw = T.w_values(z1)
            sage: n = len(currw)
            sage: epsilon = min([abs(currw[i] - currw[n-j-1]) for i in range(n) for j in range(n-i-1)])/3
            sage: T._compute_delta(z1, epsilon) # abs tol 1e-8
            0.381881307912987
        """
        if self._certification:
            if wvalues is None:
                wvalues = self.w_values(z1)
            # For computation of rho. Need the branch locus + roots of a0.
            badpoints = self.branch_locus + self._a0roots
            rho = min(abs(z1 - z) for z in badpoints) / 2
            Y = max(abs(self._fastcall_dfdz(z1, wi)/self._fastcall_dfdw(z1, wi))
                    for wi in wvalues)

            # compute M
            upperbounds = [sum(ak[k] * (abs(z1) + rho)**k
                               for k in range(ak.degree()))
                           for ak in self._aks]
            upperbounds.reverse()
            # If a0 is a constant polynomial, it is obviously bounded below.
            if not self._a0roots:
                lowerbound = self._CC(self._a0) / 2
            else:
                lowerbound = self._a0[self._a0.degree()]*prod(abs((zk - z1) - rho) for zk in self._a0roots) / 2
            M = 2 * max((upperbounds[k]/lowerbound).abs().nth_root(k+1)
                        for k in range(self.degree-1))
            return rho*(((rho*Y - epsilon)**2 + 4*epsilon*M).sqrt() - (rho*Y + epsilon))/(2*M - 2*rho*Y)
        else:
            # Instead, we just compute the minimum distance between branch
            # points and the point in question.
            return min(abs(b - z1) for b in self.branch_locus) / 2

    def homotopy_continuation(self, edge):
        r"""
        Perform homotopy continuation along an edge of the Voronoi diagram using
        Newton iteration.

        INPUT:

        - ``edge`` -- a tuple of integers indicating an edge of the Voronoi
          diagram

        OUTPUT:

        A list of complex numbers corresponding to the points which are reached
        when traversing along the direction of the edge. The ordering of these
        points indicates how they have been permuted due to the weaving of the
        curve.

        EXAMPLES:

        We check that continued values along an edge correspond (up to the
        appropriate permutation) to what is stored. Note that the permutation
        was originally computed from this data::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)
            sage: edge1 = sorted(S.edge_permutations())[0]
            sage: sigma = S.edge_permutations()[edge1]
            sage: continued_values = S.homotopy_continuation(edge1)
            sage: stored_values = S.w_values(S._vertices[edge1[1]])
            sage: all( abs(continued_values[i]-stored_values[sigma(i)]) < 1e-8 for i in range(3))
            True
        """
        i0, i1 = edge
        ZERO = self._RR.zero()
        ONE = self._RR.one()
        datastorage = []
        z_start = self._CC(self._vertices[i0])
        z_end = self._CC(self._vertices[i1])
        path_length = abs(z_end - z_start)

        def path(t):
            return z_start * (1 - t) + z_end * t
        # Primary procedure.
        T = ZERO
        currw = self.w_values(path(T))
        n = len(currw)
        epsilon = min([abs(currw[i] - currw[j]) for i in range(1,n) for j in range(i)])/3
        datastorage += [(T,currw,epsilon)]
        while T < ONE:
            delta = self._compute_delta(path(T), epsilon, wvalues=currw)/path_length
            # Move along the path by delta.
            T += delta
            # If T exceeds 1, just set it to 1 and compute.
            if T > ONE:
                delta -= (T-ONE)
                T = ONE
            while True:
                try:
                    neww = self._determine_new_w(path(T),currw,epsilon)
                except ConvergenceError:
                    delta /= 2
                    T -= delta
                else:
                    break
            currw = neww
            epsilon = min([abs(currw[i] - currw[j]) for i in range(1,n) for j in range(i)])/3
            datastorage += [(T,currw,epsilon)]
        self._L[edge] = datastorage
        return currw

    def _determine_new_w(self, z0, oldw, epsilon):
        r"""
        A procedure to Newton iterate a list of w-values simultaneously.

        Used primarily for moving along the surface for integration or
        homotopy continuation.

        INPUT:

        - ``z0`` -- a complex number

        - ``oldw`` -- a list of w-values which are presumed to be guesses of
          the w-values above ``z0``.

        - ``epsilon`` -- the minimum distance between the points of ``oldw``
          divided by 3

        OUTPUT:

        A list of points the same length as ``oldw`` corresponding to the new
        Newton iterated points.

        However, if the Newton iteration exceeds the allotted attempts,
        or exits the ``epsilon`` ball, raises a convergence error.

        EXAMPLES:

        First, a trivial example where we guess exactly what the roots are::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: z0 = S._vertices[0]
            sage: epsilon = 0.1
            sage: oldw = S.w_values(z0)
            sage: neww = S._determine_new_w(z0,oldw,epsilon); neww #abs tol 0.00000001
                [-0.934613146929672 + 2.01088055918363*I,
                 0.934613146929672 - 2.01088055918363*I]

        Which should be exactly the same as the w-values we started with.::

            sage: abs(neww[0] - oldw[0]) #abs tol 0.00000001
            0.000000000000...
            sage: abs(neww[1] - oldw[1]) #abs tol 0.00000001
            0.000000000000...

        Here is an example where we exit the ``epsilon`` bound. This approach is
        based on the homotopy continuation procedure which traverses along a
        path and attempts Newton iteration::

            sage: g = z^3*w + w^3 + z
            sage: T = RiemannSurface(g)
            sage: z0 = T._vertices[2]*(0.9) - T._vertices[15]*(0.1)
            sage: epsilon = 0.5
            sage: oldw = T.w_values(T._vertices[2])
            sage: T._determine_new_w(z0,oldw,epsilon)
            [-0.562337685361648 + 0.151166007149998*I,
             0.640201585779414 - 1.48567225836436*I,
             -0.0778639004177661 + 1.33450625121437*I]
        """
        # Tools of Newton iteration.
        F = self._fastcall_f
        dF = self._fastcall_dfdw
        neww = []
        prec = self._CC.prec()
        # Iterate over all roots.
        for i in range(len(oldw)):
            delta = F(z0, oldw[i]) / dF(z0, oldw[i])
            Ndelta = delta.norm()
            wi = oldw[i] - delta
            # it is possible in theory that Newton iteration fails to
            # converge without escaping. We catch this by capping the
            # number of iterations by 100
            for j in range(100):
                # If we exceed the epsilon bound from homotopy continuation,
                # terminate.
                if abs(wi - oldw[i]) >= epsilon:
                    raise ConvergenceError("Newton iteration escaped neighbourhood")
                new_delta = F(z0, wi) / dF(z0, wi)
                Nnew_delta = new_delta.norm()
                # If we found the root exactly, or if delta only affects half the digits and
                # stops getting smaller, we decide that we have converged.
                if (new_delta == 0) or (Nnew_delta >= Ndelta and
                            Ndelta.sign_mantissa_exponent()[2]+prec < wi.norm().sign_mantissa_exponent()[2]):
                    neww.append(wi)
                    break
                delta = new_delta
                Ndelta = Nnew_delta
                wi -= delta
            # If we run 100 iterations without a result, terminate.
            else:
                raise ConvergenceError("Newton iteration fails to converge after %s iterations" % j)
        return neww

    def _newton_iteration(self, z0, oldw, epsilon):
        r"""
        A non-vectorized Newton iteration procedure used for integration.

        INPUT:

        - ``z0`` -- a complex number.

        - ``oldw`` -- a w-value which is presumed to be a guess of one of
          the w-values above ``z0``.

        - ``epsilon`` -- the minimum distance between the w-values divided by 3.

        OUTPUT:

        A complex number, which should be a w-value above ``z0``.

        However, if the Newton iteration exceeds the allotted attempts,
        or exits the ``epsilon`` ball, raises a convergence error.

        EXAMPLES:

        First, a trivial example where we guess exactly what the root is::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: z0 = S._vertices[0]
            sage: epsilon = 0.1
            sage: oldw = S.w_values(z0)[0]
            sage: neww = S._newton_iteration(z0,oldw,epsilon); neww #abs tol 0.00000001
            -0.934613146929672 + 2.01088055918363*I

        Which should be exactly the same as the w-value we started with::

            sage: oldw - neww #abs tol 0.00000001
            0.000000000000000

        Here is an example where we exit the epsilon bound. This approach is
        based on the homotopy continuation procedure which traverses along a
        path and attempts Newton iteration::

            sage: g = z^3*w + w^3 + z
            sage: T = RiemannSurface(g)
            sage: z0 = T._vertices[2]*(0.9) - T._vertices[15]*(0.1)
            sage: epsilon = 0.5
            sage: oldw = T.w_values(T._vertices[2])[0]
            sage: T._newton_iteration(z0, oldw, epsilon)
            -0.562337685361648 + 0.151166007149998*I
        """
        F = self._fastcall_f
        dF = self._fastcall_dfdw
        prec = self._CC.prec()
        delta = F(z0, oldw) / dF(z0, oldw)
        Ndelta = delta.norm()
        neww = oldw - delta
        eps_squared = epsilon**2
        # it is possible in theory that Newton iteration fails to converge
        # without escaping. We catch this by capping the number of iterations
        # by 100
        for j in range(100):
            if (neww - oldw).norm() > eps_squared:
                raise ConvergenceError("Newton iteration escaped neighbourhood")
            new_delta = F(z0, neww) / dF(z0, neww)
            Nnew_delta = new_delta.norm()
            # If we found the root exactly, or if delta only affects half the digits and
            # stops getting smaller, we decide that we have converged.
            if (new_delta == 0) or (Nnew_delta>=Ndelta and
                    Ndelta.sign_mantissa_exponent()[2]+prec < neww.norm().sign_mantissa_exponent()[2]):
                return neww
            delta = new_delta
            Ndelta = Nnew_delta
            neww-=delta
        raise ConvergenceError("Newton iteration fails to converge")

    @cached_method
    def upstairs_edges(self):
        r"""
        Compute the edgeset of the lift of the downstairs graph onto the Riemann
        surface.

        OUTPUT:

        An edgeset between vertices (i, j), where i corresponds to the i-th
        point in the Voronoi diagram vertices, and j is the j-th w-value
        associated with that point.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 + z^3 - z^2
            sage: S = RiemannSurface(f)
            sage: edgeset = S.upstairs_edges()
            sage: len(edgeset) == S.degree*len(S.downstairs_edges())
            True
            sage: {(v[0],w[0]) for v,w in edgeset} == set(S.downstairs_edges())
            True
        """
        edgeset = []
        n = len(self._wvalues[0])
        # Lifts each edge individually.
        for e in self.downstairs_edges():
            i0, i1 = e
            # Epsilon for checking w-value later.
            epsilon = min([abs(self._wvalues[i1][i] - self._wvalues[i1][n-j-1]) for i in range(n) for j in range(n-i-1)])/3
            # Homotopy continuation along e.
            homotopycont = self.homotopy_continuation(e)
            for i in range(len(homotopycont)):
                # Checks over the w-values of the next point to check which it is.
                for j in range(len(self._wvalues[i1])):
                    if abs(homotopycont[i] - self._wvalues[i1][j]) < epsilon:
                        # Once it finds the appropriate w-value, adds the edge.
                        edgeset = edgeset + [[(i0, i), (i1, j)]]
                        continue
        return edgeset

    def _edge_permutation(self, edge):
        r"""
        Compute the permutation of the w-values above a point in the z-plane
        when moving along an edge in the Voronoi diagram.

        INPUT:

        - ``edge`` -- an edge on the Voronoi diagram

        OUTPUT:

        A permutation corresponding to how the roots interchange when moving
        along the edge.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)

        Compute the edge permutation of (1,2) on the Voronoi diagram::

            sage: S._edge_permutation((1,2))
            (0,2,1)

        This indicates that while traversing along the direction of `(5,16)`,
        the 2nd and 3rd layers of the Riemann surface are interchanging.
        """
        if edge in self.downstairs_edges():
            # find all upstairs edges that are lifts of the given
            # downstairs edge and store the corresponding indices at
            # start and end that label the branches upstairs.
            L = [(j0, j1) for ((i0, j0), (i1, j1)) in self.upstairs_edges()
                 if edge == (i0, i1)]
            # we should be finding exactly "degree" of these
            assert len(L) == self.degree
            # and as a corollary of how we construct them, the indices
            # at the start should be in order
            assert all(a == b[0] for a, b in enumerate(L))
            return self._Sn([j1 for j0, j1 in L])
        raise ValueError('edge not in Voronoi diagram')

    @cached_method
    def edge_permutations(self) -> dict:
        r"""
        Compute the permutations of branches associated to each edge.

        Over the vertices of the Voronoi decomposition around the branch locus,
        we label the fibres. By following along an edge, the lifts of the edge
        induce a permutation of that labelling.

        OUTPUT:

        A dictionary with as keys the edges of the Voronoi decomposition and as
        values the corresponding permutations.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 + z^2+1
            sage: S = RiemannSurface(f)
            sage: S.edge_permutations()
            {(0, 2): (),
             (0, 4): (),
             (1, 2): (),
             (1, 3): (0,1),
             (1, 6): (),
             (2, 0): (),
             (2, 1): (),
             (2, 5): (0,1),
             (3, 1): (0,1),
             (3, 4): (),
             (4, 0): (),
             (4, 3): (),
             (5, 2): (0,1),
             (5, 7): (),
             (6, 1): (),
             (6, 7): (),
             (7, 5): (),
             (7, 6): ()}
        """
        D = {e: self._edge_permutation(e) for e in self.downstairs_edges()}
        for (a, b), p in list(D.items()):
            D[(b, a)] = p**(-1)
        return D

    @cached_method
    def monodromy_group(self):
        r"""
        Compute local monodromy generators of the Riemann surface.

        For each branch point, the local monodromy is encoded by a permutation.
        The permutations returned correspond to positively oriented loops around
        each branch point, with a fixed base point. This means the generators
        are properly conjugated to ensure that together they generate the global
        monodromy. The list has an entry for every finite point stored in
        ``self.branch_locus``, plus an entry for the ramification above infinity.

        OUTPUT:

        A list of permutations, encoding the local monodromy at each branch
        point.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z, w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)
            sage: G = S.monodromy_group(); G
            [(0,1,2), (0,1), (0,2), (1,2), (1,2), (1,2), (0,1), (0,2), (0,2)]

        The permutations give the local monodromy generators for the branch
        points::

            sage: list(zip(S.branch_locus + [unsigned_infinity], G)) #abs tol 0.0000001
            [(0.000000000000000, (0,1,2)),
             (-1.31362670141929, (0,1)),
             (-0.819032851784253 - 1.02703471138023*I, (0,2)),
             (-0.819032851784253 + 1.02703471138023*I, (1,2)),
             (0.292309440469772 - 1.28069133740100*I, (1,2)),
             (0.292309440469772 + 1.28069133740100*I, (1,2)),
             (1.18353676202412 - 0.569961265016465*I, (0,1)),
             (1.18353676202412 + 0.569961265016465*I, (0,2)),
             (Infinity, (0,2))]

        We can check the ramification by looking at the cycle lengths and verify
        it agrees with the Riemann-Hurwitz formula::

            sage: 2*S.genus-2 == -2*S.degree + sum(e-1 for g in G for e in g.cycle_type())
            True
        """
        n = len(self.branch_locus)
        G = Graph(self.downstairs_edges())
        # we get all the regions
        loops = [self.voronoi_diagram.regions[i][:]
                 for i in self.voronoi_diagram.point_region]
        # and construct their Voronoi centers as complex numbers
        centers = self.branch_locus + [self._CC(x, y) for x, y in self.voronoi_diagram.points[n:]]
        for center, loop in zip(centers, loops):
            if -1 in loop:
                # for loops involving infinity we take the finite part of the path
                i = loop.index(-1)
                loop[:] = loop[i+1:]+loop[:i]
            else:
                # and for finite ones we close the paths
                loop.append(loop[0])
            # we make sure the loops are positively oriented wrt. their center
            v0 = self._vertices[loop[0]]
            v1 = self._vertices[loop[1]]
            M = Matrix([list(v0 - center), list(v1 - center)])
            if M.det() < 0:
                loop.reverse()

        # we stitch together the paths that are part of loops through
        # infinity. There should be a unique way of doing so.
        inf_loops = loops[n:]
        inf_path = inf_loops.pop()
        while (inf_loops):
            inf_path += (inf_loops.pop())[1:]
        assert inf_path[0] == inf_path[-1]

        loops = loops[:n]
        loops.append(inf_path)

        P0 = loops[0][0]
        monodromy_gens = []
        edge_perms = self.edge_permutations()
        SG = self._Sn
        for c in loops:
            to_loop = G.shortest_path(P0, c[0])
            to_loop_perm = SG.prod(edge_perms[(to_loop[i], to_loop[i + 1])]
                                   for i in range(len(to_loop) - 1))
            c_perm = SG.prod(edge_perms[(c[i], c[i + 1])]
                             for i in range(len(c) - 1))
            monodromy_gens.append(to_loop_perm * c_perm * ~to_loop_perm)
        return monodromy_gens

    @cached_method
    def homology_basis(self):
        r"""
        Compute the homology basis of the Riemann surface.

        OUTPUT:

        A list of paths `L = [P_1, \dots, P_n]`. Each path `P_i` is of the form
        `(k, [p_1 ... p_m, p_1])`, where `k` is the number of times to traverse
        the path (if negative, to traverse it backwards), and the `p_i` are
        vertices of the upstairs graph.

        EXAMPLES:

        In this example, there are two paths that form the homology basis::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: g = w^2 - z^4 + 1
            sage: S = RiemannSurface(g)
            sage: S.homology_basis() #random
            [[(1, [(3, 1), (5, 0), (9, 0), (10, 0), (2, 0), (4, 0),
                (7, 1), (10, 1), (3, 1)])],
             [(1, [(8, 0), (6, 0), (7, 0), (10, 0), (2, 0), (4, 0),
                (7, 1), (10, 1), (9, 1), (8, 0)])]]

        In order to check that the answer returned above is reasonable, we
        test some basic properties. We express the faces of the downstairs graph
        as ZZ-linear combinations of the edges and check that the projection
        of the homology basis upstairs projects down to independent linear
        combinations of an even number of faces::

            sage: dg = S.downstairs_graph()
            sage: edges = dg.edges()
            sage: E = ZZ^len(edges)
            sage: edge_to_E = { e[:2]: E.gen(i) for i,e in enumerate(edges)}
            sage: edge_to_E.update({ (e[1],e[0]): -E.gen(i) for i,e in enumerate(edges)})
            sage: face_span = E.submodule([sum(edge_to_E[e] for e in f) for f in dg.faces()])
            sage: def path_to_E(path):
            ....:     k,P = path
            ....:     return k*sum(edge_to_E[(P[i][0],P[i+1][0])] for i in range(len(P)-1))
            sage: hom_basis = [sum(path_to_E(p) for p in loop) for loop in S.homology_basis()]
            sage: face_span.submodule(hom_basis).rank()
            2
            sage: [sum(face_span.coordinate_vector(b))%2 for b in hom_basis]
            [0, 0]
        """
        if self.genus == 0:
            return []

        edgesu = self.upstairs_edges()
        cycles = Graph(edgesu).cycle_basis()
        # Computing the Gram matrix.
        cn = len(cycles)
        # Forming a list of lists of zeroes.
        # Later this will be converted into a matrix.
        intersectionprod = [[0] * cn for _ in cycles]

        # as it turns out, in extreme examples argument computation
        # can be quite dominant so we cache this (since we may end up
        # using these values multiple times)
        direction_cache = {}

        def direction(center, neighbour):
            k = (center, neighbour)
            if k not in direction_cache:
                theta = (self._vertices[neighbour] - self._vertices[center]).argument()
                direction_cache[k] = theta
                return theta
            else:
                return direction_cache[k]

        # This loop will start at the entry (0,1), and proceed along the row up
        # til (0,cn-1).
        # Then it will go to entry (1,2), and proceed along the row, etc.
        for i in range(1, cn):
            for j in range(i):
                # Initializing the intersection product value.
                intsum = 0
                # Intersection of the edges
                intsec = set(cycles[i]).intersection(set(cycles[j]))
                for v in intsec:
                    # Get indices of the vertex in the cycles.
                    i0 = cycles[i].index(v)
                    i1 = cycles[j].index(v)
                    # Get the complex value of the vertex v.
                    center = cycles[i][i0][0]

                    # We are in the following situation:
                    # We have two paths a_in->v->a_out and
                    # b_in->v->b_out intersecting. We say they
                    # are "positively oriented" if the a-path
                    # and the b-path are oriented as the x and y axes, i.e.,
                    # if, when we walk around v in counter-clockwise direction,
                    # we encounter a_in,b_in,a_out,b_out.

                    # we can also have that b_in and/or b_out overlaps with
                    # a_in and/or a_out. If we just score the orientation of
                    # b_in and b_out individually, we can deal with this
                    # by just ignoring the overlapping vertex. The "half"
                    # score will be appropriately complemented at one of the
                    # next vertices.

                    a_in = cycles[i][i0-1][0]
                    a_out = cycles[i][(i0+1) % len(cycles[i])][0]
                    b_in = cycles[j][i1-1][0]
                    b_out = cycles[j][(i1+1) % len(cycles[j])][0]

                    # we can get the angles (and hence the rotation order)
                    # by taking the arguments of the differences.

                    a_in_arg = direction(center, a_in)
                    a_out_arg = direction(center, a_out)
                    b_in_arg = direction(center, b_in)
                    b_out_arg = direction(center, b_out)

                    # we make sure to test overlap on the indices, so no rounding
                    # problems occur with that.

                    if (b_in != a_in) and (b_in != a_out):
                        if ((a_in_arg<b_in_arg<a_out_arg) or
                                (b_in_arg<a_out_arg<a_in_arg) or
                                (a_out_arg<a_in_arg<b_in_arg)):
                            intsum += 1
                        elif ((a_out_arg<b_in_arg<a_in_arg) or
                                (b_in_arg<a_in_arg<a_out_arg) or
                                (a_in_arg<a_out_arg<b_in_arg)):
                            intsum -= 1
                        else:
                            raise RuntimeError("impossible edge orientation")
                    if (b_out != a_in) and (b_out != a_out):
                        if ((a_in_arg<b_out_arg<a_out_arg) or
                                (b_out_arg<a_out_arg<a_in_arg) or
                                (a_out_arg<a_in_arg<b_out_arg)):
                            intsum -= 1
                        elif ((a_out_arg<b_out_arg<a_in_arg) or
                                (b_out_arg<a_in_arg<a_out_arg) or
                                (a_in_arg<a_out_arg<b_out_arg)):
                            intsum += 1
                        else:
                            raise RuntimeError("impossible edge orientation")
                assert (intsum % 2) == 0
                intsum = intsum // 2
                intersectionprod[i][j] = intsum
                # Skew Symmetry
                intersectionprod[j][i] = -intsum
        Gmatrix = Matrix(intersectionprod)
        G_normalized, P = Gmatrix.symplectic_form()
        if G_normalized.rank() != 2 * self.genus:
            raise RuntimeError("rank of homology pairing mismatches twice stored genus")
        # Define the cycle sets.
        acycles = [[] for i in range(self.genus)]
        bcycles = [[] for i in range(self.genus)]
        # There are g a and b cycles.
        for i in range(self.genus):
            # Range over the size of the Gram matrix.
            for j in range(cn):
                # Forms the acycles and bcycles. If the entry in the
                # transformation matrix is non-zero, it adds the coefficient at
                # that entry, and the corresponding cycle. (also, forms it
                # into a loop)
                if P[i][j] != 0:
                    acycles[i] += [(P[i][j], [x for x in cycles[j]] + [cycles[j][0]])]
                if P[self.genus + i][j] != 0:
                    bcycles[i] += [(P[self.genus + i][j], [x for x in cycles[j]] + [cycles[j][0]])]
        return acycles + bcycles

    def make_zw_interpolator(self, upstairs_edge):
        r"""
        Given an upstairs edge for which continuation data has been stored,
        return a function that computes `z(t),w(t)` , where `t` in `[0,1]` is a
        parametrization of the edge.

        INPUT:

        - ``upstairs_edge`` -- a pair of integer tuples indicating an edge on
          the upstairs graph of the surface

        OUTPUT:

        A tuple (g, d), where g is the function that computes the interpolation
        along the edge and d is the difference of the z-values of the end and
        start point.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: _ = S.homology_basis()
            sage: g,d = S.make_zw_interpolator([(0,0),(1,0)]);
            sage: all(f(*g(i*0.1)).abs() < 1e-13for i in range(10))
            True
            sage: abs((g(1)[0]-g(0)[0]) - d) < 1e-13
            True
        """
        eindex = tuple(u[0] for u in upstairs_edge)
        i0, i1 = eindex
        z_start = self._vertices[i0]
        z_end = self._vertices[i1]
        currL = self._L[eindex]
        windex = upstairs_edge[0][1]

        def w_interpolate(t):
            if t < 0 or t > 1:
                raise ValueError("t outside path range")
            if t == 0:
                return z_start, currL[0][1][windex]
            elif t == 1:
                return z_end, currL[-1][1][windex]
            while True:
                i = bisect(currL, t)
                t1, w1, epsilon = currL[i]
                w1 = w1[windex]
                t2, w2, _ = currL[i + 1]
                w2 = w2[windex]
                z0 = (1-t)*z_start+t*z_end
                w0 = self._CC(((t2-t)*w1+(t-t1)*w2)/(t2-t1))
                try:
                    desired_result = self._newton_iteration(z0, w0, epsilon)
                except ConvergenceError:
                    pass
                else:
                    return z0, desired_result
                # If we did not succeed, we insert a new point in our interpolation list
                tnew = t
                while True:
                    tnew = (t1 + tnew) / 2
                    znew = (1-tnew)*self._vertices[i0]+tnew*self._vertices[i1]
                    try:
                        neww1 = self._determine_new_w(znew, currL[i][1], epsilon)
                    except ConvergenceError:
                        pass
                    else:
                        # When *no* ConvergenceError is raised, we have succeeded and we can exit
                        break
                # once the loop has succeeded we insert our new value
                t1 = tnew
                self._L[eindex].insert(i + 1, (t1, neww1, epsilon))
        return w_interpolate, (z_end - z_start)

    def simple_vector_line_integral(self, upstairs_edge, differentials):
        r"""
        Perform vectorized integration along a straight path.

        INPUT:

        - ``upstairs_edge`` -- a pair of integer tuples corresponding to an edge
          of the upstairs graph.

        - ``differentials`` -- a list of polynomials; a polynomial `g`
          represents the differential `g(z,w)/(df/dw) dz` where `f(z,w)=0` is
          the equation defining the Riemann surface.

        OUTPUT:

        A complex number, the value of the line integral.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f); S
            Riemann surface defined by polynomial f = -z^4 + w^2 + 1 = 0, with 53 bits of precision

        Since we make use of data from homotopy continuation, we need to compute
        the necessary data::

            sage: M = S.riemann_matrix()
            sage: differentials = S.cohomology_basis()
            sage: S.simple_vector_line_integral([(0,0),(1,0)], differentials) #abs tol 0.00000001
            (1.14590610929717e-16 - 0.352971844594760*I)

        .. NOTE::

            Uses data that ``homology_basis`` initializes.
        """
        w_of_t, Delta_z = self.make_zw_interpolator(upstairs_edge)
        V = VectorSpace(self._CC, len(differentials))

        def integrand(t):
            zt, wt = w_of_t(t)
            dfdwt = self._fastcall_dfdw(zt, wt)
            return V([omega(zt, wt) / dfdwt for omega in differentials])

        return integrate_vector(integrand, self._prec) * Delta_z

    def cohomology_basis(self, option=1):
        r"""
        Compute the cohomology basis of this surface.

        INPUT:

        - ``option`` -- Presently, this routine uses Singular's ``adjointIdeal``
          and passes the ``option`` parameter on. Legal values are 1, 2, 3 ,4,
          where 1 is the default. See the Singular documentation for the
          meaning. The backend for this function may change, and support for
          this parameter may disappear.

        OUTPUT:

        This returns a list of polynomials `g` representing the holomorphic
        differentials `g/(df/dw) dz`, where `f(z,w)=0` is the equation
        specifying the Riemann surface.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)
            sage: S.cohomology_basis()
            [1, w, z]
        """
        if self.genus == 0:
            self._differentials = []
            return self._differentials[0]
        if self._differentials is None:
            # Computes differentials from the adjointIdeal using Singular
            # First we homogenize
            base = self.f.base_ring()
            # It's important we use a degree ordering; see below.
            R = self._R
            k = PolynomialRing(base, names="Z,W,U", order="degrevlex")
            dehom = k.Hom(R)([R.gen(0), R.gen(1), R.one()])
            fnew = self.f(k.gen(0) / k.gen(2), k.gen(1) / k.gen(2)).numerator()

            # We load the relevant functionality into singularlib
            import sage.libs.singular.function_factory
            sage.libs.singular.function_factory.lib("paraplanecurves.lib")
            adjointIdeal = sage.libs.singular.function.singular_function("adjointIdeal")
            libsing_options = sage.libs.singular.option.LibSingularVerboseOptions()

            # We compute the adjoint ideal (note we need to silence "redefine")
            redef_save = libsing_options['redefine']
            try:
                libsing_options['redefine'] = False
                J = adjointIdeal(fnew, option)
            finally:
                libsing_options['redefine'] = redef_save

            # We are interested in the (degree-3) subspace of the adjoint ideal.
            # We compute this by intersecting with (Z,W,U)^(degree-3). Then the
            # lowest degree generators are a basis of the relevant subspace.
            d = fnew.total_degree()
            J2 = k.ideal(J).intersection(k.ideal([k.gen(0), k.gen(1), k.gen(2)])**(d - 3))
            generators = [dehom(c) for c in J2.gens() if c.degree() == d - 3]
            if len(generators) != self.genus:
                raise ValueError("computed regular differentials do not match stored genus")
            self._differentials = generators
        return self._differentials

    def _bounding_data(self, differentials):
        r"""
        Compute the data required to bound a differential on a circle.

        Given a differential, one can bound it on a circular region using its
        derivative and its minimal polynomial (in the coordinate of the base).

        INPUT:

        - ``differentials`` -- The differentials for which to compute the bounding data.

        OUTPUT:

        A tuple containing the polynomial ring the minimal polynomial lies in,
        as well as a list containing the differential, its derivative, the
        minimal polynomial, and ``a0_info`` for each differential. Here
        ``a0_info`` is the leading coefficient and roots of the leading order
        term of the minimal polynomial.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - (z^3-1)
            sage: S = RiemannSurface(f)
            sage: S._bounding_data(S.cohomology_basis())
            (Multivariate Polynomial Ring in z, g over Complex Field with 53 bits of precision,
             [(1/(2*w),
               (-3*z^2*g)/(2*z^3 - 2),
               z^3*g^2 - g^2 - 0.250000000000000,
               (1.00000000000000,
                [1.00000000000000,
                 -0.500000000000000 - 0.866025403784439*I,
                 -0.500000000000000 + 0.866025403784439*I]))])

        """
        # This copies previous work by NB, outputting the zipped list required for a certified line integral. 
        RB = self._R.base_ring()
        P = PolynomialRing(RB, 'Z')
        k = P.fraction_field()
        KP = PolynomialRing(k, 'W') #W->fraction field
        fZW = self.f(P.gen(0), KP.gen(0))
        L = k.extension(fZW, 'Wb')
        dfdw_L = self._dfdw(P.gen(0), L.gen(0))
        integrand_list = [h/self._dfdw for h in differentials]
        # minpoly_univ gives the minimal polynomial for h, in variable x, with coefficients given by polynomials 
        # with coefficients in P (i.e. rational polynomials in Z).
        minpoly_univ = [(h(P.gen(0), L.gen(0))/dfdw_L).minpoly().numerator()
                        for h in differentials]
        RBzg = PolynomialRing(RB, ['z', 'g'])
        # The following line changes the variables in these minimal polynomials as Z -> z, x -> G, then evaluates
        # at G = QQzg.gens(1) ( = g )
        RBzgG = PolynomialRing(RBzg, 'G')
        minpoly_list = [RBzgG([c(RBzg.gen(0)) for c in list(h)])(RBzg.gen(1))
                        for h in minpoly_univ]
        # h(z,g)=0 --> dg/dz = - dhdz/dhdg
        dgdz_list = [-h.derivative(RBzg.gen(0))/h.derivative(RBzg.gen(1))
                     for h in minpoly_list]
        
        CCzg = PolynomialRing(self._CC, ['z','g'])
        CCminpoly_list = [CCzg(h) for h in minpoly_list]
        
        a0_list = [P(h.leading_coefficient()) for h in minpoly_univ]
        #CCa0_list = [self._CCz(a0) for a0 in a0_list]
        # Note that because the field over which the Riemann surface is defined
        # is embedded into CC, it has characteristic 0, and so we know the 
        # irreducible factors are all separable, i.e. the roots have multiplicity
        # one. 
        a0_info = [(self._CC(a0.leading_coefficient()),
                    flatten([self._CCz(F).roots(multiplicities=False)*m 
                             for F, m in a0.factor()]))
                   for a0 in a0_list]

        #a0_list = [self._CCz(h.leading_coefficient()) for h in minpoly_univ]
        #a0_info = [(a0.leading_coefficient(), flatten([[r]*n for r, n in a0.roots()])) for a0 in a0_list]
        # Note that we are assuming here that the root finder is working correctly, which will need
        # some thought.
        return CCzg, list(zip(integrand_list, dgdz_list, CCminpoly_list, a0_info))

    def rigorous_line_integral(self, upstairs_edge, differentials, bounding_data):
        r"""
        Perform vectorized integration along a straight path.

        Using the error bounds for Gauss-Legendre integration found in [Neu2018]_
        and the method for bounding an algebraic integrand on a circular domain
        developed by Gao, calculates (semi-)rigorously the integral of a list of 
        differentials along an edge of the upstairs graph

        INPUT:

        - ``upstairs_edge`` -- a pair of integer tuples corresponding to an edge
          of the upstairs graph.

        - ``differentials`` -- a list of polynomials; a polynomial `g`
          represents the differential `g(z,w)/(df/dw) dz` where `f(z,w)=0` is
          the equation defining the Riemann surface.

        - ``bounding_data`` -- tuple containing the data required for bounding
          the integrands. This should be in the form of the output from 
          ``_bounding_data``.

        OUTPUT:

        A complex number, the value of the line integral.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f); S
            Riemann surface defined by polynomial f = -z^4 + w^2 + 1 = 0, with 53 bits of precision

        Since we make use of data from homotopy continuation, we need to compute
        the necessary data::

            sage: _ = S.homology_basis()
            sage: differentials = S.cohomology_basis()
            sage: bounding_data = S._bounding_data(differentials)
            sage: S.rigorous_line_integral([(0,0), (1,0)], differentials, 
                                                bounding_data)
            (1.80277751848459e-16 - 0.352971844594760*I)

        .. NOTE::

            Uses data that ``homology_basis`` initializes. In practice is it is 
            more efficient to set ``differentials`` to a fast-callable version 
            of differentials to speed up execution. 

        .. TODO::

            Note that bounding_data contains the information of the integrands,
            so we are being inefficient by also asking for differentials as an 
            input. It would be simpler and ensure consistency if the 
            differentials were read from the bounding data.

            Moreover, the current implementation bounds along a line by 
            splitting it up into segments, each of which can be covered entirely
            by a single circle, and then placing inside that the ellipse 
            requierd to bound as per [Neu2018]_. This is reliably more efficient
            than the heuristic method, especially in poorly-conditioned cases 
            where discriminant points are close together around the edges, but
            in the case where the branch locus is well separated, it can require
            slightly more nodes than necessary. One may want to include a method
            here to transition in this regime to an algorithm that covers the 
            entire line with one ellipse, then bounds along that ellipse with 
            multiple circles. 
        """
        # Note that this, in it's current formalism, makes no check that bounding data at all corresponds to 
        # the differentials given. The onus is then on the design of other funcions which use it.
    
        # CCzg is required to be known as we need to know the ring which the minpolys lie in. 
        CCzg, bounding_data_list = bounding_data
        
        i0, _ = upstairs_edge[0]
        i1, _ = upstairs_edge[1]
        z0 = self._vertices[i0]
        z1 = self._vertices[i1]
        zwt, z1_minus_z0 = self.make_zw_interpolator(upstairs_edge)
        
        # list of (centre, radius) pairs that still need to be processed
        ball_stack = [(self._RR(0.5), self._RR(0.5))]
        alpha = self._RR(0.912) 
        # alpha set manually for scaling purposes. Basic benchmarking shows 
        # that ~0.9 is a sensible value. 
        E_global = self._RR(2)**(-self._prec+3)
        K = 2
        # The parameter K could be tuned, but basic benchmarking seems to show
        # that 2 is a sensible choice

        # Output will iteratively store the output of the integral. 
        V = VectorSpace(self._CC, len(differentials))
        output = V(0)

        while ball_stack:
            ct, rt = ball_stack.pop()
            cz = (1-ct)*z0+ct*z1 # This is the central z-value of our ball.
            distances = [(cz-b).abs() for b in self.branch_locus] # Distance to the discriminant points
            rho_z = min(distances)
            rho_t = rho_z/(z1-z0).abs()
            if rho_t > rt:
                rho_t = alpha*rho_t+(1-alpha)*rt # sqrt(rho_t*rt) could also work
                rho_z = rho_t*(z1-z0).abs()
                delta_z = (alpha*rho_t+(1-alpha)*rt)*(z1-z0).abs()
                expr = rho_t/rt+((rho_t/rt)**2-1).sqrt(self._prec) # Note this is really exp(arcosh(rho_t/rt))
                N = 3
                cw = zwt(ct)[1]
                for g, dgdz, minpoly,(a0lc,a0roots) in bounding_data_list:
                    z_1 = a0lc.abs()*prod((cz-r).abs()-rho_z for r in a0roots)
                    n = minpoly.degree(CCzg.gen(1))
                    # Note the structure of the code is currently s.t 'z' has to be the variable in
                    # the minpolys.
                    ai_new = [(minpoly.coefficient({CCzg.gen(1):i}))(z=cz+self._CCz.gen(0)) for i
                                in range(n)]
                    ai_pos = [ self._RRz([c.abs() for c in h.list()]) for h in ai_new]
                    m = [a(rho_z)/z_1 for a in ai_pos]
                    l = len(m)
                    M_tilde = 2*max((m[i].abs())**(1/self._RR(l-i)) for i in range(l))
                    cg = g(cz,cw)
                    cdgdz = dgdz(cz,cg)
                    Delta = delta_z*cdgdz.abs()+ (delta_z**2)*M_tilde/(rho_z*(rho_z-delta_z))
                    M = Delta #+ abs(cg)
                    N_required = ((64*M/(15*(1-1/expr)*E_global)).log()/(2*expr.log())).ceil()
                    N = max(N,N_required)

                N = (K*(N.sqrt(self._prec)/K).ceil())**2
                # Rounding is sensible as it allows the cache of nodes in 
                # sage.numerical.gauss_legendre to be used.
                # Quadratic rounding can be shown to be a sensible choice through the 
                # basic argument that nodes is quadratic in N 
                
                ct_minus_rt = ct-rt
                two_rt = 2*rt
                def integrand(t):
                    zt, wt = zwt(ct_minus_rt+t*two_rt)
                    dfdwt = self._fastcall_dfdw(zt, wt)
                    return V([h(zt,wt)/dfdwt for h in differentials])

                output += two_rt*integrate_vector_N(integrand,self._prec,N)
            else:
                ball_stack.append((ct-rt/2, rt/2))
                ball_stack.append((ct+rt/2, rt/2))

        return output*z1_minus_z0

    def matrix_of_integral_values(self, differentials, integration_method="heuristic"):
        r"""
        Compute the path integrals of the given differentials along the homology
        basis.

        The returned answer has a row for each differential. If the Riemann
        surface is given by the equation `f(z,w)=0`, then the differentials are
        encoded by polynomials g, signifying the differential `g(z,w)/(df/dw)
        dz`.

        INPUT:

        - ``differentials`` -- a list of polynomials.

        - ``integration_method`` -- (default: ``'heuristic'``). String specifying
          the integration method to use. The options are ``'heuristic'`` and 
          ``'rigorous'``.

        OUTPUT:

        A matrix, one row per differential, containing the values of the path
        integrals along the homology basis of the Riemann surface.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(x^3 + y^3 + 1)
            sage: B = S.cohomology_basis()
            sage: m = S.matrix_of_integral_values(B)
            sage: parent(m)
            Full MatrixSpace of 1 by 2 dense matrices over Complex Field with 53 bits of precision
            sage: (m[0,0]/m[0,1]).algdep(3).degree() #curve is CM, so the period is quadratic
            2

        """
        cycles = self.homology_basis()

        def normalize_pairs(L):
            r"""
            Return a list of edges encoded by the path in L.
            The edges are normalized to be in the direction in which
            the homotopy continuation should have been computed along them.
            """
            R = []
            for i in range(len(L) - 1):
                if L[i][0] < L[i + 1][0]:
                    R.append((L[i], L[i + 1]))
                else:
                    R.append((L[i + 1], L[i]))
            return R
        occurring_edges = set()
        occurring_edges.update(*[normalize_pairs(p[1]) for h in cycles
                                 for p in h])

        fcd = [fast_callable(omega, domain=self._CC) for omega in differentials]

        if integration_method=="heuristic":
            line_int = lambda edge: self.simple_vector_line_integral(edge,fcd)
        elif integration_method=="rigorous":
            bd = self._bounding_data(differentials)
            line_int = lambda edge: self.rigorous_line_integral(edge,fcd,bd)
        else:
            raise ValueError("Invalid integration method")

        #integral_dict = {upstairs_edge:
        #                 self.simple_vector_line_integral(upstairs_edge,
        #                                                  differentials)
        #                 for upstairs_edge in occurring_edges}
        integral_dict = dict()
        for upstairs_edge in occurring_edges:
            integral_dict[upstairs_edge] = line_int(upstairs_edge)
        
        rows = []
        for cycle in cycles:
            V = VectorSpace(self._CC, len(differentials)).zero()
            for multiplicity, loop in cycle:
                for i in range(len(loop) - 1):
                    if loop[i][0] < loop[i + 1][0]:
                        direction = 1
                        upstairs_edge = (loop[i], loop[i + 1])
                    else:
                        direction = -1
                        upstairs_edge = (loop[i + 1], loop[i])
                    V += (multiplicity * direction) * integral_dict[upstairs_edge]
            rows.append(V)
        return Matrix(rows).transpose()

    @cached_method
    def period_matrix(self):
        r"""
        Compute the period matrix of the surface.

        OUTPUT:

        A matrix of complex values.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f, prec=30)
            sage: M = S.period_matrix()

        The results are highly arbitrary, so it is hard to check if the result
        produced is correct. The closely related ``riemann_matrix`` is somewhat
        easier to test.::

            sage: parent(M)
            Full MatrixSpace of 3 by 6 dense matrices over Complex Field with 30 bits of precision
            sage: M.rank()
            3

        One can check that the two methods give similar answers::
        
            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: f = y^2 - x^3 + 1
            sage: S = RiemannSurface(f, integration_method="rigorous")
            sage: T = RiemannSurface(f)
            sage: S.riemann_matrix()
            [-0.500000000000000 + 0.866025403784439*I]
            sage: T.riemann_matrix()
            [-0.500000000000000 + 0.866025403784438*I]

        and that in this example the rigorous method accurate to more digits::

            sage: (sqrt(3)/2).n(53)
            0.866025403784439
        """
        #differentials = [fast_callable(omega, domain=self._CC)
        #                 for omega in self.cohomology_basis()]
        differentials = self.cohomology_basis()
        return self.matrix_of_integral_values(differentials, self._integration_method)

    def riemann_matrix(self):
        r"""
        Compute the Riemann matrix.

        OUTPUT:

        A matrix of complex values.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f, prec=60)
            sage: M = S.riemann_matrix()

        The Klein quartic has a Riemann matrix with values is a quadratic
        field::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2-x+2)
            sage: all(len(m.algdep(6).roots(K)) > 0 for m in M.list())
            True
        """
        PeriodMatrix = self.period_matrix()
        Am = PeriodMatrix[0:self.genus,0:self.genus]
        RM = numerical_inverse(Am)*PeriodMatrix[0:self.genus,self.genus:2*self.genus]
        return RM

    def plot_paths(self):
        r"""
        Make a graphical representation of the integration paths.

        This returns a two dimensional plot containing the branch points (in red) and
        the integration paths (obtained from the Voronoi cells of the branch
        points). The integration paths are plotted by plotting the points that
        have been computed for homotopy continuation, so the density gives an
        indication of where numerically sensitive features occur.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - x^3 - x)
            sage: S.plot_paths()
            Graphics object consisting of 2 graphics primitives
        """
        from sage.plot.point import point2d
        P = []

        # trigger the computation of the homology basis, so that self._L is present
        self.homology_basis()

        for e in self._L.keys():
            z0 = self._vertices[e[0]]
            z1 = self._vertices[e[1]]

            def path(t):
                return (1 - t) * z0 + t * z1
            T = self._L[e]
            P += [path(t[0]) for t in T]
        return point2d(P, size=1) + point2d(self.branch_locus, color="red")

    def plot_paths3d(self, thickness=0.01):
        r"""
        Return the homology basis as a graph in 3-space.

        The homology basis of the surface is constructed by taking the Voronoi
        cells around the branch points and taking the inverse image of the edges
        on the Riemann surface. If the surface is given by the equation
        `f(z,w)`, the returned object gives the image of this graph in 3-space
        with coordinates `\left(\operatorname{Re}(z), \operatorname{Im}(z),
        \operatorname{Im}(w)\right)`.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2-x^3-x)
            sage: S.plot_paths3d()
            Graphics3d Object
        """
        from sage.plot.graphics import Graphics
        from sage.plot.plot3d.shapes2 import point3d, line3d
        P = Graphics()

        # trigger the computation of the homology basis, so that
        # self._L is present
        self.homology_basis()

        for e in self._L.keys():
            z0 = self._vertices[e[0]]
            z1 = self._vertices[e[1]]

            def path(t):
                z = (1-t)*z0+t*z1
                return (z.real_part(),z.imag_part())
            T = self._L[e]
            color = "blue"
            for i in range(self.degree):
                P += line3d([path(t[0])+(t[1][i].imag_part(),) for t in T],color=color,thickness=thickness)
            for z,ws in zip(self._vertices,self._wvalues):
                for w in ws:
                    P += point3d([z.real_part(),z.imag_part(),w.imag_part()],color="purple", size=20)
        return P

    def endomorphism_basis(self, b=None, r=None):
        r"""
        Numerically compute a `\ZZ`-basis for the endomorphism ring.

        Let `\left(I | M \right)` be the normalized period matrix (`M` is the
        `g\times g` :meth:`riemann_matrix`). We consider the system of matrix
        equations `MA + C = (MB + D)M` where `A, B, C, D` are `g\times g`
        integer matrices.  We determine small integer (near) solutions using LLL
        reductions.  These solutions are returned as `2g \times 2g` integer
        matrices obtained by stacking `\left(D | B\right)` on top of `\left(C |
        A\right)`.

        INPUT:

        - ``b`` -- integer (default provided). The equation coefficients are
          scaled by `2^b` before rounding to integers.

        - ``r`` -- integer (default: ``b/4``). Solutions that have all
          coefficients smaller than `2^r` in absolute value are reported as
          actual solutions.

        OUTPUT:

        A list of `2g \times 2g` integer matrices that, for large enough ``r``
        and ``b-r``, generate the endomorphism ring.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(x^3 + y^3 + 1)
            sage: B = S.endomorphism_basis(); B #random
            [
            [1 0]  [ 0 -1]
            [0 1], [ 1  1]
            ]
            sage: sorted([b.minpoly().disc() for b in B])
            [-3, 1]

        """
        M = self.riemann_matrix()
        return integer_matrix_relations(M,M,b,r)

    def homomorphism_basis(self, other, b=None, r=None):
        r"""
        Numerically compute a `\ZZ`-basis for module of homomorphisms to a given
        complex torus.

        Given another complex torus (given as the analytic Jacobian of a Riemann
        surface), numerically compute a basis for the homomorphism module. The
        answer is returned as a list of 2g x 2g integer matrices T=(D, B; C, A)
        such that if the columns of (I|M1) generate the lattice defining the
        Jacobian of the Riemann surface and the columns of (I|M2) do this for
        the codomain, then approximately we have (I|M2)T=(D+M2C)(I|M1), i.e., up
        to a choice of basis for `\CC^g` as a complex vector space, we we
        realize (I|M1) as a sublattice of (I|M2).

        INPUT:

        - ``b`` -- integer (default provided). The equation coefficients are
          scaled by `2^b` before rounding to integers.

        - ``r`` -- integer (default: ``b/4``). Solutions that have all
          coefficients smaller than `2^r` in absolute value are reported as
          actual solutions.

        OUTPUT:

        A list of `2g \times 2g` integer matrices that, for large enough ``r``
        and ``b-r``, generate the homomorphism module.

        EXAMPLES::

            sage: S1 = EllipticCurve("11a1").riemann_surface()
            sage: S2 = EllipticCurve("11a3").riemann_surface()
            sage: [m.det() for m in S1.homomorphism_basis(S2)]
            [5]
        """
        M1 = self.riemann_matrix()
        M2 = other.riemann_matrix()
        return integer_matrix_relations(M2,M1,b,r)

    def tangent_representation_numerical(self, Rs, other=None):
        r"""
        Compute the numerical tangent representations corresponding to the
        homology representations in ``Rs``.

        The representations on homology ``Rs`` have to be given with respect to
        the symplectic homology basis of the Jacobian of ``self`` and ``other``.
        Such matrices can for example be obtained via
        :meth:`endomorphism_basis`.

        Let `P` and `Q` be the period matrices of ``self`` and ``other``. Then
        for a homology representation `R`, the corresponding tangential
        representation `T` satisfies `T P = Q R`.

        INPUT:

        - ``Rs`` -- a set of matrices on homology to be converted to their
          tangent representations.

        - ``other`` (default: ``self``) -- the codomain, another Riemann
          surface.

        OUTPUT:

        The numerical tangent representations of the matrices in ``Rs``.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: A.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - (x^6 + 2*x^4 + 4*x^2 + 8), prec = 100)
            sage: P = S.period_matrix()
            sage: Rs = S.endomorphism_basis()
            sage: Ts = S.tangent_representation_numerical(Rs)
            sage: all(((T*P - P*R).norm() < 2^(-80)) for [T, R] in zip(Ts, Rs))
            True
        """
        if not other:
            other = self
        P = self.period_matrix()
        CCP = P.base_ring()
        g = self.genus
        Q = other.period_matrix()
        Ptsubinv = numerical_inverse((P.transpose())[list(range(g))])
        Ts = []
        for R in Rs:
            QRtsub = ((Q * R).transpose())[list(range(g))]
            Tt = Ptsubinv * QRtsub
            T = Tt.transpose().change_ring(CCP)
            Ts.append(T)
        return Ts

    def tangent_representation_algebraic(self, Rs, other=None, epscomp=None):
        r"""
        Compute the algebraic tangent representations corresponding to the
        homology representations in ``Rs``.

        The representations on homology ``Rs`` have to be given with respect to
        the symplectic homology basis of the Jacobian of ``self`` and ``other``.
        Such matrices can for example be obtained via
        :meth:`endomorphism_basis`.

        Let `P` and `Q` be the period matrices of ``self`` and ``other``. Then
        for a homology representation `R`, the corresponding tangential
        representation `T` satisfies `T P = Q R`.

        INPUT:

        - ``Rs`` -- a set of matrices on homology to be converted to their
          tangent representations.

        - ``other`` (default: ``self``) -- the codomain, another Riemann
          surface.

        - ``epscomp`` -- real number (default: ``2^(-prec + 30)``). Used to
          determine whether a complex number is close enough to a root of a
          polynomial.

        OUTPUT:

        The algebraic tangent representations of the matrices in ``Rs``.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: A.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - (x^6 + 2*x^4 + 4*x^2 + 8), prec = 100)
            sage: Rs = S.endomorphism_basis()
            sage: Ts = S.tangent_representation_algebraic(Rs)
            sage: Ts[0].base_ring().maximal_order().discriminant() == 8
            True
        """
        if not epscomp:
            epscomp = 2**(-self._prec + 30)
        QQalg = QQ.algebraic_closure()

        def polynomialize_element(alpha):
            d = 1
            while True:
                d += 1
                dep = algdep(alpha, d, height_bound=10**d)
                if dep and dep(alpha) < epscomp:
                    return dep

        def algebraize_element(alpha):
            alphaPol = polynomialize_element(alpha)
            CC = alpha.parent()
            for tup in alphaPol.roots(QQalg):
                rt = tup[0]
                if (alpha - CC(rt)).abs() < epscomp:
                    return rt
            raise AssertionError('No close root found while algebraizing')

        def algebraize_matrices(Ts):
            nr = Ts[0].nrows()
            nc = Ts[0].ncols()
            TsAlg = [T.apply_map(algebraize_element) for T in Ts]
            elts = [x for TAl in TsAlg for x in TAl.list()]
            eltsAlg = number_field_elements_from_algebraics(elts)[1]
            L = eltsAlg[0].parent()
            TsAlgL = []
            for i in range(len(Ts)):
                TAlgL = [eltsAlg[j] for j in range(i*nr*nc, (i + 1)*nr*nc)]
                TsAlgL.append(Matrix(L, nr, nc, TAlgL))
            return TsAlgL

        Ts = self.tangent_representation_numerical(Rs, other=other)
        return algebraize_matrices(Ts)

    def rosati_involution(self, R):
        r"""
        Compute the Rosati involution of an endomorphism.

        The endomorphism in question should be given by its homology
        representation with respect to the symplectic basis of the Jacobian.

        INPUT:

        - ``R`` -- integral matrix.

        OUTPUT:

        The result of applying the Rosati involution to ``R``.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: A.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - (x^6 + 2*x^4 + 4*x^2 + 8), prec = 100)
            sage: Rs = S.endomorphism_basis()
            sage: S.rosati_involution(S.rosati_involution(Rs[1])) == Rs[1]
            True
        """
        def standard_symplectic_matrix(n):
            one = Matrix.identity(n)
            zero = Matrix.zero(n)
            return Matrix.block([[zero, -one], [one, zero]])
        g = self.genus
        if not(R.nrows() == 2 * g == R.ncols()):
            raise AssertionError("Matrix is not the homology representation of an endomorphism")
        J = standard_symplectic_matrix(g)
        return -J * R.transpose() * J

    def symplectic_isomorphisms(self, other=None, hom_basis=None, b=None, r=None):
        r"""
        Numerically compute symplectic isomorphisms.

        INPUT:

        - ``other`` (default: ``self``) -- the codomain, another Riemann
          surface.

        - ``hom_basis`` (default: ``None``) -- a `\ZZ`-basis of the
          homomorphisms from ``self`` to ``other``, as obtained from
          :meth:`homomorphism_basis`. If you have already calculated this
          basis, it saves time to pass it via this keyword argument. Otherwise
          the method will calculate it.

        - ``b`` -- integer (default provided): as for
          :meth:`homomorphism_basis`, and used in its invocation if
          (re)calculating said basis.

        - ``r`` -- integer (default: ``b/4``).  as for
          :meth:`homomorphism_basis`, and used in its invocation if
          (re)calculating said basis.

        OUTPUT:

        This returns the combinations of the elements of
        :meth:`homomorphism_basis` that correspond to symplectic
        isomorphisms between the Jacobians of ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: f = y^2 - (x^6 + 2*x^4 + 4*x^2 + 8)
            sage: X = RiemannSurface(f, prec=100)
            sage: P = X.period_matrix()
            sage: g = y^2 - (x^6 + x^4 + x^2 + 1)
            sage: Y = RiemannSurface(g, prec=100)
            sage: Q = Y.period_matrix()
            sage: Rs = X.symplectic_isomorphisms(Y)
            sage: Ts = X.tangent_representation_numerical(Rs, other = Y)
            sage: test1 = all(((T*P - Q*R).norm() < 2^(-80)) for [T, R] in zip(Ts, Rs))
            sage: test2 = all(det(R) == 1 for R in Rs)
            sage: test1 and test2
            True
        """
        if not other:
            other = self
        if hom_basis:
            Rs = hom_basis
        else:
            Rs = self.homomorphism_basis(other=other, b=b, r=r)
        r = len(Rs)
        g = self.genus
        A = PolynomialRing(QQ, r, 'x')
        gensA = A.gens()
        # Use that the trace is positive definite; we could also put this as an
        # extra condition when determining the endomorphism basis to speed up
        # that calculation slightly
        R = sum(gensA[i] * Rs[i].change_ring(A) for i in range(r))
        tr = (R * self.rosati_involution(R)).trace()
        # Condition tr = 2 g creates ellipsoid
        M = Matrix(ZZ, r, r, [tr.derivative(gen1).derivative(gen2)
                              for gen1 in gensA for gen2 in gensA])
        vs = M.__pari__().qfminim(4*g)[2].sage().transpose()
        vs = [v for v in vs if v * M * v == 4*g]
        vs += [-v for v in vs]
        RsIso = []
        for v in vs:
            R = sum(v[i] * Rs[i] for i in range(r))
            if R * self.rosati_involution(R) == 1:
                RsIso.append(R)
        return RsIso

    def symplectic_automorphism_group(self, endo_basis=None, b=None, r=None):
        r"""
        Numerically compute the symplectic automorphism group as a permutation
        group.

        INPUT:

        - ``endo_basis`` (default: ``None``) -- a `\ZZ`-basis of the
          endomorphisms of ``self``, as obtained from
          :meth:`endomorphism_basis`. If you have already calculated this
          basis, it saves time to pass it via this keyword argument. Otherwise
          the method will calculate it.

        - ``b`` -- integer (default provided): as for
          :meth:`homomorphism_basis`, and used in its invocation if
          (re)calculating said basis.

        - ``r`` -- integer (default: ``b/4``).  as for
          :meth:`homomorphism_basis`, and used in its invocation if
          (re)calculating said basis.

        OUTPUT:

        The symplectic automorphism group of the Jacobian of the Riemann
        surface. The automorphism group of the Riemann surface itself can be
        recovered from this; if the curve is hyperelliptic, then it is
        identical, and if not, then one divides out by the central element
        corresponding to multiplication by -1.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: A.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2 - (x^6 + 2*x^4 + 4*x^2 + 8), prec = 100)
            sage: G = S.symplectic_automorphism_group()
            sage: G.as_permutation_group().is_isomorphic(DihedralGroup(4))
            True
        """
        RsAut = self.symplectic_isomorphisms(hom_basis=endo_basis, b=b, r=r)
        return MatrixGroup(RsAut)

    def __add__(self, other):
        r"""
        Return the disjoint union of the Riemann surface and the other argument.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y> = QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S1+S1
            Riemann surface sum with period lattice of rank 4
        """
        return RiemannSurfaceSum([self, other])


def integer_matrix_relations(M1, M2, b=None, r=None):
    r"""
    Determine integer relations between complex matrices.

    Given two square matrices with complex entries of size g, h respectively,
    numerically determine an (approximate) ZZ-basis for the 2g x 2h matrices
    with integer entries of the shape (D, B; C, A) such that B+M1*A=(D+M1*C)*M2.
    By considering real and imaginary parts separately we obtain `2gh` equations
    with real coefficients in `4gh` variables. We scale the coefficients by a
    constant `2^b` and round them to integers, in order to obtain an integer
    system of equations. Standard application of LLL allows us to determine near
    solutions.

    The user can specify the parameter `b`, but by default the system will
    choose a `b` based on the size of the coefficients and the precision with
    which they are given.

    INPUT:

    - ``M1`` -- square complex valued matrix

    - ``M2`` -- square complex valued matrix of same size as M1

    - ``b`` -- integer (default provided). The equation coefficients are scaled
      by `2^b` before rounding to integers.

    - ``r`` -- integer (default: ``b/4``). The vectors found by LLL that satisfy
      the scaled equations to within `2^r` are reported as solutions.

    OUTPUT:

    A list of 2g x 2h integer matrices that, for large enough `r`, `b-r`,
    generate the ZZ-module of relevant transformations.

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import integer_matrix_relations
        sage: M1=M2=matrix(CC,2,2,[sqrt(d) for d in [2,-3,-3,-6]])
        sage: T=integer_matrix_relations(M1,M2)
        sage: id=parent(M1)(1)
        sage: M1t=[id.augment(M1) * t for t in T]
        sage: [((m[:,:2]^(-1)*m)[:,2:]-M2).norm() < 1e-13 for m in M1t]
        [True, True]
    """
    if not(M1.is_square() and M2.is_square()):
        raise ValueError("matrices need to be square")
    prec = min(M1.base_ring().precision(),M2.base_ring().precision())
    H = max(max(abs(m.real_part()) for m in M1.list() + M2.list()),
            max(abs(m.imag_part()) for m in M1.list() + M2.list()))
    if b is None:
        b = prec-5-H.log2().floor()
    if r is None:
        r = b//4
    S = 2**b
    if H*S > 2**(prec-4):
        raise ValueError("insufficient precision for b=%s" % b)
    g1 = M1.ncols()
    g2 = M2.ncols()
    CC = M1.base_ring() if (M1.base_ring().precision() <= M2.base_ring().precision()) else M2.base_ring()
    V = ["%s%s" % (n, i) for n in ["a","b","c","d"] for i in range(1,1+g1*g2)]
    R = PolynomialRing(CC, V)
    vars = R.gens()
    A = Matrix(R, g1, g2, vars[:g1*g2])
    B = Matrix(R, g1, g2, vars[g1*g2:2*g1*g2])
    C = Matrix(R, g1, g2, vars[2*g1*g2:3*g1*g2])
    D = Matrix(R, g1, g2, vars[3*g1*g2:4*g1*g2])
    W = ((M1*A+B) - (M1*C+D)*M2).list()
    vars = R.gens()
    mt = Matrix(ZZ,[[1 if i==j else 0 for j in range(4*g1*g2)] +
      [(S*w.monomial_coefficient(vars[i]).real_part()).round() for w in W] +
      [(S*w.monomial_coefficient(vars[i]).imag_part()).round() for w in W] for i in range(len(vars))])
    # we compute an LLL-reduced basis of this lattice:
    mtL = mt.LLL()

    def vectomat(v):
        A = Matrix(g1,g2,v[:g1*g2].list())
        B = Matrix(g1,g2,v[g1*g2:2*g1*g2].list())
        C = Matrix(g1,g2,v[2*g1*g2:3*g1*g2].list())
        D = Matrix(g1,g2,v[3*g1*g2:4*g1*g2].list())
        return D.augment(B).stack(C.augment(A))
    c = 2**r
    return [vectomat(v) for v in mtL if all(a.abs() <= c for a in v[g1*g2:])]


class RiemannSurfaceSum(RiemannSurface):
    r"""
    Represent the disjoint union of finitely many Riemann surfaces.

    Rudimentary class to represent disjoint unions of Riemann surfaces. Exists
    mainly (and this is the only functionality actually implemented) to
    represents direct products of the complex tori that arise as analytic
    Jacobians of Riemann surfaces.

    INPUT:

    - L -- list of RiemannSurface objects

    EXAMPLES::

        sage: _.<x> = QQ[]
        sage: SC = HyperellipticCurve(x^6-2*x^4+3*x^2-7).riemann_surface(prec=60)
        sage: S1 = HyperellipticCurve(x^3-2*x^2+3*x-7).riemann_surface(prec=60)
        sage: S2 = HyperellipticCurve(1-2*x+3*x^2-7*x^3).riemann_surface(prec=60)
        sage: len(SC.homomorphism_basis(S1+S2))
        2
    """
    def __init__(self, L):
        r"""
        TESTS::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y> = QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S2 = RiemannSurface(y^2-x^3-x-5)
            sage: S = RiemannSurfaceSum([S1,S2])
            sage: S.riemann_matrix() == S1.riemann_matrix().block_sum(S2.riemann_matrix())
            True
        """
        if not all(isinstance(l, RiemannSurface) for l in L):
            raise ValueError("summands must be RiemannSurface objects")
        prec = min(l._prec for l in L)
        self._prec = prec
        self.genus = sum(s.genus for s in L)
        it = iter(L)
        s = next(it)
        g = s.genus
        PM = s.period_matrix()
        PM1 = PM[:g, :g]
        PM2 = PM[:g, g:2*g]
        tau = s.riemann_matrix()
        for s in it:
            g = s.genus
            PM = s.period_matrix()
            PM1 = PM1.block_sum(PM[:g, :g])
            PM2 = PM2.block_sum(PM[:g, g:2*g])
            tau = tau.block_sum(s.riemann_matrix())
        self.PM = block_matrix([[PM1, PM2]], subdivide=False)
        self.tau = tau

    def period_matrix(self):
        r"""
        Return the period matrix of the surface.

        This is just the diagonal block matrix constructed from the period
        matrices of the constituents.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y> = QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S2 = RiemannSurface(y^2-x^3-x-5)
            sage: S = RiemannSurfaceSum([S1,S2])
            sage: S1S2 = S1.period_matrix().block_sum(S2.period_matrix())
            sage: S.period_matrix() == S1S2[[0,1],[0,2,1,3]]
            True
        """
        return self.PM

    def riemann_matrix(self):
        r"""
        Return the normalized period matrix of the surface.

        This is just the diagonal block matrix constructed from the Riemann
        matrices of the constituents.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y> = QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S2 = RiemannSurface(y^2-x^3-x-5)
            sage: S = RiemannSurfaceSum([S1,S2])
            sage: S.riemann_matrix() == S1.riemann_matrix().block_sum(S2.riemann_matrix())
            True
        """
        return self.tau

    def __repr__(self) -> str:
        r"""
        Return string describing Riemann surface sum.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y> = QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S2 = RiemannSurface(y^2-x^3-x-5)
            sage: RiemannSurfaceSum([S1,S2])
            Riemann surface sum with period lattice of rank 4
        """
        return "Riemann surface sum with period lattice of rank " + str(2 * self.genus)

    def __add__(self, other):
        r"""
        Return the disjoint union of the Riemann surface and the other argument.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, RiemannSurfaceSum
            sage: R.<x,y> = QQ[]
            sage: S1 = RiemannSurface(y^2-x^3-x-1)
            sage: S1+S1+S1
            Riemann surface sum with period lattice of rank 6
        """
        return RiemannSurfaceSum([self, other])
