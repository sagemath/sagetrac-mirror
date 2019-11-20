r"""
Base class for polyhedra over `\QQ`
"""
from __future__ import absolute_import

from sage.rings.all import QQ
from sage.misc.all import prod
from .base import Polyhedron_base


class Polyhedron_QQ(Polyhedron_base):
    r"""
    Base class for Polyhedra over `\QQ`

    TESTS::

        sage: p = Polyhedron([(0,0)], base_ring=QQ);  p
        A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex
        sage: TestSuite(p).run()
    """
    def _is_zero(self, x):
        """
        Test whether ``x`` is zero.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=QQ)
            sage: p._is_zero(0)
            True
            sage: p._is_zero(1/100000)
            False
        """
        return x==0

    def _is_nonneg(self, x):
        """
        Test whether ``x`` is nonnegative.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=QQ)
            sage: p._is_nonneg(1)
            True
            sage: p._is_nonneg(-1/100000)
            False
        """
        return x>=0

    def _is_positive(self, x):
        """
        Test whether ``x`` is positive.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=QQ)
            sage: p._is_positive(1)
            True
            sage: p._is_positive(0)
            False
        """
        return x>0

    _base_ring = QQ

    def integral_points_count(self, verbose=False, use_Hrepresentation=False,
                              explicit_enumeration_threshold=1000, preprocess=True, **kwds):
        r"""
        Return the number of integral points in the polyhedron.

        This method uses the optional package ``latte_int``
        if an estimate for lattice points based on bounding boxes
        exceeds ``explicit_enumeration_threshold``.

        INPUT:

        - ``verbose`` (boolean; ``False`` by default) -- whether to display
          verbose output.

        - ``use_Hrepresentation`` - (boolean; ``False`` by default) -- whether
          to send the H or V representation to LattE

        - ``preprocess`` - (boolean; ``True`` by default) -- whether, if the integral hull
          is known to lie in a coordinate hyperplane, to tighten bounds to reduce dimension

        .. SEEALSO::

            :mod:`~sage.interfaces.latte` the interface to LattE interfaces

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.integral_points_count()
            27
            sage: P.integral_points_count(explicit_enumeration_threshold=0) # optional - latte_int
            27

        We enlarge the polyhedron to force the use of the generating function methods
        implemented in LattE integrale, rather than explicit enumeration.

            sage: (1000000000*P).integral_points_count(verbose=True) # optional - latte_int
            This is LattE integrale...
            ...
            Total time:...
            8000000012000000006000000001

        We shrink the polyhedron a little bit::

            sage: Q = P*(8/9)
            sage: Q.integral_points_count()
            1
            sage: Q.integral_points_count(explicit_enumeration_threshold=0) # optional - latte_int
            1

        Unbounded polyhedra (with or without lattice points) are not supported::

            sage: P = Polyhedron(vertices=[[1/2, 1/3]], rays=[[1, 1]])
            sage: P.integral_points_count()
            Traceback (most recent call last):
            ...
            NotImplementedError: ...
            sage: P = Polyhedron(vertices=[[1, 1]], rays=[[1, 1]])
            sage: P.integral_points_count()
            Traceback (most recent call last):
            ...
            NotImplementedError: ...

        "Fibonacci" knapsacks (preprocessing helps a lot)::

            sage: def fibonacci_knapsack(d, b, backend=None):
            ....:     lp = MixedIntegerLinearProgram(base_ring=QQ)
            ....:     x = lp.new_variable(nonnegative=True)
            ....:     lp.add_constraint(lp.sum(fibonacci(i+3)*x[i] for i in range(d)) <= b)
            ....:     return lp.polyhedron(backend=backend)
            sage: fibonacci_knapsack(20, 12).integral_points_count() # does not finish with preprocess=False
            33

        TESTS:

        We check that :trac:`21491` is fixed::

            sage: P = Polyhedron(ieqs=[], eqns=[[-10,0,1],[-10,1,0]])
            sage: P.integral_points_count() # optional - latte_int
            1
            sage: P = Polyhedron(ieqs=[], eqns=[[-11,0,2],[-10,1,0]])
            sage: P.integral_points_count() # optional - latte_int
            0
        """
        if self.is_empty():
            return 0
        if self.dimension() == 0:
            return 1 if self.is_lattice_polytope() else 0
        if not self.is_compact():
            # LattE just prints the warning 'readCddExtFile:: Given polyhedron is unbounded!!!'
            # but then returns 0.
            raise NotImplementedError('Unbounded polyhedra are not supported')

        # for small bounding boxes, it is faster to naively iterate over the points of the box
        box_min, box_max = self.bounding_box(integral_hull=True)
        if box_min is None:
            return 0
        box_points = prod(max_coord-min_coord+1 for min_coord, max_coord in zip(box_min, box_max))

        if explicit_enumeration_threshold is None or box_points <= explicit_enumeration_threshold:
            return len(self.integral_points())

        p = self

        if preprocess:
            # If integral hull is known to lie in a coordinate hyperplane,
            # tighten bounds to reduce dimension.
            rat_box_min, rat_box_max = self.bounding_box(integral=False)
            if any(a == b and (ra < a or b < rb)
                   for ra, a, b, rb in zip(rat_box_min, box_min, box_max, rat_box_max)):
                lp, x = self.to_linear_program(return_variable=True)
                for i, a in enumerate(box_min):
                    lp.set_min(x[i], a)
                for i, b in enumerate(box_max):
                    lp.set_max(x[i], b)
                p = lp.polyhedron() # this recomputes the double description, which is wasteful
                if p.is_empty():
                    return 0
                if p.dimension() == 0:
                    return 1 if p.is_lattice_polytope() else 0

        # LattE does not like V-representation of lower-dimensional polytopes.
        if not p.is_full_dimensional():
            use_Hrepresentation = True

        from sage.interfaces.latte import count
        return count(
                p.cdd_Hrepresentation() if use_Hrepresentation else p.cdd_Vrepresentation(),
                cdd=True,
                verbose=verbose,
                **kwds)

    def ehrhart_polynomial(self,engine=None,variable='t',verbose=False,
            dual=None, irrational_primal=None, irrational_all_primal=None,
            maxdet=None, no_decomposition=None, compute_vertex_cones=None,
            smith_form=None, dualization=None, triangulation=None,
            triangulation_max_height=None, **kwds):
        r"""
        Return the Ehrhart polynomial of this polyhedron.

        The polyhedron must be a lattice polytope. Let `P` be a lattice
        polytope in `\RR^d` and define `L(P,t) = \# (tP\cap \ZZ^d)`.
        Then E. Ehrhart proved in 1962 that `L` coincides with a
        rational polynomial of degree `d` for integer `t`. `L` is called the
        *Ehrhart polynomial* of `P`. For more information see the
        :wikipedia:`Ehrhart_polynomial`. The Ehrhart polynomial may be computed
        using either LattE Integrale or Normaliz by setting ``engine``  to
        'latte' or 'normaliz' respectively.

        INPUT:

        - ``engine`` -- string; The backend to use. Allowed values are:

          * ``None`` (default); When no input is given the Ehrhart polynomial
            is computed using LattE Integrale (optional)
          * ``'latte'``; use LattE integrale program (optional)
          * ``'normaliz'``; use Normaliz program (optional package pynormaliz).
            The backend of ``self`` must be set to 'normaliz'.

        -  ``variable`` -- string (default: 't'); The variable in which the
           Ehrhart polynomial should be expressed.

        - When the ``engine`` is 'latte', the additional input values are:

          * ``verbose`` - boolean (default: ``False``); If ``True``, print the
            whole output of the LattE command.

          The following options are passed to the LattE command, for details
          consult `the LattE documentation
          <https://www.math.ucdavis.edu/~latte/software/packages/latte_current/>`__:

          * ``dual`` - boolean; triangulate and signed-decompose in the dual
            space
          * ``irrational_primal`` - boolean; triangulate in the dual space,
            signed-decompose in the primal space using irrationalization.
          * ``irrational_all_primal`` - boolean; triangulate and signed-decompose
            in the primal space using irrationalization.
          * ``maxdet`` -- integer; decompose down to an index (determinant) of
            ``maxdet`` instead of index 1 (unimodular cones).
          * ``no_decomposition`` -- boolean; do not signed-decompose
            simplicial cones.
          * ``compute_vertex_cones`` -- string; either 'cdd' or 'lrs' or '4ti2'
          * ``smith_form`` -- string; either 'ilio' or 'lidia'
          * ``dualization`` -- string; either 'cdd' or '4ti2'
          * ``triangulation`` - string; 'cddlib', '4ti2' or 'topcom'
          * ``triangulation_max_height`` - integer; use a uniform distribution
            of height from 1 to this number

        OUTPUT:

        A univariate polynomial in ``variable`` over a rational field.

        .. SEEALSO::

            :mod:`~sage.interfaces.latte` the interface to LattE Integrale
            `PyNormaliz <https://pypi.python.org/pypi/PyNormaliz/1.5>`_

        EXAMPLES:

        To start, we find the Ehrhart polynomial of a three-dimensional
        ``simplex``, first using ``engine='latte'``. Leaving the engine
        unspecified sets the ``engine`` to 'latte' by default::

            sage: simplex = Polyhedron(vertices=[(0,0,0),(3,3,3),(-3,2,1),(1,-1,-2)])
            sage: simplex = simplex.change_ring(QQ)
            sage: poly = simplex.ehrhart_polynomial(engine='latte')  # optional - latte_int
            sage: poly                                               # optional - latte_int
            7/2*t^3 + 2*t^2 - 1/2*t + 1
            sage: poly(1)                                            # optional - latte_int
            6
            sage: len(simplex.integral_points())                     # optional - latte_int
            6
            sage: poly(2)                                            # optional - latte_int
            36
            sage: len((2*simplex).integral_points())                 # optional - latte_int
            36

        Now we find the same Ehrhart polynomial, this time using
        ``engine='normaliz'``. To use the Normaliz engine, the ``simplex`` must
        be defined with ``backend='normaliz'``::

            sage: simplex = Polyhedron(vertices=[(0,0,0),(3,3,3),(-3,2,1),(1,-1,-2)], backend='normaliz') # optional - pynormaliz
            sage: simplex = simplex.change_ring(QQ)                                                       # optional - pynormaliz
            sage: poly = simplex.ehrhart_polynomial(engine = 'normaliz')                                  # optional - pynormaliz
            sage: poly                                                                                    # optional - pynormaliz
            7/2*t^3 + 2*t^2 - 1/2*t + 1

        If the ``engine='normaliz'``, the backend should be ``'normaliz'``, otherwise
        it returns an error::

            sage: simplex = Polyhedron(vertices=[(0,0,0),(3,3,3),(-3,2,1),(1,-1,-2)])
            sage: simplex = simplex.change_ring(QQ)
            sage: simplex.ehrhart_polynomial(engine='normaliz')  # optional - pynormaliz
            Traceback (most recent call last):
            ...
            TypeError: The backend of the polyhedron should be 'normaliz'

        The polyhedron should be compact::

            sage: C = Polyhedron(backend='normaliz',rays=[[1,2],[2,1]])  # optional - pynormaliz
            sage: C = C.change_ring(QQ)                                  # optional - pynormaliz
            sage: C.ehrhart_polynomial()                                 # optional - pynormaliz
            Traceback (most recent call last):
            ...
            ValueError: Ehrhart polynomial only defined for compact polyhedra

        The polyhedron should have integral vertices::

            sage: L = Polyhedron(vertices = [[0],[1/2]])
            sage: L.ehrhart_polynomial()
            Traceback (most recent call last):
            ...
            TypeError: the polytope has nonintegral vertices, use ehrhart_quasipolynomial with backend 'normaliz'
        """
        # check if ``self`` is compact and has vertices in ZZ
        if self.is_empty():
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            from sage.rings.rational_field import QQ
            R = PolynomialRing(QQ, variable)
            return R.zero()

        if not self.is_compact():
            raise ValueError("Ehrhart polynomial only defined for compact polyhedra")

        if any(not v.is_integral() for v in self.vertex_generator()):
            raise TypeError("the polytope has nonintegral vertices, use ehrhart_quasipolynomial with backend 'normaliz'")
        # Passes to specific latte or normaliz subfunction depending on engine
        if engine is None:
            # set default engine to latte
            engine = 'latte'
        if engine == 'latte':
            poly = self._ehrhart_polynomial_latte(verbose, dual,
            irrational_primal, irrational_all_primal, maxdet,
            no_decomposition, compute_vertex_cones, smith_form,
            dualization, triangulation, triangulation_max_height,
            **kwds)
            return poly.change_variable_name(variable)
            # TO DO: replace this change of variable by creating the appropriate
            #        polynomial ring in the latte interface.

        elif engine == 'normaliz':
            return self._ehrhart_polynomial_normaliz(variable)
        else:
            raise ValueError("engine must be 'latte' or 'normaliz'")

    def ehrhart_quasipolynomial(self, variable='t', engine=None, verbose=False,
            dual=None, irrational_primal=None, irrational_all_primal=None,
            maxdet=None, no_decomposition=None, compute_vertex_cones=None,
            smith_form=None, dualization=None, triangulation=None,
            triangulation_max_height=None, **kwds):
        r"""
        Compute the Ehrhart quasipolynomial of this polyhedron with rational
        vertices.

        If the polyhedron is a lattice polytope, returns the Ehrhart polynomial,
        a univariate polynomial in ``variable`` over a rational field.
        If the polyhedron  has rational, nonintegral vertices, returns a tuple
        of polynomials in ``variable`` over a rational field.
        The Ehrhart counting function of a polytope `P` with rational
        vertices is given by a *quasipolynomial*. That is, there exists a
        positive integer `l` and `l` polynomials
        `ehr_{P,i} \text{ for } i \in \{1,\dots,l \}` such that if `t` is
        equivalent to `i` mod `l` then `tP \cap \mathbb Z^d = ehr_{P,i}(t)`.

        INPUT:

        - ``variable`` -- string (default: 't'); The variable in which the
          Ehrhart polynomial should be expressed.

        - ``engine`` -- string; The backend to use. Allowed values are:

          * ``None`` (default); When no input is given the Ehrhart polynomial
            is computed using LattE Integrale (optional)
          * ``'latte'``; use LattE integrale program (optional)
          * ``'normaliz'``; use Normaliz program (optional package pynormaliz).
            The backend of ``self`` must be set to 'normaliz'.

        - When the ``engine`` is 'latte', the additional input values are:

          * ``verbose`` - boolean (default: ``False``); If ``True``, print the
            whole output of the LattE command.

          The following options are passed to the LattE command, for details
          consult `the LattE documentation
          <https://www.math.ucdavis.edu/~latte/software/packages/latte_current/>`__:

          * ``dual`` - boolean; triangulate and signed-decompose in the dual
            space
          * ``irrational_primal`` - boolean; triangulate in the dual space,
            signed-decompose in the primal space using irrationalization.
          * ``irrational_all_primal`` - boolean; triangulate and signed-decompose
            in the primal space using irrationalization.
          * ``maxdet`` -- integer; decompose down to an index (determinant) of
            ``maxdet`` instead of index 1 (unimodular cones).
          * ``no_decomposition`` -- boolean; do not signed-decompose
            simplicial cones.
          * ``compute_vertex_cones`` -- string; either 'cdd' or 'lrs' or '4ti2'
          * ``smith_form`` -- string; either 'ilio' or 'lidia'
          * ``dualization`` -- string; either 'cdd' or '4ti2'
          * ``triangulation`` - string; 'cddlib', '4ti2' or 'topcom'
          * ``triangulation_max_height`` - integer; use a uniform distribution of
            height from 1 to this number

        OUTPUT:

        A univariate polynomial over a rational field or a tuple of such
        polynomials.

        .. SEEALSO::

            :mod:`~sage.interfaces.latte` the interface to LattE Integrale
            `PyNormaliz <https://pypi.python.org/pypi/PyNormaliz/1.5>`_

        .. WARNING::

            If the polytope has rational, non integral vertices,
            it must have ``backend='normaliz'``.

        EXAMPLES:

        As a first example, consider the line segment [0,1/2]. If we
        dilate this line segment by an even integral factor `k`,
        then the dilated line segment will contain `k/2 +1` lattice points.
        If `k` is odd then there will be `k/2+1/2` lattice points in
        the dilated line segment. Note that it is necessary to set the
        backend of the polytope to 'normaliz'::

            sage: line_seg = Polyhedron(vertices=[[0],[1/2]],backend='normaliz') # optional - pynormaliz
            sage: line_seg                                                       # optional - pynormaliz
            A 1-dimensional polyhedron in QQ^1 defined as the convex hull of 2 vertices
            sage: line_seg.ehrhart_quasipolynomial()                             # optional - pynormaliz
            (1/2*t + 1, 1/2*t + 1/2)

        For a more exciting example, let us look at the subpolytope of the
        3 dimensional permutahedron fixed by the reflection
        across the hyperplane `x_1 = x_4`::

            sage: verts = [[3/2, 3, 4, 3/2],
            ....:  [3/2, 4, 3, 3/2],
            ....:  [5/2, 1, 4, 5/2],
            ....:  [5/2, 4, 1, 5/2],
            ....:  [7/2, 1, 2, 7/2],
            ....:  [7/2, 2, 1, 7/2]]
            sage: subpoly = Polyhedron(vertices=verts, backend='normaliz') # optional - pynormaliz
            sage: eq = subpoly.ehrhart_quasipolynomial()    # optional - pynormaliz
            sage: eq                                        # optional - pynormaliz
            (4*t^2 + 3*t + 1, 4*t^2 + 2*t)
            sage: eq = subpoly.ehrhart_quasipolynomial()    # optional - pynormaliz
            sage: eq                                        # optional - pynormaliz
            (4*t^2 + 3*t + 1, 4*t^2 + 2*t)
            sage: even_ep = eq[0]                           # optional - pynormaliz
            sage: odd_ep  = eq[1]                           # optional - pynormaliz
            sage: even_ep(2)                                # optional - pynormaliz
            23
            sage: ts = 2*subpoly                            # optional - pynormaliz
            sage: ts.integral_points_count()                # optional - pynormaliz latte_int
            23
            sage: odd_ep(1)                                 # optional - pynormaliz
            6
            sage: subpoly.integral_points_count()           # optional - pynormaliz latte_int
            6

        A polytope with rational nonintegral vertices must have
        ``backend='normaliz'``::

            sage: line_seg = Polyhedron(vertices=[[0],[1/2]])
            sage: line_seg.ehrhart_quasipolynomial()
            Traceback (most recent call last):
            ...
            TypeError: The backend of the polyhedron should be 'normaliz'

        The polyhedron should be compact::

            sage: C = Polyhedron(backend='normaliz',rays=[[1/2,2],[2,1]])  # optional - pynormaliz
            sage: C.ehrhart_quasipolynomial()                              # optional - pynormaliz
            Traceback (most recent call last):
            ...
            ValueError: Ehrhart quasipolynomial only defined for compact polyhedra

        If the polytope happens to be a lattice polytope, the Ehrhart
        polynomial is returned::

            sage: simplex = Polyhedron(vertices=[(0,0,0),(3,3,3),(-3,2,1),(1,-1,-2)], backend='normaliz') # optional - pynormaliz
            sage: simplex = simplex.change_ring(QQ)                                                       # optional - pynormaliz
            sage: poly = simplex.ehrhart_quasipolynomial(engine='normaliz')                               # optional - pynormaliz
            sage: poly                                                                                    # optional - pynormaliz
            7/2*t^3 + 2*t^2 - 1/2*t + 1
            sage: simplex.ehrhart_polynomial()                                                            # optional - pynormaliz latte_int
            7/2*t^3 + 2*t^2 - 1/2*t + 1
        """
        if self.is_empty():
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            from sage.rings.rational_field import QQ
            R = PolynomialRing(QQ, 't')
            return R.zero()

        if not self.is_compact():
            raise ValueError("Ehrhart quasipolynomial only defined for compact polyhedra")

        if engine is None:
            # setting the default to 'normaliz'
            engine = 'normaliz'
        if engine == 'normaliz':
            return self._ehrhart_quasipolynomial_normaliz(variable)
        if engine == 'latte':
            if any(not v.is_integral() for v in self.vertex_generator()):
                raise TypeError("the polytope has nonintegral vertices, the engine and backend of self should be 'normaliz'")
            poly = self._ehrhart_polynomial_latte(verbose, dual,
            irrational_primal, irrational_all_primal, maxdet,
            no_decomposition, compute_vertex_cones, smith_form,
            dualization, triangulation, triangulation_max_height,
            **kwds)
            return poly.change_variable_name(variable)
            # TO DO: replace this change of variable by creating the appropriate
            #        polynomial ring in the latte interface.
        else:
            raise TypeError("the engine should be 'latte' or 'normaliz'")

    def _ehrhart_quasipolynomial_normaliz(self, variable='t'):
        r"""
        Compute the Ehrhart quasipolynomial of a lattice or rational polytope
        using the Normaliz engine.

        INPUT:

        - ``variable`` -- string (default: 't'); The variable in which the
          Ehrhart polynomial is expressed.

        OUTPUT:

        A univariate polynomial over a rational field or a tuple of such
        polynomials.

        EXAMPLES:

        The subpolytope of the 3 dimensional permutahedron fixed by the
        reflection across the hyperplane `x_1 = x_4`::

            sage: verts = [[3/2, 3, 4, 3/2],
            ....:  [3/2, 4, 3, 3/2],
            ....:  [5/2, 1, 4, 5/2],
            ....:  [5/2, 4, 1, 5/2],
            ....:  [7/2, 1, 2, 7/2],
            ....:  [7/2, 2, 1, 7/2]]
            sage: subpoly = Polyhedron(vertices=verts, backend='normaliz')         # optional - pynormaliz
            sage: eq = subpoly._ehrhart_quasipolynomial_normaliz()  # optional - pynormaliz
            sage: eq                                                # optional - pynormaliz
            (4*t^2 + 3*t + 1, 4*t^2 + 2*t)
            sage: even_ep = eq[0]                                   # optional - pynormaliz
            sage: odd_ep  = eq[1]                                   # optional - pynormaliz
            sage: even_ep(2)                                        # optional - pynormaliz
            23
            sage: ts = 2*subpoly                                    # optional - pynormaliz
            sage: ts.integral_points_count()                        # optional - pynormaliz latte_int
            23
            sage: odd_ep(1)                                         # optional - pynormaliz
            6
            sage: subpoly.integral_points_count()                   # optional - pynormaliz latte_int
            6

        TESTS::

            sage: line_seg = Polyhedron(vertices=[[0],[1/2]])
            sage: line_seg._ehrhart_quasipolynomial_normaliz()      # optional - pynormaliz
            Traceback (most recent call last):
            ...
            TypeError: The backend of the polyhedron should be 'normaliz'
        """
        raise TypeError("The backend of the polyhedron should be 'normaliz'")

    _ehrhart_polynomial_normaliz = _ehrhart_quasipolynomial_normaliz

    def _ehrhart_polynomial_latte(self, verbose=False, dual=None,
            irrational_primal=None, irrational_all_primal=None, maxdet=None,
            no_decomposition=None, compute_vertex_cones=None, smith_form=None,
            dualization=None, triangulation=None, triangulation_max_height=None,
            **kwds):
        r"""
        Return the Ehrhart polynomial of this polyhedron using LattE integrale.

        Let `P` be a lattice polytope in `\RR^d` and define `L(P,t) = \# (tP
        \cap \ZZ^d)`. Then E. Ehrhart proved in 1962 that `L` coincides with a
        rational polynomial of degree `d` for integer `t`. `L` is called the
        *Ehrhart polynomial* of `P`. For more information see the
        :wikipedia:`Ehrhart_polynomial`.

        INPUT:

        - ``verbose`` - boolean (default: ``False``); if ``True``, print the
          whole output of the LattE command.

        The following options are passed to the LattE command, for details you
        should consult `the LattE documentation
        <https://www.math.ucdavis.edu/~latte/software/packages/latte_current/>`__:

        - ``dual`` - boolean; triangulate and signed-decompose in the dual
          space

        - ``irrational_primal`` - boolean; triangulate in the dual space,
          signed-decompose in the primal space using irrationalization.

        - ``irrational_all_primal`` - boolean; triangulate and signed-decompose
          in the primal space using irrationalization.

        - ``maxdet`` -- integer; decompose down to an index (determinant) of
          ``maxdet`` instead of index 1 (unimodular cones).

        - ``no_decomposition`` -- boolean; do not signed-decompose simplicial cones.

        - ``compute_vertex_cones`` -- string; either 'cdd' or 'lrs' or '4ti2'

        - ``smith_form`` -- string; either 'ilio' or 'lidia'

        - ``dualization`` -- string; either 'cdd' or '4ti2'

        - ``triangulation`` - string; 'cddlib', '4ti2' or 'topcom'

        - ``triangulation_max_height`` - integer; use a uniform distribution of
          height from 1 to this number

        .. NOTE::

            Any additional argument is forwarded to LattE's executable
            ``count``. All occurrences of '_' will be replaced with a '-'.

        OUTPUT:

        A univariate polynomial over a rational field.

        ALGORITHM:

        This method calls the program ``count`` from LattE integrale, a program
        for lattice point enumeration (see
        https://www.math.ucdavis.edu/~latte/).

        .. SEEALSO::

            :mod:`~sage.interfaces.latte` the interface to LattE integrale

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(0,0,0),(3,3,3),(-3,2,1),(1,-1,-2)])
            sage: p = P._ehrhart_polynomial_latte()    # optional - latte_int
            sage: p                                    # optional - latte_int
            7/2*t^3 + 2*t^2 - 1/2*t + 1
            sage: p(1)                                 # optional - latte_int
            6
            sage: len(P.integral_points())             # optional - latte_int
            6
            sage: p(2)                                 # optional - latte_int
            36
            sage: len((2*P).integral_points())         # optional - latte_int
            36

        The unit hypercubes::

            sage: from itertools import product
            sage: def hypercube(d):
            ....:     return Polyhedron(vertices=list(product([0,1],repeat=d)))
            sage: hypercube(3)._ehrhart_polynomial_latte()   # optional - latte_int
            t^3 + 3*t^2 + 3*t + 1
            sage: hypercube(4)._ehrhart_polynomial_latte()   # optional - latte_int
            t^4 + 4*t^3 + 6*t^2 + 4*t + 1
            sage: hypercube(5)._ehrhart_polynomial_latte()   # optional - latte_int
            t^5 + 5*t^4 + 10*t^3 + 10*t^2 + 5*t + 1
            sage: hypercube(6)._ehrhart_polynomial_latte()   # optional - latte_int
            t^6 + 6*t^5 + 15*t^4 + 20*t^3 + 15*t^2 + 6*t + 1

        TESTS:

        Test options::

            sage: P = Polyhedron(ieqs=[[1,-1,1,0], [-1,2,-1,0], [1,1,-2,0]], eqns=[[-1,2,-1,-3]], base_ring=ZZ)

            sage: p = P._ehrhart_polynomial_latte(maxdet=5, verbose=True)  # optional - latte_int
            This is LattE integrale ...
            ...
            Invocation: count --ehrhart-polynomial '--redundancy-check=none' --cdd '--maxdet=5' /dev/stdin
            ...
            sage: p    # optional - latte_int
            1/2*t^2 + 3/2*t + 1

            sage: p = P._ehrhart_polynomial_latte(dual=True, verbose=True)  # optional - latte_int
            This is LattE integrale ...
            ...
            Invocation: count --ehrhart-polynomial '--redundancy-check=none' --cdd --dual /dev/stdin
            ...
            sage: p   # optional - latte_int
            1/2*t^2 + 3/2*t + 1

            sage: p = P._ehrhart_polynomial_latte(irrational_primal=True, verbose=True)   # optional - latte_int
            This is LattE integrale ...
            ...
            Invocation: count --ehrhart-polynomial '--redundancy-check=none' --cdd --irrational-primal /dev/stdin
            ...
            sage: p   # optional - latte_int
            1/2*t^2 + 3/2*t + 1

            sage: p = P._ehrhart_polynomial_latte(irrational_all_primal=True, verbose=True)  # optional - latte_int
            This is LattE integrale ...
            ...
            Invocation: count --ehrhart-polynomial '--redundancy-check=none' --cdd --irrational-all-primal /dev/stdin
            ...
            sage: p   # optional - latte_int
            1/2*t^2 + 3/2*t + 1

        Test bad options::

            sage: P._ehrhart_polynomial_latte(bim_bam_boum=19)   # optional - latte_int
            Traceback (most recent call last):
            ...
            RuntimeError: LattE integrale program failed (exit code 1):
            ...
            Invocation: count --ehrhart-polynomial '--redundancy-check=none' --cdd '--bim-bam-boum=19' /dev/stdin
            Unknown command/option --bim-bam-boum=19
        """

        # note: the options below are explicitely written in the function
        # declaration in order to keep tab completion (see #18211).
        kwds.update({
            'dual'                    : dual,
            'irrational_primal'       : irrational_primal,
            'irrational_all_primal'   : irrational_all_primal,
            'maxdet'                  : maxdet,
            'no_decomposition'        : no_decomposition,
            'compute_vertex_cones'    : compute_vertex_cones,
            'smith_form'              : smith_form,
            'dualization'             : dualization,
            'triangulation'           : triangulation,
            'triangulation_max_height': triangulation_max_height})

        from sage.interfaces.latte import count
        ine = self.cdd_Hrepresentation()
        return count(ine, cdd=True, ehrhart_polynomial=True, verbose=verbose, **kwds)
