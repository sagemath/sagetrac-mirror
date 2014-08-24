r"""
Products of Projective Spaces

Products of projective spaces of varying dimension are convenient
ambient spaces for complete intersections. Group actions on them, and
the interplay with representation theory, provide many interesting
examples of algebraic varieties.

EXAMPLES::

    sage: X = ProductProjectiveSpaces([3, 3], QQ);  X
    Product of projective spaces P^3 x P^3 over Rational Field
    sage: X.dimension()
    6

The coordinate names can optionally be specified::

    sage: P2xP2 = ProductProjectiveSpaces([2, 2], QQ, names=['x', 'y'])
    sage: P2xP2.coordinate_ring().inject_variables()
    Defining x0, x1, x2, y0, y1, y2

    sage: P2xP2.<x,y> = ProductProjectiveSpaces([2, 2], QQ)
    sage: P2xP2
    Product of projective spaces P^2 x P^2 over Rational Field
    sage: P2xP2.coordinate_ring()
    Multivariate Polynomial Ring in x0, x1, x2, y0, y1, y2 over Rational Field

    sage: P2xP2.<x,y,z, r,s,t> = ProductProjectiveSpaces([2, 2], QQ)
    sage: P2xP2.coordinate_ring()
    Multivariate Polynomial Ring in x, y, z, r, s, t over Rational Field
"""

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.all import (PolynomialRing, ZZ, QQ)
from sage.rings.commutative_ring import is_CommutativeRing
from sage.misc.cachefunc import cached_method
from sage.structure.parent_gens import normalize_names
from sage.schemes.generic.ambient_space import AmbientSpace


def ProductProjectiveSpaces(dim_list, ring=None, names='x'):
    r"""
    Construct product of projective spaces.

    INPUT:

    - ``dim_list`` -- list/tuple/iterator of integers. The dimension
      of the factors.

    - ``ring`` -- ring (default: ``QQ``). The base ring.

    - ``names`` -- string or list/tuple/iterable of strings. The
      variable names.
    

    EXAMPLES::

        sage: ProductProjectiveSpaces([3, 2, 1], QQ)
        Product of projective spaces P^3 x P^2 x P^1 over Rational Field
    """
    dim_list = tuple(map(ZZ, dim_list))
    if ring is None:
        ring = QQ
    n_vars = sum(1+d for d in dim_list)
    if isinstance(names, basestring):
        names = normalize_names(n_vars, names)
    else:
        name_list = list(names)
        if len(name_list) == len(dim_list):
            names = []
            for name, dim in zip(name_list, dim_list):
                names += normalize_names(dim+1, name)
        else:
            n_vars = sum(1+d for d in dim_list)
            names = normalize_names(n_vars, name_list)
    if not is_CommutativeRing(ring):
        raise ValueError('ring must be commutative')
    return ProductProjectiveSpaces_ring(dim_list, ring, names)


class ProductProjectiveSpaces_ring(AmbientSpace):

    def __init__(self, dim_tuple, ring, names=None):
        """
        Product of Projective Spaces

        INPUT:

        - ``dim_tuple`` -- tuple of integers. The dimenson of the
          individual projective spaces.

        - ``ring`` -- a commutative ring. The base ring.

        - ``names`` -- optional string or list/tuple/iterable of
          strings. The names of homogeneous variables.

        EXAMPLES::
    
            sage: X.<x,y,z,w> = ProductProjectiveSpaces([1, 1], QQ)
            sage: X.base_scheme()
            Spectrum of Rational Field
            sage: X.base_ring()
            Rational Field
            sage: X.structure_morphism()
            Scheme morphism:
              From: Product of projective spaces P^1 x P^1 over Rational Field
              To:   Spectrum of Rational Field
              Defn: Structure map
            sage: X.coordinate_ring()
            Multivariate Polynomial Ring in x, y, z, w over Rational Field
    
        Loading and saving::

            sage: loads(X.dumps()) == X
            True
        """
        assert isinstance(dim_tuple, tuple)
        assert all(x.parent() is ZZ for x in dim_tuple)
        assert is_CommutativeRing(ring)
        AmbientSpace.__init__(self, sum(dim_tuple), ring)
        self._assign_names(names)
        self._dim = dim_tuple

    def _iter_factors(self, homogeneous_coordinates):
        """
        Split homogeneous coordinates into the factor projective spaces.

        INPUT:

        - ``homogeneous_coordinates`` -- list/tuple/iterable for
          :meth:`ngens` items.

        OUTPUT:

        Iterate over the factors and yield lists of length equal to
        the number of the homogeneous coordinates of that factor.

        EXAMPLES::

            sage: Y = ProductProjectiveSpaces([1, 2, 4, 8])
            sage: Y._iter_factors(range(19)).next()
            [0, 1]
            sage: tuple(Y._iter_factors(range(19)))
            ([0, 1], 
             [2, 3, 4], 
             [5, 6, 7, 8, 9], 
             [10, 11, 12, 13, 14, 15, 16, 17, 18])
        """
        homogeneous = list(homogeneous_coordinates)
        for d in self._dim:
            yield homogeneous[0:d+1]
            homogeneous = homogeneous[d+1:]
        if len(homogeneous) != 0:
            raise ValueError('got more items than homogenous coordinates')

    @cached_method
    def factor(self, i=None):
        """
        Return a tuple of projective space factors.

        INPUT:

        - ``i`` -- integer or ``None`` (default). The factor to
          return. 

        OUTPUT:

        If ``i`` is not specified, all component projective spaces in
        a tuple. If the index ``i`` is specified, the ``i``-th factor
        is returned.

        EXAMPLES::

            sage: Y.<a, b, c, d> = ProductProjectiveSpaces([1, 2, 4, 8])
            sage: Y.factor()
            (Projective Space of dimension 1 over Rational Field,
             Projective Space of dimension 2 over Rational Field,
             Projective Space of dimension 4 over Rational Field,
             Projective Space of dimension 8 over Rational Field)
            sage: Y.factor(0).coordinate_ring()
            Multivariate Polynomial Ring in a0, a1 over Rational Field
            sage: Y.factor(1).coordinate_ring()
            Multivariate Polynomial Ring in b0, b1, b2 over Rational Field
            sage: Y.factor(2).coordinate_ring()
            Multivariate Polynomial Ring in c0, c1, c2, c3, c4 over Rational Field
            sage: Y.factor(3).coordinate_ring()
            Multivariate Polynomial Ring in d0, d1, d2, d3, d4, d5, d6, d7, d8 over Rational Field
        """
        if i is not None:
            return self.factor()[i]
        from sage.schemes.projective.projective_space import ProjectiveSpace
        result = []
        name_iter = self._iter_factors(self.variable_names())
        for n in self.factor_dim():
            Pn = ProjectiveSpace(n, self.base_ring(), names=next(name_iter))
            result.append(Pn)
        return tuple(result)

    def factor_dim(self, i=None):
        """
        Return the dimension of the ``i``-th factor projective space

        INPUT:

        - ``i`` -- integer or ``None`` (default). The index of the
          projective space factor.

        OUTPUT:

        * If ``i`` is an integer: Integer. The dimension of the
          ``i``-th factor.

        * If ``i`` is ``None``: Tuple of integers. The dimensions of
          the factors.

        EXAMPLES::

            sage: Y = ProductProjectiveSpaces([1, 2, 4, 8])
            sage: Y.factor_dim()
            (1, 2, 4, 8)
            sage: Y.factor_dim(0)
            1
            sage: Y.factor_dim(1)
            2
            sage: Y.factor_dim(2)
            4
            sage: Y.factor_dim(3)
            8
        """
        if i is None:
            return self._dim
        else:
            return self._dim[i]

    def num_factors(self):
        """
        Return the number of factor projective space

        OUTPUT:

        Integer. The number of projective spaces that the product
        consists of.

        EXAMPLES::

            sage: Y = ProductProjectiveSpaces([1, 2, 4, 8])
            sage: Y.num_factors()
            4
        """
        return len(self._dim)

    def ngens(self):
        """
        Return the number of generators of self, i.e. the number of
        variables in the coordinate ring of self.

        EXAMPLES::

            sage: ProductProjectiveSpaces([3, 4], QQ).ngens()
            9
            sage: ProductProjectiveSpaces([1, 7], ZZ).ngens()
            10
        """
        return self.dimension_relative() + self.num_factors()

    def _check_satisfies_equations(self, v):
        """
        Check whether ``v`` defines a point.

        OUTPUT:

        ``True`` if ``v`` defines a point on the scheme self. Raise a
        TypeError otherwise.

        EXAMPLES::

            sage: P = ProductProjectiveSpaces([2, 2], ZZ)
            sage: P._check_satisfies_equations([1, 0, 0, 0, 0, 1])
            True

            sage: P = ProductProjectiveSpaces([1, 1], QQ)
            sage: P._check_satisfies_equations((1/2, 0, 1/2, 0))
            True

            sage: P = ProductProjectiveSpaces([1, 1], ZZ)
            sage: P._check_satisfies_equations([1, 1, 0, 0])
            Traceback (most recent call last):
            ...
            TypeError: The zero vector is not a point in projective space
        """
        for i, v_i in enumerate(self._iter_factors(v)):
            self.factor(i)._check_satisfies_equations(v_i)
        return True

    @cached_method
    def coordinate_ring(self):
        """
        Return the coordinate ring of this scheme.

        OUTPUT:

        Polynomial ring in the homogeneous variables.

        EXAMPLES::

            sage: ProductProjectiveSpaces([3,2]).coordinate_ring()
            Multivariate Polynomial Ring in x0, x1, x2, x3, x4, x5, x6 over Rational Field

            sage: X.<a,b,c,d> = ProductProjectiveSpaces([1, 1], GF(19^2,'alpha'))
            sage: X.coordinate_ring()
            Multivariate Polynomial Ring in a, b, c, d over Finite Field in alpha of size 19^2
        """
        return PolynomialRing(self.base_ring(),
                              self.variable_names(), 
                              self.ngens())

    def _degree(self, polynomial):
        """
        Return the homogeneous degrees.

        INPUT:

        A polynomial in :meth:`coordinate_ring`.

        OUTPUT:

        A tuple of integers, one for each projective space factor. A
        ``ValueError`` is raised if the polynomial is not homogenous.

        EXAMPLES::

            sage: P1xP1.<x,y, s,t> = ProductProjectiveSpaces([1,1], QQ)
            sage: P1xP1._degree(x^2*t + y^2*s)
            (2, 1)
            sage: P1xP1._degree(x + s)
            Traceback (most recent call last):
            ...
            ValueError: polynomial is not multi-graded homogeneous
        """
        def monomial_degree(d):
            return tuple(sum(d_i) for d_i in self._iter_factors(d))
        summand_degrees = map(monomial_degree, polynomial.exponents())
        result = summand_degrees[0]
        if not all(result == deg for deg in summand_degrees):
            raise ValueError('polynomial is not multi-graded homogeneous')
        return result

    def _validate(self, polynomials):
        """
        Validate that each polynomial is graded homogeneous.

        INPUT:

        - ``polynomials`` -- list/tuple of polynomials in the
          coordinate ring.

        OUTPUT:

        ``polynomials`` is returned unmodified if they are all graded
        homogeneous. A ``ValueError`` is raised otherwise.

        EXAMPLES::

            sage: P.<x, y, z, s, t> = ProductProjectiveSpaces([2, 1], ZZ)
            sage: P._validate([x*y*s - z^2*t, x*t])
            [x*y*s - z^2*t, x*t]

            sage: P._validate((x*y*s - z*t, x*t))
            Traceback (most recent call last):
            ...
            ValueError: polynomial is not multi-graded homogeneous

            sage: P._validate(x*y - z)
            Traceback (most recent call last):
            ...
            ValueError: argument must be list or tuple
        """
        if not isinstance(polynomials, (list, tuple)):
            raise ValueError('argument must be list or tuple')
        for f in polynomials:
            self._degree(f)    # raises ValueError if not homegeneous
        return polynomials

    def __cmp__(left, right):
        """
        Compare two products of projective spaces.

        INPUT:

        - ``right`` -- anything.

        OUTPUT:

        ``-1``, ``1``, or ``+1`` depending how ``left`` and ``right``
        compare.

        EXAMPLES::

            sage: ProductProjectiveSpaces([3, 1], QQ, ['a', 'b']) \
            ....:    == ProductProjectiveSpaces([3, 1], QQ, ['a', 'b'])
            True
            sage: ProductProjectiveSpaces([3, 1], QQ, ['a', 'b']) \
            ....:    == ProductProjectiveSpaces([3, 1], ZZ, ['a', 'b'])
            False
            sage: ProductProjectiveSpaces([3, 1], QQ, ['a', 'b']) \
            ....:    == AffineSpace(ZZ, 2, 'a')
            False
            sage: P = ProductProjectiveSpaces([3, 1], QQ, ['a', 'b'])
            sage: loads(P.dumps()) == P
            True
        """
        if not isinstance(right, ProductProjectiveSpaces_ring):
            return -1
        c = cmp(left.base_ring(), right.base_ring())
        if c != 0:
            return c
        c = cmp(left.factor_dim(), right.factor_dim())
        if c != 0:
            return c
        return cmp(left.coordinate_ring(), right.coordinate_ring())

    def _latex_(self):
        r"""
        Return a LaTeX representation of this projective space.

        EXAMPLES::

            sage: ProductProjectiveSpaces([1, 2, 3], ZZ, 'x')._latex_()
            '{\\mathbf P}_{\\Bold{Z}}^1 \\times 
             {\\mathbf P}_{\\Bold{Z}}^2 \\times 
             {\\mathbf P}_{\\Bold{Z}}^3'
        """
        return r' \times '.join(P._latex_() for P in self.factor())

    def _repr_(self):
        """
        Return a string representation of this projective space.

        EXAMPLES::

            sage: ProjectiveSpace(1, ZZ, 'x')
            Projective Space of dimension 1 over Integer Ring

        TESTS::

            sage: ProjectiveSpace(3, Zp(5), 'y')._repr_()
            'Projective Space of dimension 3 over 5-adic Ring with capped relative precision 20'
        """
        return ''.join([
            'Product of projective spaces ',
            ' x '.join(['P^{0}'.format(d) for d in self._dim]),
            ' over ',
            str(self.base_ring())])

    def _repr_generic_point(self, v=None):
        """
        Return a string representation of the point.

        INPUT:

        - ``v`` -- list/tuple/iterable of homogeneous coordinates or
          ``None`` (default). If ``None``, the representation of the
          generic point of the projective space is returned.

        EXAMPLES::

            sage: P.<x, y, s, t> = ProductProjectiveSpaces([1, 1], ZZ)
            sage: P._repr_generic_point([1, 2, 3, 4])
            '(1 : 2 | 3 : 4)'
            sage: P._repr_generic_point()
            '(x : y | s : t)'
        """
        if v is None:
            v = self.gens()
        factors = []
        for v_i in self._iter_factors(v):
            factors.append(' : '.join(map(str, v_i)))
        return '(' + ' | '.join(factors) + ')'

    def _latex_generic_point(self, v=None):
        """
        Return a LaTeX representation of the generic point
        corresponding to the list of polys on this projective space.

        If polys is None, the representation of the generic point of
        the projective space is returned.

        EXAMPLES::

            sage: P.<x, y, s, t> = ProductProjectiveSpaces([1, 1], ZZ)
            sage: P._latex_generic_point([1, 2, 3, 4])
            '\\left(1 : 2 | 3 : 4\\right)'
            sage: P._latex_generic_point()
            '\\left(x : y | s : t\\right)'
        """
        if v is None:
            v = self.gens()
        if v is None:
            v = self.gens()
        factors = []
        for v_i in self._iter_factors(v):
            factors.append(' : '.join([x._latex_() for x in v_i]))
        return r'\left(' + ' | '.join(factors) + r'\right)'

    def change_ring(self, ring):
        r"""
        Return a new product of projetive spaces with the base ring changed.

        INPUT:

        - ``ring`` -- commutative ring.

        OUTPUT:

        Product of projective spaces over ``ring``.

        EXAMPLES::

            sage: P.<x, y, s, t> = ProductProjectiveSpaces([1, 1], ZZ)
            sage: PQ = P.change_ring(QQ); PQ
            Product of projective spaces P^1 x P^1 over Rational Field
            sage: PQ.change_ring(GF(5))
            Product of projective spaces P^1 x P^1 over Finite Field of size 5
        """
        return ProductProjectiveSpaces(self.factor_dim(), ring, self.variable_names())

    def _point_homset(self, *args, **kwds):
        """
        Construct a point Hom-set.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(2, GF(3))
            sage: P2._point_homset(Spec(GF(3)), P2)
            Set of rational points of Projective Space of dimension 2 over Finite Field of size 3
        """
        from sage.schemes.product_projective.homset \
            import SchemeHomset_points_product_projective_spaces_ring
        return SchemeHomset_points_product_projective_spaces_ring(*args, **kwds)

    def _point(self, *args, **kwds):
        """
        Construct a point.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(2, GF(3))
            sage: point_homset = P2._point_homset(Spec(GF(3)), P2)
            sage: P2._point(point_homset, [1,2,3])
            (2 : 1 : 0)
        """
        from sage.schemes.product_projective.point \
            import SchemeMorphism_point_product_projective_spaces_ring
        return SchemeMorphism_point_product_projective_spaces_ring(*args, **kwds)

    def _an_element_(self):
        r"""
        Returns a (preferably typical) element of ``self``.

        This is used both for illustration and testing purposes.

        OUTPUT: a point in the projective space ``self``.

        EXAMPLES::

            sage: ProductProjectiveSpaces([1,2,3], ZZ).an_element()
            (0 : 1 | 2 : 3 : 4 | 5 : 6 : 7 : 8)
            sage: ProductProjectiveSpaces([3, 2, 1], PolynomialRing(ZZ,'y')).an_element()
            (0 : y : 2*y : 3*y | 4*y : 5*y : 6*y | 7*y : 8*y)
        """
        R = self.base_ring()
        v = [i * R.an_element() for i in range(self.ngens())]
        if all(v_i == 0 for v_i in v):
            v[0] = R.one()
        return self(v)

    def complete_intersection(self, *degrees):
        """
        Return complete intersection defined by degrees.

        INPUT:

        - ``degrees`` - list/tuple/iterable of length equals to the
          number of factors. Each entry is either an integer or a
          list/tuple/iterables of integer. The degree(s) of the
          equation(s) in the corresponding projective space factor.

        EXAMPLES::

            sage: P2xP2.<x,y> = ProductProjectiveSpaces([2, 2], QQ)
            sage: P2xP2.complete_intersection(3, 3)
            Complete intersection in Product of projective spaces P^2 x P^2 over Rational Field
              P^2   |   3
              P^2   |   3
        """
        from sage.schemes.product_projective.complete_intersection \
            import GenericCompleteIntersection
        return GenericCompleteIntersection(self, degrees)
        
