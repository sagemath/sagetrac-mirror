"""
Complete Intersections in Products of Projective Spaces

A complete intersection is a subscheme where the equations are
everywhere transverse. This assumption allows us to deduce many
properties of the subvariety without knowing the actual equations,
just from the degrees of the equations. To construct such a "generic"
complete intersection, you just have to pass the degrees::

    sage: P2xP2.<x,y> = ProductProjectiveSpaces([2, 2], QQ)
    sage: P2xP2.complete_intersection(3, 3)
    Complete intersection in Product of projective spaces P^2 x P^2 over Rational Field
      P^2   |   3
      P^2   |   3
 """

from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme



class GenericCompleteIntersection(AlgebraicScheme):

    def __init__(self, product_projective_spaces, degrees):
        """
        Generic complete intersection only defined by degrees of equations.

        INPUT:

        - ``product_projective_spaces`` -- a product of projective
          spaces. The ambient space.

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
        super(GenericCompleteIntersection, self).__init__(product_projective_spaces)
        def _normalize(d):
            try:
                return tuple([int(d)])
            except (TypeError, ValueError):
                return tuple(int(d_i) for d_i in d)
        degrees = map(_normalize, degrees)
        codim = len(degrees[0])
        if not all(codim == len(d) for d in degrees):
            raise ValueError('number of equations mismatch')
        self._codim = codim
        self._deg = degrees

    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.

        EXAMPLES::

            sage: P1xP1.<x,y> = ProductProjectiveSpaces([1, 1], QQ)
            sage: P1xP1.complete_intersection(3, 3)._repr_()
            'Complete intersection in Product of projective spaces P^1 x P^1 
             over Rational Field\n  P^1   |   3\n  P^1   |   3'
        """
        header = 'Complete intersection in {0}\n'.format(self.ambient_space())
        rows = []
        for dim, degrees in zip(self.ambient_space().factor_dim(), self._deg):
            rows.append(['P^{0}'.format(dim), '|'] + list(degrees))
        from sage.misc.table import table
        return header + str(table(rows))

    def dimension(self):
        """
        Return the dimension of the complete intersection.

        OUTPUT:

        Integer.

        EXAMPLES::

           sage: P2xP2.<x,y> = ProductProjectiveSpaces([2, 2], QQ)
           sage: P2xP2.complete_intersection(3, 3).dimension()
           3
        """
        return self.ambient_space().dimension() - self._codim

    def matrix(self):
        """
        Return the degree matrix.

        OUTPUT:

        Integer matrix.

        EXAMPLES::

           sage: P2xP3.<x,y> = ProductProjectiveSpaces([2, 3], QQ)
           sage: X = P2xP3.complete_intersection([1, 2], [3, 1]);  X
           Complete intersection in Product of projective spaces P^2 x P^3 over Rational Field
             P^2   |   1   2
             P^3   |   3   1
           sage: X.matrix()
           [1 2]
           [3 1]
           sage: matrix(X)           
           [1 2]
           [3 1]
        """
        return matrix(ZZ, self._deg)

    _matrix_ = matrix

    def degree(self, row):
        """
        Return the degree(s) of the equations for a give projective space factor.

        INPUT:

        - ``row`` -- integer. Index for one of the factor projective
          spaces. The row of :meth:`matrix`.

        OUTPUT:

        Tuple of integers. The degrees of the equation(s).

            sage: P2xP3.<x,y> = ProductProjectiveSpaces([2, 3], QQ)
            sage: X = P2xP3.complete_intersection([1, 2], [3, 1]);  X
            Complete intersection in Product of projective spaces P^2 x P^3 over Rational Field
              P^2   |   1   2
              P^3   |   3   1
            sage: X.degree(0)
            (1, 2)
            sage: X.degree(1)
            (3, 1)
        """
        return self._deg[row]

    def equation(self, col):
        """
        Return the degree(s) in the projective space factor(s) of an equation.

        INPUT:

        - ``col`` -- integer. Index for one of the equations. The
          column of :meth:`matrix`.

        OUTPUT:

        Tuple of integers. The degrees of the equation(s).

        EXAMPLES::

            sage: P2xP3.<x,y> = ProductProjectiveSpaces([2, 3], QQ)
            sage: X = P2xP3.complete_intersection([1, 2], [3, 1]);  X
            Complete intersection in Product of projective spaces P^2 x P^3 over Rational Field
              P^2   |   1   2
              P^3   |   3   1
            sage: X.equation(0)
            (1, 3)
            sage: X.equation(1)
            (2, 1)
        """
        return tuple(deg[col] for deg in self._deg)
