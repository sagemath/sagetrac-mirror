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

    def num_equations(self):
        """
        Return the number of equations.

        OUTPUT:
        
        Integer. The number of equations. Equals the codimension in
        the ambient space, since this is a complete intersection.

        EXAMPLES::

            sage: P2xP3.<x,y> = ProductProjectiveSpaces([2, 3], QQ)
            sage: X = P2xP3.complete_intersection([1, 2], [3, 1])
            sage: X.num_equations()
            2
            sage: X.ambient_space().dimension() - X.dimension()
            2
        """
        return self._codim

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

    def __cmp__(left, right):
        """
        Compare two complete intersections.

        INPUT:

        - ``right`` -- anything.

        OUTPUT:

        ``-1``, ``1``, or ``+1`` depending how ``left`` and ``right``
        compare.

        EXAMPLES::

            sage: P.<a,b> = ProductProjectiveSpaces([3, 1], QQ)
            sage: X = P.complete_intersection(2, 2)
            sage: P.complete_intersection(2, 2) == X
            True
            sage: loads(dumps(X)) == X
            True
            sage: X == P.complete_intersection(2, 3)
            False
        """
        if not isinstance(right, GenericCompleteIntersection):
            return -1
        c = cmp(left.ambient_space(), right.ambient_space())
        if c != 0:
            return c
        c = cmp(left._codim, right._codim)
        if c != 0:
            return c
        return cmp(left._deg, right._deg)

    def projective_space_permutations(self):
        """
        Return the row permutation symmetries.

        OUTPUT:

        The row permutations (that is, permutations of projective
        space factors) that can be compensated by a permutation of the
        equations (columns).
        
        EXAMPLES::

            sage: P2xP2.<x,y> = ProductProjectiveSpaces([2, 2], QQ)
            sage: X = P2xP2.complete_intersection(3, 3)
            sage: X.projective_space_permutations()
            Subgroup of (Symmetric group of order 2! as a permutation group)
            generated by [(1,2)]

            sage: four_P2 = ProductProjectiveSpaces([2] * 4, QQ)
            sage: X = four_P2.complete_intersection(
            ....:     [1,1,1,0,0], [0,0,1,1,1], [1,1,1,0,0], [0,0,1,1,1])
            sage: X.projective_space_permutations()
            Subgroup of (Symmetric group of order 4! as a permutation group)
            generated by [(1,2,3,4), (1,3)]
        """
        ambient = self.ambient_space()
        n = ambient.num_factors()
        # Necessary condition: A permutation must preserve the ambient
        # space factors, total degrees, and degrees up to permutation
        projective_spaces = tuple(
            (ambient.factor_dim(i), sum(self.degree(i)), tuple(sorted(self.degree(i))))
            for i in range(n))
        # Sufficient condition: the action on the columns (equations)
        # preserves the set of equations.
        equations = frozenset(
            self.equation(i) for i in range(self.num_equations()))
        def equation_action(permutation):
            return frozenset(permutation(eq) for eq in equations)
        # Find all permutations
        result = []
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        G = SymmetricGroup(n)
        for p in SymmetricGroup(n):
            if p.is_one():
                continue
            if p(projective_spaces) != projective_spaces:
                continue    # necessary condition
            if equation_action(p) != equations:
                continue    # sufficient condititon
            result.append(p)
        H = G.subgroup(result)
        return G.subgroup(H.gens_small())
        
                
