r"""
Points on products of projective spaces.
"""

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.misc.cachefunc import cached_method
from sage.schemes.generic.morphism import SchemeMorphism, SchemeMorphism_point


#*******************************************************************
# Projective varieties
#*******************************************************************
class SchemeMorphism_point_product_projective_spaces_ring(SchemeMorphism_point):

    def __init__(self, X, v, check=True):
        """
        A rational point of a product of projective spaces over a ring.
    
        INPUT:
    
        - ``X`` -- a homset of a subscheme of an ambient product of
           projective space over a field `K`
    
        - ``v`` -- a list or tuple of coordinates in `K`
    
        - ``check`` -- boolean (optional, default:``True``). Whether to
          check the input for consistency.
    
        EXAMPLES::
    
            sage: P = ProductProjectiveSpaces([1,2], ZZ)
            sage: P(1,2, 3,4,5)
            (1 : 2 | 3 : 4 : 5)
        """
        SchemeMorphism.__init__(self, X)
        if check:
            PP = X.codomain()
            if len(v) != PP.ngens():
                raise ValueError('coordinate list length mismatch')
            PP._check_satisfies_equations(v)
        self._coords = tuple(v)

    @cached_method
    def factor(self):
        """
        Return the point restricted to all factor projective spaces.

        OUTPUT:

        Tuple of projective points.

        EXAMPLES::

            sage: P1xP2 = ProductProjectiveSpaces([1,2], ZZ)
            sage: p = P1xP2(1,2, 3,4,5);  p
            (1 : 2 | 3 : 4 : 5)
            sage: p.factor()
            ((1 : 2), (3 : 4 : 5))
        """
        PP = self.codomain()
        return tuple(P(v) for P, v in zip(PP.factor(), PP._iter_factors(self._coords)))
        
    def __eq__(left, right):
        """
        Tests the projective equality of two points.

        INPUT:

        - ``right`` - a point on projective space

        OUTPUT:

        Boolean. ``True`` if ``self`` and ``right`` define the same
        point. ``False`` otherwise.

        Examples::

            sage: PS = ProductProjectiveSpaces([2, 1], ZZ)
            sage: P = PS([1,2,3, 1,2])
            sage: Q = PS([1,2,3, 2,4])
            sage: P == Q
            True
        """
        if not isinstance(right, SchemeMorphism_point):
            try:
                right = left.codomain()(right)
            except TypeError:
                return False
        return all(p.__eq__(q) for p, q in zip(left.factor(), right.factor()))

    def __ne__(left, right):
        """
        Tests the projective equality of two points.

        INPUT:

        - ``right`` - a point on projective space

        OUTPUT:

        Boolean. ``True`` if ``self`` and ``right`` define the same
        point. ``False`` otherwise.

        Examples::

            sage: PS = ProductProjectiveSpaces([2, 1], ZZ)
            sage: P = PS([1,2,3, 1,2])
            sage: Q = PS([1,2,3, 2,4])
            sage: P != Q
            False
        """
        if not isinstance(right, SchemeMorphism_point):
            try:
                right = left.codomain()(right)
            except TypeError:
                return True
        return any(p.__ne__(q) for p, q in zip(left.factor(), right.factor()))

