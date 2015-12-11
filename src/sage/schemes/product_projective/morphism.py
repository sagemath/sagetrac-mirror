r"""
Polynomial morphisms for products of projective spaces

This class builds on the projective space class and its point and morphism classes.

EXAMPLES::

    sage: P1xP1.<x,y, u,v> = ProductProjectiveSpaces(QQ,[1,1])
    sage: H = End(P1xP1)
    sage: H([x^2*u, y^2*v, x*v^2, y*u^2])
    Scheme endomorphism of Product of projective spaces P^1 x P^1 over Rational Field
      Defn: Defined by sending (x : y , u : v) to 
            (x^2*u : y^2*v , x*v^2 : y*u^2).
"""
#*****************************************************************************
# Copyright (C) 2014 Ben Hutz <bn4941@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
# http://www.gnu.org/licenses/
#*****************************************************************************
from sage.schemes.generic.morphism import SchemeMorphism_polynomial


class ProductProjectiveSpaces_morphism_ring(SchemeMorphism_polynomial):
    r"""
    The class of morphisms on products of projective spaces.
    The components are projective space morphisms.

    EXAMPLES::

        sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2,1],QQ)
        sage: H = T.Hom(T)
        sage: H([x^2,y^2,z^2,w^2,u^2])
        Scheme endomorphism of Product of projective spaces P^2 x P^1 over Rational Field
          Defn: Defined by sending (x : y : z , w : u) to 
                (x^2 : y^2 : z^2 , w^2 : u^2).
    """

    def __init__(self, parent, polys, check = True):
        r"""
        The Python constructor.

        INPUT:

        - ``parent`` -- Homset

        - ``polys`` -- anything that defines a point in the class

        - ``check`` -- Boolean. Whether or not to perform input checks
          (Default:`` True``)

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2,1],QQ)
            sage: H = T.Hom(T)
            sage: H([x^2*u,y^2*w,z^2*u,w^2,u^2])
            Scheme endomorphism of Product of projective spaces P^2 x P^1 over Rational Field
              Defn: Defined by sending (x : y : z , w : u) to 
                    (x^2*u : y^2*w : z^2*u , w^2 : u^2).

        ::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2,1],QQ)
            sage: H = T.Hom(T)
            sage: H([x^2*u,y^2*w,z^2*u,w^2,u*z])
            Traceback (most recent call last):
            ...
            TypeError: polys (=[x^2*u, y^2*w, z^2*u, w^2, z*u]) must be
            multi-homogeneous of the same degrees (by component)

        ::

            sage: R.<s,t> = PolynomialRing(QQ)
            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1,2],QQ)
            sage: P.<u,v,w,s,t,r> = ProductProjectiveSpaces([3,1],QQ)
            sage: H = Hom(Z,P)
            sage: f = H([a^2,b^2,a^2,a*b,a*x,b*z]); f
            Scheme morphism:
              From: Product of projective spaces P^1 x P^2 over Rational Field
              To:   Product of projective spaces P^3 x P^1 over Rational Field
              Defn: Defined by sending (a : b , x : y : z) to
                    (a^2 : b^2 : a^2 : a*b , a*x : b*z).

        ::

            sage: Z.<a,b,c,x,y,z> = ProductProjectiveSpaces([1,3],QQ)
            sage: P.<u,v,w,s,t,r> = ProductProjectiveSpaces([2,2],QQ)
            sage: H = Hom(Z,P)
            sage: f = H([a^2,b^2,c^2,x^2,y^2,z^2])
            Traceback (most recent call last):
            ...
            TypeError: polys (=[a^2, b^2, c^2, x^2, y^2, z^2]) must be
            multi-homogeneous of the same degrees (by component)
        """
        if check:
            #check multi-homogeneous
            #if self is a subscheme, we may need the lift of the polynomials
            try:
                polys[0].exponents()
            except AttributeError:
                polys = [f.lift() for f in polys]
            target = parent.codomain().ambient_space()
            dom = parent.domain().ambient_space()
            from sage.schemes.product_projective.space import is_ProductProjectiveSpaces
            if is_ProductProjectiveSpaces(target):
                splitpolys = target._factors(polys)
                for m in range(len(splitpolys)):
                    d = dom._degree(splitpolys[m][0])
                    if not all(d == dom._degree(f) for f in splitpolys[m]):
                        raise  TypeError("polys (=%s) must be multi-homogeneous of the same degrees (by component)"%polys)
            else:
                #we are mapping into some other kind of space
                target._validate(polys)

        SchemeMorphism_polynomial.__init__(self, parent, polys, check)

    def __getitem__(self, i):
        r"""
        Return the `i`-th coordinate polynomial.

        INPUT:

        - `i` -- integer

        OUTPUT:

        The (multi)-homomgeneous polynomial that is the `i`-th coordinate.

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2,1],QQ)
            sage: H = T.Hom(T)
            sage: F = H([x^2*u,y^2*w,z^2*u,w^2,u^2])
            sage: F[2]
            z^2*u
        """
        return(self._polys[i])

    def _repr_defn(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:
        
        String.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1,1], QQ)
            sage: H = Hom(P,P)
            sage: f = H([x^2,y^2,z,w])
            sage: f._repr_defn()
            'Defined by sending (x : y , z : w) to \n(x^2 : y^2 , z : w).'
        """
        s  = 'Defined by sending '
        s += self.domain().ambient_space()._repr_generic_point()
        s += ' to \n'
        s += self.codomain().ambient_space()._repr_generic_point(self._polys)
        s += '.'
        return s

    def __call__(self, P, check = True):
        r"""
        Make morphisms of products of projective spaces callable.

        INPUT:

        - ``P`` -- a point in the domain.

        - ``check`` -- Boolean - whether or not to perform the input checks
          on the image point (Default: ``True``)

        OUTPUT:

        The image point in the codomain.

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2,1],QQ)
            sage: H = T.Hom(T)
            sage: F = H([x^2*u,y^2*w,z^2*u,w^2,u^2])
            sage: F(T([2,1,3,0,1]))
            (4/9 : 0 : 1 , 0 : 1)

        ::

            sage: PP.<x,y,z,u,v,w> = ProductProjectiveSpaces(QQ,[1,1,1])
            sage: HP = End(PP)
            sage: f = HP([v*x^2,w*y^2,z^2,u^2,v^2,w^2])
            sage: Q = PP([0,1,1,1,1,1])
            sage: f(Q)
            (0 : 1 , 1 : 1 , 1 : 1)
        """
        A = self.codomain()
        Q = list(P)
        newP = [f(Q) for f in self.defining_polynomials()]
        return(A.point(newP, check))
