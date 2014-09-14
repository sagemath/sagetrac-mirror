# -*- coding: utf-8 -*-
r"""
Cohomology Ring of Product of Projective Spaces

To work with the cohomology of a product of projective spaces you
should always use the
:meth:`~sage.schemes.product_porjective.space.ProductProjectiveSpaces_ring.cohomology_class`
method. Do not create instances of the class in this module directly.

EXAMPLES::

    sage: P1xP2.<x0,x1,y0,y1,y2> = ProductProjectiveSpaces([1, 2], QQ)
    sage: HH = P1xP2.cohomology_ring()
    sage: (1 + HH(y0)) * (1 + HH(y1))
    [y2^2 + 2*y2 + 1]
"""


#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************



from sage.misc.all import cached_method, latex
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.all import PolynomialRing, QQ
from sage.rings.quotient_ring_element import QuotientRingElement
from sage.rings.quotient_ring import QuotientRing_generic


class CohomologyClass(QuotientRingElement):

    def __init__(self, cohomology_ring, representative):
        r"""
        Construct the cohomology class.

        INPUT:

        - ``cohomology_ring`` -- :class:`CohomologyRing`.

        - ``representative`` -- a polynomial in the generators of the
          cohomology ring.

        OUTPUT:

        An instance of :class:`CohomologyClass`.

        EXAMPLES::

            sage: P1xP2.<x0,x1,y0,y1,y2> = ProductProjectiveSpaces([1, 2], QQ)
            sage: H = P1xP2.cohomology_ring()
            sage: from sage.schemes.toric.variety import CohomologyClass
            sage: CohomologyClass(H, H.defining_ideal().ring().zero() )
            [0]
        """
        assert representative in cohomology_ring.defining_ideal().ring(), \
            'The given representative is not in the parent polynomial ring.'
        super(CohomologyClass, self).__init__(cohomology_ring, representative)

    def _repr_(self):
        r"""
        Return a string representation of the cohomology class.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: P1xP2.<x0,x1,y0,y1,y2> = ProductProjectiveSpaces([1, 2], QQ)
            sage: P1xP2.cohomology_ring().gen(0)._repr_()
            '[x1]'
        """
        return '['+super(CohomologyClass,self)._repr_()+']'

    def _latex_(self):
        r"""
        Return a latex representation of the cohomology class.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: P1xP2.<x0,x1,y0,y1,y2> = ProductProjectiveSpaces([1, 2], QQ)
            sage: cohomology_class = P1xP2.cohomology_ring().gen(0)
            sage: cohomology_class._latex_()
            '\\left[ x_{1} \\right]'
        """
        return r'\left[ %s \right]' % latex(self.lift())

    def deg(self):
        r"""
        The degree of the cohomology class.

        OUTPUT:

        An integer `d` such that the cohomology class is in degree
        `2d`. If the cohomology class is of mixed degree, the highest
        degree is returned.

        EXAMPLES::

            sage: P1xP2.<x0,x1,y0,y1,y2> = ProductProjectiveSpaces([1, 2], QQ)
            sage: P1xP2.cohomology_ring().gen(0).deg()
            1
            sage: P1xP2.cohomology_ring().zero().deg()
            -1
        """
        return self.lift().degree()

    def part_of_degree(self, d):
        r"""
        Project the (mixed-degree) cohomology class to the given degree.

        .. MATH::

            \mathop{pr}\nolimits_d:~ H^\bullet(\mathbbb{P},\QQ) \to H^{2d}(\mathbbb{P},\QQ)

        INPUT:

        - An integer ``d``

        OUTPUT:

        - The degree-``2d`` part of the cohomology class.

        EXAMPLES::

            sage: P1xP1.<t0,t1,y0,y1> = ProductProjectiveSpaces([1, 1], QQ)
            sage: t = P1xP1.cohomology_ring().gen(0)
            sage: y = P1xP1.cohomology_ring().gen(2)
            sage: 3*t+4*t^2*y+y+t*y+t+1
            [t1*y1 + 4*t1 + y1 + 1]
            sage: (3*t+4*t^2*y+y+t*y+t+1).part_of_degree(1)
            [4*t1 + y1]
        """
        Q = self.parent()
        # We iterate over monomials of self.lift()
        p = [x for x in self.lift() if x[1].total_degree() == d]
        if len(p)==0:
            return Q.zero()
        else:
            return Q(sum(x[0]*x[1] for x in p))

    def exp(self):
        """
        Exponentiate ``self``.

        .. NOTE::

            The exponential `\exp(x)` of a rational number `x` is
            usually not rational. Therefore, the cohomology class must
            not have a constant (degree zero) part. The coefficients
            in the Taylor series of `\exp` are rational, so any
            cohomology class without constant term can be
            exponentiated.

        OUTPUT

        The cohomology class `\exp(` ``self`` `)` if the constant part
        vanishes, otherwise a ``ValueError`` is raised.

        EXAMPLES::

            sage: P1xP2.<x0,x1,y0,y1,y2> = ProductProjectiveSpaces([1, 2], QQ)
            sage: H_class = P1xP2.cohomology_ring().gen(2)
            sage: H_class
            [y2]
            sage: H_class.exp()
            [1/2*y2^2 + y2 + 1]
        """
        if not self.part_of_degree(0).is_zero():
            raise ValueError('Must not have a constant part.')
        from sage.functions.other import factorial
        exp_x = self.parent().one()
        for d in range(1,self.parent()._variety.dimension()+1):
            exp_x += self**d / factorial(d)
        return exp_x


class CohomologyRing(QuotientRing_generic, UniqueRepresentation):

    Element = CohomologyClass

    def __init__(self, variety):
        r"""
        The cohomology ring of a product of projective spaces.

        .. NOTE::

            Projective spaces over any field of characteristic 0 are
            treated as if they were varieties over `\CC`.

        INPUT:

        - ``variety`` -- a product of projective spaces, see
          :mod:`sage.schemes.product_projective`.

        EXAMPLES::

            sage: P1xP1.<x,y> = ProductProjectiveSpaces([1, 1], QQ)
            sage: HH = P1xP1.cohomology_ring();  HH
            Rational cohomology ring of a Product of projective
            spaces P^1 x P^1 over Rational Field
            sage: HH.defining_ideal()
            Ideal (x0^2, x0 - x1, y0^2, y0 - y1) of Multivariate Polynomial
            Ring in x0, x1, y0, y1 over Rational Field
        """
        from sage.schemes.product_projective.space import is_ProductProjectiveSpaces
        if not is_ProductProjectiveSpaces(variety):
            raise ValueError('not a product of projective spaces')
        self._variety = variety
        R = PolynomialRing(QQ, variety.variable_names())
        self._polynomial_ring = R
        relations = []
        for gens in self._variety._iter_factors(R.gens()):
            x0 = gens[0]
            relations.append(x0 ** len(gens))
            for xi in gens[1:]:
                relations.append(x0 - xi)
        I = R.ideal(relations)
        super(CohomologyRing, self).__init__(R, I, names=variety.variable_names())

    def _repr_(self):
        r"""
        Return a string representation of the cohomology ring.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: P2xP2.<x,y> = ProductProjectiveSpaces([2, 2], QQ)
            sage: P2xP2.cohomology_ring()._repr_()
            'Rational cohomology ring of a Product of projective
             spaces P^2 x P^2 over Rational Field'
        """
        return 'Rational cohomology ring of a '+self._variety._repr_()

    def _latex_(self):
        r"""
        Return a latex representation of the cohomology ring.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: P2xP2.<x,y> = ProductProjectiveSpaces([2, 2], QQ)
            sage: cohomology_ring = P2xP2.cohomology_ring()
            sage: print cohomology_ring._latex_()
            H^\ast\left({\mathbf P}_{\Bold{Q}}^2 \times 
            {\mathbf P}_{\Bold{Q}}^2,\QQ\right)
        """
        return 'H^\\ast\\left('+self._variety._latex_()+',\QQ\\right)'

    def _element_constructor_(self,x):
        r"""
        Construct a :class:`CohomologyClass`.

        INPUT::

        - ``x`` -- something that defines a cohomology class. Either a
          cohomology class or something that can be converted into a
          polynomial in the homogeneous coordinates.

        OUTPUT:

        The :class:`CohomologyClass` defined by ``x``.

        EXAMPLES::

            sage: P2xP2.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: HH = P2xP2.cohomology_ring()
            sage: HH(x0*x1)
            [x2^2]

        Numbers will be converted to constants in the ring::

            sage: HH(1)
            [1]
            sage: type(HH(1))
            <class 'sage.schemes.product_projective.cohomology_ring.CohomologyRing_with_category.element_class'>
        """
        if isinstance(x, CohomologyClass) and x.parent()==self:
            return x
        if isinstance(x, QuotientRingElement):
            x = x.lift()
        else:
            x = self.cover_ring()(x)
        return self.element_class(self, x)

    # We definitely should not override __call__, but since our
    # superclass QuotientRing_generic does not adhere to the coercion
    # model we cannot either. See
    # http://trac.sagemath.org/sage_trac/ticket/9429
    def __call__(self, x, coerce=True):
        r"""
        Turn ``x`` into a ``CohomologyClass``.

        EXAMPLES::

            sage: P2xP2.<x0,x1,x2,y0,y1,y2> = ProductProjectiveSpaces([2, 2], QQ)
            sage: H = P2xP2.cohomology_ring()
            sage: H(1)
            [1]
            sage: type( H(1) )
            <class 'sage.schemes.product_projective.cohomology_ring.CohomologyRing_with_category.element_class'>
        """
        return self._element_constructor_(x)

    @cached_method
    def gens(self):
        r"""
        Return the generators of the cohomology ring.

        OUTPUT:

        A tuple of generators, one for each projective space factor.

        EXAMPLES::

            sage: P1xP2.<x0,x1,y0,y1,y2> = ProductProjectiveSpaces([1, 2], QQ)
            sage: P1xP2.cohomology_ring().gens()
            ([x1], [y2])
        """
        result = []
        for gens in self._variety._iter_factors(self.cover_ring().gens()):
            result.append(self.element_class(self, gens[0]))
        return tuple(result)

    def gen(self, i):
        r"""
        Return the generators of the cohomology ring.

        INPUT:

        - ``i`` -- integer.

        OUTPUT:

        The generator of the cohomology ring associated to the
        ``i``-th homogeneous coordinate.

        EXAMPLES::

            sage: P1xP2.<x0,x1,y0,y1,y2> = ProductProjectiveSpaces([1, 2], QQ)
            sage: HH = P1xP2.cohomology_ring()
            sage: HH.gen(0)
            [x1]
            sage: HH.gen(1)
            [x1]
            sage: HH.gen(2)
            [y2]
            sage: HH.gen(3)
            [y2]
            sage: HH.gen(4)
            [y2]
        """
        return self.element_class(self, self.cover_ring().gen(i))


#*****************************************************************
def is_CohomologyClass(x):
    r"""
    Check whether ``x`` is a cohomology class of a toric variety.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    ``True`` or ``False`` depending on whether ``x`` is an instance of
    :class:`CohomologyClass`

    EXAMPLES::

        sage: P1xP2.<x0,x1,y0,y1,y2> = ProductProjectiveSpaces([1, 2], QQ)
        sage: HH = P1xP2.cohomology_ring()
        sage: from sage.schemes.product_projective.cohomology_ring import is_CohomologyClass
        sage: is_CohomologyClass(HH.one())
        True
        sage: is_CohomologyClass('z')
        False
    """
    return isinstance(x, CohomologyClass)


