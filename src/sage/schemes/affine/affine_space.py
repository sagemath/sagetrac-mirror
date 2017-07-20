"""
Affine `n` space over a ring
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from six import integer_types

from sage.functions.orthogonal_polys import chebyshev_T, chebyshev_U
from sage.rings.all import (PolynomialRing, ZZ, Integer)
from sage.rings.ring import is_Ring
from sage.rings.rational_field import is_RationalField
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.categories.map import Map
from sage.categories.fields import Fields
_Fields = Fields()
from sage.categories.homset import End
from sage.categories.number_fields import NumberFields
from sage.misc.all import latex
from sage.structure.category_object import normalize_names
from sage.schemes.generic.scheme import AffineScheme
from sage.schemes.generic.ambient_space import AmbientSpace
from sage.schemes.affine.affine_homset import SchemeHomset_points_affine
from sage.schemes.affine.affine_morphism import (SchemeMorphism_polynomial_affine_space,
                                                 SchemeMorphism_polynomial_affine_space_field,
                                                 SchemeMorphism_polynomial_affine_space_finite_field)
from sage.schemes.affine.affine_point import (SchemeMorphism_point_affine,
                                              SchemeMorphism_point_affine_field,
                                              SchemeMorphism_point_affine_finite_field)



def is_AffineSpace(x):
    r"""
    Returns True if ``x`` is an affine space.

    EXAMPLES::

        sage: from sage.schemes.affine.affine_space import is_AffineSpace
        sage: is_AffineSpace(AffineSpace(5, names='x'))
        True
        sage: is_AffineSpace(AffineSpace(5, GF(9, 'alpha'), names='x'))
        True
        sage: is_AffineSpace(Spec(ZZ))
        False
    """
    return isinstance(x, AffineSpace_generic)

def AffineSpace(n, R=None, names='x'):
    r"""
    Return affine space of dimension ``n`` over the ring ``R``.

    EXAMPLES:

    The dimension and ring can be given in either order::

        sage: AffineSpace(3, QQ, 'x')
        Affine Space of dimension 3 over Rational Field
        sage: AffineSpace(5, QQ, 'x')
        Affine Space of dimension 5 over Rational Field
        sage: A = AffineSpace(2, QQ, names='XY'); A
        Affine Space of dimension 2 over Rational Field
        sage: A.coordinate_ring()
        Multivariate Polynomial Ring in X, Y over Rational Field

    Use the divide operator for base extension::

        sage: AffineSpace(5, names='x')/GF(17)
        Affine Space of dimension 5 over Finite Field of size 17

    The default base ring is `\ZZ`::

        sage: AffineSpace(5, names='x')
        Affine Space of dimension 5 over Integer Ring

    There is also an affine space associated to each polynomial ring::

        sage: R = GF(7)['x, y, z']
        sage: A = AffineSpace(R); A
        Affine Space of dimension 3 over Finite Field of size 7
        sage: A.coordinate_ring() is R
        True
    """
    if is_MPolynomialRing(n) and R is None:
        R = n
        A = AffineSpace(R.ngens(), R.base_ring(), R.variable_names())
        A._coordinate_ring = R
        return A
    if isinstance(R, integer_types + (Integer,)):
        n, R = R, n
    if R is None:
        R = ZZ  # default is the integers
    if names is None:
        if n == 0:
            names = ''
        else:
            raise TypeError("you must specify the variables names of the coordinate ring")
    names = normalize_names(n, names)
    if R in _Fields:
        if is_FiniteField(R):
            return AffineSpace_finite_field(n, R, names)
        else:
            return AffineSpace_field(n, R, names)
    return AffineSpace_generic(n, R, names)


class AffineSpace_generic(AmbientSpace, AffineScheme):
    """
    Affine space of dimension `n` over the ring `R`.

    EXAMPLES::

        sage: X.<x,y,z> = AffineSpace(3, QQ)
        sage: X.base_scheme()
        Spectrum of Rational Field
        sage: X.base_ring()
        Rational Field
        sage: X.category()
        Category of schemes over Rational Field
        sage: X.structure_morphism()
        Scheme morphism:
          From: Affine Space of dimension 3 over Rational Field
          To:   Spectrum of Rational Field
          Defn: Structure map

    Loading and saving::

        sage: loads(X.dumps()) == X
        True

    We create several other examples of affine spaces::

        sage: AffineSpace(5, PolynomialRing(QQ, 'z'), 'Z')
        Affine Space of dimension 5 over Univariate Polynomial Ring in z over Rational Field

        sage: AffineSpace(RealField(), 3, 'Z')
        Affine Space of dimension 3 over Real Field with 53 bits of precision

        sage: AffineSpace(Qp(7), 2, 'x')
        Affine Space of dimension 2 over 7-adic Field with capped relative precision 20

    Even 0-dimensional affine spaces are supported::

        sage: AffineSpace(0)
        Affine Space of dimension 0 over Integer Ring
    """
    def __init__(self, n, R, names):
        """
        EXAMPLES::

            sage: AffineSpace(3, Zp(5), 'y')
            Affine Space of dimension 3 over 5-adic Ring with capped relative precision 20
        """
        AmbientSpace.__init__(self, n, R)
        self._assign_names(names)
        AffineScheme.__init__(self, self.coordinate_ring(), R)

    def __iter__(self):
        """
        Return iterator over the elements of this affine space when defined over a finite field.

        EXAMPLES::

            sage: FF = FiniteField(3)
            sage: AA = AffineSpace(FF, 0)
            sage: [ x for x in AA ]
            [()]
            sage: AA = AffineSpace(FF, 1, 'Z')
            sage: [ x for x in AA ]
            [(0), (1), (2)]
            sage: AA.<z,w> = AffineSpace(FF, 2)
            sage: [ x for x in AA ]
            [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2)]

        AUTHOR:

        - David Kohel
        """
        n = self.dimension_relative()
        R = self.base_ring()
        zero = R(0)
        P = [ zero for _ in range(n) ]
        yield self(P)
        iters = [ iter(R) for _ in range(n) ]
        for x in iters: next(x) # put at zero
        i = 0
        while i < n:
            try:
                P[i] = next(iters[i])
                yield self(P)
                i = 0
            except StopIteration:
                iters[i] = iter(R)  # reset
                next(iters[i]) # put at zero
                P[i] = zero
                i += 1

    def ngens(self):
        """
        Return the number of generators of self, i.e. the number of
        variables in the coordinate ring of self.

        EXAMPLES::

            sage: AffineSpace(3, QQ).ngens()
            3
            sage: AffineSpace(7, ZZ).ngens()
            7
        """
        return self.dimension_relative()

    def rational_points(self, F=None):
        """
        Return the list of ``F``-rational points on the affine space self,
        where ``F`` is a given finite field, or the base ring of self.

        EXAMPLES::

            sage: A = AffineSpace(1, GF(3))
            sage: A.rational_points()
            [(0), (1), (2)]
            sage: A.rational_points(GF(3^2, 'b'))
            [(0), (b), (b + 1), (2*b + 1), (2), (2*b), (2*b + 2), (b + 2), (1)]

            sage: AffineSpace(2, ZZ).rational_points(GF(2))
            [(0, 0), (1, 0), (0, 1), (1, 1)]

        TESTS::

            sage: AffineSpace(2, QQ).rational_points()
            Traceback (most recent call last):
            ...
            TypeError: base ring (= Rational Field) must be a finite field
            sage: AffineSpace(1, GF(3)).rational_points(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: second argument (= Integer Ring) must be a finite field
        """
        if F is None:
            if not is_FiniteField(self.base_ring()):
                raise TypeError("base ring (= %s) must be a finite field"%self.base_ring())
            return [ P for P in self ]
        elif not is_FiniteField(F):
            raise TypeError("second argument (= %s) must be a finite field"%F)
        return [ P for P in self.base_extend(F) ]

    def __eq__(self, right):
        """
        Compare the space with ``right``.

        EXAMPLES::

            sage: AffineSpace(QQ, 3, 'a') == AffineSpace(ZZ, 3, 'a')
            False
            sage: AffineSpace(ZZ, 1, 'a') == AffineSpace(ZZ, 0, 'a')
            False
            sage: A = AffineSpace(ZZ, 1, 'x')
            sage: loads(A.dumps()) == A
            True
        """
        if not isinstance(right, AffineSpace_generic):
            return False
        return (self.dimension_relative() == right.dimension_relative() and
                self.coordinate_ring() == right.coordinate_ring())

    def __ne__(self, other):
        """
        Check whether the space is not equal to ``other``.

        EXAMPLES::

            sage: AffineSpace(QQ, 3, 'a') != AffineSpace(ZZ, 3, 'a')
            True
            sage: AffineSpace(ZZ, 1, 'a') != AffineSpace(ZZ, 0, 'a')
            True
        """
        return not (self == other)

    def _latex_(self):
        r"""
        Return a LaTeX representation of this affine space.

        EXAMPLES::

            sage: print(latex(AffineSpace(1, ZZ, 'x')))
            \mathbf{A}_{\Bold{Z}}^1

        TESTS::

            sage: AffineSpace(3, Zp(5), 'y')._latex_()
            '\\mathbf{A}_{\\ZZ_{5}}^3'
        """
        return "\\mathbf{A}_{%s}^%s"%(latex(self.base_ring()), self.dimension_relative())

    def _morphism(self, *args, **kwds):
        """
        Construct a morphism determined by action on points of this affine space.

        INPUT:

        Same as for
        :class:`~sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space`.

        OUTPUT:

        A new instance of
        :class:`~sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space`.

        EXAMPLES::

            sage: AA = AffineSpace(QQ, 3, 'a')
            sage: AA.inject_variables()
            Defining a0, a1, a2
            sage: EndAA = AA.Hom(AA)
            sage: AA._morphism(EndAA, [a0*a1, a1*a2, a0*a2])
            Scheme endomorphism of Affine Space of dimension 3 over Rational Field
              Defn: Defined on coordinates by sending (a0, a1, a2) to
                    (a0*a1, a1*a2, a0*a2)
        """
        return SchemeMorphism_polynomial_affine_space(*args, **kwds)

    def _point_homset(self, *args, **kwds):
        """
        Construct a Hom-set for this affine space.

        INPUT:

        Same as for
        :class:`~sage.schemes.affine.affine_homset.SchemeHomset_points_affine`.

        OUTPUT:

        A new instance of
        :class:`~sage.schemes.affine.affine_homset.SchemeHomset_points_affine`.

        EXAMPLES::

            sage: AA = AffineSpace(QQ, 3, 'a')
            sage: AA._point_homset(Spec(QQ), AA)
            Set of rational points of Affine Space of dimension 3 over Rational Field
        """
        return SchemeHomset_points_affine(*args, **kwds)

    def _point(self, *args, **kwds):
        r"""
        Construct a point of affine space.

        INPUT:

        Same as for
        :class:`~sage.schemes.affine.affine_point.SchemeMorphism_point_affine`.

        OUTPUT:

        A new instance of
        :class:`~sage.schemes.affine.affine_point.SchemeMorphism_point_affine`.

        TESTS::

            sage: AA = AffineSpace(QQ, 3, 'a')
            sage: AA._point(AA.point_homset(), [0, 1, 2])
            (0, 1, 2)
        """
        return SchemeMorphism_point_affine(*args, **kwds)

    def _repr_(self):
        """
        Return a string representation of this affine space.

        EXAMPLES::

            sage: AffineSpace(1, ZZ, 'x')
            Affine Space of dimension 1 over Integer Ring

        TESTS::

            sage: AffineSpace(3, Zp(5), 'y')._repr_()
            'Affine Space of dimension 3 over 5-adic Ring with capped relative precision 20'
        """
        return "Affine Space of dimension %s over %s"%(self.dimension_relative(), self.base_ring())

    def _repr_generic_point(self, polys=None):
        """
        Return a string representation of the generic point
        corresponding to the list of polys on this affine space.

        If polys is None, the representation of the generic point of
        the affine space is returned.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, ZZ)
            sage: A._repr_generic_point([y-x^2])
            '(-x^2 + y)'
            sage: A._repr_generic_point()
            '(x, y)'
        """
        if polys is None:
            polys = self.gens()
        return '(%s)'%(", ".join([str(f) for f in polys]))

    def _latex_generic_point(self, v=None):
        """
        Return a LaTeX representation of the generic point
        corresponding to the list of polys ``v`` on this affine space.

        If ``v`` is None, the representation of the generic point of
        the affine space is returned.

        EXAMPLES::

            sage: A.<x, y> = AffineSpace(2, ZZ)
            sage: A._latex_generic_point([y-x^2])
            '\\left(- x^{2} + y\\right)'
            sage: A._latex_generic_point()
            '\\left(x, y\\right)'
        """
        if v is None:
            v = self.gens()
        return '\\left(%s\\right)'%(", ".join([str(latex(f)) for f in v]))

    def _check_satisfies_equations(self, v):
        """
        Return True if ``v`` defines a point on the scheme self; raise a
        TypeError otherwise.

        EXAMPLES::

            sage: A = AffineSpace(3, ZZ)
            sage: A._check_satisfies_equations([1, 1, 0])
            True
            sage: A._check_satisfies_equations((0, 1, 0))
            True
            sage: A._check_satisfies_equations([0, 0, 0])
            True
            sage: A._check_satisfies_equations([1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            TypeError: the list v=[1, 2, 3, 4, 5] must have 3 components
            sage: A._check_satisfies_equations([1/2, 1, 1])
            Traceback (most recent call last):
            ...
            TypeError: the components of v=[1/2, 1, 1] must be elements of Integer Ring
            sage: A._check_satisfies_equations(5)
            Traceback (most recent call last):
            ...
            TypeError: the argument v=5 must be a list or tuple
        """
        if not isinstance(v, (list, tuple)):
            raise TypeError('the argument v=%s must be a list or tuple'%v)
        n = self.ngens()
        if not len(v) == n:
            raise TypeError('the list v=%s must have %s components'%(v, n))
        R = self.base_ring()
        from sage.structure.sequence import Sequence
        if not Sequence(v).universe() == R:
            raise TypeError('the components of v=%s must be elements of %s'%(v, R))
        return True

    def __pow__(self, m):
        """
        Return the Cartesian power of this space.

        INPUT:

        - ``m`` -- integer.

        OUTPUT:

        - affine ambient space.

        EXAMPLES::

            sage: A = AffineSpace(1, QQ, 'x')
            sage: A5 = A^5; A5
            Affine Space of dimension 5 over Rational Field
            sage: A5.variable_names()
            ('x0', 'x1', 'x2', 'x3', 'x4')
            sage: A2 = AffineSpace(2, QQ, "x, y")
            sage: A4 = A2^2; A4
            Affine Space of dimension 4 over Rational Field
            sage: A4.variable_names()
            ('x0', 'x1', 'x2', 'x3')

        As you see, custom variable names are not preserved by power operator,
        since there is no natural way to make new ones in general.
        """
        mm = int(m)
        if mm != m:
            raise ValueError("m must be an integer")
        return AffineSpace(self.dimension_relative() * mm, self.base_ring())

    def __mul__(self, right):
        r"""
        Create the product of affine spaces.

        INPUT:

        - ``right`` - an affine space or subscheme.

        OUTPUT: an affine space.= or subscheme.

        EXAMPLES::

            sage: A1 = AffineSpace(QQ, 1, 'x')
            sage: A2 = AffineSpace(QQ, 2, 'y')
            sage: A3 = A1*A2; A3
            Affine Space of dimension 3 over Rational Field
            sage: A3.variable_names()
            ('x', 'y0', 'y1')

            ::

            sage: A2 = AffineSpace(ZZ, 2, 't')
            sage: A3 = AffineSpace(ZZ, 3, 'x')
            sage: A3.inject_variables()
            Defining x0, x1, x2
            sage: X = A3.subscheme([x0*x2 - x1])
            sage: A2*X
            Closed subscheme of Affine Space of dimension 5 over Integer Ring defined by:
              x0*x2 - x1

        ::

            sage: S = ProjectiveSpace(QQ, 3, 'x')
            sage: T = AffineSpace(2, QQ, 'y')
            sage: T*S
            Traceback (most recent call last):
            ...
            TypeError: Projective Space of dimension 3 over Rational Field
            must be an affine space or affine subscheme
        """
        if self.base_ring() != right.base_ring():
            raise ValueError ('Must have the same base ring')

        from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme

        if isinstance(right, AffineSpace_generic):
            if self is right:
                return self.__pow__(2)
            return AffineSpace(self.dimension_relative() + right.dimension_relative(),\
                    self.base_ring(), self.variable_names() + right.variable_names())
        elif isinstance(right, AlgebraicScheme_subscheme):
            AS = self*right.ambient_space()
            CR = AS.coordinate_ring()
            n = self.ambient_space().coordinate_ring().ngens()

            phi = self.ambient_space().coordinate_ring().hom(list(CR.gens()[:n]), CR)
            psi = right.ambient_space().coordinate_ring().hom(list(CR.gens()[n:]), CR)
            return AS.subscheme([phi(t) for t in self.defining_polynomials()] + [psi(t) for t in right.defining_polynomials()])
        else:
            raise TypeError('%s must be an affine space or affine subscheme'%right)

    def change_ring(self, R):
        r"""
        Return an affine space over ring ``R`` and otherwise the same as this space.

        INPUT:

        - ``R`` -- commutative ring or morphism.

        OUTPUT:

        - affine space over ``R``.

        .. NOTE::

            There is no need to have any relation between `R` and the base ring
            of  this space, if you want to have such a relation, use
            ``self.base_extend(R)`` instead.

        EXAMPLES::

            sage: A.<x,y,z> = AffineSpace(3, ZZ)
            sage: AQ = A.change_ring(QQ); AQ
            Affine Space of dimension 3 over Rational Field
            sage: AQ.change_ring(GF(5))
            Affine Space of dimension 3 over Finite Field of size 5

        ::

            sage: K.<w> = QuadraticField(5)
            sage: A = AffineSpace(K,2,'t')
            sage: A.change_ring(K.embeddings(CC)[1])
            Affine Space of dimension 2 over Complex Field with 53 bits of precision
        """
        if isinstance(R, Map):
            return AffineSpace(self.dimension_relative(), R.codomain(), self.variable_names())
        else:
            return AffineSpace(self.dimension_relative(), R, self.variable_names())

    def coordinate_ring(self):
        """
        Return the coordinate ring of this scheme, if defined.

        EXAMPLES::

            sage: R = AffineSpace(2, GF(9,'alpha'), 'z').coordinate_ring(); R
            Multivariate Polynomial Ring in z0, z1 over Finite Field in alpha of size 3^2
            sage: AffineSpace(3, R, 'x').coordinate_ring()
            Multivariate Polynomial Ring in x0, x1, x2 over Multivariate Polynomial Ring
            in z0, z1 over Finite Field in alpha of size 3^2
        """
        try:
            return self._coordinate_ring
        except AttributeError:
            self._coordinate_ring = PolynomialRing(self.base_ring(),
            self.dimension_relative(), names=self.variable_names())
            return self._coordinate_ring

    def _validate(self, polynomials):
        """
        If ``polynomials`` is a tuple of valid polynomial functions on the affine space,
        return ``polynomials``, otherwise raise TypeError.

        Since this is an affine space, all polynomials are valid.

        INPUT:

        - ``polynomials`` -- tuple of polynomials in the coordinate ring of
          this space.

        OUTPUT:

        - tuple of polynomials in the coordinate ring of this space.

        EXAMPLES::

            sage: A.<x, y, z> = AffineSpace(3, ZZ)
            sage: A._validate((x*y - z, 1))
            (x*y - z, 1)
        """
        return polynomials

    def projective_embedding(self, i=None, PP=None):
        """
        Returns a morphism from this space into an ambient projective space
        of the same dimension.

        INPUT:


        -  ``i`` -- integer (default: dimension of self = last
           coordinate) determines which projective embedding to compute. The
           embedding is that which has a 1 in the i-th coordinate, numbered
           from 0.

        -  ``PP`` -- (default: None) ambient projective space, i.e.,
           codomain of morphism; this is constructed if it is not
           given.

        EXAMPLES::

            sage: AA = AffineSpace(2, QQ, 'x')
            sage: pi = AA.projective_embedding(0); pi
            Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1) to
                    (1 : x0 : x1)
            sage: z = AA(3, 4)
            sage: pi(z)
            (1/4 : 3/4 : 1)
            sage: pi(AA(0,2))
            (1/2 : 0 : 1)
            sage: pi = AA.projective_embedding(1); pi
            Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x0, x1) to
                    (x0 : 1 : x1)
            sage: pi(z)
            (3/4 : 1/4 : 1)
            sage: pi = AA.projective_embedding(2)
            sage: pi(z)
            (3 : 4 : 1)

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: A.projective_embedding(2).codomain().affine_patch(2) == A
            True
        """
        n = self.dimension_relative()
        if i is None:
            try:
                i = self._default_embedding_index
            except AttributeError:
                i = int(n)
        else:
            i = int(i)

        try:
            phi = self.__projective_embedding[i]
            #assume that if you've passed in a new codomain you want to override
            #the existing embedding
            if PP is None or phi.codomain() == PP:
                return(phi)
        except AttributeError:
            self.__projective_embedding = {}
        except KeyError:
            pass

        #if no i-th embedding exists, we may still be here with PP==None
        if PP is None:
            from sage.schemes.projective.projective_space import ProjectiveSpace
            PP = ProjectiveSpace(n, self.base_ring())
        elif PP.dimension_relative() != n:
            raise ValueError("projective Space must be of dimension %s"%(n))

        R = self.coordinate_ring()
        v = list(R.gens())
        if n < 0 or n >self.dimension_relative():
            raise ValueError("argument i (=%s) must be between 0 and %s, inclusive"%(i,n))
        v.insert(i, R(1))
        phi = self.hom(v, PP)
        self.__projective_embedding[i] = phi
        #make affine patch and projective embedding match
        PP.affine_patch(i,self)
        return phi

    def subscheme(self, X):
        """
        Return the closed subscheme defined by ``X``.

        INPUT:


        -  ``X`` - a list or tuple of equations.


        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: X = A.subscheme([x, y^2, x*y^2]); X
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x,
              y^2,
              x*y^2

        ::

            sage: X.defining_polynomials ()
            (x, y^2, x*y^2)
            sage: I = X.defining_ideal(); I
            Ideal (x, y^2, x*y^2) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: I.groebner_basis()
            [y^2, x]
            sage: X.dimension()
            0
            sage: X.base_ring()
            Rational Field
            sage: X.base_scheme()
            Spectrum of Rational Field
            sage: X.structure_morphism()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x,
              y^2,
              x*y^2
              To:   Spectrum of Rational Field
              Defn: Structure map
            sage: X.dimension()
            0
        """
        from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_affine
        return AlgebraicScheme_subscheme_affine(self, X)

    def _an_element_(self):
        r"""
        Returns an element of this affine space,used both for illustration and testing purposes.

        OUTPUT: A point in the affine space.

        EXAMPLES::

            sage: AffineSpace(ZZ, 2, 'x').an_element()
            (5, 4)

            sage: AffineSpace(Qp(5), 2, 'x').an_element()
            (5^2 + O(5^22), 4*5 + O(5^21))
        """
        n = self.dimension_relative()
        R = self.base_ring()
        return self([(5 - i) * R.an_element() for i in range(n)])

    def chebyshev_polynomial(self, n, kind='first'):
        """
        Generates an endomorphism of this affine line by a Chebyshev polynomial.
        
        Chebyshev polynomials are a sequence of recursively defined orthogonal
        polynomials. Chebyshev of the first kind are defined as `T_0(x) = 1`,
        `T_1(x) = x`, and `T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x)`. Chebyshev of
        the second kind are defined as `U_0(x) = 1`,
        `U_1(x) = 2x`, and `U_{n+1}(x) = 2xU_n(x) - U_{n-1}(x)`.
    
        INPUT:
    
        - ``n`` -- a non-negative integer.
    
        - ``kind`` -- ``first`` or ``second`` specifying which kind of chebyshev the user would like
          to generate. Defaults to ``first``.
    
        OUTPUT: :class:`SchemeMorphism_polynomial_affine_space`
    
        EXAMPLES::
    
            sage: A.<x> = AffineSpace(QQ, 1)
            sage: A.chebyshev_polynomial(5, 'first')
            Scheme endomorphism of Affine Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x) to
            (16*x^5 - 20*x^3 + 5*x)
    
        ::
    
            sage: A.<x> = AffineSpace(QQ, 1)
            sage: A.chebyshev_polynomial(3, 'second')
            Scheme endomorphism of Affine Space of dimension 1 over Rational Field
            Defn: Defined on coordinates by sending (x) to
            (8*x^3 - 4*x)
    
        ::
    
            sage: A.<x> = AffineSpace(QQ, 1)
            sage: A.chebyshev_polynomial(3, 2)
            Traceback (most recent call last):
            ...
            ValueError: keyword 'kind' must have a value of either 'first' or 'second'
    
        ::
    
            sage: A.<x> = AffineSpace(QQ, 1)
            sage: A.chebyshev_polynomial(-4, 'second')
            Traceback (most recent call last):
            ...
            ValueError: first parameter 'n' must be a non-negative integer

        ::

            sage: A = AffineSpace(QQ, 2, 'x')
            sage: A.chebyshev_polynomial(2)
            Traceback (most recent call last):
            ...
            TypeError: affine space must be of dimension 1
        """
        if self.dimension_relative() != 1:
            raise TypeError("affine space must be of dimension 1")
        n = ZZ(n)
        if (n < 0):
            raise ValueError("first parameter 'n' must be a non-negative integer")
        H = End(self)
        if kind == 'first':
            return H([chebyshev_T(n, self.gen(0))])
        elif kind == 'second':
            return H([chebyshev_U(n, self.gen(0))])
        else:
            raise ValueError("keyword 'kind' must have a value of either 'first' or 'second'")

class AffineSpace_field(AffineSpace_generic):
    def _point(self, *args, **kwds):
        """
        Construct a point.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = AffineSpace(3, GF(3))
            sage: point_homset = P2._point_homset(Spec(GF(3)), P2)
            sage: P2._point(point_homset, [1, 2, 3])
            (1, 2, 0)
        """
        return SchemeMorphism_point_affine_field(*args, **kwds)

    def _morphism(self, *args, **kwds):
        """
        Construct a morphism.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = AffineSpace(3, GF(3))
            sage: P2._morphism(P2.Hom(P2), [x, y, z])
            Scheme endomorphism of Affine Space of dimension 3 over Finite Field of size 3
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x, y, z)
        """
        return SchemeMorphism_polynomial_affine_space_field(*args, **kwds)

    def points_of_bounded_height(self,bound):
        r"""
        Returns an iterator of the points in this affine space of
        absolute height of at most the given bound.

        Bound check  is strict for the rational field.
        Requires this space to be affine space over a number field. Uses the
        Doyle-Krumm algorithm for computing algebraic numbers up
        to a given height [Doyle-Krumm]_.

        INPUT:

        - ``bound`` - a real number.

        OUTPUT:

        - an iterator of points in self.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: list(A.points_of_bounded_height(3))
            [(0, 0), (1, 0), (-1, 0), (1/2, 0), (-1/2, 0), (2, 0), (-2, 0), (0, 1),
            (1, 1), (-1, 1), (1/2, 1), (-1/2, 1), (2, 1), (-2, 1), (0, -1), (1, -1),
            (-1, -1), (1/2, -1), (-1/2, -1), (2, -1), (-2, -1), (0, 1/2), (1, 1/2),
            (-1, 1/2), (1/2, 1/2), (-1/2, 1/2), (2, 1/2), (-2, 1/2), (0, -1/2), (1, -1/2),
            (-1, -1/2), (1/2, -1/2), (-1/2, -1/2), (2, -1/2), (-2, -1/2), (0, 2), (1, 2),
            (-1, 2), (1/2, 2), (-1/2, 2), (2, 2), (-2, 2), (0, -2), (1, -2), (-1, -2), (1/2, -2),
            (-1/2, -2), (2, -2), (-2, -2)]

        ::

            sage: u = QQ['u'].0
            sage: A.<x,y> = AffineSpace(NumberField(u^2 - 2, 'v'), 2)
            sage: len(list(A.points_of_bounded_height(6)))
            121
        """
        if (is_RationalField(self.base_ring())):
            ftype = False # stores whether field is a number field or the rational field
        elif (self.base_ring() in NumberFields()): # true for rational field as well, so check is_RationalField first
            ftype = True
        else:
            raise NotImplementedError("self must be affine space over a number field.")

        bound = bound**(1/self.base_ring().absolute_degree()) # convert to relative height

        n = self.dimension_relative()
        R = self.base_ring()
        zero = R(0)
        P = [ zero for _ in range(n) ]
        yield self(P)
        if (ftype == False):
            iters = [ R.range_by_height(bound) for _ in range(n) ]
        else:
            iters = [ R.elements_of_bounded_height(bound) for _ in range(n) ]
        for x in iters: next(x) # put at zero
        i = 0
        while i < n:
            try:
                P[i] = next(iters[i])
                yield self(P)
                i = 0
            except StopIteration:
                if (ftype == False):
                    iters[i] = R.range_by_height(bound) # reset
                else:
                    iters[i] = R.elements_of_bounded_height(bound)
                next(iters[i]) # put at zero
                P[i] = zero
                i += 1

    def weil_restriction(self):
        r"""
        Compute the Weil restriction of this affine space over some extension
        field.

        If the field is a finite field, then this computes
        the Weil restriction to the prime subfield.

        OUTPUT: Affine space of dimension ``d * self.dimension_relative()``
                over the base field of ``self.base_ring()``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<w> = NumberField(x^5-2)
            sage: AK.<x,y> = AffineSpace(K, 2)
            sage: AK.weil_restriction()
            Affine Space of dimension 10 over Rational Field
            sage: R.<x> = K[]
            sage: L.<v> = K.extension(x^2+1)
            sage: AL.<x,y> = AffineSpace(L, 2)
            sage: AL.weil_restriction()
            Affine Space of dimension 4 over Number Field in w with defining
            polynomial x^5 - 2
        """
        try:
            X = self.__weil_restriction
        except AttributeError:
            L = self.base_ring()
            if L.is_finite():
                d = L.degree()
                K = L.prime_subfield()
            else:
                d = L.relative_degree()
                K = L.base_field()

            if d == 1:
                X = self
            else:
                X = AffineSpace(K, d*self.dimension_relative(), 'z')
            self.__weil_restriction = X
        return X

    def curve(self,F):
        r"""
        Return a curve defined by ``F`` in this affine space.

        INPUT:

        - ``F`` -- a polynomial, or a list or tuple of polynomials in the coorinate ring
          of this affine space.

        EXAMPLES::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: A.curve([y - x^4, z - y^5])
            Affine Curve over Rational Field defined by -x^4 + y, -y^5 + z
        """
        from sage.schemes.curves.constructor import Curve
        return Curve(F, self)

class AffineSpace_finite_field(AffineSpace_field):
    def _point(self, *args, **kwds):
        """
        Construct a point.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = AffineSpace(3, GF(3))
            sage: point_homset = P2._point_homset(Spec(GF(3)), P2)
            sage: P2._point(point_homset, [1, 2, 3])
            (1, 2, 0)
        """
        return SchemeMorphism_point_affine_finite_field(*args, **kwds)

    def _morphism(self, *args, **kwds):
        """
        Construct a morphism.

        For internal use only. See :mod:`morphism` for details.

        TESTS::

            sage: P2.<x,y,z> = AffineSpace(3, GF(3))
            sage: P2._morphism(P2.Hom(P2), [x, y, z])
            Scheme endomorphism of Affine Space of dimension 3 over Finite Field of size 3
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x, y, z)
        """
        return SchemeMorphism_polynomial_affine_space_finite_field(*args, **kwds)

#fix the pickles from moving affine_space.py
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.schemes.generic.affine_space',
                           'AffineSpace_generic',
                           AffineSpace_generic)
