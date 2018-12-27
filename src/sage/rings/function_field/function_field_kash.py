# -*- coding: utf-8 -*-
"""
Function Fields implemented with an expect interface to the optional
package kash.

This implementation is incomplete.  Many standard methods have not
yet been implemented.

EXAMPLES::

    sage: F.<x> = FunctionField(QQ, implementation='kash')
    sage: R.<Y> = F[]
    sage: L.<y> = F.extension(Y^2 - x^8 - 1)
    sage: O = L.maximal_order()
    sage: I = O.ideal(x, y-1)
    sage: P = I.place()
    sage: D = P.divisor()
    sage: D.basis_function_space()
    [1]
    sage: (2*D).basis_function_space()
    [1]
    sage: (3*D).basis_function_space()
    [1]
    sage: (4*D).basis_function_space()
    [1, 1/x^4*y + 1/x^4]


    sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = PolynomialRing(K)
    sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
    sage: O = F.maximal_order()
    sage: I = O.ideal(y)
    sage: I.divisor()
    2*Place (x, (1/(x^3 + x^2 + x))*y^2)
     + 2*Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)

    sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^2+Y+x+1/x)
    sage: O = L.maximal_order()
    sage: I = O.ideal(y)
    sage: I.divisor()
    -1*Place (x, x*y)
     + Place (x^2 + 1, x*y)
"""

from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten

from sage.interfaces.kash import kash, KashElement

from sage.structure.richcmp import richcmp
from sage.structure.factorization import Factorization

from sage.rings.rational_field import QQ
from sage.rings.fraction_field import FractionField
from sage.rings.finite_rings.finite_field_base import FiniteField

from .function_field import RationalFunctionField
from .function_field import FunctionFieldElement_rational
from .function_field import FunctionField_polymod
from .function_field import FunctionFieldElement_polymod
from .order import FunctionFieldMaximalOrder
from .order import FunctionFieldMaximalOrderInfinite
from .ideal import FunctionFieldIdeal
from .place import FunctionFieldPlace
from .divisor import FunctionFieldDivisor

class RationalFunctionField_kash(RationalFunctionField):
    """
    Rational function field `K(t)` in one variable, over the rationals,
    implemented using kash.

    EXAMPLES::

        sage: K.<t> = FunctionField(QQ, implementation='kash'); K
        Rational function field in t over Rational Field
        sage: K.gen()
        t
        sage: 1/t + t^3 + 5
        (t^4 + 5*t + 1)/t

    There are various ways to get at the underlying fields and rings
    associated to a rational function field::

        sage: K.<t> = FunctionField(QQ, implementation='kash')
        sage: K.base_field()
        Rational function field in t over Rational Field
        sage: K.field()
        Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: K.constant_field()
        Rational Field
        sage: K.maximal_order()
        Maximal order of Rational function field in t over Rational Field

    We define a morphism::

        sage: K.<t> = FunctionField(QQ, implementation='kash')
        sage: L = FunctionField(QQ, 'tbar') # give variable name as second input
        sage: K.hom(L.gen())
        Function Field morphism:
          From: Rational function field in t over Rational Field
          To:   Rational function field in tbar over Rational Field
          Defn: t |--> tbar
    """

    def __init__(self, constant_field, names, category=None):
        """
        Create a rational function field in one variable over the rational
        field.

        INPUT:

        - ``constant_field`` -- QQ

        - ``names`` -- string or tuple of length 1

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ, implementation='kash'); K
            Rational function field in t over Rational Field
            sage: K.category()
            Category of function fields
            sage: FunctionField(QQ[I], 'alpha', implementation='kash')   # not implemented
            Rational function field in alpha over Number Field in I with defining polynomial x^2 + 1

        Must be over a field::

            sage: FunctionField(ZZ, 't')
            Traceback (most recent call last):
            ...
            TypeError: constant_field must be a field
        """

        RationalFunctionField.__init__(self, constant_field, names, category)

        if constant_field is QQ:
            self.kash_constant_field = kash.RationalField()
        elif isinstance(constant_field, FiniteField):
            assert constant_field.order().is_prime()
            self.kash_constant_field = kash.GaloisField(constant_field.order())
        #elif isinstance(constant_field, NumberField):
        else:
            raise ValueError("The constant field must be either QQ or a finite field.")

        # we seem to need this to avoid getting variable names like '$.1' Kash's output
        self.kash_constant_field.PolynomialAlgebra().AssignNames_(['"x"'])

        self.kash = self.kash_constant_field.RationalFunctionField()

    def extension(self, f, names=None):
        """
        Create an extension `L = K[y]/(f(y))` of the rational function field.

        INPUT:

        - ``f`` -- univariate polynomial over self

        - ``names`` -- string or length-1 tuple

        OUTPUT:

        - a function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: K.extension(y^5 - x^3 - 3*x + x*y)
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        A nonintegral defining polynomial::

            sage: K.<t> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: K.extension(y^3 + (1/t)*y + t^3/(t+1))
            Function field in y defined by y^3 + 1/t*y + t^3/(t + 1)

        The defining polynomial need not be monic or integral::

            sage: K.extension(t*y^3 + (1/t)*y + t^3/(t+1))
            Function field in y defined by t*y^3 + 1/t*y + t^3/(t + 1)
        """
        from . import constructor
        return constructor.FunctionFieldExtension(f, names, implementation='kash')

    @cached_method
    def place_set(self):
        """
        Return the set of all places of the function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ, implementation='kash')
            sage: K.place_set()
            Set of places of Rational function field in t over Rational Field
        """
        from .place import PlaceSet
        return PlaceSet(self)

class FunctionField_polymod_kash(FunctionField_polymod):
    """
    Function fields defined by a univariate polynomial, as an extension of the
    base field, implemented using kash.

    EXAMPLES:

    We make a function field defined by a degree 5 polynomial over the
    rational function field over the rational numbers::

        sage: K.<x> = FunctionField(QQ, implementation='kash')
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L
        Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x

    We next make a function field over the above nontrivial function
    field L::

        sage: S.<z> = L[]
        sage: M.<z> = L.extension(z^2 + y*z + y); M
        Function field in z defined by z^2 + y*z + y
        sage: 1/z
        ((x/(-x^4 - 1))*y^4 - 2*x^2/(-x^4 - 1))*z - 1
        sage: z * (1/z)
        1

    We drill down the tower of function fields::

        sage: M.base_field()
        Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
        sage: M.base_field().base_field()
        Rational function field in x over Rational Field
        sage: M.base_field().base_field().constant_field()
        Rational Field
        sage: M.constant_base_field()
        Rational Field

    The polynomial must be irreducible::

        sage: K.<x>=FunctionField(QQ, implementation='kash')
        sage: R.<y> = K[]
        sage: L.<y>=K.extension(x^2-y^2)
        Traceback (most recent call last):
        ...
        ValueError: The polynomial must be irreducible
    """

    def __init__(self, polynomial, names, category=None):
        """
        Create a function field defined as an extension of another function
        field by adjoining a root of a univariate polynomial.

        INPUT:

        - ``polynomial`` -- univariate polynomial over a function field

        - ``names`` -- tuple of length 1 or string; variable names

        - ``category`` -- category (default: category of function fields)

        EXAMPLES:

        We create an extension of a function field::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: L = K.extension(y^5 - x^3 - 3*x + x*y); L
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        Note the type::

            sage: type(L)
            <class 'sage.rings.function_field.function_field_kash.FunctionField_polymod_kash_with_category'>

        We can set the variable name, which doesn't have to be y::

            sage: L.<w> = K.extension(y^5 - x^3 - 3*x + x*y); L
            Function field in w defined by w^5 + x*w - x^3 - 3*x

        TESTS:

        Test that :trac:`17033` is fixed::

            sage: K.<t> = FunctionField(QQ)
            sage: R.<x> = QQ[]
            sage: M.<z> = K.extension(x^7-x-t)
            sage: M(x)
            z
            sage: M('z')
            z
            sage: M('x')
            Traceback (most recent call last):
            ...
            TypeError: unable to evaluate 'x' in Fraction Field of Univariate
            Polynomial Ring in t over Rational Field
        """

        FunctionField_polymod.__init__(self, polynomial, names, category)

        assert isinstance(polynomial.base_ring(), RationalFunctionField_kash)

        kash_base_field = polynomial.base_ring().kash
        kTy = kash_base_field.PolynomialAlgebra()

        x = kash_base_field.gen(1)
        y = kTy.gen(1)

        # Goofiness because 'subs' doesn't recurse into coefficient
        # rings, so we convert our polynomial to a list of
        # coefficients, map them to Kash, and pass the list to Kash's
        # Element constructor, which builds a polynomial from a list.

        # Also, FunctionFieldElement's aren't callable, and thus the
        # extra call to element() to convert the FunctionFieldElement
        # to a FractionFieldElement

        kash_polynomial = kTy.Element(map(lambda c: c.element()(x), list(polynomial)))

        try:
            self.kash = kash_polynomial.FunctionField()
        except TypeError as ex:
            if 'Polynomial must be irreducible' in ex.args[0]:
                raise ValueError("The polynomial must be irreducible")
            raise ex

        ybar = self.kash.gen(1)

        # Set up reverse map to convert elements back from kash

        self.reverse_map = {x : polynomial.base_ring().gen(0), ybar : self.gen(0)}

    def to_kash(self, polynomial):

        # Same kind of goofiness as above.  When building an element
        # of an algebraic extension, Kash requires the list to be the
        # same length as the degree of the extension.  Also, Kash's
        # Element function doesn't accept constant polynomials - they
        # have to actually be constants.

        x = self.base_ring().kash.gen(1)
        coeffs = list(self(polynomial).element())
        coeffs += [self.base_ring().zero()] * (self.degree() - len(coeffs))
        def coeff_to_kash(c):
            c = c.element()
            # c will be an element of a fraction field
            if c.denominator() == 1 and c.numerator().is_constant():
                return c.numerator().constant_coefficient()
            else:
                return c(x)
        coeffs = map(coeff_to_kash, coeffs)
        return self.kash.Element(coeffs)

    @cached_method
    def place_set(self):
        """
        Return the set of all places of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: L.place_set()
            Set of places of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        from .place import PlaceSet
        return PlaceSet(self)

    @cached_method
    def divisor_group(self):
        """
        Return the group of divisors attached to the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2))
            sage: L.divisor_group()
            Divisor group of Function field in y defined by y^3 + (-x^3 + 1)/(x^3 - 2)
        """
        from .divisor import DivisorGroup
        return DivisorGroup(self)

    @cached_method
    def maximal_order(self):
        """
        Return the maximal order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash');
            sage: R.<t> = PolynomialRing(K);
            sage: F.<y> = K.extension(t^4 + x^12*t^2 + x^18*t + x^21 + x^18);
            sage: O = F.maximal_order()
            sage: O.basis()
            (1, 1/x^4*y, 1/x^9*y^2, 1/x^13*y^3)

            sage: K.<x> = FunctionField(GF(2), implementation='kash');
            sage: R.<t> = PolynomialRing(K);
            sage: F.<y> = K.extension(t^4 + x^12*t^2 + x^18*t + x^21 + x^18);
            sage: O = F.maximal_order()
            sage: O.basis()
            (1, 1/x^4*y, 1/x^11*y^2 + 1/x^2, 1/x^15*y^3 + 1/x^6*y)

        The basis of the maximal order *always* starts with 1. This is assumed
        in some algorithms.
        """
        return FunctionFieldMaximalOrder_kash(self)

    @cached_method
    def maximal_order_infinite(self):
        """
        Return the maximal infinite order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: F.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: L.maximal_order_infinite()
            Maximal infinite order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        return FunctionFieldMaximalOrderInfinite_kash(self)

class FunctionFieldMaximalOrder_kash(FunctionFieldMaximalOrder):
    """
    Base class of kash-implemented maximal orders of function fields.
    """

    def __init__(self, field, category=None):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ, implementation='kash'); K
            Rational function field in t over Rational Field
            sage: R = K.maximal_order(); R
            Maximal order of Rational function field in t over Rational Field
        """

        FunctionFieldMaximalOrder.__init__(self, field, category)
        self.kash = field.kash.MaximalOrderFinite()

    def _element_constructor_(self, f, check=True):
        """
        Construct an element of this order from ``f``.

        INPUT:

        - ``f`` -- element convertible to the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2-x*Y+x^2+1)
            sage: O = L.maximal_order()
            sage: y in O
            True
            sage: 1/y in O
            False
            sage: x in O
            True
            sage: 1/x in O
            False
            sage: L.<y>=K.extension(Y^2+Y+x+1/x)
            sage: O = L.maximal_order()
            sage: 1 in O
            True
            sage: y in O
            False
            sage: x*y in O
            True
            sage: x^2*y in O
            True
        """
        field = self.function_field()

        #if f.parent() is field:
        #    f = f.element()
        #f = self._field._ring(f)
        if check:
            if not kash._contains(self._field.to_kash(f).name(), self.kash.name()):
                raise TypeError("%r is not an element of %r"%(f,self))
        # return f and not field._element_class(self, f) because
        # 1. that's what FunctionFieldMaximalOrder_global does
        # 2. that's what makes "d in O" work (it tests if O(d) == d)
        # 3. but it's probably not right (O(d).parent() should report the order, not the field)
        # return field._element_class(self, f)
        return f

    def polynomial(self):
        """
        Return the defining polynomial of the function field of which this is an order.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.polynomial()
            y^4 + x*y + 4*x + 1
        """
        return self._field.polynomial()

    def basis(self):
        """
        Return the basis of the order as a module over the polynomial ring.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)
            sage: O = L.equation_order()
            sage: O.basis()
            (1, y, y^2, y^3)
        """

        return tuple(self.kash.Basis().sage(self._field.reverse_map))

    def ideal(self, *gens):
        """
        Return the fractional ideal generated by ``gens``.

        INPUT:

        - ``gens`` -- list of generators or an ideal in a ring which
          coerces to this order

        EXAMPLES::

            sage: K.<y> = FunctionField(QQ, implementation='kash')
            sage: O = K.maximal_order()
            sage: O.ideal(y)
            Ideal (y) of Maximal order of Rational function field in y over Rational Field
            sage: O.ideal([y,1/y]) == O.ideal(y,1/y) # multiple generators may be given as a list
            True

        A fractional ideal of a nontrivial extension::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2-4)
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: S = L.equation_order()
            sage: S.ideal(1/y)
            Ideal (1, (-1/(x^3 + 1))*y) of Order in Function field in y defined by y^2 - x^3 - 1
            sage: I2 = S.ideal(x^2-4); I2
            Ideal (x^2 - 4, (x^2 - 4)*y) of Order in Function field in y defined by y^2 - x^3 - 1
            sage: I2 == S.ideal(I)
            True

            sage: K.<x> = FunctionField(GF(7), implementation='kash'); R.<y> = K[]
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2-4)
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: S = L.equation_order()
            sage: S.ideal(1/y)
            Ideal (1, (6/(x^3 + 1))*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: I2 = S.ideal(x^2-4); I2
            Ideal (x^2 + 3, (x^2 + 3)*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: I2 == S.ideal(I)
            True
        """

        return FunctionFieldIdeal_kash(self, gens)

class FunctionFieldMaximalOrderInfinite_kash(FunctionFieldMaximalOrderInfinite):
    """
    Base class of kash-implemented maximal infinite orders of function fields.
    """

    def __init__(self, field, category=None):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ, implementation='kash'); K
            Rational function field in t over Rational Field
            sage: R = K.maximal_order_infinite(); R
            Maximal infinite order of Rational function field in t over Rational Field
        """

        FunctionFieldMaximalOrderInfinite.__init__(self, field, category)
        self.kash = field.kash.MaximalOrderInfinite()

    def _element_constructor_(self, f):
        """
        Make ``f`` an element of this order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x*y)
            sage: 1 in Oinf
            True
            sage: 1/x*y in Oinf
            True
            sage: x*y in Oinf
            False
            sage: 1/x in Oinf
            True
        """
        if not f.parent() is self.function_field():
            f = self.function_field()(f)

        if not kash._contains(self._field.to_kash(f).name(), self.kash.name()):
            raise TypeError("%r is not an element of %r"%(f,self))

        return f

    def basis(self):
        """
        Return a basis of this order as a module over the maximal order
        of the base function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[]
            sage: L.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x^2*y, 1/x^4*y^2)

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x*y)

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); _.<t> = K[]
            sage: L.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x^2*y, 1/x^4*y^2)

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.basis()
            (1, 1/x*y)
        """

        return tuple(self.kash.Basis().sage(self._field.reverse_map))

    def ideal(self, *gens):
        """
        Return the fractional ideal generated by ``gens``.

        INPUT:

        - ``gens`` -- tuple of elements of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(x,y); I
            Ideal (x, y) of Maximal infinite order of Function field
            in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(x,y); I
            Ideal (x, y) of Maximal infinite order of Function field
            in y defined by y^2 + y + (x^2 + 1)/x
        """

        return FunctionFieldIdeal_kash(self, gens)

    def decomposition(self):
        """
        Return prime ideal decomposition of `pO_\infty` where `p` is the unique
        prime ideal of the maximal infinite order.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: Oinf = F.maximal_order_infinite()
            sage: Oinf.decomposition()
            [(Ideal (1/x, 1/x^2*y - 1) of Maximal infinite order
             of Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2, 1, 1),
             (Ideal (1/x, 1/x^4*y^2 + 1/x^2*y + 1) of Maximal infinite order
             of Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2, 2, 1)]

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.decomposition()
            [(Ideal (1/x, 1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x, 1, 2)]
        """

        return [ ( FunctionFieldIdeal_kash(self, i), i.Degree().sage(), i.RamificationIndex().sage() )
                 for i in self.kash.Decomposition()]

class FunctionFieldIdeal_kash(FunctionFieldIdeal):
    """
    Fractional ideal of the maximal order of a kash-implemented function field.
    """

    def __init__(self, ring, gens):
        """
        Initialize.

        INPUT:

        - ``ring`` -- maximal order

        - ``gens``-- generators
        """

        FunctionFieldIdeal.__init__(self, ring)

        if isinstance(gens, KashElement):
            self.kash = gens
            self._gens = tuple(gens.Generators().sage(ring._field.reverse_map))
            #self._gens = tuple(map(lambda x : x.sage(ring._field.reverse_map), gens.Generators()))
        else:
            self._gens = tuple(flatten(gens))
            self.kash = ring.kash.Ideal(map(ring._field.to_kash, self._gens))

    def __hash__(self):
        """
        Return hash computed from the data.
        """

        return hash( (self._ring, self._gens) )

    def _richcmp_(self, other, op):
        """
        Compare the ideal with the other ideal with respect to ``op``.

        INPUT:

        - ``other`` -- ideal

        - ``op`` -- comparison operator

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: I == I + I
            True
            sage: I == I * I
            False

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: I == I + I
            True
            sage: I == I * I
            False
            sage: I < I * I
            True
            sage: I > I * I
            False
        """

        return richcmp((self.denominator(), self.gens_over_base()), (other.denominator(), other.gens_over_base()), op)

    def __repr__(self):
        """
        Return a string representation of the ideal.
        """
        return "Ideal (%s) of %s"%(', '.join([repr(g) for g in self.gens()]), self.ring())

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is in the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal([y]); I
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 - x^3 - 1
            sage: x * y in I
            True
            sage: y / x in I
            False
            sage: y^2 - 2 in I
            False

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal([y]); I
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
            sage: x * y in I
            True
            sage: y / x in I
            False
            sage: y^2 - 2 in I
            False
        """
        kashx = self._ring._field.to_kash(x)
        return kash.eval("%s in %s" % (kashx.name(), self.kash.name())) == 'TRUE'

    def __invert__(self):
        """
        Return the inverse fractional ideal of the ideal.

        EXAMPLES::
            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: ~I
            Ideal (1, (1/(x^3 + 1))*y) of Maximal order of Function field in y defined by y^2 - x^3 - 1
            sage: I^(-1)
            Ideal (1, (1/(x^3 + 1))*y) of Maximal order of Function field in y defined by y^2 - x^3 - 1
            sage: ~I * I
            Ideal (1) of Maximal order of Function field in y defined by y^2 - x^3 - 1

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: ~I
            Ideal (x, (x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: I^(-1)
            Ideal (x, (x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: ~I * I
            Ideal (1) of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """

        return FunctionFieldIdeal_kash(self._ring, 1 / self.kash)

    def _add_(self, other):
        """
        Add with other ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I + J
            Ideal (x, y) of Maximal order of Function field in y defined by y^2 - x^3*y - x

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I + J
            Ideal (1, y) of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """

        return FunctionFieldIdeal_kash(self._ring, self.kash + other.kash)

    def _mul_(self, other):
        """
        Multiply with other ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I * J
            Ideal (x^4 + x^2 - x, x*y + x^2) of Maximal order
            of Function field in y defined by y^2 - x^3*y - x

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I * J
            Ideal ((x^5 + x^3 + x^2 + 1)/x,
                   y + (1/2*x^4 + 1/2*x^3 + x^2 + 1/2*x + 1/2)/x)
                of Maximal order of Function field in y
                defined by y^2 + y + (x^2 + 1)/x

        TESTS:

        Verify the examples using Sage's standard ideal operations::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: IJ = I * J

            sage: R = PolynomialRing(QQ, K.gens() + L.gens())
            sage: Q = R.quo(L.polynomial())
            sage: I1 = ideal(Q(g*b) for g in I.gens() for b in O.basis())
            sage: J1 = ideal(Q(g*b) for g in J.gens() for b in O.basis())
            sage: IJ = ideal(Q(g*b) for g in (IJ).gens() for b in O.basis())
            sage: I1 * J1 == IJ
            True

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: IJ = I * J

            sage: Q = R.quo(L.polynomial().numerator())
            sage: I1 = ideal(Q(g*b*x) for g in I.gens() for b in O.basis())
            sage: J1 = ideal(Q(g*b) for g in J.gens() for b in O.basis())
            sage: IJ = ideal(Q(g*b*x) for g in (IJ).gens() for b in O.basis())
            sage: I1 * J1 == IJ
            True
        """

        return FunctionFieldIdeal_kash(self._ring, self.kash * other.kash)

    def denominator(self):
        """
        Return the denominator of the fractional ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y/(y+1))
            sage: d = I.denominator(); d
            x^3
            sage: d in O
            True
        """
        return self.kash.Denominator().sage(self._ring._field.reverse_map)

    def factor(self):
        """
        Return the factorization of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()
            True

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()
            True
        """
        factors = self._factor()
        return Factorization(factors, cr=True)

    def _factor(self):
        """
        Return the list of prime and multiplicity pairs of the
        factorization of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[]
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()  # indirect doctest
            True
        """
        factors = []
        for f,m in self.kash.Factorization():
            factors.append( (FunctionFieldIdeal_kash(self._ring, f), m) )
        return factors

    def gens(self):
        """
        Return the generators of the ideal.

        This provides whatever set of generators as quickly
        as possible.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: I.gens()
            (y + x,)

            sage: L.<y> = K.extension(Y^2 +Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: I.gens()
            (y + x,)
        """
        return self._gens

    def gens_over_base(self):
        """
        Return the generators of the ideal as a module over the
        maximal order of the base rational function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: I.gens_over_base()
            (x^4 + x^2 - x, y + x)

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: I.gens_over_base()
            (x^3 + 1, y + x)
        """
        return tuple(self.kash.Basis().sage(self._ring._field.reverse_map))

    def is_prime(self):
        """
        Return ``True`` if the ideal is a prime ideal.

        If checked to be a prime ideal, then the ideal can be used
        as a prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I.is_prime()
            False
            sage: [f.is_prime() for f,_ in I.factor()]
            [True, True]

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I.is_prime()
            False
            sage: [f.is_prime() for f,_ in I.factor()]
            [True, True]
        """

        return bool(self.kash.IsPrime())

    def place(self):
        """
        Return the place corresponding to the prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.place() for f,_ in I.factor()]
            [Place (x, (1/(x^3 + x^2 + x))*y^2),
             Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)]

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.place() for f,_ in I.factor()]
            [Place (x, x*y), Place (x^2 + 1, x*y)]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")
        return FunctionFieldPlace_kash(self._ring._field, self)


class FunctionFieldPlace_kash(FunctionFieldPlace):
    """
    Places of kash-implemented function field.
    """

    def __init__(self, field, prime):
        """
        Initialize the place.

        INPUT:

        - ``field`` -- function field

        - ``prime`` -- prime ideal associated with the place
        """

        FunctionFieldPlace.__init__(self, field, prime)

        self.kash = prime.kash.Place()

    def is_infinite_place(self):
        """
        Return ``True`` if the place is at infinite.
        """
        F = self.function_field()
        return self.prime_ideal().ring() == F.maximal_order_infinite()
