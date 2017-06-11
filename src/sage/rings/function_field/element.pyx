r"""
Elements

Sage provides arithmetic with elements of function fields.

EXAMPLES:

Arithmetic with rational functions::

    sage: K.<t> = FunctionField(QQ)
    sage: f = t - 1
    sage: g = t^2 - 3
    sage: h = f^2/g^3
    sage: h.valuation(t-1)
    2
    sage: h.valuation(t)
    0
    sage: h.valuation(t^2 - 3)
    -3

Derivatives of elements in separable extensions::

    sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
    sage: (y^3 + x).derivative()
    ((x^2 + 1)/x^2)*y + (x^4 + x^3 + 1)/x^3

The divisor of an element of a global function field::

    sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
    sage: y.divisor()
    -1*Place (1/x, 1/x*y)
     - Place (x, x*y)
     + 2*Place (x + 1, x*y)

AUTHORS:

- William Stein: initial version

- Robert Bradshaw (2010-05-27): cythonize function field elements

- Julian Rueth (2011-06-28): treat zero correctly

- Maarten Derickx (2011-09-11): added doctests, fixed pickling

- Kwankyu Lee (2017-04-30): added elements for global function fields

"""
#*****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#       Copyright (C) 2010 Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2011 Julian Rueth <julian.rueth@gmail.com>
#       Copyright (C) 2011 Maarten Derickx <m.derickx.student@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method

from sage.structure.element cimport FieldElement, RingElement, ModuleElement, Element
from sage.structure.richcmp cimport richcmp, richcmp_not_equal

def is_FunctionFieldElement(x):
    """
    Return True if ``x`` is any type of function field element.

    EXAMPLES::

        sage: t = FunctionField(QQ,'t').gen()
        sage: sage.rings.function_field.element.is_FunctionFieldElement(t)
        True
        sage: sage.rings.function_field.element.is_FunctionFieldElement(0)
        False
    """
    if isinstance(x, FunctionFieldElement): return True
    from function_field import is_FunctionField
    return is_FunctionField(x.parent())

def make_FunctionFieldElement(parent, element_class, representing_element):
    """
    Used for unpickling FunctionFieldElement objects (and subclasses).

    EXAMPLES::

        sage: from sage.rings.function_field.element import make_FunctionFieldElement
        sage: K.<x> = FunctionField(QQ)
        sage: make_FunctionFieldElement(K, K._element_class, (x+1)/x)
        (x + 1)/x
    """
    return element_class(parent, representing_element, reduce=False)

cdef class FunctionFieldElement(FieldElement):
    """
    Abstract base class for function field elements.

    EXAMPLES::

        sage: t = FunctionField(QQ,'t').gen()
        sage: isinstance(t, sage.rings.function_field.element.FunctionFieldElement)
        True
    """
    cdef readonly object _x
    cdef readonly object _matrix

    def __reduce__(self):
        """
        EXAMPLES::

            sage: K = FunctionField(QQ,'x')
            sage: f = K.random_element()
            sage: loads(f.dumps()) == f
            True
        """
        return (make_FunctionFieldElement,
                (self._parent, type(self), self._x))

    cdef FunctionFieldElement _new_c(self):
        cdef type t = type(self)
        cdef FunctionFieldElement x = <FunctionFieldElement>t.__new__(t)
        x._parent = self._parent
        return x

    def __pari__(self):
        r"""
        Coerce the element to PARI.

        PARI does not know about general function field elements, so this
        raises an Exception.

        TESTS:

        Check that :trac:`16369` has been resolved::

            sage: K.<a> = FunctionField(QQ)
            sage: R.<b> = K[]
            sage: L.<b> = K.extension(b^2-a)
            sage: b.__pari__()
            Traceback (most recent call last):
            ...
            NotImplementedError: PARI does not support general function field elements.

        """
        raise NotImplementedError("PARI does not support general function field elements.")

    def _latex_(self):
        """
        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: latex((t+1)/t)
            \frac{t + 1}{t}
            sage: latex((t+1)/t^67)
            \frac{t + 1}{t^{67}}
            sage: latex((t+1/2)/t^67)
            \frac{t + \frac{1}{2}}{t^{67}}
        """
        return self._x._latex_()

    def matrix(self):
        r"""
        Return the matrix of multiplication by the element, interpreted as an
        element of a vector space over its base field.

        EXAMPLES:

        A rational function field::

            sage: K.<t> = FunctionField(QQ)
            sage: t.matrix()
            [t]
            sage: (1/(t+1)).matrix()
            [1/(t + 1)]

        Now an example in a nontrivial extension of a rational function field::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: y.matrix()
            [     0      1]
            [-4*x^3      x]
            sage: y.matrix().charpoly('Z')
            Z^2 - x*Z + 4*x^3

        An example in a relative extension, where neither function
        field is rational::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: M.<T> = L[]; Z.<alpha> = L.extension(T^3 - y^2*T + x)
            sage: alpha.matrix()
            [          0           1           0]
            [          0           0           1]
            [         -x x*y - 4*x^3           0]

        We show that this matrix does indeed work as expected when making a
        vector space from a function field::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x))
            sage: V, from_V, to_V = L.vector_space()
            sage: y5 = to_V(y^5); y5
            ((x^4 + 1)/x, 2*x, 0, 0, 0)
            sage: y4y = to_V(y^4) * y.matrix(); y4y
            ((x^4 + 1)/x, 2*x, 0, 0, 0)
            sage: y5 == y4y
            True
        """
        if self._matrix is None:
            # Multiply each power of field generator on the left by this
            # element; make matrix whose rows are the coefficients of the
            # result, and transpose.
            K = self.parent()
            v = []
            x = K.gen()
            a = K(1)
            d = K.degree()
            for n in range(d):
                v += (a*self).list()
                a *= x
            k = K.base_ring()
            import sage.matrix.matrix_space
            M = sage.matrix.matrix_space.MatrixSpace(k, d)
            self._matrix = M(v)
        return self._matrix

    def trace(self):
        """
        Return the trace of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: y.trace()
            x
        """
        return self.matrix().trace()

    def norm(self):
        """
        Return the norm of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: y.norm()
            4*x^3

        The norm is relative::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3); R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y^2*z + x)
            sage: z.norm()
            -x
            sage: z.norm().parent()
            Function field in y defined by y^2 - x*y + 4*x^3
        """
        return self.matrix().determinant()

    def characteristic_polynomial(self, *args, **kwds):
        """
        Return the characteristic polynomial of the element. Give an optional
        input string to name the variable in the characteristic polynomial.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3); R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y^2*z + x)
            sage: x.characteristic_polynomial('W')
            W - x
            sage: y.characteristic_polynomial('W')
            W^2 - x*W + 4*x^3
            sage: z.characteristic_polynomial('W')
            W^3 + (-x*y + 4*x^3)*W + x
        """
        return self.matrix().characteristic_polynomial(*args, **kwds)

    charpoly = characteristic_polynomial

    def minimal_polynomial(self, *args, **kwds):
        """
        Return the minimal polynomial of the element. Give an optional input
        string to name the variable in the characteristic polynomial.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3); R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y^2*z + x)
            sage: x.minimal_polynomial('W')
            W - x
            sage: y.minimal_polynomial('W')
            W^2 - x*W + 4*x^3
            sage: z.minimal_polynomial('W')
            W^3 + (-x*y + 4*x^3)*W + x
        """
        return self.matrix().minimal_polynomial(*args, **kwds)

    minpoly = minimal_polynomial

    def is_integral(self):
        r"""
        Determine if the element is integral over the maximal order of the base field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: y.is_integral()
            True
            sage: (y/x).is_integral()
            True
            sage: (y/x)^2 - (y/x) + 4*x
            0
            sage: (y/x^2).is_integral()
            False
            sage: (y/x).minimal_polynomial('W')
            W^2 - W + 4*x
        """
        R = self.parent().base_field().maximal_order()
        return all([a in R for a in self.minimal_polynomial()])

cdef class FunctionFieldElement_polymod(FunctionFieldElement):
    """
    Elements of a finite extension of a function field.
    """
    def __init__(self, parent, x, reduce=True):
        """
        Initialize.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: x*y + 1/x^3
            x*y + 1/x^3
        """
        FieldElement.__init__(self, parent)
        if reduce:
            self._x = x % self._parent.polynomial()
        else:
            self._x = x

    def element(self):
        """
        Return the underlying polynomial that represents the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<T> = K[]
            sage: L.<y> = K.extension(T^2 - x*T + 4*x^3)
            sage: f = y/x^2 + x/(x^2+1); f
            1/x^2*y + x/(x^2 + 1)
            sage: f.element()
            1/x^2*y + x/(x^2 + 1)
        """
        return self._x

    def _repr_(self):
        """
        Return the string representation of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: y._repr_()
            'y'
        """
        return self._x._repr(name=self.parent().variable_name())

    def __nonzero__(self):
        """
        Return True if the element is not zero.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: bool(y)
            True
            sage: bool(L(0))
            False
            sage: bool(L.coerce(L.polynomial()))
            False
        """
        return not not self._x

    def __hash__(self):
        """
        Return the hash of the element.

        TESTS::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: len({hash(y^i+x^j) for i in [-2..2] for j in [-2..2]}) == 25
            True
        """
        return hash(self._x)

    cpdef _richcmp_(self, other, int op):
        """
        Do rich comparison with the other element with respect to ``op``

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: L(0) == 0
            True
            sage: y != L(2)
            True
        """
        cdef FunctionFieldElement left = <FunctionFieldElement>self
        cdef FunctionFieldElement right = <FunctionFieldElement>other
        return richcmp(left._x, right._x, op)

    cpdef _add_(self, right):
        """
        Add the element with the other element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: (2*y + x/(1+x^3))  +  (3*y + 5*x*y)         # indirect doctest
            (5*x + 5)*y + x/(x^3 + 1)
            sage: (y^2 - x*y + 4*x^3)==0                      # indirect doctest
            True
            sage: -y+y
            0
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x + (<FunctionFieldElement>right)._x
        return res

    cpdef _sub_(self, right):
        """
        Subtract the other element from the element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: (2*y + x/(1+x^3))  -  (3*y + 5*x*y)         # indirect doctest
            (-5*x - 1)*y + x/(x^3 + 1)
            sage: y-y
            0
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x - (<FunctionFieldElement>right)._x
        return res

    cpdef _mul_(self, right):
        """
        Multiply the element with the other element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: y  *  (3*y + 5*x*y)                          # indirect doctest
            (5*x^2 + 3*x)*y - 20*x^4 - 12*x^3
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = (self._x * (<FunctionFieldElement>right)._x) % self._parent.polynomial()
        return res

    cpdef _div_(self, right):
        """
        Divide the element with the other element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: (2*y + x/(1+x^3))  /  (2*y + x/(1+x^3))       # indirect doctest
            1
            sage: 1 / (y^2 - x*y + 4*x^3)                       # indirect doctest
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot invert 0
        """
        return self * ~right

    def __invert__(self):
        """
        Return the multiplicative inverse of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: a = ~(2*y + 1/x); a                           # indirect doctest
            (-x^2/(8*x^5 + x^2 + 1/2))*y + (2*x^3 + x)/(16*x^5 + 2*x^2 + 1)
            sage: a*(2*y + 1/x)
            1
        """
        if self.is_zero():
            raise ZeroDivisionError("Cannot invert 0")
        P = self._parent
        return P(self._x.xgcd(P._polynomial)[1])

    def list(self):
        """
        Return the list of the coefficients prepresenting the element.

        If the function field is `K[y]/(f(y))`, then return the coefficients of
        the reduced presentation of the element as a polynomial in `K[y]`.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: a = ~(2*y + 1/x); a
            (-x^2/(8*x^5 + x^2 + 1/2))*y + (2*x^3 + x)/(16*x^5 + 2*x^2 + 1)
            sage: a.list()
            [(2*x^3 + x)/(16*x^5 + 2*x^2 + 1), -x^2/(8*x^5 + x^2 + 1/2)]
            sage: (x*y).list()
            [0, x]
        """
        return self._x.padded_list(self.parent().degree())

cdef class FunctionFieldElement_rational(FunctionFieldElement):
    """
    Elements of a rational function field.
    """
    def __init__(self, parent, x, reduce=True):
        """
        Initialize.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ); K
            Rational function field in t over Rational Field
            sage: FunctionField(QQ,'t').gen()^3
            t^3
        """
        FieldElement.__init__(self, parent)
        self._x = x

    def __pari__(self):
        r"""
        Coerce the element to PARI.

        EXAMPLES::

            sage: K.<a> = FunctionField(QQ)
            sage: ((a+1)/(a-1)).__pari__()
            (a + 1)/(a - 1)

        """
        return self.element().__pari__()

    def element(self):
        """
        Return the underlying fraction field element that represents the element.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(7))
            sage: t.element()
            t
            sage: type(t.element())
            <type 'sage.rings.fraction_field_FpT.FpTElement'>

            sage: K.<t> = FunctionField(GF(131101))
            sage: t.element()
            t
            sage: type(t.element())
            <class 'sage.rings.fraction_field_element.FractionFieldElement_1poly_field'>
        """
        return self._x

    def list(self):
        """
        Return a list with just the element.

        The list represents the element when the rational function field is
        viewed as a (one-dimensional) vector space over itself.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t.list()
            [t]
        """
        return [self]

    def _repr_(self):
        """
        Return the string representation of the element.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t._repr_()
            't'
        """
        return repr(self._x)

    def __nonzero__(self):
        """
        Return True if the element is not zero.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: bool(t)
            True
            sage: bool(K(0))
            False
            sage: bool(K(1))
            True
        """
        return not not self._x

    def __hash__(self):
        """
        Return the hash of the element.

        TESTS:

        It would be nice if the following would produce a list of
        15 distinct hashes::

            sage: K.<t> = FunctionField(QQ)
            sage: len({hash(t^i+t^j) for i in [-2..2] for j in [i..2]})
            10
        """
        return hash(self._x)

    cpdef _richcmp_(self, other, int op):
        """
        Compare the element with the other element with respect to ``op``

        INPUT:

        - ``other`` -- element

        - ``op`` -- comparison operator

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t > 0
            True
            sage: t < t^2
            True
        """
        cdef FunctionFieldElement left
        cdef FunctionFieldElement right
        try:
            left = <FunctionFieldElement?>self
            right = <FunctionFieldElement?>other
            lp = left._parent
            rp = right._parent
            if lp != rp:
                return richcmp_not_equal(lp, rp, op)
            return richcmp(left._x, right._x, op)
        except TypeError:
            return NotImplemented

    cpdef _add_(self, right):
        """
        Add the element with the other element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t + (3*t^3)                      # indirect doctest
            3*t^3 + t
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x + (<FunctionFieldElement>right)._x
        return res

    cpdef _sub_(self, right):
        """
        Subtract the other element from the element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t - (3*t^3)                      # indirect doctest
            -3*t^3 + t
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x - (<FunctionFieldElement>right)._x
        return res

    cpdef _mul_(self, right):
        """
        Multiply the element with the other element

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: (t+1) * (t^2-1)                  # indirect doctest
            t^3 + t^2 - t - 1
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x * (<FunctionFieldElement>right)._x
        return res

    cpdef _div_(self, right):
        """
        Divide the element with the other element

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: (t+1) / (t^2 - 1)                # indirect doctest
            1/(t - 1)
        """
        cdef FunctionFieldElement res = self._new_c()
        res._parent = self._parent.fraction_field()
        res._x = self._x / (<FunctionFieldElement>right)._x
        return res

    def numerator(self):
        """
        Return the numerator of the rational function.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t+1) / (t^2 - 1/3); f
            (t + 1)/(t^2 - 1/3)
            sage: f.numerator()
            t + 1
        """
        return self._x.numerator()

    def denominator(self):
        """
        Return the denominator of the rational function.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t+1) / (t^2 - 1/3); f
            (t + 1)/(t^2 - 1/3)
            sage: f.denominator()
            t^2 - 1/3
        """
        return self._x.denominator()

    def _valuation(self, v):
        """
        Return the valuation of the rational function with respect to
        ``v``

        INPUT:

        - ``v`` -- irreducible polynomial

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t-1)^2 * (t+1) / (t^2 - 1/3)^3
            sage: f.valuation(t-1)  # indirect doctest
            2
            sage: f.valuation(t)  # indirect doctest
            0
            sage: f.valuation(t^2 - 1/3)  # indirect doctest
            -3
        """
        R = self._parent._ring
        return self._x.valuation(R(self._parent(v)._x))

    def valuation(self, place):
        """
        Return the valuation of the rational function at the place

        INPUT:

        - ``place`` -- place of the rational function field or an irreducible
          polynomial associated with a place

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t-1)^2 * (t+1) / (t^2 - 1/3)^3
            sage: p = (t-1).divisor_of_zeros().support()[0]
            sage: f.valuation(p)
            2
            sage: f.valuation(t-1)
            2
        """
        from .place import is_Place

        if not is_Place(place):
            return self._valuation(place)

        prime = place.prime_ideal()
        ideal = prime.ring().ideal(self)
        return prime.valuation(ideal)

    def is_square(self):
        """
        Return whether the element is a square.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t.is_square()
            False
            sage: (t^2/4).is_square()
            True
            sage: f = 9 * (t+1)^6 / (t^2 - 2*t + 1); f.is_square()
            True

            sage: K.<t> = FunctionField(GF(5))
            sage: (-t^2).is_square()
            True
            sage: (-t^2).sqrt()
            2*t
        """
        return self._x.is_square()

    def sqrt(self, all=False):
        """
        Return the square root of the rational function.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = t^2 - 2 + 1/t^2; f.sqrt()
            (t^2 - 1)/t
            sage: f = t^2; f.sqrt(all=True)
            [t, -t]

        TESTS::

            sage: K(4/9).sqrt()
            2/3
            sage: K(0).sqrt(all=True)
            [0]
        """
        if all:
            return [self._parent(r) for r in self._x.sqrt(all=True)]
        else:
            return self._parent(self._x.sqrt())

    def factor(self):
        """
        Factor the rational function.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t+1) / (t^2 - 1/3)
            sage: f.factor()
            (t + 1) * (t^2 - 1/3)^-1
            sage: (7*f).factor()
            (7) * (t + 1) * (t^2 - 1/3)^-1
            sage: ((7*f).factor()).unit()
            7
            sage: (f^3).factor()
            (t + 1)^3 * (t^2 - 1/3)^-3
        """
        P = self.parent()
        F = self._x.factor()
        from sage.structure.factorization import Factorization
        return Factorization([(P(a),e) for a,e in F], unit=F.unit())

    def inverse_mod(self, I):
        """
        Return an inverse of the element modulo the integral ideal `I`, if `I`
        and the element together generate the unit ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order(); I = O.ideal(x^2+1)
            sage: t = O(x+1).inverse_mod(I); t
            -1/2*x + 1/2
            sage: (t*(x+1) - 1) in I
            True
        """
        assert  len(I.gens()) == 1
        f = I.gens()[0]._x
        assert f.denominator() == 1
        assert self._x.denominator() == 1
        return self.parent()(self._x.numerator().inverse_mod(f.numerator()))

    @cached_method
    def divisor(self):
        """
        Return the divisor of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: f = 1/(x^3 + x^2 + x)
            sage: f.divisor()
            3*Place (1/x)
             - Place (x)
             - Place (x^2 + x + 1)
        """
        if self.is_zero():
            raise ValueError("divisor not defined for zero")

        F = self.parent()
        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        I = O.ideal(self)
        J = Oinf.ideal(self)
        return I.divisor() + J.divisor()

    def divisor_of_zeros(self):
        """
        Return the divisor of zeros for the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: f = 1/(x^3 + x^2 + x)
            sage: f.divisor_of_zeros()
            3*Place (1/x)
        """
        if self.is_zero():
            raise ValueError("divisor of zeros not defined for zero")

        F = self.parent()
        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        I = O.ideal(self)
        J = Oinf.ideal(self)
        return I.divisor_of_zeros() + J.divisor_of_zeros()

    def divisor_of_poles(self):
        """
        Return the divisor of poles for the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: f = 1/(x^3 + x^2 + x)
            sage: f.divisor_of_poles()
            Place (x)
             + Place (x^2 + x + 1)
        """
        if self.is_zero():
            raise ValueError("divisor of poles not defined for zero")

        F = self.parent()
        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        I = O.ideal(self)
        J = Oinf.ideal(self)
        return I.divisor_of_poles() + J.divisor_of_poles()

    def differential(self):
        """
        Return the differential `dx` where `x` is the element.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = 1 / t
            sage: f.differential()
            (-1/t^2) d(t)
        """
        from .differential import differential

        return differential(self.parent(), 1, self)

    def derivative(self):
        """
        Return the derivative of the rational function.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t + 1) / (t^2 - 1/3)
            sage: f.derivative()
            (-t^2 - 2*t - 1/3)/(t^4 - 2/3*t^2 + 1/9)
        """
        D = self.parent().derivation()
        return D(self)

    def hasse_derivative(self, i, separating_element=None):
        """
        Return the `i`-th Hasse derivative of the element with respect to the
        separating element.

        INPUT:

        - ``i`` -- nonnegative integer

        - ``separating_element`` -- separating element of the function field;
            the default is the generator of the rational function field

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(2))
            sage: f = t^2
            sage: f.hasse_derivative(2)
            1
        """
        D = self.parent().hasse_derivation()
        return D(self, i, separating_element)

cdef class FunctionFieldElement_global(FunctionFieldElement_polymod):
    """
    Elements of global function fields
    """
    cdef public dict __cached_methods

    def valuation(self, place):
        """
        Return the valuation of the element at the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_infinite()[0]
            sage: y.valuation(p)
            -1
        """
        prime = place.prime_ideal()
        ideal = prime.ring().ideal(self)
        return prime.valuation(ideal)

    # for consistency with rational function field elements
    _valuation = valuation

    @cached_method
    def divisor(self):
        """
        Return the divisor for the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: y.divisor()
            -1*Place (1/x, 1/x*y)
             - Place (x, x*y)
             + 2*Place (x + 1, x*y)
        """
        if self.is_zero():
            raise ValueError("not defined for zero")

        F = self.parent()
        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        I = O.ideal(self)
        J = Oinf.ideal(self)
        return I.divisor() + J.divisor()

    def divisor_of_zeros(self):
        """
        Return divisor of zeros for the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: (x/y).divisor_of_zeros()
            3*Place (x, x*y)
        """
        if self.is_zero():
            raise ValueError("divisor of zeros not defined for zero")

        F = self.parent()
        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        I = O.ideal(self)
        J = Oinf.ideal(self)
        return I.divisor_of_zeros() + J.divisor_of_zeros()

    def divisor_of_poles(self):
        """
        Return the divisor of poles for the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: (x/y).divisor_of_poles()
            Place (1/x, 1/x*y) + 2*Place (x + 1, x*y)
        """
        if self.is_zero():
            raise ValueError("divisor of poles not defined for zero")

        F = self.parent()
        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        I = O.ideal(self)
        J = Oinf.ideal(self)
        return I.divisor_of_poles() + J.divisor_of_poles()

    def differential(self):
        """
        Return the differential `dx` where `x` is the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x +1/x)
            sage: (y^3 + x).differential()
            (((x^2 + 1)/x^2)*y + (x^4 + x^3 + 1)/x^3) d(x)
        """
        M = self.parent().space_of_differentials()
        return M(self)

    def derivative(self):
        """
        Return the derivative of the element.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: (y^3 + x).derivative()
            ((x^2 + 1)/x^2)*y + (x^4 + x^3 + 1)/x^3
        """
        D = self.parent().derivation()
        return D(self)

    def hasse_derivative(self, i, separating_element=None):
        """
        Return the ``i``-th order Hasse derivative of the element with respect
        to the separating element.

        INPUT:

        - ``i`` -- nonnegative integer

        - ``separating_element`` -- separating element of the function field;
          the default is the generator of the base rational function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: (y^3 + x).hasse_derivative(2)
            1/x^3*y + (x^6 + x^4 + x^3 + x^2 + x + 1)/x^5
        """
        D = self.parent().hasse_derivation()
        return D(self, i, separating_element)
