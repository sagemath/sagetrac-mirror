r"""
Elements of Puiseux polynomial rings
"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.rings.integer cimport Integer
from sage.structure.element import is_Element, coerce_binop
from sage.misc.misc import union
from sage.structure.factorization import Factorization
from sage.misc.derivative import multi_derivative
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.structure.richcmp cimport richcmp, rich_to_bool
from sage.rings.rational cimport Rational
from sage.arith.all import lcm


cdef class PuiseuxPolynomial(CommutativeAlgebraElement):
    """
    Base class for Puiseux polynomials.
    """
    cdef PuiseuxPolynomial _new_c(self):
        """
        Return a new Puiseux polynomial.

        EXAMPLES::

            sage: L.<x> = PuiseuxPolynomialRing(QQ) # indirect doctest
            sage: x*x
            x*y
        """
        cdef type t = type(self)
        cdef PuiseuxPolynomial ans
        ans = t.__new__(t)
        ans._parent = self._parent
        return ans

    cpdef _add_(self, other):
        """
        Abstract addition method

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(ZZ)
            sage: from sage.rings.polynomial.laurent_polynomial import PuiseuxPolynomial
            sage: PuiseuxPolynomial._add_(x, x)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    cpdef _mul_(self, other):
        """
        Abstract multiplication method

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(ZZ)
            sage: from sage.rings.polynomial.laurent_polynomial import PuiseuxPolynomial
            sage: PuiseuxPolynomial._mul_(x, x)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    cpdef _floordiv_(self, other):
        """
        Abstract floor division method

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(ZZ)
            sage: from sage.rings.polynomial.laurent_polynomial import PuiseuxPolynomial
            sage: PuiseuxPolynomial._floordiv_(x, x)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _integer_(self, ZZ):
        r"""
        Convert this Puiseux polynomial to an integer.

        This is only possible if the Puiseux polynomial is constant.

        OUTPUT:

        An integer.

        TESTS::

            sage: L.<a> = PuiseuxPolynomialRing(QQ)
            sage: L(42)._integer_(ZZ)
            42
            sage: a._integer_(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: a is not constant
            sage: L(2/3)._integer_(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: ZZ(L(42))
            42

        ::

            sage: L.<a, b> = PuiseuxPolynomialRing(QQ)
            sage: L(42)._integer_(ZZ)
            42
            sage: a._integer_(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: a is not constant
            sage: L(2/3)._integer_(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: ZZ(L(42))
            42
        """
        if not self.is_constant():
            raise ValueError('{} is not constant'.format(self))
        return ZZ(self.constant_coefficient())

    def _rational_(self):
        r"""
        Convert this Puiseux polynomial to a rational.

        This is only possible if the Puiseux polynomial is constant.

        OUTPUT:

        A rational.

        TESTS::

            sage: L.<a> = PuiseuxPolynomialRing(QQ)
            sage: L(42)._rational_()
            42
            sage: a._rational_()
            Traceback (most recent call last):
            ...
            ValueError: a is not constant
            sage: QQ(L(2/3))
            2/3

        ::

            sage: L.<a, b> = PuiseuxPolynomialRing(QQ)
            sage: L(42)._rational_()
            42
            sage: a._rational_()
            Traceback (most recent call last):
            ...
            ValueError: a is not constant
            sage: QQ(L(2/3))
            2/3
        """
        if not self.is_constant():
            raise ValueError('{} is not constant'.format(self))
        from sage.rings.rational_field import QQ
        return QQ(self.constant_coefficient())

    def change_ring(self, R):
        """
        Return a copy of this Puiseux polynomial, with coefficients in ``R``.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: a = x^2 + 3*x^3 + 5*x^-1
            sage: a.change_ring(GF(3))
            2*x^-1 + x^2

        Check that :trac:`22277` is fixed::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: a = 2*x^2 + 3*x^3 + 4*x^-1
            sage: a.change_ring(GF(3))
            -x^2 + x^-1
        """
        return self._parent.change_ring(R)(self)

    cpdef long number_of_terms(self) except -1:
        """
        Abstract method for number of terms

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(ZZ)
            sage: from sage.rings.polynomial.laurent_polynomial import PuiseuxPolynomial
            sage: PuiseuxPolynomial.number_of_terms(x)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def hamming_weight(self):
        """
        Return the hamming weight of ``self``.

        The hamming weight is number of non-zero coefficients and
        also known as the weight or sparsity.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(ZZ)
            sage: f = x^3 - 1
            sage: f.hamming_weight()
            2
        """
        return self.number_of_terms()

    cpdef dict dict(self):
        """
        Abstract ``dict`` method.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(ZZ)
            sage: from sage.rings.polynomial.laurent_polynomial import PuiseuxPolynomial
            sage: PuiseuxPolynomial.dict(x)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


cdef class PuiseuxPolynomial_univariate(PuiseuxPolynomial):
    """
    A univariate Puiseux polynomial in `t`.

    INPUT:

    - ``parent`` -- a Puiseux polynomial ring

    - ``f`` -- a Puiseux polynomial (or something can be coerced to one)

    AUTHORS:

    - Tom Boothby (2011) copied this class almost verbatim from
      ``laurent_series_ring_element.pyx``, so most of the credit goes to
      William Stein, David Joyner, and Robert Bradshaw
    - Travis Scrimshaw (09-2013): Cleaned-up and added a few extra methods
    """
    def __init__(self, parent, f, n=None):
        r"""
        Create the Puiseux polynomial `f`.

        EXAMPLES::

            sage: R.<q> = PuiseuxPolynomialRing(ZZ)
            sage: R([1,2,3],1/2)
            1 + 2*q + 3*q^2
            sage: TestSuite(q^-3 + 3*q + 2).run()

        ::

            sage: S.<s> = PuiseuxPolynomialRing(GF(5))
            sage: T.<t> = PolynomialRing(pAdicRing(5))
            sage: S(t)
            s
            sage: parent(S(t))
            Univariate Puiseux Polynomial Ring in s over Finite Field of size 5
            sage: parent(S(t)[1])
            Finite Field of size 5

        ::

            sage: R({})
            0
        """
        CommutativeAlgebraElement.__init__(self, parent)

        if isinstance(f, PuiseuxPolynomial_univariate):
            if (<PuiseuxPolynomial_univariate>f).__u._parent is parent._R:
                F = (<PuiseuxPolynomial_univariate>f)
                self.__u = F.__u
                self.__n = F.__n
            else:
                F = (<PuiseuxPolynomial_univariate>f)
                self.__u = parent._R(F.__u)
                self.__n = F.__n
        elif isinstance(f, Polynomial):
            self.__u = f
            if n is not None:
                self.__n = n
            else:
                self.__n = 1
        elif isinstance(f, dict):
            n = lcm([exp.denominator() for exp in f])
            f = parent._R({Integer(n * i): c for i, c in f.items()})
            self.__u = f
            self.__n = n

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: R.<q> = PuiseuxPolynomialRing(ZZ)
            sage: elt = q^-3 + 2 + q
            sage: loads(dumps(elt)) == elt
            True
        """
        return PuiseuxPolynomial_univariate, (self._parent, self.__u)

    def is_unit(self):
        """
        Return ``True`` if this Puiseux polynomial is a unit in this ring.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: (2+t).is_unit()
            False
            sage: f = 2*t
            sage: f.is_unit()
            True
            sage: 1/f
            1/2*t^-1
            sage: R(0).is_unit()
            False
            sage: R.<s> = PuiseuxPolynomialRing(ZZ)
            sage: g = 2*s
            sage: g.is_unit()
            False
            sage: 1/g
            1/2*s^-1

        ALGORITHM: A Puiseux polynomial is a unit if and only if its "unit
        part" is a unit.
        """
        return self.__u.is_constant() and self.__u.coefficients()[0].is_unit()

    def is_zero(self):
        """
        Return ``1`` if ``self`` is 0, else return ``0``.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = 1/x + x + x^2 + 3*x^4
            sage: f.is_zero()
            0
            sage: z = 0*f
            sage: z.is_zero()
            1
        """
        return self.__u.is_zero()

    def __nonzero__(self):
        """
        Check if ``self`` is non-zero.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = 1/x + x + x^2 + 3*x^4
            sage: not f
            False
            sage: z = 0*f
            sage: not z
            True
        """
        return not self.__u.is_zero()

    def _im_gens_(self, codomain, im_gens):
        """
        Return the image of ``self`` under the morphism defined by
        ``im_gens`` in ``codomain``.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: H = Hom(R, QQ)
            sage: mor = H(2)
            sage: mor(t^2 + t^-2)
            17/4
            sage: 4 + 1/4
            17/4
        """
        return codomain(self(im_gens[0]))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: 2 + t^(1/2) + (2/3)*t^3
            2 + 2/3*t^3
        """
        if self.is_zero():
            return "0"
        s = " "
        v = self.__u.list()
        m = len(v)
        X = self._parent.variable_name()
        atomic_repr = self._parent.base_ring()._repr_option('element_is_atomic')
        first = True
        for n in xrange(m):
            x = v[n]
            e = Rational(n) / self.__n
            if x != 0:
                x = str(x)
                if not first:
                    s += " + "
                if not atomic_repr and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                    x = "({})".format(x)
                if e == 1:
                    var = "*{}".format(X)
                elif e == 0:
                    var = ""
                else:
                    var = "*{}^{}".format(X, e)
                s += "{}{}".format(x, var)
                first = False
        s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        return s[1:]

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = (17/2)*x^-2 + x + x^2 + 3*x^4
            sage: latex(f)
            \frac{\frac{17}{2}}{x^{2}} + x + x^{2} + 3x^{4}

        Verify that :trac:`6656` has been fixed::

            sage: R.<a,b>=PolynomialRing(QQ)
            sage: T.<x>=PuiseuxPolynomialRing(R)
            sage: y = a*x+b*x
            sage: y._latex_()
            '\\left(a + b\\right)x'
            sage: latex(y)
            \left(a + b\right)x

        TESTS::

            sage: L.<lambda2> = PuiseuxPolynomialRing(QQ)
            sage: latex(L.an_element())
            \lambda_{2}
            sage: L.<y2> = PuiseuxPolynomialRing(QQ)
            sage: latex(L.an_element())
            y_{2}
        """
        from sage.misc.latex import latex

        if self.is_zero():
            return "0"
        s = " "
        v = self.__u.list()
        valuation = self.__n
        m = len(v)
        X = self._parent.latex_variable_names()[0]
        atomic_repr = self._parent.base_ring()._repr_option('element_is_atomic')
        first = True
        for n in xrange(m):
            x = v[n]
            e = n + valuation
            x = latex(x)
            if x != '0':
                if not first:
                    s += " + "
                if not atomic_repr and e > 0 and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                    x = "\\left({}\\right)".format(x)
                if e == 1:
                    var = "|{}".format(X)
                elif e == 0:
                    var = ""
                elif e > 0:
                    var = "|{}^{{{}}}".format(X,e)
                if e >= 0:
                    s += "{}{}".format(x,var)
                else: # negative e
                    if e == -1:
                        s += "\\frac{{{}}}{{{}}}".format(x, X)
                    else:
                        s += "\\frac{{{}}}{{{}^{{{}}}}}".format(x, X,-e)
                first = False
        s = s.replace(" + -", " - ")
        s = s.replace(" 1|"," ")
        s = s.replace(" -1|", " -")
        s = s.replace("|","")

        return s[1:]

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: R = PuiseuxPolynomialRing(QQ, 't')

            sage: assert hash(R.zero()) == 0
            sage: assert hash(R.one()) == 1
            sage: assert hash(QQ['t'].gen()) == hash(R.gen())

            sage: for _ in range(20):
            ....:     p = QQ.random_element()
            ....:     assert hash(R(p)) == hash(p), "p = {}".format(p)

            sage: S.<t> = QQ[]
            sage: for _ in range(20):
            ....:     p = S.random_element()
            ....:     assert hash(R(p)) == hash(p), "p = {}".format(p)
            ....:     assert hash(R(t*p)) == hash(t*p), "p = {}".format(p)

        Check that :trac:`21272` is fixed::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: hash(R.zero()) == hash(t - t)
            True
        """
        return hash((self.__u, self.__n))

    def __getitem__(self, i):
        """
        Return the `i`-th coefficient of ``self``.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: f = -5/t^(10) + t + t^2 - 10/3*t^3; f
            -5*t^-10 + t + t^2 - 10/3*t^3
            sage: f[-10]
            -5
            sage: f[1]
            1
            sage: f[3]
            -10/3
            sage: f[-9]
            0
            sage: f = -5/t^(10) + 1/3 + t + t^2 - 10/3*t^3; f
            -5*t^-10 + 1/3 + t + t^2 - 10/3*t^3

        Slicing is deprecated::

            sage: f[-10:2]
            doctest:...: DeprecationWarning: polynomial slicing with a start index is deprecated, use list() and slice the resulting list instead
            See http://trac.sagemath.org/18940 for details.
            -5*t^-10 + 1/3 + t
            sage: f[0:]
            1/3 + t + t^2 - 10/3*t^3
            sage: f[:3]
            -5*t^-10 + 1/3 + t + t^2
            sage: f[-14:5:2]
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomial slicing with a step is not defined
        """
        cdef PuiseuxPolynomial_univariate ret
        if isinstance(i, slice):
            start = i.start - self.__n if i.start is not None else 0
            stop = i.stop - self.__n if i.stop is not None else self.__u.degree() + 1
            f = self.__u[start:stop:i.step]  # deprecation(18940)
            ret = <PuiseuxPolynomial_univariate> self._new_c()
            ret.__u = f
            ret.__n = self.__n
            return ret

        return self.__u[i - self.__n]

    cpdef long number_of_terms(self) except -1:
        """
        Return the number of non-zero coefficients of ``self``.

        Also called weight, hamming weight or sparsity.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(ZZ)
            sage: f = x^3 - 1
            sage: f.number_of_terms()
            2
            sage: R(0).number_of_terms()
            0
            sage: f = (x+1)^100
            sage: f.number_of_terms()
            101

        The method :meth:`hamming_weight` is an alias::

            sage: f.hamming_weight()
            101
        """
        return self.__u.number_of_terms()

    def __iter__(self):
        """
        Iterate through the coefficients from the first nonzero one to the
        last nonzero one.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: f = -5/t^(2) + t + t^2 - 10/3*t^3; f
            -5*t^-2 + t + t^2 - 10/3*t^3
            sage: for a in f: print(a)
            -5
            0
            0
            1
            1
            -10/3
        """
        return iter(self.__u)

    def _symbolic_(self, R):
        """
        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = x^3 + 2/x
            sage: g = f._symbolic_(SR); g
            (x^4 + 2)/x
            sage: g(x=2)
            9

            sage: g = SR(f)
            sage: g(x=2)
            9

        Since :trac:`24072` the symbolic ring does not accept positive
        characteristic::

            sage: R.<w> = PuiseuxPolynomialRing(GF(7))
            sage: SR(2*w^3 + 1)
            Traceback (most recent call last):
            ...
            TypeError: positive characteristic not allowed in symbolic computations
        """
        d = {repr(g): R.var(g) for g in self._parent.gens()}
        return self.subs(**d)

    cpdef dict dict(self):
        """
        Return a dictionary representing ``self``.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: Q.<t> = PuiseuxPolynomialRing(R)
            sage: f = (x^3 + x*t^(3/2))^3 + t^2; f
            x^9 + 3*x^7*t^3/2 + t^2 + 3*x^5*t^3 + x^3*t^9/2
            sage: f.dict()
            {0: x^9, 3/2: 3*x^7, 2: 1, 3: 3*x^5, 9/2: x^3}
        """
        cdef dict d = self.__u.dict()
        return {Rational(k) / self.__n: d[k] for k in d}

    def coefficients(self):
        """
        Return the nonzero coefficients of ``self``.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: f = -5*t^(7/2) + t + t^2 - 10/3*t^3
            sage: f.coefficients()
            [1, 1, -10/3, -5]
        """
        return self.__u.coefficients()

    def exponents(self):
        """
        Return the exponents appearing in ``self`` with nonzero coefficients.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: f = -5*t^(7/2) + t + t^2 - 10/3*t^3
            sage: f.exponents()
            [1, 2, 3, 7/2]
        """
        return [Rational(i) / self.__n for i in self.__u.exponents()]

    def __setitem__(self, n, value):
        """
        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: f = t^2 + t^-3
            sage: f[2] = 5
            Traceback (most recent call last):
            ...
            IndexError: Puiseux polynomials are immutable
        """
        raise IndexError("Puiseux polynomials are immutable")

    cpdef _unsafe_mutate(self, i, value):
        r"""
        Sage assumes throughout that commutative ring elements are
        immutable. This is relevant for caching, etc. But sometimes you
        need to change a Puiseux polynomial and you really know what you're
        doing. That's when this function is for you.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: f = t^2 + t^-3
            sage: f._unsafe_mutate(2, 3)
            sage: f
            t^-3 + 3*t^2
        """
        j = i - self.__n
        if j >= 0:
            self.__u._unsafe_mutate(j, value)
        else: # off to the left
            if value != 0:
                self.__n = self.__n + j
                R = self._parent.base_ring()
                coeffs = [value] + [R.zero() for _ in range(1,-j)] + self.__u.list()
                self.__u = self.__u._parent(coeffs)

    cpdef _add_(self, right_m):
        """
        Add two Puiseux polynomials with the same parent.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: t + t
            2*t
            sage: f = 1 + t^2 + t^(3/2) - 17/3 * t^4
            sage: g = 2 + t^(3/5)
            sage: f + g
            3 + t^3/5 + t^3/2 + t^2 - 17/3*t^4
            sage: f + 0
            1 + t^3/2 + t^2 - 17/3*t^4
            sage: 0 + f
            1 + t^3/2 + t^2 - 17/3*t^4
            sage: R(0) + R(0)
            0
            sage: t^(1/2) + t^(1/3)
            t^1/3 + t^1/2

        ALGORITHM: reduction to common denominator of exponents
        """
        cdef PuiseuxPolynomial_univariate right = <PuiseuxPolynomial_univariate>right_m
        cdef long m
        cdef PuiseuxPolynomial_univariate ret

        # 1. Special case when one or the other is 0.
        if not right:
            return self
        if not self:
            return right

        m = lcm(self.__n, right.__n)

        t = self.__u.parent().gen()
        f1 = self.__u.subs({t: t**(m // self.__n)})
        f2 = right.__u.subs({t: t**(m // right.__n)})

        # 3. Add
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> (f1 + f2)
        ret.__n = m
        return ret

    cpdef _sub_(self, right_m):
        """
        Subtract two Puiseux polynomials with the same parent.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: t - t
            0
            sage: t^(3/5) + 2 * t^(4/7)
            2*t^4/7 + t^3/5

        ALGORITHM: reduction to common denominator of exponents
        """
        cdef PuiseuxPolynomial_univariate right = <PuiseuxPolynomial_univariate>right_m
        cdef long m
        cdef PuiseuxPolynomial_univariate ret

        # 1. Special case when one or the other is 0.
        if not right:
            return self
        if not self:
            return -right

        m = lcm(self.__n, right.__n)

        t = self.__u.parent().gen()
        f1 = self.__u.subs({t: t**(m // self.__n)})
        f2 = right.__u.subs({t: t**(m // right.__n)})

        # 3. Subtract
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> (f1 - f2)
        ret.__n = m
        return ret

    def degree(self):
        """
        Return the degree of ``self``.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(ZZ)
            sage: g = x^2 - x^4
            sage: g.degree()
            4
            sage: g = -10*x^(121/5) + x^2 - x^7
            sage: g.degree()
            121/5
        """
        return self.__u.degree() / self.__n

    def __neg__(self):
        """
        Return the negative of ``self``.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(ZZ)
            sage: -(1+t^5)
            -1 - t^5
        """
        cdef PuiseuxPolynomial_univariate ret
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> -self.__u
        ret.__n = self.__n
        # No need to normalize
        return ret

    cpdef _mul_(self, right_r):
        """
        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(GF(2))
            sage: f = 1 + x + x^2 + 3*x^4
            sage: g = 1 - x + x^2 - x^(11/3)
            sage: f*g
            1 + x^2 + x^11/3 + x^14/3 + x^5 + x^17/3 + x^6 + x^23/3
        """
        cdef PuiseuxPolynomial_univariate right = <PuiseuxPolynomial_univariate>right_r

        m = lcm(self.__n, right.__n)

        t = self.__u.parent().gen()
        f1 = self.__u.subs({t: t**(m // self.__n)})
        f2 = right.__u.subs({t: t**(m // right.__n)})

        cdef PuiseuxPolynomial_univariate ret
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> (f1 * f2)
        ret.__n = m
        return ret

    cpdef _rmul_(self, Element c):
        """
        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(ZZ)
            sage: f = 1 + x + x^(4/3) + 3*x^4
            sage: 3 * f
            3 + 3*x + 3*x^4/3 + 9*x^4
        """
        cdef PuiseuxPolynomial_univariate ret
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> self.__u._rmul_(c)
        ret.__n = self.__n
        return ret

    cpdef _lmul_(self, Element c):
        """
        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(ZZ)
            sage: f = 1 + x + x^(4/3) + 3*x^4
            sage: f * 3
            3 + 3*x + 3*x^4/3 + 9*x^4
        """
        cdef PuiseuxPolynomial_univariate ret
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> self.__u._lmul_(c)
        ret.__n = self.__n
        return ret

    def is_monomial(self):
        r"""
        Return ``True`` if ``self`` is a monomial; that is, if ``self``
        is `x^n` for some integer `n`.

        EXAMPLES::

            sage: k.<z> = PuiseuxPolynomialRing(QQ)
            sage: z.is_monomial()
            True
            sage: k(1).is_monomial()
            True
            sage: (z+1).is_monomial()
            False
        """
        return self.__u.is_monomial()

    def __pow__(_self, r, dummy):
        """
        EXAMPLES::

            sage: x = PuiseuxPolynomialRing(QQ,'x').0
            sage: f = x + x^2 + 3*x^(1/2)
            sage: g = 1 - x^(1/3)
            sage: f^3
            27*x^3/2 + 27*x^2 + 9*x^5/2 + 28*x^3 + 18*x^7/2 + 3*x^4 + 9*x^9/2 + 3*x^5 + x^6

            sage: g^4
            1 - 4*x^1/3 + 6*x^2/3 - 4*x + x^4/3

            sage: (x**(1/3))**(11/7)
            x^11/21

            sage: (1+x)**(1/2)
            Traceback (most recent call last):
            ...
            ValueError: exponent must be an integer
        """
        cdef PuiseuxPolynomial_univariate self = _self
        cdef Rational right = Rational(r)
        num, denom = right.numerator(), right.denominator()
        if right != r:
            raise ValueError("exponent must be a rational")
        if self.__u.is_term():
            return self._parent.element_class(self._parent, self.__u**num,
                                              self.__n * denom)
        else:
            if denom != 1:
                raise ValueError("exponent must be an integer")
            return self._parent.element_class(self._parent, self.__u**num,
                                              self.__n)

    cpdef _floordiv_(self, rhs):
        """
        Perform division with remainder and return the quotient.

        EXAMPLES::

            sage: L.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = x^3 + x^-3
            sage: g = x^-1 + x
            sage: f // g
            x^-2 - 1 + x^2
            sage: g * (f // g) == f
            True
            sage: f // 1
            x^-3 + x^3
            sage: 1 // f
            0
        """
        cdef PuiseuxPolynomial_univariate right = <PuiseuxPolynomial_univariate> rhs
        cdef PuiseuxPolynomial_univariate ret
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> (self.__u // right.__u)
        ret.__n = self.__n - right.__n
        return ret

    def shift(self, k):
        r"""
        Return this Puiseux polynomial multiplied by the power `t^k`.

        Does not change this polynomial.

        TO DO : k could be any rational..

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ['y'])
            sage: f = (t+t^(1/2))^4; f
            t^2 + 4*t^5/2 + 6*t^3 + 4*t^7/2 + t^4

            sage: f.shift(10)
            t^6 + 4*t^8 + 6*t^10 + 4*t^12 + t^14
            sage: f >> 10
            t^-14 + 4*t^-12 + 6*t^-10 + 4*t^-8 + t^-6
            sage: f << 4
            1 + 4*t^2 + 6*t^4 + 4*t^6 + t^8
        """
        cdef PuiseuxPolynomial_univariate ret
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        t = self.__u.parent().gen()
        f1 = self.__u.subs({t: t**1})
        ret.__u = f1
        ret.__n = self.__n
        # No need to normalize
        return ret

    def __lshift__(PuiseuxPolynomial_univariate self, k):
        """
        Return the left shift of ``self``.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: f = (t+t^-1)^4; f
            t^-4 + 4*t^-2 + 6 + 4*t^2 + t^4
            sage: f << 4
            1 + 4*t^2 + 6*t^4 + 4*t^6 + t^8
        """
        cdef PuiseuxPolynomial_univariate ret
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = self.__u
        ret.__n = self.__n + k
        # No need to normalize
        return ret

    def __rshift__(PuiseuxPolynomial_univariate self, k):
        """
        Return the right shift of ``self``.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: f = (t+t^-1)^4; f
            t^-4 + 4*t^-2 + 6 + 4*t^2 + t^4
            sage: f >> 10
            t^-14 + 4*t^-12 + 6*t^-10 + 4*t^-8 + t^-6
        """
        cdef PuiseuxPolynomial_univariate ret
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = self.__u
        ret.__n = self.__n - k
        # No need to normalize
        return ret

    cpdef _div_(self, rhs):
        """
        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = x + x^2 + 3*x^4
            sage: g = 1/x^7 - x + x^2 - x^4
            sage: f / x
            1 + x + 3*x^3
            sage: f / g
            (3*x^11 + x^9 + x^8)/(-x^11 + x^9 - x^8 + 1)
            sage: (x^-2 + x)*(x^-2 + 1) / ((x^5 + x^8)*(x + 2))
            (x^2 + 1)/(x^10 + 2*x^9)
            sage: (x^-2 + x)*(x^-2 + 1) / ((x^-5 + x^-8)*(x + 2))
            (x^6 + x^4)/(x + 2)
        """
        cdef PuiseuxPolynomial_univariate right = <PuiseuxPolynomial_univariate> rhs
        if right.__u.is_zero():
            raise ZeroDivisionError
        return self * ~right

    def __invert__(self):
        """
        Return the inverse of ``self``.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: i = ~(t^-2); i
            t^2
            sage: i.parent() is R
            True
            sage: i = ~(2*t^2); i
            1/2*t^-2
            sage: i.parent() is R
            True
            sage: i = ~(t^-2 + 2 + t^2); i
            t^2/(t^4 + 2*t^2 + 1)
            sage: i.parent()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        cdef PuiseuxPolynomial_univariate ret
        if self.__u.is_constant(): # this has a single term c*x^n
            ret = <PuiseuxPolynomial_univariate> self._new_c()
            if self.__u.is_unit():
                ret.__u = self.__u.inverse_of_unit()
                ret.__n = -self.__n
                return ret
            # Enlarge the ring so we can divide by the coefficient
            R = self._parent.base_ring().fraction_field()
            P = self._parent.change_ring(R)
            return P.element_class(P, ~R(self.__u), -self.__n)
        P = self._parent._R
        if self.__n < 0:
            return P.gen()**-self.__n / self.__u
        return P.one() / (P.gen()**self.__n * self.__u)

    def inverse_of_unit(self):
        """
        Return the inverse of ``self`` if a unit.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: (t^-2).inverse_of_unit()
            t^2
            sage: (t + 2).inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: element is not a unit
        """
        if self.is_unit():
            return ~self
        raise ArithmeticError("element is not a unit")

    def _fraction_pair(self):
        """
        Return one representation of ``self`` as a pair
        ``(numerator, denominator)``.

        Here both the numerator and the denominator are polynomials.

        This is used for coercion into the fraction field.

        EXAMPLES::

            sage: L.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = 4*x^-7 + 3*x^3 + 1 + 2*x^4 + x^6
            sage: f._fraction_pair()
            (x^13 + 2*x^11 + 3*x^10 + x^7 + 4, x^7)
        """
        P = self._parent._R
        numer = self.__u
        denom = P.one()
        return (numer, denom)

    def gcd(self, right):
        """
        Return the gcd of ``self`` with ``right`` where the common divisor
        ``d`` makes both ``self`` and ``right`` into polynomials with
        the lowest possible degree.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: t.gcd(2)
            1
            sage: gcd(t^-2 + 1, t^-4 + 3*t^-1)
            t^-4
            sage: gcd((t^-2 + t)*(t + t^-1), (t^5 + t^8)*(1 + t^-2))
            t^-3 + t^-1 + 1 + t^2
        """
        b = <PuiseuxPolynomial_univariate> self._parent(right)
        cdef PuiseuxPolynomial_univariate ret
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = self.__u.gcd(b.__u)
        ret.__n = min(self.__n, b.__n)
        return ret

    @coerce_binop
    def quo_rem(self, right_r):
        """
        Attempts to divide ``self`` by ``right`` and returns a quotient and
        a remainder.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: (t^-3 - t^3).quo_rem(t^-1 - t)
            (t^-2 + 1 + t^2, 0)
            sage: (t^-2 + 3 + t).quo_rem(t^-4)
            (t^2 + 3*t^4 + t^5, 0)
            sage: (t^-2 + 3 + t).quo_rem(t^-4 + t)
            (0, 1 + 3*t^2 + t^3)
        """
        cdef PuiseuxPolynomial_univariate right = <PuiseuxPolynomial_univariate> right_r
        q,r = self.__u.quo_rem(right.__u)
        cdef PuiseuxPolynomial_univariate ql, qr
        ql = <PuiseuxPolynomial_univariate> self._new_c()
        ql.__u = <ModuleElement> q
        ql.__n = self.__n - right.__n
        qr = <PuiseuxPolynomial_univariate> self._new_c()
        qr.__u = <ModuleElement> r
        qr.__n = 0
        return (ql, qr)

    cpdef _richcmp_(self, right_r, int op):
        r"""
        Comparison of ``self`` and ``right_r``.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = x^(1/2) + 1 + x
            sage: g = x^(1/2) + 1
            sage: f == g
            False

        ::

            sage: f = x^(1/2) + 1 + x
            sage: g = x^(1/3) + 2
            sage: f == g
            False
            sage: f != g
            True
            sage: f < g
            True
            sage: f <= g
            True
            sage: f > g
            False
            sage: f >= g
            False

        ::

            sage: f = x^(1/4) + 1 + x
            sage: g = x^(1/2) + 2
            sage: f == g
            False
            sage: f < g
            False
            sage: f > g
            True
        """
        cdef PuiseuxPolynomial_univariate right = <PuiseuxPolynomial_univariate> right_r

        zero = self._parent.base_ring().zero()

        if not self and not right:
            return rich_to_bool(op, 0)

        # zero pad coefficients on the left, to line them up for comparison
        cdef long n = min(self.__n, right.__n)
        x = [zero] * (self.__n - n) + self.__u.list()
        y = [zero] * (right.__n - n) + right.__u.list()

        # zero pad on right to make the lists the same length
        # (this is necessary since the power series list() function just
        # returns the coefficients of the underlying polynomial, which may
        # have zeroes in the high coefficients)
        if len(x) < len(y):
            x.extend([zero] * (len(y) - len(x)))
        elif len(y) < len(x):
            y.extend([zero] * (len(x) - len(y)))

        return richcmp(x, y, op)

    def valuation(self, p=None):
        """
        Return the valuation of ``self``.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(ZZ)
            sage: f = x**(1/2) + x^2 + 3*x^4
            sage: g = 1 - x + x^(3/7) - x^4
            sage: f.valuation()
            1/2
            sage: g.valuation()
            0
        """
        return self.__u.valuation(p) / self.__n

    def truncate(self, n):
        """
        Return a polynomial with degree at most `n-1` whose `j`-th coefficients
        agree with ``self`` for all `j < n`.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = 1/x^12 + x^3 + x^5 + x^9
            sage: f.truncate(10)
            x^-12 + x^3 + x^5 + x^9
            sage: f.truncate(5)
            x^-12 + x^3
            sage: f.truncate(-16)
            0
        """
        if n <= self.valuation():
            return self._parent.zero()
        cdef PuiseuxPolynomial_univariate ret
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> self.__u.truncate(n - self.__n)
        ret.__n = self.__n
        return ret

    def variable_name(self):
        """
        Return the name of variable of ``self`` as a string.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = 1/x + x^2 + 3*x^4
            sage: f.variable_name()
            'x'
        """
        return self._parent.variable_name()

    def variables(self):
        """
        Return the tuple of variables occuring in this Puiseux polynomial.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = 1/x + x^2 + 3*x^4
            sage: f.variables()
            (x,)
            sage: R.one().variables()
            ()
        """
        if self.is_constant():
            return ()
        return self._parent.gens()

    def polynomial_construction(self):
        """
        Return the polynomial and the denominator of power used to construct the
        Puiseux polynomial `u(t^{1/n})`.

        OUTPUT:

        A tuple ``(u, n)`` where ``u`` is the underlying polynomial and ``n``
        is the power of the exponent denominator.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = 1 + x^2 + 3*x^(7/11)
            sage: f.polynomial_construction()
            (x^22 + 3*x^7 + 1, 11)
        """
        return (self.__u, self.__n)

    def is_constant(self):
        """
        Return whether this Puiseux polynomial is constant.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: x.is_constant()
            False
            sage: R.one().is_constant()
            True
            sage: (x^(2/3)).is_constant()
            False
            sage: (x^2).is_constant()
            False
            sage: (x^2 + 2).is_constant()
            False
            sage: R(0).is_constant()
            True
            sage: R(42).is_constant()
            True
            sage: x.is_constant()
            False
        """
        return self.__u.is_constant()

    def __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = 1 + x^2 + 3*x^(4/3)
            sage: cf = copy(f)
            sage: cf == f
            True
            sage: cf is not f
            True
        """
        from copy import copy
        cdef PuiseuxPolynomial_univariate ret
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = copy(self.__u)
        ret.__n = self.__n
        # No need to normalize
        return ret

    def derivative(self, *args):
        """
        The formal derivative of this Puiseux polynomial, with respect
        to variables supplied in args.

        Multiple variables and iteration counts may be supplied. See
        documentation for the global :func:`derivative` function for more
        details.

        .. SEEALSO::

           :meth:`_derivative`

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: g = x^10 - x + x^2 - x^(4/5)
            sage: g.derivative()
            -10*x^-11 - 1 + 2*x - 4*x^3
            sage: g.derivative(x)
            -10*x^-11 - 1 + 2*x - 4*x^3

        ::

            sage: R.<t> = PolynomialRing(ZZ)
            sage: S.<x> = PuiseuxPolynomialRing(R)
            sage: f = 2*t + (3*t^2 + 6*t)*x**(1/2)
            sage: f.derivative()
            -2*t*x^-2 + (3*t^2 + 6*t)
            sage: f.derivative(x)
            -2*t*x^-2 + (3*t^2 + 6*t)
            sage: f.derivative(t)
            2*x^-1 + (6*t + 6)*x
        """
        return multi_derivative(self, args)

    def _derivative(self, var=None):
        """
        The formal derivative of this Puiseux series with respect to ``var``.

        If ``var`` is ``None`` or the generator of this ring, it's the formal
        derivative as expected. Otherwise, ``_derivative(var)`` gets called
        recursively on each coefficient.

        .. SEEALSO::

           :meth:`derivative`

        EXAMPLES::

            sage: R.<x> = PuiseuxPolynomialRing(QQ)
            sage: f = x^2 + 3*x^(6/5)
            sage: f._derivative()
            18/5*x^1/5 + 2*x
            sage: f._derivative(x)
            18/5*x^1/5 + 2*x

        Differentiating with respect to something other than the generator
        gets recursed into the base ring::

            sage: R.<t> = PolynomialRing(ZZ)
            sage: S.<x> = PuiseuxPolynomialRing(R)
            sage: f = 2*t + (3*t^2 + 6*t)*x
            sage: f._derivative(t)
            2 + (6*t + 6)*x
        """
        cdef PuiseuxPolynomial_univariate ret
        if var is not None and var is not self._parent.gen():
            # call _derivative() recursively on coefficients
            u = [coeff._derivative(var) for coeff in self.__u.list(copy=False)]
            ret = <PuiseuxPolynomial_univariate> self._new_c()
            ret.__u = <ModuleElement> self._parent._R(u)
            ret.__n = self.__n
            return ret

        # compute formal derivative with respect to generator
        if self.is_zero():
            return self  # this is already 0
        cdef long n = self.__n
        cdef dict a = self.__u.dict()
        cdef dict b = {}
        for exp, coeff in a.items():
            b[exp - n] = coeff * exp / n
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> self._parent._R(b)
        ret.__n = self.__n
        return ret

    def integral(self):
        r"""
        The formal integral of this Puiseux series with 0 constant term.

        EXAMPLES:

        The integral may or may not be defined if the base ring
        is not a field.

        ::

            sage: t = PuiseuxPolynomialRing(QQ, 't').0
            sage: f = 2*t^(3/4) + 3*t^2
            sage: f.integral()
            8/7*t^7/4 + t^3

        ::

            sage: f = t^3
            sage: f.integral()
            Traceback (most recent call last):
            ...
            ArithmeticError: coefficients of integral cannot be coerced into the base ring

        The integral of `1/t` is `\log(t)`, which is not given by a
        Puiseux polynomial::

            sage: t = PuiseuxPolynomialRing(ZZ,'t').0
            sage: f = -1/t^3 - 31/t
            sage: f.integral()
            Traceback (most recent call last):
            ...
            ArithmeticError: the integral of is not a Puiseux polynomial, since t^-1 has nonzero coefficient

        Another example with just one negative coefficient::

            sage: A.<t> = PuiseuxPolynomialRing(QQ)
            sage: f = -2*t^(-4)
            sage: f.integral()
            2/3*t^-3
            sage: f.integral().derivative() == f
            True
        """
        cdef long i, n = self.__n
        cdef PuiseuxPolynomial_univariate ret
        cdef dict a = self.__u.dict()
        cdef dict b = {}
        for exp, coeff in a.items():
            b[exp + n] = coeff * n / (exp + n)
        try:
            u = self._parent._R(b)
        except TypeError:
            raise ArithmeticError("coefficients of integral cannot be coerced into the base ring")
        ret = <PuiseuxPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> u
        ret.__n = n
        return ret

    def __call__(self, *x, **kwds):
        """
        Compute value of this Puiseux polynomial at ``x``.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(ZZ)
            sage: f = t^(-2) + t^2
            sage: f(2)
            17/4
            sage: f(-1)
            2
            sage: f(1/3)
            82/9
            sage: f(t=-1)
            2
            sage: f(x=-1)
            t^-2 + t^2
            sage: f()
            t^-2 + t^2
            sage: f(1,2)
            Traceback (most recent call last):
            ...
            TypeError: number of arguments does not match number of
             variables in parent
        """
        if kwds:
            f = self.subs(**kwds)
            if x: # If there are non-keyword arguments
                return f(*x)
            else:
                return f

        if not x:
            return self
        if len(x) != 1:
            raise TypeError("number of arguments does not match number"
                            " of variables in parent")
        if isinstance(x[0], tuple):
            x = x[0]
        return self.__u(x) * (x[0]**self.__n)

    def factor(self):
        """
        Return a Puiseux monomial (the unit part of the factorization) and
        a factored polynomial.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(ZZ)
            sage: f = 4*t^-7 + 3*t^3 + 2*t^4 + t^-6
            sage: f.factor()
            (t^-7) * (4 + t + 3*t^10 + 2*t^11)
        """
        cdef PuiseuxPolynomial_univariate u, d
        pf = self.__u.factor()
        u = <PuiseuxPolynomial_univariate> self._new_c()
        u.__u = pf.unit()
        u.__n = self.__n

        f = []
        for t in pf:
            d = <PuiseuxPolynomial_univariate> self._new_c()
            d.__u = t[0]
            d.__n = 0
            if d.is_unit():
                u *= d ** t[1]
            else:
                f.append((d, t[1]))

        return Factorization(f, unit=u)

    def residue(self):
        """
        Return the residue of ``self``.

        The residue is the coefficient of `t^-1`.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: f = 3*t^-2 - t^-1 + 3 + t^2
            sage: f.residue()
            -1
            sage: g = -2*t^-2 + 4 + 3*t
            sage: g.residue()
            0
            sage: f.residue().parent()
            Rational Field
        """
        return self.__u[-1 - self.__n]

    def constant_coefficient(self):
        """
        Return the coefficient of the constant term of ``self``.

        EXAMPLES::

            sage: R.<t> = PuiseuxPolynomialRing(QQ)
            sage: f = 3*t^-2 - t^-1 + 3 + t^2
            sage: f.constant_coefficient()
            3
            sage: g = -2*t^-2 + t^-1 + 3*t
            sage: g.constant_coefficient()
            0
        """
        return self.__u[-self.__n]



