"""
Power series implemented using PARI

AUTHORS:

- Peter Bruin (December 2013): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013 Peter Bruin <P.Bruin@warwick.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/stdsage.pxi"
include "sage/libs/pari/decl.pxi"

from sage.libs.pari.gen cimport gen as pari_gen
from sage.libs.pari.pari_instance cimport PariInstance
import sage.libs.pari.all
cdef PariInstance pari = sage.libs.pari.all.pari

from sage.rings.polynomial.polynomial_element cimport Polynomial
from power_series_ring_element cimport PowerSeries
from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.structure.parent cimport Parent
from infinity import infinity

cdef class PowerSeries_pari(PowerSeries):
    """
    A power series.

    """
    def __init__(self, parent, f=0, prec=infinity, check=True):
        """
        TEST::

            sage: R.<q> = PowerSeriesRing(CC, implementation='pari')
            sage: loads(q.dumps()) == q
            True
            sage: TestSuite(q).run()

        """
        cdef Parent f_parent
        cdef pari_gen g
        cdef long t
        cdef str v = parent.variable_name()

        if not check and PY_TYPE_CHECK(f, pari_gen):
            # Fast construction for PARI objects of suitable type
            # (series, polynomials, scalars and rational functions).
            # We ignore the precision argument and use the precision
            # of f, or the default precision if f is a rational
            # function.
            g = <pari_gen>f
            t = typ(g.g)
            if t == t_SER and varn(g.g) == pari.get_var(v):
                prec = lg(g.g) - 2 + valp(g.g)
            elif t == t_RFRAC:
                prec = parent.default_prec()
                g = g.Ser(v, prec - g.valuation(v))
            else:
                prec = infinity
            self.g = g
            PowerSeries.__init__(self, parent, prec)
            return

        R = parent.base_ring()
        P = parent._poly_ring()

        if PY_TYPE_CHECK(f, PowerSeries):  # not only PowerSeries_pari
            f_parent = (<PowerSeries>f)._parent
            if f_parent is parent:
                if prec is infinity:
                    prec = (<PowerSeries>f)._prec
                g = f._pari_()
            elif R.has_coerce_map_from(f_parent):
                g = R.coerce(f)._pari_()
            else:
                if prec is infinity:
                    prec = f.prec()
                g = f.polynomial().change_ring(R)._pari_()
        elif PY_TYPE_CHECK(f, Polynomial):
            f_parent = (<Polynomial>f)._parent
            if f_parent is P:
                g = f._pari_()
            elif R.has_coerce_map_from(f_parent):
                g = R.coerce(f)._pari_()
            else:
                g = P.coerce(f)._pari_()
        elif PY_TYPE_CHECK(f, pari_gen):
            g = f
            t = typ(g.g)
            if t == t_POL:
                g = P(g)._pari_()
            elif t == t_SER and varn(g.g) == pari.get_var(v):
                if valp(g.g) < 0:
                    raise ValueError('series has negative valuation')
                if prec is infinity:
                    prec = lg(g.g) - 2 + valp(g.g)
                g = P(g.Pol(v))._pari_()
            elif t == t_RFRAC:
                if prec is infinity:
                    prec = parent.default_prec()
                g = P.fraction_field()(g)._pari_()
                g = g.Ser(v, prec - g.valuation(v))
            elif t == t_VEC:
                g = P(g.Polrev(v))._pari_()
            else:
                g = R(g)._pari_()
        elif PY_TYPE_CHECK(f, list) or PY_TYPE_CHECK(f, tuple):
            g = pari([R.coerce(x) for x in f]).Polrev(v)
        else:
            g = R.coerce(f)._pari_()

        if prec is infinity:
            self.g = g
        else:
            if not g:
                self.g = g.Ser(v, prec)
            else:
                self.g = g.Ser(v, prec - g.valuation(v))

        PowerSeries.__init__(self, parent, prec)

    def __hash__(self):
        """
        Return a hash of ``self``.

        TEST::

            sage: R.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: hash(t^2 + 1) == hash(pari(t^2 + 1))
            True

        """
        return hash(self.g)

    def __reduce__(self):
        """
        Used for pickling.

        EXAMPLES::

            sage: A.<z> = PowerSeriesRing(RR, implementation='pari')
            sage: f = z - z^3 + O(z^10)
            sage: f == loads(dumps(f)) # indirect doctest
            True

        """
        return PowerSeries_pari, (self._parent, self.g, self._prec, False)

    def __richcmp__(left, right, int op):
       """
       Used for comparing power series.

       EXAMPLES::

           sage: R.<t> = PowerSeriesRing(ZZ, implementation='pari')
           sage: f = 1 + t + t^7 - 5*t^10
           sage: g = 1 + t + t^7 - 5*t^10 + O(t^15)
           sage: f == f
           True
           sage: f < g
           False
           sage: f == g
           True

       """
       return (<Element>left)._richcmp(right, op)

    def _pari_(self):
        """
        Convert ``self`` to a PARI object.

        TEST::

            sage: R.<t> = PowerSeriesRing(GF(7), implementation='pari')
            sage: (3 - t^3 + O(t^5))._pari_()
            Mod(3, 7) + Mod(6, 7)*t^3 + O(t^5)

        """
        return self.g

    def polynomial(self):
        """
        Convert ``self`` to a polynomial.

        EXAMPLE::

            sage: R.<t> = PowerSeriesRing(GF(7), implementation='pari')
            sage: f = 3 - t^3 + O(t^5)
            sage: f.polynomial()
            6*t^3 + 3

        """
        return self._parent._poly_ring()(self.list())

    def valuation(self):
        """
        Return the valuation of ``self``.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ, implementation='pari')
            sage: (5 - t^8 + O(t^11)).valuation()
            0
            sage: (-t^8 + O(t^11)).valuation()
            8
            sage: O(t^7).valuation()
            7
            sage: R(0).valuation()
            +Infinity

        """
        if not self.g:
            return self._prec
        return self.g.valuation(self._parent.variable_name())

    def __nonzero__(self):
        """
        Return True if ``self`` is nonzero, and False otherwise.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(GF(11), implementation='pari')
            sage: (1 + t + O(t^18)).__nonzero__()
            True
            sage: R(0).__nonzero__()
            False
            sage: O(t^18).__nonzero__()
            False

        """
        return self.g.__nonzero__()

    def __call__(self, *x, **kwds):
        """
        Evaluate ``self`` at `x = a`.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: f = t^2 + t^3 + O(t^6)
            sage: f(t^3)
            t^6 + t^9 + O(t^18)
            sage: f(t=t^3)
            t^6 + t^9 + O(t^18)
            sage: f(f)
            t^4 + 2*t^5 + 2*t^6 + 3*t^7 + O(t^8)
            sage: f(f)(f) == f(f(f))
            True

        The following demonstrates that the problems raised in :trac:`3979`
        and :trac:`5367` are solved::

            sage: [f(t^2 + O(t^n)) for n in [9, 10, 11]]
            [t^4 + t^6 + O(t^11), t^4 + t^6 + O(t^12), t^4 + t^6 + O(t^12)]
            sage: f(t^2)
            t^4 + t^6 + O(t^12)

        It is possible to substitute a series for which only the precision
        is defined::

            sage: f(O(t^5))
            O(t^10)

        or to substitute a polynomial (the result belonging to the power
        series ring over the same base ring)::

            sage: P.<z> = ZZ[]
            sage: g = f(z + z^3); g
            z^2 + z^3 + 2*z^4 + 3*z^5 + O(z^6)
            sage: g.parent()
            Power Series Ring in z over Integer Ring

        A series defined over another ring can be substituted::

            sage: S.<u> = PowerSeriesRing(GF(7), implementation='pari')
            sage: f(2*u + u^3 + O(u^5))
            4*u^2 + u^3 + 4*u^4 + 5*u^5 + O(u^6)

        Substituting `p`-adic numbers::

            sage: f(100 + O(5^7))
            5^4 + 3*5^5 + 4*5^6 + 2*5^7 + 2*5^8 + O(5^9)

            sage: ff = PowerSeriesRing(pAdicRing(5), 't', implementation='pari')(f)
            sage: ff
            (1 + O(5^20))*t^2 + (1 + O(5^20))*t^3 + O(t^6)

            sage: ff(100 + O(5^7))
            5^4 + 3*5^5 + 4*5^6 + 2*5^7 + 2*5^8 + O(5^9)

            sage: ff(100 + O(2^7))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: '5-adic Ring with capped relative precision 20' and '2-adic Ring with capped relative precision 20'

        The argument must have valuation at least 1, unless the series
        is actually a polynomial::

            sage: f(0)
            0
            sage: f(1 + t)
            Traceback (most recent call last):
            ...
            PariError: non positive valuation in a series substitution

            sage: f(t^-2)
            Traceback (most recent call last):
            ...
            PariError: non positive valuation in a series substitution

            sage: f(2 + O(5^3))
            Traceback (most recent call last):
            ...
            ValueError: can only substitute elements of positive valuation

            sage: g = t^2 + t^3
            sage: g(1 + t + O(t^2))
            2 + 5*t + O(t^2)
            sage: g(3)
            36

        Substitution of variables belonging to the base ring can be
        done using keywords::

            sage: P.<a> = GF(5)[]
            sage: Q.<x> = PowerSeriesRing(P, implementation='pari')
            sage: h = (1 - a*x)^-1 + O(x^7); h
            1 + a*x + a^2*x^2 + a^3*x^3 + a^4*x^4 + a^5*x^5 + a^6*x^6 + O(x^7)
            sage: h(x^2, a=3)
            1 + 3*x^2 + 4*x^4 + 2*x^6 + x^8 + 3*x^10 + 4*x^12 + O(x^14)

        """
        if len(kwds) >= 1:
            name = self._parent.variable_name()
            if kwds.has_key(name):  # the series variable is specified by a keyword
                if len(x) > 0:
                    raise ValueError("must not specify %s keyword and positional argument" % name)
                x = [kwds[name]]
                del kwds[name]

        if len(x) != 1:
            raise ValueError("must specify exactly one positional argument")

        a = x[0]

        s = self._prec
        if s is infinity:
            return self.polynomial()(a)

        # Determine the parent of the result.
        P = self._parent
        Q = a.parent()
        if not Q.has_coerce_map_from(P.base_ring()):
            from sage.structure.element import canonical_coercion
            a = canonical_coercion(P.base_ring()(0), a)[1]
            Q = a.parent()

        # The result is defined if the ring Q is complete with respect
        # to an ideal I, and the element a lies in I.  Here we only
        # implement a few special cases.
        from sage.rings.padics.padic_generic import pAdicGeneric
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        from sage.rings.power_series_ring import PowerSeriesRing_generic
        from sage.rings.laurent_series_ring import LaurentSeriesRing_generic
        if isinstance(Q, pAdicGeneric):
            # Substitution of p-adic numbers in power series is
            # currently not implemented in PARI (2.5.5).
            t = a.valuation()
            if t <= 0:
                raise ValueError("can only substitute elements of positive valuation")
            return Q(self.polynomial()(a)).add_bigoh(t * self._prec)
        elif isinstance(Q, (PowerSeriesRing_generic, LaurentSeriesRing_generic)):
            pass
        elif isinstance(Q, PolynomialRing_general):
            Q = Q.completion(Q.gen())
        elif Q.is_exact() and not a:
            pass
        else:
            raise ValueError('cannot substitute %s in %s' % (a, self))

        if len(kwds) == 0:
            return Q(self.g(a))
        else:
            kwds[P.variable_name()] = a
            return Q(self.g(**kwds))


    def __getitem__(self, n):
        """
        Return the `n`-th coefficient of self.

        If `n` is a slice object, this returns a power series of the
        same precision, whose coefficients are the same as ``self``
        for those indices in the slice, and 0 otherwise.

        Returns 0 for negative coefficients.  Raises an ``IndexError``
        if trying to access beyond known coefficients.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = 3/2 - 17/5*t^3 + O(t^5)
            sage: f[3]
            -17/5
            sage: f[-2]
            0
            sage: f[4]
            0
            sage: f[5]
            Traceback (most recent call last):
            ...
            IndexError: index out of bounds
            sage: f[1:4]
            -17/5*t^3 + O(t^5)

            sage: R.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: f = (2-t)^5; f
            32 - 80*t + 80*t^2 - 40*t^3 + 10*t^4 - t^5
            sage: f[2:4]
            80*t^2 - 40*t^3
            sage: f[5:9]
            -t^5
            sage: f[2:7:2]
            80*t^2 + 10*t^4
            sage: f[10:20]
            0
            sage: f[10:]
            0
            sage: f[:4]
            32 - 80*t + 80*t^2 - 40*t^3

            sage: f = 1 + t^3 - 4*t^4 + O(t^7) ; f
            1 + t^3 - 4*t^4 + O(t^7)
            sage: f[2:4]
            t^3 + O(t^7)
            sage: f[4:9]
            -4*t^4 + O(t^7)
            sage: f[2:7:2]
            -4*t^4 + O(t^7)
            sage: f[10:20]
            O(t^7)
            sage: f[10:]
            O(t^7)
            sage: f[:4]
            1 + t^3 + O(t^7)

        """
        cdef long t
        if isinstance(n, slice):
            # get values from slice object
            start = n.start if n.start is not None else 0
            stop = self._prec if n.stop is None else n.stop
            if stop is infinity:
                stop = self.polynomial().degree() + 1
            step = 1 if n.step is None else n.step

            # find corresponding polynomial
            poly = self.polynomial()[start:stop]
            if step is not None:
                coeffs = poly.padded_list(stop)
                for i in range(start, stop):
                    if (i - start) % step:
                        coeffs[i] = 0
                poly = self._parent._poly_ring()(coeffs)

            # return the power series
            return PowerSeries_pari(self._parent, poly,
                                    prec=self._prec, check=False)
        elif n < 0:
            return self.base_ring()(0)
        else:
            t = typ(self.g.g)
            if t == t_POL or t == t_SER:
                h = self.g[n]
            else:
                h = self.g
            return self.base_ring()(h)

    def __invert__(self):
        """
        Return the multiplicative inverse of ``self``.

        TEST::

            sage: R.<t> = PowerSeriesRing(QQ, default_prec=6, implementation='pari')
            sage: ~(R(1-t))
            1 + t + t^2 + t^3 + t^4 + t^5 + O(t^6)

        """
        cdef pari_gen h = ~self.g
        if h.valuation(self._parent.variable_name()) < 0:
            return self._parent.laurent_series_ring()(h)
        return PowerSeries_pari(self._parent, h, check=False)

    def __neg__(self):
        """
        Return the negative of ``self``.

        TEST::

            sage: R.<t> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = t + 17/5*t^3 + 2*t^4 + O(t^5)
            sage: -f
            -t - 17/5*t^3 - 2*t^4 + O(t^5)

        """
        return PowerSeries_pari(self._parent, -self.g, check=False)

    def __pow__(PowerSeries_pari self, n, m):
        """
        Exponentiation of power series.

        TESTS::

            sage: R.<t> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = 3 - t^3 + O(t^5)
            sage: a = f^3; a
            27 - 27*t^3 + O(t^5)
            sage: b = f^-3; b
            1/27 + 1/27*t^3 + O(t^5)

        """
        cdef pari_gen h = self.g ** n
        if h.valuation(self._parent.variable_name()) < 0:
            return self._parent.laurent_series_ring()(h)
        return PowerSeries_pari(self._parent, h, check=False)

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Addition of power series.

        TESTS::

            sage: R.<x> = PowerSeriesRing(ZZ, implementation='pari')
            sage: f = x^4 + O(x^5); f
            x^4 + O(x^5)
            sage: g = x^2 + O(x^3); g
            x^2 + O(x^3)
            sage: f+g
            x^2 + O(x^3)

        """
        return PowerSeries_pari(self._parent, self.g + (<PowerSeries_pari>right).g, check=False)

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Subtraction of power series.

        TEST::

            sage: k.<w> = ZZ[]
            sage: R.<t> = PowerSeriesRing(k, implementation='pari')
            sage: w*t^2 -w*t +13 - (w*t^2 + w*t)
            13 - 2*w*t

        """
        return PowerSeries_pari(self._parent, self.g - (<PowerSeries_pari>right).g, check=False)

    cpdef RingElement _mul_(self, RingElement right):
        """
        Multiplication of power series.

        TEST::

            sage: k.<w> = PowerSeriesRing(ZZ, implementation='pari')
            sage: (1+17*w+15*w^3+O(w^5))*(19*w^10+O(w^12))
            19*w^10 + 323*w^11 + O(w^12)

        """
        return PowerSeries_pari(self._parent, self.g * (<PowerSeries_pari>right).g, check=False)

    cpdef ModuleElement _rmul_(self, RingElement c):
        """
        Right multiplication by a scalar.

        TEST::

            sage: R.<t> = PowerSeriesRing(GF(7), implementation='pari')
            sage: f = t + 3*t^4 + O(t^11)
            sage: f * GF(7)(3)
            3*t + 2*t^4 + O(t^11)

        """
        return PowerSeries_pari(self._parent, self.g * c, check=False)

    cpdef ModuleElement _lmul_(self, RingElement c):
        """
        Left multiplication by a scalar.

        TEST::

            sage: R.<t> = PowerSeriesRing(GF(11), implementation='pari')
            sage: f = 1 + 3*t^4 + O(t^120)
            sage: 2 * f
            2 + 6*t^4 + O(t^120)

        """
        return PowerSeries_pari(self._parent, c * self.g, check=False)

    cpdef RingElement _div_(self, RingElement right):
        """
        Division of power series.

        TEST::

            sage: R.<t> = PowerSeriesRing(GF(11), default_prec=8, implementation='pari')
            sage: f = t/(1 - t); f
            t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + O(t^8)
            sage: f.parent()
            Power Series Ring in t over Finite Field of size 11
            sage: g = (1 - t)/t; g
            t^-1 + 10
            sage: g.parent()
            Laurent Series Ring in t over Finite Field of size 11

        """
        cdef pari_gen h = self.g / (<PowerSeries_pari>right).g
        if h.valuation(self._parent.variable_name()) < 0:
            return self._parent.laurent_series_ring()(h)
        return PowerSeries_pari(self._parent, h, check=False)

    def list(self):
        """
        Return the list of known coefficients for ``self``.

        This is just the list of coefficients of the underlying
        polynomial; it need not have length equal to ``self.prec()``.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: f = 1 - 5*t^3 + t^5 + O(t^7)
            sage: f.list()
            [1, 0, 0, -5, 0, 1]

            sage: S.<u> = PowerSeriesRing(pAdicRing(5), implementation='pari')
            sage: (2 + u).list()
            [2 + O(5^20), 1 + O(5^20)]

        """
        cdef pari_gen g = self.g
        cdef long vn = pari.get_var(self._parent.variable_name())
        R = self.base_ring()
        if typ(g.g) == t_SER and varn(g.g) == vn:
            g = g.truncate()
        if typ(g.g) == t_POL and varn(g.g) == vn:
            # t_POL has 2 codewords.  Use new_ref instead of g[i] for speed.
            return [R(pari.new_ref(gel(g.g, i), g)) for i in xrange(2, lg(g.g))]
        else:
            return [R(g)]

    def padded_list(self, n=None):
        """
        Return a list of coefficients of ``self`` up to (but not
        including) `q^n`.

        The list is padded with zeroes on the right so that it has
        length `n`.

        INPUT:

        - ``n`` - (optional) a non-negative integer.  If `n` is not
           given, it will be taken to be the precision of self, unless
           this is +Infinity, in which case we just return
           ``self.list()``.

        EXAMPLES::

            sage: R.<q> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = 1 - 17*q + 13*q^2 + 10*q^4 + O(q^7)
            sage: f.list()
            [1, -17, 13, 0, 10]
            sage: f.padded_list(7)
            [1, -17, 13, 0, 10, 0, 0]
            sage: f.padded_list(10)
            [1, -17, 13, 0, 10, 0, 0, 0, 0, 0]
            sage: f.padded_list(3)
            [1, -17, 13]
            sage: f.padded_list()
            [1, -17, 13, 0, 10, 0, 0]
            sage: g = 1 - 17*q + 13*q^2 + 10*q^4
            sage: g.list()
            [1, -17, 13, 0, 10]
            sage: g.padded_list()
            [1, -17, 13, 0, 10]
            sage: g.padded_list(10)
            [1, -17, 13, 0, 10, 0, 0, 0, 0, 0]

        """
        if n is None:
            if self._prec is infinity:
                return self.list()
            else:
                n = self._prec
        if not n:
            return []

        cdef pari_gen g = self.g
        cdef long l, m

        R = self.base_ring()
        if typ(g.g) == t_POL and varn(g.g) == pari.get_var(self._parent.variable_name()):
            l = lg(g.g) - 2  # t_POL has 2 codewords
            if n <= l:
                return [R(pari.new_ref(gel(g.g, i + 2), g)) for i in xrange(n)]
            else:
                return ([R(pari.new_ref(gel(g.g, i + 2), g)) for i in xrange(l)]
                        + [R.zero_element()] * (n - l))
        elif typ(g.g) == t_SER and varn(g.g) == pari.get_var(self._parent.variable_name()):
            l = lg(g.g) - 2  # t_SER has 2 codewords
            m = valp(g.g)
            if n <= m:
                return [R.zero_element()] * n
            elif n <= l + m:
                return ([R.zero_element()] * m
                        + [R(pari.new_ref(gel(g.g, i + 2), g)) for i in xrange(n - m)])
            else:
                return ([R.zero_element()] * m
                        + [R(pari.new_ref(gel(g.g, i + 2), g)) for i in xrange(l)]
                        + [R.zero_element()] * (n - l - m))
        else:
            return [R(g)] + [R.zero_element()] * (n - 1)

    def dict(self):
        """
        Return a dictionary of coefficients for ``self``.

        This is simply a dict for the underlying polynomial; it need
        not have keys corresponding to every number smaller than
        ``self.prec()``.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: f = 1 + t^10 + O(t^12)
            sage: f.dict()
            {0: 1, 10: 1}

        """
        return self.polynomial().dict()

    def _derivative(self, var=None):
        """
        Return the derivative of ``self`` with respect to the variable
        ``var``.

        If ``var`` is ``None``, the variable of the power series ring
        is used.

        .. seealso::

            :meth:`derivative()`

        EXAMPLES::

            sage: R.<w> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = 2 + 3*w^2 + w^10 + O(w^100); f
            2 + 3*w^2 + w^10 + O(w^100)
            sage: f._derivative()
            6*w + 10*w^9 + O(w^99)
            sage: f._derivative(w)
            6*w + 10*w^9 + O(w^99)

            sage: R.<t> = PolynomialRing(ZZ)
            sage: S.<x> = PowerSeriesRing(R, implementation='pari')
            sage: f = t^3*x^4 + O(x^5)
            sage: f._derivative()
            4*t^3*x^3 + O(x^4)
            sage: f._derivative(x)
            4*t^3*x^3 + O(x^4)
            sage: f._derivative(t)
            3*t^2*x^4 + O(x^5)

        """
        if var is None:
            var = self._parent.variable_name()
        return PowerSeries_pari(self._parent, self.g.deriv(var), check=False)

    def integral(self, var=None):
        """
        Return the formal integral of ``self``.

        By default, the integration variable is the variable of the
        power series.  Otherwise, the integration variable is the
        optional parameter ``var``.

        .. NOTE::

            The integral is always chosen so the constant term is 0.

        EXAMPLES::

            sage: k.<w> = PowerSeriesRing(QQ, implementation='pari')
            sage: (1+17*w+15*w^3+O(w^5)).integral()
            w + 17/2*w^2 + 15/4*w^4 + O(w^6)
            sage: (w^3 + 4*w^4 + O(w^7)).integral()
            1/4*w^4 + 4/5*w^5 + O(w^8)
            sage: (3*w^2).integral()
            w^3

        TESTS::

            sage: t = PowerSeriesRing(QQ, 't', implementation='pari').gen()
            sage: f = t + 5*t^2 + 21*t^3
            sage: g = f.integral() ; g
            1/2*t^2 + 5/3*t^3 + 21/4*t^4
            sage: g.parent()
            Power Series Ring in t over Rational Field

            sage: R.<a> = QQ[]
            sage: t = PowerSeriesRing(R, 't', implementation='pari').gen()
            sage: f = a*t +5*t^2
            sage: f.integral()
            1/2*a*t^2 + 5/3*t^3
            sage: f.integral(a)
            1/2*a^2*t + 5*a*t^2

        """
        if var is None:
            var = self._parent.variable_name()
        return PowerSeries_pari(self._parent, self.g.intformal(var), check=False)

    def reversion(self, precision=None):
        """
        Return the reversion of ``self``.

        The reversion of a power series `f` is the power series `g`
        such that `g(f(x)) = x`.  This exists if and only if the
        valuation of ``self`` is exactly 1 and the coefficient of `x`
        is a unit.

        If the optional argument ``precision`` is given, the reversion
        is returned with this precision.  If ``f`` has infinite
        precision and the argument ``precision`` is not given, then
        the reversion is returned with the default precision of
        ``f.parent()``.

        EXAMPLES::

            sage: R.<x> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = 2*x + 3*x^2 - x^4 + O(x^5)
            sage: g = f.reversion()
            sage: g
            1/2*x - 3/8*x^2 + 9/16*x^3 - 131/128*x^4 + O(x^5)
            sage: f(g)
            x + O(x^5)
            sage: g(f)
            x + O(x^5)

            sage: A.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: a = t - t^2 - 2*t^4 + t^5 + O(t^6)
            sage: b = a.reversion(); b
            t + t^2 + 2*t^3 + 7*t^4 + 25*t^5 + O(t^6)
            sage: a(b)
            t + O(t^6)
            sage: b(a)
            t + O(t^6)

            sage: B.<b,c> = PolynomialRing(ZZ)
            sage: A.<t> = PowerSeriesRing(B, implementation='pari')
            sage: f = t + b*t^2 + c*t^3 + O(t^4)
            sage: g = f.reversion(); g
            t - b*t^2 + (2*b^2 - c)*t^3 + O(t^4)
            sage: f(g)
            t + O(t^4)
            sage: g(f)
            t + O(t^4)

            sage: A.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: B.<x> = PowerSeriesRing(A, implementation='pari')
            sage: f = (1 - 3*t + 4*t^3 + O(t^4))*x + (2 + t + t^2 + O(t^3))*x^2 + O(x^3)
            sage: g = f.reversion(); g
            (1 + 3*t + 9*t^2 + 23*t^3 + O(t^4))*x + (-2 - 19*t - 118*t^2 + O(t^3))*x^2 + O(x^3)

        The optional argument ``precision`` sets the precision of the output::

            sage: R.<x> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = 2*x + 3*x^2 - 7*x^3 + x^4 + O(x^5)
            sage: g = f.reversion(precision=3); g
            1/2*x - 3/8*x^2 + O(x^3)
            sage: f(g)
            x + O(x^3)
            sage: g(f)
            x + O(x^3)

        If the input series has infinite precision, the precision of the
        output is automatically set to the default precision of the parent
        ring::

            sage: R.<x> = PowerSeriesRing(QQ, default_prec=20, implementation='pari')
            sage: (x - x^2).reversion()  # get some Catalan numbers
            x + x^2 + 2*x^3 + 5*x^4 + 14*x^5 + 42*x^6 + 132*x^7 + 429*x^8 + 1430*x^9 + 4862*x^10 + 16796*x^11 + 58786*x^12 + 208012*x^13 + 742900*x^14 + 2674440*x^15 + 9694845*x^16 + 35357670*x^17 + 129644790*x^18 + 477638700*x^19 + O(x^20)
            sage: (x - x^2).reversion(precision=3)
            x + x^2 + O(x^3)

        TESTS::

            sage: R.<x> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = 1 + 2*x + 3*x^2 - x^4 + O(x^5)
            sage: f.reversion()
            Traceback (most recent call last):
            ...
            PariError: valuation not equal to 1 in serreverse

        """
        cdef PowerSeries_pari f
        if self._prec is infinity:
            if precision is None:
                precision = self._parent.default_prec()
            f = self.add_bigoh(precision)
        else:
            if precision is None:
                precision = self._prec
            f = self
        return PowerSeries_pari(self._parent, f.g.serreverse(), precision, check=True)
