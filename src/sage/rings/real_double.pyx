r"""
Double Precision Real Numbers

EXAMPLES:

We create the real double vector space of dimension `3`::

    sage: V = RDF^3; V
    Vector space of dimension 3 over Real Double Field

Notice that this space is unique::

    sage: V is RDF^3
    True
    sage: V is FreeModule(RDF, 3)
    True
    sage: V is VectorSpace(RDF, 3)
    True

Also, you can instantly create a space of large dimension::

    sage: V = RDF^10000

TESTS:

Test NumPy conversions::

    sage: RDF(1).__array_interface__
    {'typestr': '=f8'}
    sage: import numpy
    sage: numpy.array([RDF.pi()]).dtype
    dtype('float64')
"""

# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport libc.math
from libc.string cimport memcpy
from cpython.object cimport *
from cpython.float cimport *

from cysignals.signals cimport sig_on, sig_off

from sage.ext.stdsage cimport PY_NEW
from sage.cpython.python_debug cimport if_Py_TRACE_REFS_then_PyObject_INIT

from sage.libs.gsl.all cimport *

gsl_set_error_handler_off()

import math, operator

import sage.rings.integer
import sage.rings.rational

from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ

from sage.categories.morphism cimport Morphism
from sage.structure.coerce cimport is_numpy_type
from sage.misc.randstate cimport randstate, current_randstate
from sage.structure.richcmp cimport rich_to_bool
from sage.arith.constants cimport *

cimport gmpy2


new_gen_from_real_double_element = None


def is_RealDoubleField(x):
    """
    Returns ``True`` if ``x`` is the field of real double precision numbers.

    This function is deprecated. Use :func:`isinstance` with
    :class:`~sage.rings.abc.RealDoubleField` instead.

    EXAMPLES::

        sage: from sage.rings.real_double import is_RealDoubleField
        sage: is_RealDoubleField(RDF)
        doctest:warning...
        DeprecationWarning: is_RealDoubleField is deprecated;
        use isinstance(..., sage.rings.abc.RealDoubleField) instead
        See https://trac.sagemath.org/32610 for details.
        True
        sage: is_RealDoubleField(RealField(53))
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(32610, 'is_RealDoubleField is deprecated; use isinstance(..., sage.rings.abc.RealDoubleField) instead')
    return isinstance(x, RealDoubleField_class)


cdef class RealDoubleField_class(sage.rings.abc.RealDoubleField):
    """
    An approximation to the field of real numbers using double
    precision floating point numbers. Answers derived from calculations
    in this approximation may differ from what they would be if those
    calculations were performed in the true field of real numbers. This
    is due to the rounding errors inherent to finite precision
    calculations.

    EXAMPLES::

        sage: RR == RDF
        False
        sage: RDF == RealDoubleField()    # RDF is the shorthand
        True

    ::

        sage: RDF(1)
        1.0
        sage: RDF(2/3)
        0.6666666666666666

    A ``TypeError`` is raised if the coercion doesn't make sense::

        sage: RDF(QQ['x'].0)
        Traceback (most recent call last):
        ...
        TypeError: cannot convert nonconstant polynomial
        sage: RDF(QQ['x'](3))
        3.0

    One can convert back and forth between double precision real
    numbers and higher-precision ones, though of course there may be
    loss of precision::

        sage: a = RealField(200)(2).sqrt(); a
        1.4142135623730950488016887242096980785696718753769480731767
        sage: b = RDF(a); b
        1.4142135623730951
        sage: a.parent()(b)
        1.4142135623730951454746218587388284504413604736328125000000
        sage: a.parent()(b) == b
        True
        sage: b == RR(a)
        True

    TESTS::

        sage: RDF.is_finite()
        False
    """
    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: R = RealDoubleField()
            sage: TestSuite(R).run()
        """
        from sage.categories.fields import Fields
        Field.__init__(self, self, category=Fields().Infinite().Metric().Complete())
        self._populate_coercion_lists_(init_no_parent=True,
                                       convert_method_name='_real_double_')

    _element_constructor_ = RealDoubleElement

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: loads(dumps(RDF)) is RDF
            True
        """
        return RealDoubleField, ()

    cpdef bint is_exact(self) except -2:
        """
        Returns ``False``, because doubles are not exact.

        EXAMPLES::

            sage: RDF.is_exact()
            False
        """
        return False

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(RDF) # indirect doctest
            \Bold{R}
        """
        return "\\Bold{R}"

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(RDF, verify=True)
            # Verified
            RDF
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: RDF._sage_input_(SageInputBuilder(), False)
            {atomic:RDF}
        """
        return sib.name('RDF')

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RealDoubleField() # indirect doctest
            Real Double Field
            sage: RDF
            Real Double Field
        """
        return "Real Double Field"

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: RDF._repr_option('element_is_atomic')
            True
        """
        if key == 'element_is_atomic':
            return True
        return super(RealDoubleField_class, self)._repr_option(key)

    def __richcmp__(self, x, op):
        """
        Compare ``self`` to ``x``.

        EXAMPLES::

            sage: RDF == 5
            False
            sage: loads(dumps(RDF)) == RDF
            True
        """
        if isinstance(x, RealDoubleField_class):
            return rich_to_bool(op, 0)
        if op == Py_NE:
            return True
        return NotImplemented

    def construction(self):
        r"""
        Returns the functorial construction of ``self``, namely, completion of
        the rational numbers with respect to the prime at `\infty`.

        Also preserves other information that makes this field unique (i.e.
        the Real Double Field).

        EXAMPLES::

            sage: c, S = RDF.construction(); S
            Rational Field
            sage: RDF == c(S)
            True
        """
        from sage.categories.pushout import CompletionFunctor
        return (CompletionFunctor(sage.rings.infinity.Infinity,
                                  53,
                                  {'type': 'RDF'}),
               sage.rings.rational_field.QQ)

    def complex_field(self):
        """
        Return the complex field with the same precision as ``self``, i.e.,
        the complex double field.

        EXAMPLES::

            sage: RDF.complex_field()
            Complex Double Field
        """
        from sage.rings.complex_double import CDF
        return CDF

    def algebraic_closure(self):
        """
        Return the algebraic closure of ``self``, i.e., the complex double
        field.

        EXAMPLES::

            sage: RDF.algebraic_closure()
            Complex Double Field
        """
        from sage.rings.complex_double import CDF
        return CDF

    cpdef _coerce_map_from_(self, S):
        """
        Canonical coercion of ``S`` to the real double field.

        The rings that canonically coerce to the real double field are:

        - the real double field itself
        - int, long, integer, and rational rings
        - numpy integers and floatings
        - the real lazy field
        - the MPFR real field with at least 53 bits of precision

        EXAMPLES::

            sage: RDF.coerce(5) # indirect doctest
            5.0
            sage: RDF.coerce(9499294r)
            9499294.0
            sage: RDF.coerce(61/3)
            20.333333333333332
            sage: parent(RDF(3) + CDF(5))
            Complex Double Field
            sage: parent(CDF(5) + RDF(3))
            Complex Double Field
            sage: CDF.gen(0) + 5.0
            5.0 + 1.0*I
            sage: RLF(2/3) + RDF(1)
            1.6666666666666665

            sage: import numpy
            sage: RDF.coerce(numpy.int8('1'))
            1.0
            sage: RDF.coerce(numpy.float64('1'))
            1.0

            sage: RDF.coerce(pi)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Symbolic Ring to Real Double Field

        Test that :trac:`15695` is fixed (see also :trac:`18076`)::

            sage: 1j + numpy.float64(2)
            2.00000000000000 + 1.00000000000000*I
            sage: parent(_)
            Complex Field with 53 bits of precision
        """
        if S is int or S is float:
            return ToRDF(S)

        from .rational_field import QQ
        from .real_lazy import RLF
        if S is ZZ or S is QQ or S is RLF:
            return ToRDF(S)

        from .real_mpfr import RR, RealField_class
        if isinstance(S, RealField_class):
            if S.prec() >= 53:
                return ToRDF(S)
            else:
                return None
        elif is_numpy_type(S):
            import numpy
            if issubclass(S, numpy.integer) or issubclass(S, numpy.floating):
                return ToRDF(S)
            else:
                return None

        connecting = RR._internal_coerce_map_from(S)
        if connecting is not None:
            return ToRDF(RR) * connecting

    def _magma_init_(self, magma):
        r"""
        Return a string representation of ``self`` in the Magma language.

        EXAMPLES:

        Magma handles precision in decimal digits, so we lose a bit::

            sage: magma(RDF) # optional - magma # indirect doctest
            Real field of precision 15
            sage: 10^15 < 2^53 < 10^16
            True

        When we convert back from Magma, we convert to a generic real field
        that has 53 bits of precision::

            sage: magma(RDF).sage() # optional - magma
            Real Field with 53 bits of precision
        """
        return "RealField(%s : Bits := true)" % self.prec()

    def _fricas_init_(self):
        r"""
        Return the FriCAS representation of the real double field.

        EXAMPLES::

            sage: fricas(RDF)    # indirect doctest, optional - fricas
            DoubleFloat
        """
        return "DoubleFloat"

    def _polymake_init_(self):
        r"""
        Return the polymake representation of the real double field.

        EXAMPLES::

            sage: polymake(RDF)    #optional - polymake # indirect doctest
            Float
        """
        return '"Float"'

    def precision(self):
        """
        Return the precision of this real double field in bits.

        Always returns 53.

        EXAMPLES::

            sage: RDF.precision()
            53
        """
        return 53

    prec = precision

    def to_prec(self, prec):
        """
        Return the real field to the specified precision. As doubles have
        fixed precision, this will only return a real double field if ``prec``
        is exactly 53.

        EXAMPLES::

            sage: RDF.to_prec(52)
            Real Field with 52 bits of precision
            sage: RDF.to_prec(53)
            Real Double Field
        """
        if prec == 53:
            return self
        else:
            from .real_mpfr import RealField
            return RealField(prec)


    def gen(self, n=0):
        """
        Return the generator of the real double field.

        EXAMPLES::

            sage: RDF.0
            1.0
            sage: RDF.gens()
            (1.0,)
        """
        if n != 0:
            raise ValueError("only 1 generator")
        return RealDoubleElement(1)

    def ngens(self):
        """
        Return the number of generators which is always 1.

        EXAMPLES::

            sage: RDF.ngens()
            1
        """
        return 1

    def characteristic(self):
        """
        Returns 0, since the field of real numbers has characteristic 0.

        EXAMPLES::

            sage: RDF.characteristic()
            0
        """
        return Integer(0)

    cdef _new_c(self, double value):
        cdef RealDoubleElement x
        x = PY_NEW(RealDoubleElement)
        x._value = value
        return x

    def random_element(self, double min=-1, double max=1):
        """
        Return a random element of this real double field in the interval
        ``[min, max]``.

        EXAMPLES::

            sage: RDF.random_element().parent() is RDF
            True
            sage: -1 <= RDF.random_element() <= 1
            True
            sage: 100 <= RDF.random_element(min=100, max=110) <= 110
            True
        """
        cdef randstate rstate = current_randstate()

        return self._new_c((max-min)*rstate.c_rand_double() + min)

    def name(self):
        """
        The name of ``self``.

        EXAMPLES::

            sage: RDF.name()
            'RealDoubleField'
        """
        return "RealDoubleField"

    def __hash__(self):
        """
        Return the hash value of ``self``.

        This class is intended for use as a singleton so any instance
        of it should be equivalent from a hashing perspective.

        TESTS::

            sage: from sage.rings.real_double import RealDoubleField_class
            sage: hash(RDF) == hash(RealDoubleField_class())
            True
        """
        return 1157042230

    def pi(self):
        r"""
        Returns `\pi` to double-precision.

        EXAMPLES::

            sage: RDF.pi()
            3.141592653589793
            sage: RDF.pi().sqrt()/2
            0.8862269254527579
        """
        return self(M_PI)

    def euler_constant(self):
        """
        Return Euler's gamma constant to double precision.

        EXAMPLES::

            sage: RDF.euler_constant()
            0.5772156649015329
        """
        return self(M_EULER)

    def log2(self):
        r"""
        Return `\log(2)` to the precision of this field.

        EXAMPLES::

            sage: RDF.log2()
            0.6931471805599453
            sage: RDF(2).log()
            0.6931471805599453
        """
        return self(M_LN2)

    def factorial(self, int n):
        """
        Return the factorial of the integer `n` as a real number.

        EXAMPLES::

            sage: RDF.factorial(100)
            9.332621544394415e+157
        """
        if n < 0:
            raise ArithmeticError("n must be nonnegative")
        return self(gsl_sf_fact(n))

    def zeta(self, n=2):
        """
        Return an `n`-th root of unity in the real field, if one
        exists, or raise a ``ValueError`` otherwise.

        EXAMPLES::

            sage: RDF.zeta()
            -1.0
            sage: RDF.zeta(1)
            1.0
            sage: RDF.zeta(5)
            Traceback (most recent call last):
            ...
            ValueError: No 5th root of unity in self
        """
        if n == 1:
            return self(1)
        elif n == 2:
            return self(-1)
        raise ValueError("No %sth root of unity in self" % n)

    def NaN(self):
        """
        Return Not-a-Number ``NaN``.

        EXAMPLES::

            sage: RDF.NaN()
            NaN
        """
        return self(0)/self(0)

    nan = NaN

    def _factor_univariate_polynomial(self, f):
        """
        Factor the univariate polynomial ``f``.

        INPUT:

        - ``f`` -- a univariate polynomial defined over the double precision
          real numbers

        OUTPUT:

        - A factorization of ``f`` over the double precision real numbers
          into a unit and monic irreducible factors

        .. NOTE::

            This is a helper method for
            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.factor`.

        TESTS::

            sage: R.<x> = RDF[]
            sage: RDF._factor_univariate_polynomial(x)
            x
            sage: RDF._factor_univariate_polynomial(2*x)
            (2.0) * x
            sage: RDF._factor_univariate_polynomial(x^2)
            x^2
            sage: RDF._factor_univariate_polynomial(x^2 + 1)
            x^2 + 1.0
            sage: RDF._factor_univariate_polynomial(x^2 - 1)
            (x - 1.0) * (x + 1.0)

        The implementation relies on the ``roots()`` method which often reports
        roots not to be real even though they are::

            sage: f = (x-1)^3
            sage: f.roots(ring=CDF)  # abs tol 2e-5
            [(1.0000065719436413, 1),
             (0.9999967140281792 - 5.691454546815028e-06*I, 1),
             (0.9999967140281792 + 5.691454546815028e-06*I, 1)]

        This leads to the following incorrect factorization::

            sage: f.factor()  # abs tol 2e-5
            (x - 1.0000065719436413) * (x^2 - 1.9999934280563585*x + 0.9999934280995487)
        """
        roots = f.roots(sage.rings.complex_double.CDF)

        # collect real roots and conjugate pairs of non-real roots
        real_roots = [(r, e) for r, e in roots if r.imag().is_zero()]
        non_real_roots = {r: e for r, e in roots if not r.imag().is_zero()}
        assert all(non_real_roots[r.conj()] == e for r, e in non_real_roots.items()), "Bug in root finding code over RDF - roots must always come in conjugate pairs"
        non_real_roots = [(r, e) for r, e in non_real_roots.items() if r.imag() > 0]

        # turn the roots into irreducible factors
        x = f.parent().gen()
        real_factors = [(x - r.real(), e) for r, e in real_roots]
        non_real_factors = [(x**2 - (r + r.conj()).real()*x + (r*r.conj()).real(), e) for r, e in non_real_roots]

        # make the factors monic
        from sage.structure.factorization import Factorization
        return Factorization([(g.monic(), e) for g, e in real_factors + non_real_factors], f.leading_coefficient())


cdef class RealDoubleElement(FieldElement):
    """
    An approximation to a real number using double precision floating
    point numbers. Answers derived from calculations with such
    approximations may differ from what they would be if those
    calculations were performed with true real numbers. This is due to
    the rounding errors inherent to finite precision calculations.
    """

    __array_interface__ = {'typestr': '=f8'}

    def __cinit__(self):
        """
        Initialize ``self`` for cython.

        EXAMPLES::

            sage: RDF(2.3) # indirect doctest
            2.3
        """
        (<Element>self)._parent = _RDF

    def __init__(self, x):
        """
        Create a new ``RealDoubleElement`` with value ``x``.

        EXAMPLES::

            sage: RDF(10^100)
            1e+100

        TESTS::

            sage: from gmpy2 import *
            sage: RDF(mpz(42))
            42.0
            sage: RDF(mpq(3/4))
            0.75
            sage: RDF(mpq('4.1'))
            4.1
        """
        self._value = float(x)

    def _magma_init_(self, magma):
        r"""
        Return a string representation of ``self`` in the Magma language.

        EXAMPLES::

            sage: RDF(10.5)
            10.5
            sage: magma(RDF(10.5)) # optional - magma # indirect doctest
            10.5000000000000
        """
        return "%s!%s" % (self.parent()._magma_init_(magma), self)

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: a = RDF(-2.7)
            sage: loads(dumps(a)) == a
            True
        """
        return RealDoubleElement, (self._value, )

    cdef _new_c(self, double value):
        cdef RealDoubleElement x
        x = PY_NEW(RealDoubleElement)
        x._value = value
        return x

    def prec(self):
        """
        Return the precision of this number in bits.

        Always returns 53.

        EXAMPLES::

            sage: RDF(0).prec()
            53
        """

        return 53

    def ulp(self):
        """
        Returns the unit of least precision of ``self``, which is the
        weight of the least significant bit of ``self``. This is always
        a strictly positive number. It is also the gap between this
        number and the closest number with larger absolute value that
        can be represented.

        EXAMPLES::

            sage: a = RDF(pi)
            sage: a.ulp()
            4.440892098500626e-16
            sage: b = a + a.ulp()

        Adding or subtracting an ulp always gives a different number::

            sage: a + a.ulp() == a
            False
            sage: a - a.ulp() == a
            False
            sage: b + b.ulp() == b
            False
            sage: b - b.ulp() == b
            False

        Since the default rounding mode is round-to-nearest, adding or
        subtracting something less than half an ulp always gives the
        same number, unless the result has a smaller ulp. The latter
        can only happen if the input number is (up to sign) exactly a
        power of 2::

            sage: a - a.ulp()/3 == a
            True
            sage: a + a.ulp()/3 == a
            True
            sage: b - b.ulp()/3 == b
            True
            sage: b + b.ulp()/3 == b
            True
            sage: c = RDF(1)
            sage: c - c.ulp()/3 == c
            False
            sage: c.ulp()
            2.220446049250313e-16
            sage: (c - c.ulp()).ulp()
            1.1102230246251565e-16

        The ulp is always positive::

            sage: RDF(-1).ulp()
            2.220446049250313e-16

        The ulp of zero is the smallest positive number in RDF::

            sage: RDF(0).ulp()
            5e-324
            sage: RDF(0).ulp()/2
            0.0

        Some special values::

            sage: a = RDF(1)/RDF(0); a
            +infinity
            sage: a.ulp()
            +infinity
            sage: (-a).ulp()
            +infinity
            sage: a = RDF('nan')
            sage: a.ulp() is a
            True

        The ulp method works correctly with small numbers::

            sage: u = RDF(0).ulp()
            sage: u.ulp() == u
            True
            sage: x = u * (2^52-1)  # largest denormal number
            sage: x.ulp() == u
            True
            sage: x = u * 2^52  # smallest normal number
            sage: x.ulp() == u
            True

        """
        # First, check special values
        if self._value == 0:
            return RealDoubleElement(libc.math.ldexp(1.0, -1074))
        if gsl_isnan(self._value):
            return self
        if gsl_isinf(self._value):
            return self.abs()

        # Normal case
        cdef int e
        libc.math.frexp(self._value, &e)
        e -= 53
        # Correction for denormals
        if e < -1074:
            e = -1074
        return RealDoubleElement(libc.math.ldexp(1.0, e))

    def real(self):
        """
        Return ``self`` - we are already real.

        EXAMPLES::

            sage: a = RDF(3)
            sage: a.real()
            3.0
        """
        return self

    def imag(self):
        """
        Return the imaginary part of this number, which is zero.

        EXAMPLES::

            sage: a = RDF(3)
            sage: a.imag()
            0.0
        """
        return RealDoubleElement(0)

    def __complex__(self):
        """
        Return ``self`` as a python complex number.

        EXAMPLES::

            sage: a = 2303
            sage: RDF(a)
            2303.0
            sage: complex(RDF(a))
            (2303+0j)
        """
        return complex(self._value,0)

    def _integer_(self, ZZ=None):
        """
        If this floating-point number is actually an integer, return
        that integer.  Otherwise, raise an exception.

        EXAMPLES::

            sage: ZZ(RDF(237.0)) # indirect doctest
            237
            sage: ZZ(RDF(0.0/0.0))
            Traceback (most recent call last):
            ...
            ValueError: cannot convert float NaN to integer
            sage: ZZ(RDF(1.0/0.0))
            Traceback (most recent call last):
            ...
            OverflowError: cannot convert float infinity to integer
            sage: ZZ(RDF(-123456789.0))
            -123456789
            sage: ZZ(RDF((2.0))^290)
            1989292945639146568621528992587283360401824603189390869761855907572637988050133502132224
            sage: ZZ(RDF(-2345.67))
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert non-integral float to integer
        """
        return Integer(self._value)

    def __mpfr__(self):
        """
        Convert Sage ``RealDoubleElement`` to gmpy2 ``mpfr``.

        EXAMPLES::

            sage: RDF(42.2).__mpfr__()
            mpfr('42.200000000000003')
            sage: from gmpy2 import mpfr
            sage: mpfr(RDF(5.1))
            mpfr('5.0999999999999996')

        TESTS::

            sage: RDF().__mpfr__(); raise NotImplementedError("gmpy2 is not installed")
            Traceback (most recent call last):
            ...
            NotImplementedError: gmpy2 is not installed
        """
        return gmpy2.mpfr(self._value)

    def _interface_init_(self, I=None):
        """
        Return ``self`` formatted as a string, suitable as input to another
        computer algebra system. (This is the default function used for
        exporting to other computer algebra systems.)

        EXAMPLES::

            sage: s1 = RDF(sin(1)); s1
            0.8414709848078965
            sage: s1._interface_init_()
            '0.8414709848078965'
            sage: s1 == RDF(gp(s1))
            True
        """
        return repr(self._value)

    def _mathematica_init_(self):
        """
        TESTS:

        Check that :trac:`28814` is fixed::

            sage: mathematica(RDF(1e25))   # optional - mathematica
            1.*^25
            sage: mathematica(RDF(1e-25))  # optional - mathematica
            1.*^-25
        """
        from .real_mpfr import RR
        return RR(self._value)._mathematica_init_()

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(RDF(NaN))
            RDF(NaN)
            sage: sage_input(RDF(-infinity), verify=True)
            # Verified
            -RDF(infinity)
            sage: sage_input(RDF(-infinity)*polygen(RDF))
            R.<x> = RDF[]
            -RDF(infinity)*x + RDF(NaN)
            sage: sage_input(RDF(pi), verify=True)
            # Verified
            RDF(3.1415926535897931)
            sage: sage_input(RDF(-e), verify=True, preparse=False)
            # Verified
            -RDF(2.718281828459045...)
            sage: sage_input(RDF(pi)*polygen(RDF), verify=True, preparse=None)
            # Verified
            R = RDF['x']
            x = R.gen()
            3.1415926535897931*x
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: RDF(22/7)._sage_input_(sib, True)
            {atomic:3.1428571428571428}
            sage: RDF(22/7)._sage_input_(sib, False)
            {call: {atomic:RDF}({atomic:3.1428571428571428})}
        """
        cdef int isinf = gsl_isinf(self._value)
        cdef bint isnan = gsl_isnan(self._value)
        if isinf or isnan:
            if isnan:
                v = sib.name('NaN')
            else:
                v = sib.name('infinity')
            v = sib(self.parent())(v)
            if isinf < 0:
                v = -v
            return v

        from sage.rings.all import ZZ, RR

        cdef bint negative = self._value < 0
        if negative:
            self = -self

        # There are five possibilities for printing this floating-point
        # number, ordered from prettiest to ugliest (IMHO).
        # 1) An integer: 42
        # 2) A simple literal: 3.14159
        # 3) A coerced integer: RDF(42)
        # 4) A coerced literal: RDF(3.14159)
        # 5) A coerced RR value: RDF(RR('3.14159'))

        # str(self) works via libc, which we don't necessarily trust
        # to produce the best possible answer.  So this function prints
        # via RR/MPFR.  Without the preparser, input works via libc as
        # well, but we don't have a choice about that.

        # To use choice 1 or choice 3, this number must be an integer.
        cdef bint can_use_int_literal = \
            self.abs() < (Integer(1) << self.prec()) and self in ZZ

        self_str = RR(self._value).str(truncate=False, skip_zeroes=True)

        # To use choice 2 or choice 4, we must be able to read
        # numbers of this precision as a literal.
        cdef bint can_use_float_literal = \
            (sib.preparse() or float(self_str) == self)

        if can_use_int_literal or can_use_float_literal:
            if can_use_int_literal:
                v = sib.int(self._integer_())
            else:
                v = sib.float_str(self_str)
        else:
            v = sib(RR(self))
        if not coerced:
            v = sib(self.parent())(v)

        if negative:
            v = -v

        return v

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: RDF(-2/3)
            -0.6666666666666666
            sage: a = RDF(2); a
            2.0
            sage: a^2
            4.0
            sage: RDF("nan")
            NaN
            sage: RDF(1)/RDF(0)
            +infinity
            sage: RDF(-1)/RDF(0)
            -infinity
        """
        return double_repr(self._value)

    def __format__(self, format_spec):
        """
        Return a formatted string representation of this real number.

        EXAMPLES::

            sage: format(RDF(32/3), '.4f')
            '10.6667'
            sage: '{:.4e}'.format(RDF(2/3))
            '6.6667e-01'
        """
        return format(float(self), format_spec)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(RDF(3.4)) # indirect doctest
            3.4
            sage: latex(RDF(2e-100)) # indirect doctest
            2 \times 10^{-100}
        """
        s = self.str()
        parts = s.split('e')
        if len(parts) > 1:
            # scientific notation
            if parts[1][0] == '+':
                parts[1] = parts[1][1:]
            s = "%s \\times 10^{%s}" % (parts[0], parts[1])
        return s

    def __hash__(self):
        """
        Return the hash of ``self``, which coincides with the python float
        (and often int) type.

        EXAMPLES::

            sage: hash(RDF(1.2)) == hash(1.2r)
            True
        """
        return hash(self._value)

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of ``self`` under the homomorphism from the rational
        field to ``codomain``.

        This always just returns ``self`` coerced into the ``codomain``.

        EXAMPLES::

            sage: RDF(2.1)._im_gens_(RR, [RR(1)])
            2.10000000000000
            sage: R = RealField(20)
            sage: RDF(2.1)._im_gens_(R, [R(1)])
            2.1000
        """
        return codomain(self) # since 1 |--> 1

    def str(self):
        """
        Return the informal string representation of ``self``.

        EXAMPLES::

            sage: a = RDF('4.5'); a.str()
            '4.5'
            sage: a = RDF('49203480923840.2923904823048'); a.str()
            '49203480923840.29'
            sage: a = RDF(1)/RDF(0); a.str()
            '+infinity'
            sage: a = -RDF(1)/RDF(0); a.str()
            '-infinity'
            sage: a = RDF(0)/RDF(0); a.str()
            'NaN'

        We verify consistency with ``RR`` (mpfr reals)::

            sage: str(RR(RDF(1)/RDF(0))) == str(RDF(1)/RDF(0))
            True
            sage: str(RR(-RDF(1)/RDF(0))) == str(-RDF(1)/RDF(0))
            True
            sage: str(RR(RDF(0)/RDF(0))) == str(RDF(0)/RDF(0))
            True
        """
        return double_repr(self._value)

    def __copy__(self):
        """
        Return copy of ``self``, which since ``self`` is immutable, is just
        ``self``.

        EXAMPLES::

            sage: r = RDF('-1.6')
            sage: r.__copy__() is r
            True
        """
        return self

    def __deepcopy__(self, memo):
        """
        EXAMPLES::

            sage: r = RDF('-1.6')
            sage: deepcopy(r) is r
            True
        """
        return self

    def integer_part(self):
        """
        If in decimal this number is written ``n.defg``, returns ``n``.

        EXAMPLES::

            sage: r = RDF('-1.6')
            sage: a = r.integer_part(); a
            -1
            sage: type(a)
            <class 'sage.rings.integer.Integer'>
            sage: r = RDF(0.0/0.0)
            sage: a = r.integer_part()
            Traceback (most recent call last):
            ...
            TypeError: Attempt to get integer part of NaN
        """
        if gsl_isnan(self._value):
            raise TypeError("Attempt to get integer part of NaN")
        else:
            return Integer(int(self._value))

    def sign_mantissa_exponent(self):
        r"""
        Return the sign, mantissa, and exponent of ``self``.

        In Sage (as in MPFR), floating-point numbers of precision `p`
        are of the form `s m 2^{e-p}`, where `s \in \{-1, 1\}`,
        `2^{p-1} \leq m < 2^p`, and `-2^{30} + 1 \leq e \leq 2^{30} -
        1`; plus the special values ``+0``, ``-0``, ``+infinity``,
        ``-infinity``, and ``NaN`` (which stands for Not-a-Number).

        This function returns `s`, `m`, and `e-p`.  For the special values:

        - ``+0`` returns ``(1, 0, 0)``
        - ``-0`` returns ``(-1, 0, 0)``
        - the return values for ``+infinity``, ``-infinity``, and ``NaN`` are
          not specified.

        EXAMPLES::

            sage: a = RDF(exp(1.0)); a
            2.718281828459045
            sage: sign,mantissa,exponent = RDF(exp(1.0)).sign_mantissa_exponent()
            sage: sign,mantissa,exponent
            (1, 6121026514868073, -51)
            sage: sign*mantissa*(2**exponent) == a
            True

        The mantissa is always a nonnegative number::

            sage: RDF(-1).sign_mantissa_exponent()
            (-1, 4503599627370496, -52)

        TESTS::

            sage: RDF('+0').sign_mantissa_exponent()
            (1, 0, 0)
            sage: RDF('-0').sign_mantissa_exponent()
            (-1, 0, 0)
        """
        from sage.rings.all import RR
        return RR(self._value).sign_mantissa_exponent()

    def as_integer_ratio(self):
        """
        Return a coprime pair of integers ``(a, b)`` such that ``self``
        equals ``a / b`` exactly.

        EXAMPLES::

            sage: RDF(0).as_integer_ratio()
            (0, 1)
            sage: RDF(1/3).as_integer_ratio()
            (6004799503160661, 18014398509481984)
            sage: RDF(37/16).as_integer_ratio()
            (37, 16)
            sage: RDF(3^60).as_integer_ratio()
            (42391158275216203520420085760, 1)
        """
        nd = float.as_integer_ratio(self._value)
        return (Integer(nd[0]), Integer(nd[1]))

    ########################
    #   Basic Arithmetic
    ########################
    def __invert__(self):
        """
        Compute the multiplicative inverse of ``self``.

        EXAMPLES::

            sage: a = RDF(-1.5)*RDF(2.5)
            sage: a.__invert__()
            -0.26666666666666666
            sage: ~a
            -0.26666666666666666
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = 1.0 / self._value
        return x

    cpdef _add_(self, right):
        """
        Add two real numbers with the same parent.

        EXAMPLES::

            sage: RDF('-1.5') + RDF('2.5') # indirect doctest
            1.0
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value + (<RealDoubleElement>right)._value
        return x

    cpdef _sub_(self, right):
        """
        Subtract two real numbers with the same parent.

        EXAMPLES::

            sage: RDF('-1.5') - RDF('2.5') # indirect doctest
            -4.0
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value - (<RealDoubleElement>right)._value
        return x

    cpdef _mul_(self, right):
        """
        Multiply two real numbers with the same parent.

        EXAMPLES::

            sage: RDF('-1.5') * RDF('2.5') # indirect doctest
            -3.75
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value * (<RealDoubleElement>right)._value
        return x

    cpdef _div_(self, right):
        """
        Divide ``self`` by ``right``.

        EXAMPLES::

            sage: RDF('-1.5') / RDF('2.5') # indirect doctest
            -0.6
            sage: RDF(1)/RDF(0)
            +infinity
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value / (<RealDoubleElement>right)._value
        return x

    def __neg__(self):
        """
        Negates ``self``.

        EXAMPLES::

            sage: -RDF('-1.5')
            1.5
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = -self._value
        return x

    def conjugate(self):
        r"""
        Returns the complex conjugate of this real number, which is
        the real number itself.

        EXAMPLES::

            sage: RDF(4).conjugate()
            4.0
        """
        return self

    def __abs__(self):
        """
        Returns the absolute value of ``self``.

        EXAMPLES::

            sage: abs(RDF(1.5))
            1.5
            sage: abs(RDF(-1.5))
            1.5
            sage: abs(RDF(0.0))
            0.0
            sage: abs(RDF(-0.0))
            0.0
        """
        # Use signbit instead of >= to handle -0.0 correctly
        if not libc.math.signbit(self._value):
            return self
        else:
            return self._new_c(-self._value)

    cpdef RealDoubleElement abs(RealDoubleElement self):
        """
        Returns the absolute value of ``self``.

        EXAMPLES::

            sage: RDF(1e10).abs()
            10000000000.0
            sage: RDF(-1e10).abs()
            10000000000.0
        """
        if self._value >= 0:
            return self
        else:
            return self._new_c(-self._value)

    def __lshift__(x, y):
        """
        LShifting a double is not supported; nor is lshifting a
        :class:`RealDoubleElement`.

        TESTS::

            sage: RDF(2) << 3
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for <<
        """
        raise TypeError("unsupported operand type(s) for <<")

    def __rshift__(x, y):
        """
        RShifting a double is not supported; nor is rshifting a
        :class:`RealDoubleElement`.

        TESTS::

            sage: RDF(2) >> 3
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for >>
        """
        raise TypeError("unsupported operand type(s) for >>")

    def multiplicative_order(self):
        r"""
        Returns `n` such that ``self^n == 1``.

        Only `\pm 1` have finite multiplicative order.

        EXAMPLES::

            sage: RDF(1).multiplicative_order()
            1
            sage: RDF(-1).multiplicative_order()
            2
            sage: RDF(3).multiplicative_order()
            +Infinity
        """
        if self._value == 1:
            return 1
        elif self._value == -1:
            return 2
        return sage.rings.infinity.infinity

    def sign(self):
        """
        Returns -1,0, or 1 if ``self`` is negative, zero, or positive;
        respectively.

        EXAMPLES::

            sage: RDF(-1.5).sign()
            -1
            sage: RDF(0).sign()
            0
            sage: RDF(2.5).sign()
            1
        """
        if not self._value:
            return 0
        if self._value > 0:
            return 1
        return -1


    ###################
    # Rounding etc
    ###################

    def round(self):
        """
        Given real number `x`, rounds up if fractional part is greater than
        `0.5`, rounds down if fractional part is less than `0.5`.

        EXAMPLES::

            sage: RDF(0.49).round()
            0
            sage: a=RDF(0.51).round(); a
            1
        """
        return Integer(round(self._value))

    def floor(self):
        """
        Return the floor of ``self``.

        EXAMPLES::

            sage: RDF(2.99).floor()
            2
            sage: RDF(2.00).floor()
            2
            sage: RDF(-5/2).floor()
            -3
        """
        return Integer(math.floor(self._value))

    def ceil(self):
        """
        Return the ceiling of ``self``.

        EXAMPLES::

            sage: RDF(2.99).ceil()
            3
            sage: RDF(2.00).ceil()
            2
            sage: RDF(-5/2).ceil()
            -2
        """
        return Integer(math.ceil(self._value))

    ceiling = ceil

    def trunc(self):
        """
        Truncates this number (returns integer part).

        EXAMPLES::

            sage: RDF(2.99).trunc()
            2
            sage: RDF(-2.00).trunc()
            -2
            sage: RDF(0.00).trunc()
            0
        """
        return Integer(int(self._value))

    def frac(self):
        """
        Return a real number in `(-1, 1)`. It satisfies the relation:
        ``x = x.trunc() + x.frac()``

        EXAMPLES::

            sage: RDF(2.99).frac()
            0.9900000000000002
            sage: RDF(2.50).frac()
            0.5
            sage: RDF(-2.79).frac()
            -0.79
        """
        return self._new_c(self._value - int(self._value))

    ###########################################
    # Conversions
    ###########################################

    def __float__(self):
        """
        Return ``self`` as a python float.

        EXAMPLES::

            sage: float(RDF(1.5))
            1.5
            sage: type(float(RDF(1.5)))
            <... 'float'>
        """
        return self._value

    def _rpy_(self):
        """
        Returns ``self.__float__()`` for rpy to convert into the
        appropriate R object.

        EXAMPLES::

            sage: n = RDF(2.0)
            sage: n._rpy_()
            2.0
            sage: type(n._rpy_())
            <... 'float'>
        """
        return self.__float__()

    def __int__(self):
        """
        Return integer truncation of this real number.

        EXAMPLES::

            sage: int(RDF(2.99))
            2
            sage: int(RDF(-2.99))
            -2
        """
        return int(self._value)

    def _complex_mpfr_field_(self, CC):
        """
        EXAMPLES::

            sage: a = RDF(1/3)
            sage: CC(a)
            0.333333333333333
            sage: a._complex_mpfr_field_(CC)
            0.333333333333333

        If we coerce to a higher-precision field the extra bits appear
        random; they are actually 0's in base 2.

        ::

            sage: a._complex_mpfr_field_(ComplexField(100))
            0.33333333333333331482961625625
            sage: a._complex_mpfr_field_(ComplexField(100)).str(2)
            '0.01010101010101010101010101010101010101010101010101010100000000000000000000000000000000000000000000000'
        """
        return CC(self._value)

    def _complex_double_(self, CDF):
        """
        Return ``self`` as a complex double.

        EXAMPLES::

            sage: CDF(RDF(1/3)) # indirect doctest
            0.3333333333333333
        """
        return CDF(self._value)

    def __pari__(self):
        """
        Return a PARI representation of ``self``.

        EXAMPLES::

            sage: RDF(1.5).__pari__()
            1.50000000000000
        """
        global new_gen_from_real_double_element
        if new_gen_from_real_double_element is None:
            from sage.libs.pari.convert_sage_real_double import new_gen_from_real_double_element
        return new_gen_from_real_double_element(self)


    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    def is_NaN(self):
        """
        Check if ``self`` is ``NaN``.

        EXAMPLES::

            sage: RDF(1).is_NaN()
            False
            sage: a = RDF(0)/RDF(0)
            sage: a.is_NaN()
            True
        """
        return gsl_isnan(self._value)

    def is_positive_infinity(self):
        r"""
        Check if ``self`` is `+\infty`.

        EXAMPLES::

            sage: a = RDF(1)/RDF(0)
            sage: a.is_positive_infinity()
            True
            sage: a = RDF(-1)/RDF(0)
            sage: a.is_positive_infinity()
            False
        """
        return gsl_isinf(self._value) > 0

    def is_negative_infinity(self):
        r"""
        Check if ``self`` is `-\infty`.

        EXAMPLES::

            sage: a = RDF(2)/RDF(0)
            sage: a.is_negative_infinity()
            False
            sage: a = RDF(-3)/RDF(0)
            sage: a.is_negative_infinity()
            True
        """
        return gsl_isinf(self._value) < 0

    def is_infinity(self):
        r"""
        Check if ``self`` is `\infty`.

        EXAMPLES::

            sage: a = RDF(2); b = RDF(0)
            sage: (a/b).is_infinity()
            True
            sage: (b/a).is_infinity()
            False
        """
        return gsl_isinf(self._value)

    cpdef _richcmp_(left, right, int op):
        """
        Rich comparison of ``left`` and ``right``.

        EXAMPLES::

            sage: RDF(2) < RDF(0)
            False
            sage: RDF(2) == RDF(4/2)
            True
            sage: RDF(-2) > RDF(-4)
            True

        TESTS:

        Check comparisons with ``NaN`` (:trac:`16515`)::

            sage: n = RDF('NaN')
            sage: n == n
            False
            sage: n == RDF(1)
            False
        """
        # We really need to use the correct operators, to deal
        # correctly with NaNs.
        cdef double x = (<RealDoubleElement>left)._value
        cdef double y = (<RealDoubleElement>right)._value
        if op == Py_LT:
            return x < y
        elif op == Py_LE:
            return x <= y
        elif op == Py_EQ:
            return x == y
        elif op == Py_NE:
            return x != y
        elif op == Py_GT:
            return x > y
        else:
            return x >= y


    ############################
    # Special Functions
    ############################

    def NaN(self):
        """
        Return Not-a-Number ``NaN``.

        EXAMPLES::

            sage: RDF.NaN()
            NaN
        """
        return self(0)/self(0)

    nan = NaN

    def sqrt(self, extend=True, all=False):
        """
        The square root function.

        INPUT:

        -  ``extend`` -- bool (default: ``True``); if ``True``, return a
           square root in a complex field if necessary if ``self`` is negative;
           otherwise raise a ``ValueError``.

        -  ``all`` -- bool (default: ``False``); if ``True``, return a
           list of all square roots.

        EXAMPLES::

            sage: r = RDF(4.0)
            sage: r.sqrt()
            2.0
            sage: r.sqrt()^2 == r
            True

        ::

            sage: r = RDF(4344)
            sage: r.sqrt()
            65.90902821313632
            sage: r.sqrt()^2 - r
            0.0

        ::

            sage: r = RDF(-2.0)
            sage: r.sqrt()
            1.4142135623730951*I

        ::

            sage: RDF(2).sqrt(all=True)
            [1.4142135623730951, -1.4142135623730951]
            sage: RDF(0).sqrt(all=True)
            [0.0]
            sage: RDF(-2).sqrt(all=True)
            [1.4142135623730951*I, -1.4142135623730951*I]
        """
        if self._value >= 0:
            x = self._new_c(libc.math.sqrt(self._value))
            if all:
                if x.is_zero():
                    return [x]
                else:
                    return [x, -x]
            else:
                return x
        if not extend:
            raise ValueError("negative number %s does not have a square root in the real field" % self)
        import sage.rings.complex_double
        return self._complex_double_(sage.rings.complex_double.CDF).sqrt(all=all)

    def is_square(self):
        """
        Return whether or not this number is a square in this field. For
        the real numbers, this is ``True`` if and only if ``self`` is
        non-negative.

        EXAMPLES::

            sage: RDF(3.5).is_square()
            True
            sage: RDF(0).is_square()
            True
            sage: RDF(-4).is_square()
            False
        """
        return self._value >= 0

    def is_integer(self):
        """
        Return True if this number is a integer

        EXAMPLES::

            sage: RDF(3.5).is_integer()
            False
            sage: RDF(3).is_integer()
            True
        """
        return self._value in ZZ


    def cube_root(self):
        """
        Return the cubic root (defined over the real numbers) of ``self``.

        EXAMPLES::

            sage: r = RDF(125.0); r.cube_root()
            5.000000000000001
            sage: r = RDF(-119.0)
            sage: r.cube_root()^3 - r  # rel tol 1
            -1.4210854715202004e-14
        """
        return self.nth_root(3)


    def nth_root(self, int n):
        """
        Return the `n^{th}` root of ``self``.

        INPUT:

        -  ``n`` -- an integer

        OUTPUT:

        The output is a complex double if ``self`` is negative and `n` is even,
        otherwise it is a real double.

        EXAMPLES::

            sage: r = RDF(-125.0); r.nth_root(3)
            -5.000000000000001
            sage: r.nth_root(5)
            -2.6265278044037674
            sage: RDF(-2).nth_root(5)^5  # rel tol 1e-15
            -2.000000000000001
            sage: RDF(-1).nth_root(5)^5
            -1.0
            sage: RDF(3).nth_root(10)^10
            2.9999999999999982
            sage: RDF(-1).nth_root(2)
            6.123233995736757e-17 + 1.0*I
            sage: RDF(-1).nth_root(4)
            0.7071067811865476 + 0.7071067811865475*I
        """
        if n == 0:
            return RealDoubleElement(float('nan'))
        if self._value < 0:
            if GSL_IS_EVEN(n):
                return self._complex_double_(sage.rings.complex_double.CDF).nth_root(n)
            else:
                return - ( (-self) ** (float(1)/n) )
        else:
            return self ** (float(1)/n)

    cdef __pow_double(self, double exponent, double sign):
        """
        If ``sign == 1`` or ``self >= 0``, return ``self ^ exponent``.
        If ``sign == -1`` and ``self < 0``, return ``- abs(self) ^ exponent``.
        """
        cdef double v = self._value
        if v >= 0:
            if v == 1:
                return self
            elif exponent == 0:
                return self._new_c(1.0)
            elif v == 0:
                if exponent < 0:
                    raise ZeroDivisionError("0.0 cannot be raised to a negative power")
                return self
            sign = 1.0
        else:  # v < 0
            expmod2 = libc.math.fmod(exponent, 2.0)
            if expmod2 == 0.0:
                pass
            elif expmod2 == 1.0:
                sign = -1.0
            else:
                raise ValueError("negative number cannot be raised to a fractional power")
            v = -v
        return self._new_c(sign * gsl_sf_exp(gsl_sf_log(v) * exponent))

    cpdef _pow_(self, other):
        """
        Return ``self`` raised to the real double power ``other``.

        EXAMPLES::

            sage: a = RDF('1.23456')
            sage: a^a
            1.2971114817819216

        TESTS::

            sage: RDF(0) ^ RDF(0.5)
            0.0
            sage: RDF(0) ^ (1/2)
            0.0
            sage: RDF(0) ^ RDF(0)
            1.0
            sage: RDF(0) ^ RDF(-1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: 0.0 cannot be raised to a negative power
            sage: RDF(-1) ^ RDF(0)
            1.0
            sage: RDF(-1) ^ RDF(1)
            -1.0
            sage: RDF(-1) ^ RDF(0.5)
            Traceback (most recent call last):
            ...
            ValueError: negative number cannot be raised to a fractional power
        """
        return self.__pow_double((<RealDoubleElement>other)._value, 1)

    cpdef _pow_int(self, n):
        """
        Return ``self`` raised to the integer power ``n``.

        TESTS::

            sage: RDF(1) ^ (2^1000)
            1.0
            sage: RDF(1) ^ (2^1000 + 1)
            1.0
            sage: RDF(1) ^ (-2^1000)
            1.0
            sage: RDF(1) ^ (-2^1000 + 1)
            1.0
            sage: RDF(-1) ^ (2^1000)
            1.0
            sage: RDF(-1) ^ (2^1000 + 1)
            -1.0
            sage: RDF(-1) ^ (-2^1000)
            1.0
            sage: RDF(-1) ^ (-2^1000 + 1)
            -1.0

        ::

            sage: base = RDF(1.0000000000000002)
            sage: base._pow_int(0)
            1.0
            sage: base._pow_int(1)
            1.0000000000000002
            sage: base._pow_int(2)
            1.0000000000000004
            sage: base._pow_int(3)
            1.0000000000000007
            sage: base._pow_int(2^57)
            78962960182680.42
            sage: base._pow_int(2^57 + 1)
            78962960182680.42

        ::

            sage: base = RDF(-1.0000000000000002)
            sage: base._pow_int(0)
            1.0
            sage: base._pow_int(1)
            -1.0000000000000002
            sage: base._pow_int(2)
            1.0000000000000004
            sage: base._pow_int(3)
            -1.0000000000000007
            sage: base._pow_int(2^57)
            78962960182680.42
            sage: base._pow_int(2^57 + 1)
            -78962960182680.42
        """
        return self.__pow_double(n, -1.0 if (n & 1) else 1.0)

    cdef _pow_long(self, long n):
        """
        Compute ``self`` raised to the power ``n``.

        EXAMPLES::

            sage: RDF('1.23456') ^ 20
            67.64629770385...
            sage: RDF(3) ^ 32
            1853020188851841.0
            sage: RDF(2)^(-1024)
            5.562684646268003e-309

        TESTS::

            sage: base = RDF(1.0000000000000002)
            sage: base ^ RDF(2^31)
            1.000000476837272
            sage: base ^ (2^57)
            78962960182680.42
            sage: base ^ RDF(2^57)
            78962960182680.42
        """
        if -2048 <= n <= 2048:
            # For small exponents, it is possible that the powering
            # is exact either because the base is a power of 2
            # (e.g. 2.0^1000) or because the exact result has few
            # significant digits (e.g. 3.0^10). Here, we use the
            # square-and-multiply algorithm by GSL.
            return self._new_c(gsl_pow_int(self._value, <int>n))
        # If the exponent is sufficiently large in absolute value, the
        # result cannot be exact (except if the base is -1.0, 0.0 or
        # 1.0 but those cases are handled by __pow_double too). The
        # log-and-exp algorithm from __pow_double will be more precise
        # than square-and-multiply.

        # We do need to take care of the sign since the conversion
        # of n to double might change an odd number to an even number.
        return self.__pow_double(<double>n, -1.0 if (n & 1) else 1.0)

    cdef _log_base(self, double log_of_base):
        if self._value == 0:
            return RDF(-1)/RDF(0)
        elif self._value < 0:
            return RDF.NaN()
        sig_on()
        a = self._new_c(gsl_sf_log(self._value) / log_of_base)
        sig_off()
        return a

    def log(self, base=None):
        """
        Return the logarithm.

        INPUT:

        - ``base`` -- integer or ``None`` (default). The base of the
          logarithm. If ``None`` is specified, the base is `e` (the so-called
          natural logarithm).

        OUTPUT:

        The logarithm of ``self``.  If ``self`` is positive, a double
        floating point number. Infinity if ``self`` is zero. A
        imaginary complex floating point number if ``self`` is
        negative.

        EXAMPLES::

            sage: RDF(2).log()
            0.6931471805599453
            sage: RDF(2).log(2)
            1.0
            sage: RDF(2).log(pi)
            0.6055115613982801
            sage: RDF(2).log(10)
            0.30102999566398114
            sage: RDF(2).log(1.5)
            1.7095112913514547
            sage: RDF(0).log()
            -infinity
            sage: RDF(-1).log()
            3.141592653589793*I
            sage: RDF(-1).log(2)  # rel tol 1e-15
            4.532360141827194*I

        TESTS:

        Make sure that we can take the log of small numbers accurately
        and the fix doesn't break preexisting values (:trac:`12557`)::

            sage: R = RealField(128)
            sage: def check_error(x):
            ....:   x = RDF(x)
            ....:   log_RDF = x.log()
            ....:   log_RR = R(x).log()
            ....:   diff = R(log_RDF) - log_RR
            ....:   if abs(diff) < log_RDF.ulp():
            ....:       return True
            ....:   print("logarithm check failed for %s (diff = %s ulp)"% \
            ....:       (x, diff/log_RDF.ulp()))
            ....:   return False
            sage: all( check_error(2^x) for x in range(-100,100) )
            True
            sage: all( check_error(x) for x in sxrange(0.01, 2.00, 0.01) )
            True
            sage: all( check_error(x) for x in sxrange(0.99, 1.01, 0.001) )
            True
            sage: RDF(1.000000001).log()
            1.000000082240371e-09
            sage: RDF(1e-17).log()
            -39.14394658089878
            sage: RDF(1e-50).log()
            -115.12925464970229
        """
        if self < 0:
            from sage.rings.complex_double import CDF
            return CDF(self).log(base)
        if base is None:
            return self._log_base(1)
        else:
            if isinstance(base, RealDoubleElement):
                return self._log_base(base._log_base(1))
            else:
                return self._log_base(gsl_sf_log(float(base)))

    def log2(self):
        """
        Return log to the base 2 of ``self``.

        EXAMPLES::

            sage: r = RDF(16.0)
            sage: r.log2()
            4.0

        ::

            sage: r = RDF(31.9); r.log2()
            4.9954845188775066
        """
        if self < 0:
            from sage.rings.complex_double import CDF
            return CDF(self).log(2)
        sig_on()
        a = self._new_c(gsl_sf_log(self._value) * M_1_LN2)
        sig_off()
        return a


    def log10(self):
        """
        Return log to the base 10 of ``self``.

        EXAMPLES::

            sage: r = RDF('16.0'); r.log10()
            1.2041199826559248
            sage: r.log() / RDF(log(10))
            1.2041199826559246
            sage: r = RDF('39.9'); r.log10()
            1.6009728956867482
        """
        if self < 0:
            from sage.rings.complex_double import CDF
            return CDF(self).log(10)
        sig_on()
        a = self._new_c(gsl_sf_log(self._value) * M_1_LN10)
        sig_off()
        return a

    def logpi(self):
        r"""
        Return log to the base `\pi` of ``self``.

        EXAMPLES::

            sage: r = RDF(16); r.logpi()
            2.4220462455931204
            sage: r.log() / RDF(log(pi))
            2.4220462455931204
            sage: r = RDF('39.9'); r.logpi()
            3.2203023346075152
        """
        if self < 0:
            from sage.rings.complex_double import CDF
            return CDF(self).log(M_PI)
        sig_on()
        a = self._new_c(gsl_sf_log(self._value) * M_1_LNPI)
        sig_off()
        return a

    def exp(self):
        r"""
        Return `e^\mathtt{self}`.

        EXAMPLES::

            sage: r = RDF(0.0)
            sage: r.exp()
            1.0

        ::

            sage: r = RDF('32.3')
            sage: a = r.exp(); a
            106588847274864.47
            sage: a.log()
            32.3

        ::

            sage: r = RDF('-32.3')
            sage: r.exp()
            9.381844588498685e-15

        ::

            sage: RDF(1000).exp()
            +infinity
        """
        sig_on()
        a = self._new_c(gsl_sf_exp(self._value))
        sig_off()
        return a

    def exp2(self):
        """
        Return `2^\mathtt{self}`.

        EXAMPLES::

            sage: r = RDF(0.0)
            sage: r.exp2()
            1.0

        ::

            sage: r = RDF(32.0)
            sage: r.exp2()
            4294967295.9999967

        ::

            sage: r = RDF(-32.3)
            sage: r.exp2()
            1.8911724825302065e-10
        """
        sig_on()
        a = self._new_c(gsl_sf_exp(self._value * M_LN2))
        sig_off()
        return a

    def exp10(self):
        r"""
        Return `10^\mathtt{self}`.

        EXAMPLES::

            sage: r = RDF(0.0)
            sage: r.exp10()
            1.0

        ::

            sage: r = RDF(32.0)
            sage: r.exp10()
            1.0000000000000069e+32

        ::

            sage: r = RDF(-32.3)
            sage: r.exp10()
            5.011872336272702e-33
        """
        sig_on()
        a = self._new_c(gsl_sf_exp(self._value * M_LN10))
        sig_off()
        return a

    def cos(self):
        """
        Return the cosine of ``self``.

        EXAMPLES::

            sage: t=RDF.pi()/2
            sage: t.cos()
            6.123233995736757e-17
        """
        return self._new_c(gsl_sf_cos(self._value))

    def sin(self):
        """
        Return the sine of ``self``.

        EXAMPLES::

            sage: RDF(2).sin()
            0.9092974268256817
        """
        return self._new_c(gsl_sf_sin(self._value))

    def dilog(self):
        r"""
        Return the dilogarithm of ``self``.

        This is defined by the
        series `\sum_n x^n/n^2` for `|x| \le 1`. When the absolute
        value of ``self`` is greater than 1, the returned value is the
        real part of (the analytic continuation to `\CC` of) the
        dilogarithm of ``self``.

        EXAMPLES::

            sage: RDF(1).dilog()  # rel tol 1.0e-13
            1.6449340668482264
            sage: RDF(2).dilog()  # rel tol 1.0e-13
            2.46740110027234
        """
        return self._new_c(gsl_sf_dilog(self._value))

    def restrict_angle(self):
        r"""
        Return a number congruent to ``self`` mod `2\pi` that lies in
        the interval `(-\pi, \pi]`.

        Specifically, it is the unique `x \in (-\pi, \pi]` such
        that ```self`` `= x + 2\pi n` for some `n \in \ZZ`.

        EXAMPLES::

            sage: RDF(pi).restrict_angle()
            3.141592653589793
            sage: RDF(pi + 1e-10).restrict_angle()
            -3.1415926534897936
            sage: RDF(1+10^10*pi).restrict_angle()
            0.9999977606...
        """
        return self._new_c(gsl_sf_angle_restrict_symm(self._value))

    def tan(self):
        """
        Return the tangent of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/3
            sage: q.tan()
            1.7320508075688767
            sage: q = RDF.pi()/6
            sage: q.tan()
            0.5773502691896256
        """
        cdef double denom
        cos = gsl_sf_cos(self._value)
        a = self._new_c(gsl_sf_sin(self._value) / cos)
        return a

    def sincos(self):
        """
        Return a pair consisting of the sine and cosine of ``self``.

        EXAMPLES::

            sage: t = RDF.pi()/6
            sage: t.sincos()
            (0.49999999999999994, 0.8660254037844387)
        """
        return self.sin(), self.cos()

    def hypot(self, other):
        r"""
        Computes the value `\sqrt{s^2 + o^2}` where `s` is ``self`` and `o`
        is ``other`` in such a way as to avoid overflow.

        EXAMPLES::

            sage: x = RDF(4e300); y = RDF(3e300)
            sage: x.hypot(y)
            5e+300
            sage: sqrt(x^2+y^2) # overflow
            +infinity
        """
        sig_on()
        a = self._new_c(gsl_sf_hypot(self._value, float(other)))
        sig_off()
        return a

    def arccos(self):
        """
        Return the inverse cosine of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/3
            sage: i = q.cos()
            sage: i.arccos() == q
            True
        """
        return self._new_c(libc.math.acos(self._value))

    def arcsin(self):
        """
        Return the inverse sine of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/5
            sage: i = q.sin()
            sage: i.arcsin() == q
            True
        """
        return self._new_c(libc.math.asin(self._value))

    def arctan(self):
        """
        Return the inverse tangent of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/5
            sage: i = q.tan()
            sage: i.arctan() == q
            True
        """
        return self._new_c(libc.math.atan(self._value))


    def cosh(self):
        """
        Return the hyperbolic cosine of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/12
            sage: q.cosh()
            1.0344656400955106
        """
        return self._new_c(gsl_ldexp( gsl_sf_exp(self._value) + gsl_sf_exp(-self._value), -1)) # (e^x + e^-x)/2

    def sinh(self):
        """
        Return the hyperbolic sine of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/12
            sage: q.sinh()
            0.26480022760227073
        """
        return self._new_c(gsl_ldexp( gsl_sf_expm1(self._value) - gsl_sf_expm1(-self._value), -1)) # (e^x - e^-x)/2

    def tanh(self):
        """
        Return the hyperbolic tangent of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/12
            sage: q.tanh()
            0.25597778924568454
        """
        return self.sinh() / self.cosh()

    def acosh(self):
        """
        Return the hyperbolic inverse cosine of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/2
            sage: i = q.cosh(); i
            2.5091784786580567
            sage: abs(i.acosh()-q) < 1e-15
            True
        """
        return self._new_c(gsl_acosh(self._value))

    def arcsinh(self):
        """
        Return the hyperbolic inverse sine of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/2
            sage: i = q.sinh(); i
            2.3012989023072947
            sage: abs(i.arcsinh()-q) < 1e-15
            True
        """
        return self._new_c(gsl_asinh(self._value))

    def arctanh(self):
        """
        Return the hyperbolic inverse tangent of ``self``.

        EXAMPLES::

            sage: q = RDF.pi()/2
            sage: i = q.tanh(); i
            0.9171523356672744
            sage: i.arctanh() - q  # rel tol 1
            4.440892098500626e-16
        """
        return self._new_c(gsl_atanh(self._value))

    def sech(self):
        r"""
        Return the hyperbolic secant of ``self``.

        EXAMPLES::

            sage: RDF(pi).sech()
            0.08626673833405443
            sage: CDF(pi).sech()
            0.08626673833405443
        """
        return 1/self.cosh()

    def csch(self):
        r"""
        Return the hyperbolic cosecant of ``self``.

        EXAMPLES::

            sage: RDF(pi).csch()
            0.08658953753004694
            sage: CDF(pi).csch()  # rel tol 1e-15
            0.08658953753004696
        """
        return 1/self.sinh()

    def coth(self):
        r"""
        Return the hyperbolic cotangent of ``self``.

        EXAMPLES::

            sage: RDF(pi).coth()
            1.003741873197321
            sage: CDF(pi).coth()
            1.0037418731973213
        """
        return self.cosh() / self.sinh()

    def agm(self, other):
        r"""
        Return the arithmetic-geometric mean of ``self`` and ``other``. The
        arithmetic-geometric mean is the common limit of the sequences
        `u_n` and `v_n`, where `u_0` is ``self``,
        `v_0` is other, `u_{n+1}` is the arithmetic mean
        of `u_n` and `v_n`, and `v_{n+1}` is the
        geometric mean of `u_n` and `v_n`. If any operand is negative, the
        return value is ``NaN``.

        EXAMPLES::

            sage: a = RDF(1.5)
            sage: b = RDF(2.3)
            sage: a.agm(b)
            1.8786484558146697

        The arithmetic-geometric mean always lies between the geometric and
        arithmetic mean::

            sage: sqrt(a*b) < a.agm(b) < (a+b)/2
            True
        """
        cdef double a = self._value
        cdef double b = other
        cdef double eps = 2.0**-51
        if a < 0 or b < 0:
            return self._parent.nan()
        while True:
            a1 = (a+b)/2
            b1 = libc.math.sqrt(a*b)
            if abs((b1/a1)-1) < eps: return self._new_c(a1)
            a, b = a1, b1

    def erf(self):
        """
        Return the value of the error function on ``self``.

        EXAMPLES::

            sage: RDF(6).erf()
            1.0
        """
        return self._new_c(gsl_sf_erf(self._value))

    def gamma(self):
        """
        Return the value of the Euler gamma function on ``self``.

        EXAMPLES::

            sage: RDF(6).gamma()
            120.0
            sage: RDF(1.5).gamma()  # rel tol 1e-15
            0.8862269254527584
        """
        sig_on()
        a = self._new_c(gsl_sf_gamma(self._value))
        sig_off()
        return a

    def zeta(self):
        r"""
        Return the Riemann zeta function evaluated at this real number.

        .. NOTE::

           PARI is vastly more efficient at computing the Riemann zeta
           function. See the example below for how to use it.

        EXAMPLES::

            sage: RDF(2).zeta()  # rel tol 1e-15
            1.6449340668482269
            sage: RDF.pi()^2/6
            1.6449340668482264
            sage: RDF(-2).zeta()
            0.0
            sage: RDF(1).zeta()
            +infinity
        """
        if self._value == 1:
            return self._new_c(1)/self._new_c(0)
        return self._new_c(gsl_sf_zeta(self._value))

    def algebraic_dependency(self, n):
        """
        Return a polynomial of degree at most `n` which is
        approximately satisfied by this number.

        .. NOTE::

            The resulting polynomial need not be irreducible, and indeed
            usually won't be if this number is a good approximation to an
            algebraic number of degree less than `n`.

        ALGORITHM:

        Uses the PARI C-library :pari:`algdep` command.

        EXAMPLES::

            sage: r = sqrt(RDF(2)); r
            1.4142135623730951
            sage: r.algebraic_dependency(5)
            x^2 - 2
        """
        return sage.arith.all.algdep(self,n)

    algdep = algebraic_dependency

cdef class ToRDF(Morphism):
    def __init__(self, R):
        """
        Fast morphism from anything with a ``__float__`` method to an ``RDF``
        element.

        EXAMPLES::

            sage: f = RDF.coerce_map_from(ZZ); f
            Native morphism:
              From: Integer Ring
              To:   Real Double Field
            sage: f(4)
            4.0
            sage: f = RDF.coerce_map_from(QQ); f
            Native morphism:
              From: Rational Field
              To:   Real Double Field
            sage: f(1/2)
            0.5
            sage: f = RDF.coerce_map_from(int); f
            Native morphism:
              From: Set of Python objects of class 'int'
              To:   Real Double Field
            sage: f(3r)
            3.0
            sage: f = RDF.coerce_map_from(float); f
            Native morphism:
              From: Set of Python objects of class 'float'
              To:   Real Double Field
            sage: f(3.5)
            3.5
        """
        from sage.categories.homset import Hom
        if isinstance(R, type):
            from sage.sets.pythonclass import Set_PythonType
            R = Set_PythonType(R)
        Morphism.__init__(self, Hom(R, RDF))

    cpdef Element _call_(self, x):
        """
        Send ``x`` to the image under this map.

        EXAMPLES::

            sage: f = RDF.coerce_map_from(float)
            sage: f(3.5) # indirect doctest
            3.5
        """
        cdef RealDoubleElement r = <RealDoubleElement>PY_NEW(RealDoubleElement)
        r._value = PyFloat_AsDouble(x)
        return r

    def _repr_type(self):
        """
        Return the representation type of ``self``.

        EXAMPLES::

            sage: RDF.coerce_map_from(float)._repr_type()
            'Native'
        """
        return "Native"


#####################################################
# unique objects
#####################################################
cdef RealDoubleField_class _RDF
_RDF = RealDoubleField_class()

RDF = _RDF   # external interface

def RealDoubleField():
    """
    Return the unique instance of the
    :class:`real double field<RealDoubleField_class>`.

    EXAMPLES::

        sage: RealDoubleField() is RealDoubleField()
        True
    """
    global _RDF
    return _RDF

def is_RealDoubleElement(x):
    """
    Check if ``x`` is an element of the real double field.

    EXAMPLES::

        sage: from sage.rings.real_double import is_RealDoubleElement
        sage: is_RealDoubleElement(RDF(3))
        True
        sage: is_RealDoubleElement(RIF(3))
        False
    """
    return isinstance(x, RealDoubleElement)


################# FAST CREATION CODE ######################
########### Based on fast integer creation code   #########
######## There is nothing to see here, move along   #######

# We use a global element to steal all the references
# from.  DO NOT INITIALIZE IT AGAIN and DO NOT REFERENCE IT!
cdef RealDoubleElement global_dummy_element
global_dummy_element = RealDoubleElement(0)

# A global pool for performance when elements are rapidly created and destroyed.
# It operates on the following principles:
#
# - The pool starts out empty.
# - When an new element is needed, one from the pool is returned
#   if available, otherwise a new RealDoubleElement object is created
# - When an element is collected, it will add it to the pool
#   if there is room, otherwise it will be deallocated.
DEF element_pool_size = 50

cdef PyObject* element_pool[element_pool_size]
cdef int element_pool_count = 0

# used for profiling the pool
cdef int total_alloc = 0
cdef int use_pool = 0


cdef PyObject* fast_tp_new(type t, args, kwds):
    global element_pool, element_pool_count, total_alloc, use_pool

    cdef PyObject* new

    # for profiling pool usage
    # total_alloc += 1

    # If there is a ready real double in the pool, we will
    # decrement the counter and return that.

    if element_pool_count > 0:

        # for profiling pool usage
        # use_pool += 1

        element_pool_count -= 1
        new = <PyObject *> element_pool[element_pool_count]

    # Otherwise, we have to create one.

    else:

        # allocate enough room for the RealDoubleElement,
        # sizeof_RealDoubleElement is sizeof(RealDoubleElement).
        # The use of PyObject_Malloc directly assumes
        # that RealDoubleElements are not garbage collected, i.e.
        # they do not possess references to other Python
        # objects (As indicated by the Py_TPFLAGS_HAVE_GC flag).
        # See below for a more detailed description.

        new = <PyObject*>PyObject_Malloc( sizeof(RealDoubleElement) )

        # Now set every member as set in z, the global dummy RealDoubleElement
        # created before this tp_new started to operate.

        memcpy(new, (<void*>global_dummy_element), sizeof(RealDoubleElement) )

    # This line is only needed if Python is compiled in debugging mode
    # './configure --with-pydebug' or SAGE_DEBUG=yes. If that is the
    # case a Python object has a bunch of debugging fields which are
    # initialized with this macro.
    if_Py_TRACE_REFS_then_PyObject_INIT(
            new, Py_TYPE(global_dummy_element))

    # The global_dummy_element may have a reference count larger than
    # one, but it is expected that newly created objects have a
    # reference count of one. This is potentially unneeded if
    # everybody plays nice, because the gobal_dummy_element has only
    # one reference in that case.

    # Objects from the pool have reference count zero, so this
    # needs to be set in this case.

    new.ob_refcnt = 1

    return new

cdef void fast_tp_dealloc(PyObject* o):

    # If there is room in the pool for a used integer object,
    # then put it in rather than deallocating it.

    global element_pool, element_pool_count

    if element_pool_count < element_pool_size:

        # And add it to the pool.
        element_pool[element_pool_count] = o
        element_pool_count += 1
        return

    # Free the object. This assumes that Py_TPFLAGS_HAVE_GC is not
    # set. If it was set another free function would need to be
    # called.

    PyObject_Free(o)


from sage.misc.allocator cimport hook_tp_functions
hook_tp_functions(global_dummy_element, <newfunc>(&fast_tp_new), <destructor>(&fast_tp_dealloc), False)


cdef double_repr(double x):
    """
    Convert a double to a string with maximum precision.
    """
    if gsl_finite(x):
        return repr(x)
    cdef int v = gsl_isinf(x)
    if v > 0:
        return "+infinity"
    if v < 0:
        return "-infinity"
    return "NaN"
