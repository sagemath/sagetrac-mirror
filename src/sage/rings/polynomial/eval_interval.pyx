# cython: boundscheck=False, wraparound=False
"""
Evaluation of polynomials in intervals
"""
from __future__ import division

include "sage/ext/interrupt.pxi"
from .polynomial_element cimport Polynomial_generic_dense
from sage.structure.element cimport parent_c as parent, coercion_model
from sage.rings.real_mpfr cimport RealNumber, RealField_class
from sage.rings.real_mpfi import is_RealIntervalField
from sage.rings.complex_interval_field import is_ComplexIntervalField
from sage.libs.mpfr cimport *
from sage.libs.flint.flint cimport FLINT_BIT_COUNT


cpdef eval_polynomial_interval(pol, v):
    """
    Evaluate the polynomial ``pol`` in the interval ``v`` by separating
    the interval in a "center" and "error" part.

    EXAMPLES::

        sage: R.<x> = RIF[]
        sage: pol = x - x^2
        sage: v = RIF(0,1)
        sage: (pol(v)).endpoints()
        (0.000000000000000, 0.250000000000000)

    Calling this function directly::

        sage: from sage.rings.polynomial.eval_interval import eval_polynomial_interval
        sage: eval_polynomial_interval(pol, v).endpoints()
        (0.000000000000000, 0.250000000000000)

    Compare with the naive calculations::

        sage: (v*(1 - v)).endpoints()
        (0.000000000000000, 1.00000000000000)
        sage: (v - v^2).endpoints()
        (-1.00000000000000, 1.00000000000000)

    Test small and big intervals::

        sage: pol(RIF(0.9,1)).endpoints()
        (-6.54858112181245e-17, 0.0925000000000001)
        sage: pol(RIF(-1,2)).endpoints()
        (-2.00000000000000, 0.250000000000000)

    This function is still far from optimal. The Chebyshev polynomials
    for example take values in `[-1,1]` on the interval `[-1,1]`::

        sage: f = chebyshev_T(20, polygen(RIF))
        sage: f(RIF(-1,1)).endpoints()
        (-2.26195350000000e7, 2.26195370000000e7)
        sage: f(RIF(0.9,1)).endpoints()
        (-10.7011727397009, 11.7391595902137)
        sage: f(RIF(0.99,1)).endpoints()
        (-1.39294943457185, 1.00000002117272)
        sage: f(RIF(0.999,1)).endpoints()
        (0.619516553371075, 1.00000001839622)

    Coercion works properly::

        sage: R.<x> = RR[]
        sage: pol = x^3 + 1
        sage: r = pol(RealIntervalField(24)(1,2))
        sage: r.endpoints()
        (0.874999, 9.00001)
        sage: parent(r)
        Real Interval Field with 24 bits of precision
        sage: r = pol(RealIntervalField(64)(1,2))
        sage: r.endpoints()
        (0.874999999999999, 9.00000000000001)
        sage: parent(r)
        Real Interval Field with 53 bits of precision
    """
    V = parent(v)
    F = coercion_model.common_parent(V, pol.base_ring())

    if is_RealIntervalField(F):
        Freal = F
    elif is_ComplexIntervalField(F):
        Freal = F._real_field()
    else:
        raise TypeError("cannot evaluate polynomial with base ring {} at interval".format(pol.base_ring()))

    cdef list coeffs
    try:
        coeffs = (<Polynomial_generic_dense?>pol).__coeffs
    except TypeError:
        coeffs = pol.list()

    cdef Py_ssize_t deg = len(coeffs) - 1
    if deg < 0:
        return F.zero()

    # Convert coefficients and given interval to F if needed
    if parent(coeffs[0]) is not F:
        coeffs = [F(c) for c in coeffs]

    if deg <= 1:
        if deg == 0:
            return coeffs[0]
        else:  # deg == 1
            return (v * coeffs[1]) + coeffs[0]

    # Split interval in center c + error e
    if V is F:
        c = F(v.center())
        e = v - c
    else:
        # Since F is the result after coercion, we know that the
        # precision of F is at most the precision of V.
        # We compute the center with the smaller precision of F, such
        # that the conversions V(c) and F(c) are certainly exact.
        c = F(v).center()
        e = F(v - V(c))
        c = F(c)

    zero = F.zero()
    one = F.one()

    # Compute pol(c)
    cdef Py_ssize_t i
    value = coeffs[deg]
    for i from deg > i >= 0:
        value = value * c + coeffs[i]

    # Error term
    err = eval_polynomial_error(coeffs, c, e)
    # Make err into a box [-err, err]
    unit = Freal(-1, 1)
    if F is not Freal:  # F is complex
        unit = F(unit, unit)
    value += err * unit

    return value


def geom_sum_diff(RealNumber x, RealNumber e, unsigned long n):
    r"""
    Define `(q)[n] = 1 + q + q^2 + ... + q^{n-1}`. Return an
    approximation of `(x+e)[n] - (x)[n]` assuming `x > 0` and `e \geq 0`.

    The approximation is done with the rounding mode of `x`. There is
    no guarantee of exact rounding, but the rounding will be in the
    correct direction if the rounding mode is directed (for example,
    with ``RNDU`` rounding, the returned value is guaranteed to be an
    upper bound).

    ALGORITHM: we use

    .. MATH::

       \frac{(x+e)^n - 1}{(x+e) - 1} - \frac{x^n - 1}{x - 1} = \left(\frac{(x+e)^n - 1}{(x+e) - 1} - n\right) - \left(\frac{x^n - 1}{x - 1} - n\right).

    To evaluate both terms occuring in this formula, we use a
    direct computation if `x \leq 1 - 1/n`. For larger `x`, we use

    .. MATH::

        \frac{x^n - 1}{x - 1} - n = \left(\frac{1}{x-1} + n\right) \left(\exp\left(n \log(x) - \log(1 + n(x-1))\right) - 1 \right).

    TESTS::

        sage: bits = ZZ(sys.maxint).ndigits(2) + 1
        sage: fields = []
        sage: for prec in [20, 40, 80, 160, 240, 320]:
        ....:     for rnd in ["RNDN", "RNDD", "RNDU"]:
        ....:         fields.append(RealField(prec, rnd=rnd))
        sage: RRR = RealField(1200)
        sage: from sage.rings.polynomial.eval_interval import geom_sum_diff
        sage: def test_geom_sum_diff(x0, e0, n):
        ....:     x = RRR(x0)
        ....:     e = RRR(e0)
        ....:     y = x + e
        ....:     b1 = (y^n-1)/(y-1) if y != 1 else n
        ....:     b2 = (x^n-1)/(x-1) if x != 1 else n
        ....:     b = b1 - b2
        ....:     if b.is_infinity() or b.is_NaN():
        ....:         return
        ....:
        ....:     a = geom_sum_diff(x0, e0, n)
        ....:     u = RR((RRR(a) - b)/a.ulp())
        ....:     rnd = parent(x0).rounding_mode()
        ....:     if abs(u) >= 2 or (rnd == "RNDU" and u < 0) or (rnd == "RNDD" and u > 0) or (rnd == "RNDN" and abs(u) >= 1):
        ....:         print("geom_sum_diff({}, {}, {}) with rounding {} wrong by {} ulp".format(x0.str(truncate=False), e0.str(truncate=False), n, rnd, u))
        sage: for i in range(100000):  # long time
        ....:     R1 = fields[randrange(len(fields))]
        ....:     R2 = fields[randrange(len(fields))]
        ....:     x = R1.random_element(-100, 100).exp()
        ....:     e = R2.random_element(-100, 100).exp()
        ....:     k = randrange(10)
        ....:     if k == 0:         # numbers close to 1 are bad cases
        ....:         x = abs(1 - x)
        ....:     elif k == 1:
        ....:         e = abs(1 - x)
        ....:     elif k == 2:
        ....:         x += 1
        ....:     elif k == 3:
        ....:         e += 1
        ....:     n = randrange(5) + 2^randrange(bits)
        ....:     test_geom_sum_diff(x, e, n)
    """
    cdef RealField_class R = <RealField_class>(x._parent)

    # Trivial cases
    if mpfr_zero_p(e.value) or n <= 1:
        return R.zero()
    if n == 2:
        return R(e)

    # Handle NaN/infinity
    if not mpfr_number_p(x.value) or not mpfr_number_p(e.value):
        if mpfr_nan_p(x.value):
            return R(x)
        if not mpfr_number_p(e.value):
            return R(e)
        return R(x)

    cdef RealNumber ret = R._new()
    cdef mp_prec_t prec = R.__prec
    cdef mpfr_rnd_t rnd, rndopp
    cdef mp_exp_t exp_1x

    sig_on()

    # Compute exp_1x = largest exponent of 1 and x
    exp_1x = mpfr_get_exp(x.value)
    if exp_1x < 0:
        exp_1x = 0

    # We need enough precision for max(1,x,e) + e
    cdef mp_prec_t prec2 = exp_1x - mpfr_get_exp(e.value)
    if prec2 > 0:
        prec += prec2

    # Given that x and e are positive, the answer also must be
    # positive. We compute the rounding mode and its opposite.
    if R.rnd == MPFR_RNDN:
        rnd = MPFR_RNDN
        rndopp = MPFR_RNDN
    if R.rnd == MPFR_RNDD or R.rnd == MPFR_RNDZ:
        rnd = MPFR_RNDD
        rndopp = MPFR_RNDU
    else:
        rnd = MPFR_RNDU
        rndopp = MPFR_RNDD

    prec += FLINT_BIT_COUNT(n)
    # Additional headroom of 10 bits
    prec += 10

    cdef mpfr_t r1, r2, t, u
    mpfr_init2(r1, prec)
    mpfr_init2(r2, prec)
    mpfr_init2(t, prec)
    mpfr_init2(u, prec)

    mpfr_add(u, x.value, e.value, rnd)  # u = x + e
    mpfr_geom_sum_minus_n(r1, u, t, u, n, rnd, rndopp)
    mpfr_geom_sum_minus_n(r2, x.value, t, u, n, rndopp, rnd)
    mpfr_sub(ret.value, r1, r2, R.rnd)

    mpfr_clear(r1)
    mpfr_clear(r2)
    mpfr_clear(t)
    mpfr_clear(u)
    sig_off()

    return ret


cdef inline void mpfr_geom_sum_minus_n(mpfr_t rop, mpfr_t x, mpfr_t t, mpfr_t u, unsigned long n, mpfr_rnd_t rnd, mpfr_rnd_t rndopp):
    r"""
    Compute `(x^n - 1)/(x - 1) - n`.
    The result is stored in ``rop``, using ``t`` and ``u`` as scratch
    variables. None of the ``mpfr_t`` variables may alias another,
    except for ``x`` and ``u``.

    The result is rounded in the direction given by ``rnd``, which must
    be ``MPFR_RNDN``, ``MPFR_RNDD`` or ``MPFR_RNDU``.
    The argument ``rndopp`` denotes the opposite rounding mode.
    """
    # Let rop = x - 1 (with "rndopp" rounding)
    cdef inex = mpfr_sub_ui(rop, x, 1, rndopp)

    if mpfr_zero_p(rop):  # x == 1 and rop == 0
        return

    # First compute n (x - 1)
    mpfr_mul_ui(t, rop, n, rndopp)

    # If this is <= -1 (then in particular x < 1), use the formula
    # (x^n - 1)/(x - 1) - n
    # Mind the rounding because x < 1
    if mpfr_cmp_si(t, -1) <= 0:
        mpfr_pow_ui(t, x, n, rndopp)
        mpfr_sub_ui(t, t, 1, rndopp)
        if inex:  # We need x-1 with rounding rnd
            mpfr_sub_ui(rop, x, 1, rnd)
        mpfr_div(u, t, rop, rnd)
        mpfr_sub_ui(rop, u, n, rnd)
    else:
        # The result is A * B where
        # A = 1/(x-1) + n
        # B = exp(n log(x) - log(1+n(x-1))) - 1
        # where both A and B are > 0
        mpfr_log(u, x, rnd)
        mpfr_mul_ui(u, u, n, rnd)
        mpfr_log1p(t, t, rndopp)
        mpfr_sub(u, u, t, rnd)
        mpfr_expm1(u, u, rnd)
        mpfr_ui_div(t, 1, rop, rnd)
        mpfr_add_ui(t, t, n, rnd)
        mpfr_mul(rop, u, t, rnd)


cpdef eval_polynomial_error(list coeffs, x, e):
    r"""
    Return a bound on `|f(x+e) - f(x)|`.

    INPUT:

    - ``coeffs`` -- a list of intervals representing the coefficients
      of `f`

    - ``x`` -- an interval, assumed to be an exact number

    - ``e`` -- an interval for the error term

    EXAMPLES::

        sage: from sage.rings.polynomial.eval_interval import eval_polynomial_error
        sage: R.<x> = RIF[]
        sage: pol = x^3
        sage: eval_polynomial_error(list(pol), RIF(10), RIF(1))
        331.000000012233
        sage: pol(11) - pol(10)
        331

    If the given polynomial is a geometric series, the bound should be
    sharp (up to rounding errors)::

        sage: pol = x^5 + x^4 + x^3 + x^2 + x + 1
        sage: eval_polynomial_error(list(pol), RIF(1), RIF(1))
        57.0000000000001
        sage: pol(2) - pol(1)
        57
        sage: eval_polynomial_error(list(pol), RIF(2), RIF(-1,1))
        301.000000000001
        sage: pol(3) - pol(2)
        301
        sage: eval_polynomial_error(list(pol), RIF(1/10), RIF(2))
        75.9490000000001
        sage: pol(1/10+2) - pol(1/10)
        75.949000000000?
        sage: eval_polynomial_error(list(pol), RIF(10), RIF(0))
        0.000000000000000

    ::

        sage: pol = (x^3 + x^2 + x + 1)(x*28); pol
        21952*x^3 + 784*x^2 + 28*x + 1
        sage: eval_polynomial_error(list(pol), RIF(3), RIF(1))
        817740.000000002
        sage: pol(4)-pol(3)
        817740
        sage: pol = (x^3 + x^2 + x + 1)(x/28); pol
        0.0000455539358600583?*x^3 + 0.001275510204081633?*x^2 + 0.03571428571428572?*x + 1
        sage: eval_polynomial_error(list(pol), RIF(3), RIF(1))
        0.0463283527696794
        sage: pol(4)-pol(3)
        0.046328352769679?

    Some trivial cases::

        sage: eval_polynomial_error([], RIF(1234), RIF(10))
        0.000000000000000
        sage: eval_polynomial_error([RIF(1)], RIF(1234), RIF(10))
        0.000000000000000
        sage: eval_polynomial_error([RIF(1), RIF(1)], RIF(1234), RIF(10))
        10
    """
    cdef Py_ssize_t i, deg = len(coeffs)-1
    emag = e.magnitude()
    R = parent(emag)

    # Trivial cases
    if deg <= 0 or not emag:
        return R.zero()
    if deg == 1:
        return coeffs[1] * e.magnitude()

    # Special case x == 0: just compute f(e) - f(0)
    if x.contains_zero():
        epow = e
        res = coeffs[1] * epow
        for i in range(2, deg+1):
            epow *= e
            res += coeffs[i] * epow
        return res.magnitude()

    # List of log |c_i| (i >= 1)
    cdef list logc = []
    for i in range(1, deg+1):
        logc.append(coeffs[i].magnitude().log())

    # Compute A and B such that |c_i| <= B * A^i (for i >= 1)
    # In other words: logc_i <= logB + (i+1) * logA

    # Smallest allowed value for log(A) to avoid under/overflow
    logmin = R(-1073741824.0)   # -2^30

    # We first try A = 1/|x|
    logA = -x.magnitude().log()
    if logA < logmin:
        logA = logmin
    logB = logc[0] - logA
    cdef Py_ssize_t i0 = 0   # Index where equality was reached
    for i in range(1, deg):
        t = logc[i] + R(-(i+1)) * logA
        if t > logB:
            logB = t
            i0 = i

    # If the best index was the first or the last, we can improve our
    # bound.
    if i0 == 0:
        # C = B * A
        # We need logc_i <= logC + i * logA
        logC = logc[0]
        logA = logc[1] - logC
        for i in range(2, deg):
            t = (logc[i] - logC) / R(i)
            if t > logA:
                logA = t
        if logA < logmin:
            logA = logmin
        logB = logC - logA
    elif i0 == deg - 1:
        # C = B * A^deg
        # We need logc_i <= logC + (deg-(i+1)) * (-logA)
        logC = logc[deg-1]
        mlogA = logc[deg-2] - logC  # -log(A)
        for i in range(0, deg-2):
            t = (logc[i] - logC) / R(deg-(i+1))
            if t > mlogA:
                mlogA = t
        if mlogA < logmin:
            mlogA = logmin
        logA = -mlogA
        logB = logC + R(-deg) * logA

    # Now an upper bound for the error is given by the error
    # for the function x -> B * ((Ax)^(deg+1) - 1)/(Ax - 1)
    A = logA.exp()
    return logB.exp() * geom_sum_diff(A*x.magnitude(), A*e.magnitude(), deg+1)
