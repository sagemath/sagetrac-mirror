r"""
Univariate rational functions over `\QQ` implemented via FLINT
"""

include 'sage/ext/interrupt.pxi'

from cpython.object cimport Py_LE, Py_EQ, Py_GE, Py_NE

from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport mpq_numref, mpq_denref

from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.fmpq_poly cimport fmpq_poly_get_numerator, fmpq_poly_denref

from sage.rings.rational cimport Rational
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_rational_flint cimport Polynomial_rational_flint
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint

from sage.structure.element cimport Element, ModuleElement, RingElement, FieldElement

cdef int _init_fmpz_poly(mpz_t lcm, fmpz_poly_t ret, op) except -1:
    cdef long i
    cdef mpz_t temp

    if isinstance(op, (list, tuple)):
        if len(op) == 0:
            fmpz_poly_zero(ret)
            mpz_set_ui(lcm, 1u)
            return 0
        elif len(op) == 1:
            return _init_fmpz_poly(lcm, ret, op[0])

        op = [coeff if isinstance(coeff, (Rational, Integer, int))
                else Rational(coeff) for coeff in op]

        sig_on()

        mpz_set_ui(lcm, 1u)
        for coeff in op:
            if isinstance(coeff, Rational):
                mpz_lcm(lcm, lcm, mpq_denref((<Rational>coeff).value))

        sig_off()

        mpz_init(temp)

        if not sig_on_no_except():
            mpz_clear(temp)
            cython_check_exception()

        fmpz_poly_realloc(ret, len(op))

        for i, coeff in enumerate(op):

            if isinstance(coeff, Rational):
                mpz_divexact(temp, lcm, mpq_denref((<Rational>coeff).value))
                mpz_mul(temp, temp, mpq_numref((<Rational>coeff).value))
            elif isinstance(coeff, Integer):
                mpz_mul(temp, lcm, (<Integer>coeff).value)
            else: # isinstance(coeff, int)
                mpz_mul_si(temp, lcm, coeff)

            fmpz_poly_set_coeff_mpz(ret, i, temp)

        sig_off()

        mpz_clear(temp)

    elif isinstance(op, Rational):
        fmpz_poly_set_mpz(ret, mpq_numref((<Rational>op).value))
        mpz_set(lcm, mpq_denref((<Rational>op).value))
    elif isinstance(op, Polynomial_rational_flint):
        sig_on()
        fmpq_poly_get_numerator(ret, (<Polynomial_rational_flint>op).__poly)
        fmpz_get_mpz(lcm, fmpq_poly_denref((<Polynomial_rational_flint>op).__poly))
        sig_off()
    else:
        if isinstance(op, int):
            fmpz_poly_set_si(ret, op)
        elif isinstance(op, Integer):
            fmpz_poly_set_mpz(ret, (<Integer>op).value)
        elif isinstance(op, Polynomial_integer_dense_flint):
            sig_on()
            fmpz_poly_set(ret, (<Polynomial_integer_dense_flint>op).__poly)
            sig_off()
        else:
            raise ValueError
        mpz_set_ui(lcm, 1u)

cdef class RationalFunctionQQ(FieldElement):
    """
    Univariate rational functions over the rationals, implemented via FLINT.
    """

    ###########################################################################
    # Allocation & initialisation                                             #
    ###########################################################################

    def __cinit__(self):
        fmpz_poly_q_init(self._func)

    def __dealloc__(self):
        fmpz_poly_q_clear(self._func)

    cdef RationalFunctionQQ _new(self):
        cdef RationalFunctionQQ res
        res = RationalFunctionQQ.__new__(RationalFunctionQQ)
        res._parent = self._parent
        return res

    def __init__(self, parent, numerator, denominator = 1):

        if denominator == 0:
            raise ValueError

        cdef mpz_t d_num, d_den, g

        mpz_init(d_num)
        mpz_init(d_den)
        mpz_init(g)

        if not sig_on_no_except():
            mpz_clear(d_num)
            mpz_clear(d_den)
            mpz_clear(g)
            cython_check_exception()

        _init_fmpz_poly(d_num, fmpz_poly_q_numref(self._func), numerator)
        _init_fmpz_poly(d_den, fmpz_poly_q_denref(self._func), denominator)

        mpz_gcd(g, d_num, d_den)
        mpz_divexact(d_num, d_num, g)
        mpz_divexact(d_den, d_den, g)

        fmpz_poly_q_scalar_mul_mpz(self._func, self._func, d_den)
        fmpz_poly_q_scalar_div_mpz(self._func, self._func, d_num)

        fmpz_poly_q_canonicalise(self._func)

        sig_off()

        mpz_clear(d_num)
        mpz_clear(d_den)
        mpz_clear(g)

        self._parent = parent

    def _repr_(self):
        return fmpz_poly_q_get_python_str(
                self._func, name=self._parent.variable_name())

    def _latex_(self):
        return fmpz_poly_q_get_python_str(
                self._func, name=self._parent.variable_name(), latex=True)

    def __copy__(self):
        cdef RationalFunctionQQ res = self._new()
        sig_on()
        fmpz_poly_q_set(res._func, self._func)
        sig_off()
        return res

    ###########################################################################
    # Basis access                                                            #
    ###########################################################################

    def degree(self):
        cdef Integer deg = Integer.__new__(Integer)
        mpz_set_si(deg.value,
                fmpz_poly_degree(fmpz_poly_q_numref(self._func))
                    - fmpz_poly_degree(fmpz_poly_q_denref(self._func))
                )
        return deg

    def __call__(self, value):
        cdef Rational qq_res
        cdef RationalFunctionQQ func_res

        try:
            value = Rational(value)
        except TypeError:
            pass

        if isinstance(value, Rational):
            qq_res = Rational.__new__(Rational)
            sig_on()
            fmpz_poly_q_evaluate(
                    qq_res.value, self._func, (<Rational> value).value)
            sig_off()
            return qq_res

        if isinstance(value, RationalFunctionQQ):
            # TODO: composition
            pass

        # TODO: general coercion
        raise ValueError

    ###########################################################################
    # Comparisons                                                             #
    ###########################################################################

    def is_zero(self):
        return bool(fmpz_poly_q_is_zero(self._func))

    def __nonzero__(self):
        return not fmpz_poly_q_is_zero(self._func)

    cdef _richcmp_c_impl(left, Element right, int op):
        if op in (Py_LE, Py_EQ, Py_GE):
            return bool(fmpz_poly_q_equal(
                    left._func, (<RationalFunctionQQ>right)._func))
        if op == Py_NE:
            return not fmpz_poly_q_equal(
                    left._func, (<RationalFunctionQQ>right)._func)
        return False

    ###########################################################################
    # Arithmetic                                                              #
    ###########################################################################

    cpdef ModuleElement _add_(left, ModuleElement right):
        cdef RationalFunctionQQ res = left._new()

        sig_on()
        fmpz_poly_q_add(res._func,
                left._func, (<RationalFunctionQQ> right)._func)
        sig_off()

        return res

    cpdef ModuleElement _iadd_(self, ModuleElement right):
        sig_on()
        fmpz_poly_q_add(self._func,
                self._func, (<RationalFunctionQQ> right)._func)
        sig_off()

        return self

    cpdef ModuleElement _sub_(left, ModuleElement right):
        cdef RationalFunctionQQ res = left._new()

        sig_on()
        fmpz_poly_q_sub(res._func,
                left._func, (<RationalFunctionQQ> right)._func)
        sig_off()

        return res

    cpdef ModuleElement _isub_(self, ModuleElement right):
        sig_on()
        fmpz_poly_q_sub(self._func,
                self._func, (<RationalFunctionQQ> right)._func)
        sig_off()

        return self

    cpdef ModuleElement _neg_(self):
        cdef RationalFunctionQQ res = self._new()

        sig_on()
        fmpz_poly_q_neg(res._func, self._func)
        sig_off()

        return res

    cpdef RingElement _mul_(left, RingElement right):
        cdef RationalFunctionQQ res = left._new()

        sig_on()
        fmpz_poly_q_mul(res._func,
                left._func, (<RationalFunctionQQ> right)._func)
        sig_off()

        return res

    cpdef RingElement _imul_(self, RingElement right):
        sig_on()
        fmpz_poly_q_mul(self._func,
                self._func, (<RationalFunctionQQ> right)._func)
        sig_off()

        return self

    cpdef RingElement _div_(left, RingElement right):
        if not right:
            raise ZeroDivisionError("division by zero")

        cdef RationalFunctionQQ res = left._new()

        sig_on()
        fmpz_poly_q_div(res._func,
                left._func, (<RationalFunctionQQ> right)._func)
        sig_off()

        return res

    cpdef RingElement _idiv_(self, RingElement right):
        if not right:
            raise ZeroDivisionError("division by zero")

        sig_on()
        fmpz_poly_q_div(self._func,
                self._func, (<RationalFunctionQQ> right)._func)
        sig_off()

        return self

    def __invert__(self):
        if not self:
            raise ZeroDivisionError("division by zero")

        cdef RationalFunctionQQ res = self._new()

        sig_on()
        fmpz_poly_q_inv(res._func, self._func)
        sig_off()

        return res

    cpdef ModuleElement _rmul_(self, RingElement left):
        cdef RationalFunctionQQ res = self._new()

        sig_on()
        fmpz_poly_q_scalar_mul_mpq(res._func,
                self._func, (<Rational> left).value)
        sig_off()

        return res

    cpdef ModuleElement _lmul_(self, RingElement right):
        return self._rmul_(right)

    def __pow__(RationalFunctionQQ self, exp, ignored):
        cdef RationalFunctionQQ res
        cdef fmpz_poly_struct *poly
        cdef int shift

        res = self._new()

        if exp == 0:
            fmpz_poly_q_one(res._func)
            return res
        elif exp < 0:
            if not self:
                raise ZeroDivisionError("negative exponent in power of zero")

            sig_on()
            fmpz_poly_q_inv(res._func, self._func)
            sig_off()

            exp = -exp
        else:
            sig_on()
            fmpz_poly_q_set(res._func, self._func)
            sig_off()

        if fmpz_poly_q_is_zero(res._func) or fmpz_poly_q_is_one(res._func):
            return res

        sig_on()

        # optimize common exponentiations (e.g. x^n, (1/x)^n)
        poly = fmpz_poly_q_numref(res._func)
        for shift in range(fmpz_poly_length(poly)):
            if not fmpz_is_zero(fmpz_poly_get_coeff_ptr(poly, shift)):
                break

        if not shift:
            poly = fmpz_poly_q_denref(res._func)
            for shift in range(fmpz_poly_length(poly)):
                if not fmpz_is_zero(fmpz_poly_get_coeff_ptr(poly, shift)):
                    break

        if shift:
            fmpz_poly_shift_right(poly, poly, shift)

        if not fmpz_poly_q_is_one(res._func):
            fmpz_poly_q_pow(res._func, res._func, <int> exp)

        if shift:
            fmpz_poly_shift_left(poly, poly, <int>(shift*exp))

        sig_off()

        return res
