"""
Finite field elements implemented via FLINT's fq_t type

AUTHORS:

- Jean-Pierre Flori (July 2014): initial version, based on
  element_pari_ffelt.py by Peter Bruin.
"""

#*****************************************************************************
#      Copyright (C) 2013 Peter Bruin <peter.bruin@math.uzh.ch>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"
include "sage/libs/pari/decl.pxi"
include "sage/libs/pari/pari_err.pxi"

from sage.libs.flint.types cimport *
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.fmpz_mod_poly cimport *
from sage.libs.flint.fq cimport *

from element_base cimport FinitePolyExtElement
from integer_mod import IntegerMod_abstract

import sage.rings.integer
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.finite_field_flint_fq cimport FiniteField_flint_fq
from sage.rings.finite_rings.element_pari_ffelt cimport FiniteFieldElement_pari_ffelt
from sage.libs.pari.pari_instance cimport PariInstance
cdef PariInstance pari = sage.libs.pari.pari_instance.pari
from sage.libs.pari.gen cimport gen as pari_gen
from sage.interfaces.gap import is_GapElement
from sage.modules.free_module_element import FreeModuleElement
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.multi_polynomial_element import MPolynomial
from sage.rings.rational import Rational
from sage.structure.element cimport Element, ModuleElement, RingElement

cdef class FiniteFieldElement_flint_fq(FinitePolyExtElement):
    """
    An element of a finite field.

    EXAMPLE::

        sage: K = FiniteField(10007^10, 'a', impl='flint_fq')
        sage: a = K.gen(); a
        a
        sage: type(a)
        <type 'sage.rings.finite_rings.element_flint_fq.FiniteFieldElement_flint_fq'>

    TESTS::

        sage: n = 63
        sage: m = 3;
        sage: K.<a> = GF(2^n, impl='flint_fq')
        sage: f = conway_polynomial(2, n)
        sage: f(a) == 0
        True
        sage: e = (2^n - 1) / (2^m - 1)
        sage: conway_polynomial(2, m)(a^e) == 0
        True

        sage: K.<a> = FiniteField(2^16, impl='flint_fq')
        sage: K(0).is_zero()
        True
        sage: (a - a).is_zero()
        True
        sage: a - a
        0
        sage: a == a
        True
        sage: a - a == 0
        True
        sage: a - a == K(0)
        True
        sage: TestSuite(a).run()

    Test creating elements from basic Python types::

        sage: K.<a> = FiniteField(7^20, impl='flint_fq')
        sage: K(int(8))
        1
        sage: K(long(-2^300))
        6
    """

    def __init__(FiniteFieldElement_flint_fq self, object parent, object x):
        """
        Initialise ``self`` with the given ``parent`` and value
        converted from ``x``.

        This is called when constructing elements from Python.

        TEST::

            sage: from sage.rings.finite_rings.element_flint_fq import FiniteFieldElement_flint_fq
            sage: K = FiniteField(101^2, 'a', impl='flint_fq')
            sage: x = FiniteFieldElement_flint_fq(K, 'a + 1')
            sage: x
            a + 1
        """
        # FinitePolyExtElement.__init__(self, parent)
        self._parent = parent
        self._cparent = (<FiniteField_flint_fq>parent)._ctx
        # cannot do that in cinit as we need cparent
        fq_init(self.val, self._cparent)
        self.initialized = 1
        self.construct_from(x)

    def __dealloc__(FiniteFieldElement_flint_fq self):
        """
        Cython deconstructor.
        """
        if self.initialized:
            fq_clear(self.val, self._cparent)

    cdef FiniteFieldElement_flint_fq _new(FiniteFieldElement_flint_fq self):
        """
        Create an empty element with the same parent as ``self``.
        """
        cdef FiniteFieldElement_flint_fq x
        x = FiniteFieldElement_flint_fq.__new__(FiniteFieldElement_flint_fq)
        x._parent = self._parent
        x._cparent = self._cparent
        fq_init(x.val, x._cparent)
        self.initialized = 1
        return x

    cdef void set_from_fq(FiniteFieldElement_flint_fq self, fq_t val) except *:
        """
        Sets self from an fq_t.
        """
        fq_set(self.val, val, self._cparent)

    cdef void construct_from(FiniteFieldElement_flint_fq self, object x) except *:
        """
        Initialise ``self`` from the Sage object `x`.
        """
        cdef Integer x_INT
        cdef GEN x_GEN
        cdef fmpz_t x_flint

        if isinstance(x, FiniteFieldElement_flint_fq):
            if self._parent is (<FiniteFieldElement_flint_fq>x)._parent:
                fq_set(self.val, (<FiniteFieldElement_flint_fq>x).val, self._cparent)
            else:
                raise TypeError("no coercion defined")

        elif isinstance(x, Integer):
            fmpz_init_set_readonly(x_flint, (<Integer>x).value)
            fq_set_fmpz(self.val, x_flint, self._cparent)
            fmpz_clear_readonly(x_flint)

        elif isinstance(x, int) or isinstance(x, long):
            x_INT = Integer(x)
            self.construct_from(x_INT)

        elif isinstance(x, IntegerMod_abstract):
            if self._parent.characteristic().divides(x.modulus()):
                x_INT = Integer(x)
                fmpz_init_set_readonly(x_flint, x_INT.value)
                fq_set_fmpz(self.val, x_flint, self._cparent)
                fmpz_clear_readonly(x_flint)
            else:
                raise TypeError("no coercion defined")

        elif x is None:
            fq_zero(self.val, self._cparent)

        #elif isinstance(x, FiniteFieldElement_base):
        #    raise NotImplementedError

        elif isinstance(x, pari_gen):
            x_GEN = (<pari_gen>x).g

            pari_catch_sig_on()
            if gequal0(x_GEN):
                fq_zero(self.val, self._cparent)
                pari_catch_sig_off()
            elif gequal1(x_GEN):
                fq_one(self.val, self._cparent)
                pari_catch_sig_off()
            else:
                t = typ(x_GEN)
                p = self._parent.characteristic()
                f = self._parent.degree()
                modulus = self._parent.modulus().change_ring(ZZ)
                if t == t_FFELT:
                    # The following calls pari_catch_sig_off()
                    x = (<PariInstance>pari).new_gen(FF_p(x_GEN))
                    pari_catch_sig_on()
                    y = FF_f(x_GEN)
                    R = modulus.parent()
                    pari_catch_sig_on()
                    z = R((<PariInstance>pari).new_gen(FF_mod(x_GEN)))
                    if x == p and y == f and z == modulus:
                        # The following calls pari_catch_sig_off()
                        x = (<PariInstance>pari).new_gen(FF_to_FpXQ(x_GEN))
                        self.construct_from(R(x))
                        return
                elif t == t_INT:
                    pari_catch_sig_off()
                    self.construct_from(Integer(x))
                    return
                elif t == t_INTMOD:
                    # The following calls pari_catch_sig_off()
                    x = (<PariInstance>pari).new_gen(gel(x_GEN, 1))
                    if x % p == 0:
                        self.construct_from(Integer(x))
                        return
                elif t == t_FRAC:
                    # The following calls pari_catch_sig_off()
                    x = (<PariInstance>pari).new_gen(gel(x_GEN, 2))
                    if x % p != 0:
                        self.construct_from(Rational(x))
                        return
                else:
                    pari_catch_sig_off()
                raise TypeError("no coercion defined")

        elif (isinstance(x, FreeModuleElement)
              and x.parent() is self._parent.vector_space()):
            n = len(x)
            while n > 0 and x[n - 1] == 0:
                n -= 1
            if n == 0:
                fq_zero(self.val, self._cparent)
            else:
                fmpz_poly_zero(self.val)
                for i in xrange(n):
                    x_INT = Integer(x[i])
                    fmpz_init_set_readonly(x_flint, x_INT.value)
                    fmpz_poly_set_coeff_fmpz(self.val, i, x_flint)
                    fmpz_clear_readonly(x_flint)
                fq_reduce(self.val, self._cparent)

        elif isinstance(x, Rational):
            self.construct_from(x % self._parent.characteristic())

        elif isinstance(x, Polynomial):
            if x.base_ring() is not self._parent.base_ring():
                x = x.change_ring(self._parent.base_ring())
            self.construct_from(x.substitute(self._parent.gen()))

        elif isinstance(x, MPolynomial) and x.is_constant():
            self.construct_from(x.constant_coefficient())

        elif isinstance(x, list):
            if len(x) == self._parent.degree():
                self.construct_from(self._parent.vector_space()(x))
            else:
                Fp = self._parent.base_ring()
                self.construct_from(self._parent.polynomial_ring()([Fp(y) for y in x]))

        elif isinstance(x, str):
            self.construct_from(self._parent.polynomial_ring()(x))

        elif is_GapElement(x):
            from sage.interfaces.gap import gfq_gap_to_sage
            try:
                self.construct_from(gfq_gap_to_sage(x, self._parent))
            except (ValueError, IndexError, TypeError):
                raise TypeError("no coercion defined")

        else:
            raise TypeError("no coercion defined")

    def _repr_(FiniteFieldElement_flint_fq self):
        """
        Return the string representation of ``self``.

        EXAMPLE::

            sage: k.<c> = GF(3^17, impl='flint_fq')
            sage: c^20  # indirect doctest
            c^4 + 2*c^3
        """
        return (<bytes> fq_get_str_pretty(self.val, self._cparent)).replace('+', ' + ')

    def __hash__(FiniteFieldElement_flint_fq self):
        """
        Return the hash of ``self``.  This is by definition equal to
        the hash of ``self.polynomial()``.

        EXAMPLE::

            sage: k.<a> = GF(3^15, impl='flint_fq')
            sage: R = GF(3)['a']; aa = R.gen()
            sage: hash(a^2 + 1) == hash(aa^2 + 1)
            True
        """
        return hash(self.polynomial())

    def __reduce__(FiniteFieldElement_flint_fq self):
        """
        For pickling.

        TEST::

            sage: K.<a> = FiniteField(10007^10, impl='flint_fq')
            sage: loads(a.dumps()) == a
            True
        """
        return unpickle_FiniteFieldElement_flint_fq, (self._parent, str(self))

    def __copy__(FiniteFieldElement_flint_fq self):
        """
        Return a copy of ``self``.

        TESTS::

            sage: k.<a> = FiniteField(3^3, impl='flint_fq')
            sage: a
            a
            sage: b = copy(a); b
            a
            sage: a == b
            True
            sage: a is b
            False
        """
        cdef FiniteFieldElement_flint_fq x = self._new()
        fq_set(x.val, self.val, x._cparent)
        return x

    cdef int _cmp_c_impl(FiniteFieldElement_flint_fq self, Element other) except -2:
        """
        Comparison of finite field elements.

        TESTS::

            sage: k.<a> = FiniteField(3^3, impl='flint_fq')
            sage: a == 1
            False
            sage: a^0 == 1
            True
            sage: a == a
            True
            sage: a < a^2
            True
            sage: a > a^2
            False
        """
        if fq_equal(self.val, (<FiniteFieldElement_flint_fq>other).val, self._cparent):
            return 0
        else:
            r = cmp(fmpz_poly_degree(self.val), fmpz_poly_degree((<FiniteFieldElement_flint_fq>other).val))
            if r:
                return r
            return cmp(self.polynomial(), other.polynomial())

    def __richcmp__(FiniteFieldElement_flint_fq left, object right, int op):
        """
        Rich comparison of finite field elements.

        EXAMPLE::

            sage: k.<a> = GF(2^20, impl='flint_fq')
            sage: e = k.random_element()
            sage: f = loads(dumps(e))
            sage: e is f
            False
            sage: e == f
            True
            sage: e != (e + 1)
            True

        .. NOTE::

            Finite fields are unordered.  However, for the purpose of
            this function, we adopt the lexicographic ordering on the
            representing polynomials.

        EXAMPLE::

            sage: K.<a> = GF(2^100, impl='flint_fq')
            sage: a < a^2
            True
            sage: a > a^2
            False
            sage: a+1 > a^2
            False
            sage: a+1 < a^2
            True
            sage: a+1 < a
            False
            sage: a+1 == a
            False
            sage: a == a
            True
        """
        return (<Element>left)._richcmp(right, op)

    cpdef ModuleElement _add_(FiniteFieldElement_flint_fq self, ModuleElement right):
        """
        Addition.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='flint_fq')
            sage: a + a^2 # indirect doctest
            a^2 + a
        """
        cdef FiniteFieldElement_flint_fq x = self._new()
        fq_add(x.val, (<FiniteFieldElement_flint_fq>self).val,
                      (<FiniteFieldElement_flint_fq>right).val,
                      self._cparent)
        return x

    cpdef ModuleElement _sub_(FiniteFieldElement_flint_fq self, ModuleElement right):
        """
        Subtraction.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='flint_fq')
            sage: a - a # indirect doctest
            0
        """
        cdef FiniteFieldElement_flint_fq x = self._new()
        fq_sub(x.val, (<FiniteFieldElement_flint_fq>self).val,
                      (<FiniteFieldElement_flint_fq>right).val,
                      self._cparent)
        return x

    cpdef RingElement _mul_(FiniteFieldElement_flint_fq self, RingElement right):
        """
        Multiplication.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='flint_fq')
            sage: (a^12 + 1)*(a^15 - 1) # indirect doctest
            a^15 + 2*a^12 + a^11 + 2*a^10 + 2
        """
        cdef FiniteFieldElement_flint_fq x = self._new()
        fq_mul(x.val, (<FiniteFieldElement_flint_fq>self).val,
                      (<FiniteFieldElement_flint_fq>right).val,
                      self._cparent)
        return x

    cpdef RingElement _div_(FiniteFieldElement_flint_fq self, RingElement right):
        """
        Division.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='flint_fq')
            sage: (a - 1) / (a + 1) # indirect doctest
            2*a^16 + a^15 + 2*a^14 + a^13 + 2*a^12 + a^11 + 2*a^10 + a^9 + 2*a^8 + a^7 + 2*a^6 + a^5 + 2*a^4 + a^3 + 2*a^2 + a + 1
        """
        if fq_is_zero((<FiniteFieldElement_flint_fq>right).val, self._cparent):
            raise ZeroDivisionError
        cdef FiniteFieldElement_flint_fq x = self._new()
        cdef fq_t rinv
        fq_init(rinv, self._cparent)
        fq_inv(rinv, (<FiniteFieldElement_flint_fq>right).val, self._cparent)
        fq_mul(x.val, (<FiniteFieldElement_flint_fq>self).val,
                      rinv, self._cparent)
        fq_clear(rinv, self._cparent)
        # Not in the 2.4 branch
        #fq_div(x.val, (<FiniteFieldElement_flint_fq>self).val,
        #              (<FiniteFieldElement_flint_fq>right).val,
        #              self._cparent)
        return x

    def is_zero(FiniteFieldElement_flint_fq self):
        """
        Return ``True`` if ``self`` equals 0.

        EXAMPLE::

            sage: F.<a> = FiniteField(5^3, impl='flint_fq')
            sage: a.is_zero()
            False
            sage: (a - a).is_zero()
            True
        """
        return bool(fq_is_zero(self.val, self._cparent))

    def is_one(FiniteFieldElement_flint_fq self):
        """
        Return ``True`` if ``self`` equals 1.

        EXAMPLE::

            sage: F.<a> = FiniteField(5^3, impl='flint_fq')
            sage: a.is_one()
            False
            sage: (a/a).is_one()
            True
        """
        return bool(fq_is_one(self.val, self._cparent))

    def is_unit(FiniteFieldElement_flint_fq self):
        """
        Return ``True`` if ``self`` is non-zero.

        EXAMPLE::

            sage: F.<a> = FiniteField(5^3, impl='flint_fq')
            sage: a.is_unit()
            True
        """
        return not bool(fq_is_zero(self.val, self._cparent))

    __nonzero__ = is_unit

    def __pos__(FiniteFieldElement_flint_fq self):
        """
        Unitary positive operator...

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='flint_fq')
            sage: +a
            a
        """
        return self

    def __neg__(FiniteFieldElement_flint_fq self):
        """
        Negation.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='flint_fq')
            sage: -a
            2*a
        """
        cdef FiniteFieldElement_flint_fq x = self._new()
        fq_neg(x.val, (<FiniteFieldElement_flint_fq>self).val, self._cparent)
        return x

    def __invert__(FiniteFieldElement_flint_fq self):
        """
        Return the multiplicative inverse of ``self``.

        EXAMPLE::

            sage: k.<a> = FiniteField(3^2, impl='flint_fq')
            sage: ~a
            a + 2
            sage: (a+1)*a
            2*a + 1
            sage: ~((2*a)/a)
            2
        """
        if fq_is_zero(self.val, self._cparent):
            raise ZeroDivisionError
        cdef FiniteFieldElement_flint_fq x = self._new()
        fq_inv(x.val, (<FiniteFieldElement_flint_fq>self).val, self._cparent)
        return x

    def __pow__(FiniteFieldElement_flint_fq self, object exp, object other):
        """
        Exponentiation.

        TESTS::

            sage: K.<a> = GF(5^10, impl='flint_fq')
            sage: n = (2*a)/a
            sage: n^-15
            2

        Large exponents are not a problem::

            sage: e = 3^10000
            sage: a^e
            2*a^9 + a^5 + 4*a^4 + 4*a^3 + a^2 + 3*a
            sage: a^(e % (5^10 - 1))
            2*a^9 + a^5 + 4*a^4 + 4*a^3 + a^2 + 3*a
        """
        cdef Integer exp_INT
        cdef fmpz_t exp_fmpz
        inv = 0
        if exp == 0:
            return self._parent.one_element()
        if exp < 0:
            if fq_is_zero(self.val, self._cparent):
                raise ZeroDivisionError
            inv = 1
            exp = -exp

        cdef FiniteFieldElement_flint_fq x = self._new()
        if isinstance(exp, int) or isinstance(exp, long):
            sig_on()
            fq_pow_ui(x.val, self.val, exp, self._cparent)
            sig_off()
        else:
            exp_INT = Integer(exp)
            sig_on()
            fmpz_init_set_readonly(exp_fmpz, exp_INT.value)
            fq_pow(x.val, self.val, exp_fmpz, self._cparent)
            fmpz_clear_readonly(exp_fmpz)
            sig_off()

        if inv:
            sig_on()
            fq_inv(x.val, x.val, self._cparent)
            sig_off()

        return x

    def polynomial(FiniteFieldElement_flint_fq self):
        """
        Return the unique representative of ``self`` as a polynomial
        over the prime field whose degree is less than the degree of
        the finite field over its prime field.

        EXAMPLES::

            sage: k.<a> = FiniteField(3^2, impl='flint_fq')
            sage: pol = a.polynomial()
            sage: pol
            a
            sage: parent(pol)
            Univariate Polynomial Ring in a over Finite Field of size 3

        ::

            sage: k = FiniteField(3^4, 'alpha', impl='flint_fq')
            sage: a = k.gen()
            sage: a.polynomial()
            alpha
            sage: (a**2 + 1).polynomial()
            alpha^2 + 1
            sage: (a**2 + 1).polynomial().parent()
            Univariate Polynomial Ring in alpha over Finite Field of size 3
        """
        cdef fmpz_t cflint
        cdef mpz_t cgmp
        cdef Integer cint
        cdef int i, n
        cdef FiniteField_flint_fq K
        cdef list clist
        fmpz_init(cflint)
        cint = Integer.__new__(Integer)
        K = self._parent
        n = K.degree()
        clist = []
        for i in xrange(n):
            fmpz_poly_get_coeff_fmpz(cflint, self.val, i)
            flint_mpz_init_set_readonly(cgmp, cflint)
            cint.set_from_mpz(cgmp)
            flint_mpz_clear_readonly(cgmp)
            clist.append(Integer(cint))
        fmpz_clear(cflint)
        return K.polynomial_ring()(clist)

    # No FLINT implementation yet
    #def charpoly(FiniteFieldElement_flint_fq self, object var='x'):
    #def is_square(FiniteFieldElement_flint_fq self):
    #def sqrt(FiniteFieldElement_flint_fq self, extend=False, all=False):
    #def multiplicative_order(FiniteFieldElement_flint_fq self):

    def norm(self):
        """
        Return the norm of self down to the prime subfield.

        This is the product of the Galois conjugates of self.

        EXAMPLES::

            sage: S.<b> = GF(5^2, impl="flint_fq"); S
            Finite Field in b of size 5^2
            sage: b.norm()
            2

        Next we consider a cubic extension::

            sage: S.<a> = GF(5^3, impl="flint_fq"); S
            Finite Field in a of size 5^3
            sage: a.norm()
            2
            sage: a * a^5 * (a^25)
            2
        """
        cdef fmpz_t t
        cdef mpz_t tgmp
        cdef Integer tint
        fmpz_init(t)
        tint = Integer.__new__(Integer)
        fq_norm(t, self.val, self._cparent)
        flint_mpz_init_set_readonly(tgmp, t)
        tint.set_from_mpz(tgmp)
        flint_mpz_clear_readonly(tgmp)
        fmpz_clear(t)
        return self.parent().prime_subfield()(tint)

    def trace(self):
        """
        Return the trace of this element, which is the sum of the
        Galois conjugates.

        EXAMPLES::

            sage: S.<a> = GF(5^3, impl="flint_fq"); S
            Finite Field in a of size 5^3
            sage: a.trace()
            0
            sage: a + a^5 + a^25
            0
            sage: z = a^2 + a + 1
            sage: z.trace()
            2
            sage: z + z^5 + z^25
            2
        """
        cdef fmpz_t t
        cdef mpz_t tgmp
        cdef Integer tint
        fmpz_init(t)
        tint = Integer.__new__(Integer)
        fq_trace(t, self.val, self._cparent)
        flint_mpz_init_set_readonly(tgmp, t)
        tint.set_from_mpz(tgmp)
        flint_mpz_clear_readonly(tgmp)
        fmpz_clear(t)
        return self.parent().prime_subfield()(tint)

    # JPF: the following should definitely go into element_base.pyx.
    def log(FiniteFieldElement_flint_fq self, FiniteFieldElement_flint_fq base):
        """
        Return `x` such that `b^x = a`, where `x`
        is `a` and `b` is the base.

        INPUT:

        -  ``b`` -- finite field element that generates the
           multiplicative group.

        OUTPUT:

        Integer `x` such that `a^x = b`, if it exists. Raises a
        ``ValueError`` exception if no such `x` exists.

        EXAMPLES::

            sage: F = GF(17)
            sage: F(3^11).log(F(3))
            11
            sage: F = GF(113)
            sage: F(3^19).log(F(3))
            19
            sage: F = GF(next_prime(10000))
            sage: F(23^997).log(F(23))
            997

        ::

            sage: F = FiniteField(2^10, 'a', impl='flint_fq')
            sage: g = F.gen()
            sage: b = g; a = g^37
            sage: a.log(b)
            37
            sage: b^37; a
            a^8 + a^7 + a^4 + a + 1
            a^8 + a^7 + a^4 + a + 1

        AUTHORS:

        - David Joyner and William Stein (2005-11)
        """
        from  sage.groups.generic import discrete_log

        b = self.parent()(base)
        # TODO: This function is TERRIBLE!
        return discrete_log(self, b)

    def lift(FiniteFieldElement_flint_fq self):
        """
        If ``self`` is an element of the prime field, return a lift of
        this element to an integer.

        EXAMPLE::

            sage: k = FiniteField(next_prime(10^10)^2, 'u', impl='flint_fq')
            sage: a = k(17)/k(19)
            sage: b = a.lift(); b
            7894736858
            sage: b.parent()
            Integer Ring
        """
        if fq_is_zero(self.val, self._cparent):
            return Integer(0)
        f = self.polynomial()
        if f.degree() == 0:
            return f.constant_coefficient().lift()
        else:
            raise ValueError("element is not in the prime field")

    def _integer_(self, ZZ=None):
        """
        Lift to a Sage integer, if possible.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='flint_fq')
            sage: b = k(2)
            sage: b._integer_()
            2
            sage: a._integer_()
            Traceback (most recent call last):
            ...
            ValueError: element is not in the prime field
        """
        return self.lift()

    def __int__(self):
        """
        Lift to a python int, if possible.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='flint_fq')
            sage: b = k(2)
            sage: int(b)
            2
            sage: int(a)
            Traceback (most recent call last):
            ...
            ValueError: element is not in the prime field
        """
        return int(self.lift())

    def __long__(self):
        """
        Lift to a python long, if possible.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='flint_fq')
            sage: b = k(2)
            sage: long(b)
            2L
        """
        return long(self.lift())

    def __float__(self):
        """
        Lift to a python float, if possible.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='flint_fq')
            sage: b = k(2)
            sage: float(b)
            2.0
        """
        return float(self.lift())

    def _pari_(self, var=None):
        """
        Return a PARI object representing ``self``.

        INPUT:

        - var -- ignored

        EXAMPLE::

            sage: k.<a> = FiniteField(3^3, impl='flint_fq')
            sage: b = a**2 + 2*a + 1
            sage: b._pari_()
            a^2 + 2*a + 1
        """
        return pari(self._pari_init_())

    def _pari_init_(self):
        """
        Return a string representing ``self`` in PARI.

        EXAMPLE::

            sage: k.<a> = GF(3^17, impl='flint_fq')
            sage: a._pari_init_()
            'subst(a+3*a,a,ffgen(Mod(1, 3)*x^17 + Mod(2, 3)*x + Mod(1, 3),a))'
            sage: k(1)._pari_init_()
            'subst(1+3*a,a,ffgen(Mod(1, 3)*x^17 + Mod(2, 3)*x + Mod(1, 3),a))'

        This is used for conversion to GP. The element is displayed
        as "a" but has correct arithmetic::

            sage: gp(a)
            a
            sage: gp(a).type()
            t_FFELT
            sage: gp(a)^100
            2*a^16 + 2*a^15 + a^4 + a + 1
            sage: gp(a^100)
            2*a^16 + 2*a^15 + a^4 + a + 1
            sage: gp(k(0))
            0
            sage: gp(k(0)).type()
            t_FFELT
        """
        ffgen = "ffgen(%s,%s)" % (self._parent.modulus()._pari_init_(), self._parent.variable_name())
        # Add this "zero" to ensure that the polynomial is not constant
        zero = "%s*%s" % (self._parent.characteristic(), self._parent.variable_name())
        return "subst(%s+%s,%s,%s)" % (self, zero, self._parent.variable_name(), ffgen)

    def _magma_init_(self, magma):
        """
        Return a string representing ``self`` in Magma.

        EXAMPLE::

            sage: GF(7)(3)._magma_init_(magma)            # optional - magma
            'GF(7)!3'
        """
        k = self._parent
        km = magma(k)
        return str(self).replace(k.variable_name(), km.gen(1).name())

    def _gap_init_(self):
        r"""
        Return the a string representing ``self`` in GAP.

        .. NOTE::

           The order of the parent field must be `\leq 65536`.  This
           function can be slow since elements of non-prime finite
           fields are represented in GAP as powers of a generator for
           the multiplicative group, so a discrete logarithm must be
           computed.

        EXAMPLE::

            sage: F = FiniteField(2^3, 'a', impl='flint_fq')
            sage: a = F.multiplicative_generator()
            sage: gap(a) # indirect doctest
            Z(2^3)
            sage: b = F.multiplicative_generator()
            sage: a = b^3
            sage: gap(a)
            Z(2^3)^3
            sage: gap(a^3)
            Z(2^3)^2

        You can specify the instance of the Gap interpreter that is used::

            sage: F = FiniteField(next_prime(200)^2, 'a', impl='flint_fq')
            sage: a = F.multiplicative_generator ()
            sage: a._gap_ (gap)
            Z(211^2)
            sage: (a^20)._gap_(gap)
            Z(211^2)^20

        Gap only supports relatively small finite fields::

            sage: F = FiniteField(next_prime(1000)^2, 'a', impl='flint_fq')
            sage: a = F.multiplicative_generator ()
            sage: gap._coerce_(a)
            Traceback (most recent call last):
            ...
            TypeError: order must be at most 65536
        """
        F = self._parent
        if F.order() > 65536:
            raise TypeError("order must be at most 65536")

        if self == 0:
            return '0*Z(%s)'%F.order()
        assert F.degree() > 1
        g = F.multiplicative_generator()
        n = self.log(g)
        return 'Z(%s)^%s'%(F.order(), n)


def unpickle_FiniteFieldElement_flint_fq(parent, elem):
    """
    EXAMPLE::

        sage: k.<a> = GF(2^20, impl='flint_fq')
        sage: e = k.random_element()
        sage: f = loads(dumps(e)) # indirect doctest
        sage: e == f
        True
    """
    return parent(elem)
