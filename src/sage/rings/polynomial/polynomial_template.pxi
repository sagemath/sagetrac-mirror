"""
Polynomial Template for C/C++ Library Interfaces
"""
#*****************************************************************************
#       Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#       Copyright (C) 2008 Robert Bradshaw
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.structure.element cimport ModuleElement, Element, RingElement
from sage.structure.element import coerce_binop, bin_op
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.integer cimport Integer
from sage.libs.all import pari_gen

import operator

from sage.interfaces.all import singular as singular_default

def make_element(parent, args):
    return parent(*args)

cdef inline Polynomial_template element_shift(self, int n):
     if not PY_TYPE_CHECK(self, Polynomial_template):
         if n > 0:
             error_msg = "Cannot shift %s << %n."%(self, n)
         else:
             error_msg = "Cannot shift %s >> %n."%(self, n)
         raise TypeError(error_msg)

     if n == 0:
         return self

     cdef celement *gen = celement_new((<Polynomial_template>self)._cparent)
     celement_gen(gen, 0, (<Polynomial_template>self)._cparent)
     celement_pow(gen, gen, abs(n), NULL, (<Polynomial_template>self)._cparent)
     cdef type T = self.__class__
     cdef Polynomial_template r = <Polynomial_template>T.__new__(T)
     celement_construct(&r.x, (<Polynomial_template>self)._cparent)
     r._parent = (<Polynomial_template>self)._parent
     r._cparent = (<Polynomial_template>self)._cparent

     if n > 0:
         celement_mul(&r.x, &(<Polynomial_template>self).x, gen, (<Polynomial_template>self)._cparent)
     else:
         celement_floordiv(&r.x, &(<Polynomial_template>self).x, gen, (<Polynomial_template>self)._cparent)

     celement_delete(gen, (<Polynomial_template>self)._cparent)
     return r

cdef class Polynomial_template(Polynomial):
    r"""
    Template for interfacing to external C / C++ libraries for implementations of polynomials.

    AUTHORS:

    - Robert Bradshaw (2008-10): original idea for templating
    - Martin Albrecht (2008-10): initial implementation

    This file implements a simple templating engine for linking univariate
    polynomials to their C/C++ library implementations. It requires a
    'linkage' file which implements the ``celement_`` functions (see
    :mod:`sage.libs.ntl.ntl_GF2X_linkage` for an example). Both parts are
    then plugged together by inclusion of the linkage file when inheriting from
    this class. See :mod:`sage.rings.polynomial.polynomial_gf2x` for an
    example.

    We illustrate the generic glueing using univariate polynomials over
    `\mathop{\mathrm{GF}}(2)`.

    .. note::

        Implementations using this template MUST implement coercion from base
        ring elements and ``__getitem__``. See
        :class:`~sage.rings.polynomial.polynomial_gf2x.Polynomial_GF2X` for an
        example.
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: P(0)
            0
            sage: P(GF(2)(1))
            1
            sage: P(3)
            1
            sage: P([1,0,1])
            x^2 + 1
            sage: P(map(GF(2),[1,0,1]))
            x^2 + 1
        """
        cdef celement *gen
        cdef celement *monomial
        cdef Py_ssize_t deg

        Polynomial.__init__(self, parent, is_gen=is_gen)

        (<Polynomial_template>self)._cparent = get_cparent(self._parent)

        if is_gen:
            celement_construct(&self.x, (<Polynomial_template>self)._cparent)
            celement_gen(&self.x, 0, (<Polynomial_template>self)._cparent)

        elif PY_TYPE_CHECK(x, Polynomial_template):
            try:
                celement_construct(&self.x, (<Polynomial_template>self)._cparent)
                celement_set(&self.x, &(<Polynomial_template>x).x, (<Polynomial_template>self)._cparent)
            except NotImplementedError:
                raise TypeError("%s not understood."%x)

        elif PY_TYPE_CHECK(x, int) or PY_TYPE_CHECK(x, Integer):
            try:
                celement_construct(&self.x, (<Polynomial_template>self)._cparent)
                celement_set_si(&self.x, int(x), (<Polynomial_template>self)._cparent)
            except NotImplementedError:
                raise TypeError("%s not understood."%x)

        elif PY_TYPE_CHECK(x, list) or PY_TYPE_CHECK(x, tuple):
            celement_construct(&self.x, (<Polynomial_template>self)._cparent)
            gen = celement_new((<Polynomial_template>self)._cparent)
            monomial = celement_new((<Polynomial_template>self)._cparent)

            celement_set_si(&self.x, 0, (<Polynomial_template>self)._cparent)
            celement_gen(gen, 0, (<Polynomial_template>self)._cparent)

            deg = 0
            for e in x:
                # r += parent(e)*power
                celement_pow(monomial, gen, deg, NULL, (<Polynomial_template>self)._cparent)
                celement_mul(monomial, &(<Polynomial_template>self.__class__(parent, e)).x, monomial, (<Polynomial_template>self)._cparent)
                celement_add(&self.x, &self.x, monomial, (<Polynomial_template>self)._cparent)
                deg += 1

            celement_delete(gen, (<Polynomial_template>self)._cparent)
            celement_delete(monomial, (<Polynomial_template>self)._cparent)

        elif PY_TYPE_CHECK(x, dict):
            celement_construct(&self.x, (<Polynomial_template>self)._cparent)
            gen = celement_new((<Polynomial_template>self)._cparent)
            monomial = celement_new((<Polynomial_template>self)._cparent)

            celement_set_si(&self.x, 0, (<Polynomial_template>self)._cparent)
            celement_gen(gen, 0, (<Polynomial_template>self)._cparent)

            for deg, coef in x.iteritems():
                celement_pow(monomial, gen, deg, NULL, (<Polynomial_template>self)._cparent)
                celement_mul(monomial, &(<Polynomial_template>self.__class__(parent, coef)).x, monomial, (<Polynomial_template>self)._cparent)
                celement_add(&self.x, &self.x, monomial, (<Polynomial_template>self)._cparent)

            celement_delete(gen, (<Polynomial_template>self)._cparent)
            celement_delete(monomial, (<Polynomial_template>self)._cparent)

        elif PY_TYPE_CHECK(x, pari_gen):
            k = (<Polynomial_template>self)._parent.base_ring()
            x = [k(w) for w in x.list()]
            self.__class__.__init__(self, parent, x, check=True, is_gen=False, construct=construct)
        elif PY_TYPE_CHECK(x, Polynomial):
            k = (<Polynomial_template>self)._parent.base_ring()
            x = [k(w) for w in list(x)]
            Polynomial_template.__init__(self, parent, x, check=True, is_gen=False, construct=construct)
        elif PY_TYPE_CHECK(x, FractionFieldElement) and (x.parent().base() is parent or x.parent().base() == parent) and x.denominator() == 1:
            x = x.numerator()
            self.__class__.__init__(self, parent, x, check=check, is_gen=is_gen, construct=construct)
        else:
            x = parent.base_ring()(x)
            self.__class__.__init__(self, parent, x, check=check, is_gen=is_gen, construct=construct)

    def get_cparent(self):
        return <long> self._cparent

    def __reduce__(self):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: loads(dumps(x)) == x
            True
        """
        return make_element, ((<Polynomial_template>self)._parent, (self.list(), False, self.is_gen()))

    def list(self):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: x.list()
            [0, 1]
            sage: list(x)
            [0, 1]
        """
        cdef Py_ssize_t i
        return [self[i] for i in range(celement_len(&self.x, (<Polynomial_template>self)._cparent))]

    def __dealloc__(self):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: del x

        TEST:

        The following has been a problem in a preliminary version of
        :trac:`12313`::

            sage: K.<z> = GF(4)
            sage: P.<x> = K[]
            sage: del P
            sage: del x
            sage: import gc
            sage: _ = gc.collect()
        """
        celement_destruct(&self.x, (<Polynomial_template>self)._cparent)

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: x + 1
            x + 1
        """
        cdef type T = self.__class__
        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)

        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        celement_add(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, (<Polynomial_template>self)._cparent)
        #assert(r._parent(pari(self) + pari(right)) == r)
        return r

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: x - 1
            x + 1
        """
        cdef type T = self.__class__
        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        celement_sub(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, (<Polynomial_template>self)._cparent)
        #assert(r._parent(pari(self) - pari(right)) == r)
        return r

    def __neg__(self):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: -x
            x
        """
        cdef type T = self.__class__
        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        celement_neg(&r.x, &self.x, (<Polynomial_template>self)._cparent)
        #assert(r._parent(-pari(self)) == r)
        return r

    cpdef ModuleElement _rmul_(self, RingElement left):
        """
        EXAMPLES::

            sage: P.<x> = GF(2)[]
            sage: t = x^2 + x + 1
            sage: 0*t
            0
            sage: 1*t
            x^2 + x + 1

            sage: R.<y> = GF(5)[]
            sage: u = y^2 + y + 1
            sage: 3*u
            3*y^2 + 3*y + 3
            sage: 5*u
            0
            sage: (2^81)*u
            2*y^2 + 2*y + 2
            sage: (-2^81)*u
            3*y^2 + 3*y + 3
        """
        cdef type T = self.__class__
        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        celement_mul_scalar(&r.x, &(<Polynomial_template>self).x, left, (<Polynomial_template>self)._cparent)
        return r

    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        EXAMPLES::

            sage: P.<x> = GF(2)[]
            sage: t = x^2 + x + 1
            sage: t*0
            0
            sage: t*1
            x^2 + x + 1

            sage: R.<y> = GF(5)[]
            sage: u = y^2 + y + 1
            sage: u*3
            3*y^2 + 3*y + 3
            sage: u*5
            0
        """
        # all currently implemented rings are commutative
        return self._rmul_(right)

    cpdef RingElement _mul_(self, RingElement right):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: x*(x+1)
            x^2 + x
        """
        cdef type T = self.__class__
        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        celement_mul(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, (<Polynomial_template>self)._cparent)
        #assert(r._parent(pari(self) * pari(right)) == r)
        return r

    @coerce_binop
    def gcd(self, Polynomial_template other):
        """
        Return the greatest common divisor of self and other.

        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: f = x*(x+1)
            sage: f.gcd(x+1)
            x + 1
            sage: f.gcd(x^2)
            x
        """
        if(celement_is_zero(&self.x, (<Polynomial_template>self)._cparent)):
            return other
        if(celement_is_zero(&other.x, (<Polynomial_template>self)._cparent)):
            return self

        cdef type T = self.__class__
        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        celement_gcd(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>other).x, (<Polynomial_template>self)._cparent)
        #assert(r._parent(pari(self).gcd(pari(other))) == r)
        return r

    @coerce_binop
    def xgcd(self, Polynomial_template other):
        """
        Computes extended gcd of self and other.

        EXAMPLE::

            sage: P.<x> = GF(7)[]
            sage: f = x*(x+1)
            sage: f.xgcd(x+1)
            (x + 1, 0, 1)
            sage: f.xgcd(x^2)
            (x, 1, 6)
        """
        if(celement_is_zero(&self.x, (<Polynomial_template>self)._cparent)):
            return other, self._parent(0), self._parent(1)
        if(celement_is_zero(&other.x, (<Polynomial_template>self)._cparent)):
            return self, self._parent(1), self._parent(0)

        cdef type T = self.__class__
        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent

        cdef Polynomial_template s = <Polynomial_template>T.__new__(T)
        celement_construct(&s.x, (<Polynomial_template>self)._cparent)
        s._parent = (<Polynomial_template>self)._parent
        s._cparent = (<Polynomial_template>self)._cparent

        cdef Polynomial_template t = <Polynomial_template>T.__new__(T)
        celement_construct(&t.x, (<Polynomial_template>self)._cparent)
        t._parent = (<Polynomial_template>self)._parent
        t._cparent = (<Polynomial_template>self)._cparent

        celement_xgcd(&r.x, &s.x, &t.x, &(<Polynomial_template>self).x, &(<Polynomial_template>other).x, (<Polynomial_template>self)._cparent)
        #rp, sp, tp = pari(self).xgcd(pari(other))
        #assert(r._parent(rp) == r)
        #assert(s._parent(sp) == s)
        #assert(t._parent(tp) == t)
        return r,s,t

    def __floordiv__(self, right):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: x//(x + 1)
            1
            sage: (x + 1)//x
            1
        """
        # We can't use @coerce_binop for operators in cython classes,
        # so we use sage.structure.element.bin_op to handle coercion.
        if type(self) is not type(right) or \
                (<Polynomial_template>self)._parent is not (<Polynomial_template>right)._parent:
            return bin_op(self, right, operator.mod)
        cdef Polynomial_template _right = <Polynomial_template>right

        if celement_is_zero(&_right.x, (<Polynomial_template>self)._cparent):
            raise ZeroDivisionError
        cdef type T = self.__class__
        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        #assert(r._parent(pari(self) // pari(right)) == r)
        celement_floordiv(&r.x, &(<Polynomial_template>self).x, &(<Polynomial_template>right).x, (<Polynomial_template>self)._cparent)
        return r

    def __mod__(self, other):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: (x^2 + 1) % x^2
            1


        TESTS::

        We test that #10578 is fixed::

            sage: P.<x> = GF(2)[]
            sage: x % 1r
            0

        """
        # We can't use @coerce_binop for operators in cython classes,
        # so we use sage.structure.element.bin_op to handle coercion.
        if type(self) is not type(other) or \
                (<Polynomial_template>self)._parent is not (<Polynomial_template>other)._parent:
            return bin_op(self, other, operator.mod)
        cdef Polynomial_template _other = <Polynomial_template>other

        if celement_is_zero(&_other.x, (<Polynomial_template>self)._cparent):
            raise ZeroDivisionError

        cdef type T = self.__class__
        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        celement_mod(&r.x, &(<Polynomial_template>self).x, &_other.x, (<Polynomial_template>self)._cparent)
        #assert(r._parent(pari(self) % pari(other)) == r)
        return r

    @coerce_binop
    def quo_rem(self, Polynomial_template right):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: f = x^2 + x + 1
            sage: f.quo_rem(x + 1)
            (x, 1)
        """
        if celement_is_zero(&right.x, (<Polynomial_template>self)._cparent):
            raise ZeroDivisionError

        cdef type T = self.__class__
        cdef Polynomial_template q = <Polynomial_template>T.__new__(T)
        celement_construct(&q.x, (<Polynomial_template>self)._cparent)
        q._parent = (<Polynomial_template>self)._parent
        q._cparent = (<Polynomial_template>self)._cparent

        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent

        celement_quorem(&q.x, &r.x, &(<Polynomial_template>self).x, &right.x, (<Polynomial_template>self)._cparent)
        return q,r

    def __long__(self):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: int(x)
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce polynomial with degree 1 to integer.

            sage: int(P(1))
            1
        """
        if celement_len(&self.x, (<Polynomial_template>self)._cparent) > 1:
            raise ValueError("Cannot coerce polynomial with degree %d to integer."%(self.degree()))
        return int(self[0])

    def __nonzero__(self):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: bool(x), x.is_zero()
            (True, False)
            sage: bool(P(0)), P(0).is_zero()
            (False, True)
        """
        return not celement_is_zero(&self.x, (<Polynomial_template>self)._cparent)

    def __richcmp__(left, right, int op):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: x != 1
            True
            sage: x < 1
            False
            sage: x > 1
            True
        """
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
        """
        return celement_cmp(&(<Polynomial_template>left).x, &(<Polynomial_template>right).x, (<Polynomial_template>left)._cparent)

    def __hash__(self):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: {x:1}
            {x: 1}
        """
        cdef long result = 0 # store it in a c-int and just let the overflowing additions wrap
        cdef long result_mon
        cdef long c_hash
        cdef long var_name_hash
        cdef int i
        for i from 0<= i <= self.degree():
            if i == 1:
                # we delay the hashing until now to not waste it one a constant poly
                var_name_hash = hash(self.variable_name())
            #  I'm assuming (incorrectly) that hashes of zero indicate that the element is 0.
            # This assumption is not true, but I think it is true enough for the purposes and it
            # it allows us to write fast code that omits terms with 0 coefficients.  This is
            # important if we want to maintain the '==' relationship with sparse polys.
            c_hash = hash(self[i])
            if c_hash != 0:
                if i == 0:
                    result += c_hash
                else:
                    # Hash (self[i], generator, i) as a tuple according to the algorithm.
                    result_mon = c_hash
                    result_mon = (1000003 * result_mon) ^ var_name_hash
                    result_mon = (1000003 * result_mon) ^ i
                    result += result_mon
        if result == -1:
            return -2
        return result


    def __pow__(self, ee, modulus):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: P.<x> = GF(2)[]
            sage: x^1000
            x^1000
            sage: (x+1)^2
            x^2 + 1
            sage: (x+1)^(-2)
            1/(x^2 + 1)
            sage: f = x^9 + x^7 + x^6 + x^5 + x^4 + x^2 + x
            sage: h = x^10 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + 1
            sage: (f^2) % h
            x^9 + x^8 + x^7 + x^5 + x^3
            sage: pow(f, 2, h)
            x^9 + x^8 + x^7 + x^5 + x^3
        """
        if not PY_TYPE_CHECK(self, Polynomial_template):
            raise NotImplementedError("%s^%s not defined."%(ee,self))
        cdef bint recip = 0, do_sig

        cdef long e
        try:
            e = ee
        except OverflowError:
            return Polynomial.__pow__(self, ee, modulus)
        if e != ee:
            raise TypeError("Only integral powers defined.")
        elif e < 0:
            recip = 1 # delay because powering frac field elements is slow
            e = -e
        if not self:
            if e == 0:
                raise ArithmeticError("0^0 is undefined.")
        cdef type T = self.__class__
        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)

        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        parent = (<Polynomial_template>self)._parent
        r._parent = parent
        r._cparent = (<Polynomial_template>self)._cparent

        if modulus is None:
            celement_pow(&r.x, &(<Polynomial_template>self).x, e, NULL, (<Polynomial_template>self)._cparent)
        else:
            if parent is not (<Polynomial_template>modulus)._parent and parent != (<Polynomial_template>modulus)._parent:
                modulus = parent._coerce_(modulus)
            celement_pow(&r.x, &(<Polynomial_template>self).x, e, &(<Polynomial_template>modulus).x, (<Polynomial_template>self)._cparent)

        #assert(r._parent(pari(self)**ee) == r)
        if recip:
            return ~r
        else:
            return r

    def __copy__(self):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: copy(x) is x
            False
            sage: copy(x) == x
            True
        """
        cdef type T = self.__class__
        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        celement_set(&r.x, &self.x, (<Polynomial_template>self)._cparent)
        return r

    def is_gen(self):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: x.is_gen()
            True
            sage: (x+1).is_gen()
            False
        """
        cdef celement *gen = celement_new((<Polynomial_template>self)._cparent)
        celement_gen(gen, 0, (<Polynomial_template>self)._cparent)
        cdef bint r = celement_equal(&self.x, gen, (<Polynomial_template>self)._cparent)
        celement_delete(gen, (<Polynomial_template>self)._cparent)
        return r

    def shift(self, int n):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: f = x^3 + x^2 + 1
            sage: f.shift(1)
            x^4 + x^3 + x
            sage: f.shift(-1)
            x^2 + x
        """
        return element_shift(self, n)

    def __lshift__(self, int n):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: f = x^3 + x^2 + 1
            sage: f << 1
            x^4 + x^3 + x
            sage: f << -1
            x^2 + x
        """
        return element_shift(self, n)

    def __rshift__(self, int n):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: x>>1
            1
            sage: (x^2 + x)>>1
            x + 1
            sage: (x^2 + x) >> -1
            x^3 + x^2
        """
        return element_shift(self, -n)

    def is_zero(self):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: x.is_zero()
            False
        """
        return celement_is_zero(&self.x, (<Polynomial_template>self)._cparent)

    def is_one(self):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: P(1).is_one()
            True
        """
        return celement_is_one(&self.x, (<Polynomial_template>self)._cparent)

    def degree(self):
        """
        EXAMPLE::

            sage: P.<x> = GF(2)[]
            sage: x.degree()
            1
            sage: P(1).degree()
            0
            sage: P(0).degree()
            -1
        """
        return Integer(celement_len(&self.x, (<Polynomial_template>self)._cparent)-1)

    cpdef Polynomial truncate(self, long n):
        r"""
        Returns this polynomial mod `x^n`.

        EXAMPLES::

            sage: R.<x> =GF(2)[]
            sage: f = sum(x^n for n in range(10)); f
            x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
            sage: f.truncate(6)
            x^5 + x^4 + x^3 + x^2 + x + 1
        """
        cdef type T = self.__class__
        cdef Polynomial_template r = <Polynomial_template>T.__new__(T)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent

        if n <= 0:
            return r

        cdef celement *gen = celement_new((<Polynomial_template>self)._cparent)
        celement_gen(gen, 0, (<Polynomial_template>self)._cparent)
        celement_pow(gen, gen, n, NULL, (<Polynomial_template>self)._cparent)

        celement_mod(&r.x, &self.x, gen, (<Polynomial_template>self)._cparent)
        celement_delete(gen, (<Polynomial_template>self)._cparent)
        return r

    def _singular_(self, singular=singular_default, have_ring=False):
        r"""
        Return Singular representation of this polynomial

        INPUT:

        - ``singular`` -- Singular interpreter (default: default interpreter)
        - ``have_ring`` -- set to True if the ring was already set in Singular

        EXAMPLE::

            sage: P.<x> = PolynomialRing(GF(7))
            sage: f = 3*x^2 + 2*x + 5
            sage: singular(f)
            3*x^2+2*x-2
        """
        if not have_ring:
            self.parent()._singular_(singular).set_ring() #this is expensive
        return singular(self._singular_init_())

    def _derivative(self, var=None):
        r"""
        Returns the formal derivative of self with respect to var.

        var must be either the generator of the polynomial ring to which
        this polynomial belongs, or None (either way the behaviour is the
        same).

        .. seealso:: :meth:`.derivative`

        EXAMPLES::

            sage: R.<x> = Integers(77)[]
            sage: f = x^4 - x - 1
            sage: f._derivative()
            4*x^3 + 76
            sage: f._derivative(None)
            4*x^3 + 76

            sage: f._derivative(2*x)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to 2*x

            sage: y = var("y")
            sage: f._derivative(y)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to y
        """
        if var is not None and var is not self._parent.gen():
            raise ValueError("cannot differentiate with respect to %s" % var)

        P = self.parent()
        x = P.gen()
        res = P(0)
        for i,c in enumerate(self.list()[1:]):
            res += (i+1)*c*x**i
        return res
