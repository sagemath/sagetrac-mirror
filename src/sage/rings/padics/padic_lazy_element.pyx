# ****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************


from sage.libs.flint.types cimport *
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *

cdef extern from "sage/rings/padics/padic_lazy_element_helper.c":
    cdef void flint_randinit(flint_rand_t state)
    cdef void flint_randclear(flint_rand_t state)
    cdef fmpz* get_coeff(fmpz_poly_t poly, slong i)
    cdef void get_slice(fmpz_poly_t slice, fmpz_poly_t poly, slong start, slong length)
    cdef void iadd_coeff(fmpz_poly_t poly, fmpz_t summand, slong i)
    cdef void isub_coeff(fmpz_poly_t poly, fmpz_t summand, slong i)
    cdef void iadd_shifted(fmpz_poly_t poly, fmpz_poly_t summand, slong shift)
    cdef void reduce_coeff(fmpz_poly_t poly, slong i, fmpz_t modulus)

from sage.ext.stdsage cimport PY_NEW
from sage.structure.element import coerce_binop

from sage.rings.all import ZZ
from sage.rings.integer cimport Integer
from sage.rings.infinity import Infinity

from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.padics.precision_error import PrecisionError

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1
MAXORDP = ZZ(maxordp)


cpdef lazy_sum(parent, summands):
    n = len(summands)
    signs = n * [True]
    return pAdicLazyElement_add(parent, summands, signs)


cdef class pAdicLazyElement(pAdicGenericElement):
    def __cinit__(self):
        fmpz_poly_init(self._digits)

    def __dealloc__(self):
        fmpz_poly_clear(self._digits)

    def __init__(self, parent):
        pAdicGenericElement.__init__(self, parent)
        self._valuation = 0
        self._precrel = 0
        # should find something better (prime_pow?)
        cdef Integer p = self._parent.prime()
        fmpz_set_mpz(self._prime, p.value)

    def prime(self):
        return self._parent.prime()

    cdef int _next_c(self):
        raise NotImplementedError("must be implemented in subclasses")

    cdef int _jump_c(self, slong prec):
        cdef int error
        cdef slong i, stop
        stop = min(prec, maxordp) - self._precrel - self._valuation
        for i in range(stop):
            error = self._next_c()
            if error:
                return error
        return prec >= maxordp

    def jump(self, prec, quiet=True):
        if prec not in ZZ:
            raise ValueError("precision must be an integer")
        prec = ZZ(prec)
        if prec >= maxordp:
            raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
        error = self._jump_c(long(prec))
        if not quiet and error == 1:
            raise PrecisionError
        if error == 2:
            raise RecursionError("definition looks circular")

    def expansion(self, n=None):
        if n is None:
            return ExpansionIter(self, self._valuation, self._valuation + self._precrel)
        if isinstance(n, slice):
            if n.step is not None:
                raise NotImplementedError
            if n.start is None:
                start = 0
            elif n.start in ZZ:
                start = ZZ(n.start)
            else:
                raise ValueError
            if n.stop is None:
                stop = None
            elif n.stop in ZZ:
                start = ZZ(n.stop)
            else:
                raise ValueError
            return ExpansionIter(self, start, stop)
        else:
            if n not in ZZ:
                raise IndexError("index must be an integer")
            n = ZZ(n)
            if n >= maxordp:
                raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
            return self._digit(long(n))

    cdef Integer _digit(self, slong i):
        self.jump(i+1, quiet=False)
        cdef Integer ans = PY_NEW(Integer)
        cdef fmpz* coeff = get_coeff(self._digits, i - self._valuation)
        fmpz_get_mpz(ans.value, coeff)
        return ans

    def _repr_(self):
        self.jump(self._parent.default_prec())
        return pAdicGenericElement._repr_(self)

    cdef bint _is_equal(self, pAdicLazyElement right, slong prec):
        cdef slong i
        self._jump_c(prec)
        right._jump_c(prec)
        if self._valuation >= prec:
            return right._valuation >= prec
        if right._valuation >= prec:
            return False
        if self._valuation != right._valuation:
            return False
        for i in range(prec - self._valuation):
            if not fmpz_equal(get_coeff(self._digits, i), get_coeff(right._digits, i)):
                return False
        return True

    @coerce_binop
    def is_equal_at_precision(self, other, prec):
        return self._is_equal(other, long(prec))

    def __eq__(self, other):
        prec = min(self.precision_absolute(), other.precision_absolute())
        prec = max(prec, self._parent.default_prec())
        return self.is_equal_at_precision(other, prec)

    cpdef bint _is_exact_zero(self) except -1:
        return isinstance(self, pAdicLazyElement_zero)

    cpdef bint _is_inexact_zero(self) except -1:
        return self._precrel == 0

    def is_exact(self):
        return self._is_exact_zero()

    def precision_absolute(self):
        if self._is_exact_zero():
            return Infinity
        return ZZ(self._precrel + self._valuation)

    def precision_relative(self):
        return ZZ(self._precrel)

    def valuation(self):
        if self._is_exact_zero():
            return Infinity
        return ZZ(self._valuation)

    def unit_part(self):
        if self._precrel == 0:
            raise PrecisionError
        return self >> self._valuation

    def __rshift__(self, s):
        cdef slong shift = long(s)
        if shift:
            return pAdicLazyElement_shift((<pAdicLazyElement>self)._parent, self, shift, not (<pAdicLazyElement>self)._parent.is_field())
        else:
            return self

    def __lshift__(self, s):
        return self.__rshift__(-s)

    cpdef _add_(self, other):
        cdef list summands
        cdef list signs
        if isinstance(self, pAdicLazyElement_zero):
            return other
        if isinstance(other, pAdicLazyElement_zero):
            return self
        if isinstance(self, pAdicLazyElement_add):
            summands = list((<pAdicLazyElement_add>self)._summands)
            signs = list((<pAdicLazyElement_add>self)._signs)
        else:
            summands = [self]
            signs = [True]
        if isinstance(other, pAdicLazyElement_add):
            summands.extend((<pAdicLazyElement_add>other)._summands)
            signs.extend((<pAdicLazyElement_add>other)._signs)
        else:
            summands.append(other)
            signs.append(True)
        return pAdicLazyElement_add(self._parent, summands, signs)

    cpdef _sub_(self, other):
        cdef list summands
        cdef list signs
        if isinstance(self, pAdicLazyElement_zero):
            return -other
        if isinstance(other, pAdicLazyElement_zero):
            return self
        if isinstance(self, pAdicLazyElement_add):
            summands = list((<pAdicLazyElement_add>self)._summands)
            signs = list((<pAdicLazyElement_add>self)._signs)
        else:
            summands = [self]
            signs = [True]
        if isinstance(other, pAdicLazyElement_add):
            summands.extend((<pAdicLazyElement_add>other)._summands)
            signs.extend([ not sign for sign in (<pAdicLazyElement_add>other)._signs ])
        else:
            summands.append(other)
            signs.append(False)
        return pAdicLazyElement_add(self._parent, summands, signs)

    cpdef _neg_(self):
        cdef list summands
        cdef list signs
        if isinstance(self, pAdicLazyElement_add):
            summands = list((<pAdicLazyElement_add>self)._summands)
            signs = [ not sign for sign in (<pAdicLazyElement_add>self)._signs ]
        else:
            summands = [self]
            signs = [False]
        return pAdicLazyElement_add(self._parent, summands, signs)

    cpdef _mul_(self, other):
        if isinstance(self, pAdicLazyElement_zero) or isinstance(other, pAdicLazyElement_one):
            return self
        if isinstance(self, pAdicLazyElement_one) or isinstance(other, pAdicLazyElement_zero):
            return other
        return pAdicLazyElement_mul(self._parent, self, <pAdicLazyElement>other)

    cpdef _div_(self, other):
        raise NotImplementedError("division is not implemented")


# Assignations
##############

# Zero

cdef class pAdicLazyElement_zero(pAdicLazyElement):
    def __init__(self, parent):
        pAdicLazyElement.__init__(self, parent)
        self._valuation = maxordp

    cdef int _jump_c(self, slong prec):
        return 0

    cdef int _next_c(self):
        return 0


# One

cdef class pAdicLazyElement_one(pAdicLazyElement):
    def __init__(self, parent):
        pAdicLazyElement.__init__(self, parent)
        fmpz_poly_set_ui(self._digits, 1)

    cdef int _jump_c(self, slong prec):
        if self._precrel < prec:
            self._precrel = prec
        return 0

    cdef int _next_c(self):
        self._precrel += 1
        return 0


# Value

cdef class pAdicLazyElement_value(pAdicLazyElement):
    def __init__(self, parent, value, prec=None, maxprec=None):
        pAdicLazyElement.__init__(self, parent)
        cdef Integer x = ZZ(value)
        fmpz_poly_set_mpz(self._digits, x.value)
        self._maxprec = maxprec

    cdef int _next_c(self):
        if (self._maxprec is not None) and (self._valuation + self._precrel >= self._maxprec):
            return 1
        reduce_coeff(self._digits, self._precrel, self._prime)
        if self._precrel == 0 and fmpz_is_zero(get_coeff(self._digits, 0)):
            self._valuation += 1
            fmpz_poly_shift_right(self._digits, self._digits, 1)
        else:
            self._precrel += 1
        return 0


# Random

cdef class pAdicLazyElement_random(pAdicLazyElement):
    def __cinit__(self):
        flint_randinit(self._randstate)

    def __dealloc__(self):
        flint_randclear(self._randstate)

    def __init__(self, parent):
        pAdicLazyElement.__init__(self, parent)

    cdef int _next_c(self):
        cdef fmpz_t r;
        fmpz_randm(r, self._randstate, self._prime)
        if self._precrel == 0 and fmpz_is_zero(r):
            self._valuation += 1
        else:
            fmpz_poly_set_coeff_fmpz(self._digits, self._precrel, r);
            self._precrel += 1
        return 0


# Operations
############

# Shift

cdef class pAdicLazyElement_shift(pAdicLazyElement):
    def __init__(self, parent, pAdicLazyElement x, slong s, bint truncate):
        cdef slong start = 0
        cdef fmpz_poly_t digits;
        pAdicLazyElement.__init__(self, parent)
        self._x = x
        self._shift = s
        if truncate and x._valuation < s:
            start = s - x._valuation
            self._valuation = 0
            self._precrel = x._precrel - start
        else:
            self._valuation = x._valuation - s
            self._precrel = x._precrel
        if self._precrel < 0:
            self._precrel = 0
        else:
            get_slice(digits, x._digits, start, self._precrel)
            fmpz_poly_set(self._digits, digits)

    cdef int _next_c(self):
        cdef slong n = self._valuation + self._precrel
        cdef pAdicLazyElement x = self._x
        cdef fmpz* digit
        cdef int error

        error = x._jump_c(n + self._shift + 1)
        if error:
            return error
        digit = get_coeff(x._digits, n + self._shift - x._valuation)
        fmpz_poly_set_coeff_fmpz(self._digits, self._precrel, digit)
        if self._precrel == 0 and fmpz_is_zero(get_coeff(self._digits, 0)):
            self._valuation += 1
            fmpz_poly_shift_right(self._digits, self._digits, 1)
        else:
            self._precrel += 1
        return 0


# Addition

cdef class pAdicLazyElement_add(pAdicLazyElement):
    def __init__(self, parent, summands, signs):
        pAdicLazyElement.__init__(self, parent)
        self._valuation = min((<pAdicLazyElement>summand)._valuation for summand in summands)
        self._summands = summands
        self._signs = signs

    cdef int _next_c(self):
        cdef pAdicLazyElement summand
        cdef slong n = self._valuation + self._precrel
        cdef fmpz* coeff
        cdef int error

        for summand in self._summands:
            error = summand._jump_c(n+1)
            if error:
                return error
        cdef slong i
        for i in range(len(self._summands)):
            summand = self._summands[i]
            coeff = get_coeff(summand._digits, n - summand._valuation)
            if self._signs[i]:
                iadd_coeff(self._digits, coeff, n - self._valuation)
            else:
                isub_coeff(self._digits, coeff, n - self._valuation)
        reduce_coeff(self._digits, self._precrel, self._prime)
        if self._precrel == 0 and fmpz_is_zero(get_coeff(self._digits, 0)):
            self._valuation += 1
            fmpz_poly_shift_right(self._digits, self._digits, 1)
        else:
            self._precrel += 1
        return 0


# Multiplication

cdef class pAdicLazyElement_mul(pAdicLazyElement):
    def __cinit__(self):
        fmpz_init(self._tmp_coeff)
        fmpz_poly_init(self._tmp_poly)

    def __dealloc__(self):
        fmpz_clear(self._tmp_coeff)
        fmpz_poly_clear(self._tmp_poly)

    def __init__(self, parent, x, y):
        pAdicLazyElement.__init__(self, parent)
        self._x = <pAdicLazyElement?>x
        self._y = <pAdicLazyElement?>y
        self._valuation = self._x._valuation + self._y._valuation

    cdef int _next_c(self):
        cdef pAdicLazyElement x = self._x
        cdef pAdicLazyElement y = self._y
        cdef slong n = self._valuation + self._precrel

        cdef int errorx = x._jump_c(n - y._valuation + 1)
        cdef int errory = y._jump_c(n - x._valuation + 1)
        if self._precrel == 0:
            self._valuation = x._valuation + y._valuation
            if self._valuation > n:
                return 0
            if self._valuation < n or x._precrel == 0 or y._precrel == 0:
                if errorx and errory:
                    return min(errorx, errory)
                else:
                    return max(errorx, errory)
        elif errorx or errory:
            return max(errorx, errory)

        n = self._precrel
        fmpz_mul(self._tmp_coeff, get_coeff(x._digits, 0), get_coeff(y._digits, n))
        iadd_coeff(self._digits, self._tmp_coeff, n)
        if n:
            fmpz_mul(self._tmp_coeff, get_coeff(x._digits, n), get_coeff(y._digits, 0))
            iadd_coeff(self._digits, self._tmp_coeff, n)

        cdef slong m = n + 2
        cdef slong len = 1
        cdef fmpz_poly_t slicex, slicey
        while (m & 1 == 0) and (m > 3):
            m >>= 1
            len <<= 1
            get_slice(slicex, x._digits, len - 1, len)
            get_slice(slicey, y._digits, (m-1)*len - 1, len)
            fmpz_poly_mul(self._tmp_poly, slicex, slicey)
            iadd_shifted(self._digits, self._tmp_poly, n)
            if m > 2:
                get_slice(slicex, x._digits, (m-1)*len - 1, len)
                get_slice(slicey, y._digits, len - 1, len)
                fmpz_poly_mul(self._tmp_poly, slicex, slicey)
                iadd_shifted(self._digits, self._tmp_poly, n)

        reduce_coeff(self._digits, self._precrel, self._prime)
        self._precrel += 1
        return 0


#cdef class pAdicLazyElement_lmul(pAdicLazyElement):


# Division

#cdef class pAdicLazyElement_div(pAdicLazyElement):


# Square root

#cdef class pAdicLazyElement_sqrt(pAdicLazyElement):



# Self-referent definitions
###########################

cdef class pAdicLazyElement_selfref(pAdicLazyElement):
    def __init__(self, parent, valuation):
        pAdicLazyElement.__init__(self, parent)
        self._valuation = valuation
        self._definition = None
        self._next = False

    cpdef set(self, pAdicLazyElement definition):
        if self._definition is not None:
            raise ValueError("this self-referent number is already defined")
        self._definition = <pAdicLazyElement>definition
        cdef slong valsve = self._valuation
        cdef int error = self._jump_c(self._valuation + self._parent.default_prec())
        if error == 2:
            self._definition = None
            self._valuation = valsve
            self._precrel = 0
            self._next = False
            raise RecursionError("definition looks circular")

    def __eq__(self, other):
        if self._definition is None:
            self.set(other)
        return pAdicLazyElement.__eq__(self, other)

    cdef int _next_c(self):
        cdef pAdicLazyElement definition = self._definition
        cdef fmpz* digit
        cdef slong diffval
        cdef int error

        if definition is None:
            return 1
        if self._next:
            return 2

        self._next = True
        error = definition._jump_c(self._valuation + self._precrel + 1)
        if not error:
            diffval = self._valuation - definition._valuation
            if diffval < 0:
                self._valuation = definition._valuation
            else:
                digit = get_coeff(definition._digits, self._precrel + diffval)
                fmpz_poly_set_coeff_fmpz(self._digits, self._precrel, digit);
                if self._precrel == 0 and fmpz_is_zero(digit):
                    self._valuation += 1
                    fmpz_poly_shift_right(self._digits, self._digits, 1)
                else:
                    self._precrel += 1
        self._next = False
        return error


# Expansion
###########

cdef class ExpansionIter(object):
    cdef pAdicLazyElement elt
    cdef slong start
    cdef slong stop
    cdef slong current

    def __init__(self, pAdicLazyElement elt, start, stop):
        self.elt = elt
        if start is None:
            self.start = 0
        else:
            if start >= maxordp:
                raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
            self.start = long(start)
        if stop is None:
            self.stop = self.start - 1
        else: 
            if stop >= maxordp:
                raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
            self.stop = long(stop)
            if self.stop <= self.start:
                self.stop = self.start
        self.current = self.start

    def __repr__(self):
        return "%s-adic expansion of %s" % (self.elt._parent.prime(), self.elt)

    def __len__(self):
        if self.stop < self.start:
            raise NotImplementedError("infinite sequence")
        return ZZ(self.stop - self.start)

    def __iter__(self):
        return self

    def __next__(self):
        if self.stop < self.start or self.current < self.stop:
            digit = self.elt._digit(self.current)
            self.current += 1
            return digit
        else:
            raise StopIteration
