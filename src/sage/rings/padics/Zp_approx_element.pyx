from libc.stdint cimport uintptr_t

from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.types cimport mpz_t

from sage.rings.integer cimport Integer


cdef extern from "pow_prime.h":
    cdef cppclass PowPrime:
        PowPrime() except +
        void init(mpz_t prime, const int cache_loglength)
        mpz_t* power(const long n)

cdef extern from "Zp_approximated_element.h":
    cdef cppclass ZpApproximatedElement:
        ZpApproximatedElement() except +
        ZpApproximatedElement(mpz_t value, PowPrime*) except +
        char* repr()
        long valuation()
        void ireduce(const long, bool)
        ZpApproximatedElement add(ZpApproximatedElement&, const long, bool)
        void iadd(ZpApproximatedElement&, const long, bool)
        void operator+=(ZpApproximatedElement&)
        ZpApproximatedElement operator+(ZpApproximatedElement&)



prime = Integer(3)
cdef PowPrime pow
pow.init(prime.value, 7)

cdef class ZpApproxElement:
    cdef ZpApproximatedElement _value

    def __init__(self, v):
        cdef Integer vv = Integer(v)
        self._value = ZpApproximatedElement(vv.value, &pow)

    cpdef ZpApproxElement _add_(ZpApproxElement self, ZpApproxElement other):
        cdef ZpApproxElement ans = PY_NEW(ZpApproxElement)
        ans._value = self._value + other._value
        ans._value.ireduce(20, True)
        return ans

    cpdef ZpApproxElement _add2_(ZpApproxElement self, ZpApproxElement other):
        cdef ZpApproxElement ans = PY_NEW(ZpApproxElement)
        ans._value = self._value
        ans._value.iadd(other._value, 20, 1)
        return ans

    def __add__(self, other):
        return self._add_(other)

    def __repr__(self):
        return str(self._value.repr())
