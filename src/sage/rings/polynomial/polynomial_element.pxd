from sage.structure.element import Element, CommutativeAlgebraElement
from sage.structure.element cimport Element, CommutativeAlgebraElement, ModuleElement
from sage.structure.parent cimport Parent
from sage.rings.integer cimport Integer
from .polynomial_compiled cimport CompiledPolynomialFunction


cdef class Polynomial(CommutativeAlgebraElement):
    cdef Polynomial _new_generic(self, list coeffs)
    cdef char _is_gen
    cdef CompiledPolynomialFunction _compiled
    cdef Polynomial truncate(self, long n)
    cdef Polynomial inverse_series_trunc(self, long prec)
    cdef long _hash_c(self) except -1
    cdef constant_coefficient(self)
    cdef Polynomial _new_constant_poly(self, a, Parent P)
    cdef list list(self, bint copy=*)
    cdef _mul_generic(self, right)
    cdef _square_generic(self)

    cdef bint is_zero(self) except -1
    cdef bint is_one(self) except -1
    cdef bint is_term(self) except -1

    cdef dict _mpoly_dict_recursive(self, tuple variables=*, base_ring=*)

    cdef _add_(self, other)
    cdef _mul_(self, other)
    cdef _floordiv_(self, right)
    cdef Polynomial _mul_trunc_(self, Polynomial right, long n)
    cdef Polynomial _power_trunc(self, unsigned long n, long prec)
    cdef Polynomial _mul_term(self, Polynomial term, bint term_on_right)

    # UNSAFE, only call from an inplace operator
    # may return a new element if not possible to modify inplace
    cdef _inplace_truncate(self, long n)

    cdef get_coeff_c(self, Py_ssize_t i)
    cdef get_unsafe(self, Py_ssize_t i)
    cdef long number_of_terms(self)

    # See 23227
    cdef _add_(self, right)
    cdef _mul_(self, right)
    cdef _floordiv_(self, right)

    cdef public dict __cached_methods

cdef class Polynomial_generic_dense(Polynomial):
    cdef Polynomial_generic_dense _new_c(self, list coeffs, Parent P)
    cdef list __coeffs
    cdef int __normalize(self) except -1
    cdef list list(self, bint copy=*)

cdef class Polynomial_generic_dense_inexact(Polynomial_generic_dense):
    pass

cdef is_Polynomial(f)
cdef Polynomial generic_power_trunc(Polynomial p, Integer n, long prec)
cdef list _dict_to_list(dict x, zero)

cdef bint polynomial_is_variable(x)

