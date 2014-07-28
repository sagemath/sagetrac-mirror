from sage.libs.flint.fq cimport *

from finite_field_base cimport FiniteField

cdef class FiniteField_flint_fq(FiniteField):
    cdef fq_ctx_struct *_ctx
    #cdef fq_ctx_t _ctx
    cdef int _ctx_initialized

    cdef public object _modulus
    cdef public object _degree
    cdef public object _gen
    cdef public int __hash
