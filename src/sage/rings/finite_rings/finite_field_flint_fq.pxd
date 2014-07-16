include "sage/libs/flint/fq.pxi"

from finite_field_base cimport FiniteField

cdef class FiniteField_flint_fq(FiniteField):
    cdef fq_ctx_struct *_ctx
    #cdef fq_ctx_t _ctx

    cdef public object _modulus
    cdef public object _degree
    cdef public object _gen
    cdef public int __hash
