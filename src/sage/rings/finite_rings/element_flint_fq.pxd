include "sage/libs/flint/fq.pxi"

from sage.rings.finite_rings.element_base cimport FinitePolyExtElement

cdef class FiniteFieldElement_flint_fq(FinitePolyExtElement):
    cdef fq_t val
    cdef fq_ctx_struct *_cparent
    #cdef fq_ctx_t _cparent
    cdef FiniteFieldElement_flint_fq _new(FiniteFieldElement_flint_fq self)
    cdef void construct_from(FiniteFieldElement_flint_fq self, object x) except *
