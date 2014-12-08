from sage.libs.flint.types cimport fq_t, fq_ctx_struct

from sage.rings.finite_rings.element_base cimport FinitePolyExtElement

cdef class FiniteFieldElement_flint_fq(FinitePolyExtElement):
    cdef fq_t val
    cdef int initialized
    cdef fq_ctx_struct *_cparent
    cdef FiniteFieldElement_flint_fq _new(FiniteFieldElement_flint_fq self)
    cdef void set_from_fq(FiniteFieldElement_flint_fq self, fq_t val) except *
    cdef void construct_from(FiniteFieldElement_flint_fq self, object x) except *
