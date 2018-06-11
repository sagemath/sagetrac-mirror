cimport sage.rings.ring

cdef class MPolynomialRing_base(sage.rings.ring.CommutativeRing):
    cdef object __ngens
    cdef object __term_order
    cdef public object _has_singular
    cdef public object _magma_gens
    cdef public dict _magma_cache
