from sage.structure.element cimport AlgebraElement


cdef class BoundedSeries(AlgebraElement):
    cdef char __is_gen
    cdef _prec
    cdef _valuation_final_terms
