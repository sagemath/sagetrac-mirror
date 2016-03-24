from sage.structure.element cimport FieldElement

cdef enum ReductionType:
    Unreduced = 0,
    ReductionFailed,
    Reduced,
    NormalForm

cdef class FractionFieldElement(FieldElement):
    cdef object __numerator
    cdef object __denominator
    cdef ReductionType __reduction

