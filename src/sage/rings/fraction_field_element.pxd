from sage.structure.element cimport FieldElement

cdef class FractionFieldElement(FieldElement):

    cdef object __numerator
    cdef object __denominator
    cdef bint _is_reduced

    cpdef bint is_zero(self)
    cpdef bint is_one(self)

    cpdef reduce(self)
    cdef normalize_unit_denominator(self)
