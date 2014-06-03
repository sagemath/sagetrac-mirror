from sage.structure.element cimport FieldElement

from sage.libs.flint.fmpz_poly_q cimport *

cdef class RationalFunctionQQ(FieldElement):
    cdef fmpz_poly_q_t _func

    cdef RationalFunctionQQ _new(self)
