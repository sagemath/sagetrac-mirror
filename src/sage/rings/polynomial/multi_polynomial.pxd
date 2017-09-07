from sage.structure.element cimport CommutativeRingElement

cdef class MPolynomial(CommutativeRingElement):
    cdef long _hash_c(self) except -1
    cpdef _mod_(self, right)
