from sage.structure.element cimport CommutativeRingElement

cdef class MPolynomial(CommutativeRingElement):
    cdef long _hash_c(self) except -1
    cdef _mod_(self, right)
    cdef dict _mpoly_dict_recursive(self, tuple vars=*, base_ring=*)

