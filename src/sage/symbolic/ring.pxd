from sage.symbolic.expression cimport Expression
from sage.rings.ring cimport CommutativeRing

cdef class SymbolicRing(CommutativeRing):
    cdef public dict _symbols

    cpdef Expression symbol(self, name=*, latex_name=*, domain=*)
