from sage.rings.ring cimport CommutativeRing
from sage.rings.padics.pool cimport Pool

cdef class LocalGeneric(CommutativeRing):
    cdef Pool _pool_empty
    cdef Pool _pool
