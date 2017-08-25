from sage.structure.element cimport CommutativeRingElement
from sage.rings.padics.pool cimport Pool

cdef class LocalGenericElement(CommutativeRingElement):
    cdef Pool _pool
