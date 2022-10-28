from sage.structure.element cimport Element
from sage.categories.map cimport Map


cdef class CCtoCDF(Map):

    cdef Element _call_(self, x)
