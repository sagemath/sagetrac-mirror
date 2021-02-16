from sage.structure.element cimport MultiplicativeGroupElement
from gappy.gapobj cimport GapObj


cdef class ElementLibGAP(MultiplicativeGroupElement):
    cdef GapObj _libgap
    cpdef GapObj gap(self)
    cpdef _mul_(self, other)
