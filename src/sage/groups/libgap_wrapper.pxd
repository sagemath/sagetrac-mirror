from sage.structure.element cimport MultiplicativeGroupElement
from sage.libs.gap.element cimport GapElement


cdef class ElementLibGAP(MultiplicativeGroupElement):
    cdef GapElement _libgap
    cdef GapElement gap(self)
    cdef _mul_(self, other)
