from sage.structure.element cimport MultiplicativeGroupElement, Element, MonoidElement, Matrix
from sage.groups.libgap_wrapper cimport ElementLibGAP

cdef is_MatrixGroupElement(x)

cdef class MatrixGroupElement_generic(MultiplicativeGroupElement):
    cdef public Matrix _matrix

    cdef _act_on_(self, x, bint self_on_left)
    cdef _mul_(self, other)
    cdef list list(self)

cdef class MatrixGroupElement_gap(ElementLibGAP):
    cdef _act_on_(self, x, bint self_on_left)
    cdef list list(self)

