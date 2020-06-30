from sage.structure.element cimport MonoidElement
from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.structure.parent cimport Parent

cdef class FreeAbelianMonoidElement(MonoidElement):
    cdef Vector_integer_dense _element_vector

