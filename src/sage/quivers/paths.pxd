from sage.structure.element cimport MonoidElement, Element
from sage.data_structures.bounded_integer_sequences cimport biseq_t

cdef class QuiverPath(MonoidElement):
    cdef biseq_t _path
    cdef int _start, _end
    cdef QuiverPath _new_(self, int start, int end)
    cdef _mul_(self, other)
    cdef _mod_(self, right)
    cdef tuple complement(self, QuiverPath subpath)
    cdef bint has_subpath(self, QuiverPath subpath) except -1
    cdef bint has_prefix(self, QuiverPath subpath) except -1
