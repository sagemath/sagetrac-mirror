from sage.structure.element cimport MonoidElement, Element
from sage.data_structures.bounded_integer_sequences cimport *
from sage.libs.gmp.types cimport *
from sage.libs.gmp.mpn cimport mpn_cmp

include "sage/ext/python.pxi"
include "sage/ext/cdefs.pxi"
include "sage/ext/stdsage.pxi"
include "sage/libs/ntl/decl.pxi"
include "sage/ext/interrupt.pxi"

cdef extern from "Python.h":
    bint PySlice_Check(PyObject* ob)

cdef class QuiverPath(MonoidElement):
    cdef biseq_t _path
    cdef int _start, _end
    cdef QuiverPath _new_(self, int start, int end)
    cpdef bint has_subpath(self, QuiverPath subpath) except -1
    cpdef bint has_initial_segment(self, QuiverPath subpath) except -1

cpdef QuiverPath NewQuiverPath(Q, start, end, data, itembitsize, length)
