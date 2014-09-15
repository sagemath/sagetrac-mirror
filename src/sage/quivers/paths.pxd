from sage.structure.element cimport MonoidElement, Element
from sage.misc.bounded_integer_sequences cimport biseq, biseq_t, dealloc_biseq, getitem_biseq, concat_biseq, startswith_biseq, contains_biseq, max_overlap_biseq, slice_biseq, list_to_biseq, biseq_to_list

include "sage/ext/python.pxi"
include "sage/ext/cdefs.pxi"
include "sage/ext/stdsage.pxi"
include "sage/libs/ntl/decl.pxi"
include "sage/ext/interrupt.pxi"
cdef extern from "gmp.h":
    cdef int mp_bits_per_limb

cdef extern from "Python.h":
    bint PySlice_Check(PyObject* ob)

cdef extern from "mpz_pylong.h":
    cdef long mpz_pythonhash(mpz_t src)

cdef class QuiverPath(MonoidElement):
    cdef biseq_t _path
    cdef int _start, _end
    cdef QuiverPath _new_(self, int start, int end, biseq_t data)
    cpdef int has_subpath(self, QuiverPath subpath) except -1
    cpdef int has_initial_segment(self, QuiverPath subpath) except -1

cpdef QuiverPath NewQuiverPath(Q, start, end, data, bitsize, itembitsize, length)
