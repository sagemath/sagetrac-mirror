from .matrix_dense cimport Matrix_dense
from sage.libs.gap.element cimport GapElement
from sage.groups.libgap_wrapper cimport ElementLibGAP

cdef class Matrix_gap(Matrix_dense):
    cdef GapElement _libgap

    cpdef GapElement gap(self)
    cdef Matrix_gap _new(self, Py_ssize_t nrows, Py_ssize_t ncols)

