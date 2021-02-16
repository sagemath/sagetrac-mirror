from .matrix_dense cimport Matrix_dense
from gappy.gapobj cimport GapObj

cdef class Matrix_gap(Matrix_dense):
    cdef GapObj _libgap

    cpdef GapObj gap(self)
    cdef Matrix_gap _new(self, Py_ssize_t nrows, Py_ssize_t ncols)

