from .matrix cimport Matrix

cdef class Matrix_dense(Matrix):
    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value)
    cdef void set_unsafe_ui(self, Py_ssize_t i, Py_ssize_t j, unsigned long value)
    cdef unsigned long get_unsafe_ui(self, Py_ssize_t i, Py_ssize_t j)

