cdef class Matrix_mod_dense(matrix_dense.Matrix_dense):
    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value):
        raise NotImplementedError

    cdef int get_unsafe_int(self, Py_ssize_t i, Py_ssize_t j):
        raise NotImplementedError
