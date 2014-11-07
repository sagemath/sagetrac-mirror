include 'sage/modules/vector_modn_sparse_h.pxi'

cimport matrix_sparse
cimport matrix_dense
from matrix_mod_dense cimport Matrix_mod_dense

cdef class Matrix_modn_sparse(matrix_sparse.Matrix_sparse):
    cdef c_vector_modint* rows
    cdef public int p
    cdef swap_rows_c(self, Py_ssize_t n1, Py_ssize_t n2)
    cdef set_block_unsafe(self, Py_ssize_t row, Py_ssize_t col, Matrix_mod_dense block)

    cdef _init_linbox(self)
