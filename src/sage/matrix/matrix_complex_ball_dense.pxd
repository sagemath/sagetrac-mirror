from sage.libs.arb.acb_mat cimport acb_mat_t
cimport matrix_dense
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.structure.parent cimport Parent

cdef void matrix_to_acb_mat(acb_mat_t target, source)
cdef Matrix_generic_dense acb_mat_to_matrix(
    acb_mat_t source, Parent CIF)

cdef class Matrix_complex_ball_dense(matrix_dense.Matrix_dense):
    cdef acb_mat_t value
    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x)
    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j)
