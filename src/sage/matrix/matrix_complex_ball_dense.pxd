from sage.libs.arb.types cimport acb_mat_t
cimport matrix_dense

cdef class Matrix_complex_ball_dense(matrix_dense.Matrix_dense):
    cdef acb_mat_t value
