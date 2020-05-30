from sage.libs.arb.types cimport arb_mat_t
from .matrix_dense cimport Matrix_dense

cdef class Matrix_real_ball_dense(Matrix_dense):
    cdef arb_mat_t value
