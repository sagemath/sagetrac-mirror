
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense

cdef class Matrix_test_sage(Matrix_generic_dense):
    cpdef test_func(self)
    cpdef row_echelon_test(self, pivot_choice =* ,  l =*, compute_pivots =*)
    cdef alter_precision(self,i,j,col_ind)
    cpdef row_echelon_test2(self,   l =*, compute_pivots =*)