include "sage/ext/cdefs.pxi"

cimport matrix_dense

cdef class Matrix_integer_2x2(matrix_dense.Matrix_dense):
    cdef __mpz_struct *a
    cdef __mpz_struct *b
    cdef __mpz_struct *c
    cdef __mpz_struct *d
    cdef mpz_t _entries[4]

    cdef Matrix_integer_2x2 _new_c(self)


