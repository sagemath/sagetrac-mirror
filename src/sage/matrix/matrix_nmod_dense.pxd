from sage.libs.gmp.types cimport *
from sage.libs.flint.types cimport nmod_mat_t

from .matrix_dense cimport Matrix_dense
from  sage.rings.finite_rings.integer_mod cimport NativeIntStruct

cdef class Matrix_nmod_dense(Matrix_dense):
    cdef nmod_mat_t _matrix
    cdef NativeIntStruct _modulus
    cdef void set_unsafe_si(self, Py_ssize_t i, Py_ssize_t j, long value)
    cdef long get_unsafe_si(self, Py_ssize_t i, Py_ssize_t j)

    cdef Matrix_nmod_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols)
    cdef void strong_echelon_form(self)
    cdef void howell_form(self)


