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
    cpdef _shift_mod(self, mp_limb_t modulus, mp_limb_t shift=*, bint mul=*, bint domod=*)
