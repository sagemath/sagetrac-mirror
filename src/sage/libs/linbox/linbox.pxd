# distutils: extra_compile_args = LINBOX_CFLAGS
# distutils: libraries = LINBOX_LIBRARIES
# distutils: library_dirs = LINBOX_LIBDIR
# distutils: language = c++ 

from sage.libs.gmp.types cimport mpz_t
from sage.modules.vector_modn_sparse cimport c_vector_modint

from sage.matrix.matrix_integer_dense cimport mod_int

ctypedef enum SparseAlgorithm:
    SparseElimination
    Wiedemann

cdef class Linbox_matrix_modn_sparse:
    cdef size_t nrows, ncols
    cdef void *_M

    cdef set(self, int modulus, size_t nrows, size_t ncols, c_vector_modint *rows)
    cdef size_t rank(self, SparseAlgorithm algorithm)
    cdef void solve(self, c_vector_modint **x, c_vector_modint *b, SparseAlgorithm algorithm)

cdef class Linbox_integer_dense:
    cdef void* _M
    cdef size_t nrows, ncols

    cdef set(self, size_t nrows, size_t ncols, mpz_t ** matrix)
    cdef size_t rank(self)
    cdef det(self)
    cdef minpoly(self)
    cdef charpoly(self)

