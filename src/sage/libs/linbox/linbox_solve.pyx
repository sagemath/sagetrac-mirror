r"""
Deprecated linbox interface

This module is a work around the compilation problems encountered
in :trac:`24544`.
"""

## NOTE: The sig_on()/sig_off() stuff can't go in here -- it has to be in the
## code that calls these functions.  Otherwise strangely objects get left
## in an incorrect state.

from sage.libs.gmp.mpz cimport *
from sage.rings.integer cimport Integer
from sage.misc.misc import verbose, get_verbose
from sage.modules.vector_modn_sparse cimport set_entry

##########################################################################
## Sparse matrices modulo p.
##########################################################################

cdef extern from "linbox/linbox-sage.h":
    ctypedef struct vector_uint "std::vector<unsigned int>":
        void (*push_back)(unsigned int)
        int (*get "operator[]") (size_t i)
        int (*size)()


    int linbox_modn_sparse_matrix_rank(unsigned long modulus, size_t nrows, size_t ncols, void
* rows, int reorder) #, int **pivots)
    vector_uint linbox_modn_sparse_matrix_solve(unsigned long modulus, size_t numrows, size_t numcols, void *a,  void *b, int method)

cdef class Linbox_modn_sparse:
    cdef set(self, int modulus, size_t nrows, size_t ncols, c_vector_modint *rows):
        self.rows = rows
        self.nrows = nrows
        self.ncols = ncols
        self.modulus = modulus

    cdef void solve(self, c_vector_modint **x, c_vector_modint *b, int method):
        """
        """
        cdef vector_uint X
        X = linbox_modn_sparse_matrix_solve(self.modulus, self.nrows, self.ncols, self.rows, b, method)

        for i from 0 <= i < X.size():
            set_entry(x[0], i, X.get(i))
