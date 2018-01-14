from sage.libs.linbox.givaro cimport Modular_uint64

from sage.libs.linbox.linbox cimport (SparseMatrix_Modular_uint64,
        DenseVector_Modular_uint64)

from sage.modules.vector_modn_dense cimport Vector_modn_dense

cdef class LinBoxSolver_modn_sparse(object):
    r"""
    Solving linear equations in Sage using LinBox.
    """
    cdef object Vin
    cdef object Vout
    cdef size_t nrows
    cdef SparseMatrix_Modular_uint64 * A
    cdef DenseVector_Modular_uint64 * b
    cdef DenseVector_Modular_uint64 * res
    cdef Vector_modn_dense res_sage
    cdef Modular_uint64 * F
