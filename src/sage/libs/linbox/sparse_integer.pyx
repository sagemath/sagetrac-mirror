# distutils: extra_compile_args = LINBOX_CFLAGS
# distutils: libraries = LINBOX_LIBRARIES
# distutils: library_dirs = LINBOX_LIBDIR
# distutils: language = c++

r"""
Interface between sparse matrices and linbox
"""

from sage.libs.gmp.types cimport mpz_t, mpz_srcptr, mpz_ptr
from sage.libs.gmp.mpz cimport mpz_set
from sage.libs.flint.types cimport fmpz, fmpz_t
from sage.libs.flint.fmpz cimport fmpz_get_mpz, fmpz_set_mpz
from sage.rings.integer cimport Integer

##########################################################################
## Sparse matrices over ZZ
##########################################################################

cdef extern from "gmp++/gmp++.h":
    cdef cppclass GivaroInteger "Givaro::Integer":
        mpz_ptr get_mpz()
        mpz_srcptr get_mpz_const()

cdef extern from "givaro/zring.h":
    cdef cppclass GivaroIntegerRing "Givaro::ZRing<Givaro::Integer>":
        ctypedef GivaroInteger Element

cdef extern from "linbox/matrix/sparse-matrix.h":
    ## template <class _Field>
    ## using DenseMatrix = BlasMatrix<_Field> ;
    ##
    # indirect include from linbox/matrix/densematrix/blas-matrix.h
    ##
    ## template <class _Field, class _Storage>
    ## class BlasMatrix
    cdef cppclass LinBoxIntegerSparseMatrix "LinBox::SparseMatrix<Givaro::ZRing<Givaro::Integer>>":
        ctypedef GivaroIntegerRing Field
        ctypedef GivaroInteger Element
        size_t rowdim()
        size_t coldim()
        LinBoxIntegerSparseMatrix(Field &F, size_t m, size_t n)
        void setEntry(size_t i, size_t j, Element &a)
        Element &getEntry(size_t i, size_t j)

cdef extern from "linbox/solutions/rank.h":
    ## template <class Blackbox, class Method, class DomainCategory>
    ## inline unsigned long &rank (unsigned long &r, const Blackbox &A,
    ##                             const DomainCategory &tag, const Method &M);
    unsigned long & LinBoxIntegerSparse_rank "LinBox::rank" (unsigned long &, LinBoxIntegerSparseMatrix)

# return the rank of A
cdef unsigned long linbox_sparse_mat_rank(dict A, size_t nrows, size_t ncols):
    cdef GivaroIntegerRing ZZ
    cdef LinBoxIntegerSparseMatrix * LBA
    cdef unsigned long r = 0

    LBA = new LinBoxIntegerSparseMatrix(ZZ, nrows, ncols)
    cdef size_t i,j
    cdef GivaroInteger t

    for i,j in A:
        mpz_set(t.get_mpz(), (<Integer> A[i,j]).value)
        LBA.setEntry(i, j, t)
    LinBoxIntegerSparse_rank(r, LBA[0])

    del LBA

    return r

