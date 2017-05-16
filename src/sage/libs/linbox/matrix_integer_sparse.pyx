# distutils: extra_compile_args = LINBOX_CFLAGS
# distutils: libraries = LINBOX_LIBRARIES
# distutils: library_dirs = LINBOX_LIBDIR
# distutils: language = c++
r"""
LinBox interface to Sage sparse integer matrices
"""
#*****************************************************************************
#       Copyright (C) 2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function

from libcpp.vector cimport vector as cppvector
from libcpp.string cimport string as cppstring
from libcpp.typeinfo cimport type_info
from cython.operator cimport typeid

from sage.ext.stdsage cimport PY_NEW

from sage.libs.gmp.types cimport mpz_srcptr, mpz_ptr
from sage.libs.gmp.mpz cimport mpz_set

from .methods cimport LinBoxMethod

from sage.rings.integer cimport Integer

from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.modules.vector_integer_sparse cimport mpz_vector,  mpz_vector_get_entry, mpz_vector_set_entry
from sage.matrix.matrix_integer_sparse cimport Matrix_integer_sparse


cdef extern from "<sstream>" namespace "std":
    cdef cppclass ostringstream:
        cppstring str()


cdef extern from "givaro/givconfig.h":
    pass
cdef extern from "linbox/linbox-config.h":
    pass


cdef extern from "gmp++/gmp++.h":
    cdef cppclass GivaroInteger "Givaro::Integer":
        mpz_ptr get_mpz()
        mpz_srcptr get_mpz_const()


cdef extern from "givaro/zring.h":
    cdef cppclass GivaroIntegerRing "Givaro::ZRing<Givaro::Integer>":
        ctypedef GivaroInteger Element


cdef extern from "linbox/matrix/sparse-matrix.h":
    cdef cppclass LinBoxIntegerSparseMatrix "LinBox::SparseMatrix<Givaro::ZRing<Givaro::Integer>>":
        ctypedef GivaroIntegerRing Field
        ctypedef GivaroInteger Element
        LinBoxIntegerSparseMatrix(Field &F, size_t m, size_t n)
        size_t rowdim()
        size_t coldim()
        void setEntry(size_t i, size_t j, Element &a)
        Element &getEntry(size_t i, size_t j)
        Field& field()


cdef extern from "linbox/vector/vector.h":
    cdef cppclass LinBoxIntegerDenseVector "LinBox::DenseVector<Givaro::ZRing<Givaro::Integer>>":
        ctypedef GivaroIntegerRing Field
        ctypedef GivaroInteger Element
        cppclass iterator:
            GivaroInteger operator*()
            iterator operator++()
            bint operator==()
            bint operator!=()
        iterator begin()
        iterator end()
        LinBoxIntegerDenseVector (Field &F)
        LinBoxIntegerDenseVector (Field &F, long& m)
        LinBoxIntegerDenseVector (Field &F, cppvector[GivaroInteger]&)
        cppvector[GivaroInteger]& refRep()


cdef extern from "wrap.h":
    void fprint_linbox_vector "WRAP_OUT" (ostringstream, LinBoxIntegerDenseVector)
    void fprint_linbox_matrix "WRAP_OUT" (ostringstream, LinBoxIntegerSparseMatrix)
    void fprint_string "WRAP_OUT" (ostringstream, cppstring)
    void fprint_linbox_integer "WRAP_OUT" (ostringstream, GivaroInteger)


cdef extern from "linbox/solutions/solve.h":
    LinBoxIntegerDenseVector& LinBoxIntegerSparse_solve "LinBox::solve" (LinBoxIntegerDenseVector &,
                                                     GivaroInteger &,
                                                     LinBoxIntegerSparseMatrix &,
                                                     LinBoxIntegerDenseVector &,
                                                     LinBoxMethod.BlasElimination) except +RuntimeError

    LinBoxIntegerDenseVector& LinBoxIntegerSparse_solve "LinBox::solve" (LinBoxIntegerDenseVector &,
                                                     GivaroInteger &,
                                                     LinBoxIntegerSparseMatrix &,
                                                     LinBoxIntegerDenseVector &,
                                                     LinBoxMethod.SparseElimination) except +RuntimeError

    LinBoxIntegerDenseVector& LinBoxIntegerSparse_solve "LinBox::solve" (LinBoxIntegerDenseVector &,
                                                     GivaroInteger &,
                                                     LinBoxIntegerSparseMatrix &,
                                                     LinBoxIntegerDenseVector &,
                                                     LinBoxMethod.Wiedemann) except +RuntimeError

###############################################################################
# Begining of code                                                            #
###############################################################################

cdef LinBoxIntegerSparseMatrix * new_linbox_matrix_from_sage_integer_sparse(Matrix_integer_sparse m):
    cdef GivaroIntegerRing ZZ
    cdef LinBoxIntegerSparseMatrix * A
    A = new LinBoxIntegerSparseMatrix(ZZ, m._nrows, m._ncols)
    cdef size_t i,j,k
    cdef mpz_vector * v
    cdef GivaroInteger t
    for i in range(m._nrows):
        v = m._matrix + i  # a mpz_vector
        for k in range(v.num_nonzero):
            j = v.positions[k]
            mpz_set(t.get_mpz(), v.entries[k])
            A.setEntry(i, j, t)
    return A

cdef LinBoxIntegerDenseVector * new_linbox_vector_from_sage_integer_dense(Vector_integer_dense v):
    cdef cppvector[GivaroInteger] * vec = new cppvector[GivaroInteger](v._degree)
    cdef size_t i
    for i in range(v._degree):
        mpz_set(vec[0][i].get_mpz(), v._entries[i])

    cdef GivaroIntegerRing ZZ
    cdef LinBoxIntegerDenseVector * V = new LinBoxIntegerDenseVector(ZZ, vec[0])
    del vec
    return V

def linbox_integer_sparse_solve(Matrix_integer_sparse m, Vector_integer_dense b, method):
    r"""
    Return a pair ``(a, d)`` so that ``d * b = m * a``

    INPUT:

    - ``m`` -- a sparse integer matrix

    - ``b`` -- a dense integer vector

    - ``method`` -- either ``'blas_elimination'``,
      ``'sparse_elimination'`` or ``'wiedemann'``

    EXAMPLES::

        sage: from sage.libs.linbox.matrix_integer_sparse import linbox_integer_sparse_solve
        sage: m = matrix(ZZ, 4, sparse=True)
        sage: m[0,0] = m[1,2] = m[2,0] = m[3,3] = 2
        sage: m[0,2] = m[1,1] = -1
        sage: m[2,3] = m[3,0] = -3

        sage: b0 = vector((1,1,1,1))
        sage: linbox_integer_sparse_solve(m, b0, 'blas_elimination')
        ((-1, -7, -3, -1), 1)
        sage: linbox_integer_sparse_solve(m, b0, 'sparse_elimination')
        ((-1, -7, -3, -1), 1)
        sage: linbox_integer_sparse_solve(m, b0, 'wiedemann')
        ((-1, -7, -3, -1), 1)

        sage: b1 = vector((1,2,3,4))
        sage: linbox_integer_sparse_solve(m, b1, 'blas_elimination')
        ((-18, -92, -41, -17), 5)
        sage: linbox_integer_sparse_solve(m, b1, 'sparse_elimination')
        ((-18, -92, -41, -17), 5)
        sage: linbox_integer_sparse_solve(m, b1, 'wiedemann')
        ((-18, -92, -41, -17), 5)

        sage: a1, d1 = linbox_integer_sparse_solve(m, b1)
        sage: d1 * b1 == m * a1
        True

    Rectangular systems are weird::

        sage: m = matrix(ZZ, 1, 2, [1, 1], sparse=True)
        sage: b = vector(ZZ, [1])
        sage: linbox_integer_sparse_solve(m, b, 'blas_elimination')
        ((1, 0), 1)
        sage: linbox_integer_sparse_solve(m, b, 'sparse_elimination')  # bug!!!
        ((0, 0), 0)
        sage: linbox_integer_sparse_solve(m, b, 'wiedemann')           # bug!!!
        ((0, 0), 0)

        sage: m = matrix(ZZ, 2, 1, [1, 1], sparse=True)
        sage: b = vector(ZZ, [1,1])
        sage: linbox_integer_sparse_solve(m, b, 'blas_elimination')
        ((1), 1)

        sage: m = matrix(ZZ, 2, 1, [1, 1], sparse=True)
        sage: linbox_integer_sparse_solve(m, b, 'sparse_elimination') # bug!!!
        ((0), 0)

        sage: m = matrix(ZZ, 2, 1, [1, 1], sparse=True)
        sage: linbox_integer_sparse_solve(m, b, 'wiedemann')          # bug!!!
        ((0), 0)
    """
    if b._degree != m._nrows:
        raise ValueError("dimension mismatch, the number of rows of the matrix must be the degree of the vector")

    cdef GivaroIntegerRing ZZ
    cdef LinBoxIntegerDenseVector * A = new LinBoxIntegerDenseVector(ZZ, m._ncols)
    cdef LinBoxIntegerDenseVector * B = new_linbox_vector_from_sage_integer_dense(b)
    cdef LinBoxIntegerSparseMatrix * M = new_linbox_matrix_from_sage_integer_sparse(m)
    cdef GivaroInteger D

    # debug info
    print("type(ZZ) is %s" % typeid(ZZ).name())
    print("type(M) is %s" % typeid(M[0]).name()) 
    print("type(A) is %s" %typeid(A[0]).name())
    print("type(D) is %s" %typeid(D).name())
    print("type(B) is %s" %typeid(B[0]).name())
    print("type(BlasElim) is %s" %typeid(LinBoxMethod.BlasElimination()).name())
    print("type(SparseElim) is %s" %typeid(LinBoxMethod.SparseElimination()).name())
    print("type(Wiedemann) is %s" %typeid(LinBoxMethod.Wiedemann()).name())

    # solve
    if method == 'blas_elimination':
        print('solve using BlasElimination')
        LinBoxIntegerSparse_solve(A[0], D, M[0], B[0], LinBoxMethod.BlasElimination())
    elif method == 'sparse_elimination':
        print('solve using SparseElimination')
        LinBoxIntegerSparse_solve(A[0], D, M[0], B[0], LinBoxMethod.SparseElimination())
    elif method == 'wiedemann':
        print('solve using Wiedemann')
        LinBoxIntegerSparse_solve(A[0], D, M[0], B[0], LinBoxMethod.Wiedemann())
    else:
        print('Error')
        raise ValueError("method '{}' not available".format(method))

    # print information on stdout
    cdef ostringstream out
    fprint_string(out, "M=\n")
    fprint_linbox_matrix(out, M[0])
    fprint_string(out, "\nB=")
    fprint_linbox_vector(out, B[0])
    fprint_string(out, "\nA=")
    fprint_linbox_vector(out, A[0])
    fprint_string(out, " / ")
    fprint_linbox_integer(out, D)
    print(out.str())
