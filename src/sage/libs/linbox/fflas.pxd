# distutils: extra_compile_args = FFLASFFPACK_CFLAGS
# distutils: libraries = FFLASFFPACK_LIBRARIES
# distutils: library_dirs = FFLASFFPACK_LIBDIR
# distutils: language = c++ 

from .givaro cimport Modular_double, Modular_float, Dense, Sparse
from .givaro cimport givvector, Poly1Dom
from libcpp.vector cimport vector

ctypedef Poly1Dom[Modular_double, Dense] PolynomialRing_Modular_double
ctypedef Poly1Dom[Modular_float, Dense] PolynomialRing_Modular_float

ctypedef givvector[Modular_double.Element] ModDoubleDensePolynomial
ctypedef givvector[Modular_float.Element] ModFloatDensePolynomial

cdef extern from "fflas-ffpack/paladin/blockcuts.inl" namespace "FFLAS::CuttingStrategy":
    ctypedef struct Single:
        pass
    ctypedef struct Row:
        pass
    ctypedef struct Column:
        pass
    ctypedef struct Block:
        pass
    ctypedef struct Recursive:
        pass

cdef extern from "fflas-ffpack/paladin/blockcuts.inl" namespace "FFLAS::StrategyParameter":
    ctypedef struct Fixed:
        pass
    ctypedef struct Threads:
        pass
    ctypedef struct Grain:
        pass
    ctypedef struct TwoD:
        pass
    ctypedef struct TwoDAdaptive:
        pass
    ctypedef struct ThreeD:
        pass
    ctypedef struct ThreeDInPlace:
        pass
    ctypedef struct ThreeDAdaptive:
        pass

cdef extern from "fflas-ffpack/paladin/blockcuts.inl" namespace "FFLAS::ParSeqHelper":
    cdef cppclass Parallel[C,P]:
        Parallel(size_t n)
        size_t& set_numthreads(size_t n)

    cdef cppclass Sequential:
        Sequential()

cdef extern from "fflas-ffpack/fflas/fflas_helpers.h" namespace "FFLAS::MMHelperAlgo":
    ctypedef struct Auto:
        pass
    ctypedef struct Classic:
        pass
    ctypedef struct Winograd:
        pass
    ctypedef struct WinogradPar:
        pass
    ctypedef struct Bini:
        pass

cdef extern from "fflas-ffpack/fflas-ffpack.h" namespace "FFLAS":
    ctypedef enum FFLAS_TRANSPOSE:
        FflasNoTrans
        FflasTrans

    ctypedef enum FFLAS_SIDE:
        FflasRight

    # double
    void fgemv (Modular_double F, FFLAS_TRANSPOSE transA,
             size_t nrows, size_t ncols,
             Modular_double.Element alpha, Modular_double.Element* A,
             size_t lda, Modular_double.Element* X, size_t incX,
             Modular_double.Element beta, Modular_double.Element* Y,
             size_t incY)

    Modular_double.Element* fgemm (Modular_double F,
             FFLAS_TRANSPOSE transA, FFLAS_TRANSPOSE transB,
             size_t nrowsA, size_t ncolsB, size_t ncolsA,
             Modular_double.Element alpha, Modular_double.Element* A,
             size_t A_stride, Modular_double.Element* B, int B_stride,
             Modular_double.Element beta, Modular_double.Element* C,
             size_t C_stride)

    Modular_double.Element* fgemm[C,P] (Modular_double F,
             FFLAS_TRANSPOSE transA, FFLAS_TRANSPOSE transB,
             size_t nrowsA, size_t ncolsB, size_t ncolsA,
             Modular_double.Element alpha, Modular_double.Element* A,
             size_t A_stride, Modular_double.Element* B, int B_stride,
             Modular_double.Element beta, Modular_double.Element* C,
             size_t C_stride, Parallel[C,P]& H)


    # float
    void fgemv (Modular_float F, FFLAS_TRANSPOSE transA,
             size_t nrows, size_t ncols,
             Modular_float.Element alpha, Modular_float.Element* A,
             size_t lda, Modular_float.Element* X, size_t incX,
             Modular_float.Element beta, Modular_float.Element* Y,
             size_t incY)

    Modular_float.Element* fgemm (Modular_float F,
             FFLAS_TRANSPOSE transA, FFLAS_TRANSPOSE transB,
             size_t nrowsA, size_t ncolsB, size_t ncolsA,
             Modular_float.Element alpha, Modular_float.Element* A,
             size_t A_stride, Modular_float.Element* B, int B_stride,
             Modular_float.Element beta, Modular_float.Element* C,
             size_t C_stride)

    Modular_float.Element* fgemm[C,P] (Modular_float F,
             FFLAS_TRANSPOSE transA, FFLAS_TRANSPOSE transB,
             size_t nrowsA, size_t ncolsB, size_t ncolsA,
             Modular_float.Element alpha, Modular_float.Element* A,
             size_t A_stride, Modular_float.Element* B, int B_stride,
             Modular_float.Element beta, Modular_float.Element* C,
             size_t C_stride, Parallel[C,P]& H)



cdef extern from "fflas-ffpack/fflas-ffpack.h" namespace "FFPACK":
    # double
    bint IsSingular (Modular_double F,
                     size_t nrows, size_t ncols, Modular_double.Element* A,
                     size_t A_stride)

    Modular_double.Element* Invert (Modular_double F, size_t order,
                                    Modular_double.Element* A, size_t A_stride, int nullity)

    Modular_double.Element Det (Modular_double F,
                                size_t nrows, size_t ncols,
                                Modular_double.Element* A, size_t A_stride)

    int Rank (Modular_double,
              size_t nrows, size_t ncols,
              Modular_double.Element *A, size_t lda)

    size_t ReducedRowEchelonForm (Modular_double F, size_t a, size_t b,
                                  Modular_double.Element* matrix,
                                  size_t s, size_t* P, size_t* Q)

    void applyP (Modular_double F,
                 FFLAS_SIDE s, FFLAS_TRANSPOSE tr,
                 size_t nr, size_t foo, size_t r,
                 Modular_double.Element* matrix, size_t nc, size_t* Q)

    void MinPoly ( Modular_double& F,
                   vector[Modular_double.Element] minP, size_t N,
                   Modular_double.Element*A, size_t lda)

    void CharPoly ( PolynomialRing_Modular_double& R,
                    ModDoubleDensePolynomial& charp, size_t N,
                    Modular_double.Element* A, size_t lda)

    # float

    bint IsSingular (Modular_float F,
                     size_t nrows, size_t ncols, Modular_float.Element* A,
                     size_t A_stride)

    Modular_float.Element* Invert (Modular_float F, size_t order,
                                   Modular_float.Element* A, size_t A_stride, int nullity)

    Modular_float.Element Det (Modular_float F,
                               size_t nrows, size_t ncols,
                               Modular_float.Element* A, size_t A_stride)

    int Rank (Modular_float,
              size_t nrows, size_t ncols,
              Modular_float.Element *A, size_t lda)

    size_t ReducedRowEchelonForm (Modular_float F, size_t a, size_t b,
                                  Modular_float.Element* matrix,
                                  size_t s, size_t* P, size_t* Q)

    void applyP (Modular_float F,
                 FFLAS_SIDE s, FFLAS_TRANSPOSE tr,
                 size_t nr, size_t foo, size_t r,
                 Modular_float.Element* matrix, size_t nc, size_t* Q)

    void MinPoly ( Modular_float F,
                   vector[Modular_float.Element] minP, size_t N,
                   Modular_float.Element* A, size_t lda)

    void CharPoly ( PolynomialRing_Modular_float& F,
                    ModFloatDensePolynomial& charp, size_t N,
                    Modular_float.Element* A, size_t lda )
