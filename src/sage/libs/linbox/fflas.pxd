# distutils: extra_compile_args = FFLASFFPACK_CFLAGS
# distutils: libraries = FFLASFFPACK_LIBRARIES
# distutils: library_dirs = FFLASFFPACK_LIBDIR
# distutils: language = c++ 

from .givaro cimport (Modular_double as ModDoubleField,
        Modular_float as ModFloatField, Dense, Sparse)
from .givaro cimport givvector, Poly1Dom
from libcpp.vector cimport vector


ctypedef Poly1Dom[ModDoubleField, Dense] PolynomialRing_Modular_double
ctypedef Poly1Dom[ModFloatField, Dense] PolynomialRing_Modular_float

# bug in Cython
# Modular.Element can not be used to instantiate templates
# givvector below
cdef extern from *:
    cdef cppclass ModDoubleFieldElement "Givaro::Modular<double>::Element":
        pass
    cdef cppclass ModFloatFieldElement "Givaro::Modular<float>::Element":
        pass

ctypedef givvector[ModDoubleFieldElement] ModDoubleDensePolynomial
ctypedef givvector[ModFloatFieldElement] ModFloatDensePolynomial

cdef extern from "fflas-ffpack/fflas-ffpack.h" namespace "FFLAS":
    ctypedef enum FFLAS_TRANSPOSE:
        FflasNoTrans
        FflasTrans

    ctypedef enum FFLAS_SIDE:
        FflasRight

    # double
    void fgemv \
            (ModDoubleField F, FFLAS_TRANSPOSE transA,
             size_t nrows, size_t ncols,
             ModDoubleField.Element alpha, ModDoubleField.Element* A,
             size_t lda, ModDoubleField.Element* X, size_t incX,
             ModDoubleField.Element beta, ModDoubleField.Element* Y,
             size_t incY)

    ModDoubleField.Element* fgemm \
            (ModDoubleField F,
             FFLAS_TRANSPOSE transA, FFLAS_TRANSPOSE transB,
             size_t nrowsA, size_t ncolsB, size_t ncolsA,
             ModDoubleField.Element alpha, ModDoubleField.Element* A,
             size_t A_stride, ModDoubleField.Element* B, int B_stride,
             ModDoubleField.Element beta, ModDoubleField.Element* C,
             size_t C_stride)


    # float
    void fgemv \
            (ModFloatField F, FFLAS_TRANSPOSE transA,
             size_t nrows, size_t ncols,
             ModFloatField.Element alpha, ModFloatField.Element* A,
             size_t lda, ModFloatField.Element* X, size_t incX,
             ModFloatField.Element beta, ModFloatField.Element* Y,
             size_t incY)

    ModFloatField.Element* fgemm \
            (ModFloatField F,
             FFLAS_TRANSPOSE transA, FFLAS_TRANSPOSE transB,
             size_t nrowsA, size_t ncolsB, size_t ncolsA,
             ModFloatField.Element alpha, ModFloatField.Element* A,
             size_t A_stride, ModFloatField.Element* B, int B_stride,
             ModFloatField.Element beta, ModFloatField.Element* C,
             size_t C_stride)

cdef extern from "fflas-ffpack/fflas-ffpack.h" namespace "FFPACK":
    # double
    bint IsSingular (ModDoubleField F,
                     size_t nrows, size_t ncols, ModDoubleField.Element* A,
                     size_t A_stride)

    ModDoubleField.Element* Invert (ModDoubleField F, size_t order,
                                    ModDoubleField.Element* A, size_t A_stride, int nullity)

    ModDoubleField.Element Det (ModDoubleField F,
                                size_t nrows, size_t ncols,
                                ModDoubleField.Element* A, size_t A_stride)

    int Rank (ModDoubleField,
              size_t nrows, size_t ncols,
              ModDoubleField.Element *A, size_t lda)

    size_t ReducedRowEchelonForm (ModDoubleField F, size_t a, size_t b,
                                  ModDoubleField.Element* matrix,
                                  size_t s, size_t* P, size_t* Q)

    void applyP (ModDoubleField F,
                 FFLAS_SIDE s, FFLAS_TRANSPOSE tr,
                 size_t nr, size_t foo, size_t r,
                 ModDoubleField.Element* matrix, size_t nc, size_t* Q)

    void MinPoly ( ModDoubleField& F,
                   vector[ModDoubleField.Element] minP, size_t N,
                   ModDoubleField.Element*A, size_t lda)

    void CharPoly ( PolynomialRing_Modular_double& R,
                    ModDoubleDensePolynomial& charp, size_t N,
                    ModDoubleField.Element* A, size_t lda)

    # float

    bint IsSingular (ModFloatField F,
                     size_t nrows, size_t ncols, ModFloatField.Element* A,
                     size_t A_stride)

    ModFloatField.Element* Invert (ModFloatField F, size_t order,
                                   ModFloatField.Element* A, size_t A_stride, int nullity)

    ModFloatField.Element Det (ModFloatField F,
                               size_t nrows, size_t ncols,
                               ModFloatField.Element* A, size_t A_stride)

    int Rank (ModFloatField,
              size_t nrows, size_t ncols,
              ModFloatField.Element *A, size_t lda)

    size_t ReducedRowEchelonForm (ModFloatField F, size_t a, size_t b,
                                  ModFloatField.Element* matrix,
                                  size_t s, size_t* P, size_t* Q)

    void applyP (ModFloatField F,
                 FFLAS_SIDE s, FFLAS_TRANSPOSE tr,
                 size_t nr, size_t foo, size_t r,
                 ModFloatField.Element* matrix, size_t nc, size_t* Q)

    void MinPoly ( ModFloatField F,
                   vector[ModFloatField.Element] minP, size_t N,
                   ModFloatField.Element* A, size_t lda)

    void CharPoly ( PolynomialRing_Modular_float& F,
                    ModFloatDensePolynomial& charp, size_t N,
                    ModFloatField.Element* A, size_t lda )
