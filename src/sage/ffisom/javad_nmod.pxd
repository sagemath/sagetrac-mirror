from sage.libs.flint.nmod_poly cimport *

cdef extern from "javad_nmod/ff_isom.h":
    cdef cppclass FFIsomorphism:
        FFIsomorphism(nmod_poly_t, nmod_poly_t) except +
        nmod_poly_t x_image

cdef extern from "javad_nmod/ff_isom_base_change.h":
    cdef cppclass FFIsomBaseChange:
        FFIsomBaseChange() except +
        void change_basis(nmod_poly_t result, const nmod_poly_t f, const nmod_poly_t g, const nmod_poly_t modulus)
