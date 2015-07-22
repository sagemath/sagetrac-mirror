from sage.libs.flint.nmod_poly cimport *

cdef extern from "javad_nmod/ff_isom.h":
    cdef cppclass FFIsomorphism:
        FFIsomorphism(nmod_poly_t, nmod_poly_t) except +
        void compute_generators(nmod_poly_t g1, nmod_poly_t g2)
        void build_isomorphism(nmod_poly_t g1, nmod_poly_t g2)
        void compute_image_using_modcomp(nmod_poly_t image, const nmod_poly_t f)
        void get_x_image(nmod_poly_t x)

cdef extern from "javad_nmod/ff_isom_base_change.h":
    cdef cppclass FFIsomBaseChange:
        FFIsomBaseChange() except +
        void change_basis(nmod_poly_t result, const nmod_poly_t f, const nmod_poly_t g, const nmod_poly_t modulus)
