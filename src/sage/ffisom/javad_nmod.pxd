from sage.libs.flint.nmod_poly cimport *

cdef extern from "javad_nmod/ff_isom.h":
    cdef cppclass FFIsomorphism:
        FFIsomorphism(nmod_poly_t, nmod_poly_t) except +
