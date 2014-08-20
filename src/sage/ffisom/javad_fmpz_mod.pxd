from sage.libs.flint.fmpz_mod_poly cimport *

cdef extern from "javad_fmpz_mod/ff_isom.h":
    cdef cppclass FFIsomorphism:
        FFIsomorphism(fmpz_mod_poly_t, fmpz_mod_poly_t) except +
        fmpz_mod_poly_t x_image
