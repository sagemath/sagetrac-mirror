from sage.libs.flint.fmpz_mod_poly cimport *

cdef extern from "javad_fmpz_mod/ff_isom.h":
    cdef cppclass FFIsomorphism:
        FFIsomorphism(fmpz_mod_poly_t, fmpz_mod_poly_t) except +
        fmpz_mod_poly_t x_image

cdef extern from "javad_fmpz_mod/ff_isom_base_change.h":
    cdef cppclass FFIsomBaseChange:
        FFIsomBaseChange() except +
        void change_basis(fmpz_mod_poly_t result, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const fmpz_mod_poly_t modulus)
