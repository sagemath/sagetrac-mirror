# distutils: libraries = flint

from sage.libs.flint.types cimport fmpz_mod_poly_t, fmpz_t, slong

cdef extern from "flint/fmpz_mod_poly.h":
    void fmpz_mod_poly_init(fmpz_mod_poly_t poly, const fmpz_t p)
    void fmpz_mod_poly_clear(fmpz_mod_poly_t poly)

    void fmpz_mod_poly_set_coeff_fmpz(fmpz_mod_poly_t poly, slong n, const fmpz_t x)
    void fmpz_mod_poly_get_coeff_fmpz(fmpz_t x, const fmpz_mod_poly_t poly, slong n)
