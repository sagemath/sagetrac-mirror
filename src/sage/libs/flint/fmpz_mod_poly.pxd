from sage.libs.flint.fmpz cimport fmpz_t

cdef extern from "flint/fmpz_mod_poly.h":
    ctypedef struct fmpz_mod_poly_struct:
        pass
    ctypedef fmpz_mod_poly_struct fmpz_mod_poly_t[1]

    void fmpz_mod_poly_init(fmpz_mod_poly_t poly, const fmpz_t p)
    void fmpz_mod_poly_clear(fmpz_mod_poly_t poly)

    void fmpz_mod_poly_set_coeff_fmpz(fmpz_mod_poly_t poly, long n, const fmpz_t x)
    void fmpz_mod_poly_get_coeff_fmpz(fmpz_t x, const fmpz_mod_poly_t poly, long n)
