include "sage/libs/ntl/decl.pxi"

cdef extern from "flint/fmpz.h":

    ctypedef long fmpz
    ctypedef long * fmpz_t
    ctypedef void * mpz_t

    void fmpz_init(fmpz_t x)

    void fmpz_set_ui(fmpz_t res, unsigned long x)
    void fmpz_set_si(fmpz_t res, long x)

    void fmpz_clear(fmpz_t f)
    void fmpz_print(fmpz_t f)
    int fmpz_is_one(fmpz_t f)
    int fmpz_sgn(fmpz_t f)

    void fmpz_get_mpz(mpz_t rop, fmpz_t op)
    void fmpz_set_mpz(fmpz_t rop, mpz_t op)
    void fmpz_neg(fmpz_t rop, fmpz_t op)

    void fmpz_init_set_readonly(fmpz_t rop, mpz_t op)
    void fmpz_clear_readonly(fmpz_t op)
    void flint_mpz_init_set_readonly(mpz_t rop, fmpz_t op)
    void flint_mpz_clear_readonly(mpz_t rop)

    void fmpz_add_ui(fmpz_t f, fmpz_t g, unsigned long c)
