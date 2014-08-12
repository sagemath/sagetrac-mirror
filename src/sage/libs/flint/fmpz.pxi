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

    void fmpz_neg(fmpz_t rop, fmpz_t op)
    void fmpz_get_mpz(mpz_t rop, const fmpz_t op)
    void fmpz_set_mpz(fmpz_t rop, const mpz_t op)

    void fmpz_add_ui(fmpz_t f, fmpz_t g, unsigned long c)

    void fmpz_init_set(fmpz_t f, const fmpz_t g)
    void fmpz_init_set_ui(fmpz_t f, const unsigned long g)

    void fmpz_abs(fmpz_t f1, const fmpz_t f2)
    void fmpz_set(fmpz_t f, const fmpz_t val)
    void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h)
    int fmpz_sgn(const fmpz_t f)
    double fmpz_get_d(const fmpz_t f)

    void fmpz_fdiv_r(fmpz_t f, const fmpz_t g, const fmpz_t h)
    unsigned long fmpz_fdiv_ui(const fmpz_t g, unsigned long x)
    void fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_lcm(fmpz_t f, const fmpz_t g, const fmpz_t h)

    void fmpz_fdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)

    size_t fmpz_sizeinbase(const fmpz_t f, int b)
    char * fmpz_get_str(char * str, int b, const fmpz_t f)
    int fmpz_set_str(fmpz_t f, const char * str, int b)
    int fmpz_cmp(const fmpz_t f, const fmpz_t g)
    int fmpz_cmp_si(const fmpz_t f, int g)
    int fmpz_cmp_ui(const fmpz_t f, unsigned int g)

    void fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h)
