from sage.libs.flint.fmpz cimport fmpz_t
from sage.libs.flint.nmod_poly cimport *

cdef extern from "flint/fq_nmod.h":
    ctypedef struct fq_nmod_ctx_struct:
        nmod_poly_t modulus
    ctypedef fq_nmod_ctx_struct fq_nmod_ctx_t[1]
    ctypedef nmod_poly_struct fq_nmod_struct
    ctypedef nmod_poly_t fq_nmod_t

    void fq_nmod_ctx_init(fq_nmod_ctx_t ctx, const fmpz_t p, long d, const char *var)
    void fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, const fmpz_t p, long d, const char *var)
    void fq_nmod_ctx_init_modulus(fq_nmod_ctx_t ctx,
                         nmod_poly_t modulus,
                         const char *var)

    void fq_nmod_ctx_clear(fq_nmod_ctx_t ctx)

    fmpz_t fq_nmod_ctx_prime(fq_nmod_ctx_t ctx)
    long fq_nmod_ctx_degree(const fq_nmod_ctx_t ctx)
    void fq_nmod_ctx_order(fmpz_t f, const fq_nmod_ctx_t ctx)

    void fq_nmod_ctx_print(const fq_nmod_ctx_t ctx)

    #  Memory managment 

    void fq_nmod_init(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
    void fq_nmod_clear(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
    void fq_nmod_reduce(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

    #  Basic arithmetic 

    void fq_nmod_add(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
    void fq_nmod_sub(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
    void fq_nmod_sub_one(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)
    void fq_nmod_neg(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)
    void fq_nmod_mul(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
    void fq_nmod_mul_fmpz(fq_nmod_t rop, const fq_nmod_t op, const fmpz_t x, const fq_nmod_ctx_t ctx)
    void fq_nmod_mul_si(fq_nmod_t rop, const fq_nmod_t op, long x, const fq_nmod_ctx_t ctx)
    void fq_nmod_mul_ui(fq_nmod_t rop, const fq_nmod_t op, unsigned long x, const fq_nmod_ctx_t ctx)
    void fq_nmod_sqr(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    void fq_nmod_inv(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)
    void fq_nmod_gcdinv(fq_nmod_t rop, fq_nmod_t inv, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    void fq_nmod_pow(fq_nmod_t rop, const fq_nmod_t op1, const fmpz_t e, const fq_nmod_ctx_t ctx)
    void fq_nmod_pow_ui(fq_nmod_t rop, const fq_nmod_t op, const unsigned long e, const fq_nmod_ctx_t ctx)
    void fq_nmod_pth_root(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)

    #  Comparison 

    int fq_nmod_equal(const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
    int fq_nmod_is_zero(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    int fq_nmod_is_one(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    #  Assignments and conversions 

    void fq_nmod_set(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    void fq_nmod_set_fmpz(fq_nmod_t rop, const fmpz_t x, const fq_nmod_ctx_t ctx)
    void fq_nmod_set_ui(fq_nmod_t rop, const unsigned long x, const fq_nmod_ctx_t ctx)
    void fq_nmod_set_si(fq_nmod_t rop, const long x, const fq_nmod_ctx_t ctx)
    void fq_nmod_set_coeff_fmpz(fq_nmod_t rop, const fmpz_t x, const unsigned long n, const fq_nmod_ctx_t ctx)
    void fq_nmod_swap(fq_nmod_t op1, fq_nmod_t op2, const fq_nmod_ctx_t ctx)
    void fq_nmod_zero(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
    void fq_nmod_one(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
    void fq_nmod_gen(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

    #  Output 

    void fq_nmod_print(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    int fq_nmod_print_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    char * fq_nmod_get_str(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    char * fq_nmod_get_str_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    #  Special functions 

    void fq_nmod_trace(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    void fq_nmod_frobenius(fq_nmod_t rop, const fq_nmod_t op, long e, const fq_nmod_ctx_t ctx)
    void fq_nmod_norm(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    # Templated functions

    # Not in 2.4 branch
    #void fq_nmod_div(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
