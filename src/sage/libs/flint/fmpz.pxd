include "sage/libs/ntl/decl.pxi"

from sage.libs.gmp.types cimport mpz_t
from flint cimport ulong, slong

cdef extern from "flint/fmpz.h":

    ctypedef slong fmpz
    ctypedef fmpz fmpz_t[1]

    # Memory management
    void fmpz_init(fmpz_t)
    void fmpz_init2(fmpz_t, ulong limbs)

    void fmpz_clear(fmpz_t)

    void fmpz_init_set(fmpz_t, fmpz_t)
    void fmpz_init_set_ui(fmpz_t, ulong)

    # Conversion
    void fmpz_set_ui(fmpz_t, ulong)
    void fmpz_neg_ui(fmpz_t, ulong)
    ulong fmpz_get_ui(fmpz_t)

    void fmpz_set_si(fmpz_t, slong)
    slong fmpz_get_si(fmpz_t)

    void fmpz_set_d(fmpz_t, double)
    double fmpz_get_d(fmpz_t)
    double fmpz_get_d_2exp(slong *, fmpz_t)

    void fmpz_set_mpz(fmpz_t, mpz_t)
    void fmpz_get_mpz(mpz_t, fmpz_t)

    int fmpz_set_str(fmpz_t, char *, int)
    char *fmpz_get_str(char *, int, fmpz_t)

    void fmpz_set_uiui(fmpz_t, mp_limb_t, mp_limb_t)
    void fmpz_neg_uiui(fmpz_t, mp_limb_t, mp_limb_t)
    void fmpz_set_ui_smod(fmpz_t, mp_limb_t, mp_limb_t)

    void flint_mpz_init_set_readonly(mpz_t, fmpz_t)
    void flint_mpz_clear_readonly(mpz_t)

    void fmpz_init_set_readonly(fmpz_t, mpz_t)
    void fmpz_clear_readonly(fmpz_t)

    # Input and output
    int fmpz_read(fmpz_t)
    #int fmpz_fread(FILE *, fmpz_t)
    #size_t fmpz_inp_raw(fmpz_t, FILE *)

    int fmpz_print(fmpz_t)
    #int fmpz_fprint(FILE *, fmpz_t)
    #size_t fmpz_out_raw(FILE *, fmpz_t)

    # Comparison
    int fmpz_cmp(fmpz_t, fmpz_t)
    int fmpz_cmp_ui(fmpz_t, ulong)
    int fmpz_cmp_si(fmpz_t, slong)

    int fmpz_cmpabs(fmpz_t, fmpz_t)

    int fmpz_equal(fmpz_t, fmpz_t)
    int fmpz_equal_ui(fmpz_t, ulong)
    int fmpz_equal_si(fmpz_t, slong)

    int fmpz_is_zero(fmpz_t)
    int fmpz_is_one(fmpz_t)
    int fmpz_is_pm1(fmpz_t)
    int fmpz_is_even(fmpz_t)
    int fmpz_is_odd(fmpz_t)

    # Basic arithmetic
    void fmpz_neg(fmpz_t, fmpz_t)
    void fmpz_abs(fmpz_t, fmpz_t)

    void fmpz_add(fmpz_t, fmpz_t, fmpz_t)
    void fmpz_add_ui(fmpz_t, fmpz_t, ulong)

    void fmpz_sub(fmpz_t, fmpz_t, fmpz_t)
    void fmpz_sub_ui(fmpz_t, fmpz_t, ulong)

    void fmpz_mul(fmpz_t, fmpz_t, fmpz_t)
    void fmpz_mul_ui(fmpz_t, fmpz_t, ulong)
    void fmpz_mul_si(fmpz_t, fmpz_t, slong)
    void fmpz_mul2_uiui(fmpz_t, fmpz_t, ulong, ulong)
    void fmpz_mul_2exp(fmpz_t, fmpz_t, ulong)

    void fmpz_addmul(fmpz_t, fmpz_t, fmpz_t)
    void fmpz_addmul_ui(fmpz_t, fmpz_t, ulong)

    void fmpz_submul(fmpz_t, fmpz_t, fmpz_t)
    void fmpz_submul_ui(fmpz_t, fmpz_t, ulong)

    void fmpz_cdiv_q(fmpz_t, fmpz_t, fmpz_t)
    void fmpz_cdiv_q_ui(fmpz_t, fmpz_t, ulong)
    void fmpz_cdiv_q_si(fmpz_t, fmpz_t, slong)
    void fmpz_cdiv_q_2exp(fmpz_t, fmpz_t, ulong)

    void fmpz_fdiv_q(fmpz_t, fmpz_t, fmpz_t)
    void fmpz_fdiv_q_ui(fmpz_t, fmpz_t, ulong)
    void fmpz_fdiv_q_si(fmpz_t, fmpz_t, slong)
    void fmpz_fdiv_q_2exp(fmpz_t, fmpz_t, ulong)

    ulong fmpz_fdiv_ui(fmpz_t, ulong)

    void fmpz_fdiv_r(fmpz_t, fmpz_t, fmpz_t)
    void fmpz_fdiv_r_2exp(fmpz_t, fmpz_t, ulong)

    void fmpz_fdiv_qr(fmpz_t, fmpz_t, fmpz_t, fmpz_t)

    void fmpz_tdiv_q(fmpz_t, fmpz_t, fmpz_t)
    void fmpz_tdiv_q_ui(fmpz_t, fmpz_t, ulong)
    void fmpz_tdiv_q_si(fmpz_t, fmpz_t, slong)
    void fmpz_tdiv_q_2exp(fmpz_t, fmpz_t, ulong)

    ulong fmpz_tdiv_ui(fmpz_t, ulong)

    void fmpz_tdiv_qr(fmpz_t, fmpz_t, fmpz_t, fmpz_t)

    void fmpz_divexact(fmpz_t, fmpz_t, fmpz_t)
    void fmpz_divexact_ui(fmpz_t, fmpz_t, ulong)
    void fmpz_divexact_si(fmpz_t, fmpz_t, slong)
    void fmpz_divexact2_uiui(fmpz_t, fmpz_t, ulong, ulong)

    void fmpz_mul_tdiv_q_2exp(fmpz_t, fmpz_t, fmpz_t, ulong)
    void fmpz_mul_si_tdiv_q_2exp(fmpz_t, fmpz_t, slong, ulong)

    int fmpz_divisible(fmpz_t, fmpz_t)
    int fmpz_divisible_si(fmpz_t, slong)

    void fmpz_mod(fmpz_t, fmpz_t, fmpz_t)
    void fmpz_mod_ui(fmpz_t, fmpz_t, ulong)

    #void fmpz_preinvn_init(fmpz_preinvn_t, fmpz_t)
    #void fmpz_preinvn_clear(fmpz_preinvn_t)

    #void fmpz_fdiv_preinvn(fmpz_t, fmpz_t, fmpz_t, fmpz_t, fmpz_preinvn_t)

    void fmpz_pow_ui(fmpz_t, fmpz_t, ulong)

    void fmpz_powm(fmpz_t, fmpz_t, fmpz_t, fmpz_t)
    void fmpz_powm_ui(fmpz_t, fmpz_t, ulong, ulong)

    slong fmpz_clog(fmpz_t, fmpz_t)
    slong fmpz_clog_ui(fmpz_t, ulong)

    slong fmpz_flog(fmpz_t, fmpz_t)
    slong fmpz_flog_ui(fmpz_t, ulong)

    double fmpz_dlog(fmpz_t)

    int fmpz_sqrtmod(fmpz_t, fmpz_t, fmpz_t)

    int fmpz_is_square(fmpz_t)
    void fmpz_sqrt(fmpz_t, fmpz_t)
    void fmpz_sqrtrem(fmpz_t, fmpz_t, fmpz_t)

    void fmpz_root(fmpz_t, fmpz_t, slong)

    void fmpz_fac_ui(fmpz_t, ulong)
    void fmpz_fib_ui(fmpz_t, ulong)
    void fmpz_bin_uiui(fmpz_t, ulong, ulong)
    void fmpz_rfac_ui(fmpz_t, fmpz_t, ulong)
    void fmpz_rfac_uiui(fmpz_t, ulong, ulong)
