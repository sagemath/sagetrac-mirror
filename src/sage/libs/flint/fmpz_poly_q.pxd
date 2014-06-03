from sage.libs.gmp.types cimport mpz_t, mpq_t

from flint cimport *
from fmpz_poly cimport fmpz_poly_struct

cdef extern from "flint/fmpz_poly_q.h":
    ctypedef struct fmpz_poly_q_struct:
        pass

    ctypedef fmpz_poly_q_struct fmpz_poly_q_t[1]

    # Accessing numerator and denominator
    fmpz_poly_struct *fmpz_poly_q_numref(const fmpz_poly_q_t)
    fmpz_poly_struct *fmpz_poly_q_denref(const fmpz_poly_q_t)
    void fmpz_poly_q_canonicalise(fmpz_poly_q_t)
    int fmpz_poly_q_is_canonical(const fmpz_poly_q_t)

    # Memory management
    void fmpz_poly_q_init(fmpz_poly_q_t)
    void fmpz_poly_q_clear(fmpz_poly_q_t)

    # Randomization
    #void fmpz_poly_q_randtest(fmpz_poly_q_t, flint_rand_t, slong, mp_bitcnt_t, slong, mp_bitcnt_t)
    #void fmpz_poly_q_randtest_non_zero(fmpz_poly_q_t, flint_rand_t, slong, mp_bitcnt_t, slong, mp_bitcnt_t)

    # Assignment
    void fmpz_poly_q_set(fmpz_poly_q_t, const fmpz_poly_q_t)
    void fmpz_poly_q_set_si(fmpz_poly_q_t, slong)
    void fmpz_poly_q_swap(fmpz_poly_q_t, fmpz_poly_q_t)
    void fmpz_poly_q_zero(fmpz_poly_q_t)
    void fmpz_poly_q_one(fmpz_poly_q_t)
    void fmpz_poly_q_neg(fmpz_poly_q_t, const fmpz_poly_q_t)
    void fmpz_poly_q_inv(fmpz_poly_q_t, const fmpz_poly_q_t)

    # Comparison
    int fmpz_poly_q_is_zero(const fmpz_poly_q_t)
    int fmpz_poly_q_is_one(const fmpz_poly_q_t)
    int fmpz_poly_q_equal(const fmpz_poly_q_t, const fmpz_poly_q_t)

    # Addition and subtraction
    void fmpz_poly_q_add_in_place(fmpz_poly_q_t, const fmpz_poly_q_t)
    void fmpz_poly_q_sub_in_place(fmpz_poly_q_t, const fmpz_poly_q_t)
    void fmpz_poly_q_add(
            fmpz_poly_q_t, const fmpz_poly_q_t, const fmpz_poly_q_t)
    void fmpz_poly_q_sub(
            fmpz_poly_q_t, const fmpz_poly_q_t, const fmpz_poly_q_t)
    void fmpz_poly_q_addmul(
            fmpz_poly_q_t, const fmpz_poly_q_t, const fmpz_poly_q_t)
    void fmpz_poly_q_submul(
            fmpz_poly_q_t, const fmpz_poly_q_t, const fmpz_poly_q_t)

    # Scalar multiplication and division
    void fmpz_poly_q_scalar_mul_si(fmpz_poly_q_t, const fmpz_poly_q_t, slong)
    void fmpz_poly_q_scalar_mul_mpz(
            fmpz_poly_q_t, const fmpz_poly_q_t, const mpz_t)
    void fmpz_poly_q_scalar_mul_mpq(
            fmpz_poly_q_t, const fmpz_poly_q_t, const mpq_t)
    void fmpz_poly_q_scalar_div_si(fmpz_poly_q_t, const fmpz_poly_q_t, slong)
    void fmpz_poly_q_scalar_div_mpz(
            fmpz_poly_q_t, const fmpz_poly_q_t, const mpz_t)
    void fmpz_poly_q_scalar_div_mpq(
            fmpz_poly_q_t, const fmpz_poly_q_t, const mpq_t)

    # Multiplation and division
    void fmpz_poly_q_mul(
            fmpz_poly_q_t, const fmpz_poly_q_t, const fmpz_poly_q_t)
    void fmpz_poly_q_div(
            fmpz_poly_q_t, const fmpz_poly_q_t, const fmpz_poly_q_t)

    # Powering
    void fmpz_poly_q_pow(fmpz_poly_q_t, const fmpz_poly_q_t, ulong)

    # Derivative
    void fmpz_poly_q_derivative(fmpz_poly_q_t, const fmpz_poly_q_t)

    # Evaluation
    void fmpz_poly_q_evaluate(mpq_t, const fmpz_poly_q_t, const mpq_t)

    # Input and output
    int fmpz_poly_q_set_str(fmpz_poly_q_t, const char *)
    char *fmpz_poly_q_get_str(const fmpz_poly_q_t)
    char *fmpz_poly_q_get_str_pretty(const fmpz_poly_q_t, const char *)
    int fmpz_poly_q_print(const fmpz_poly_q_t)
    int fmpz_poly_q_print_pretty(const fmpz_poly_q_t, const char *)

cdef str fmpz_poly_q_get_python_str(
        const fmpz_poly_q_t func, str name=?, bint latex=?, int base=?)
