from sage.libs.gmp.types cimport mp_limb_t, mp_limb_signed_t

cdef extern from "flint/flint.h":
    cdef long FLINT_BITS
    cdef long FLINT_D_BITS

    cdef unsigned long FLINT_BIT_COUNT(unsigned long)

cdef extern from "flint/fmpz.h":
    void _fmpz_cleanup()
    void _fmpz_cleanup_mpz_content()
