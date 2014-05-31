### Type Declarations ###

cdef extern from "gmp.h":

    # Underlying typedefs
    ctypedef unsigned long int mp_limb_t
    ctypedef long int mp_limb_signed_t

    ctypedef mp_limb_t mp_bitcnt_t

    ctypedef long int mp_size_t
    ctypedef long int mp_exp_t

    ctypedef mp_limb_t * mp_ptr
    ctypedef const mp_limb_t * mp_srcptr

    # The internal structure is not guaranteed to stay the same with future
    # releases of gmp, so we treat them as a black boxes:
    ctypedef struct __mpz_struct:
        pass
    ctypedef struct __mpq_struct:
        pass
    ctypedef struct __mpf_struct:
        pass
    ctypedef struct __gmp_randstate_struct:
        pass

    # User facing types
    ctypedef __mpz_struct mpz_t[1]
    ctypedef __mpq_struct mpq_t[1]
    ctypedef __mpf_struct mpf_t[1]
    ctypedef __gmp_randstate_struct gmp_randstate_t[1]
