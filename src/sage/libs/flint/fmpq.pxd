from sage.libs.gmp.types cimport mpq_t

cdef extern from "flint/fmpq.h":
    ctypedef struct fmpq:
        pass

    ctypedef fmpq fmpq_t[1]

    # Memory management
    void fmpq_init(fmpq_t)
    void fmpq_clear(fmpq_t)

    # Conversion
    void fmpq_set_mpq(fmpq_t, const mpq_t)
    void fmpq_get_mpq(mpq_t, const fmpq_t)
