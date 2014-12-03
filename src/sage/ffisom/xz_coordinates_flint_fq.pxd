from sage.libs.flint.fmpz cimport fmpz_t
from sage.libs.flint.fq cimport fq_t, fq_ctx_t

cdef extern from "ellmul/fq_weierstrass_xz.h":
    ctypedef struct fq_weierstrass_xz_struct:
         pass
    ctypedef fq_weierstrass_xz_struct fq_weierstrass_xz_t[1]

    void fq_weierstrass_xz_set_ui(fq_weierstrass_xz_t, unsigned long, unsigned long, fq_ctx_t)
    void fq_weierstrass_xz_set_fq(fq_weierstrass_xz_t, fq_t, fq_t, fq_ctx_t)

    void fq_weierstrass_xz_init(fq_weierstrass_xz_t, fq_ctx_t)
    void fq_weierstrass_xz_clear(fq_weierstrass_xz_t, fq_ctx_t)

    void fq_weierstrass_xz_dbl(fq_t, fq_t, fq_t, fq_t, fq_weierstrass_xz_t, fq_ctx_t)
    void fq_weierstrass_xz_dadd(fq_t, fq_t, fq_t, fq_t, fq_t, fq_t, fq_t, fq_t, fq_weierstrass_xz_t, fq_ctx_t)
    void fq_weierstrass_xz_ladd(fq_t, fq_t, fq_t, fq_t, fq_t, fq_t, fq_t, fq_t, fq_t, fq_t, fq_weierstrass_xz_t, fq_ctx_t)

    void fq_weierstrass_xz_mul_ltr(fq_t, fq_t, fq_t, fq_t, fmpz_t, fq_weierstrass_xz_t, fq_ctx_t)

