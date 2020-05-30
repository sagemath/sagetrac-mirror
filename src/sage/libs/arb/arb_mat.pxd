# distutils: depends = arb_mat.h

from sage.libs.arb.types cimport arb_t, arb_mat_t

# arb_mat.h
cdef extern from "arb_wrap.h":
    arb_t arb_mat_entry(arb_mat_t mat, unsigned long i, unsigned long j)
    void arb_mat_init(arb_mat_t mat, long r, long c)
    void arb_mat_clear(arb_mat_t mat)
    bint arb_mat_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)
