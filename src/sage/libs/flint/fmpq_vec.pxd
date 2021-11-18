# distutils: libraries = flint
# distutils: depends = flint/fmpq_vec.h

from sage.libs.flint.types cimport fmpq, slong, ulong, fmpq_t

# flint/fmpq_vec.h
cdef extern from "flint_wrap.h":
    fmpq * _fmpq_vec_init(slong)
    void _fmpq_vec_clear(fmpq *, slong)
    void _fmpq_vec_dot(fmpq_t, const fmpq *, const fmpq *, slong)
