# distutils: libraries = gmp flint arb
# distutils: depends = bernoulli.h

from ..flint.types cimport fmpq_t, ulong

# bernoulli.h
cdef extern from "arb_wrap.h":
    void bernoulli_fmpq_ui(fmpq_t b, ulong n)
