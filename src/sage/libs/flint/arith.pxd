from sage.libs.flint.fmpz cimport fmpz_t
from sage.libs.flint.fmpq cimport fmpq_t

cdef extern from "flint/arith.h":
    void arith_number_of_partitions(fmpz_t x, unsigned long n)
    void arith_dedekind_sum(fmpq_t, fmpz_t, fmpz_t)
