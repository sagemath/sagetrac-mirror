from sage.data_structures.bitset cimport bitset_t
from .matroid cimport Matroid
from .basis_exchange_matroid cimport BasisExchangeMatroid
from .set_system cimport SetSystem

cdef class BasisMatroid(BasisExchangeMatroid):
    cdef bitset_t _bb
    cdef bitset_t _b
    cdef SetSystem _nonbases
    cdef _bases_invariant_var
    cdef SetSystem _bases_partition_var
    cdef _bases_invariant2_var
    cdef SetSystem _bases_partition2_var
    cdef _bases_invariant3_var
    cdef SetSystem _bases_partition3_var

    cdef reset_current_basis(self)

    cdef _is_basis(self, X)

    cdef bases_count(self)
    cdef bases(self)
    cdef nonbases(self)

    cdef truncation(self)
    cdef _extension(self, e, H)
    cdef _with_coloop(self, e)
    cdef relabel(self, l)

    cdef _bases_invariant(self)
    cdef _bases_partition(self)
    cdef _bases_invariant2(self)
    cdef _bases_partition2(self)
    cdef _bases_invariant3(self)
    cdef _bases_partition3(self)
    cdef _reset_invariants(self)
    cdef  bint is_distinguished(self, e)
    cdef _is_relaxation(self, M, morphism)
    cdef _is_isomorphism(self, M, morphism)
    cdef _isomorphism(self, other)
    cdef _is_isomorphic(self, other, certificate=*)


cdef  binom_init(long n, long k)
cdef  long set_to_index(bitset_t S)
cdef  index_to_set(bitset_t, long, long, long)
