from sage.data_structures.bitset cimport *
from sage.data_structures.bitset_base cimport bitset_t, bitset_s

from .matroid cimport Matroid
from .set_system cimport SetSystem

cdef class BasisExchangeMatroid(Matroid):
    cdef long _groundset_size, _matroid_rank, _bitset_size
    cdef bitset_t _current_basis, _inside, _outside, _input, _input2, _output, _temp
    cdef tuple _E
    cdef dict _idx
    cdef frozenset _groundset

    cdef _bcount
    cdef _weak_invariant_var, _strong_invariant_var, _heuristic_invariant_var
    cdef SetSystem _weak_partition_var, _strong_partition_var, _heuristic_partition_var

    cdef __relabel(self, l)

    cdef __pack(self, bitset_t, X)
    cdef __unpack(self, bitset_t)
    cdef bint __is_exchange_pair(self, long x, long y) except -1
    cdef int __exchange(self, long x, long y) except -1
    cdef int __move(self, bitset_t X, bitset_t Y) except -1
    cdef __fundamental_cocircuit(self, bitset_t, long x)
    cdef __fundamental_circuit(self, bitset_t, long y)

    cdef __max_independent(self, bitset_t, bitset_t)
    cdef __circuit(self, bitset_t, bitset_t)
    cdef __closure(self, bitset_t, bitset_t)
    cdef __max_coindependent(self, bitset_t, bitset_t)
    cdef __cocircuit(self, bitset_t, bitset_t)
    cdef __coclosure(self, bitset_t, bitset_t)

    cdef __augment(self, bitset_t, bitset_t, bitset_t)
    cdef bint __is_independent(self, bitset_t F) except -1
    cdef __move_current_basis(self, bitset_t, bitset_t)

    cdef bint _set_current_basis(self, F)

    cdef groundset(self)
    cdef groundset_list(self)
    cdef full_rank(self)
    cdef full_corank(self)

    cdef basis(self)
    cdef _move_current_basis(self, X, Y)

    cdef _max_independent(self, F)
    cdef _rank(self, F)
    cdef _circuit(self, F)
    cdef _fundamental_circuit(self, B, e)
    cdef _closure(self, F)

    cdef _max_coindependent(self, F)
    cdef _corank(self, F)
    cdef _cocircuit(self, F)
    cdef _fundamental_cocircuit(self, B, e)
    cdef _coclosure(self, F)

    cdef _augment(self, X, Y)
    cdef _is_independent(self, F)

    cdef f_vector(self)
    cdef  _f_vector_rec(self, object f_vec, bitset_t* flats, bitset_t* todo, long elt, long rnk)
    cdef flats(self, R)
    cdef  _flats_rec(self, SetSystem Rflats, long R, bitset_t* flats, bitset_t* todo, long elt, long rnk)
    cdef coflats(self, R)
    cdef  _coflats_rec(self, SetSystem Rcoflats, long R, bitset_t* coflats, bitset_t* todo, long elt, long cornk)
    cdef _flat_element_inv(self, long k)
    cdef  _flat_element_inv_rec(self, object f_inc, long R, bitset_t* flats, bitset_t* todo, long elt, long i)

    cdef bases_count(self)
    cdef independent_r_sets(self, long r)
    cdef bases(self)
    cdef dependent_r_sets(self, long r)
    cdef nonbases(self)

    cdef nonspanning_circuits(self)
    cdef cocircuits(self)
    cdef circuits(self)

    cdef _characteristic_setsystem(self)
    cdef _weak_invariant(self)
    cdef _weak_partition(self)
    cdef _strong_invariant(self)
    cdef _strong_partition(self)
    cdef _heuristic_invariant(self)
    cdef _heuristic_partition(self)
    cdef _flush(self)

    cdef _equitable_partition(self, P=*)
    cdef _is_isomorphic(self, other, certificate=*)
    cdef _isomorphism(self, other)
    cdef _is_isomorphism(self, other, morphism)
    cdef bint __is_isomorphism(self, BasisExchangeMatroid other, morphism)

    cdef is_valid(self)

cdef bint nxksrd(bitset_s *b, long n, long k, bint succ)
