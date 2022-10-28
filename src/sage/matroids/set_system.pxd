from sage.data_structures.bitset cimport bitset_t

cdef class SetSystem:
    cdef long _groundset_size, _bitset_size
    cdef tuple _groundset
    cdef dict _idx
    cdef bitset_t* _subsets
    cdef long _len, _capacity
    cdef bitset_t _temp

    cdef copy(self)
    cdef _relabel(self, l)
    cdef _complements(self)

    cdef resize(self, k=*)
    cdef _append(self, bitset_t X)
    cdef append(self, X)
    cdef _subset(self, long k)
    cdef subset(self, k)
    cdef _get_groundset(self)

    cdef list _incidence_count(self, E)
    cdef SetSystem _groundset_partition(self, SetSystem P, list cnt)
    cdef long subset_characteristic(self, SetSystem P, long e)
    cdef subsets_partition(self, SetSystem P=*, E=*)
    cdef _distinguish(self, Py_ssize_t v)
    cdef is_connected(self)

    cdef initial_partition(self, SetSystem P=*, E=*)
    cdef _equitable_partition(self, SetSystem P=*, EP=*)
    cdef _heuristic_partition(self, SetSystem P=*, EP=*)
    cdef _isomorphism(self, SetSystem other, SetSystem SP=*, SetSystem OP=*)
    cdef _equivalence(self, is_equiv, SetSystem other, SetSystem SP=*, SetSystem OP=*)

cdef class SetSystemIterator:
    cdef SetSystem _H
    cdef long _pointer, _len
