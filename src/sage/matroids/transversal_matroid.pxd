from sage.data_structures.bitset cimport bitset_t

from .matroid cimport Matroid
from .basis_exchange_matroid cimport BasisExchangeMatroid

cpdef graph_from_buckets(groundset, buckets)

cdef class TransversalMatroid(BasisExchangeMatroid):
    cdef list _matching, _buckets
    cdef object _D
