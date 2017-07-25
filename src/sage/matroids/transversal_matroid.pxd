from sage.data_structures.bitset cimport bitset_t

from .matroid cimport Matroid
from .basis_exchange_matroid cimport BasisExchangeMatroid

cdef class TransversalMatroid(BasisExchangeMatroid):
    cdef list _matching, _buckets
    cdef object _D
