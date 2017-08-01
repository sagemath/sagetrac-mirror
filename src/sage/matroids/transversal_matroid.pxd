from sage.data_structures.bitset cimport bitset_t

from .matroid cimport Matroid
from .basis_exchange_matroid cimport BasisExchangeMatroid

cpdef graph_from_buckets(buckets, groundset=*)

cdef class TransversalMatroid(BasisExchangeMatroid):
    cdef set _matching
    cdef frozenset _buckets
    cdef object _D

    cpdef transversal_extension(self, element=*, newset=*, sets=*)
    cpdef sets(self)
