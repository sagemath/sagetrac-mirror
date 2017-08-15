from sage.data_structures.bitset cimport bitset_t

from .matroid cimport Matroid
from .basis_exchange_matroid cimport BasisExchangeMatroid

cpdef graph_from_buckets(buckets, groundset)

cdef class TransversalMatroid(BasisExchangeMatroid):
    cdef dict _matching
    cdef object sets
    cdef object _D
    cdef list buckets # temporary

    cpdef sets(self)
    cpdef transversal_extension(self, element=*, newset=*, sets=*)
    cpdef reduce_presentation(self)
