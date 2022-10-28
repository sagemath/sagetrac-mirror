from .matroid cimport Matroid

cdef class CircuitClosuresMatroid(Matroid):
    cdef frozenset _groundset  # _E
    cdef dict _circuit_closures  # _CC
    cdef int _matroid_rank  # _R
    cdef groundset(self)
    cdef _rank(self, X)
    cdef full_rank(self)
    cdef _is_independent(self, F)
    cdef _max_independent(self, F)
    cdef _circuit(self, F)
    cdef circuit_closures(self)
    cdef _is_isomorphic(self, other, certificate=*)
