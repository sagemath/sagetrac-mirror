from sage.structure.element cimport Element, ModuleElement

cdef class IndexedFreeModuleElement(ModuleElement):
    cdef public dict _monomial_coefficients
    cdef long _hash
    cdef bint _hash_set

    cdef _add_(self, other)
    cdef _sub_(self, other)
    cdef _neg_(self)

    cdef dict monomial_coefficients(self, bint copy=*)
