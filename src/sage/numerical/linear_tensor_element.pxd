from sage.structure.element cimport Element, ModuleElement

cdef class LinearTensor(ModuleElement):
    cdef dict _f
    cdef _add_(self, other)
