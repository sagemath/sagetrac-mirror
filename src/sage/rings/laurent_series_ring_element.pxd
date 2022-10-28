from sage.structure.element cimport AlgebraElement, ModuleElement

cdef class LaurentSeries(AlgebraElement):
    cdef ModuleElement __u
    cdef long __n

    cdef __normalize(self)
    cdef _add_(self, other)
    cdef _mul_(self, other)

