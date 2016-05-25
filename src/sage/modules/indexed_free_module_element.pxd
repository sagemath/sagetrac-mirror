from sage.structure.element cimport Element

cdef class IndexedFreeModuleElement(Element):
    cdef public dict _monomial_coefficients
    cdef _hash

    cpdef _richcmp_(left, Element right, int op)
    cpdef int _cmp_(left, Element right) except -2
    cpdef _acted_upon_(self, scalar, bint self_on_left)

cpdef _divide_if_possible(x, y)
