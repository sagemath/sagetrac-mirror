
from sage.structure.element cimport Element

cdef class ElementWrapper(Element):
    cdef public object value

    cdef _richcmp_(left, right, int op)
    cdef bint _lt_by_value(self, other)

cdef class ElementWrapperCheckWrappedClass(ElementWrapper):
    pass

