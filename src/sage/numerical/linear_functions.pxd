from sage.structure.parent cimport Parent, Parent_richcmp_element_without_coercion
from sage.structure.element cimport ModuleElement, RingElement, Element

cdef is_LinearFunction(x)

cdef class LinearFunctionOrConstraint(ModuleElement):
    pass

cdef class LinearFunctionsParent_class(Parent):
    cdef _element_constructor_(self, x)
    cdef _coerce_map_from_(self, R)
    cdef public _multiplication_symbol

cdef class LinearFunction(LinearFunctionOrConstraint):
    cdef dict _f
    cdef _add_(self, other)
    cdef iteritems(self)
    cdef _acted_upon_(self, x, bint self_on_left)
    cdef is_zero(self)
    cdef equals(LinearFunction left, LinearFunction right)

cdef class LinearConstraintsParent_class(Parent):
    cdef LinearFunctionsParent_class _LF
    cdef _element_constructor_(self, left, right=?, equality=?)
    cdef _coerce_map_from_(self, R)

cdef class LinearConstraint(LinearFunctionOrConstraint):
    cdef bint equality
    cdef list constraints
    cdef equals(LinearConstraint left, LinearConstraint right)
