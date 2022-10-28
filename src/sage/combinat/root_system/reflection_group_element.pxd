from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement

cdef class ComplexReflectionGroupElement(PermutationGroupElement):
    cdef action(self, vec, on_space=*)
    cdef action_on_root_indices(self, i)

cdef class RealReflectionGroupElement(ComplexReflectionGroupElement):
    cdef bint has_left_descent(self, i)
    cdef bint has_descent(self, i, side=*, positive=*)
    cdef action(self, vec, side=*, on_space=*)
    cdef action_on_root_indices(self, i, side=*)
