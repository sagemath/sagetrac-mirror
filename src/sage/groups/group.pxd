from sage.monoids.monoid cimport Monoid

cdef class Group(Monoid):
    pass

cdef class AbelianGroup(Group):
    pass

cdef class FiniteGroup(Group):
    pass

cdef class AlgebraicGroup(Group):
    pass

