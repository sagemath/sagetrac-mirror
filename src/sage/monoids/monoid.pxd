from sage.structure.parent cimport Parent

cdef class Monoid(Parent):
    pass

cdef class AbelianMonoid(Monoid):
    pass

cdef class FiniteMonoid(Monoid):
    pass

