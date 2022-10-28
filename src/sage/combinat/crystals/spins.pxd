from sage.structure.element cimport Element

cdef class Spin(Element):
    cdef bint* _value
    cdef int _n
    cdef long _hash

    cdef Spin _new_c(self, bint* value)

cdef class Spin_crystal_type_B_element(Spin):
    cdef Spin e(self, int i)
    cdef Spin f(self, int i)
    cdef int epsilon(self, int i)
    cdef int phi(self, int i)

cdef class Spin_crystal_type_D_element(Spin):
    cdef Spin e(self, int i)
    cdef Spin f(self, int i)
    cdef int epsilon(self, int i)
    cdef int phi(self, int i)

