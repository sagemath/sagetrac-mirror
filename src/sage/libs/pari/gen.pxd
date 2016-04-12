from .types cimport *
from sage.structure.element cimport RingElement
from .pari_instance cimport PariInstance
cimport cython


cdef class gen_auto(RingElement):
    cdef GEN g
    cdef pari_sp b
    cdef dict refers_to
    cdef PariInstance _pari

@cython.final
cdef class gen(gen_auto):
    pass

cpdef gen objtogen(PariInstance pari, s)
