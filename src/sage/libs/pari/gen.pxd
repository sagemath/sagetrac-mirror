include 'decl.pxi'

from sage.structure.element cimport RingElement

cdef class gen(RingElement):
    cdef GEN g
    cdef pari_sp b
    cdef dict _refers_to
