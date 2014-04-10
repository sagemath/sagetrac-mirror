cimport borie
from sage.structure.sage_object cimport SageObject

# cdef class PermList(object)

cpdef Perm one()

cdef class Perm(SageObject):
    cdef borie.SGroup_type _p
    cpdef list as_list(self)

    cpdef Perm mult(self, Perm b)

from libcpp.vector cimport vector

cdef class PermList(object):
    cdef vector[borie.SGroup_type] _v
    cpdef append(self, Perm p)

cdef class PermListList(object):
    cdef vector[vector[borie.SGroup_type]] _v
    cpdef append(self, PermList p)

cpdef bint is_canonical(PermListList a, Perm v)
