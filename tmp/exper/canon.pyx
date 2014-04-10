import cython
from sage.rings.integer cimport Integer
from sage.structure.sage_object cimport SageObject

include 'sage/ext/interrupt.pxi'

cimport borie

N = borie.N

cdef class Perm(SageObject):
    def __init__(self, list l):
         for i in range(len(l)):
             self._p[i] = l[i]
         for i in range(len(l), borie.SGroup_vect_len):
             self._p[i] = i

    cpdef list as_list(self):
         cdef list res = []
         cdef unsigned int i
         for i in range(borie.SGroup_N):
             res.append(self._p[i])
         return res

    cpdef Perm mult(self, Perm b):
         cdef Perm res
         res = Perm.__new__(Perm)
         borie.SGroup_mult(res._p, self._p, b._p)
         return res

    def __repr__(self):
        return "".join(str(i) for i in self.as_list())

    def equal(self, Perm other):
        return self._p == other._p

cpdef Perm one():
    cdef Perm res
    res = Perm.__new__(Perm)
    res._p = borie.SGroup_one()
    return res

cdef class PermList(object):
    def __init__(self):
        pass
    cpdef append(self, Perm p):
        self._v.push_back(p._p)
    def __len__(self):
        return self._v.size()
    def __getitem__(self, unsigned int i):
        cdef Perm res
        if i >= self._v.size():
            raise IndexError
        res = Perm.__new__(Perm)
        res._p = self._v.at(i)
        return res

cdef class PermListList(object):
    def __init__(self):
        pass
    cpdef append(self, PermList p):
        self._v.push_back(p._v)
    def __len__(self):
        return self._v.size()
    def __getitem__(self, int i):
        cdef PermList res = PermList.__new__(PermList)
        res._v = self._v.at(i)
        return res


cpdef bint is_canonical(PermListList sgs, Perm v):
    return borie.is_canonical(sgs._v, v._p)

from libcpp.list cimport list as stl_list

cpdef PermList elements_of_depth(int depth, PermListList sgs):
    cdef stl_list[borie.SGroup_type] cppres
    sig_on()
    cppres = borie.elements_of_depth(depth, sgs._v)
    sig_off()
    cdef PermList res = PermList.__new__(PermList)
    cdef vector[borie.SGroup_type] resvect
    cdef borie.SGroup_type tmp
    for tmp in cppres:
        resvect.push_back(tmp)
    res._v = resvect
    return res

cpdef int elements_of_depth_number(int depth, PermListList sgs):
    return borie.elements_of_depth(depth, sgs._v).size()

class CilkError(Exception):
    pass

def cilk_start(int nproc):
    np = bytes(str(nproc))
    if not borie.cilk_start(np):
        raise CilkError, "Unable to set the number of Cilk workers"
