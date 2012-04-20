from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport ModuleElement
from sage.categories.action cimport Action
from sage.rings.padics.pow_computer cimport PowComputer_long

#cdef extern from "../../../ext/multi_modular.h":
#    ctypedef unsigned long mod_int
#    mod_int MOD_INT_MAX



cdef class Dist(ModuleElement):
    cpdef normalize(self)

cdef class Dist_vector(Dist):
    cdef public moments
    cdef Dist_vector _new_c(self)

#cdef class Dist2(Dist): # only works on 64-bit....
#    cdef long[60] moments
#    cdef int prec
#    cdef public PowComputer_long prime_pow
#    cdef Dist2 _new_c(self)

cdef class Dist_long(Dist):
    cdef long[60] moments # 38 once 2 is special-cased
    cdef int prec
    cdef public PowComputer_long prime_pow
    cdef int quasi_normalize(self) except -1
    cdef Dist_long _new_c(self)

cdef class WeightKAction(Action):
    cdef public _k
    cdef public _character
    cdef public _tuplegen
    cdef public _p
    cdef public _Np
    cdef public _actmat
    cdef public _maxprecs

    cpdef _check_mat(self, a, b, c, d)
    cpdef acting_matrix(self, g, M)
    cpdef _compute_acting_matrix(self, g, M)

cdef class WeightKAction_vector(WeightKAction):
    pass

cdef class SimpleMat(SageObject):
    cdef long* _mat
    cdef long M
    cdef bint _inited

cdef class WeightKAction_long(WeightKAction):
    pass

cdef class iScale(Action):
    pass
