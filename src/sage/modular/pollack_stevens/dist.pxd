from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport ModuleElement
from sage.categories.action cimport Action
from sage.rings.padics.pow_computer cimport PowComputer_class


cdef class Dist(ModuleElement):
    cdef normalize(self, include_zeroth_moment=*)
    cdef long ordp
    cdef long _ord_p(self)
    cdef long _relprec(self)
    cdef _unscaled_moment(self, long i)

cdef class Dist_vector(Dist):
    cdef public _moments
    cdef Dist_vector _new_c(self)
    cdef Dist_vector _addsub(self, Dist_vector right, bint negate)
    cdef _add_(self, other)


cdef class WeightKAction(Action):
    cdef public _k
    cdef public _character
    cdef public _adjuster
    cdef public _p
    cdef public _Np
    cdef public _actmat
    cdef public _maxprecs
    cdef public _symk
    cdef public _dettwist
    cdef public _Sigma0


    cdef acting_matrix(self, g, M)
    cdef _compute_acting_matrix(self, g, M)

cdef class WeightKAction_vector(WeightKAction):
    pass
