from object cimport SCIPObject
from decl cimport SCIP_CONS

cdef class Constraint(SCIPObject):
    cdef SCIP_CONS *_cons
    cdef object _name

cdef class LinearConstraint(Constraint):
    pass

cdef class QuadraticConstraint(Constraint):
    pass

cdef class XORConstraint(Constraint):
    pass


cdef class ANDConstraint(Constraint):
    pass

cdef class ORConstraint(Constraint):
    pass

cdef class LogicORConstraint(Constraint):
    pass

cdef class SetPPCConstraint(Constraint):
    pass
