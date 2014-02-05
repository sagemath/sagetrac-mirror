from decl cimport SCIP_VAR
from object cimport SCIPObject

cdef class Variable(SCIPObject):
    cdef SCIP_VAR *_var
