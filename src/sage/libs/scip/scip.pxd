from sage.numerical.backends.generic_backend cimport GenericBackend
from decl cimport SCIP_c

cdef class SCIP(GenericBackend):
    cdef SCIP_c *_scip
    cdef list _variables
    cdef list _constraints

