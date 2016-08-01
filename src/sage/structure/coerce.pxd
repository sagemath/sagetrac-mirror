from .element cimport CoercionModel
from .parent cimport Parent
from .coerce_dict cimport TripleDict

cpdef py_scalar_parent(py_type)
cpdef py_scalar_to_element(py)
cpdef bint is_numpy_type(t)

cdef class CoercionModel_cache_maps(CoercionModel):
    # This MUST be a mapping to tuples, where each
    # tuple contains at least two elements that are either
    # None or of type Morphism.
    cdef TripleDict _coercion_maps

    # This MUST be a mapping to actions.
    cdef TripleDict _action_maps

    cpdef coercion_maps(self, R, S)
    cpdef discover_coercion(self, R, S)
    cpdef verify_coercion_maps(self, R, S, homs, bint fix=*)
    cpdef verify_action(self, action, R, S, op, bint fix=*)

    cpdef get_action(self, R, S, op, r=*, s=*)
    cpdef discover_action(self, R, S, op, r=*, s=*)

    cdef bint _record_exceptions
    cpdef _record_exception(self)
    cdef readonly list _exception_stack
    cdef bint _exceptions_cleared

    cdef TripleDict _division_parents
    cpdef analyse(self, xp, yp, op=*)
    cpdef Parent division_parent(self, Parent parent)
