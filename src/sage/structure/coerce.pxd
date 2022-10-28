from .parent cimport Parent
from .coerce_dict cimport TripleDict

cdef py_scalar_parent(py_type)
cdef py_scalar_to_element(py)
cdef bint parent_is_integers(P) except -1
cdef bint is_numpy_type(t)
cdef bint is_mpmath_type(t)


cdef class CoercionModel:
    # This MUST be a mapping to tuples, where each
    # tuple contains at least two elements that are either
    # None or of type Morphism.
    cdef readonly TripleDict _coercion_maps

    # This MUST be a mapping to actions.
    cdef readonly TripleDict _action_maps

    cdef canonical_coercion(self, x, y)
    cdef bin_op(self, x, y, op)
    cdef richcmp(self, x, y, int op)

    cdef coercion_maps(self, R, S)
    cdef discover_coercion(self, R, S)
    cdef verify_coercion_maps(self, R, S, homs, bint fix=*)
    cdef verify_action(self, action, R, S, op, bint fix=*)

    cdef get_action(self, R, S, op=*, r=*, s=*)
    cdef discover_action(self, R, S, op, r=*, s=*)

    cdef bint _record_exceptions
    cdef _record_exception(self)
    cdef readonly list _exception_stack
    cdef bint _exceptions_cleared

    cdef TripleDict _division_parents
    cdef analyse(self, xp, yp, op=*)
    cdef division_parent(self, Parent P)


# Unique global coercion_model instance
cdef CoercionModel coercion_model
