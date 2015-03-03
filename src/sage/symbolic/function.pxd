from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport Element

include "sage/ext/stdsage.pxi"

# In many applications, such as plotting, these functions are called many times
# repeatedly. This method is slightly faster than sage.structure.coerce.parent
# The only difference is the PyNumber_Check clause.
from cpython.number cimport *
cdef inline parent_c(x):
    if PY_TYPE_CHECK(x, Element):
        return (<Element>x)._parent
    elif PyNumber_Check(x):
        return <object>PY_TYPE(x)
    elif hasattr(x, 'parent'):
        return x.parent()
    return <object>PY_TYPE(x)

cdef class Function(SageObject):
    cdef unsigned int _serial
    cdef int _nargs
    cdef object _name
    cdef object _alt_name
    cdef object _latex_name
    cdef object _conversions
    cdef object _evalf_params_first
    cdef _is_registered(self)
    cdef _register_function(self)

cdef class BuiltinFunction(Function):
    cdef _is_registered(self)

cdef class GinacFunction(BuiltinFunction):
    cdef object _ginac_name
    cdef _is_registered(self)
    cdef _register_function(self)

cdef class SymbolicFunction(Function):
    # cache hash value
    cdef long _hash_(self) except -1
    cdef bint __hinit
    cdef long __hcache
    cdef _is_registered(self)
