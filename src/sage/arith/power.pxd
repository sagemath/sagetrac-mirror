from cpython.number cimport PyNumber_TrueDivide
from sage.structure.element cimport Element


ctypedef fused ulong_or_object:
    unsigned long
    object


cpdef generic_power(a, n, mul=?, inv=?)
cdef generic_power_long(a, long n, mul, inv)
cdef generic_power_pos(a, ulong_or_object n)  # n > 0
cdef generic_power_pos2(a, ulong_or_object n, mul)



cdef inline invert(a):
    """
    Return ``a^(-1)``.
    """
    if isinstance(a, Element):
        return ~a
    return PyNumber_TrueDivide(type(a)(1), a)


cdef inline one(a):
    """
    Return ``a^0``.
    """
    if isinstance(a, Element):
        return (<Element>a)._parent.one()
    return type(a)(1)
