cdef class _LazyString():
    cdef func
    cdef args
    cdef kwargs
    cdef val(self)
    cdef update_lazy_string(self, args, kwds)
