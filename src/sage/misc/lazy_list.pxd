cdef class lazy_list_generic(object):
    cdef list cache                  # the cache
    cdef lazy_list_generic master   # a reference if self is a slice
    cdef Py_ssize_t start, stop, step

    cdef str name, separator, more, opening_delimiter, closing_delimiter
    cdef int preview
    cdef object cls                    # the class when making a copy

    cpdef get(self, Py_ssize_t i)
    cpdef int _fit(self, Py_ssize_t n) except -1
    cdef int update_cache_up_to(self, Py_ssize_t i) except -1
    cpdef list _get_cache_(self)

cdef class lazy_list_from_iterator(lazy_list_generic):
    cdef object iterator

cdef class lazy_list_from_function(lazy_list_generic):
    cdef object callable

cdef class lazy_list_from_update_function(lazy_list_generic):
    cdef object update_function

cdef class lazy_list_takewhile(lazy_list_generic):
    cdef object predicate
    cdef bint taking

cdef class lazy_list_dropwhile(lazy_list_generic):
    cdef object predicate
    cdef bint dropping
