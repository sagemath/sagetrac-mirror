cdef class lazy_list_generic(object):
    cdef list cache                    # the cache
    cdef lazy_list_generic master      # a reference if self is a slice or other sublist
    cdef Py_ssize_t start, stop, step  # for slicing
    cdef Py_ssize_t _start_master      # to track changes in master

    cdef str name, separator, more     # used for representation
    cdef str opening_delimiter, closing_delimiter
    cdef int preview                   # number of elements to show in repr
    cdef object cls                    # the class when making a copy
    cdef dict cls_kwds                 # additional keyword arguments

    cpdef get(self, Py_ssize_t i)
    cpdef int _fit(self, Py_ssize_t n) except -1
    cpdef int update_cache_up_to(self, Py_ssize_t i) except -1
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
    cdef Py_ssize_t to_test

cdef class lazy_list_dropwhile(lazy_list_generic):
    cdef object predicate
    cdef bint dropping
