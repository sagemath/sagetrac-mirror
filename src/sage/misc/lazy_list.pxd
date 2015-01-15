from cpython.object cimport *

cdef class lazy_list(object):
    cdef Py_ssize_t start, stop, step

cdef class lazy_list_with_cache(lazy_list):
    cdef list cache

    cdef int update_cache_up_to(self, Py_ssize_t i) except *
    
cdef class lazy_list_from_iterator(lazy_list_with_cache):
    cdef object iterator

    cdef int update_cache_up_to(self, Py_ssize_t i) except -1

cdef class lazy_list_from_function(lazy_list_with_cache):
    cdef object function

    cdef int update_cache_up_to(self, Py_ssize_t i) except -1

cdef class lazy_list_explicit(lazy_list_with_cache):
    cdef object function

    cdef int update_cache_up_to(self, Py_ssize_t i) except -1

cdef class lazy_list_periodic(lazy_list):
    cdef list pre_period, period

cdef class lazy_list_iterator(object):
    cdef lazy_list l
    cdef Py_ssize_t pos, step

cdef class lazy_list_with_cache_iterator(object):
    cdef lazy_list_with_cache l
    cdef Py_ssize_t pos, step

cdef class stopped_lazy_list_with_cache_iterator(lazy_list_with_cache_iterator):
    cdef Py_ssize_t stop

cdef class lazy_list_explicit_iterator(object):
    cdef lazy_list_explicit l
    cdef Py_ssize_t pos, step
    
cdef class stopped_lazy_list_explicit_iterator(lazy_list_explicit_iterator):
    cdef Py_ssize_t stop
    
cdef class lazy_list_periodic_iterator(object):
    cdef lazy_list_periodic l
    cdef Py_ssize_t pos, step

cdef class stopped_lazy_list_periodic_iterator(lazy_list_periodic_iterator):
    cdef Py_ssize_t stop

