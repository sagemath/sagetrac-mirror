from cpython.object cimport *

cdef class abstract_lazy_list(object):
    cdef Py_ssize_t start, stop, step

cdef class lazy_list_with_cache(abstract_lazy_list):
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

cdef class lazy_list_periodic(abstract_lazy_list):
    cdef list pre_period, period
