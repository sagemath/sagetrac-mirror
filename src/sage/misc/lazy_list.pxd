from cpython.object cimport *

cdef class lazy_list(object):
    cdef list cache
    cdef Py_ssize_t start, stop, step

    cdef int update_cache_up_to(self, Py_ssize_t i) except *

cdef class lazy_list_from_iterator(lazy_list):
    cdef object iterator

    cdef int update_cache_up_to(self, Py_ssize_t i) except -1

cdef class lazy_list_from_fun(lazy_list):
    cdef object fun

    cdef int update_cache_up_to(self, Py_ssize_t i) except -1
    
cdef class lazy_list_iterator(object):
    cdef lazy_list l
    cdef Py_ssize_t pos, step

cdef class stopped_lazy_list_iterator(object):
    cdef lazy_list l
    cdef Py_ssize_t pos, step, stop


