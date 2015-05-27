# distutils: extra_compile_args = -std=c++11

from ginac cimport GConstant

cdef class PynacConstant:
    cdef GConstant* pointer
    cdef GConstant object
    cdef object _name
