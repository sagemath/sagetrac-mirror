# distutils: language = c++
# distutils: libraries = normaliz


cimport cython
from libc.stdint                cimport uint64_t

cdef extern from "check_bad_face.cc":
    cdef int check_bad_face(size_t **PolyIneq, size_t n_coatoms, size_t m, uint64_t LHS, size_t *Hrep, size_t n_Hrep, size_t e) nogil
