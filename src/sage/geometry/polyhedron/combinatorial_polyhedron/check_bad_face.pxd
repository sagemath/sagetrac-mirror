# distutils: language = c++
# distutils: libraries = normaliz


cimport cython
from libc.stdint                cimport uint64_t, uint8_t

cdef check_from_file2(size_t m, size_t n_coatoms, size_t **PolyIneq, path, size_t n_threads)

cdef extern from "check_bad_face.cc":
    cdef int check_bad_face(size_t **PolyIneq, size_t n_coatoms, size_t m, uint64_t LHS, size_t *Hrep, size_t n_Hrep, size_t e) nogil
    cdef int check_bad_face_uint8_t "check_bad_face" (size_t **PolyIneq, size_t n_coatoms, size_t m, uint64_t LHS, uint8_t *Hrep, size_t n_Hrep, size_t e) nogil
    cdef int check_bad_faces(size_t **PolyIneq, size_t n_coatoms, size_t m,
            uint64_t *LHS, size_t n_bad_faces, uint8_t **bad_faces) nogil
