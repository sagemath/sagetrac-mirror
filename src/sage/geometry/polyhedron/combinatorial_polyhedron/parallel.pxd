# distutils: language = c++
# distutils: extra_compile_args = "-march=native"
from .face_iterator cimport iter_struct

cdef int parallel_f_vector(iter_struct **face_iter, size_t *f_vector) nogil except -1

#cdef int partial_f(iter_struct *face_iter, size_t *f_vector, size_t i) nogil except -1
