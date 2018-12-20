cdef extern from "hasse_diagram.h":
    ctypedef void* CombinatorialType_ptr;
    cdef CombinatorialType_ptr init_CombinatorialType(tuple py_tuple, unsigned int nr_vertices)
    cdef CombinatorialType_ptr init_CombinatorialType(tuple py_tuple)
    cdef unsigned int dimension(CombinatorialType_ptr C)
    cdef tuple f_vector(CombinatorialType_ptr C)
    cdef tuple edges(CombinatorialType_ptr C)
    cdef tuple ridges(CombinatorialType_ptr C)
    cdef void delete_CombinatorialType(CombinatorialType_ptr)
