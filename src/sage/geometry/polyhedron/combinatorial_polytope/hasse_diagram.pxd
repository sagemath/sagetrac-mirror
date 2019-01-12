cdef extern from "hasse_diagram.h":
    ctypedef void* CombinatorialPolytope_ptr;
    cdef CombinatorialPolytope_ptr init_CombinatorialPolytope(tuple py_tuple, unsigned int nr_vertices)
    cdef CombinatorialPolytope_ptr init_CombinatorialPolytope(tuple py_tuple)
    cdef unsigned int dimension(CombinatorialPolytope_ptr C)
    cdef tuple f_vector(CombinatorialPolytope_ptr C)
    cdef tuple edges(CombinatorialPolytope_ptr C)
    cdef tuple ridges(CombinatorialPolytope_ptr C)
    cdef tuple incidences(CombinatorialPolytope_ptr C, int dimension_one, int dimension_two)
    cdef void record_all_faces(CombinatorialPolytope_ptr C)
    cdef tuple get_faces(CombinatorialPolytope_ptr C, int dimension, unsigned int facet_repr)
    cdef unsigned long  get_flag(CombinatorialPolytope_ptr C, tuple  py_tuple)
    cdef void delete_CombinatorialPolytope(CombinatorialPolytope_ptr)
