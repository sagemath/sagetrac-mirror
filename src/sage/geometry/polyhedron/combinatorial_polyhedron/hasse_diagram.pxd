cdef extern from "hasse_diagram.h":
    ctypedef void* CombinatorialPolyhedron_ptr;
    cdef CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(tuple py_tuple, unsigned int nr_vertices)
    cdef CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(tuple py_tuple)
    cdef unsigned int dimension(CombinatorialPolyhedron_ptr C)
    cdef tuple f_vector(CombinatorialPolyhedron_ptr C)
    cdef tuple edges(CombinatorialPolyhedron_ptr C)
    cdef tuple ridges(CombinatorialPolyhedron_ptr C)
    cdef tuple incidences(CombinatorialPolyhedron_ptr C, int dimension_one, int dimension_two)
    cdef void record_all_faces(CombinatorialPolyhedron_ptr C)
    cdef tuple get_faces(CombinatorialPolyhedron_ptr C, int dimension, unsigned int facet_repr)
    cdef unsigned long  get_flag(CombinatorialPolyhedron_ptr C, tuple  py_tuple)
    cdef void delete_CombinatorialPolyhedron(CombinatorialPolyhedron_ptr)
