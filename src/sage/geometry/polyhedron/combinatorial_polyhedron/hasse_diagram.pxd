cdef extern from "hasse_diagram.h":
    ctypedef void* CombinatorialPolyhedron_ptr;
    cdef CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(tuple py_tuple, unsigned int nr_vertices, int is_unbounded, unsigned int nr_lines)
    cdef CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(tuple py_tuple, int is_unbounded, unsigned int nr_lines)
    cdef unsigned int dimension(CombinatorialPolyhedron_ptr C)
    cdef void f_vector(CombinatorialPolyhedron_ptr C, unsigned long *vector)
    cdef unsigned int ** edges(CombinatorialPolyhedron_ptr C)
    cdef unsigned int ** ridges(CombinatorialPolyhedron_ptr C)
    cdef tuple incidences(CombinatorialPolyhedron_ptr C, int dimension_one, int dimension_two)
    cdef void record_all_faces(CombinatorialPolyhedron_ptr C)
    cdef tuple get_faces(CombinatorialPolyhedron_ptr C, int dimension, unsigned int facet_repr)
    cdef unsigned long  get_flag(CombinatorialPolyhedron_ptr C, unsigned int *flagarray, unsigned int length)
    cdef void delete_CombinatorialPolyhedron(CombinatorialPolyhedron_ptr)
    cdef unsigned long get_maxnumberedges()
