cimport cython
from libc.stdint cimport uint64_t
from sage.ext.memory_allocator cimport MemoryAllocator

cdef int char_from_tuple(tuple tup, uint64_t *output,
                         size_t face_length) except 0

cdef int char_from_incidence(tuple incidences, uint64_t *output,
                             size_t face_length) except 0

cdef size_t vertex_repr_from_bitrep(uint64_t *face, size_t *output,
                                    size_t face_length) except? 0


cdef ListOfFaces get_facets_from_incidence_matrix(tuple matrix)

cdef ListOfFaces get_vertices_from_incidence_matrix(tuple matrix)

cdef ListOfFaces get_facets_bitrep_from_facets_tuple(tuple facets_input,
                                                     size_t nr_vertices)

cdef ListOfFaces get_vertices_bitrep_from_facets_tuple(tuple facets_input,
                                                       size_t nr_vertices)

@cython.final
cdef class ListOfFaces:
    cdef uint64_t **data
    cdef MemoryAllocator _mem
    cdef size_t nr_faces, face_length, nr_vertices
    cdef int copy_face(self, size_t index, uint64_t *face) except 0

cdef int calculate_dimension(ListOfFaces faces) except -2

@cython.final
cdef class FaceIterator:
    cdef uint64_t *face
    cdef int current_dimension, dimension, record_dimension, lowest_dimension
    cdef MemoryAllocator _mem
    cdef uint64_t ***newfaces2
    cdef tuple newfaces_lists
    cdef uint64_t ***newfaces
    cdef uint64_t **forbidden
    cdef size_t *nr_faces
    cdef size_t *nr_forbidden
    cdef int *first_time
    cdef size_t yet_to_yield, face_length, nr_facets, nr_vertices
    cdef size_t *output1
    cdef size_t *output2
    cdef int nr_lines

    cdef void set_record_dimension(self, int dim)
    cdef inline int next_face(self) except -1
    cdef inline int next_face_loop(self) except -1
    cdef size_t length_vertex_repr(self) except? 0
    cdef size_t facet_repr(self, size_t *output) except? 0
    cdef size_t vertex_repr(self, size_t *output) except? 0
    cdef size_t *get_output1_array(self) except NULL
    cdef size_t *get_output2_array(self) except NULL
