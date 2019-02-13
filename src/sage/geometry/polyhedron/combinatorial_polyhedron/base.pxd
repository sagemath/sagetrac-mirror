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
    cdef uint64_t ** data
    cdef MemoryAllocator _mem
    cdef size_t nr_faces, face_length, nr_vertices
    cdef int copy_face(self, size_t index, uint64_t * face) except 0
