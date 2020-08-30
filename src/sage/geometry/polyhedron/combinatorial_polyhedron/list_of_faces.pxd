cimport cython
from libc.stdint cimport uint64_t
from sage.ext.memory_allocator cimport MemoryAllocator

@cython.final
cdef class ListOfFaces:
    cdef MemoryAllocator _mem

    # ``face_length`` is the length of each face in terms of ``uint64_t``.
    cdef readonly size_t n_faces, face_length, n_atoms

    # ``data`` points to the raw data.
    # It will be of "type" ``uint64_t[n_faces][face_length]``
    cdef uint64_t **data

    cpdef ListOfFaces __copy__(self)

    cpdef int compute_dimension(self) except -2
    cdef int compute_dimension_loop(self, uint64_t **faces, size_t n_faces,
                                      size_t face_length) except -2

    cpdef ListOfFaces pyramid(self)

    cdef ListOfFaces delete_atoms_unsafe(self, uint64_t *face, int *delete)  # not in place
    cdef void delete_faces_unsafe(self, uint64_t *face, int *delete)  # in place

    cdef void get_not_inclusion_maximal_unsafe(self, int *not_inclusion_maximal)
    cdef void get_faces_all_set_unsafe(self, int *all_set)

cdef tuple face_as_combinatorial_polyhedron(ListOfFaces facets, ListOfFaces Vrep, uint64_t *face, uint64_t *coface)
