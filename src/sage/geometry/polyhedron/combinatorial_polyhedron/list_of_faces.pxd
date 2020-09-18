cimport cython
from sage.ext.memory_allocator cimport MemoryAllocator
from .face_data_structure   cimport *
from sage.data_structures.alternative_bitset cimport mp_bitcnt_t

cdef struct face_list_s:
    face_t* faces
    size_t n_faces
    size_t total_n_faces
    size_t n_atoms
    size_t n_coatoms
    bint polyhedron_is_simple
    bint* is_not_new_face

ctypedef face_list_s face_list_t[1]

@cython.final
cdef class ListOfFaces:
    cdef MemoryAllocator _mem

    # ``data`` points to the raw data.
    # It will be of "type" ``uint64_t[n_faces][face_length]``
    cdef face_list_t data

    cpdef int compute_dimension(self) except -2

    cdef inline size_t n_faces(self):
        return self.data.n_faces
    cdef inline size_t n_atoms(self):
        return self.data.n_atoms
    cdef inline size_t n_coatoms(self):
        return self.data.n_coatoms

    cpdef ListOfFaces pyramid(self)

cdef size_t get_next_level(
        face_list_t faces,
        face_list_t new_faces,
        face_list_t visited_all) nogil except -1
    # Already raps itself in ``sig_on/sig_off``.

# From list_of_faces.pxi
cdef void sort_faces_list(face_list_t faces)
cdef size_t find_face(face_t face, face_list_t faces)
cdef int face_list_init(face_list_t faces, size_t n_faces, size_t n_atoms, size_t n_coatoms) except -1
cdef int face_list_shallow_init(face_list_t faces, size_t n_faces, size_t n_atoms, size_t n_coatoms, MemoryAllocator mem) except -1
cdef void face_list_free(face_list_t face)
cdef int face_list_shallow_copy(face_list_t dst, face_list_t src) except -1
cdef int add_face_shallow(face_list_t faces, face_t face) nogil except -1
cdef int add_face_deep(face_list_t faces, face_t face) except -1

# From face.pxi
cdef void face_intersection(face_t dest, face_t A, face_t B) nogil
cdef long face_len_atoms(face_t face) nogil
cdef size_t bit_rep_to_coatom_rep(face_t face, face_list_t coatoms, size_t *output)
cdef bint face_init(face_t face, mp_bitcnt_t n_atoms, mp_bitcnt_t n_coatoms) except -1
cdef void face_free(face_t face)
cdef void face_copy(face_t dst, face_t src)
cdef long face_next_atom(face_t face, mp_bitcnt_t n)
cdef int face_add_atom_safe(face_t face, mp_bitcnt_t n) except -1
cdef void face_add_atom(face_t face, mp_bitcnt_t n)
cdef void facet_set_coatom(face_t face, mp_bitcnt_t n)
cdef void face_clear(face_t face)
