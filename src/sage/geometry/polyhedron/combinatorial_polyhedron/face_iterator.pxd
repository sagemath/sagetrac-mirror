# distutils: language = c++
# distutils: extra_compile_args = "-march=native"
# distutils: sources = bit_vector_operations.cc
cimport cython
from libc.stdint                cimport uint64_t
from sage.ext.memory_allocator  cimport MemoryAllocator
from sage.structure.sage_object cimport SageObject
from .list_of_faces             cimport ListOfFaces
from .combinatorial_face        cimport CombinatorialFace

from .bit_vector_operations               cimport iter_struct


cdef inline int next_face_loop(iter_struct *face_iter) nogil except -1
cdef inline int ignore_subfaces(iter_struct *face_iter) nogil except -1
cdef inline int ignore_supfaces(iter_struct *face_iter) nogil except -1
cdef inline size_t length_atom_repr(iter_struct *face_iter) except -1
cdef inline size_t set_coatom_repr(iter_struct *face_iter) except -1
cdef inline size_t set_atom_repr(iter_struct *face_iter) except -1
cdef inline int next_dimension(iter_struct *face_iter) nogil except -1
cdef inline int prepare_partial_iter(iter_struct *face_iter, size_t i) nogil except -1

@cython.final
cdef class FaceIterator(SageObject):
    cdef readonly bint dual         # if 1, then iterate over dual Polyhedron
    cdef iter_struct structure
    cdef MemoryAllocator _mem
    cdef tuple newfaces_lists       # tuple to hold the ListOfFaces corresponding to maybe_newfaces
    cdef tuple _V, _H, _equalities  # some copies from ``CombinatorialPolyhedron``

    # Atoms and coatoms are the vertices/facets of the Polyedron.
    # If ``dual == 0``, then coatoms are facets, atoms vertices and vice versa.
    cdef ListOfFaces atoms, coatoms

    cdef inline CombinatorialFace next_face(self)
    cdef inline int next_dimension(self) except -1
    cdef size_t length_atom_repr(self) except -1
    cdef size_t set_coatom_repr(self) except -1
    cdef size_t set_atom_repr(self) except -1
