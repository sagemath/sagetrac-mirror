# distutils: sources = sage/geometry/polyhedron/combinatorial_polyhedron/bitsets.cpp sage/geometry/polyhedron/combinatorial_polyhedron/face.cpp sage/geometry/polyhedron/combinatorial_polyhedron/face_list.cpp
# distutils: depends = sage/geometry/polyhedron/combinatorial_polyhedron/bitsets.h sage/geometry/polyhedron/combinatorial_polyhedron/face.h sage/geometry/polyhedron/combinatorial_polyhedron/face_list.h
# distutils: include_dirs = sage/geometry/polyhedron/combinatorial_polyhedron
# distutils: language = c++
# distutils: extra_compile_args = -std=c++11
# distutils: libraries = gmp

cimport cython
from sage.ext.memory_allocator cimport MemoryAllocator

cdef extern from "face.h":
    cdef struct face_struct:
        pass

cdef extern from "face_list.h":
    cdef struct face_list_struct:
        face_struct* faces
        size_t n_faces
        size_t n_atoms
        int polyhedron_is_simple

cdef void allocate_one_face(face_struct& face, size_t n_faces, size_t n_atoms, MemoryAllocator mem)

@cython.final
cdef class ListOfFaces:
    cdef MemoryAllocator _mem

    # ``data`` points to the raw data.
    cdef face_list_struct data

    cdef inline size_t n_faces(self):
        return self.data.n_faces
    cdef inline size_t n_atoms(self):
        return self.data.n_atoms

    cpdef int compute_dimension(self) except -2

    cpdef ListOfFaces pyramid(self)
