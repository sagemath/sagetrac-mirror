from sage.plot.plot3d.base cimport PrimitiveObject

from transform cimport point_c, face_c, color_c

cdef class IndexFaceSet(PrimitiveObject):
    cdef bint enclosed
    cdef bint global_texture
    cdef Py_ssize_t vcount, fcount, icount
    cdef realloc(self, vcount, fcount, icount)
    # array of {x,y,z}
    cdef point_c* vs
    # array of array of fcount indices into points, each ending with -1
    cdef int* face_indices
    # pointers into face_indices marking the begining of each face
    cdef face_c* _faces

cdef class Vertex(object):
    cdef readonly double x, y, z

cdef class VertexSequence(object):
    cdef IndexFaceSet set

cdef class Face(object):
    cdef IndexFaceSet set
    cdef face_c *face

cdef class FaceSequence(object):
    cdef IndexFaceSet set

cdef class EdgeIter(object):
    cdef Py_ssize_t i, j
    cdef object seen
    cdef IndexFaceSet set
