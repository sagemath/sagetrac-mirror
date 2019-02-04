# The actual algorithm is implemented in the C++ file partitions_c.cc
# which requires the GMP, MPFR and NTL libraries.

#distutils: language = c++
#distutils: extra_compile_args = -march=native
r"""
Several algorithms working implicitly with the hasse_diagram
(of a polytope), including calculating the f_vector, dimension,
flags, level-sets and even the face lattice.

This is a wrapper for the functions in hasse_diagram.cc.

This computes implicitely a finite atomic and coatomic lattices,
where every interval of length two has at exactly 4 elements.
(Knwon as the diamond property).
In particular this module calculates quickly the f_vector of polytopes.
The input must be a tuple of coatoms given each by a tuple of atoms.
The atoms must be labeled 0,...,n.


AUTHOR:

- Jonathan Kliem (2019-01)


"""


#*****************************************************************************
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import
from sage.rings.integer import Integer
from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from sage.combinat.posets.lattices import FiniteLatticePoset
from sage.geometry.polyhedron.base import is_Polyhedron
from sage.geometry.lattice_polytope import is_LatticePolytope

from cpython cimport array
import array
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from sage.structure.sage_object cimport SageObject
from cysignals.memory cimport sig_malloc, sig_free
from cysignals.signals cimport sig_check, sig_on, sig_off, sig_on_no_except


cdef void * aligned_malloc(size_t size, size_t align):
    #taken from https://github.com/xjw/cpp/blob/master/cpp/memory_alignment.cpp
    #they are a workaround in case that C11 is not available

    # alignment could not be less than 0
    if (size<0):
        return NULL
    # allocate necessary memory for
    # alignment +
    # area to store the address of memory returned by malloc
    cdef void *p = sig_malloc(size + align-1 + sizeof(void *))
    if (p == NULL):
        raise MemoryError()

    # address of the aligned memory according to the align parameter
    cdef void *ptr = <void *> ((<unsigned long>p + sizeof(void *) + align-1) & ~(align-1))

    # store th address of mallc() above at the beginning of our total memory area
    (<void **>ptr -1)[0] = p

    # return the address of aligned memory
    return ptr

cdef void aligned_free(void *p):
    #taken from https://github.com/xjw/cpp/blob/master/cpp/memory_alignment.cpp

    # Get address of the memory from start of total memory area
    sig_free (( <void **>p - 1)[0] )

cdef class list_of_pointers:
    cdef void ** _data
    cdef void ** _memory
    # just to be on the safe side, we keep a pointer to the memory here
    # if the user messes with `data`, this should not hurt
    cdef size_t _length

    def __init__(self, size_t length):
        self._length = length
        self._memory = <void **> sig_malloc(length * sizeof(void *))
        self._data = self._memory

    def __dealloc__(self):
        sig_free(self._memory)

    cdef void ** data(self):
        return self._data

    cdef size_t length(self):
        return self._length

cdef class list_of_list_of_pointers:
    cdef void *** _data
    cdef void *** _memory
    cdef void *** _memory1
    # just to be on the safe side, we keep a pointer to the memory here
    # if the user messes with `data`, this should not hurt
    cdef size_t _length1
    cdef size_t _length2

    def __init__(self, size_t length1, size_t length2):
        cdef size_t i
        self._length1 = length1
        self._length2 = length2
        self._memory = <void ***> sig_malloc(length1 * sizeof(void **))
        self._memory1 = <void ***> sig_malloc(length1 * sizeof(void **))
        for i in range(length1):
            self._memory1[i] = <void **> sig_malloc(length2 * sizeof(void *))
            self._memory[i] = self._memory1[i]
        self._data = self._memory

    def __dealloc(self):
        cdef size_t i
        for i in range(self._length1):
            sig_free(self._memory1[i])
        sig_free(self._memory)
        sig_free(self._memory1)

    cdef void *** data(self):
        return self._data

    cdef size_t length1(self):
        return self._length1

    cdef size_t length2(self):
        return self._length2


cdef class list_of_faces:
    cdef void ** _data
    cdef void ** _memory
    cdef void ** _memory1
    cdef size_t _chunksize
    cdef size_t _nr_faces
    cdef size_t _length_of_face
    cdef size_t _max_nr_faces
    # `_data` is not supposed to change, if memory has been allocated
    # to be on the safe side, we will store the allocated memory in
    # in `_memory`

    def __init__ (self, size_t nr_faces, size_t length, size_t chunksize):
        cdef size_t i
        self._nr_faces = nr_faces
        self._length_of_face = ((length - 1)/chunksize + 1)
        # `_length_of_face` is the length in terms chunks of size `chunksize`
        self._memory = <void **> sig_malloc(nr_faces * sizeof(void *))
        self._data = self._memory
        self._chunksize = chunksize
        self._memory1 = <void **> PyMem_Malloc(nr_faces * sizeof(void *))
        for i in range(nr_faces):
            self._memory1[i] = aligned_malloc(
                self._length_of_face*chunksize/8, chunksize/8)
            self._memory[i] = self._memory1[i]

    def __dealloc__(self):
        cdef size_t i
        for i in range(self._nr_faces):
            aligned_free(self._memory1[i])
        sig_free(self._memory1)
        sig_free(self._memory)

    cdef make_copy(self):
        return list_of_faces(
            self._nr_faces, self.length(),
            self._chunksize)

    cdef void ** data(self):
        return self._data

    cdef size_t chunksize(self):
        return self._chunksize

    cdef size_t nr_faces(self):
        return self._nr_faces

    cdef size_t length(self):
        return self._length_of_face*self._chunksize

    cdef size_t length_of_face(self):
        return self._length_of_face

    cdef size_t max_nr_faces(self):
        return self._max_nr_faces

cdef class list_of_list_of_faces:
    cdef void *** _data
    cdef void *** _memory
    cdef tuple _lists
    cdef size_t _nr_lists

    def __init__(self, size_t nr_lists, size_t nr_faces, size_t length, size_t chunksize):
        cdef size_t i
        cdef list_of_faces saver
        self._nr_lists = nr_lists
        self._memory = <void ***> sig_malloc(nr_lists * sizeof(void **))
        self._data = self._memory
        self._lists = tuple(list_of_faces(nr_faces, length, chunksize) for i in range(nr_lists))
        for i in range(nr_lists):
            saver = self._lists[i]
            self._memory[i] = saver.data()

    def __dealloc__(self):
        sig_free(self._memory)

    cdef void *** data(self):
        return self._data

    cdef size_t nr_lists(self):
        return self._nr_lists

    cdef list_of_faces get_list(self, size_t index):
        return self._lists[index]


    """
    def face_iter_method(self, dimension):
        facets = self
        nr_facets = facets.nr_faces
        face_length = facets.length_of_face
        newfaces = tuple(facets.make_copy() for _ in range(dimension - 1))
        newfaces2 = tuple(facets.make_copy(allocation=False) for _ in range(dimension - 1))
        cdef list_of_faces forbidden = facets.make_copy()
        cdef list_of_faces saver
        cdef void *** newfacespointer = <void ***> PyMem_Malloc((dimension - 1) * sizeof(void **))
        cdef void *** newfaces2pointer = <void ***> PyMem_Malloc((dimension - 1) * sizeof(void **))
        cdef FaceIterator_ptr it
        cdef unsigned int face_dimension
        for i in range(dimension -1):
            saver = newfaces[i]
            newfacespointer[i] = saver.data
            saver = newfaces2[i]
            newfaces2pointer[i] = saver.data
        saver = facets
        yield 2
        it = FaceIterator_init(saver.data, newfacespointer, newfaces2pointer, forbidden.data, dimension, nr_facets, face_length)
        yield 2
        face_dimension = FaceIterator_call(it)
        while (face_dimension < dimension):
            yield face_dimension
            face_dimension = FaceIterator_call(it)
        FaceIterator_del(it)
        PyMem_Free(newfacespointer)
        PyMem_Free(newfaces2pointer)
    """


cdef int calculate_dimension(list_of_faces faces):
    cdef size_t nr_faces
    cdef int dim
    nr_faces = faces.nr_faces()
    if nr_faces == 0:
            raise TypeError('Wrong usage of `calculate_dimension`, at least one face needed.')

    return calculate_dimension_loop(faces.data(), nr_faces,
                                    faces.length_of_face(), faces.chunksize())


cdef int calculate_dimension_loop(void ** facesdata, size_t nr_faces, size_t face_length, size_t chunksize):
    cdef size_t bitcount, new_nr_faces
    cdef list_of_faces newfaces
    cdef void ** newfacesdata
    cdef list_of_pointers newfaces2
    cdef void ** newfaces2data
    cdef int dim

    if nr_faces == 0:
            raise TypeError('Wrong usage of `calculate_dimension_loop`, at least one face needed.')

    if nr_faces == 1:
        # we expect the face to be a vertex
        # possibly it contains more than one vertex/rays/lines
        # the dimension of the Polyhedron with this face as facet is
        # `bitcount + 1`
        bitcount = CountFaceBits(facesdata[0], face_length)
        return int(bitcount)

    newfaces = list_of_faces(nr_faces, face_length*chunksize, chunksize)
    newfaces2 = list_of_pointers(nr_faces)
    newfacesdata = newfaces.data()
    newfaces2data = newfaces2.data()
    new_nr_faces = get_next_level(facesdata, nr_faces, newfacesdata,
                         newfaces2data, newfacesdata, 0, face_length)
    return calculate_dimension_loop(newfaces2data, new_nr_faces,
                                    face_length, chunksize) + 1

cdef class FaceIterator:
    cdef void * face
    cdef int current_dimension, dimension, record_dimension, lowest_dimension
    cdef list_of_list_of_pointers newfaces2_mem
    cdef list_of_list_of_faces newfaces_mem
    cdef list_of_pointers forbidden_mem
    cdef void *** newfaces2
    cdef void *** newfaces
    cdef void ** forbidden
    cdef size_t * nr_faces
    cdef size_t * nr_forbidden
    cdef size_t yet_to_yield
    cdef int * first_time
    cdef size_t length_of_face

    def __init__(self, list_of_faces facets, int dimension, int nr_lines):
        if dimension <= 0:
            raise TypeError('FaceIterator expects positive dimensions')
        if facets.nr_faces() < 0:
            raise TypeError('FaceIterator expects non-empty `list_of_faces`')
        self.dimension = dimension
        self.current_dimension = dimension - 1
        self.length_of_face = facets.length_of_face()
        self.nr_faces = <size_t *> sig_malloc(dimension * sizeof(size_t))
        self.nr_faces[dimension - 1] = facets.nr_faces()
        self.nr_forbidden = <size_t *> sig_malloc(dimension * sizeof(size_t))
        self.nr_forbidden[dimension -1] = 0
        self.newfaces2_mem = list_of_list_of_pointers(<size_t> dimension,
                                                      facets.nr_faces())
        self.newfaces2 = self.newfaces2_mem.data()
        self.newfaces2[dimension - 1] = facets.data()
        self.newfaces_mem = \
            list_of_list_of_faces(dimension -1, facets.nr_faces(),
                                  facets.length(), facets.chunksize())
        self.newfaces = self.newfaces_mem.data()
        self.forbidden_mem = list_of_pointers(facets.nr_faces())
        self.forbidden = self.forbidden_mem.data()
        self.yet_to_yield = facets.nr_faces()
        self.record_dimension = -2
        self.lowest_dimension = nr_lines
        self.first_time = <int *> sig_malloc(dimension * sizeof(int))
        self.first_time[dimension - 1] = 1


    def __dealloc__(self):
        sig_free(self.nr_faces)
        sig_free(self.nr_forbidden)
        sig_free(self.first_time)

    cdef void set_record_dimension(self, int dim):
        self.record_dimension = dim
        self.lowest_dimension = max(0,dim)

    cdef int get_f_vector(self, size_t * f_vector):
        cdef dim = self.dimension
        f_vector[0] = 1
        f_vector[dim + 1] = 0 # the full-dimensional face will be added
        try:
            while (self.current_dimension < dim):
                while (self.next_face_loop() == 0) and (self.current_dimension < dim):
                    sig_check()
                f_vector[self.current_dimension + 1] += 1
        except:
            return 0
        return 1


    cdef inline int next_face(self):
        # this calls next_face_loop until it sets a new face
        # or until its consumed and `_current_dimension` is `_dimension`
        # **** Messing with the face_iterator *****
        # suppose face_iterator returns `face` and you do not want
        # to visit and farther faces of `face` you can do the following:
        # forbidden[face_iterator_nr_forbidden] = face;
        # face_iterator_nr_forbidden++;
        # This prevents any faces of `face` of appearing in the face iterator
        self.face = NULL
        cdef int dim = self.dimension
        while (not self.next_face_loop()) and (self.current_dimension < dim):
            pass
        return self.current_dimension

    cdef inline int next_face_loop(self):
        # sets `self._face` to the next face
        # might not do so, if so and
        # `self._current_dimension == self.dimension`
        # then there are no more faces
        cdef size_t nr_faces
        cdef size_t nr_forbidden
        cdef void **faces
        cdef size_t i, newfacescounter
        if (self.current_dimension == self.dimension):
            # the function is not supposed to be called in this case
            # just to prevent it from crashing
            return 0

        nr_faces = self.nr_faces[self.current_dimension]
        nr_forbidden = self.nr_forbidden[self.current_dimension]
        faces = self.newfaces2[self.current_dimension]
        if (self.record_dimension > -2) and \
                (self.record_dimension != self.current_dimension):
            # if we are not in dimension `record_dimension`,
            # then we should not yield any faces
            # (in case `face_iterator_dimension == -2` we yield all faces)
            self.yet_to_yield = 0

        if self.yet_to_yield:
            # return the next face
            self.yet_to_yield -= 1
            self.face = faces[self.yet_to_yield]
            return 1

        if self.current_dimension <= self.record_dimension:
            # if we do not want to yield lower dimensional faces,
            # than we should go up one dimension again to look for more faces
            # (act as if we had visited all faces in lower dimensions already)
            self.current_dimension += 1
            return 0

        if self.current_dimension == self.lowest_dimension:
            # we will not yield the empty face
            # we will not yield below what is wanted
            self.current_dimension += 1
            return 0

        if nr_faces <= 1:
            # there will be no more faces from intersections
            self.current_dimension += 1
            return 0

        i = nr_faces - 1;
        self.nr_faces[self.current_dimension] -= 1
        if not self.first_time[self.current_dimension]:
            # if there exists faces[i+1], we have visited all its faces already
            # hence we should not visit any of them again
            self.forbidden[nr_forbidden] = faces[i+1]
            self.nr_forbidden[self.current_dimension] += 1
            nr_forbidden = self.nr_forbidden[self.current_dimension]

        else:
            # we will visit all the faces of faces[nr_faces]
            # once we have done so, we want to add this face to forbidden
            self.first_time[self.current_dimension] = 0

        # get the facets contained in faces[i] but not in any of the forbidden
        newfacescounter = \
            get_next_level(faces, i+1, self.newfaces[self.current_dimension-1],
                           self.newfaces2[self.current_dimension-1],
                           self.forbidden, nr_forbidden, self.length_of_face)

        if newfacescounter:
            # if there are new faces contained in faces[i+1],
            # we will set up the variables to correctly visit them on the next
            # call of `next_face_loop`
            self.current_dimension -= 1
            self.first_time[self.current_dimension] = 1
            self.nr_faces[self.current_dimension] = newfacescounter
            self.nr_forbidden[self.current_dimension] = nr_forbidden
            self.yet_to_yield = newfacescounter
            return 0

        #else:
            # if there are no faces in lower dimension,
            # then there is no need to add the face to forbidden
            # this might become important when calculating simpliness
            # and simpliciality, where we will mess with the iterator
            # and add some faces to forbidden in order to not consider subfaces
            # self.first_time[self.current_dimension] = 1

        return 0

    cdef size_t nr_vertices(self):
        if self.face:
            return CountFaceBits(self.face, self.length_of_face)

        # the face was not initialized properly
        raise LookupError('The `FaceIterator` does not point to a face.')

    cdef tuple which_ridge(self):
        r"""
        Assumes that self.face is a ridge. Returns the tuple of the facets.
        """
        cdef size_t nr_facets = self.nr_faces[self.dimension - 1]
        cdef void ** facets = self.newfaces2[self.dimension - 1]
        cdef size_t length_of_face = self.length_of_face
        cdef size_t one = nr_facets
        cdef size_t two = nr_facets
        cdef size_t counter = 0
        cdef size_t counter2 = 0
        while counter < nr_facets:
            if is_subset(self.face, facets[counter], length_of_face):
                one = counter
                counter2 = counter + 1
                counter = nr_facets
            counter += 1
        while counter2 < nr_facets:
            if is_subset(self.face, facets[counter2], length_of_face):
                two = counter2
                counter2 = nr_facets
            counter2 += 1
        return (one, two)

    cdef list ridges(self):
        self.set_record_dimension(self.dimension - 2)
        cdef list ridges
        cdef size_t counter = 0
        cdef size_t length = 1024
        ridges = [0] * 1024
        while self.next_face() != self.dimension:
            if counter + 1 > length:
                ridges += [0] * length
                length *= 2
            ridges[counter] = self.which_ridge()
            counter += 1
        return ridges[:counter]






cdef extern from "hasse_diagram.cc":
    cdef const unsigned int chunksize;
    # A pointer to the underlying C++ class of CombinatorialPolyhedron.
    ctypedef void * CombinatorialPolyhedron_ptr

#    r"""
#    CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(
#       unsigned int ** facets_pointer, unsigned int nr_facets,
#       unsigned int *len_facets, unsigned int nr_vertices,
#       int is_unbounded)

#    Initalizes the C++ class of CombinatorialPolyhedron
#    and returns a pointer to it.

#    INPUT:
#     - ``facets_pointer`` to a list of facets a lists of vertices/rays/lines.
#       The vertices/rays/lines need to be labeled 0,...,``nr_vertices - 1``.
#     - ``nr_facets`` the number of facets/the length of ``facets_pointer.
#     - ``len_facets`` the length of each face, i.e.
#       ``facetspointer[i]`` should be of length ``len_facets[i]``.
#     - ``is_unbounded`` needs to be 0, if the Polyhedron is unbounded,
#       otherwise 1 + ``nr_lines``.
#       ``is_unbounded = 0`` and ``is_unbounded = 1`` give the same result
#       for bounded Polyhedra, but ``is_unbounded = 0`` might be faster.

#    WARNING::

#        - ``nr_facets`` needs to be at least two.
#    """
    cdef CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(
        unsigned int ** facets_pointer, unsigned int nr_facets,
        unsigned int * len_facets, unsigned int nr_vertices, int is_unbounded)

#    r"""
#    CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(
#       unsigned int ** incidence_matrix, unsigned int nr_facets,
#       unsigned int nr_vertices, int is_unbounded)

#    Initalizes the C++ class of CombinatorialPolyhedron and returns a
#    pointer to it.

#    INPUT:
#     - An ``incidence_matrix`` given by a list of lenght ``nr_facets`` of
#       a list of ``nr_vertices`` of 0s and 1s.
#     - ``nr_facets`` the number of facets.
#     - ``nr_vertices`` the number of vertices/rays/lines.
#     - ``is_unbounded`` needs to be 0, if the Polyhedron is unbounded,
#       otherwise 1 + ``nr_lines``.
#       ``is_unbounded = 0`` and ``is_unbounded = 1`` give the same result
#       for bounded Polyhedra, but ``is_unbounded = 0`` might be faster.

#    WARNING::

#        - ``nr_facets`` needs to be at least two.
#    """
    cdef CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(
        unsigned int ** incidence_matrix, unsigned int nr_facets,
        unsigned int nr_vertices, int is_unbounded)

#    r"""
#    delete_CombinatorialPolyhedron(CombinatorialPolyhedron_ptr)

#    Needs to be called to deallocate the memory used by the
#    CombinatorialPolyhedron object in C++
#    """
    cdef void delete_CombinatorialPolyhedron(CombinatorialPolyhedron_ptr)

#    r"""
#    unsigned int dimension(CombinatorialPolyhedron_ptr C)

#    Returns the dimension of the Polyhedron.

#    This is implicitly calculated and stored for most other functions.
#    """
    cdef unsigned int dimension(CombinatorialPolyhedron_ptr C)

#    r"""
#    void f_vector(CombinatorialPolyhedron_ptr C, unsigned long *vector)

#    Stores the f-vector of the Polyhedron in ``vector``. Vector needs to
#    be of length ``dimension + 2``
#    """
    cdef void f_vector(CombinatorialPolyhedron_ptr C, unsigned long * vector)

#    r"""
#    unsigned int ** edges(CombinatorialPolyhedron_ptr C)

#    Calculates the edges of the Polyhedron and returns a pointer to them.

#    The result ``edges`` will be a pointer of pointers to unsigned integers.
#    Each integer represents a vertex/ray/line in the order of the input at
#    initialization.

#    After getting the pointer, you can get the nr of edges by the
#    corresponding entry in the f-vector.

#    ``edges[0]`` stores ``maxnumberedges`` edges.
#    (see ``get_maxnumberedges()``)
#    The first edge will have vertices ``edges[0][0]`` and ``edges[0][1]``,
#    the second edges ``edges[0][2]`` and ``edges[0][3]`` and so on.
#    The ``maxnumberedges + 1``st edge will have vertices
#    ``edges[1][0]`` and ``edges[1][1]``.

#    WARNING:

#        - ``edges`` will have at most length ``maxnumberedges``, so there will
#          be at most ``maxnumberedges^2`` edges stored.

#        - ``edges[1]`` is not initialized if the nr of edges is at most
#        ``maxnumberedges``.
#        - ``edges[i]`` is not initialized if the nr of edges is
#          too small.

#    NOTE:
#        Might implicitly calculate the f-vector, but NOT vice versa.

#    """
    cdef unsigned int ** edges(CombinatorialPolyhedron_ptr C)

#    r"""
#    unsigned int ** ridges(CombinatorialPolyhedron_ptr C)

#    See ** edges(CombinatorialPolyhedron_ptr C).
#    """
    cdef unsigned int ** ridges(CombinatorialPolyhedron_ptr C)

#    r"""
#    unsigned long ** incidences(CombinatorialPolyhedron_ptr C,
#                                int dimension_one, int dimension_two,
#                                unsigned long * nr_incidences,
#                                unsigned int * twisted)

#    Calculates the incidences of ``dimension_one``-faces with
#    ``dimension-two``-faces.

#    Return a pointer to ``incidences`` to the incedences similar to ``edges``.

#    INPUT:

#        - pointer to the C++-class of CombinatorialPolyhedron
#        - ``dimension_one`` and ``dimension_two`` the dimensions
#          of the faces of which the incidences will be given
#        - ``nr_incidences`` an ``unsigned long[1]`` in which the nr_of
#          incidences will be stored
#        - ``twisted`` an ``unsigned int[1]``, which will indicate if
#          the result is twisted/flipped

#    OUTPUT:
#        - ``incidences`` of type
#          ``unsigned long[maxnumberincidences][maxnumberincidences]``.
#          ``incidences[0]`` stores ``maxnumberincidences`` incidenes.
#          (see ``get_maxnumberincidences()``)
#          incidences[1] stores another ``maxnumberincidences``
#          of incidences etc.
#        - ``twisted[0]`` indicates the output corresponding to an input
#          where ``dimension_one`` and ``dimension_two`` are interchanged.
#        - ``nr_incidences[0]`` the number of incidences

#    The first incidences will have faces ``incidences[0][0]``
#    and ``incidences[0][1]``,
#    where the first is of dimension ``dimension_one`` and the second is
#    of length ``dimension_two``.
#    The faces are given as their index in the list of all faces as in
#    ``get_faces``.

#    WARNING:

#        - ``incidences`` will have at most length ``maxnumberincidences``,
#          so there will be at most ``maxnumberincidences^2`` edges stored.
#        - ``incidences[1]`` is not initialized if the nr of incidences
#          at most ``maxnumberincidences``.
#        - ``incidences[i]`` is not initialized if the nr of incidences is
#          to small.

#    NOTE:
#        Implicetely calculates a list of all faces of ``dimension_one``
#        and ``dimension_two``. AFTER calling ``incidences`` the order in
#        ``get_faces`` is the same as in ``face_iterator``.

#    """
    cdef unsigned long ** incidences(CombinatorialPolyhedron_ptr C,
                                     int dimension_one, int dimension_two,
                                     unsigned long * nr_incidences,
                                     unsigned int * twisted)

#    """
#    ``record_all_faces`` generates a list of all faces in C++ as to speed
#    up different calculations later.
#    """
    cdef void record_all_faces(CombinatorialPolyhedron_ptr C)

#    """
#    get_faces(CombinatorialPolyhedron_ptr C, int dimension,
#              unsigned int facet_repr, unsigned int **faces_to_return,
#              unsigned int *length_of_faces)

#    Gets a list of all faces.

#    INPUT:

#        - ``dimension`` of the faces to return.
#          Should be at least -1 and at most ``dimension.
#        - ``facet_repr`` -- 0 for the V-represenation of the faces.
#          1 for the H-representation.
#        - ``faces_to_return`` -- an array of type
#          ``unsigned int[nr_faces][nr_vertices]`` if ``facet_repr == 0``
#          and
#          ``unsigned int[nr_faces][nr_facets]`` if ``facet_repr == 1`]
#          NOTE:
#            - get the correct nr of faces by the entry in the f-vector
#            - nr_vertices is the nr of vertices/rays/lines, i.e. the
#              length of the V-representation given on initialization
#            - If you really know what you are doing, you can initialize
#              ``faces_to_return[i]`` just long enough to make the longest
#              face fit.
#        - ``length_of_faces`` -- an array of typ ``unsigned int[nr_faces]

#    OUTPUT:

#        - faces_to_return will contain the faces as list of
#          vertices, if ``facet_repr == 0``, or
#          facets, if ``facet_repr == 1``.
#        - ``faces_to_return[i]`` will contain a face of length
#          ``length_of_faces[i]``.
#    """
    void get_faces(CombinatorialPolyhedron_ptr C, int dimension,
                   unsigned int facet_repr, unsigned int ** faces_to_return,
                   unsigned int * length_of_faces)

#    """
#    face_iterator_init(CombinatorialPolyhedron_ptr C, int dimension,
#                       unsigned int vertex_repr, unsigned int facet_repr)

#    Will initalize a face_iterator.

#    INPUT:

#        - ``dimension`` of the faces to iterator over, if
#          ``dimension == -2`` then the iterator will use al faces
#        - ``vertex_repr`` needs to be set to 1 if vertex_repr is wanted,
#          0 otherwise
#        - ``facet_repr`` needs to beset to 1 if facet_repr is wanted,
#          both facet_repr and vertex_repr are possible

#    WARNING:

#        - There can be at most one working iterator around. The iterator
#          should be used up before calling different functions, especially,
#            - `faces`
#            - `vertices`
#            - `edges` (a second call is fine)
#            - `ridges` (a second call is fine)
#            - `f_vector` (a second call is fine)
#            - `record_all_faces`
#            - `incidences`
#            - `face_lattice`
#            - `flag`
#            - `k-simplicial` #not implemented yet
#            - `k-simple`    # not implemented yet
#    """
    cdef void face_iterator_init(CombinatorialPolyhedron_ptr C,
                                 int dimension, unsigned int vertex_repr,
                                 unsigned int facet_repr)

#    """
#    unsigned int face_iterator(
#       CombinatorialPolyhedron_ptr C, unsigned int *Vface_to_return,
#       unsigned int *Vlength, unsigned int *Hface_to_return,
#       unsigned int *Hlength)

#    WARNING:

#        Do not call the face_iterator unless it has been initalized prior
#        by ``face_iterator_init``.

#        There can only be one face_iterator around.

#    INPUT:

#        - Pointer to CombinatorialPolyhedron in C.
#        - ``Vface_to_return`` of type ``unsigned int[nr_vertices]``.
#        - ``Vlength`` of type ``unsigned int[1]``
#        - ``Hface_to_return`` of type ``unsigned int[nr_facets]``.
#        - ``Hlength`` of type ``unsigned int[1]``.

#    OUTPUT:

#        - returns 0 if there was no next face, 1 otherwise
#        - if ``vertex_repr == 1`` (initalization)
#          ``Vface_to_return`` will contain ``Vlength[0]`` integers
#          corresponding to the vertices of the face
#        - if ``facet_repr == 1`` (initalization)
#          ``Hface_to_return`` will contain ``Hlength[0]`` integers
#          corresponding to the vertices of the face
#    """
    cdef unsigned int face_iterator(
        CombinatorialPolyhedron_ptr C, unsigned int * Vface_to_return,
        unsigned int * Vlength, unsigned int * Hface_to_return,
        unsigned int * Hlength)

#    """
#    get_flag(CombinatorialPolyhedron_ptr C, unsigned int *flagarray,
#             unsigned int length)

#    Returns the number of flags of type ``flagarray``.
#    ``Flagarray`` needs to be a sorted array of length ``length``,
#    ``Flagarray contain integers from 0 up to ``dimension - 1``, where
#    ``dimensions`` is the dimension of the CombinatorialPolyhedron.
#    ``dimension = dimension(CombinatorialPolyhedron_ptr C)``
#    """
    cdef unsigned long get_flag(CombinatorialPolyhedron_ptr C,
                                unsigned int * flagarray, unsigned int length)


#    """
#    ``get_maxnumberedges`` is a helper function for
#    ``edges(CombinatorialPolyhedron_ptr C)`` and
#    ``ridges(CombinatorialPolyhedron_ptr C)``

#    Edges and ridges are stored in an array of form
#    ``edges = unsigned int[maxnumberedges][maxnumberedges*2]``.
#    In order get the edges from the pointer, one needs this number.
#    ``edges`` does not store more than ``maxnumberedges^2`` edges.
#    """
    cdef unsigned long get_maxnumberedges()

#    """
#    ``get_maxnumberincidences`` is a helper function for
#    ``incidences(CombinatorialPolyhedron_ptr C, int dimension_one,
#                 int dimension_two, unsigned long * nr_incidences,
#                 unsigned int * twisted)``

#    Incidences are stored in an array of form
#    ``incidences = \
#       unsigned long[maxnumberincidences][maxnumberincidences*2]``.
#    In order get the incidences from the pointer, one needs this number.
#    ``incidences`` does not store more than ``maxnumberincidences^2``
#    incidences
#    """
    cdef unsigned long get_maxnumberincidences()


cdef extern from "helper.cc":
    cdef inline void intersection(void *A1, void *B1, void *C1, size_t face_length)

    cdef inline int is_subset(void *A1, void *B1, size_t face_length)

    cdef size_t get_next_level(void **faces, size_t lenfaces, void **nextfaces, \
                               void **nextfaces2, void **forbidden, \
                               size_t nr_forbidden, size_t face_length)
    # intersects the first `lenfaces - 1` faces of `faces` with'faces[lenfaces-1]`
    # stores the faces in `newfaces`
    # determines which ones are exactly of one dimension less
    # by considering containment
    # newfaces2 will point at those of exactly one dimension less
    # which are not contained in any of the faces in 'forbidden'
    # returns the number of those faces

    cdef size_t CountFaceBits(void* A1, size_t face_length)

    cdef void get_facets_from_incidence_matrix(unsigned int **incidence_matrix, void **facets, \
                                      size_t nr_vertices, size_t nr_facets)

    cdef void get_vertices_from_incidence_matrix(unsigned int **incidence_matrix, void **vertices, \
                                        size_t nr_vertices, size_t nr_facets)

    cdef void get_vertices_bitrep_from_facets_pointer( \
        unsigned int ** facets_input, unsigned int *len_facets, \
        void ** vertices_output, size_t nr_vertices, size_t nr_facets)

    cdef void get_facets_bitrep_from_facets_pointer( \
        unsigned int ** facets_input, unsigned int *len_facets, \
        void ** facets_output, size_t nr_vertices, size_t nr_facets)


cdef class CombinatorialPolyhedron(SageObject):
    r"""
    The class of the Combinatorial Type of a Polyehdron, a Polytope.

    INPUT:

    - ``data`` -- a ``Polyhedron``, i.e. an instance of
      :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`.

    or

    - ``data`` -- a ``LatticePolytope``, i.e. an instance of
      :class:`~sage.geometry.lattice_polytope.LatticePolytopeClass`.

    or

    - ``data`` -- an ``incidence_matrix`` as in
      :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`
      of :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`.

      * ``vertices`` -- a list of ``[vertices, rays, lines]``, if
        the rows in the incidence_matrix should correspond to names.

      * ``facets`` -- a list of facets, if
        the columns in the incidence_matrix should correspond to names.

      * ``nr_lines`` -- for bounded Polyhedra, this should be the
      default: ``None``. For unbounded Polyhedra, this needs to be set
      to the correct nr of lines, i.e. the the maximum nr of lines with
      linearly independent directions in the Polyehdron.

    or

    - ``data`` -- a list of facets,
      each facet given as a list of ``[vertices, rays, lines]``.
      If the Polyhedron is unbounded, then rays and lines are required.
      If the Polyehdron contains no lines the rays can be thought of as
      the vertices of the facets deleted from a bounded Polyhedron. See
      :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`
      on how to use rays and lines.

      * ``facets`` -- a list of names of the facets, if
        the facets given should correspond to names.

      * ``unbounded`` -- for bounded Polyhedra, this should be the
      default ``None``. For unbounded Polyhedra, this needs to be set
      to the correct nr of lines, i.e. the the maximum nr of lines with
      linearly independent directions in the Polyehdron.


    EXAMPLES:

    Input is Polyhedron::

        sage: P = polytopes.cube()
        sage: CombinatorialPolyhedron(P)
        Combinatorial Type of a Polyhedron of dimension 3 with 8 vertices

    Input is a LatticePolytope::

        sage: points = [(1,0,0), (0,1,0), (0,0,1),
        ....: (-1,0,0), (0,-1,0), (0,0,-1)]
        sage: L = LatticePolytope(points)
        sage: CombinatorialPolyhedron(L)
        Combinatorial Type of a Polyhedron of dimension 3 with 6 vertices

    Input is an incidence matrix::

        sage: data = Polyhedron(rays=[[0,1]]).incidence_matrix()
        sage: CombinatorialPolyhedron(data, nr_lines=0)
        Combinatorial Type of a half-space of dimension 1
        sage: C = CombinatorialPolyhedron(data, vertices=['myvertex'],
        ....: facets=['myfacet'], nr_lines=0)
        sage: C.Vrepresentation()
        ('myvertex',)
        sage: C.Hrepresentation()
        ('myfacet',)

    You can also give the facets explicitely::

        sage: CombinatorialPolyhedron(((1,2,3),(1,2,4),(1,3,4),(2,3,4)))
        Combinatorial Type of a Polyhedron of dimension 3 with 4 vertices
        sage: facetnames = ['facet0', 'facet1', 'facet2', 'myfacet3']
        sage: facetinc = ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        sage: C = CombinatorialPolyhedron(facetinc, facets=facetnames)
        sage: C.Vrepresentation()
        (1, 2, 3, 4)
        sage: C.Hrepresentation()
        ('facet0', 'facet1', 'facet2', 'myfacet3')

    Specifying the nr of lines is important::

        sage: P = Polyhedron(ieqs=[[1,-1,0],[1,1,0]])
        sage: C = CombinatorialPolyhedron(P) #this works fine
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices
        sage: data = P.incidence_matrix()
        sage: C = CombinatorialPolyhedron(data) #wrong usage!
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 1 vertices
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C = CombinatorialPolyhedron(data, nr_lines=1) #correct
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices
        sage: C.f_vector()
        (1, 0, 2, 1)



    """
    cdef CombinatorialPolyhedron_ptr _C
    cdef tuple _V
    cdef tuple _H
    cdef tuple _equalities
    cdef dict _Hinv
    cdef dict _Vinv
    cdef int is_trivial  # in some instances the polyhedron might not
    # have facets or otherwise produce errors in the C function
    cdef int _dimension  # manually set, if `is_trivial`
    cdef unsigned int _length_Hrep
    cdef unsigned int _length_Vrep
    cdef int _unbounded  # set to 0, if the Polyhedron is bounded,
    # otherwise it is 1 + nr_lines
    cdef int _nr_lines
    cdef list_of_faces bitrep_facets
    cdef list_of_faces bitrep_vertices
    cdef int _polar
    cdef unsigned int chunksize
    cdef size_t * _f_vector

    def __init__(self, data, vertices=None, facets=None, nr_lines=None):
        r"""
        Initializes the combinatorial polyhedron.

        See :class:`CombinatorialPolyhedron`.

        TESTS::

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],
            ....: [0,2,3],[1,2,3]])    # indirect doctests
        """
        self.chunksize = chunksize
        cdef unsigned int ** incidence_matrix
        cdef unsigned int ** facets_pointer
        cdef unsigned int * len_facets
        self._dimension = -2
        self._f_vector = NULL
        self._polar = 0
        self._nr_lines = 0
        if nr_lines is None:
            self._unbounded = 0
        else:
            self._unbounded = 1 + int(nr_lines)
            self._nr_lines = int(nr_lines)
        self.is_trivial = 0
        self._equalities = ()
        if is_Polyhedron(data):
            if data.is_empty():
                self.is_trivial = 1
                self._dimension = -1
                return
            vertices = data.Vrepresentation()
            facets = tuple(inequality for inequality in data.Hrepresentation())
            # in this case the Polyhedron does not have facets:
            if (len(vertices) == data.n_lines() + 1) and (data.n_lines > 0):
                self.is_trivial = 1
                self._dimension = data.n_lines()
                self._V = tuple(vertices)
                return
            self._unbounded
            if not data.is_compact():
                self._unbounded = 1 + data.n_lines()
                self._nr_lines = int(data.n_lines())
            data = data.incidence_matrix()
            self._length_Hrep = data.ncols()
            self._length_Vrep = data.nrows()

        if is_LatticePolytope(data):
            if data.npoints() == 0:
                self.is_trivial = 1
                self._dimension = -1
                return
            if data.npoints() == 1:
                self.is_trivial = 1
                self._dimension = 0
                self._V = data.vertices()
                return
            self._unbounded = 0
            vertices = data.vertices()
            self._length_Vrep = len(vertices)
            facets = data.facets()
            self._length_Hrep = len(facets)
            data = tuple(tuple(vert for vert in facet.vertices())
                         for facet in facets)

        if vertices:
            self._V = tuple(vertices)
            self._Vinv = {v: i for i,v in enumerate(self._V)}
        else:
            self._V = None
            self._Vinv = None

        if facets:
            facets = tuple(facets)
            test = [0] * len(facets)
            for i in range(len(facets)):
                test[i] = 1
                if hasattr(facets[i], "is_inequality"):
                    if not facets[i].is_inequality():
                        test[i] = 0
            self._H = \
                tuple(facets[i] for i in range(len(facets)) if test[i])
            # only keeping those that are actual inequalities
            self._equalities = \
                tuple(facets[i] for i in range(len(facets)) if not test[i])
            # the inequalities are saved here
            self._Hinv = {v: i for i,v in enumerate(self._H)}
        else:
            self._H = None
            self._Hinv = None

        if hasattr(data, "incidence_matrix"):
            # TODO: Better check for incidence_matrix
            data = data.incidence_matrix()

        if hasattr(data, "nrows"):  # TODO: Better check for matrix
            self._length_Hrep = data.ncols()
            self._length_Vrep = data.nrows()
            rg = range(data.nrows())
            tup = tuple(tuple(data[i,j] for i in rg)
                        for j in range(data.ncols())
                        if not all(data[i,j] for i in rg))
            # transpose and get rid of equalities (which all vertices satisfie)

            if len(tup) == 0:  # the case of the empty Polyhedron
                self.is_trivial = 1
                self._dimension = -1 + data.nrows()
                # the elements in Vrep are assumed to be one vertex
                # and otherwise lines
                return

            if len(tup) == 1:  # the case of a half space
                self.is_trivial = 2
                self._dimension = -1 + data.nrows()
                # the elements in Vrep are assumed to be one vertex
                # and otherwise lines
                return

            incidence_matrix = \
                <unsigned int**> PyMem_Malloc(len(tup)*
                                              sizeof(unsigned int *))
            for i in range(len(tup)):
                incidence_matrix[i] = \
                    <unsigned int*> PyMem_Malloc(self._length_Vrep *
                                                 sizeof(unsigned int))
                for j in range(self._length_Vrep):
                    incidence_matrix[i][j] = tup[i][j]
            self._C = init_CombinatorialPolyhedron(
                incidence_matrix, len(tup), self._length_Vrep, self._unbounded)

            nr_vertices = self._length_Vrep
            nr_facets = len(tup)
            if self._unbounded or nr_facets <= nr_vertices:
                self._polar = 0
                # initializing the facets as BitVectors
                self.bitrep_facets = \
                    list_of_faces(nr_facets, nr_vertices,
                                  self.chunksize)
                get_facets_from_incidence_matrix(
                    incidence_matrix, self.bitrep_facets.data(),
                    nr_vertices, nr_facets)
                # initializing the vertices as BitVectors
                self.bitrep_vertices = \
                    list_of_faces(nr_vertices, nr_facets,
                                  self.chunksize)
                get_vertices_from_incidence_matrix(
                    incidence_matrix, self.bitrep_vertices.data(),
                    nr_vertices, nr_facets)
            else:
                self._polar = 1
                # initializing the vertices as BitVectors
                self.bitrep_vertices = \
                    list_of_faces(nr_facets, nr_vertices,
                                  self.chunksize)
                get_facets_from_incidence_matrix(
                    incidence_matrix, self.bitrep_vertices.data(),
                    nr_vertices, nr_facets)
                # initializing the facets as BitVectors
                self.bitrep_facets = \
                    list_of_faces(nr_vertices, nr_facets,
                                  self.chunksize)
                get_vertices_from_incidence_matrix(
                    incidence_matrix, self.bitrep_facets.data(),
                    nr_vertices, nr_facets)

            # cleanup
            for i in range(len(tup)):
                PyMem_Free(incidence_matrix[i])
            PyMem_Free(incidence_matrix)

        elif isinstance(data, Integer):  # intput for a trivial Polyhedron
            if data < -1:
                TypeError("A polyhedron must have dimension at least -1")
            self.is_trivial = 1
            self._dimension = data

        else:  # assumes the facets are given as a list of vertices/rays/lines

            if len(data) == 0:
                    self.is_trivial = 1
                    self._dimension = -1
                    return

            if len(data) == 1:
                    self.is_trivial = 2
                    self._dimension = len(data[0]) - 1
                    if self._dimension <= 0:
                        self.is_trivial = 1
                        # we are treating a polyhedron equal to its affine hull
                    return

            if self._V is None:
                vertices = sorted(set.union(*map(set, data)))
                nr_vertices = len(vertices)
                if vertices != range(len(vertices)):
                    self._V = tuple(vertices)
                    self._Vinv = {v: i for i,v in enumerate(self._V)}
            else:
                nr_vertices = len(self._V)
            self._length_Vrep = nr_vertices

            if self._V is not None:
                def f(v): return self._Vinv[v]
            else:
                def f(v): return int(v)

            facets = tuple(tuple(f(i) for i in j) for j in data)
            self._length_Hrep = len(facets)
            facets_pointer = \
                <unsigned int**> PyMem_Malloc(len(facets) * sizeof(unsigned int *))
            len_facets = \
                <unsigned int*> PyMem_Malloc(len(facets) *
                                             sizeof(unsigned int))
            for i in range(len(facets)):
                len_facets[i] = len(facets[i])
                facets_pointer[i] = \
                    <unsigned int*> PyMem_Malloc(len_facets[i] *
                                                 sizeof(unsigned int))
                for j in range(len_facets[i]):
                    facets_pointer[i][j] = facets[i][j]
            self._C = init_CombinatorialPolyhedron(facets_pointer, len(facets),
                                                   len_facets, nr_vertices,
                                                   self._unbounded)

            if self._unbounded or len(facets) <= nr_vertices:
                self._polar = 0
                # initializing the facets as BitVectors
                self.bitrep_facets = \
                    list_of_faces(len(facets), nr_vertices,
                                  self.chunksize)
                get_facets_bitrep_from_facets_pointer(
                    facets_pointer, len_facets, self.bitrep_facets.data(),
                    nr_vertices, len(facets))
                # initializing the vertices as BitVectors
                self.bitrep_vertices = \
                    list_of_faces(nr_vertices, len(facets),
                                  self.chunksize)
                get_vertices_bitrep_from_facets_pointer(
                    facets_pointer, len_facets, self.bitrep_vertices.data(),
                    nr_vertices, len(facets))
            else:
                self._polar = 1
                # initializing the vertices as BitVectors
                self.bitrep_vertices = \
                    list_of_faces(len(facets), nr_vertices,
                                  self.chunksize)
                get_facets_bitrep_from_facets_pointer(
                    facets_pointer, len_facets, self.bitrep_vertices.data(),
                    nr_vertices, len(facets))
                # initializing the facets as BitVectors
                self.bitrep_facets = \
                    list_of_faces(nr_vertices, len(facets),
                                  self.chunksize)
                get_vertices_bitrep_from_facets_pointer(
                    facets_pointer, len_facets, self.bitrep_facets.data(),
                    nr_vertices, len(facets))

            # cleanup
            for i in range(len(facets)):
                PyMem_Free(facets_pointer[i])
            PyMem_Free(facets_pointer)
            PyMem_Free(len_facets)

    def __dealloc__(self):
        r"""
        This function deallocates all the memory used by underlying C++-class
        """
        if self.is_trivial > 0:
            return
        delete_CombinatorialPolyhedron(self._C)
        if self._f_vector:
            sig_free(self._f_vector)

    def _repr_(self):
        r"""
        Returns a description of the Combinatorial Polyhedron.

        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a Polyhedron of dimension 3 with 4 vertices'

            sage: P = Polyhedron(vertices=[])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of the empty Polyhedron'

            sage: P = Polyhedron(vertices=[[0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of the Polyhedron with one vertex'

            sage: P = Polyhedron(lines=[[0,0,1],[0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a trivial Polyhedron of dimension 2'

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[-1,0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a half-space of dimension 2'
        """
        if self.is_trivial == 1:
            if self._dimension == 0:
                return "Combinatorial Type of the Polyhedron with one vertex"

            if self._dimension > 0:
                return "Combinatorial Type of a trivial Polyhedron of \
                    dimension %s" % self._dimension

            return "Combinatorial Type of the empty Polyhedron"

        if self.is_trivial == 2:
            return "Combinatorial Type of a half-space of dimension %s"\
                % self._dimension

        return "Combinatorial Type of a Polyhedron of dimension %s with %s vertices" \
                % (self.dimension(), len(self.vertices()))

    def __reduce__(self):
        r"""
        Override __reduce__ to correctly pickle/unpickle.

        TESTS::

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: tuple(i for i in C.face_iter(facet_repr=True)) \
            ....: == tuple(i for i in C1.face_iter(facet_repr=True))
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: tuple(i for i in C.face_iter(facet_repr=True)) \
            ....: == tuple(i for i in C1.face_iter(facet_repr=True))
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0],
            ....:                      [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: C.f_vector() == C1.f_vector()
            True
        """
        unbounded = None
        if self._unbounded:
            unbounded = Integer(unbounded) - 1

        if self.is_trivial == 1:
            return (CombinatorialPolyhedron, (Integer(self._dimension),
                    self.Vrepresentation(), self.Hrepresentation(), unbounded))

        if self.is_trivial == 2:
            pickletuple = tuple(0 for _ in range(self._dimension + 1))
            return (CombinatorialPolyhedron, ((pickletuple,),
                    self.Vrepresentation(), self.Hrepresentation(), unbounded))
            # if `data` is a `tuple` containining exactly one tuple
            # of some length, than __init__ will figure it is a halfspace

        return (CombinatorialPolyhedron, (self.facets(),
                self.Vrepresentation(), self.Hrepresentation(), unbounded))

    def Vrepresentation(self):
        r"""
        Return a list of names of ``[vertices, rays, lines]``.

        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0,0], [0,1,0], \
            ....:                      [0,0,1],[0,0,-1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.Vrepresentation()
            (A line in the direction (0, 0, 1),
             A ray in the direction (1, 0, 0),
             A vertex at (0, 0, 0),
             A ray in the direction (0, 1, 0))

        """
        if self._V is not None:
            return self._V
        else:
            return tuple(Integer(i) for i in range(self._length_Vrep))

    def Hrepresentation(self):
        r"""
        Returns a list of names of facets and possibly some equalities.

        EXAMPLES::

            sage: P = polytopes.permutahedron(3)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.Hrepresentation()
            (An equation (1, 1, 1) x - 6 == 0,
             An inequality (0, -1, -1) x + 5 >= 0,
             An inequality (0, 0, -1) x + 3 >= 0,
             An inequality (0, -1, 0) x + 3 >= 0,
             An inequality (0, 1, 0) x - 1 >= 0,
             An inequality (0, 1, 1) x - 3 >= 0,
             An inequality (0, 0, 1) x - 1 >= 0)

        """
        if self._H is not None:
            return self._equalities + self._H
        else:
            return tuple(Integer(i) for i in range(self._length_Hrep))

    def vertices(self, names=True):
        r"""
        Returns the elements in the ``Vrepresentation`` that are vertices.

        In the case of an unbounded Polyhedron, there might be lines and
        rays in the Vrepresentation.

        If ``names`` is set to ``False``, then the vertices are given by
        their indices in the Vrepresentation.

        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[0,0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.vertices()
            (A vertex at (0, 0, 0),)
            sage: C.Vrepresentation()
            (A vertex at (0, 0, 0),
             A ray in the direction (0, 0, 1),
             A ray in the direction (0, 1, 0),
             A ray in the direction (1, 0, 0))
            sage: P = polytopes.cross_polytope(3)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.vertices()
            (A vertex at (-1, 0, 0),
             A vertex at (0, -1, 0),
             A vertex at (0, 0, -1),
             A vertex at (0, 0, 1),
             A vertex at (0, 1, 0),
             A vertex at (1, 0, 0))
            sage: C.vertices(names=False)
            (0, 1, 2, 3, 4, 5)

            sage: points = [(1,0,0), (0,1,0), (0,0,1),
            ....:           (-1,0,0), (0,-1,0), (0,0,-1)]
            sage: L = LatticePolytope(points)
            sage: C = CombinatorialPolyhedron(L)
            sage: C.vertices()
            (M(1, 0, 0), M(0, 1, 0), M(0, 0, 1), M(-1, 0, 0), M(0, -1, 0), M(0, 0, -1))
            sage: C.vertices(names=False)
            (0, 1, 2, 3, 4, 5)

        """

        return tuple(i[0] for i in self.faces(0, names=names))

    def facets(self, names=True):
        r"""
        Returns the facets as lists of ``[vertices, rays, lines]``.

        If ``names`` is ``False``, then the vertices in the facets
        are given by their indices in the Vrepresentation.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.facets()
            ((A vertex at (-1, -1, 1),
              A vertex at (-1, 1, 1),
              A vertex at (1, -1, 1),
              A vertex at (1, 1, 1)),
             (A vertex at (-1, 1, -1),
              A vertex at (-1, 1, 1),
              A vertex at (1, 1, -1),
              A vertex at (1, 1, 1)),
             (A vertex at (1, -1, -1),
              A vertex at (1, -1, 1),
              A vertex at (1, 1, -1),
              A vertex at (1, 1, 1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, -1, 1),
              A vertex at (-1, 1, -1),
              A vertex at (-1, 1, 1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, 1, -1),
              A vertex at (1, -1, -1),
              A vertex at (1, 1, -1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, -1, 1),
              A vertex at (1, -1, -1),
              A vertex at (1, -1, 1)))
            sage: C.facets(names=False)
            ((1, 3, 5, 7),
             (2, 3, 6, 7),
             (4, 5, 6, 7),
             (0, 1, 2, 3),
             (0, 2, 4, 6),
             (0, 1, 4, 5))


        """
        return self.faces(self.dimension()-1, names=names)

    def edges(self, names=True):
        r"""
        Returns the edges of the CombinatorialPolyhedron,
        i.e. the rank 1 faces, which contain 2 vertices.

        If ``names`` is set to ``False``, then the vertices in the edges
        are given by their indices in the Vrepresentation.

        If you want to compute all faces of dimension 1,
        use :meth:`CombinatorialPolyhedron.faces` instead.

        .. NOTE::

            To compute edges and f_vector, first compute edges.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(3,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ((A vertex at (0, 0, 0), A vertex at (4, 16, 64)),
             (A vertex at (1, 1, 1), A vertex at (4, 16, 64)),
             (A vertex at (2, 4, 8), A vertex at (4, 16, 64)),
             (A vertex at (3, 9, 27), A vertex at (4, 16, 64)),
             (A vertex at (0, 0, 0), A vertex at (3, 9, 27)),
             (A vertex at (1, 1, 1), A vertex at (3, 9, 27)),
             (A vertex at (0, 0, 0), A vertex at (2, 4, 8)),
             (A vertex at (1, 1, 1), A vertex at (2, 4, 8)),
             (A vertex at (0, 0, 0), A vertex at (1, 1, 1)))
            sage: C.edges(names=False)
            ((0, 4), (1, 4), (2, 4), (3, 4), (0, 3), (1, 3), (0, 2), (1, 2), (0, 1))

            sage: P = Polyhedron(rays=[[-1,0],[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ()
            sage: C.faces(1)
            ((A line in the direction (1, 0), A vertex at (0, 0)),)

            sage: P = Polyhedron(vertices=[[0,0],[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ((A vertex at (0, 0), A vertex at (1, 0)),)

        """
        if self.is_trivial > 0:
            return ()

        cdef unsigned int ** edgepointer = edges(self._C)
        cdef unsigned long maxnumberedges = get_maxnumberedges()
        # the edges are being saved in a list basically
        # with the first entry the first vertex of the first edges,
        # the second entry the second vertex of that edge

        # possibly there are many edges,
        # hence they are stored in an array of arrays,
        # with each array containing maxnumberedges of edges

        if self.dimension() <= 0:
            return ()

        cdef unsigned long nr_edges = self.f_vector()[2]
        if self.f_vector()[1] < 2:
            # in case of a cone or a Polyhedron there are no edges
            return ()

        if self.dimension() == 1:
            # the C-module assumes at least dimension 2 to give edges
            return self.faces(1, names=names)

        if nr_edges > maxnumberedges*maxnumberedges:
            raise ValueError("Cannot calculate %s edges" % nr_edges)

        if self._V is not None and names is True:
            def f(i): return self._V[i]
        else:
            def f(i): return Integer(i)

        def vertex_one(i):
            return f(edgepointer[i / maxnumberedges][2*(i % maxnumberedges)])

        def vertex_two(i):
            return f(edgepointer[i / maxnumberedges][2*(i % maxnumberedges)+1])
        return tuple((vertex_one(i), vertex_two(i)) for i in range(nr_edges))

        # TODO: In the unbounded case, one should check if the vertices
        # of edges are vertices of the Polyhedron
        #
        # P = Polyhedron(vertices=[[0,0],[1,0]], rays=[[0,1]])

    def edge_graph(self, names=True):
        r"""
        Returns the edge graph.

        If ``names`` is set to ``False``, the vertices will carry names
        according to the indexing of the Vrepresentation.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(3,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edge_graph()
            Graph on 5 vertices
            sage: G = C.edge_graph()
            sage: G.degree()
            [3, 4, 4, 3, 4]

        """

        return Graph(self.edges(names=names), format="list_of_edges")

    def dimension(self):
        r"""
        Returns the dimension of the ``CombinatorialPolyehdron``.

        EXAMPLES::

            sage: C = CombinatorialPolyhedron([(1,2,3), (1,2,4),
            ....:                              (1,3,4), (2,3,4)])
            sage: C.dimension()
            3

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
            sage: CombinatorialPolyhedron(P).dimension()
            3


        """
        if self._dimension == -2:
            self._dimension = calculate_dimension(self.bitrep_facets)
        return Integer(self._dimension)

    def ridges(self, add_equalities=False, names=True):
        r"""
        Returns the ridges.

        The ridges of the CombinatorialPolyhedron are the faces
        contained in exactly two facets.

        If you want to compute all faces of codimension 1,
        use :meth:`CombinatorialPolyhedron.faces` instead.

        The ridges will be given by the facets, they are contained in.

        - If ``add_equalities`` is ``True``, then equalities the entire
          Polyhedron satisfies, are added.

        - If ``names`` is ``False``, then the facets in the ridges are
          given by their indices in the Hrepresentation.

        .. NOTE::

            To compute ridges and f_vector, compute ridges first.

        EXAMPLES::

            sage: P = polytopes.permutahedron(2)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.ridges()
            ((An inequality (0, -1) x + 2 >= 0, An inequality (0, 1) x - 1 >= 0),)
            sage: C.ridges(add_equalities=True)
            (((An equation (1, 1) x - 3 == 0, An inequality (0, -1) x + 2 >= 0),
              (An equation (1, 1) x - 3 == 0, An inequality (0, 1) x - 1 >= 0)),)

            sage: P = polytopes.cyclic_polytope(4,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.ridges()
            ((An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-12, 19, -8, 1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (24, -26, 9, -1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (-12, 19, -8, 1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (-12, 19, -8, 1) x + 0 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (-12, 19, -8, 1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (8, -14, 7, -1) x + 0 >= 0))
            sage: C.ridges(names=False)
            ((0, 4),
             (1, 4),
             (2, 4),
             (3, 4),
             (0, 3),
             (1, 3),
             (2, 3),
             (0, 2),
             (1, 2),
             (0, 1))

            sage: P = Polyhedron(rays=[[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C
            Combinatorial Type of a half-space of dimension 1
            sage: C.ridges()
            ()
            sage: C.faces(0, facet_repr=True)
            ((An equation (0, 1) x + 0 == 0, An inequality (1, 0) x + 0 >= 0),)


        """

        if self.is_trivial > 1:
            return ()
        if self.f_vector()[-2] <= 1:
            # in some instances there are no ridges as intersections of facets,
            # i.e. if we treat a ray, or a halfspace
            return ()

        # The new approach is going to look somewhat like this.
        """
        cdef FaceIterator face_iter
        cdef list_of_faces facets = self.bitrep_facets
        cdef int dim = self.dimension()
        cdef list ridges
        face_iter = FaceIterator(facets, dim, self._nr_lines)
        ridges = face_iter.ridges()
        if self._H is not None and names is True:
            def f(i): return self._H[i]
        else:
            def f(i): return Integer(i)
        if add_equalities:
            return tuple(self._equalities + (f(i[0]),f(i[1])) for i in ridges)
        else:
            return tuple((f(i[0]),f(i[1])) for i in ridges)
        """

        # old version
        cdef unsigned int ** ridgepointer = ridges(self._C)
        cdef unsigned long maxnumberedges = get_maxnumberedges()
        # the ridges are being saved in a list basically,
        # with the first entry the first facet of the first ridge,
        # the second entry the second facet of that ridges
        # as there might be many ridges they are saved in an array of arrays,
        # with each array containing maxnumberedges of ridges
        if self.dimension() <= 0:
            return ()
        cdef unsigned long nr_ridges = self.f_vector()[-3]
        if nr_ridges > maxnumberedges*maxnumberedges:
            raise ValueError("Cannot calculate %s ridges" % nr_ridges)

        if self._H is not None and names is True:
            def f(i): return self._H[i]
        else:
            def f(i): return Integer(i)

        def facet_one(i):
            return f(ridgepointer[i/maxnumberedges][2*(i % maxnumberedges)])

        def facet_two(i):
            return f(ridgepointer[i/maxnumberedges][2*(i % maxnumberedges)+1])

        if add_equalities:
            return tuple(
                ((self._equalities + (facet_one(i),)),
                 (self._equalities + (facet_two(i),))) for i in range(nr_ridges))
        else:
            return tuple((facet_one(i), facet_two(i)) for i in range(nr_ridges))

    def ridge_graph(self, names=True):
        r"""
        Returns the ridge graph.

        The ridge graph of the CombinatorialPolyhedron consists of
        ridges as edges and facets as vertices.

        If ``names`` is ``False``, the ``vertices`` of the graph  will
        be the incidences of the facets in the Hrepresentation.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(4,6)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.ridge_graph()
            Graph on 9 vertices

        """
        return Graph(self.ridges(names=names), format="list_of_edges")

    def f_vector(self):
        r"""
        Calculates the ``f_vector`` of the CombinatorialPolyhedron.

        The ``f_vector`` contains the number of faces of dimension ``k``
        for each ``k`` in ``range(-1, self.dimension() + 1)``.

        .. NOTE::

            If you also want to compute edges and/or ridges, do so
            first.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.f_vector()
            (1, 120, 240, 150, 30, 1)

            sage: P = polytopes.cyclic_polytope(6,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.f_vector()
            (1, 10, 45, 120, 185, 150, 50, 1)

        TESTS::

            sage: from itertools import combinations
            sage: N = combinations(range(20),19)
            sage: C = CombinatorialPolyhedron(list(N))
            sage: try:
            ....:     alarm(0.0001)
            ....:     C.f_vector()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: C.f_vector()
            (1,
             20,
             190,
             1140,
             4845,
             15504,
             38760,
             77520,
             125970,
             167960,
             184756,
             167960,
             125970,
             77520,
             38760,
             15504,
             4845,
             1140,
             190,
             20,
             1)


        ALGORITHM:

        The number of facets is assumed to be at least two here.
        The algorithm to visit all proper faces exactly once is roughly
        equivalent to:

        facets = [set(facet) for facet in self.facets()]
        ComputeNextStep(facets, [])

        #this algorithm assumes at each step to receive all facets of
        #some face except those contained in a face of forbidden
        def ComputeNextStep(faces, forbidden):

            for face in faces:
                pass #do something here with that face

            while len(faces) > 1:
                one_face = faces.pop()
                newfaces = [one_face.intersection(face) for face in faces]
                #newfaces contains all intersection

                newfaces2 = []
                for face1 in newfaces:
                    # face1 is a facet of one_face iff
                    # it is not contained in another facet
                    if all(not face1 < face2 for face2 in newfaces):
                        newfaces2.append(face1)
                #newfaces2 contains all facets of one_face not contained
                #in any one of forbidden and maybe some that are
                #contained in one of forbidden

                newfaces3 = []
                for face1 in newfaces2:
                    if all(not face1 < face2 for face2 in forbidden):
                        newfaces3.append(face1)
                #newfaces3 contains exactly all facets of one_face but
                #those contained in one face of forbidden

                #visit all faces in one_face that are not contained in
                #one of forbidden
                ComputeNextStep(newfaces3, forbidden)

                #we have visited all faces in one_face, so we should not
                #visit one ever again

                forbidden.append(one_face)

            return


        """
        cdef int dim = self.dimension()
        if self._f_vector is NULL:
            self._calculate_f_vector()
        if self._f_vector is NULL:
            raise KeyboardInterrupt('Interrupt on user intput')
        if self._polar:
            f = tuple(Integer(self._f_vector[dim+1-i]) for i in range(dim + 2))
        else:
            f = tuple(Integer(self._f_vector[i]) for i in range(dim + 2))
        return f

    cdef void _calculate_f_vector(self):
        cdef FaceIterator face_iter
        cdef list_of_faces facets
        cdef int dim = self.dimension()
        cdef int sucess = 0
        self._f_vector = <size_t *> sig_malloc((dim + 2)*sizeof(size_t))
        for i in range(dim + 2):
            self._f_vector[i] = 0
        self._f_vector[0] = 1
        self._f_vector[dim + 1] = 1

        if self.is_trivial == 1:
            return
        if self.is_trivial == 2:
            self._f_vector[dim] = 1
            return

        facets = self.bitrep_facets
        face_iter = FaceIterator(facets, dim, self._nr_lines)
        #if not sig_on_no_except():
        if not face_iter.get_f_vector(self._f_vector):
            sig_free(self._f_vector)
            self._f_vector = NULL

    def _record_all_faces(self):
        r"""
        Records all faces of the Polyhedron. For quicker acces later.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C._record_all_faces()
        """
        if self.is_trivial >= 1:
            return
        record_all_faces(self._C)

    def faces(self, dimension, facet_repr=False, names=True):
        r"""
        Gets all k-faces for specified dimenion k.

        By default faces are given as tuple of vertices.

        If ``facet_repr`` is set to ``True``, then vertices are given in
        as tuple of facets.

        If ``names`` is set to ``False``, then the vertices are given by
        their indices in the Vrepresentation.


        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.faces(2)
            ((A vertex at (0, 0, 0, 1),
              A vertex at (0, 0, 1, 0),
              A vertex at (0, 1, 0, 0)),
             (A vertex at (0, 0, 0, 1),
              A vertex at (0, 0, 1, 0),
              A vertex at (1, 0, 0, 0)),
             (A vertex at (0, 0, 0, 1),
              A vertex at (0, 1, 0, 0),
              A vertex at (1, 0, 0, 0)),
             (A vertex at (0, 0, 1, 0),
              A vertex at (0, 1, 0, 0),
              A vertex at (1, 0, 0, 0)))
            sage: C.faces(2, facet_repr=True)
            ((An equation (1, 1, 1, 1) x - 1 == 0,
              An inequality (0, -1, -1, -1) x + 1 >= 0),
             (An equation (1, 1, 1, 1) x - 1 == 0, An inequality (0, 1, 0, 0) x + 0 >= 0),
             (An equation (1, 1, 1, 1) x - 1 == 0, An inequality (0, 0, 1, 0) x + 0 >= 0),
             (An equation (1, 1, 1, 1) x - 1 == 0, An inequality (0, 0, 0, 1) x + 0 >= 0))

            sage: P = Polyhedron(rays=[[1,0],[0,1]], vertices=[[1,0],[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.faces(1)
            ((A ray in the direction (0, 1), A vertex at (0, 1)),
             (A ray in the direction (1, 0), A vertex at (1, 0)),
             (A vertex at (0, 1), A vertex at (1, 0)))
            sage: C.faces(1, names=False)
            ((0, 1), (2, 3), (1, 3))
            sage: C.faces(1, facet_repr=True)
            ((An inequality (1, 0) x + 0 >= 0,),
             (An inequality (0, 1) x + 0 >= 0,),
             (An inequality (1, 1) x - 1 >= 0,))

        ALGORITHM:

        See :meth:`f_vector` for a description
        on how all faces are visited.
        """
        if self.is_trivial == 1:
            # the case of a Polyhedron equal to its affine hull
            if dimension == -1:
                return ((),)
            if dimension == self._dimension:
                if facet_repr is True:
                    return ((),)
                if self._V is not None and names is True:
                    return (tuple(self._V),)
                else:
                    if dimension == 0:
                        return (Integer(0),)
                    return ((),)
            return ()

        if self.is_trivial == 2:  # the case of a half-space
            if dimension == -1:
                if facet_repr is True:
                    if self._H is not None and names is True:
                        return (self._equalities + tuple(self._H),)
                    else:
                        return ((Integer(0),),)
                return ((),)
            if dimension == self._dimension - 1:
                if facet_repr is True:
                    if self._H is not None and names is True:
                        return (self._equalities + tuple(self._H),)
                    else:
                        return ((Integer(0),),)
                return ((),)
            if dimension == self._dimension:
                return ((),)
            return ()

        dim = self.dimension()
        if dimension not in range(-1, dim + 1):
            return ()

        # creating an array of arrays for the C function to store the faces in
        cdef unsigned long number_of_faces = self.f_vector()[dimension + 1]
        cdef unsigned int length_of_face
        if facet_repr is False:
            size_of_face = self._length_Vrep
        else:
            size_of_face = self._length_Hrep
        cdef unsigned int ** faces_to_return = \
            <unsigned int**> PyMem_Malloc(number_of_faces *
                                          sizeof(unsigned int *))
        for i in range(number_of_faces):
            faces_to_return[i] = \
                <unsigned int*> PyMem_Malloc(size_of_face *
                                             sizeof(unsigned int))
        cdef unsigned int * length_of_faces = \
            <unsigned int*> PyMem_Malloc(number_of_faces
                                         * sizeof(unsigned int))

        # filling the array
        get_faces(self._C, dimension, facet_repr,
                  faces_to_return, length_of_faces)

        # translating the result to the desired representation
        addtuple = ()
        if (facet_repr):
            facet_repr = 1
            if self._H is not None and names is True:
                def f(i): return self._H[i]
            else:
                def f(i): return Integer(i)
            addtuple = self._equalities
        else:
            facet_repr = 0
            if self._V is not None and names is True:
                def f(i): return self._V[i]
            else:
                def f(i): return Integer(i)

        returntuple = \
            tuple(addtuple + tuple(f(faces_to_return[facecounter][j])
                  for j in range(length_of_faces[facecounter]))
                  for facecounter in range(number_of_faces))

        # cleaning up
        for i in range(number_of_faces):
            PyMem_Free(faces_to_return[i])
        PyMem_Free(faces_to_return)
        PyMem_Free(length_of_faces)

        return returntuple

    def face_iter(self, dimension=None, vertex_repr=True,
                  facet_repr=False, names=True):
        r"""
        Iterator over all faces of specified dimension.

        If ``dimension`` is not specified then iterate over all faces.

        If ``vertex_repr`` is ``True``, then give the faces as lists of
        elements in ``Vrepresentation``.

        If ``facet_repr`` is ``True``, then give the faces as lists of
        elements in ``Hrepresentation``.

        - Both ``vertex_repr`` and ``facet_repr`` can be set to ``True,
          this will give each face as tuple of the form
          (``vertex_repr``, ``facet_rerpr``).

        If ``names`` is ``False`` then vertices and facets are labeled
        by their indexes.

        .. WARNING::

            There can only be one face iterator around. The second
            one will invalid the first one. Even worse, if you then call
            the first one again, the output will be incorrect and the
            second one will not yield all faces.

            This shares resources with other methods. Do not call other
            methods while iteration, s.t.

            - :meth:`faces`
            - :meth:`vertices`
            - :meth:`edges` (a second call is fine)
            - :meth:`ridges` (a second call is fine)
            - :meth:`f_vector` (a second call is fine)
            - :meth:`_record_all_faces`
            - :meth:`incidences`
            - :meth:`face_lattice`
            - :meth:`flag`
            - :meth:`k-simplicial`
            - :meth:`k-simple`

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(dimension=2)
            sage: next(it)
            (A vertex at (4, 1, 5, 2, 3),
             A vertex at (4, 2, 5, 1, 3),
             A vertex at (5, 1, 4, 2, 3),
             A vertex at (5, 2, 4, 1, 3))
            sage: next(it)
            (A vertex at (4, 1, 5, 2, 3),
             A vertex at (4, 1, 5, 3, 2),
             A vertex at (5, 1, 4, 2, 3),
             A vertex at (5, 1, 4, 3, 2))
            sage: it = C.face_iter(dimension=2, names=False)
            sage: next(it)
            (76, 82, 100, 106)
            sage: next(it)
            (76, 77, 100, 101)
            sage: it = C.face_iter(dimension=2, facet_repr=True)
            sage: next(it)
            ((A vertex at (4, 1, 5, 2, 3),
              A vertex at (4, 2, 5, 1, 3),
              A vertex at (5, 1, 4, 2, 3),
              A vertex at (5, 2, 4, 1, 3)),
             (An equation (1, 1, 1, 1, 1) x - 15 == 0,
              An inequality (0, 1, 0, 1, 0) x - 3 >= 0,
              An inequality (0, 1, 0, 1, 1) x - 6 >= 0))
            sage: it = C.face_iter(dimension=2, vertex_repr=False,
            ....:                  facet_repr=True, names=False)
            sage: next(it)
            (28, 29)
            sage: next(it)
            (25, 29)


        TESTS::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: f = C.f_vector()
            sage: altf = tuple(len(tuple(C.face_iter(i)))
            ....:              for i in range(-1,C.dimension()+1))
            sage: altf == f
            True
            sage: allfaces = tuple(tuple(C.face_iter(i))
            ....:                  for i in range(-1,C.dimension()+1))
            sage: allfaces2 = tuple(tuple(C.face_iter(i))
            ....:                   for i in range(-1,C.dimension()+1))
            sage: all(sorted(sorted(j) for j in allfaces[i]) == \
            ....:   sorted(sorted(j) for j in allfaces2[i])
            ....:          for i in range(-1, C.dimension()+1))
            True

            sage: P = polytopes.cyclic_polytope(5,20)
            sage: C = CombinatorialPolyhedron(P)
            sage: f = C.f_vector()
            sage: altf = tuple(len(tuple(C.face_iter(i)))
            ....:              for i in range(-1,C.dimension()+1))
            sage: altf == f
            True
            sage: allfaces = tuple(tuple(C.face_iter(i))
            ....:                  for i in range(-1,C.dimension()+1))
            sage: C._record_all_faces()
            sage: allfaces2 = tuple(tuple(C.face_iter(i))
            ....:                   for i in range(-1,C.dimension()+1))
            sage: all(sorted(sorted(j) for j in allfaces[i]) == \
            ....: sorted(sorted(j) for j in allfaces2[i])
            ....:        for i in range(-1, C.dimension()+1))
            True
            sage: C = CombinatorialPolyhedron(P)
            sage: allfaces = tuple(tuple(C.face_iter(i, vertex_repr=False,
            ....:                                    facet_repr=True)
            ....:                  for i in range(-1, C.dimension()+1)))
            sage: C._record_all_faces()
            sage: allfaces2 = tuple(tuple(C.face_iter(i, vertex_repr=False,
            ....:                                     facet_repr=True)
            ....:                         for i in range(-1, C.dimension()+1)))
            sage: all(sorted(sorted(j) for j in allfaces[i]) == \
            ....: sorted(sorted(j) for j in allfaces2[i])
            ....:        for i in range(-1, C.dimension()+1))
            True



        ALGORITHM:

        See :meth:`f_vector` for a description
        on how all faces are visited.
        """

        if dimension is not None:
            dimension = int(dimension)
            dimensionrange = (dimension,)
        else:
            dimensionrange = range(-1, self.dimension()+1)
            dimension = -2

        if not facet_repr:
            vertex_repr = True

        if self.is_trivial > 0:  # taking care of the trivial polynomial
            for dim in dimensionrange:
                if vertex_repr and facet_repr:
                    vert = self.faces(dim, names=names)
                    fac = self.faces(dim, facet_repr=True, names=names)
                    for i in range(len(vert)):
                        yield (vert[i], fac[i])
                elif vertex_repr:
                    vert = self.faces(dim, names=names)
                    for i in range(len(vert)):
                        yield vert[i]
                else:
                    fac = self.faces(dim, facet_repr=True, names=names)
                    for i in range(len(fac)):
                        yield fac[i]
            return

        #return face_iter(self.bitrep_facets, self.dimension())

        if 0 == dimension:
            if vertex_repr and facet_repr:
                vert = self.faces(0, names=names)
                fac = self.faces(0, facet_repr=True, names=names)
                for i in range(len(vert)):
                    yield (vert[i], fac[i])
            elif vertex_repr:
                vert = self.faces(0, names=names)
                for i in range(len(vert)):
                    yield vert[i]
            else:
                fac = self.faces(0, facet_repr=True, names=names)
                for i in range(len(fac)):
                    yield fac[i]
            return

        if dimension == self.dimension() - 1:
            if vertex_repr and facet_repr:
                vert = self.faces(dimension, names=names)
                fac = self.faces(dimension, facet_repr=True, names=names)
                for i in range(len(vert)):
                    yield (vert[i], fac[i])
            elif vertex_repr:
                vert = self.faces(dimension, names=names)
                for i in range(len(vert)):
                    yield vert[i]
            else:
                fac = self.faces(dimension, facet_repr=True, names=names)
                for i in range(len(fac)):
                    yield fac[i]
            return

        if -1 in dimensionrange:
            if vertex_repr and facet_repr:
                vert = self.faces(-1, names=names)
                fac = self.faces(-1, facet_repr=True, names=names)
                for i in range(len(vert)):
                    yield (vert[i], fac[i])
            elif vertex_repr:
                vert = self.faces(-1, names=names)
                for i in range(len(vert)):
                    yield vert[i]
            else:
                fac = self.faces(-1, facet_repr=True, names=names)
                for i in range(len(fac)):
                    yield fac[i]
            if -1 == dimension:
                return

        dim = self.dimension()
        if dim in dimensionrange:
            if vertex_repr and facet_repr:
                vert = self.faces(dim, names=names)
                fac = self.faces(dim, facet_repr=True, names=names)
                for i in range(len(vert)):
                    yield (vert[i], fac[i])
            elif vertex_repr:
                vert = self.faces(dim, names=names)
                for i in range(len(vert)):
                    yield vert[i]
            else:
                fac = self.faces(dim, facet_repr=True, names=names)
                for i in range(len(fac)):
                    yield fac[i]
            if dim == dimension:
                return

        dim = self.dimension()
        if dimension not in range(-2, dim + 1):
            return ()

        # creating two arrays for the C function to store the faces in
        cdef unsigned int successful
        size_of_faceV = self._length_Vrep
        size_of_faceH = self._length_Hrep
        cdef unsigned int VSize[1]
        VSize[0] = 0
        cdef unsigned int HSize[1]
        HSize[0] = 0
        cdef unsigned int * Vface_to_return = \
            <unsigned int*> PyMem_Malloc(size_of_faceV * sizeof(unsigned int))
        cdef unsigned int * Hface_to_return = \
            <unsigned int*> PyMem_Malloc(size_of_faceH * sizeof(unsigned int))

        # translating the result to the desired representation
        addtuple = self._equalities
        if self._H is not None and names is True:
            def h(i): return self._H[i]
        else:
            def h(i): return Integer(i)

        if self._V is not None and names is True:
            def v(i): return self._V[i]
        else:
            def v(i): return Integer(i)

        addtuple = ()
        if names:
            addtuple = self._equalities

        # init
        face_iterator_init(self._C, dimension, int(vertex_repr),
                           int(facet_repr))

        # filling the array
        successful = face_iterator(self._C, Vface_to_return, VSize,
                                   Hface_to_return, HSize)
        while successful:
            if vertex_repr and facet_repr:
                yield (tuple(v(Vface_to_return[i]) for i in range(VSize[0])),
                       addtuple + tuple(h(Hface_to_return[i])
                       for i in range(HSize[0])))
            elif vertex_repr:
                yield tuple(v(Vface_to_return[i]) for i in range(VSize[0]))
            else:
                yield addtuple + tuple(h(Hface_to_return[i])
                                       for i in range(HSize[0]))
            successful = face_iterator(self._C, Vface_to_return, VSize,
                                       Hface_to_return, HSize)

        # cleaning up
        PyMem_Free(Vface_to_return)
        PyMem_Free(Hface_to_return)

    def incidences(self, dimension_one, dimension_two):
        r"""
        Returns all incidences between ``dimension_one``-faces and
        ``dimension_two``-faces.

        Incidences are given as tuple of integers, where the integer
        corresponds to the enumeration according to :meth:`faces`

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.incidences(0,1)
            ((0, 0),
             (4, 0),
             (0, 1),
             (1, 1),
             (4, 2),
             (5, 2),
             (1, 3),
             (5, 3),
             (0, 4),
             (2, 4),
             (4, 5),
             (6, 5),
             (2, 6),
             (6, 6),
             (2, 7),
             (3, 7),
             (1, 8),
             (3, 8),
             (6, 9),
             (7, 9),
             (5, 10),
             (7, 10),
             (3, 11),
             (7, 11))
            sage: C.incidences(1,2)
            ((3, 0),
             (8, 0),
             (10, 0),
             (11, 0),
             (6, 1),
             (7, 1),
             (9, 1),
             (11, 1),
             (2, 2),
             (5, 2),
             (9, 2),
             (10, 2),
             (1, 3),
             (4, 3),
             (7, 3),
             (8, 3),
             (0, 4),
             (4, 4),
             (5, 4),
             (6, 4),
             (0, 5),
             (1, 5),
             (2, 5),
             (3, 5))
            sage: C.incidences(0,2)
            ((1, 0),
             (3, 0),
             (5, 0),
             (7, 0),
             (2, 1),
             (3, 1),
             (6, 1),
             (7, 1),
             (4, 2),
             (5, 2),
             (6, 2),
             (7, 2),
             (0, 3),
             (1, 3),
             (2, 3),
             (3, 3),
             (0, 4),
             (2, 4),
             (4, 4),
             (6, 4),
             (0, 5),
             (1, 5),
             (4, 5),
             (5, 5))
        """
        if self.is_trivial == 1:
            if dimension_one in (-1, self._dimension) \
                    and dimension_two in (-1, self._dimension):
                return ((0, 0),)
            return ()
        if self.is_trivial == 2:
            if dimension_one in (-1, self._dimension, self._dimension - 1) \
                and dimension_two in (-1, self._dimension,
                                      self._dimension - 1):
                return ((0, 0),)
            return ()

        cdef unsigned long maxnumberincidences = get_maxnumberincidences()
        cdef unsigned long nr_incidences[1]
        nr_incidences[:] = [0]
        cdef unsigned int twisted[1]
        twisted[:] = [0]
        cdef unsigned long ** incidencepointer \
            = incidences(self._C, dimension_one, dimension_two,
                         nr_incidences, twisted)

        if nr_incidences[0] > maxnumberincidences*maxnumberincidences:
            raise ValueError("Cannot calculate %s incidences"
                             % nr_incidences[0])

        if twisted[0] == 0:
            def incidence_one(i):
                return Integer(incidencepointer[i / maxnumberincidences]
                                               [2*(i % maxnumberincidences)])

            def incidence_two(i):
                return Integer(incidencepointer[i / maxnumberincidences]
                                               [2*(i % maxnumberincidences)+1])
        else:
            def incidence_two(i):
                return Integer(incidencepointer[i / maxnumberincidences]
                                               [2*(i % maxnumberincidences)])

            def incidence_one(i):
                return Integer(incidencepointer[i / maxnumberincidences]
                                               [2*(i % maxnumberincidences)+1])
        return tuple((incidence_one(i), incidence_two(i))
                     for i in range(nr_incidences[0]))

    def new_face_iter(self):
        cdef FaceIterator face_iter
        cdef list_of_faces facets = self.bitrep_facets
        cdef int dim = self.dimension()
        cdef int next_dim
        face_iter = FaceIterator(facets, dim, self._nr_lines)
        next_dim = face_iter.next_face()
        while next_dim is not dim:
            yield next_dim
            next_dim = face_iter.next_face()

    def face_lattice(self, vertices=False, facets=False, names=False):
        r"""
        Generates the face-lattice.

        INPUT:

        - ``vertices`` -- if set to ``True`` the elements in the lattice
          will be named according to the vertices they contain.

        - ``facets`` -- if set to ``True`` the elements in the lattice
          will be named according to the facets they are contained in.

        - ``names`` -- if set to ``False``, facets and vertices will be
          named according to their index.

        Both vertices and facets can be set to ``True``.

        In the case of a trivial Polyhedron, which is equal to its own
        affine hull, ``facets`` will be set to ``False``, as the
        elements need distinct names.

        In the case of a half-space ``vertices`` and ``facets`` will be
        set to ``False``.

        OUTPUT:

        - a :class:'~sage.combinat.posets.lattices.FiniteLatticePoset'

        .. NOTE::

            As :class:'~sage.combinat.posets.lattices.FiniteLatticePoset'
            is awfully slow with elements having meaningful labels,
            the default of this function is to not do so.


        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0],[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.incidences(0,1)
            ((0, 0), (0, 1))
            sage: C.face_lattice()
            Finite lattice containing 5 elements
            sage: C.face_lattice(vertices=True).atoms()
            [(A vertex at (0, 0),)]

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: P1 = Polyhedron(rays=[[1,0], [-1,0]])
            sage: C1 = CombinatorialPolyhedron(P1)
            sage: C.face_lattice().is_isomorphic(C1.face_lattice())
            True

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice()
            Finite lattice containing 542 elements
        """

        # we must ignore part of the input to ensure an injective relabeling
        if self.is_trivial == 1:
            facets = False
        if self.is_trivial == 2:
            vertices = False
            facets = False

        f_vector = self.f_vector()
        self._record_all_faces()
        dimension = self.dimension()
        dic = {}
        range_f_vector = [k for k in range(-1, dimension + 1)
                          if f_vector[k+1] > 0]
        range_f_vector1 = range_f_vector[1:-1]
        for k in range_f_vector:
            faces = (self.faces(k), self.faces(k, facet_repr=True))
            dic[k] = tuple((faces[0][i], faces[1][i])
                           for i in range(f_vector[k+1]))
        edges0 = ()
        if dimension >= 0:
            edges0 = tuple((i[0], i[1] + 1) for i in
                           self.incidences(-1, range_f_vector[1]))
            # we must take care of the case, when there are no vertices

        edges = edges0 + tuple((i[0] + sum(f_vector[:k+1]),
                                i[1] + sum(f_vector[:k+2]))
                               for k in range_f_vector1
                               for i in self.incidences(k, k+1))
        all_faces = tuple(i for k in range_f_vector for i in dic[k])
        if vertices and facets:
            def f(i): return all_faces[i]
        elif vertices:
            def f(i): return all_faces[i][0]
        elif facets:
            def f(i): return all_faces[i][1]
        else:
            def f(i): return i
        V = tuple(range(sum(f_vector)))
        D = DiGraph([V, edges], format='vertices_and_edges')
        D.relabel(f)
        return FiniteLatticePoset(D)

    def flag(self, *flag):
        r"""
        Returns the number of flags of given type.

        flag(i) is equivalent to f_vector(i).

        flag(i_1,...,i_n) will count the number of tuples
        (face_1,...,face_n), where each face_j is an i_j face and
        face_1 is contained in face_2 is contained in face_3 ...

        The implementation sorts the input arguments.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(8,15)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.flag(1,3,6)
            378000

            sage: P = polytopes.permutahedron(6)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.flag(1,4)
            7200
            sage: C.flag(1,4,4)
            7200
            sage: C.flag(1,3,4,4)
            21600
            sage: C.flag(1,4,3)
            21600
            sage: C.flag(1,2,3,4)
            43200
        """
        for number in flag:
            if not isinstance(number, Integer):
                return TypeError("All arguments of \
                    combinatorialPolyhedron.flag() must be integers.")
        dim = self.dimension()
        flag = set(number for number in flag if number in range(-1, dim+1))
        if flag == set():
            return 0
        if flag <= set([-1, dim]):
            return 1
        if self.is_trivial == 1:
            return 0
        if self.is_trivial == 2:
            if set(flag) <= set([-1, dim-1, dim]):
                return 1
            return 0
        cdef array.array flagarray = \
            array.array('I', sorted(number for number in flag
                                    if number in range(0, dim)))
        return Integer(get_flag(self._C, flagarray.data.as_uints,
                                len(flagarray)))


# Error checking on intput!
# check for containments, shouldn't take long but is very nice to the user


# add k-simple, k-simplicial
# Example:
# sage: for i in Combinations(6,3):
# ....:     x.append(list(Integer(j in i) for j in range(6)))
# P = Polyhedron(vertices=x)
# this is 2-simplicial and 6-2 simple 6-1 dimensional polyhedron
# taken from lecture notes Guenter M. Ziegler
