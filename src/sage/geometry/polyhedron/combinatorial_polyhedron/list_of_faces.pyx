# distutils: depends = sage/geometry/polyhedron/combinatorial_polyhedron/bit_vector_operations.cc
# distutils: language = c++
# distutils: extra_compile_args = -std=c++11

r"""
List of faces

This module provides a class to store faces of a polyhedron in Bit-representation.

This class allocates memory to store the faces in.
A face will be stored as vertex-incidences, where each Bit represents an incidence.
In :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.conversions` there a methods to actually convert facets of a polyhedron
to bit-representations of vertices stored in :class:`ListOfFaces`.

Moreover, :class:`ListOfFaces` calculates the dimension of a polyhedron, assuming the
faces are the facets of this polyhedron.

Each face is stored over-aligned according to the ``chunktype``.

.. SEEALSO::

    :mod:`sage.geometry.polyhedron.combinatorial_polyhedron.base`.

EXAMPLES:

Provide enough space to store `20` faces as incidences to `60` vertices::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
    ....: import ListOfFaces
    sage: face_list = ListOfFaces(20, 60)

Obtain the facets of a polyhedron::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import incidence_matrix_to_bit_rep_of_facets
    sage: P = polytopes.cube()
    sage: face_list = incidence_matrix_to_bit_rep_of_facets(P.incidence_matrix())
    sage: face_list = incidence_matrix_to_bit_rep_of_facets(P.incidence_matrix())
    sage: face_list.compute_dimension()
    3

Obtain the Vrepresentation of a polyhedron as facet-incidences::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import incidence_matrix_to_bit_rep_of_Vrep
    sage: P = polytopes.associahedron(['A',3])
    sage: face_list = incidence_matrix_to_bit_rep_of_Vrep(P.incidence_matrix())
    sage: face_list.compute_dimension()
    3

Obtain the facets of a polyhedron as :class:`ListOfFaces` from a facet list::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import facets_tuple_to_bit_rep_of_facets
    sage: facets = ((0,1,2), (0,1,3), (0,2,3), (1,2,3))
    sage: face_list = facets_tuple_to_bit_rep_of_facets(facets, 4)

Likewise for the Vrepresenatives as facet-incidences::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import facets_tuple_to_bit_rep_of_Vrep
    sage: facets = ((0,1,2), (0,1,3), (0,2,3), (1,2,3))
    sage: face_list = facets_tuple_to_bit_rep_of_Vrep(facets, 4)

Obtain the matrix of a list of faces::

    sage: face_list.matrix()
    [1 1 1 0]
    [1 1 0 1]
    [1 0 1 1]
    [0 1 1 1]

.. SEEALSO::

    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.base`,
    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator`,
    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.conversions`,
    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.polyhedron_faces_lattice`.

AUTHOR:

- Jonathan Kliem (2019-04)
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

from sage.structure.element import is_Matrix

from cysignals.signals      cimport sig_on, sig_off
from libc.string            cimport memcpy, memset
from .conversions           cimport vertex_to_bit_dictionary
from sage.matrix.matrix_integer_dense  cimport Matrix_integer_dense

cdef extern from "bit_vector_operations.cc":
    # Any Bit-representation is assumed to be `chunksize`-Bit aligned.
    cdef const size_t chunksize

    cdef inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length)
#       Return ``A & ~B == 0``.
#       A is not subset of B, iff there is a vertex in A, which is not in B.
#       ``face_length`` is the length of A and B in terms of uint64_t.

    cdef size_t get_next_level(
        uint64_t **faces, const size_t n_faces, uint64_t **nextfaces,
        uint64_t **nextfaces2, uint64_t **visited_all,
        size_t n_visited_all, size_t face_length)
#        Set ``newfaces`` to be the facets of ``faces[n_faces -1]``
#        that are not contained in a face of ``visited_all``.

#        INPUT:

#        - ``maybe_newfaces`` -- quasi of type ``uint64_t[n_faces -1][face_length]``,
#          needs to be ``chunksize``-Bit aligned
#        - ``newfaces`` -- quasi of type ``*uint64_t[n_faces -1]
#        - ``visited_all`` -- quasi of type ``*uint64_t[n_visited_all]
#        - ``face_length`` -- length of the faces

#        OUTPUT:

#        - return number of ``newfaces``
#        - set ``newfaces`` to point to the new faces

#        ALGORITHM:

#        To get all facets of ``faces[n_faces-1]``, we would have to:
#        - Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
#        - Add all the intersection of ``visited_all`` with the last face
#        - Out of both the inclusion-maximal ones are of codimension 1, i.e. facets.

#        As we have visited all faces of ``visited_all``, we alter the algorithm
#        to not revisit:
#        Step 1: Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
#        Step 2: Out of thosse the inclusion-maximal ones are some of the facets.
#                At least we obtain all of those, that we have not already visited.
#                Maybe, we get some more.
#        Step 3: Only keep those that we have not already visited.
#                We obtain exactly the facets of ``faces[n_faces-1]`` that we have
#                not visited yet.

    cdef size_t count_atoms(uint64_t *A, size_t face_length)
#        Return the number of atoms/vertices in A.
#        This is the number of set bits in A.
#        ``face_length`` is the length of A in terms of uint64_t.

    cdef int is_contained_in_one(uint64_t *face, uint64_t **faces, size_t n_faces, size_t face_length, int* skip)
#       Return whether ``face`` is contained in one of ``faces``.
#
#       Skip all ``faces[i]`` where ``skip[i]`` is ``True``.

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef class ListOfFaces:
    r"""
    A class to store the Bit-representation of faces in.

    This class will allocate the memory for the faces.

    INPUT:

    - ``n_faces`` -- the number of faces to be stored
    - ``n_atoms`` -- the total number of atoms the faces contain

    .. SEEALSO::

        :meth:`incidence_matrix_to_bit_rep_of_facets`,
        :meth:`incidence_matrix_to_bit_rep_of_Vrep`,
        :meth:`facets_tuple_to_bit_rep_of_facets`,
        :meth:`facets_tuple_to_bit_rep_of_Vrep`,
        :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator`,
        :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
        ....:     import ListOfFaces
        sage: facets = ListOfFaces(5, 13)
        sage: facets.face_length in (1, 2, 4)
        True
        sage: facets.n_atoms
        13
        sage: facets.n_faces
        5
    """
    def __init__(self, size_t n_faces, size_t n_atoms):
        r"""
        Initialize :class:`ListOfFaces`.

        See :class:`ListOfFaces`.

        TESTS:

        Checking for correct alignment of the data::

            sage: cython('''
            ....: from libc.stdint cimport uint64_t
            ....: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
            ....: cimport ListOfFaces
            ....:
            ....: cdef ListOfFaces facets
            ....: cdef size_t address
            ....: cdef size_t required_alignment
            ....:
            ....: facets = ListOfFaces(10, 13)
            ....: required_alignment = facets.face_length*8
            ....: for i in range(10):
            ....:     address = <size_t> facets.data[i]
            ....:     if not address == address & ~(required_alignment - 1):
            ....:         print('Alignment not correct')
            ....: ''')

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces).run()
        """
        self.n_faces = n_faces
        self.n_atoms = n_atoms
        self._mem = MemoryAllocator()

        # ``data`` will point to the faces as ``*uint64_t``.
        self.data = <uint64_t **> self._mem.allocarray(n_faces, sizeof(uint64_t *))

        # ``face_length`` is the length in terms of ``uint64_t``
        # NOTE: This needs to be divisible by 2, if chunksize is 128
        #       and divisible by 4, if chunksize is 256.
        self.face_length = ((n_atoms - 1)//chunksize + 1)*chunksize//64


        cdef size_t i
        for i in range(n_faces):
            # Allocate the memory for the i-th face.
            # We must allocate the memory for ListOfFaces overaligned:
            # - must be 16-byte aligned if chunksize = 128
            # - must be 32-byte aligned if chunksize = 256
            self.data[i] = <uint64_t *> \
                self._mem.aligned_malloc(chunksize//8, self.face_length*8)

    cpdef ListOfFaces __copy__(self):
        r"""
        Return a copy of self.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
            ....:     import ListOfFaces
            sage: facets = ListOfFaces(5, 13)
            sage: copy(facets).n_atoms
            13
            sage: copy(facets).n_faces
            5

        TESTS::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:     import facets_tuple_to_bit_rep_of_facets
            sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
            ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            sage: facets = facets_tuple_to_bit_rep_of_facets(bi_pyr, 6)
            sage: facets.compute_dimension()
            3
            sage: copy(facets).compute_dimension()
            3
            sage: facets.matrix() == copy(facets).matrix()
            True
            sage: copy(facets) is facets
            False
        """
        cdef ListOfFaces copy = ListOfFaces(self.n_faces, self.n_atoms)
        cdef size_t i
        for i in range(self.n_faces):
            memcpy(copy.data[i], self.data[i], self.face_length*8)
        return copy

    cpdef int compute_dimension(self) except -2:
        r"""
        Compute the dimension of a polyhedron by its facets.

        This assumes that ``self`` is the list of facets of a polyhedron.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:     import facets_tuple_to_bit_rep_of_facets, \
            ....:            facets_tuple_to_bit_rep_of_Vrep
            sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
            ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            sage: facets = facets_tuple_to_bit_rep_of_facets(bi_pyr, 6)
            sage: Vrep = facets_tuple_to_bit_rep_of_Vrep(bi_pyr, 6)
            sage: facets.compute_dimension()
            3
            sage: Vrep.compute_dimension()
            3

        ALGORITHM:

        This is done by iteration:

        Computes the facets of one of the facets (i.e. the ridges contained in
        one of the facets). Then computes the dimension of the facet, by
        considering its facets.

        Repeats until a face has only one facet. Usually this is a vertex.

        However, in the unbounded case, this might be different. The face with only
        one facet might be a ray or a line. So the correct dimension of a
        polyhedron with one facet is the number of ``[lines, rays, vertices]``
        that the facet contains.

        Hence, we know the dimension of a face, which has only one facet and
        iteratively we know the dimension of entire polyhedron we started from.

        TESTS::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
            ....:     import ListOfFaces
            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:     import incidence_matrix_to_bit_rep_of_facets, \
            ....:            incidence_matrix_to_bit_rep_of_Vrep
            sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
            ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            sage: for _ in range(10):
            ....:     points = tuple(tuple(randint(-1000,1000) for _ in range(10))
            ....:                    for _ in range(randint(3,15)))
            ....:     P = Polyhedron(vertices=points)
            ....:     inc = P.incidence_matrix()
            ....:     mod_inc = inc.delete_columns([i for i,V in enumerate(P.Hrepresentation()) if V.is_equation()])
            ....:     facets = incidence_matrix_to_bit_rep_of_facets(mod_inc)
            ....:     vertices = incidence_matrix_to_bit_rep_of_Vrep(mod_inc)
            ....:     d1 = P.dimension()
            ....:     if d1 == 0:
            ....:         continue
            ....:     d2 = facets.compute_dimension()
            ....:     d3 = vertices.compute_dimension()
            ....:     if not d1 == d2 == d3:
            ....:         print('calculation_dimension() seems to be incorrect')
        """
        if self.n_faces == 0:
            raise TypeError("at least one face needed")
        return self.compute_dimension_loop(self.data, self.n_faces, self.face_length)

    cdef int compute_dimension_loop(self, uint64_t **faces, size_t n_faces,
                                      size_t face_length) except -2:
        r"""
        Compute the dimension of a polyhedron by its facets.

        INPUT:

        - ``faces`` -- facets in Bit-representation
        - ``n_faces`` -- length of facesdata
        - ``face_length`` -- the length of each face in terms of ``uint64_t``

        OUTPUT:

        - dimension of the polyhedron

        .. SEEALSO::

            :meth:`compute_dimension`
        """
        if n_faces == 0:
            raise TypeError("wrong usage of ``compute_dimension_loop``,\n" +
                            "at least one face needed.")

        if n_faces == 1:
            # We expect the face to be the empty polyhedron.
            # Possibly it contains more than one vertex/rays/lines.
            # The dimension of a polyhedron with this face as only facet is
            # the number of atoms it contains.
            return count_atoms(faces[0], face_length)

        # ``maybe_newfaces`` are all intersection of ``faces[n_faces -1]`` with previous faces.
        # It needs to be allocated to store those faces.
        cdef ListOfFaces maybe_newfaces_mem = ListOfFaces(n_faces, face_length*64)
        cdef uint64_t **maybe_newfaces = maybe_newfaces_mem.data

        # ``newfaces`` point to the actual facets of ``faces[n_faces -1]``.
        cdef MemoryAllocator newfaces_mem = MemoryAllocator()
        cdef uint64_t **newfaces = <uint64_t **> newfaces_mem.allocarray(n_faces, sizeof(uint64_t *))

        # Calculating ``maybe_newfaces`` and ``newfaces``
        # such that ``newfaces`` points to all facets of ``faces[n_faces -1]``.
        cdef size_t new_n_faces
        sig_on()
        new_n_faces = get_next_level(faces, n_faces, maybe_newfaces,
                                      newfaces, NULL, 0, face_length)
        sig_off()

        # compute the dimension of the polyhedron,
        # by calculating dimension of one of its faces.
        return self.compute_dimension_loop(newfaces, new_n_faces, face_length) + 1

    def _test_trailing_bits(self):
        """
        Check that all trailing bits are zero.

        TESTS::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:     import facets_tuple_to_bit_rep_of_facets, \
            ....:            facets_tuple_to_bit_rep_of_Vrep
            sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
            ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            sage: facets = facets_tuple_to_bit_rep_of_facets(bi_pyr, 6)
            sage: facets._test_trailing_bits()
        """
        M = self.matrix()
        for i in range(self.n_faces):
            assert sum(M.row(i)) == count_atoms(self.data[i], self.face_length)

    cpdef ListOfFaces pyramid(self):
        r"""
        Return the list of faces of the pyramid.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:         import facets_tuple_to_bit_rep_of_facets
            sage: facets = ((0,1,2), (0,1,3), (0,2,3), (1,2,3))
            sage: face_list = facets_tuple_to_bit_rep_of_facets(facets, 4)
            sage: face_list.matrix()
            [1 1 1 0]
            [1 1 0 1]
            [1 0 1 1]
            [0 1 1 1]
            sage: face_list.pyramid().matrix()
            [1 1 1 0 1]
            [1 1 0 1 1]
            [1 0 1 1 1]
            [0 1 1 1 1]
            [1 1 1 1 0]

        Incorrect facets that illustrate how this method works::

            sage: facets = ((0,1,2,3), (0,1,2,3), (0,1,2,3), (0,1,2,3))
            sage: face_list = facets_tuple_to_bit_rep_of_facets(facets, 4)
            sage: face_list.matrix()
            [1 1 1 1]
            [1 1 1 1]
            [1 1 1 1]
            [1 1 1 1]
            sage: face_list.pyramid().matrix()
            [1 1 1 1 1]
            [1 1 1 1 1]
            [1 1 1 1 1]
            [1 1 1 1 1]
            [1 1 1 1 0]

        ::

            sage: facets = ((), (), (), ())
            sage: face_list = facets_tuple_to_bit_rep_of_facets(facets, 4)
            sage: face_list.matrix()
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            sage: face_list.pyramid().matrix()
            [0 0 0 0 1]
            [0 0 0 0 1]
            [0 0 0 0 1]
            [0 0 0 0 1]
            [1 1 1 1 0]
        """
        cdef ListOfFaces copy
        cdef size_t i, j

        # ``copy`` has a new atom and a new coatom.
        copy = ListOfFaces(self.n_faces + 1, self.n_atoms + 1)

        for i in range(self.n_faces):
            for j in range(self.face_length, copy.face_length):
                copy.data[i][j] = 0

            # All old coatoms contain their respective old atoms.
            # Also all of them contain the new atom.
            memcpy(copy.data[i], self.data[i], self.face_length*8)
            copy.data[i][self.n_atoms//64] |= vertex_to_bit_dictionary(self.n_atoms % 64)

        # The new coatom contains all atoms, but the new atom.
        for i in range(copy.face_length):
            copy.data[self.n_faces][i] = 0
        for i in range(self.n_atoms):
            copy.data[self.n_faces][i//64] |= vertex_to_bit_dictionary(i % 64)

        return copy

    cdef ListOfFaces delete_atoms_unsafe(self, uint64_t *face, int *delete):
        r"""
        Return a copy of ``self`` with each bit not set in ``face`` removed.

        Alternatively, bit ``i`` is removed if and only if ``delete[i]``.

        The bits are removed NOT unset. Thus the output will usually have
        less atoms.

        .. WARNING::

            ``face`` is assumed to be of length ``self.face_length`` or NULL.
            ``delete`` is assumed to be of length ``self.n_atoms`` or NULL.

            Exactly one of ``face`` or ``delete`` must be non-NULL.
        """

        cdef output_n_atoms
        cdef size_t i, j
        if face is not NULL:
            output_n_atoms = count_atoms(face, self.face_length)
        else:
            output_n_atoms = self.n_atoms
            for i in range(self.n_atoms):
                if delete[i]:
                    output_n_atoms -= 1
        cdef ListOfFaces output = ListOfFaces(self.n_faces, output_n_atoms)

        cdef size_t counter = 0
        cdef size_t self_pos, output_pos, self_bit, output_bit
        cdef uint64_t new_value
        cdef bint prev_was_removed = False

        # We lazily shift ``self`` to the left to copy it.
        for i in range(self.n_atoms):
            if ((face is not NULL and not face[i//64] & vertex_to_bit_dictionary(i % 64)) or
                    (delete is not NULL and delete[i])):
                # The atom will be removed.
                prev_was_removed = True
            elif (counter % 64 == 0) or prev_was_removed:
                # Either the previous atom was removed
                # or we are in the next ``uint64_t``.
                self_pos = i//64
                self_bit = i % 64
                output_pos = counter//64
                output_bit = counter % 64

                for j in range(self.n_faces):
                    # We set output.data[j][output_pos]
                    # as if no further atoms are removed.

                    # Set ``new_value`` to contain the bits
                    # i,...,i+63 of self.data[j].

                    new_value = (self.data[j][self_pos] << self_bit)
                    if self_bit and self_pos < self.face_length -1:
                        new_value |= (self.data[j][self_pos + 1] >> (64-self_bit))

                    if not output_bit:
                        output.data[j][output_pos] = new_value
                    else:
                        # Clear all bits after and including ``output_bit``.
                        output.data[j][output_pos] &= ((<uint64_t> -1) << (64 - output_bit))

                        # Copy ``new_value`` to the bits after and including ``output_bit``.
                        output.data[j][output_pos] |= (new_value >> output_bit)

                prev_was_removed = False
                counter += 1
            else:
                counter += 1

        # Finally unset all bits after ``output.n_atoms``.

        output_pos = output.n_atoms//64
        output_bit = output.n_atoms % 64

        if output_bit:
            for j in range(self.n_faces):
                # Clear all bits after and including ``output_bit``.
                output.data[j][output_pos] &= ((<uint64_t> -1) << (64 - output_bit))

        for i in range((output.n_atoms -1)//64 + 1, output.face_length):
            for j in range(self.n_faces):
                output.data[j][i] = 0

        return output

    cdef void delete_faces_unsafe(self, uint64_t *face, int *delete):
        r"""
        Deletes face ``i`` if and only if ``delete[i]``.

        Alternatively, deletes all faces such that the ``i``-th bit in ``face`` is not set.

        This will modify ``self``.

        .. WARNING::

            ``face`` is assumed to contain ``self.n_faces`` atoms or NULL.
            ``delete`` is assumed to be of length ``self.n_faces`` or NULL.

            Exactly one of ``face`` or ``delete`` must be non-NULL.
        """
        cdef size_t n_newfaces = 0
        for i in range(self.n_faces):
            if ((delete is not NULL and not delete[i]) or
                    (face is not NULL and face[i//64] & vertex_to_bit_dictionary(i % 64))):
                self.data[n_newfaces] = self.data[i]
                n_newfaces += 1

        self.n_faces = n_newfaces

    cdef void get_not_inclusion_maximal_unsafe(self, int *not_inclusion_maximal):
        r"""
        Get all faces that are not inclusion maximal.

        Set ``not_inclusion_maximal[i]`` to one if ``self.data[i]`` is not
        an inclusion-maximal face, otherwise to zero.

        If there are duplicates, all but the last duplicate will be marked as
        not inclusion maximal.

        .. WARNING::

            ``not_inclusion_maximal`` is assumed to be at least of length ``self.n_atoms`` or NULL.
        """
        cdef size_t i
        memset(not_inclusion_maximal, 0, sizeof(int)*self.n_faces)
        for i in range(self.n_faces):
            not_inclusion_maximal[i] = 1  # mark to not check inclusion of the face itself
            not_inclusion_maximal[i] = is_contained_in_one(self.data[i], self.data, self.n_faces, self.face_length, not_inclusion_maximal)

    cdef void get_faces_all_set_unsafe(self, int *all_set):
        r"""
        Get the faces that have all ``bits`` set.

        Set ``all_set[i]`` to one if ``self.data[i]``
        has all bits set, otherwise to zero.

        .. WARNING::

            ``all_set`` is assumed to be at least of length ``self.n_atoms`` or NULL.
        """
        cdef size_t i, j
        cdef uint64_t minone = <uint64_t> -1
        for i in range(self.n_faces):
            for j in range(self.n_atoms//64):
                if self.data[i][j] != minone:
                    all_set[i] = 0
                    break
            else:
                if self.n_atoms % 64 == 0:
                    all_set[i] = 1
                elif count_atoms(self.data[i] + self.n_atoms//64, 1) == self.n_atoms % 64:
                    all_set[i] = 1
                else:
                    all_set[i] = 0

    def matrix(self):
        r"""
        Obtain the matrix of self.

        Each row represents a face and each column an atom.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:     import facets_tuple_to_bit_rep_of_facets, \
            ....:     facets_tuple_to_bit_rep_of_Vrep
            sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4), (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            sage: facets = facets_tuple_to_bit_rep_of_facets(bi_pyr, 6)
            sage: Vrep = facets_tuple_to_bit_rep_of_Vrep(bi_pyr, 6)
            sage: facets.matrix()
            [1 1 0 0 1 0]
            [0 1 1 0 1 0]
            [0 0 1 1 1 0]
            [1 0 0 1 1 0]
            [1 1 0 0 0 1]
            [0 1 1 0 0 1]
            [0 0 1 1 0 1]
            [1 0 0 1 0 1]
            sage: facets.matrix().transpose() == Vrep.matrix()
            True
        """
        from sage.rings.all import ZZ
        from sage.matrix.constructor import matrix
        cdef Matrix_integer_dense M = matrix(
                ZZ, self.n_faces, self.n_atoms, 0)

        cdef size_t i,j
        for i in range(self.n_faces):
            for j in range(self.n_atoms):
                if self.data[i][j//64] & vertex_to_bit_dictionary(j % 64):
                    M.set_unsafe_si(i, j, 1)

        M.set_immutable()
        return M

cdef tuple face_as_combinatorial_polyhedron(ListOfFaces facets, ListOfFaces Vrep, uint64_t *face, uint64_t *coface):
    r"""
    Obtain facets and Vrepresentation of ``face`` as new combinatorial polyhedron.

    INPUT:

    - ``facets`` -- facets of the polyhedron
    - ``Vrep`` -- Vrepresentation of the polyhedron
    - ``face`` -- face in Vrepresentation or ``NULL``
    - ``coface`` -- face in facet-represetnation or ``NULL``

    Exactly one of ``face`` or ``coface`` must be ``NULL``.

    OUTPUT: A tuple of new facets and new Vrepresentation as :class:`ListOfFaces`.

    .. WARNING::

        ``face`` is assumed to be of length ``old_facets.face_length`` or NULL.
        ``coface`` is assumed to be of length ``old_Vrep.face_length`` or NULL.
    """
    cdef ListOfFaces new_facets, new_Vrep
    cdef int* delete
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef size_t i

    #assert facets.matrix() == Vrep.matrix().transpose()
    #print(facets.matrix(), '1')

    # Delete all atoms not in the face.
    if face is not NULL:
        new_facets = facets.delete_atoms_unsafe(face, NULL)
        new_Vrep = Vrep.__copy__()
        new_Vrep.delete_faces_unsafe(face, NULL)

        delete = <int*> mem.allocarray(new_facets.n_faces, sizeof(int))
    else:
        delete = <int*> mem.allocarray(max(facets.n_faces, facets.n_atoms), sizeof(int))

        # Set ``delete[i]`` to one if ``i`` is not an vertex of ``coface``.
        for i in range(Vrep.n_faces):
            if is_subset(coface, Vrep.data[i], Vrep.face_length):
                delete[i] = 0
            else:
                delete[i] = 1

        new_facets = facets.delete_atoms_unsafe(NULL, delete)
        new_Vrep = Vrep.__copy__()
        new_Vrep.delete_faces_unsafe(NULL, delete)

    #assert new_facets.matrix() == new_Vrep.matrix().transpose()
    #print(new_facets.matrix(), '2')

    # Delete all facets that define the face.
    new_facets.get_faces_all_set_unsafe(delete)
    new_facets.delete_faces_unsafe(NULL, delete)
    new_Vrep = new_Vrep.delete_atoms_unsafe(NULL, delete)

    #assert new_facets.matrix() == new_Vrep.matrix().transpose()
    #print(new_facets.matrix(), '3')

    # Now delete all facets that are not inclusion maximal.
    # the last copy of each duplicate will remain.
    new_facets.get_not_inclusion_maximal_unsafe(delete)
    new_facets.delete_faces_unsafe(NULL, delete)
    new_Vrep = new_Vrep.delete_atoms_unsafe(NULL, delete)

    #assert new_facets.matrix() == new_Vrep.matrix().transpose()
    #print(new_facets.matrix(), '4')

    return (new_facets, new_Vrep)
