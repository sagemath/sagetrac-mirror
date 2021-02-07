"""
Inline cython methods for lists of faces.
"""
# ****************************************************************************
#       Copyright (C) 2020 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .face_data_structure             cimport *
from libc.string                      cimport memset
from cysignals.signals                cimport sig_on, sig_off

cdef struct face_list_s:
    face_t* faces
    size_t n_faces
    size_t total_n_faces
    size_t n_atoms
    size_t n_coatoms
    bint polyhedron_is_simple
    bint* is_not_new_face

ctypedef face_list_s face_list_t[1]

#############################################################################
# Face List Initalization
#############################################################################

cdef inline int face_list_init(face_list_t faces, size_t n_faces, size_t n_atoms, size_t n_coatoms, MemoryAllocator mem) except -1:
    """
    Sets the initial values for a list of faces with given number of faces
    and number of atoms.
    """
    face_list_shallow_init(faces, n_faces, n_atoms, n_coatoms, mem)
    cdef size_t i
    for i in range(n_faces):
        face_init(faces.faces[i], n_atoms, n_coatoms, mem)

cdef inline int face_list_shallow_init(face_list_t faces, size_t n_faces, size_t n_atoms, size_t n_coatoms, MemoryAllocator mem) except -1:
    """
    Initialize ``faces`` completely, but only set up memory for the pointers to the faces.
    """
    faces.n_faces = n_faces
    faces.total_n_faces = n_faces
    faces.n_atoms = n_atoms
    faces.n_coatoms = n_coatoms
    faces.faces = <face_t *> mem.allocarray(n_faces, sizeof(face_t))
    faces.is_not_new_face = <bint *> mem.allocarray(n_faces, sizeof(bint))
    faces.polyhedron_is_simple = False

cdef inline int face_list_copy(face_list_t dst, face_list_t src) except -1:
    """
    This is a deep copy. All the data for the faces is copied.

    Asserts that ``dst`` is allocated and fits everything.
    """
    assert dst.total_n_faces >= src.n_faces
    assert dst.n_atoms >= src.n_atoms
    assert dst.n_coatoms >= src.n_coatoms

    dst.n_faces = src.n_faces
    dst.polyhedron_is_simple = src.polyhedron_is_simple

    cdef size_t i
    for i in range(src.n_faces):
        face_copy(dst.faces[i], src.faces[i])

cdef inline int face_list_shallow_copy(face_list_t dst, face_list_t src) except -1:
    """
    Copy the pointers to the faces.

    Asserts that ``dst`` contains enough space for the pointers.
    """
    assert dst.total_n_faces >= src.n_faces
    dst.n_atoms = src.n_atoms
    dst.n_coatoms = src.n_coatoms
    dst.polyhedron_is_simple = src.polyhedron_is_simple
    dst.n_faces = src.n_faces

    cdef size_t i
    for i in range(src.n_faces):
        dst.faces[i] = src.faces[i]

cdef inline int add_face_shallow(face_list_t faces, face_t face) nogil except -1:
    """
    Add a face to faces.
    """
    if not faces.total_n_faces >= faces.n_faces + 1:
        with gil:
            raise AssertionError
    faces.faces[faces.n_faces][0] = face[0]
    faces.n_faces += 1

cdef inline int add_face_deep(face_list_t faces, face_t face) except -1:
    """
    Add a face to faces.
    """
    assert faces.total_n_faces >= faces.n_faces + 1
    face_copy(faces.faces[faces.n_faces], face)
    faces.n_faces += 1

#############################################################################
# Face Comparison
#############################################################################

cdef void sort_faces_list(face_list_t faces)

cdef inline size_t find_face(face_t face, face_list_t faces):
    r"""
    Return the index of ``face`` in ``faces``.

    Return ``-1`` if the ``face`` is not contained.

    .. NOTE::

        Assumes that ``faces`` are sorted.
    """
    cdef size_t start = 0
    cdef size_t middle
    cdef size_t n_faces = faces.n_faces
    cdef face_t* faces_pt = faces.faces
    cdef int val


    while (n_faces > 1):
        # In each iteration step, we will look for ``face`` in
        # ``faces_pt[start:start+n_faces]``.
        middle = n_faces//2
        val = face_cmp(face, faces_pt[middle + start])
        if val < 0:
            # If face is in the list, then in the lower half.
            # Look for face in ``faces[start : start + middle]`` in next step.
            n_faces = middle
        elif val > 0:
            # If face is in the list, then in the upper half.
            # Look for face in ``faces[start+middle:start+n_faces]``, i.e.
            # ``faces[start + middle : (start + middle) + n_faces - middle]``.
            n_faces -= middle
            start += middle
        else:
            return middle + start
    if face_cmp(face, faces_pt[start]) == 0:
        return start
    else:
        return -1

cdef inline bint is_contained_in_one_fused(face_t face, face_list_t faces, algorithm_variant algorithm) nogil:
    """
    Return whether ``face`` is contained in one of ``faces``.
    """
    cdef size_t i
    for i in range(faces.n_faces):
        if face_issubset_fused(face, faces.faces[i], algorithm):
            return True
    return False

cdef inline bint is_not_maximal_fused(face_list_t faces, size_t j, algorithm_variant algorithm, bint* is_not_new_face) nogil:
    """
    Return whether face ``j`` is not maximal in ``faces``.
    """
    cdef size_t i
    if algorithm_variant is standard:
        for i in range(j):
            if (not is_not_new_face[i]) and face_issubset_fused(faces.faces[j], faces.faces[i], algorithm):
                # It suffices to check those faces, that are maximal.
                # This way, if multiple identical faces are maximal,
                # exactly the last one is considered maximal.
                return True
        for i in range(j+1, faces.n_faces):
            if face_issubset_fused(faces.faces[j], faces.faces[i], algorithm):
                return True
        return False
    else:
        # For simple polytopes an intersection of facets is of codimension 2,
        # if and only if it contains a coatom.
        return face_isempty(faces.faces[j])

#############################################################################
# Arithmetic
#############################################################################

cdef inline int face_list_intersection_fused(face_list_t dest, face_list_t A, face_t b, algorithm_variant algorithm) nogil except -1:
    """
    Set ``dest`` to be the intersection of each face of ``A`` with ``b``.
    """
    if not dest.total_n_faces >= A.n_faces:
        with gil:
            raise AssertionError
    if not dest.n_atoms >= A.n_atoms:
        with gil:
            raise AssertionError
    dest.n_faces = A.n_faces
    dest.polyhedron_is_simple = A.polyhedron_is_simple

    cdef size_t i
    for i in range(A.n_faces):
        face_intersection_fused(dest.faces[i], A.faces[i], b, algorithm)


cdef inline size_t get_next_level_fused(
        face_list_t faces,
        face_list_t new_faces,
        face_list_t visited_all, algorithm_variant algorithm) nogil except -1:
    """
    Set ``new_faces`` to be the facets of ``faces.faces[face.n_faces-1]``
    that are not contained in a face of ``visited_all``.

    Reduce the number of faces in ``faces`` by one.

    INPUT:

    - ``faces`` -- containing at least one face
    - ``new_faces`` -- needs to be of same size as ``faces``
    - ``visited_all`` -- the faces which have been visited before

    OUTPUT:

    - set ``new_faces`` to point to the new faces

    ALGORITHM:

    To get all facets of ``faces.faces[faces.n_faces-1]``, we would have to:
    - Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
    - Add all the intersection of ``visited_all`` with the last face
    - Out of both the inclusion-maximal ones are of codimension one, i.e. facets.

    As we have visited all faces of ``visited_all``, we alter the algorithm
    to not revisit:
    Step 1: Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
    Step 2: Out of those the inclusion-maximal ones are some of the facets.
            At least we obtain all of those, that we have not already visited.
            Maybe, we get some more.
    Step 3: Only keep those that we have not already visited.
            We obtain exactly the facets of ``faces[n_faces-1]`` that we have
            not visited yet.
    """
    # We keep track, which face in ``new_faces`` is a new face.
    cdef size_t n_faces = faces.n_faces
    cdef bint* is_not_new_face = new_faces.is_not_new_face
    memset(is_not_new_face, 0, n_faces*sizeof(bint))

    # Step 1:
    n_faces -= 1
    faces.n_faces -= 1
    face_list_intersection_fused(new_faces, faces, faces.faces[n_faces], algorithm)

    cdef size_t j
    for j in range(n_faces):
        if (is_not_maximal_fused(new_faces, j, algorithm, is_not_new_face) or  # Step 2
                is_contained_in_one_fused(new_faces.faces[j], visited_all, algorithm)):  # Step 3
            is_not_new_face[j] = True

    # Set ``new_faces`` to point to the correct ones.
    cdef size_t n_new_faces = 0
    for j in range(n_faces):
        if is_not_new_face[j]:
            continue
        # It is a new face of codimension 1.
        # Either ``faces.n_new_faces == j`` or ``new_faces.faces[n_new_faces]`` is not
        # a new face.

        swap_faces(new_faces.faces[j], new_faces.faces[n_new_faces])

        n_new_faces += 1

    new_faces.n_faces = n_new_faces
    return n_new_faces

cdef inline size_t get_next_level(
        face_list_t faces,
        face_list_t new_faces,
        face_list_t visited_all) nogil except -1:

    cdef size_t output
    sig_on()
    if faces.polyhedron_is_simple:
        output = get_next_level_fused(faces, new_faces, visited_all, <simple> 0)
    else:
        output = get_next_level_fused(faces, new_faces, visited_all, <standard> 0)
    sig_off()
    return output

cdef inline size_t bit_rep_to_coatom_rep(face_t face, face_list_t coatoms, size_t *output):
    """
    Write the coatom-representation of face in output. Return length.
    ``face_length`` is the length of ``face`` and ``coatoms[i]``
    in terms of uint64_t.
    ``n_coatoms`` length of ``coatoms``.
    """
    cdef size_t count_length = 0
    cdef size_t i
    for i in range(coatoms.n_faces):
        if face_issubset(face, coatoms.faces[i]):
            output[count_length] = i
            count_length += 1
    return count_length

cdef inline bint face_list_check_alignment(face_list_t faces):
    """
    Return whether all faces in ``faces`` are aligned correctly.
    """
    cdef size_t i
    for i in range(faces.n_faces):
       if not face_check_alignment(faces.faces[i]):
           return False
    return True
