r"""
Face iterator for polyhedra

This iterator in principle works on every graded lattice, where
every interval of length two has exactly 4 elements (diamond property).

It also works on unbounded polyhedra, as those satisfy the diamond property,
except for intervals including the empty face.
A (slightly generalized) description of the algorithm can be found in [KS2019]_.

Terminology in this module:

- Coatoms               -- the faces from which all others are constructed in
                           the face iterator. This will be facets or Vrepr.
                           In non-dual mode, faces are constructed as
                           intersections of the facets. In dual mode, the are
                           constructed theoretically as joins of vertices.
                           The coatoms are reprsented as incidences with the
                           atoms they contain.
- Atoms                 -- facets or Vrepr depending on application of algorithm.
                           Atoms are reprsented as incidences of coatoms they
                           are contained in.

.. SEEALSO::

    :mod:`sage.geometry.polyhedron.combinatorial_polyhedron.base`.

EXAMPLES:

Construct a face iterator::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator \
    ....:         import FaceIterator
    sage: P = polytopes.octahedron()
    sage: C = CombinatorialPolyhedron(P)

    sage: FaceIterator(C, False)
    Iterator over the proper faces of a 3-dimensional combinatorial polyhedron
    sage: FaceIterator(C, False, output_dimension=2)
    Iterator over the 2-faces of a 3-dimensional combinatorial polyhedron

Iterator in the non-dual mode starts with facets::

    sage: it = FaceIterator(C, False)
    sage: [next(it) for _ in range(9)]
    [A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 1-dimensional face of a 3-dimensional combinatorial polyhedron]

Iterator in the dual-mode starts with vertices::

    sage: it = FaceIterator(C, True)
    sage: [next(it) for _ in range(7)]
    [A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 1-dimensional face of a 3-dimensional combinatorial polyhedron]

Obtain the Vrepresentation::

    sage: it = FaceIterator(C, False)
    sage: face = next(it)
    sage: face.Vrepr()
    (A vertex at (0, -1, 0), A vertex at (0, 0, -1), A vertex at (1, 0, 0))
    sage: face.length_Vrepr()
    3

Obtain the facet-representation::

    sage: it = FaceIterator(C, True)
    sage: face = next(it)
    sage: face.Hrepr()
    (An inequality (-1, -1, 1) x + 1 >= 0,
     An inequality (-1, -1, -1) x + 1 >= 0,
      An inequality (-1, 1, -1) x + 1 >= 0,
       An inequality (-1, 1, 1) x + 1 >= 0)
    sage: face.Hrepr(names=False)
    (4, 5, 6, 7)
    sage: face.length_Hrepr()
    4

In non-dual mode one can ignore all faces contained in the current face::

    sage: it = FaceIterator(C, False)
    sage: face = next(it)
    sage: face.Hrepr(names=False)
    (7,)
    sage: it.ignore_subfaces()
    sage: [face.Hrepr(names=False) for face in it]
    [(6,),
    (5,),
    (4,),
    (3,),
    (2,),
    (1,),
    (0,),
    (5, 6),
    (1, 6),
    (0, 1, 5, 6),
    (4, 5),
    (0, 5),
    (0, 3, 4, 5),
    (3, 4),
    (2, 3),
    (0, 3),
    (0, 1, 2, 3),
    (1, 2),
    (0, 1)]

In dual mode one can ignore all faces that contain the current face::

    sage: it = FaceIterator(C, True)
    sage: face = next(it)
    sage: face.Vrepr(names=False)
    (5,)
    sage: it.ignore_supfaces()
    sage: [face.Vrepr(names=False) for face in it]
    [(4,),
    (3,),
    (2,),
    (1,),
    (0,),
    (3, 4),
    (2, 4),
    (0, 4),
    (0, 3, 4),
    (0, 2, 4),
    (1, 3),
    (0, 3),
    (0, 1, 3),
    (1, 2),
    (0, 2),
    (0, 1, 2),
    (0, 1)]

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

from sage.rings.integer     cimport smallInteger
from cysignals.signals      cimport sig_check, sig_on, sig_off, sig_check_no_except
from .conversions           cimport bit_repr_to_Vrepr_list
from .base                  cimport CombinatorialPolyhedron, KunzCone
from .bit_vector_operations cimport get_next_level, \
                                    count_atoms, \
                                    bit_repr_to_coatom_repr, \
                                    is_bad_face_cc, \
                                    chunksize, \
                                    get_next_level_with_nonzero
from cysignals.memory       cimport sig_malloc, sig_realloc, sig_free, sig_calloc
from cython.parallel        cimport prange, threadid
from .check_bad_face        cimport check_bad_face, check_bad_faces

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef inline int next_dimension(iter_struct *structure) nogil except -1:
    r"""
    Set attribute ``face`` to the next face and return the dimension.

    Will return the dimension of the polyhedron on failure.

    The function calls :meth:`FaceIterator.next_face_loop` until a new
    face is set or until the iterator is consumed.

    .. NOTE::

        The face_iterator can be prevented from visiting any subfaces
        (or supfaces in dual mode) as in :meth:`FaceIterator.ignore_subfaces`
        and :meth`FaceIterator.ignore_supfaces`.

        Those methods add the current face to ``visited_all`` before
        visiting sub-/supfaces instead of after. One cannot arbitrarily
        add faces to ``visited_all``, as visited_all has a maximal length.
    """
    cdef int dim = structure[0].max_dimension
    while (not next_face_loop(structure)) and (structure[0].current_dimension < dim):
        sig_check()
    structure[0]._index += 1
    return structure[0].current_dimension

cdef inline int next_face_loop(iter_struct *structure) nogil except -1:
    r"""
    Set attribute ``face`` to the next face. On success return `1`.
    Otherwise `0`. Needs to be recalled then.

    If ``self.current_dimension == self.dimension``, then the iterator is
    consumed.
    """
    # Getting ``[faces, n_faces, n_visited_all]`` according to dimension.
    cdef uint64_t **faces = structure[0].newfaces[structure[0].current_dimension]
    cdef size_t n_faces = structure[0].n_newfaces[structure[0].current_dimension]
    cdef size_t n_visited_all = structure[0].n_visited_all[structure[0].current_dimension]

    if (structure[0].output_dimension > -2) and (structure[0].output_dimension != structure[0].current_dimension):
        # If only a specific dimension was requested (i.e. ``output_dimension > -2``),
        # then we will not yield faces in other dimension.
        structure[0].yet_to_visit = 0

    if structure[0].yet_to_visit:
        # Set ``face`` to the next face.
        structure[0].yet_to_visit -= 1
        structure[0].face = faces[structure[0].yet_to_visit]
        if structure[0].current_dimension < structure[0].dimension -1:
            structure[0].nonzero_face = structure[0].nonzero_newfaces[structure[0].current_dimension][structure[0].yet_to_visit]
        if structure[0].LHS:
            structure.current_LHS = &structure[0].LHS[structure[0].current_dimension][structure[0].yet_to_visit]
            structure.current_RHS = &structure[0].RHS[structure[0].current_dimension][structure[0].yet_to_visit]
        return 1

    if structure[0].current_dimension <= structure[0].lowest_dimension:
        # We will not yield the empty face.
        # We will not yield below requested dimension.
        structure[0].current_dimension += 1
        return 0

    if n_faces <= 1:
        # There will be no more faces from intersections.
        structure[0].current_dimension += 1
        return 0

    # We will visit the last face now.
    structure[0].n_newfaces[structure[0].current_dimension] -= 1
    n_faces -= 1

    if not structure[0].first_time[structure[0].current_dimension]:
        # In this case there exists ``faces[n_faces + 1]``, of which we
        # have visited all faces, but which was not added to
        # ``visited_all`` yet.
        structure[0].visited_all[n_visited_all] = faces[n_faces + 1]
        structure[0].n_visited_all[structure[0].current_dimension] += 1
        n_visited_all = structure[0].n_visited_all[structure[0].current_dimension]
    else:
        # Once we have visited all faces of ``faces[n_faces]``, we want
        # to add it to ``visited_all``.
        structure[0].first_time[structure[0].current_dimension] = False

    # Get the faces of codimension 1 contained in ``faces[n_faces]``,
    # which we have not yet visited.
    cdef size_t newfacescounter

    cdef uint64_t **LHS = NULL
    cdef uint64_t **RHS = NULL
    if structure[0].LHS:
        LHS = &structure[0].LHS[structure[0].current_dimension-1]
        RHS = &structure[0].RHS[structure[0].current_dimension-1]

    newfacescounter = get_next_level_with_nonzero(
        faces, n_faces + 1, structure[0].maybe_newfaces[structure[0].current_dimension-1],
        structure[0].newfaces[structure[0].current_dimension-1],
        structure[0].visited_all, n_visited_all, structure[0].face_length, structure[0].is_not_newface, LHS, RHS,
        structure[0].nonzero_maybe_newfaces[structure[0].current_dimension-1],
        structure[0].nonzero_newfaces[structure[0].current_dimension-1])

    if newfacescounter:
        # ``faces[n_faces]`` contains new faces.
        # We will visted them on next call, starting with codimension 1.

        # Setting the variables correclty for next call of ``next_face_loop``.
        structure[0].current_dimension -= 1
        structure[0].first_time[structure[0].current_dimension] = True
        structure[0].n_newfaces[structure[0].current_dimension] = newfacescounter
        structure[0].n_visited_all[structure[0].current_dimension] = n_visited_all
        structure[0].yet_to_visit = newfacescounter
        return 0
    else:
        # ``faces[n_faces]`` contains no new faces.
        # Hence there is no need to add it to ``visited_all``.
        # NOTE:
        #     For the methods ``ignore_subfaces`` and ``ignore_supfaces``
        #     this step needs to be done, as ``faces[n_faces]`` might
        #     have been added manually to ``visited_all``.
        #     So this step is required to respect boundaries of ``visited_all``.
        structure[0].first_time[structure[0].current_dimension] = True
        return 0

cdef inline size_t myPow(size_t x, size_t p) nogil:
    if (p == 0): return 1
    if (p == 1):  return x

    cdef size_t tmp = myPow(x, p/2)
    if (p%2 == 0) : return tmp * tmp
    else:  return x * tmp * tmp

cdef inline int prepare_partial_iter(iter_struct *face_iter, size_t job_ID, size_t *f_vector, size_t rec_depth, bint bad) nogil except -1:
    r"""
    Prepares the face iterator to not visit the first faces of codimension ``rec_depth``.

    More precisly, we divide the job into parts as:
    Visit the `i_1`-th facet, of that the `i_2` facet, of that the `i_3` facet.
    The ``job_ID`` of that task `i_3 + i_2*r + i_1*r^2`, where `r` denotes recursion depth.
    """
    cdef size_t j, n_visited_all, n_newfaces
    cdef int prep_dim

    # Figuring out the facet by the job_ID.
    cdef size_t facet_nr = 0
    if (rec_depth > 0):
        facet_nr = job_ID/(myPow(face_iter[0].n_coatoms, (rec_depth - 1)))

    face_iter[0].face = NULL

    cdef size_t codim # We do recursion on codimension.

    cdef int d
    cdef int dimension = face_iter[0].dimension
    cdef bint only_preparing_last_level = True

    if((facet_nr != face_iter[0].current_stadium[0])):
        # The face iterator is in a completely different state.
        # We start from scratch.
        # Basically we initialize the face iterator again.
        face_iter[0].current_dimension = dimension -1
        face_iter[0].lowest_dimension = 0

        face_iter[0].n_visited_all[dimension -1] = 0
        if not face_iter[0].bounded:
            # Adding the far face as visited_all.
            face_iter[0].n_visited_all[dimension - 1] = 1

        face_iter[0].n_newfaces[dimension - 1] = face_iter[0].n_coatoms
        face_iter[0].current_stadium[0] = 0
        face_iter[0].first_time[dimension - 1] = 1
        only_preparing_last_level = False

    cdef int codim2
    cdef size_t missing_faces
    cdef size_t old_number

    for codim in range(rec_depth):
        codim2 = codim
        prep_dim = dimension - codim2 - 1

        # facet_nr and job_ID will go down recursively.
        facet_nr = job_ID/(myPow(face_iter[0].n_coatoms, (rec_depth - codim - 1)))
        job_ID = job_ID%(myPow(face_iter[0].n_coatoms, (rec_depth - codim - 1)))

        if((facet_nr != face_iter[0].current_stadium[codim]) \
                or (face_iter[0].current_dimension >= prep_dim) \
                or (codim == rec_depth - 1)):
            # In case the iterator in wrong state with respect to the current codim
            # OR the iterator was already in the wrong state with respect to higher codim
            # OR we are in the last instance of our recursion and must setup the iterator properly.

            if((codim == rec_depth - 1) and (only_preparing_last_level) and (codim > 0)):
                # We reset the iterator to the state, its claiming to be in.
                # More precisly, the iterator has visited the last of the newfaces
                # and reduced the corresponding number.
                # We undo this.
                face_iter[0].n_newfaces[prep_dim] += 1

            if((face_iter[0].current_dimension > prep_dim) \
                    or (facet_nr >= face_iter[0].n_newfaces[prep_dim] + face_iter[0].current_stadium[codim])):
                # If the dimension of the face iterator has gone up,
                # this is because there where no faces left in that dimension.
                # OR
                # The ``facet_nr`` is larger, than the number of newfaces.
                return 0

            face_iter[0].current_dimension = prep_dim
            if job_ID == 0 and not bad:
                # This is the first time we visit this face
                # and we increase the f_vector accordingly.

                # In the bad_vector case, we allow only recursion depth,
                # such that we know, that this face has high embedding dimension
                # and is therefore not bad.
                f_vector[prep_dim + 1] += 1

            # We mark faces as visited, to match the requested statium.
            n_visited_all = face_iter[0].n_visited_all[prep_dim]
            n_newfaces    = face_iter[0].n_newfaces[prep_dim]
            missing_faces = facet_nr - face_iter[0].current_stadium[codim]
            for j in range(n_newfaces - missing_faces, n_newfaces):
               face_iter[0].visited_all[n_visited_all] \
                       = face_iter[0].newfaces[prep_dim][j]
               n_visited_all += 1

            face_iter[0].n_visited_all[prep_dim]  = n_visited_all
            face_iter[0].current_stadium[codim]   = facet_nr
            face_iter[0].n_newfaces[prep_dim]    -= missing_faces

            # The face iterator is set up to
            # visit the first face of codimension `codim + 1`.
            face_iter[0].current_stadium[codim+1] = 0
            face_iter[0].first_time[prep_dim] = 1
            face_iter[0].yet_to_visit = 0

            # Make sure the iterator stops, once it has
            # done its part.
            face_iter[0].max_dimension = prep_dim
            only_preparing_last_level = False

            if codim < rec_depth -1:
                # In case we are not done preparing.

                # We keep the last face in ``newfaces``,
                # such that the iterator can prepared for different states later.
                old_number = face_iter[0].n_newfaces[dimension-codim-1]
                d = next_dimension(face_iter)
                face_iter[0].n_newfaces[prep_dim] = old_number

                face_iter[0].first_time[prep_dim] = 1
                face_iter[0].yet_to_visit = 0

    return 1

cdef int parallel_bad_vector(
        iter_struct **face_iter, size_t *bad_vector,
        size_t n_threads, size_t rec_depth,
        bint orbit_only,
        size_t start, size_t end, FILE **fp) except -1:
    """
    Distribute computation of the bad faces to threads.

    Sum up the partial results.
    """
    cdef size_t i
    cdef int j
    cdef MemoryAllocator mem = MemoryAllocator()

    if rec_depth > 3:
        raise ValueError("recursion depth can be at most 3 (as ridges are not checked for recursion depth 3)")

    cdef size_t **shared_bad = <size_t **> mem.calloc(n_threads, sizeof(size_t *))
    cdef int dimension = face_iter[0][0].dimension

    for i in range(n_threads):
        shared_bad[i] = <size_t *> mem.calloc(dimension + 2, sizeof(size_t))
    cdef size_t n_faces = face_iter[0][0].n_coatoms

    if myPow(n_faces, rec_depth) < end:
        end = myPow(n_faces, rec_depth)

    cdef size_t l, ID
    for l in prange(start, end, nogil=True, num_threads=n_threads, schedule='dynamic', chunksize=1):
        partial_bad(face_iter, shared_bad, l, rec_depth, orbit_only, fp)

    # Summing up the results..
    for i in range(n_threads):
        for j in range(dimension + 2):
            bad_vector[j] += shared_bad[i][j]

cdef int parallel_f_vector(iter_struct **face_iter, size_t *f_vector, size_t n_threads, size_t recursion_depth) except -1:
    cdef size_t i
    cdef int j
    cdef size_t rec_depth = recursion_depth
    #omp_set_num_threads(n_threads);
    cdef size_t **shared_f = <size_t **> sig_calloc(n_threads, sizeof(size_t *))
    cdef int dimension = face_iter[0][0].dimension

    for i in range(n_threads):
        shared_f[i] = <size_t *> sig_calloc(dimension + 2, sizeof(size_t))
    cdef size_t n_faces = face_iter[0][0].n_coatoms

    #iter_struct *my_iter;
    #size_t *my_f;
    #for l in prange(0, n_faces ** rec_depth, schedule='dynamic', chunksize=1):
    #pragma omp parallel for shared(face_iter, shared_f) schedule(dynamic, 1)
    cdef size_t l, ID
    for l in prange(myPow(n_faces, rec_depth), nogil=True, num_threads=n_threads, schedule='dynamic', chunksize=1):
        #partial_f(face_iter[openmp.omp_get_thread_num()], shared_f[openmp.omp_get_thread_num()], l)
        #partial_f(face_iter[omp_get_thread_num()], shared_f[omp_get_thread_num()], l);
        with gil:
            sig_check()
        partial_f(face_iter, shared_f, l, rec_depth)

    for i in range(n_threads):
        for j in range(dimension + 2):
            f_vector[j] += shared_f[i][j]

        sig_free(shared_f[i])

    sig_free(shared_f)

cdef inline int partial_f(iter_struct **face_iter_all, size_t **f_vector_all, size_t i, size_t rec_depth) nogil except -1:
    cdef size_t ID
    with gil:
        ID = threadid()
    cdef iter_struct * face_iter = face_iter_all[ID]
    cdef size_t * f_vector = f_vector_all[ID]
    cdef int rec_depth2
    cdef int d, dimension
    cdef size_t j
    if (prepare_partial_iter(face_iter, i, f_vector, rec_depth, False)):
        dimension = face_iter[0].dimension
        d = next_dimension(face_iter)
        rec_depth2 = rec_depth
        while (d < dimension -rec_depth2):
            f_vector[d + 1] += 1
            d = next_dimension(face_iter)

cdef inline int partial_bad(
        iter_struct **face_iter_all, size_t **bad_vector_all,
        size_t job_ID, size_t rec_depth, bint orbit_only, FILE **fps) nogil except -1:
    cdef size_t ID
    with gil:
        ID = threadid()
    cdef iter_struct * face_iter = face_iter_all[ID]
    cdef FILE * fp = fps[ID]

    # Figuring out, wether we can safely skip this, as we have visited this orbit already.
    cdef size_t n = face_iter.n_first_orbit_facets
    cdef size_t *to_do = face_iter.first_orbit_facets
    cdef size_t j
    cdef size_t n_facets = face_iter.n_coatoms

    cdef size_t facet_nr = 0
    if (rec_depth > 0):
        facet_nr = job_ID/(myPow(face_iter[0].n_coatoms, (rec_depth - 1)))

    facet_nr = n_facets - facet_nr

    # See if the face_nr corresponds to the last representative of an orbit.
    # (We start with the last facet and work our way to the beginning!)
    cdef bint leave = True
    for j in range(n):
        if to_do[j] == facet_nr:
            leave = False
    if facet_nr == n_facets:
        leave = False

    if leave and orbit_only:
        # This is not the first facet of an orbit.
        return 0

    cdef size_t * bad_vector = bad_vector_all[ID]
    cdef int rec_depth2
    cdef int d, dimension
    if (prepare_partial_iter(face_iter, job_ID, bad_vector, rec_depth, True)):
        dimension = face_iter[0].dimension
        d = next_dimension(face_iter)
        rec_depth2 = rec_depth
        while (d < dimension -rec_depth2):
            if is_bad_face(face_iter, fp):
                bad_vector[d + 1] += 1
            d = next_dimension(face_iter)

    '''
    # Now checking the bad faces so far.
    cdef size_t l = 0
    cdef int test = check_bad_faces(
            face_iter[l].PolyIneq, face_iter[l].n_coatoms, face_iter[l].dimension + 2,
            face_iter[l].bad_faces_LHS,
            face_iter[l].n_bad_faces,
            face_iter[l].bad_faces_pt)

    face_iter[0].n_bad_faces = 0
    face_iter[0].len_bad_faces = 0
    '''
    write_bad_faces(face_iter, fp)

cdef inline int is_bad_face(iter_struct *face_iter, FILE *fp) nogil except -1:
    """
    Return 1, if the current face is a bad face.

    If the corresponding polyhedron is non-empty, print a message.

    If it contains a point (e.g. counterexample to Wilf conjecture),
    raise an error.
    """
    # Obtain values from the iterator structure.
    cdef int dimension          = face_iter[0].dimension
    cdef uint64_t *face         = face_iter[0].face
    cdef uint64_t ** coatoms    = face_iter[0].newfaces[dimension-1]
    cdef size_t n_coatoms       = face_iter[0].n_coatoms
    cdef size_t face_length     = face_iter[0].face_length
    cdef uint64_t *LHS          = face_iter[0].LHS[dimension-1]
    cdef uint64_t *RHS          = face_iter[0].RHS[dimension-1]
    cdef uint32_t *nonzero_face = face_iter[0].nonzero_face
    cdef uint8_t *bad_faces     = face_iter[0].bad_faces
    cdef size_t len_bad_faces   = face_iter[0].len_bad_faces
    cdef uint64_t *bad_faces_LHS= face_iter[0].bad_faces_LHS
    cdef size_t n_bad_faces     = face_iter[0].n_bad_faces

    cdef size_t e  # The value for e is stored here, as ``check_bad_face`` needs it.

    cdef size_t output
    output = is_bad_face_cc(face, nonzero_face, dimension, coatoms, n_coatoms,
                            face_length, LHS, RHS, face_iter[0].current_LHS, face_iter[0].current_RHS,
                            &e, bad_faces + len_bad_faces)
    if output == 0:
        return 0

    face_iter[0].bad_faces_pt[n_bad_faces] = bad_faces + len_bad_faces
    face_iter[0].len_bad_faces += output + 2
    bad_faces_LHS[n_bad_faces] = face_iter[0].current_LHS[0]
    face_iter[0].n_bad_faces += 1

    if unlikely(face_iter[0].len_bad_faces + 200 > 10000000):
        write_bad_faces(face_iter, fp)

    return 1


    cdef int test = check_bad_face(
            face_iter[0].PolyIneq, n_coatoms, dimension + 2,
            face_iter[0].current_LHS[0],
            coatom_repr, output, e)
    return test
    if unlikely(test == 1):
        with gil:
            print("Discovered non-empty face")
    elif unlikely(test == 2):
        with gil:
            raise ValueError("Wilf conjection does not hold, counterexample: {}".format(tuple(coatom_repr[i] for i in range(output))))
    return 1

cdef inline int write_bad_faces(iter_struct *face_iter, FILE *fp) nogil except -1:
    # Obtain values from the iterator structure.
    cdef uint8_t *bad_faces     = face_iter[0].bad_faces
    cdef size_t len_bad_faces   = face_iter[0].len_bad_faces
    cdef uint64_t *bad_faces_LHS= face_iter[0].bad_faces_LHS
    cdef size_t n_bad_faces     = face_iter[0].n_bad_faces
    if len_bad_faces == 0:
        return 0
    fwrite(&len_bad_faces, 1, sizeof(size_t), fp)
    fwrite(&n_bad_faces, 1, sizeof(size_t), fp)
    fwrite(bad_faces, len_bad_faces, sizeof(uint8_t), fp)
    fwrite(bad_faces_LHS, n_bad_faces, sizeof(uint64_t), fp)
    face_iter[0].n_bad_faces = 0
    face_iter[0].len_bad_faces = 0

cdef class FaceIterator(SageObject):
    r"""
    A class to iterate over all faces of a polyhedron.

    Constructs all proper faces from the facets. In dual mode, constructs all proper
    faces from the vertices. Dual will be faster for less vertices than facets.

    INPUT:

    - ``C`` -- a :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`
    - ``dual`` -- if ``True``, then dual polyhedron is used for iteration
      (only possible for bounded Polyhedra)
    - ``output_dimension`` -- if not ``None``, then the FaceIterator will only yield
      faces of this dimension

    .. SEEALSO::

        :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`.

    EXAMPLES:

    Construct a FaceIterator::

        sage: P = polytopes.cuboctahedron()
        sage: C = CombinatorialPolyhedron(P)
        sage: it = C.face_iter()
        sage: next(it)
        A 0-dimensional face of a 3-dimensional combinatorial polyhedron

    Construct faces by the dual or not::

        sage: it = C.face_iter(dual=False)
        sage: next(it).dimension()
        2

        sage: it = C.face_iter(dual=True)
        sage: next(it).dimension()
        0

    For unbounded polyhedra only non-dual iteration is possible::

        sage: P = Polyhedron(rays=[[0,0,1], [0,1,0], [1,0,0]])
        sage: C = CombinatorialPolyhedron(P)
        sage: it = C.face_iter()
        sage: [face.Vrepr() for face in it]
        [(A vertex at (0, 0, 0),
          A ray in the direction (0, 1, 0),
          A ray in the direction (1, 0, 0)),
         (A vertex at (0, 0, 0),
          A ray in the direction (0, 0, 1),
          A ray in the direction (1, 0, 0)),
         (A vertex at (0, 0, 0),
          A ray in the direction (0, 0, 1),
          A ray in the direction (0, 1, 0)),
         (A vertex at (0, 0, 0), A ray in the direction (1, 0, 0)),
         (A vertex at (0, 0, 0), A ray in the direction (0, 1, 0)),
         (A vertex at (0, 0, 0),),
         (A vertex at (0, 0, 0), A ray in the direction (0, 0, 1))]
        sage: it = C.face_iter(dual=True)
        Traceback (most recent call last):
        ...
        ValueError: cannot iterate over dual of unbounded Polyedron

    Construct a FaceIterator only yielding dimension `2` faces::

        sage: P = polytopes.permutahedron(5)
        sage: C = CombinatorialPolyhedron(P)
        sage: it = C.face_iter(dimension=2)
        sage: counter = 0
        sage: for _ in it: counter += 1
        sage: print ('permutahedron(5) has', counter,
        ....:        'faces of dimension 2')
        permutahedron(5) has 150 faces of dimension 2
        sage: C.f_vector()
        (1, 120, 240, 150, 30, 1)

    In non-dual mode one can ignore all faces contained in the current face::

        sage: P = polytopes.cube()
        sage: C = CombinatorialPolyhedron(P)
        sage: it = C.face_iter(dual=False)
        sage: face = next(it)
        sage: face.Hrepr(names=False)
        (5,)
        sage: it.ignore_subfaces()
        sage: [face.Hrepr(names=False) for face in it]
        [(4,),
         (3,),
         (2,),
         (1,),
         (0,),
         (3, 4),
         (2, 4),
         (1, 4),
         (1, 3, 4),
         (1, 2, 4),
         (1, 3),
         (0, 3),
         (0, 1, 3),
         (1, 2),
         (0, 2),
         (0, 1, 2),
         (0, 1)]

        sage: it = C.face_iter(dual=True)
        sage: next(it)
        A 0-dimensional face of a 3-dimensional combinatorial polyhedron
        sage: it.ignore_subfaces()
        Traceback (most recent call last):
        ...
        ValueError: only possible when not in dual mode

    In dual mode one can ignore all faces that contain the current face::

        sage: it = C.face_iter(dual=True)
        sage: next(it)
        A 0-dimensional face of a 3-dimensional combinatorial polyhedron
        sage: face = next(it)
        sage: face.Vrepr(names=False)
        (6,)
        sage: [face.Vrepr(names=False) for face in it]
        [(5,),
         (4,),
         (3,),
         (2,),
         (1,),
         (0,),
         (6, 7),
         (5, 7),
         (3, 7),
         (4, 5, 6, 7),
         (2, 3, 6, 7),
         (1, 3, 5, 7),
         (4, 6),
         (2, 6),
         (0, 2, 4, 6),
         (4, 5),
         (1, 5),
         (0, 1, 4, 5),
         (0, 4),
         (2, 3),
         (1, 3),
         (0, 1, 2, 3),
         (0, 2),
         (0, 1)]

        sage: it = C.face_iter(dual=False)
        sage: next(it)
        A 2-dimensional face of a 3-dimensional combinatorial polyhedron
        sage: it.ignore_supfaces()
        Traceback (most recent call last):
        ...
        ValueError: only possible when in dual mode

    ALGORITHM:

    A (slightly generalized) description of the algorithm can be found in [KS2019]_.

    The algorithm to visit all proper faces exactly once is roughly
    equivalent to::

        faces = [set(facet) for facet in P.facets()]
        face_iterator(faces, [])

        def face_iterator(faces, visited_all):
            # Visit all faces of a polyhedron `P`, except those contained in
            # any of the visited all.

            # Assumes ``faces`` to be exactly those facets of `P`
            # that are not contained in any of the ``visited_all``.

            # Assumes ``visited_all`` to be some list of faces of
            # a polyhedron `P_2`, which contains `P` as one of its faces.

            while len(facets) > 0:
                one_face = faces.pop()
                maybe_newfaces = [one_face.intersection(face) for face in faces]

                # ``maybe_newfaces`` contains all facets of ``one_face``,
                # which we have not visited before.
                # Proof: Let `F` be a facet of ``one_face``.
                # We have a chain:
                # `P` ⊃ ``one_face`` ⊃ `F`.
                # By diamond property there exists ``second_face`` with:
                # `P` ⊃ ``second_face`` ⊃ `F`.

                # Either ``second_face`` is not an element of ``faces``:
                #     Hence ``second_face`` is contained in one of ``visited_all``.
                #     In particular, `F` is contained in ``visited_all``.
                # Or ``second_face`` is an element of ``faces``:
                #     Then, intersecting ``one_face`` with ``second_face`` gives
                #     ``F``. ∎

                # If an element in ``maybe_newfaces`` is inclusion maximal
                # and not contained in any of the ``visited_all``,
                # it is a facet of ``one_face``.
                # Any facet in ``maybe_newfaces`` of ``one_face``
                # is inlcusion maximal.
                maybe_newfaces2 = []
                for face1 in maybe_newfaces:
                    # ``face1`` is a facet of ``one_face``,
                    # iff it is not contained in another facet.
                    if all(not face1 < face2 for face2 in maybe_newfaces):
                        maybe_newfaces2.append(face1)

                # ``maybe_newfaces2`` contains only facets of ``one_face``
                # and some faces contained in any of ``visited_all``.
                # It also contains all the facets not contained in any of ``visited_all``.
                # Let ``newfaces`` be the list of all facets of ``one_face``
                # not contained in any of ``visited_all``.
                newfaces = []
                for face1 in maybe_newfaces2:
                    if all(not face1 < face2 for face2 in visited_all):
                        newfaces.append(face1)

                # By induction we can apply the algorithm, to visit all
                # faces of ``one_face`` not contained in ``visited_all``:
                face_iterator(newfaces, visited_all)

                # Finally visit ``one_face`` and add it to ``visited_all``:
                visit(one_face)
                visited_all.append(one_face)

                # Note: At this point, we have visited exactly those faces,
                # contained in any of the ``visited_all``.
    """
    def __init__(self, E, bint dual, output_dimension=None):
        r"""
        Initialize :class:`FaceIterator`.

        See :class:`FaceIterator`.

        EXAMPLES::

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: f_vector = [1, 0, 0, 0, 1]
            sage: for face in it: f_vector[face.dimension()+1] += 1
            sage: print ('f_vector of permutahedron(4): ', f_vector)
            f_vector of permutahedron(4):  [1, 24, 36, 14, 1]

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator).run()
        """
        assert(isinstance(E, CombinatorialPolyhedron), "wrong type")
        cdef KunzCone D
        cdef CombinatorialPolyhedron C = E

        if dual and not C.is_bounded():
            raise ValueError("cannot iterate over dual of unbounded Polyedron")
        cdef int i
        cdef ListOfFaces some_list  # make Cython aware of type

        self.dual = dual
        self.structure.dual = dual
        self.structure.face = NULL
        self.structure.dimension = C.dimension()
        self.structure.current_dimension = self.structure.dimension -1
        self._mem = MemoryAllocator()

        # We will not yield the empty face.
        # If there are `n` lines, than there
        # are no faces below dimension `n`.
        # The dimension of the level-sets in the face lattice jumps from `n` to `-1`.
        self.structure.lowest_dimension = 0

        if output_dimension is not None:
            if not output_dimension in range(0,self.structure.dimension):
                raise ValueError("``output_dimension`` must be the dimension of proper faces")
            if self.dual:
                # In dual mode, the dimensions are reversed.
                self.structure.output_dimension = self.structure.dimension - 1 - output_dimension
            else:
                self.structure.output_dimension = output_dimension
            self.structure.lowest_dimension = max(0, self.structure.output_dimension)
        else:
            self.structure.output_dimension = -2

        if dual:
            self.atoms = C.bitrep_facets().__copy__()
            self.coatoms = C.bitrep_Vrepr().__copy__()
        else:
            self.coatoms = C.bitrep_facets().__copy__()
            self.atoms = C.bitrep_Vrepr().__copy__()
        self.structure.face_length = self.coatoms.face_length
        self._V = C.V()
        self._H = C.H()
        self._equalities = C.equalities()
        self.structure.n_coatoms = self.coatoms.n_faces

        self.structure.atom_repr = <size_t *> self._mem.allocarray(self.coatoms.n_atoms, sizeof(size_t))
        self.structure.coatom_repr = <size_t *> self._mem.allocarray(self.coatoms.n_faces, sizeof(size_t))

        if self.structure.dimension == 0 or self.coatoms.n_faces == 0:
            # As we will only yield proper faces,
            # there is nothing to yield in those cases.
            # We have to discontinue initialization,
            # as it assumes ``self.dimension > 0`` and ``self.n_faces > 0``.
            self.structure.current_dimension = self.structure.dimension
            return
        # We may assume ``dimension > 0`` and ``n_faces > 0``.

        # Initialize ``maybe_newfaces``,
        # the place where the new faces are being stored.
        self.newfaces_lists = tuple(ListOfFaces(self.coatoms.n_faces, self.coatoms.n_atoms)
                                    for i in range(self.structure.dimension -1))
        self.structure.maybe_newfaces = <uint64_t ***> self._mem.allocarray((self.structure.dimension -1), sizeof(uint64_t **))
        for i in range(self.structure.dimension -1):
            some_list = self.newfaces_lists[i]
            self.structure.maybe_newfaces[i] = some_list.data

        # Initialize ``visited_all``.
        self.structure.visited_all = <uint64_t **> self._mem.allocarray(self.coatoms.n_faces, sizeof(uint64_t *))
        self.structure.n_visited_all = <size_t *> self._mem.allocarray(self.structure.dimension, sizeof(size_t))
        self.structure.n_visited_all[self.structure.dimension -1] = 0
        self.structure.bounded = C.is_bounded()
        if not C.is_bounded():
            # Treating the far face as if we had visited all its elements.
            # Hence we will visit all intersections of facets unless contained in the far face.

            # Regarding the length of ``self.visited_all``:
            # The last facet will not yield any new faces thus the length of ``visited_all``
            # needs to be at most ``n_facets - 1``.
            # Hence it is fine to use the first entry already for the far face,
            # as ``self.visited_all`` holds ``n_facets`` pointers.
            some_list = C.far_face()
            self.structure.visited_all[0] = some_list.data[0]
            self.structure.n_visited_all[self.structure.dimension -1] = 1

        # Initialize ``newfaces``, which will point to the new faces of codimension 1,
        # which have not been visited yet.
        self.structure.newfaces = <uint64_t ***> self._mem.allocarray(self.structure.dimension, sizeof(uint64_t **))
        for i in range(self.structure.dimension - 1):
            self.structure.newfaces[i] = <uint64_t **> self._mem.allocarray(self.coatoms.n_faces, sizeof(uint64_t *))
        self.structure.newfaces[self.structure.dimension - 1] = self.coatoms.data  # we start with coatoms

        # Initialize ``n_newfaces``.
        self.structure.n_newfaces = <size_t *> self._mem.allocarray(self.structure.dimension, sizeof(size_t))
        self.structure.n_newfaces[self.structure.dimension - 1] = self.coatoms.n_faces

        cdef size_t n,j
        # Initialize pointers to non-zero entries of newfaces etc.
        self.structure.nonzero_newfaces = <uint32_t ***> self._mem.allocarray(self.structure.dimension, sizeof(uint32_t **))
        self.structure.nonzero_maybe_newfaces = <uint32_t ***> self._mem.allocarray(self.structure.dimension, sizeof(uint32_t **))
        for i in range(self.structure.dimension - 1):
            self.structure.nonzero_newfaces[i] = <uint32_t **> self._mem.allocarray(self.coatoms.n_faces, sizeof(uint32_t *))
            self.structure.nonzero_maybe_newfaces[i] = <uint32_t **> self._mem.allocarray(self.coatoms.n_faces, sizeof(uint32_t *))
            for j in range(self.coatoms.n_faces):
                self.structure.nonzero_maybe_newfaces[i][j] = <uint32_t *> self._mem.allocarray(self.structure.face_length*64/chunksize+1, sizeof(uint32_t))

        # Initialize ``first_time``.
        self.structure.first_time = <bint *> self._mem.allocarray(self.structure.dimension, sizeof(bint))
        self.structure.first_time[self.structure.dimension - 1] = True

        self.structure.yet_to_visit = self.coatoms.n_faces
        self.structure._index = 0

        self.structure.max_dimension = self.structure.dimension
        self.structure.current_stadium = <size_t*> self._mem.calloc(self.structure.dimension, sizeof(size_t))
        self.structure.is_not_newface = <int*> self._mem.allocarray(self.coatoms.n_faces, sizeof(int))

        if isinstance(E, KunzCone):
            D = E
            n = len(E._orbit_first_element)
            self.structure.n_first_orbit_facets = n
            self.structure.first_orbit_facets = <size_t *> self._mem.calloc(n, sizeof(size_t))
            for j in range(n):
                self.structure.first_orbit_facets[j] = D._orbit_first_element[j]
            self.structure.LHS = <uint64_t **> self._mem.calloc(self.structure.dimension, sizeof(uint64_t*))
            self.structure.RHS = <uint64_t **> self._mem.calloc(self.structure.dimension, sizeof(uint64_t*))
            for i in range(self.structure.dimension):
                self.structure.LHS[i] = <uint64_t *> self._mem.calloc(self.structure.n_coatoms, sizeof(uint64_t))
                self.structure.RHS[i] = <uint64_t *> self._mem.calloc(self.structure.n_coatoms, sizeof(uint64_t))

            for j in range(self.structure.n_coatoms):
                self.structure.LHS[self.structure.dimension-1][j] = D.facets_LHS[j]
                self.structure.RHS[self.structure.dimension-1][j] = D.facets_RHS[j]
            self.structure.PolyIneq = <size_t **> self._mem.calloc(self.structure.n_coatoms, sizeof(size_t*))
            for j in range(self.structure.n_coatoms):
                self.structure.PolyIneq[j] = <size_t *> self._mem.calloc(3, sizeof(size_t))
                self.structure.PolyIneq[j][0] = D.PolyIneq[j][0]
                self.structure.PolyIneq[j][1] = D.PolyIneq[j][1]
                self.structure.PolyIneq[j][2] = D.PolyIneq[j][2]

            # Space for the bad faces.
            self.structure.bad_faces = <uint8_t *> self._mem.calloc(10000000, sizeof(uint8_t))
            self.structure.bad_faces_pt = <uint8_t **> self._mem.calloc(1000000, sizeof(uint8_t*))
            self.structure.len_bad_faces = 0
            self.structure.bad_faces_LHS = <uint64_t *> self._mem.calloc(1000000, sizeof(uint64_t))
            self.structure.n_bad_faces = 0
        else:
            self.structure.LHS = NULL
            self.structure.RHS = NULL

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: P = polytopes.associahedron(['A',3])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_iter()
            Iterator over the proper faces of a 3-dimensional combinatorial polyhedron

            sage: C.face_iter(1)
            Iterator over the 1-faces of a 3-dimensional combinatorial polyhedron
        """
        if self.structure.output_dimension != -2:
            if self.dual:
                # ouput_dimension is stored with respect to the dual
                intended_dimension = self.structure.dimension - 1 - self.structure.output_dimension
            else:
                intended_dimension = self.structure.output_dimension
            output = "Iterator over the {}-faces".format(intended_dimension)
        else:
            output = "Iterator over the proper faces"
        return output + " of a {}-dimensional combinatorial polyhedron".format(self.structure.dimension)

    def __next__(self):
        r"""
        Return the next face.

        EXAMPLES::
            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: [next(it) for _ in range(7)]
            [A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron]
        """
        cdef CombinatorialFace face = self.next_face()
        if unlikely(self.structure.current_dimension == self.structure.dimension):
            raise StopIteration

        return face

    next = __next__

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: [d for d in it]
            [A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron]
        """
        return self

    def __reduce__(self):
        r"""
        Override __reduce__ to indicate that pickle/unpickle will not work.

        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: it1 = loads(it.dumps())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def ignore_subfaces(self):
        r"""
        :class:`FaceIterator` will not visit any faces of the current face.

        Only possible when not in dual mode.

        EXAMPLES::

            sage: P = polytopes.Gosset_3_21()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(dual=False)
            sage: n_non_simplex_faces = 1
            sage: for face in it:
            ....:     if face.length_Vrepr() > face.dimension() + 1:
            ....:         n_non_simplex_faces += 1
            ....:     else:
            ....:         it.ignore_subfaces()
            ....:
            sage: n_non_simplex_faces
            127
        """
        if unlikely(self.dual):
            raise ValueError("only possible when not in dual mode")
        if unlikely(self.structure.face is NULL):
            raise ValueError("iterator not set to a face yet")

        # The current face is added to ``visited_all``.
        # This will make the iterator skip those faces.
        # Also, this face will not be added a second time to ``visited_all``,
        # as there are no new faces.
        self.structure.visited_all[self.structure.n_visited_all[self.structure.current_dimension]] = self.structure.face
        self.structure.n_visited_all[self.structure.current_dimension] += 1

    def ignore_supfaces(self):
        r"""
        :class:`FaceIterator` will not visit any faces of the current face.

        Only possible when not in dual mode.

        EXAMPLES::

            sage: P = polytopes.Gosset_3_21()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(dual=True)
            sage: n_faces_with_non_simplex_quotient = 1
            sage: for face in it:
            ....:     if face.length_Hrepr() > C.dimension() - face.dimension() + 1:
            ....:         n_faces_with_non_simplex_quotient += 1
            ....:     else:
            ....:         it.ignore_supfaces()
            ....:
            sage: n_faces_with_non_simplex_quotient
            4845
        """
        if unlikely(not self.dual):
            raise ValueError("only possible when in dual mode")
        if unlikely(self.structure.face is NULL):
            raise ValueError("iterator not set to a face yet")

        # The current face is added to ``visited_all``.
        # This will make the iterator skip those faces.
        # Also, this face will not be added a second time to ``visited_all``,
        # as there are no new faces.
        self.structure.visited_all[self.structure.n_visited_all[self.structure.current_dimension]] = self.structure.face
        self.structure.n_visited_all[self.structure.current_dimension] += 1

    cdef inline CombinatorialFace next_face(self):
        r"""
        Set attribute ``face`` to the next face and return it as
        :class:`sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace`.
        """
        self.next_dimension()
        if unlikely(self.structure.current_dimension == self.structure.dimension):
            return None
        return CombinatorialFace(self)

    cdef inline int next_dimension(self) except -1:
        r"""
        Set attribute ``face`` to the next face and return the dimension.

        Will return the dimension of the polyhedron on failure.

        The function calls :meth:`FaceIterator.next_face_loop` until a new
        face is set or until the iterator is consumed.

        .. NOTE::

            The face_iterator can be prevented from visiting any subfaces
            (or supfaces in dual mode) as in :meth:`FaceIterator.ignore_subfaces`
            and :meth`FaceIterator.ignore_supfaces`.

            Those methods add the current face to ``visited_all`` before
            visiting sub-/supfaces instead of after. One cannot arbitrarily
            add faces to ``visited_all``, as visited_all has a maximal length.
        """
        cdef int dim = self.structure.dimension
        while (not self.next_face_loop()) and (self.structure.current_dimension < dim):
            pass
        self.structure._index += 1
        return self.structure.current_dimension

    cdef inline int next_face_loop(self) except -1:
        r"""
        Set attribute ``face`` to the next face. On success return `1`.
        Otherwise `0`. Needs to be recalled then.

        If ``self.current_dimension == self.dimension``, then the iterator is
        consumed.
        """
        if unlikely(self.structure.current_dimension == self.structure.dimension):
            # The function is not supposed to be called,
            # just prevent it from crashing.
            raise StopIteration
        return next_face_loop(&self.structure)

    cdef size_t length_atom_repr(self) except -1:
        r"""
        Compute the number of atoms in the current face by counting the
        number of set bits.

        This is a shortcut of :class:`sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace.length_atom_repr`
        """
        if self.structure.face:
            return count_atoms(self.structure.face, self.structure.face_length)

        # The face was not initialize properly.
        raise LookupError("``FaceIterator`` does not point to a face")

    cdef size_t set_coatom_repr(self) except -1:
        r"""
        Set ``coatom_repr`` to be the coatom-representation of the current face.
        Return its length.

        This is a shortcut of :class:`sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace.set_coatom_repr`
        """
        cdef size_t n_coatoms = self.coatoms.n_faces
        cdef uint64_t **coatoms = self.coatoms.data
        cdef size_t face_length = self.structure.face_length
        return bit_repr_to_coatom_repr(self.structure.face, coatoms, n_coatoms,
                                       face_length, self.structure.coatom_repr)

    cdef size_t set_atom_repr(self) except -1:
        r"""
        Set ``atom_repr`` to be the atom-representation of the current face.
        Return its length.

        This is a shortcut of :class:`sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace.set_atom_repr`
        """
        cdef size_t face_length = self.structure.face_length
        return bit_repr_to_Vrepr_list(self.structure.face, self.structure.atom_repr, face_length)
