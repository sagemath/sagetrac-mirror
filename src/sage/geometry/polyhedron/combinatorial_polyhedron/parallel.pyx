# distutils: language = c++
# distutils: extra_compile_args = "-march=native"
from libc.stdlib cimport calloc, free, malloc
from cython.parallel cimport threadid, prange, parallel
from libc.stdint                cimport uint64_t
from .bit_vector_operations cimport get_next_level
from cysignals.signals      cimport sig_check, sig_on, sig_off

cimport openmp

cdef size_t rec_depth = 2

#cdef extern from "<omp.h>":
#    extern int omp_get_thread_num() nogil

cdef int parallel_f_vector(iter_struct **face_iter, size_t *f_vector) nogil except -1:
    cdef size_t **shared_f = <size_t **> calloc(8, sizeof(size_t*))
    cdef int dimension = face_iter[0][0].dimension
    cdef size_t i
    cdef int j

    for i in range(8):
        shared_f[i] = <size_t *> calloc(dimension + 2, sizeof(size_t))
    cdef size_t n_faces = face_iter[0][0].n_coatoms

    cdef int * local_buf
    cdef int n = 0
    cdef Py_ssize_t l = 0
    cdef iter_struct *my_iter
    cdef size_t *my_f
    cdef int chunksize = 1
    #local_buf[0] = 0#openmp.omp_get_thread_num()
    #my_iter = <iter_struct*> malloc(sizeof(iter_struct) * 1)
    #my_f = <size_t *> malloc(sizeof(size_t*) * 1)
    #with gil:
    #    my_iter = face_iter[local_buf[0]]
    #    my_f = shared_f[local_buf[0]]
    with nogil:
        for l in prange(0, n_faces ** rec_depth, num_threads=8, schedule='dynamic', chunksize=chunksize):
            partial_f(face_iter[openmp.omp_get_thread_num()], shared_f[openmp.omp_get_thread_num()], l)
    #free(my_f)
    #free(my_iter)

    for i in range(8):
        for j in range(0,dimension + 2):
            f_vector[j] += shared_f[i][j]
        free(shared_f[i])
    free(shared_f)


cdef int partial_f(iter_struct *face_iter, size_t *f_vector, size_t i) nogil except -1:
    prepare_partial_iter(face_iter, i, f_vector)
    cdef int d, dimension = face_iter[0].dimension
    cdef size_t j
    d = next_dimension(face_iter)
    while d < dimension -rec_depth:
        f_vector[d + 1] += 1
        d = next_dimension(face_iter)


cdef inline int next_face_loop(iter_struct *face_iter) nogil except -1:
    r"""
    Set attribute ``face`` to the next face. On success return `1`.
    Otherwise `0`. Needs to be recalled then.

    If ``face_iter[0].current_dimension == face_iter[0].dimension``, then the iterator is
    consumed.
    """

    # Getting ``[faces, n_faces, n_visited_all]`` according to dimension.
    cdef uint64_t **faces = face_iter[0].newfaces[face_iter[0].current_dimension]
    cdef size_t n_faces = face_iter[0].n_newfaces[face_iter[0].current_dimension]
    cdef size_t n_visited_all = face_iter[0].n_visited_all[face_iter[0].current_dimension]

    if (face_iter[0].output_dimension > -2) and (face_iter[0].output_dimension != face_iter[0].current_dimension):
        # If only a specific dimension was requested (i.e. ``face_iter[0].output_dimension > -2``),
        # then we will not yield faces in other dimension.
        face_iter[0].yet_to_visit = 0

    if face_iter[0].yet_to_visit:
        # Set ``face`` to the next face.
        face_iter[0].yet_to_visit -= 1
        face_iter[0].face = faces[face_iter[0].yet_to_visit]
        return 1

    if face_iter[0].current_dimension <= face_iter[0].lowest_dimension:
        # We will not yield the empty face.
        # We will not yield below requested dimension.
        face_iter[0].current_dimension += 1
        return 0

    if n_faces <= 1:
        # There will be no more faces from intersections.
        face_iter[0].current_dimension += 1
        return 0

    # We will visit the last face now.
    face_iter[0].n_newfaces[face_iter[0].current_dimension] -= 1
    n_faces -= 1

    if not face_iter[0].first_time[face_iter[0].current_dimension]:
        # In this case there exists ``faces[n_faces + 1]``, of which we
        # have visited all faces, but which was not added to
        # ``visited_all`` yet.
        face_iter[0].visited_all[n_visited_all] = faces[n_faces + 1]
        face_iter[0].n_visited_all[face_iter[0].current_dimension] += 1
        n_visited_all = face_iter[0].n_visited_all[face_iter[0].current_dimension]
    else:
        # Once we have visited all faces of ``faces[n_faces]``, we want
        # to add it to ``visited_all``.
        face_iter[0].first_time[face_iter[0].current_dimension] = False

    # Get the faces of codimension 1 contained in ``faces[n_faces]``,
    # which we have not yet visited.
    cdef size_t newfacescounter

    newfacescounter = get_next_level(
        faces, n_faces + 1, face_iter[0].maybe_newfaces[face_iter[0].current_dimension-1],
        face_iter[0].newfaces[face_iter[0].current_dimension-1],
        face_iter[0].visited_all, n_visited_all, face_iter[0].face_length, face_iter[0].is_not_newface)

    if newfacescounter:
        # ``faces[n_faces]`` contains new faces.
        # We will visted them on next call, starting with codimension 1.

        # Setting the variables correclty for next call of ``next_face_loop``.
        face_iter[0].current_dimension -= 1
        face_iter[0].first_time[face_iter[0].current_dimension] = True
        face_iter[0].n_newfaces[face_iter[0].current_dimension] = newfacescounter
        face_iter[0].n_visited_all[face_iter[0].current_dimension] = n_visited_all
        face_iter[0].yet_to_visit = newfacescounter
        return 0
    else:
        # ``faces[n_faces]`` contains no new faces.
        # Hence there is no need to add it to ``visited_all``.
        # NOTE:
        #     For the methods ``ignore_subfaces`` and ``ignore_supfaces``
        #     this step needs to be done, as ``faces[n_faces]`` might
        #     have been added manually to ``visited_all``.
        #     So this step is required to respect boundaries of ``visited_all``.
        face_iter[0].first_time[face_iter[0].current_dimension] = True
        return 0


cdef inline int next_dimension(iter_struct *face_iter) nogil except -1:
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
        visiting sub-/supfaces instead of after. One cannot arbitralily
        add faces to ``visited_all``, as visited_all has a maximal length.
    """
    cdef int dim = face_iter[0].max_dimension
    while (not next_face_loop(&face_iter[0])) and (face_iter[0].current_dimension < dim):
        #sig_check()  # too slow, comment out for benchmarking
        pass
    face_iter[0]._index += 1
    return face_iter[0].current_dimension

cdef inline int prepare_partial_iter(iter_struct *face_iter, size_t i, size_t *f_vector) nogil except -1:
    r"""
    Prepares the face iterator to not visit the fatets 0,...,-1
    """
    cdef size_t j, k
    cdef int dimension = face_iter[0].dimension
    face_iter[0].face = NULL
    face_iter[0].current_dimension = dimension -1
    face_iter[0].lowest_dimension = face_iter[0].n_lines


    face_iter[0].n_visited_all[dimension -1] = 0
    face_iter[0].n_newfaces[dimension - 1] = face_iter[0].n_coatoms
    cdef size_t rec
    cdef size_t current_i
    cdef int d
    for rec in range(rec_depth):
        #with gil:
        #    print i
        current_i = i//(face_iter[0].n_coatoms ** (rec_depth - rec - 1))
        i = i%(face_iter[0].n_coatoms ** (rec_depth - rec - 1))
        if current_i >= face_iter[0].n_newfaces[dimension-rec-1]:
            face_iter[0].current_dimension = face_iter[0].dimension -1
            face_iter[0].n_newfaces[dimension - rec - 1] = 0
            return 0
        if i == 0:
            f_vector[dimension - rec] += 1

        k = face_iter[0].n_visited_all[dimension - rec-1]
        for j in range(face_iter[0].n_newfaces[dimension-rec-1]-current_i,face_iter[0].n_newfaces[dimension-rec-1]):
           face_iter[0].visited_all[k] = face_iter[0].newfaces[dimension -rec-1][j]
           k += 1
        face_iter[0].n_visited_all[dimension - rec-1] += current_i

        face_iter[0].n_newfaces[dimension - rec-1] -= current_i
        face_iter[0].first_time[dimension - rec-1] = True
        face_iter[0].yet_to_visit = 0
        face_iter[0].max_dimension = dimension - rec -1
        if rec < rec_depth -1:
            d = next_dimension(face_iter)
            face_iter[0].yet_to_visit = 0
