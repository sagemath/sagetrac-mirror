# distutils: language = c++
# distutils: extra_compile_args = "-march=native"

cimport cython
from libc.stdint                cimport uint64_t
from libcpp.pair cimport pair as cpair
from libcpp.utility             cimport pair

cdef extern from "bit_vector_operations.cc":
    # Any Bit-representation is assumed to be `chunksize`-Bit aligned.
    cdef const size_t chunksize;
    cdef void intersection(uint64_t *A, uint64_t *B, uint64_t *C, size_t face_length)
#    Return ``A & ~B == 0``.
#    A is not subset of B, iff there is a vertex in A, which is not in B.
#    ``face_length`` is the length of A and B in terms of uint64_t.


    ctypedef cpair[int,size_t] mypair

    cdef size_t get_next_level(
        uint64_t **faces, const size_t n_faces, uint64_t **nextfaces,
        uint64_t **nextfaces2, uint64_t **visited_all,
        size_t n_visited_all, size_t face_length, int *is_not_newface, mypair *sorting_array) nogil
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

    cdef size_t bit_repr_to_coatom_repr( \
            uint64_t *face, uint64_t **coatoms, size_t n_coatoms, \
            size_t face_length, size_t *output) nogil
#        Write the coatom-representation of face in output. Return length.
#        ``face_length`` is the length of ``face`` and ``coatoms[i]``
#        in terms of uint64_t.
#        ``n_coatoms`` length of ``coatoms``.
    cdef void parallel_f_vector(iter_struct **face_iter, size_t *f_vector, size_t n_threads, size_t recursion_depth) nogil
    cdef struct iter_struct:
        bint dual                  # if 1, then iterate over dual Polyhedron
        uint64_t *face             # the current face of the iterator
        size_t *atom_repr          # a place where atom-representaion of face will be stored
        size_t *coatom_repr        # a place where coatom-representaion of face will be stored
        int current_dimension      # dimension of current face, dual dimension if ``dual``
        int dimension              # dimension of the polyhedron
        int n_lines                # ``_n_lines`` of ``CombinatorialPolyhedron``
        int output_dimension       # only faces of this (dual?) dimension are considered
        int lowest_dimension       # don't consider faces below this (dual?) dimension
        int max_dimension
        size_t face_length         # stores length of the faces in terms of uint64_t
        size_t _index, n_coatoms

        # ``visited_all`` points to faces, of which we have visited all faces already.
        # The number of faces in ``visited_all` might depend on the current dimension:
        #     Consider we visit the facets A,B of some face F.
        #     We will first visit all faces of A and then add A to visited_all.
        #     Then we visit all faces of B and add B to visited_all.
        #     Then we have visited F completely.
        #     Instead of having A and B in ``visited_all`` we will point to F.
        #     In this way, we will append ``visited_all`` in lower dimension, but
        #     will ignore those changes when going up in dimension again.
        #     This is why the number of faces in ``visited_all``depends on dimension.
        uint64_t **visited_all
        size_t *n_visited_all

        # ``maybe_newfaces`` is where all possible facets of a face are stored.
        # In dimension ``dim`` when visiting all faces of some face,
        # the intersections with other faces are stored in ``newfaces2[dim]``.
        uint64_t ***maybe_newfaces

        # ``newfaces`` will point to those faces in ``maybe_newfaces``
        # that are of codimension 1 and not already visited.
        uint64_t ***newfaces
        size_t *n_newfaces  # number of newfaces for each dimension

        # After having visited a face completely, we want to add it to ``visited_all``.
        # ``first_dim[i]`` will indicate, wether there is one more face in
        # ``newfaces[i]`` then ``n_newfaces[i]`` suggests
        # that has to be added to ``visited_all``.
        # If ``first_time[i] == False``, we still need to
        # add ``newfaces[i][n_newfaces[i]]`` to ``visited_all``.
        bint *first_time

        # The number of elements in newfaces[current_dimension],
        # that have not been visited yet.
        size_t yet_to_visit
        int *is_not_newface
        mypair *sorting_array
        size_t *current_stadium
