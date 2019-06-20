/*
#*****************************************************************************
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
*/
#include <cstdint>
#include <cstdio>
#include <omp.h>
#include <algorithm>    // std::sort
using namespace std;

#if __AVX__
    #include <immintrin.h>
#elif __SSE4_1__
    #include <emmintrin.h>
    #include <smmintrin.h>
#endif

#if __POPCNT__
    #include <immintrin.h>
#endif

size_t rec_depth = 2;

// As of now, 512bit does not have something like _mm256_testc_si256,
// which is the bottle neck of the algorithm,
// so it does not make sense to implement it.


// Determine the best way to check for inclusion/subset.
#if __AVX__
    // 256-bit commands, those operations are equivalent to the operations
    // defined in `#else`
    // (intrinsics defined in immintrin.h)
    const size_t chunksize = 256;
    inline bool is_subset(const uint64_t *A, const uint64_t *B, const size_t face_length){
        // A & ~B == 0
        // Return 1 if A is a subset of B, otherwise returns 0.

        // This done by checking if there is an element in A,
        // which is not in B
        // `face_length` is the length of A and B in terms of uint64_t
        // Note that A,B need to be 32-Byte-aligned
        for (size_t i = 0; i < face_length; i += 4){
            __m256i a = _mm256_load_si256((const __m256i*)&A[i]);
            __m256i b = _mm256_load_si256((const __m256i*)&B[i]);
            if (!_mm256_testc_si256(b, a)){ //need to be opposite order !!
                return false;
            }
        }
        return true;
    }

#elif __SSE4_1__
    // 128-bit commands, those operations are equivalent to the operations
    // defined in `#else`
    // (intrinsics defined in smmintrin.h and emmintrin.h)
    const size_t chunksize = 128;
    inline bool is_subset(uint64_t *A, uint64_t *B, size_t face_length){
        // A & ~B == 0
        // Return 1 if A is a subset of B, otherwise returns 0.

        // This done by checking if there is an element in A,
        // which is not in B
        // `face_length` is the length of A and B in terms of uint64_t
        // Note that A,B need to be 16-Byte-aligned
        size_t i;
        for (i = 0; i < face_length; i += 2){
            __m128i a = _mm_load_si128((const __m128i*)&A[i]);
            __m128i b = _mm_load_si128((const __m128i*)&B[i]);
            if (!_mm_testc_si128(b, a)){ //need to be opposite order !!
                return false;
            }
        }
        return true;
    }

#else
    // no intrinsics
    const size_t chunksize = 64;
    inline bool is_subset(uint64_t *A, uint64_t *B, size_t face_length){
        // A & ~B == 0
        // Return 1 if A is a subset of B, otherwise returns 0.

        // This done by checking if there is an element in A,
        // which is not in B
        // `face_length` is the length of A and B in terms of uint64_t
        size_t i;
        for (i = 0; i < face_length; i++){
            if (A[i] & ~B[i]){
                return false;
            }
        }
        return true;
    }

#endif

// Determine the best way to calculate intersection.
#if __AVX2__
    // 256-bit commands, those operations are equivalent to the operations
    // defined in `#else`
    // intrinsics defined in immintrin.h
    inline void intersection(const uint64_t *A, const uint64_t *B, uint64_t *C, \
                             const size_t face_length){
        // C = A & B
        // Set C to be the intersection of A and B.
        // ``face_length`` is the length of A, B and C in terms of ``uint64_t``.
        // Note that A,B,C need to be 32-Byte-aligned
        for (size_t i = 0; i < face_length; i += 4){
            __m256i a = _mm256_load_si256((const __m256i*)&A[i]);
            __m256i b = _mm256_load_si256((const __m256i*)&B[i]);
            __m256i c = _mm256_and_si256(a, b);
            _mm256_store_si256((__m256i*)&C[i],c);
        }
    }

#elif __SSE4_1__
    // actually SSE2 would be fine, but we don't want to force greater chunks,
    // because of intersection, which is not the bottleneck
    // 128-bit commands, those operations are equivalent to the operations
    // defined in `#else`
    // intrinsics defined in emmintrin.h
    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             size_t face_length){
        // C = A & B
        // Set C to be the intersection of A and B.
        // ``face_length`` is the length of A, B and C in terms of ``uint64_t``.
        // Note that A,B,C need to be 16-Byte-aligned
        size_t i;
        for (i = 0; i < face_length; i += 2){
            __m128i a = _mm_load_si128((const __m128i*)&A[i]);
            __m128i b = _mm_load_si128((const __m128i*)&B[i]);
            __m128i c = _mm_and_si128(a, b);
            _mm_store_si128((__m128i*)&C[i],c);
        }
    }

#else
    // commands, without intrinsics
    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             size_t face_length){
        // C = A & B
        // Set C to be the intersection of A and B.
        // ``face_length`` is the length of A, B and C in terms of ``uint64_t``.
        size_t i;
        for (i = 0; i < face_length; i++){
            C[i] = A[i] & B[i];
        }
    }

#endif

// Determine the best way to count the set bits in ``uint64_t*``.
#if (__POPCNT__) && (INTPTR_MAX == INT64_MAX) // 64-bit and popcnt
    inline size_t count_atoms(uint64_t* A, size_t face_length) {
        // Count the number of vertices in a face by counting bits set to 1.
        // ``face_length`` is the length of ``A`` in terms of ``uint64_t``.
        size_t i;
        unsigned int count = 0;
        for (i=0; i<face_length; i++){
            count += (size_t) _mm_popcnt_u64(A[i]);
        }
        return count;
    }

#else // popcount without intrinsics
    inline size_t count_atoms(uint64_t* A, size_t face_length) {
        // Count the number of vertices in a face by counting bits set to 1.
        // ``face_length`` is the length of ``A`` in terms of ``uint64_t``.
        size_t i;
        unsigned int count = 0;
        for (i=0; i<face_length; i++){
            uint64_t a = A[i];
            while (a){
                count += a & 1;
                a >>= 1;
            }
        }
        return count;
    }

#endif

typedef std::pair<int,size_t> mypair;

bool compare_pairs ( const mypair& l, const mypair& r)
       { return l.first > r.first; }

inline size_t get_next_level(\
        uint64_t **faces, const size_t n_faces, uint64_t **maybe_newfaces, \
        uint64_t **newfaces, uint64_t **visited_all, \
        size_t n_visited_all, size_t face_length, int *is_not_newface, mypair *sorting_array){
    /*
    Set ``newfaces`` to be the facets of ``faces[n_faces -1]``
    that are not contained in a face of ``visited_all``.

    INPUT:

    - ``newfaces`` -- quasi of type ``*uint64_t[n_faces -1]
    - ``visited_all`` -- quasi of type ``*uint64_t[n_visited_all]
    - ``face_length`` -- length of the faces

    OUTPUT:

    - return number of ``newfaces``
    - set ``newfaces`` to point to the new faces

    ALGORITHM:

    To get all facets of ``faces[n_faces-1]``, we would have to:
    - Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
    - Add all the intersection of ``visited_all`` with the last face
    - Out of both the inclusion-maximal ones are of codimension one, i.e. facets.

    As we have visited all faces of ``visited_all``, we alter the algorithm
    to not revisit:
    Step 1: Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
    Step 2: Out of thosse the inclusion-maximal ones are some of the facets.
            At least we obtain all of those, that we have not already visited.
            Maybe, we get some more.
    Step 3: Only keep those that we have not already visited.
            We obtain exactly the facets of ``faces[n_faces-1]`` that we have
            not visited yet.
    */

    int *is_newface = is_not_newface;
    // Step 1:
    uint64_t *face_to_intersect = faces[n_faces -1];
    for (size_t j = 0; j < n_faces - 1; j++){
        intersection(faces[j], face_to_intersect, maybe_newfaces[j], face_length);
    }

    // We keep track, which face in ``maybe_newfaces`` is a new face.
    //bool *is_not_newface = new bool[n_faces -1]();

    // For each face we will Step 2 and Step 3.
    for (size_t j = 0; j < n_faces-1; j++){
        is_newface[j] = 1;
        // Step 2a:
        //size_t k = minimal_test[j];
        for(size_t k = 0; k < j; k++){
            // Testing if maybe_newfaces[j] is contained in different nextface.
            if((is_newface[k]) && is_subset(maybe_newfaces[j], maybe_newfaces[k],face_length)){
                // If so, it is not inclusion-maximal and hence not of codimension 1.
                is_newface[j] = 0;
                break;
                }
            }
        if (!is_newface[j]) {
            // No further tests needed, if it is not of codimension 1.
            continue;
        }

        // Step 2b:
        //size_t k = minimal_test[j];
        for(size_t k = j+1; k < n_faces-1; k++){
            // Testing if maybe_newfaces[j] is contained in different nextface.
            if(is_subset(maybe_newfaces[j], maybe_newfaces[k],face_length)){
                // If so, it is not inclusion-maximal and hence not of codimension 1.
                is_newface[j] = 0;
                break;
                }
            }
        if (!is_newface[j]) {
            // No further tests needed, if it is not of codimension 1.
            continue;
        }

        // Step 3:
        for (size_t k = 0; k < n_visited_all; k++){
            // Testing if maybe_newfaces[j] is contained in one,
            // we have already completely visited.
            if(is_subset(maybe_newfaces[j], visited_all[k], face_length)){
                // If so, we don't want to revisit.
                is_newface[j] = 0;
                break;
            }
        }
    }

    // Set ``newfaces`` to point to the correct ones.
    size_t n_newfaces = 0;  // length of newfaces2
    mypair *idx = sorting_array;
    //mypair *idx = new mypair[n_faces-1];
    for (size_t j = 0; j < n_faces -1; j++){
        if (is_newface[j]) {
            is_newface[j] = count_atoms(maybe_newfaces[j], face_length) + 1;
        }
        idx[j].first =  is_newface[j];
        idx[j].second = j;
    }
    sort(idx, idx + n_faces - 1, compare_pairs);

    //for (size_t j = 0; j < n_faces-1; j++){
    //for (size_t index = n_faces - 2; index >= 0; index--){
    //for (size_t index = 0; index< n_faces-1; index++){
    for (size_t j = 0; j < n_faces-1; j++){
        newfaces[j] = maybe_newfaces[idx[j].second];
        is_newface[j] = idx[j].first;
        if (idx[j].first)
            n_newfaces++;
    }

    for (size_t j = 0; j < n_faces-1; j++){
        maybe_newfaces[j] = newfaces[j];
    }

    /*
    n_newfaces = 0;
    //for (size_t j = 0; j < n_faces-1; j++){
    for (std::vector<mypair>::iterator it=idx.begin(); it!=idx.end(); ++it){
        size_t j = (*it).second;
        if (!is_newface[j]) {
            // Not a new face of codimension 1.
            continue;
        }
        // It is a new face of codimension 1.
        newfaces[n_newfaces] = maybe_newfaces[j];
        n_newfaces++;
    }
    */


    return n_newfaces;
}


inline size_t get_next_level_simplex(\
        uint64_t **faces, const size_t n_faces, uint64_t **maybe_newfaces, \
        uint64_t **newfaces, uint64_t **visited_all, \
        size_t n_visited_all, size_t face_length, int *is_not_newface, mypair *sorting_array){
    /*
    Set ``newfaces`` to be the facets of ``faces[n_faces -1]``
    that are not contained in a face of ``visited_all``.

    INPUT:

    - ``newfaces`` -- quasi of type ``*uint64_t[n_faces -1]
    - ``visited_all`` -- quasi of type ``*uint64_t[n_visited_all]
    - ``face_length`` -- length of the faces

    OUTPUT:

    - return number of ``newfaces``
    - set ``newfaces`` to point to the new faces

    ALGORITHM:

    To get all facets of ``faces[n_faces-1]``, we would have to:
    - Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
    - Add all the intersection of ``visited_all`` with the last face
    - Out of both the inclusion-maximal ones are of codimension one, i.e. facets.

    As we have visited all faces of ``visited_all``, we alter the algorithm
    to not revisit:
    Step 1: Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
    Step 2: Out of thosse the inclusion-maximal ones are some of the facets.
            At least we obtain all of those, that we have not already visited.
            Maybe, we get some more.
    Step 3: Only keep those that we have not already visited.
            We obtain exactly the facets of ``faces[n_faces-1]`` that we have
            not visited yet.
    */

    int *is_newface = is_not_newface;
    // Step 1:
    for (size_t j = 0; j < n_faces - 1; j++){
        intersection(faces[j], faces[n_faces - 1], maybe_newfaces[j], face_length);
    }

    // We keep track, which face in ``maybe_newfaces`` is a new face.
    //bool *is_not_newface = new bool[n_faces -1]();

    // For each face we will Step 2 and Step 3.
    for (size_t j = 0; j < n_faces-1; j++){
        is_newface[j] = 1;
        // Step 2a:
        //size_t k = minimal_test[j];

        // If the common parent is simplex, it can never happen, that anyone is not inclusion maximal
        /*
        for(size_t k = 0; k < j; k++){
            // Testing if maybe_newfaces[j] is contained in different nextface.
            if((is_newface[k]) && is_subset(maybe_newfaces[j], maybe_newfaces[k],face_length)){
                // If so, it is not inclusion-maximal and hence not of codimension 1.
                is_newface[j] = 0;
                break;
                }
            }
        if (!is_newface[j]) {
            // No further tests needed, if it is not of codimension 1.
            continue;
        }

        // Step 2b:
        //size_t k = minimal_test[j];
        for(size_t k = j+1; k < n_faces-1; k++){
            // Testing if maybe_newfaces[j] is contained in different nextface.
            if(is_subset(maybe_newfaces[j], maybe_newfaces[k],face_length)){
                // If so, it is not inclusion-maximal and hence not of codimension 1.
                is_newface[j] = 0;
                break;
                }
            }
        if (!is_newface[j]) {
            // No further tests needed, if it is not of codimension 1.
            continue;
        }
        */

        // Step 3:
        for (size_t k = 0; k < n_visited_all; k++){
            // Testing if maybe_newfaces[j] is contained in one,
            // we have already completely visited.
            if(is_subset(maybe_newfaces[j], visited_all[k], face_length)){
                // If so, we don't want to revisit.
                is_newface[j] = 0;
                break;
            }
        }
    }

    /*
    // Set ``newfaces`` to point to the correct ones.
    size_t n_newfaces = 0;  // length of newfaces2
    mypair *idx = sorting_array;
    //mypair *idx = new mypair[n_faces-1];
    for (size_t j = 0; j < n_faces -1; j++){
        if (is_newface[j]) {
            is_newface[j] = count_atoms(maybe_newfaces[j], face_length) + 1;
        }
        idx[j].first =  is_newface[j];
        idx[j].second = j;
    }
    sort(idx, idx + n_faces - 1, compare_pairs);

    //for (size_t j = 0; j < n_faces-1; j++){
    //for (size_t index = n_faces - 2; index >= 0; index--){
    //for (size_t index = 0; index< n_faces-1; index++){
    for (size_t j = 0; j < n_faces-1; j++){
        newfaces[j] = maybe_newfaces[idx[j].second];
        if (idx[j].first)
            n_newfaces++;
    }

    for (size_t j = 0; j < n_faces-1; j++){
        maybe_newfaces[j] = newfaces[j];
    }

    */
    size_t n_newfaces = 0;
    for (size_t j = 0; j < n_faces-1; j++){
    //for (std::vector<mypair>::iterator it=idx.begin(); it!=idx.end(); ++it){
        //size_t j = (*it).second;
        if (!is_newface[j]) {
            // Not a new face of codimension 1.
            continue;
        }
        // It is a new face of codimension 1.
        newfaces[n_newfaces] = maybe_newfaces[j];
        n_newfaces++;
    }


    return n_newfaces;
}

size_t bit_repr_to_coatom_repr(uint64_t *face, uint64_t **coatoms, \
                               size_t n_coatoms, size_t face_length, \
                               size_t *output){
    /*
    Write the coatom-representation of face in output. Return length.
    ``face_length`` is the length of ``face`` and ``coatoms[i]``
    in terms of uint64_t.
    ``n_coatoms`` length of ``coatoms``.
    */
    size_t count_length = 0;
    for (size_t i = 0; i < n_coatoms; i++){
        if (is_subset(face, coatoms[i], face_length)){
            // ``face`` is contain in ``coatoms[i]``,
            // then ``i`` is an element in the coatom-represention.
            output[count_length] = i;
            count_length++;
        }
    }
    return count_length;
}
struct iter_struct{
    int dual;               // # if 1, then iterate over dual Polyhedron
    uint64_t *face;         //     # the current face of the iterator
    size_t *atom_repr;      //    # a place where atom-representaion of face will be stored
    size_t *coatom_repr;    //    # a place where coatom-representaion of face will be stored
    int current_dimension;  //    # dimension of current face, dual dimension if ``dual``
    int dimension;          //     # dimension of the polyhedron
    int n_lines;            //     # ``_n_lines`` of ``CombinatorialPolyhedron``
    int output_dimension;   //     # only faces of this (dual?) dimension are considered
    int lowest_dimension;   //     # don't consider faces below this (dual?) dimension
    int max_dimension;
    size_t face_length;     //     # stores length of the faces in terms of uint64_t
    size_t _index, n_coatoms;

    uint64_t **visited_all;
    size_t *n_visited_all;

    uint64_t ***maybe_newfaces;

    uint64_t ***newfaces;
    size_t *n_newfaces; // # number of newfaces for each dimension
    int *first_time;
    size_t yet_to_visit;
    int *is_not_newface;
    mypair *sorting_array;
    size_t *current_stadium;
};

inline size_t myPow(size_t x, size_t p){
      if (p == 0) return 1;
        if (p == 1) return x;

          size_t tmp = myPow(x, p/2);
            if (p%2 == 0) return tmp * tmp;
              else return x * tmp * tmp;
}

inline int next_face_loop(iter_struct *face_iter){
    //Getting ``[faces, n_faces, n_visited_all]`` according to dimension.
    uint64_t **faces = face_iter[0].newfaces[face_iter[0].current_dimension];
    size_t n_faces = face_iter[0].n_newfaces[face_iter[0].current_dimension];
    size_t n_visited_all = face_iter[0].n_visited_all[face_iter[0].current_dimension];

    if((face_iter[0].output_dimension > -2) && (face_iter[0].output_dimension != face_iter[0].current_dimension)){
        //If only a specific dimension was requested (i.e. ``face_iter[0].output_dimension > -2``),
        //then we will not yield faces in other dimension.
        face_iter[0].yet_to_visit = 0;
    }

    if (face_iter[0].yet_to_visit){
        // Set ``face`` to the next face.
        face_iter[0].yet_to_visit -= 1;
        face_iter[0].face = faces[face_iter[0].yet_to_visit];
        return 1;
    }

    if (face_iter[0].current_dimension <= face_iter[0].lowest_dimension){
        // We will not yield the empty face.
        // We will not yield below requested dimension.
        face_iter[0].current_dimension += 1;
        return 0;
    }

    if (n_faces <= 1){
        // There will be no more faces from intersections.
        face_iter[0].current_dimension += 1;
        return 0;
    }

    // We will visit the last face now.
    face_iter[0].n_newfaces[face_iter[0].current_dimension] -= 1;
    n_faces -= 1;

    if (!face_iter[0].first_time[face_iter[0].current_dimension]){
        // In this case there exists ``faces[n_faces + 1]``, of which we
        // have visited all faces, but which was not added to
        // ``visited_all`` yet.
        face_iter[0].visited_all[n_visited_all] = faces[n_faces + 1];
        face_iter[0].n_visited_all[face_iter[0].current_dimension] += 1;
        n_visited_all = face_iter[0].n_visited_all[face_iter[0].current_dimension];
    } else {
        // Once we have visited all faces of ``faces[n_faces]``, we want
        // to add it to ``visited_all``.
        face_iter[0].first_time[face_iter[0].current_dimension] = 0;
    }

    // Get the faces of codimension 1 contained in ``faces[n_faces]``,
    // which we have not yet visited.
    size_t newfacescounter;

    newfacescounter = get_next_level( \
        faces, n_faces + 1, face_iter[0].maybe_newfaces[face_iter[0].current_dimension-1], \
        face_iter[0].newfaces[face_iter[0].current_dimension-1],\
        face_iter[0].visited_all, n_visited_all, face_iter[0].face_length, face_iter[0].is_not_newface, face_iter[0].sorting_array);
    /*
    if (face_iter[0].is_simplex[face_iter[0].current_dimension][n_faces] < 2){
        newfacescounter = get_next_level( \
            faces, n_faces + 1, face_iter[0].maybe_newfaces[face_iter[0].current_dimension-1], \
            face_iter[0].newfaces[face_iter[0].current_dimension-1],\
            face_iter[0].visited_all, n_visited_all, face_iter[0].face_length, face_iter[0].is_not_newface, face_iter[0].sorting_array);
    } else {

        // all the newfaces are actually children of a simplex, so no need to check for inclusion maximal
        newfacescounter = get_next_level_simplex( \
            faces, n_faces + 1, face_iter[0].maybe_newfaces[face_iter[0].current_dimension-1], \
            face_iter[0].newfaces[face_iter[0].current_dimension-1],\
            face_iter[0].visited_all, n_visited_all, face_iter[0].face_length, face_iter[0].is_not_newface, face_iter[0].sorting_array);
    }

    int previous = face_iter[0].is_simplex[face_iter[0].current_dimension][n_faces] == 2;
    if (previous == 0){
        for(size_t i=0; i< newfacescounter; i++){
            int n_atoms = face_iter[0].is_not_newface[i] - 1;
            if (n_atoms != face_iter[0].current_dimension){
                // the face is not simplex
                face_iter[0].is_simplex[face_iter[0].current_dimension -1][i] = 0;
            }
            else {
                face_iter[0].is_simplex[face_iter[0].current_dimension -1][i] = 1;
            }
        }
    } else {
        for(size_t i=0; i< newfacescounter; i++){
            face_iter[0].is_simplex[face_iter[0].current_dimension -1][i] = previous + 1;
        }
    }
    */

    if (newfacescounter){
        // ``faces[n_faces]`` contains new faces.
        // We will visted them on next call, starting with codimension 1.

        // Setting the variables correclty for next call of ``next_face_loop``.
        face_iter[0].current_dimension -= 1;
        face_iter[0].first_time[face_iter[0].current_dimension] = 1;
        face_iter[0].n_newfaces[face_iter[0].current_dimension] = newfacescounter;
        face_iter[0].n_visited_all[face_iter[0].current_dimension] = n_visited_all;
        face_iter[0].yet_to_visit = newfacescounter;
        return 0;
    } else {
        //  ``faces[n_faces]`` contains no new faces.
        //  Hence there is no need to add it to ``visited_all``.
        //  NOTE:
        //      For the methods ``ignore_subfaces`` and ``ignore_supfaces``
        //      this step needs to be done, as ``faces[n_faces]`` might
        //      have been added manually to ``visited_all``.
        //      So this step is required to respect boundaries of ``visited_all``.
        face_iter[0].first_time[face_iter[0].current_dimension] = 1;
        return 0;
    }
}

inline int next_dimension(iter_struct *face_iter){
    // Set attribute ``face`` to the next face and return the dimension.

    // Will return the dimension of the polyhedron on failure.

    // The function calls :meth:`FaceIterator.next_face_loop` until a new
    // face is set or until the iterator is consumed.

    // .. NOTE::

    //     The face_iterator can be prevented from visiting any subfaces
    //     (or supfaces in dual mode) as in :meth:`FaceIterator.ignore_subfaces`
    //     and :meth`FaceIterator.ignore_supfaces`.

    //     Those methods add the current face to ``visited_all`` before
    //     visiting sub-/supfaces instead of after. One cannot arbitralily
    //     add faces to ``visited_all``, as visited_all has a maximal length.
    int dim = face_iter[0].max_dimension;
    while (!next_face_loop(&face_iter[0]) && (face_iter[0].current_dimension < dim)){}
    face_iter[0]._index += 1;
    return face_iter[0].current_dimension;
}
//#cdef extern from "<omp.h>":
//#    extern int omp_get_thread_num() nogil



inline int prepare_partial_iter(iter_struct *face_iter, size_t i, size_t *f_vector){
    // r"""
    // Prepares the face iterator to not visit the fatets 0,...,-1
    // """
    size_t j, k;
    size_t current_i = 0;
    if(rec_depth > 0){
        current_i = i/(myPow(face_iter[0].n_coatoms, (rec_depth - 1)));
    }
    size_t rec;
    int d;
    int dimension = face_iter[0].dimension;
    int add_a_face = 1;
    face_iter[0].face = NULL;
    if((current_i != face_iter[0].current_stadium[0])){
        face_iter[0].face = NULL;
        face_iter[0].current_dimension = dimension -1;
        face_iter[0].lowest_dimension = face_iter[0].n_lines;


        face_iter[0].n_visited_all[dimension -1] = 0;
        face_iter[0].n_newfaces[dimension - 1] = face_iter[0].n_coatoms;
        face_iter[0].current_stadium[0] = 0;
        face_iter[0].first_time[dimension - 1] = 1;
        add_a_face = 0;
    }
    for(size_t rec = 0; rec < rec_depth; rec++){
        current_i = i/(myPow(face_iter[0].n_coatoms, (rec_depth - rec - 1)));
        i = i%(myPow(face_iter[0].n_coatoms, (rec_depth - rec - 1)));
        int rec2 = rec;
        if((current_i != face_iter[0].current_stadium[rec]) || (face_iter[0].current_dimension >= dimension - rec2 - 1) || (rec == rec_depth -1)){
            if((rec == rec_depth - 1) && (add_a_face) && (rec > 0)){
                face_iter[0].n_newfaces[dimension-rec-1] += 1;
            }
            if((face_iter[0].current_dimension > dimension - rec2 - 1) || (current_i >= face_iter[0].n_newfaces[dimension-rec-1] + face_iter[0].current_stadium[rec])){
                return 0;
            }
            face_iter[0].current_dimension = dimension - rec2 - 1;
            if(i == 0)
                f_vector[dimension - rec] += 1;

            k = face_iter[0].n_visited_all[dimension - rec-1];
            size_t missing_faces = current_i - face_iter[0].current_stadium[rec];
            for(j = face_iter[0].n_newfaces[dimension-rec-1]- missing_faces; j < face_iter[0].n_newfaces[dimension-rec-1]; j++){
               face_iter[0].visited_all[k] = face_iter[0].newfaces[dimension -rec-1][j];
               k += 1;
            }
            face_iter[0].n_visited_all[dimension - rec-1] += missing_faces;
            face_iter[0].current_stadium[rec] = current_i;
            face_iter[0].current_stadium[rec+1] = 0;

            face_iter[0].n_newfaces[dimension - rec-1] -= missing_faces;
            face_iter[0].first_time[dimension - rec-1] = 1;
            face_iter[0].yet_to_visit = 0;
            face_iter[0].max_dimension = dimension - rec -1;
            add_a_face = 0;
            if(rec < rec_depth -1){
                size_t old_number = face_iter[0].n_newfaces[dimension-rec-1];
                d = next_dimension(face_iter);
                face_iter[0].n_newfaces[dimension-rec-1] = old_number;
                face_iter[0].first_time[dimension - rec-1] = 1;
                face_iter[0].yet_to_visit = 0;
            }
        }
    }
    return 1;
}

inline void partial_f(iter_struct *face_iter, size_t *f_vector, size_t i){
    if (prepare_partial_iter(face_iter, i, f_vector)){
        int d, dimension = face_iter[0].dimension;
        size_t j;
        d = next_dimension(face_iter);
        int rec_depth2 = rec_depth;
        while (d < dimension -rec_depth2){
            f_vector[d + 1] += 1;
            d = next_dimension(face_iter);
        }
    }
}

void parallel_f_vector(iter_struct **face_iter, size_t *f_vector, size_t n_threads, size_t recursion_depth){
    rec_depth = recursion_depth;
    omp_set_num_threads(n_threads);
    size_t **shared_f = new size_t *[n_threads]();
    int dimension = face_iter[0][0].dimension;
    int j;

    for(size_t i = 0; i < n_threads; i++){
        shared_f[i] = new size_t[dimension + 2]();
    }
    size_t n_faces = face_iter[0][0].n_coatoms;

    //iter_struct *my_iter;
    //size_t *my_f;
    //for l in prange(0, n_faces ** rec_depth, schedule='dynamic', chunksize=1):
    #pragma omp parallel for shared(face_iter, shared_f) schedule(dynamic, 1)
    for(size_t l = 0; l < myPow(n_faces,rec_depth); l++){
        //partial_f(face_iter[openmp.omp_get_thread_num()], shared_f[openmp.omp_get_thread_num()], l)
        partial_f(face_iter[omp_get_thread_num()], shared_f[omp_get_thread_num()], l);
    }

    for(size_t i = 0; i < n_threads; i ++){
        for(int j = 0; j < dimension + 2; j++){
            f_vector[j] += shared_f[i][j];
        }
        delete[] shared_f[i];
    }
    delete[] shared_f;
}
