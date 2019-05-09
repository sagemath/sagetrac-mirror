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


// As of now, 512bit does not have something like _mm256_testc_si256,
// which is the bottle neck of the algorithm,
// so it does not make sense to implement it.


// Determine the best way to check for inclusion/subset.
#if __AVX__
    // 256-bit commands, those operations are equivalent to the operations
    // defined in `#else`
    // (intrinsics defined in immintrin.h)
    const size_t chunksize = 256;
    inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
        // A & ~B == 0
        // Return 1 if A is a subset of B, otherwise returns 0.

        // This done by checking if there is an element in A,
        // which is not in B
        // `face_length` is the length of A and B in terms of uint64_t
        // Note that A,B need to be 32-Byte-aligned
        size_t i;
        for (i = 0; i < face_length; i += 4){
            __m256i a = _mm256_load_si256((const __m256i*)&A[i]);
            __m256i b = _mm256_load_si256((const __m256i*)&B[i]);
            if (!_mm256_testc_si256(b, a)){ //need to be opposite order !!
                return 0;
            }
        }
        return 1;
    }

#elif __SSE4_1__
    // 128-bit commands, those operations are equivalent to the operations
    // defined in `#else`
    // (intrinsics defined in smmintrin.h and emmintrin.h)
    const size_t chunksize = 128;
    inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
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
                return 0;
            }
        }
        return 1;
    }

#else
    // no intrinsics
    const size_t chunksize = 64;
    inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
        // A & ~B == 0
        // Return 1 if A is a subset of B, otherwise returns 0.

        // This done by checking if there is an element in A,
        // which is not in B
        // `face_length` is the length of A and B in terms of uint64_t
        size_t i;
        for (i = 0; i < face_length; i++){
            if (A[i] & ~B[i]){
                return 0;
            }
        }
        return 1;
    }

#endif

// Determine the best way to calculate intersection.
#if __AVX2__
    // 256-bit commands, those operations are equivalent to the operations
    // defined in `#else`
    // intrinsics defined in immintrin.h
    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             size_t face_length){
        // C = A & B
        // Set C to be the intersection of A and B.
        // ``face_length`` is the length of A, B and C in terms of ``uint64_t``.
        // Note that A,B,C need to be 32-Byte-aligned
        size_t i;
        for (i = 0; i < face_length; i += 4){
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

size_t get_next_level(\
        uint64_t **faces, const size_t n_faces, uint64_t **maybe_newfaces, \
        uint64_t **newfaces, uint64_t **visited_all, \
        size_t n_visited_all, size_t face_length){
    /*
    Set ``newfaces`` to be the facets of ``faces[n_faces -1]``
    that are not contained in a face of ``visited_all``.

    INPUT:

    - ``maybe_newfaces`` -- quasi of type ``uint64_t[n_faces -1][face_length]``,
      needs to be ``chunksize``-Bit aligned
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

    // Step 1:
    for (size_t j = 0; j < n_faces - 1; j++){
        intersection(faces[j], faces[n_faces - 1], maybe_newfaces[j], face_length);
    }

    // We keep track, which face in ``maybe_newfaces`` is a new face.
    int *is_not_newface = new int[n_faces -1]();

    // For each face we will Step 2 and Step 3.
    for (size_t j = 0; j < n_faces-1; j++){
        // Step 2a:
        for(size_t k = 0; k < j; k++){
            // Testing if maybe_newfaces[j] is contained in different nextface.
            if(is_subset(maybe_newfaces[j], maybe_newfaces[k],face_length)){
                // If so, it is not inclusion-maximal and hence not of codimension 1.
                is_not_newface[j] = 1;
                break;
                }
            }
        if (is_not_newface[j]) {
            // No further tests needed, if it is not of codimension 1.
            continue;
        }

        // Step 2b:
        for(size_t k = j+1; k < n_faces-1; k++){
            // Testing if maybe_newfaces[j] is contained in different nextface.
            if(is_subset(maybe_newfaces[j],maybe_newfaces[k], face_length)){
                // If so, it is not inclusion-maximal and hence not of codimension 1.
                is_not_newface[j] = 1;
            break;
            }
        }
        if (is_not_newface[j]) {
            // No further tests needed, if it is not of codimension 1.
            continue;
        }

        // Step 3:
        for (size_t k = 0; k < n_visited_all; k++){
            // Testing if maybe_newfaces[j] is contained in one,
            // we have already completely visited.
            if(is_subset(maybe_newfaces[j], visited_all[k], face_length)){
                // If so, we don't want to revisit.
                is_not_newface[j] = 1;
                break;
            }
        }
    }

    // Set ``newfaces`` to point to the correct ones.
    size_t n_newfaces = 0;  // length of newfaces2
    for (size_t j = 0; j < n_faces -1; j++){
        if (is_not_newface[j]) {
            // Not a new face of codimension 1.
            continue;
        }
        // It is a new face of codimension 1.
        newfaces[n_newfaces] = maybe_newfaces[j];
        n_newfaces++;
    }
    delete[] is_not_newface;
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
