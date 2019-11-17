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

/*
WARNING: Any Bit-representation is assumed to be `chunksize`-bit aligned.

As of now, there is no 512-bit equivalent to `_mm256_testc_si256`,
which is the bottle neck of this function.
So we do not consider anything newer than AVX2.
*/

// inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length)
#if __AVX__
    /*
    256-bit commands.
    See `#else` for a non-intrinsic version.

    Requires header `immintrin.h`.
    */
    const size_t chunksize = 256;
    inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
        /*
        Return ``A & ~B == 0``.
        A is not subset of B, iff there is a vertex in A, which is not in B.
        ``face_length`` is the length of A and B in terms of uint64_t.

        Note that A,B are assumed to be 32-byte-aligned.
        */
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

    inline int is_subset(uint64_t *A, uint64_t *B, uint32_t* nonzero){
        /*
        Return ``A & ~B == 0``.
        A is not subset of B, iff there is a vertex in A, which is not in B.
        ``face_length`` is the length of A and B in terms of uint64_t.

        Note that A,B are assumed to be 32-byte-aligned.
        */
        uint32_t i,j;
        uint32_t face_length = nonzero[0];
        for (j = 1; j < face_length; j++){
            i = nonzero[j];
            __m256i a = _mm256_load_si256((const __m256i*)&A[i]);
            __m256i b = _mm256_load_si256((const __m256i*)&B[i]);
            if (!_mm256_testc_si256(b, a)){ //need to be opposite order !!
                return 0;
            }
        }
        return 1;
    }

#elif __SSE4_1__
    /*
    128-bit commands.
    See `#else` for a non-intrinsic version.

    Requires headers `smmintrin.h` and `emmintrin.h`.
    */
    const size_t chunksize = 128;
    inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
        /*
        Return ``A & ~B == 0``.
        A is not subset of B, iff there is a vertex in A, which is not in B.
        ``face_length`` is the length of A and B in terms of uint64_t.

        Note that A,B are assumed to be 16-byte-aligned.
        */
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

    inline int is_subset(uint64_t *A, uint64_t *B, uint32_t* nonzero){
        /*
        Return ``A & ~B == 0``.
        A is not subset of B, iff there is a vertex in A, which is not in B.
        ``face_length`` is the length of A and B in terms of uint64_t.

        Note that A,B are assumed to be 32-byte-aligned.
        */
        uint32_t i,j;
        uint32_t face_length = nonzero[0];
        for (j = 1; j < face_length; j++){
            i = nonzero[j];
            __m128i a = _mm_load_si128((const __m128i*)&A[i]);
            __m128i b = _mm_load_si128((const __m128i*)&B[i]);
            if (!_mm_testc_si128(b, a)){ //need to be opposite order !!
                return 0;
            }
        }
        return 1;
    }

#else
    // No intrinsics.
    const size_t chunksize = 64;
    inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
        /*
        Return ``A & ~B == 0``.
        A is not subset of B, iff there is a vertex in A, which is not in B.
        ``face_length`` is the length of A and B in terms of uint64_t.
        */
        size_t i;
        for (i = 0; i < face_length; i++){
            if (A[i] & ~B[i]){
                return 0;
            }
        }
        return 1;
    }

    inline int is_subset(uint64_t *A, uint64_t *B, uint32_t* nonzero){
        /*
        Return ``A & ~B == 0``.
        A is not subset of B, iff there is a vertex in A, which is not in B.
        ``face_length`` is the length of A and B in terms of uint64_t.

        Note that A,B are assumed to be 32-byte-aligned.
        */
        uint32_t i,j;
        uint32_t face_length = nonzero[0];
        for (j = 1; j < face_length; j++){
            i = nonzero[j];
            if (A[i] & ~B[i]){
                return 0;
            }
        }
        return 1;
    }

#endif

// inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C,
//                          size_t face_length)
#if __AVX2__
    /*
    256-bit commands.
    See `#else` for a non-intrinsic version.

    Requires header `immintrin.h`.
    */
    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             size_t face_length){
        /*
        Set ``C = A & B``, i.e. C is the intersection of A and B.
        ``face_length`` is the length of A, B and C in terms of uint64_t.

        Note that A,B,C are assumed to be 32-byte-aligned.
        */
        size_t i;
        for (i = 0; i < face_length; i += 4){
            __m256i a = _mm256_load_si256((const __m256i*)&A[i]);
            __m256i b = _mm256_load_si256((const __m256i*)&B[i]);
            __m256i c = _mm256_and_si256(a, b);
            _mm256_store_si256((__m256i*)&C[i],c);
        }
    }

    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             uint32_t * nonzero_C, size_t face_length){
        /*
        Set ``C = A & B``, i.e. C is the intersection of A and B.
        ``face_length`` is the length of A, B and C in terms of uint64_t.

        Note that A,B,C are assumed to be 32-byte-aligned.
        */
        uint32_t i;
        uint32_t counter = 1;
        for (i = 0; i < face_length; i += 4){
            __m256i a = _mm256_load_si256((const __m256i*)&A[i]);
            __m256i b = _mm256_load_si256((const __m256i*)&B[i]);
            __m256i c = _mm256_and_si256(a, b);
            if (!_mm256_testz_si256(c,c)){
                nonzero_C[counter] = i;
                counter ++;
            }
            _mm256_store_si256((__m256i*)&C[i],c);
        }
        nonzero_C[0] = counter;
    }

#elif __AVX__
    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             size_t face_length){
        /*
        Set ``C = A & B``, i.e. C is the intersection of A and B.
        ``face_length`` is the length of A, B and C in terms of uint64_t.

        Note that A,B,C are assumed to be 16-byte-aligned.
        */
        size_t i;
        for (i = 0; i < face_length; i += 2){
            __m128i a = _mm_load_si128((const __m128i*)&A[i]);
            __m128i b = _mm_load_si128((const __m128i*)&B[i]);
            __m128i c = _mm_and_si128(a, b);
            _mm_store_si128((__m128i*)&C[i],c);
        }
    }

    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             uint32_t * nonzero_C, size_t face_length){
        /*
        Set ``C = A & B``, i.e. C is the intersection of A and B.
        ``face_length`` is the length of A, B and C in terms of uint64_t.

        Note that A,B,C are assumed to be 32-byte-aligned.
        */
        uint32_t i;
        uint32_t counter = 1;
        for (i = 0; i < face_length; i += 4){
            __m128i a = _mm_load_si128((const __m128i*)&A[i]);
            __m128i b = _mm_load_si128((const __m128i*)&B[i]);
            __m128i c = _mm_and_si128(a, b);
            int not_zero = 1 - _mm_testz_si128(c, c);
            _mm_store_si128((__m128i*)&C[i],c);

            a = _mm_load_si128((const __m128i*)&A[i+2]);
            b = _mm_load_si128((const __m128i*)&B[i+2]);
            c = _mm_and_si128(a, b);
            not_zero += 1 - _mm_testz_si128(c, c);
            _mm_store_si128((__m128i*)&C[i+2],c);

            if (not_zero){
                nonzero_C[counter] = i;
                counter ++;
            }
        }
        nonzero_C[0] = counter;
    }

#elif __SSE4_1__
    /*
    128-bit commands.
    See `#else` for a non-intrinsic version.

    Requires header `emmintrin.h`.

    Actually this is already defined for SSE2,
    but we have chosen chunksize according to `is_subset`.

    So 16-byte-alignment is only given for SSE4.1 and greater.
    */
    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             size_t face_length){
        /*
        Set ``C = A & B``, i.e. C is the intersection of A and B.
        ``face_length`` is the length of A, B and C in terms of uint64_t.

        Note that A,B,C are assumed to be 16-byte-aligned.
        */
        size_t i;
        for (i = 0; i < face_length; i += 2){
            __m128i a = _mm_load_si128((const __m128i*)&A[i]);
            __m128i b = _mm_load_si128((const __m128i*)&B[i]);
            __m128i c = _mm_and_si128(a, b);
            _mm_store_si128((__m128i*)&C[i],c);
        }
    }

    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             uint32_t * nonzero_C, size_t face_length){
        /*
        Set ``C = A & B``, i.e. C is the intersection of A and B.
        ``face_length`` is the length of A, B and C in terms of uint64_t.

        Note that A,B,C are assumed to be 32-byte-aligned.
        */
        uint32_t i;
        uint32_t counter = 1;
        for (i = 0; i < face_length; i += 2){
            __m128i a = _mm_load_si128((const __m128i*)&A[i]);
            __m128i b = _mm_load_si128((const __m128i*)&B[i]);
            __m128i c = _mm_and_si128(a, b);
            _mm_store_si128((__m128i*)&C[i],c);
            if (!_mm_testz_si128(c, c)){
                nonzero_C[counter] = i;
                counter ++;
            }
            _mm_store_si128((__m128i*)&C[i],c);
        }
        nonzero_C[0] = counter;
    }

#else
    // No intrinsics.
    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             size_t face_length){
        /*
        Set ``C = A & B``, i.e. C is the intersection of A and B.
        ``face_length`` is the length of A, B and C in terms of uint64_t.
        */
        size_t i;
        for (i = 0; i < face_length; i++){
            C[i] = A[i] & B[i];
        }
    }

    inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                             uint32_t * nonzero_C, size_t face_length){
        /*
        Set ``C = A & B``, i.e. C is the intersection of A and B.
        ``face_length`` is the length of A, B and C in terms of uint64_t.

        Note that A,B,C are assumed to be 32-byte-aligned.
        */
        uint32_t i;
        uint32_t counter = 1;
        for (i = 0; i < face_length; i++){
            C[i] = A[i] & B[i];
            if (C[i]){
                nonzero_C[counter] = i;
                counter ++;
            }
        }
        nonzero_C[0] = counter;
    }

#endif


// inline size_t CountFaceBits(uint64_t* A, size_t face_length)
#if (__POPCNT__) && (INTPTR_MAX == INT64_MAX)
    /*
    Assumes a 64-bit machine and POPCNT

    Requires header `immintrin.h`.
    */
    inline size_t count_atoms(uint64_t* A, size_t face_length) {
        /*
        Return the number of atoms/vertices in A.
        This is the number of set bits in A.
        ``face_length`` is the length of A in terms of uint64_t.
        */
        size_t i;
        unsigned int count = 0;
        for (i=0; i<face_length; i++){
            count += (size_t) _mm_popcnt_u64(A[i]);
        }
        return count;
    }

#else // popcount without intrinsics
    inline size_t count_atoms(uint64_t* A, size_t face_length) {
        /*
        Return the number of atoms/vertices in A.
        This is the number of set bits in A.
        ``face_length`` is the length of A in terms of uint64_t.
        */
        size_t i;
        unsigned int count = 0;
        for (i=0; i<face_length; i++){
            uint64_t a = A[i];
            a = a - ((a >> 1) & 0x5555555555555555ULL);
            a = (a & 0x3333333333333333ULL) + ((a >> 2) & 0x3333333333333333ULL);
            count += ( ((a + (a >> 4)) & 0x0f0f0f0f0f0f0f0fULL) * 0x0101010101010101ULL ) >> 56;
        }
        return count;
    }

#endif

size_t get_next_level(\
        uint64_t **faces, const size_t n_faces, uint64_t **maybe_newfaces, \
        uint64_t **newfaces, uint64_t **visited_all, \
        size_t n_visited_all, size_t face_length, int *is_not_newface, uint64_t **LHS, uint64_t **RHS){
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

    uint64_t *current_LHS = NULL;
    uint64_t *current_RHS = NULL;
    uint64_t *next_LHS = NULL;
    uint64_t *next_RHS = NULL;
    if (LHS){
        current_LHS = LHS[1];
        current_RHS = RHS[1];
        next_LHS = LHS[0];
        next_RHS = RHS[0];
    }


    // Step 1:
    for (size_t j = 0; j < n_faces - 1; j++){
        intersection(faces[j], faces[n_faces - 1], maybe_newfaces[j], face_length);
        if (LHS){
            next_LHS[j] = current_LHS[j] | current_LHS[n_faces -1];
            next_RHS[j] = current_RHS[j] | current_RHS[n_faces -1];
        }
    }

    // For each face we will Step 2 and Step 3.
    for (size_t j = 0; j < n_faces-1; j++){
        is_not_newface[j] = 0;
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
        if (LHS){
            next_LHS[n_newfaces] = next_LHS[j];
            next_RHS[n_newfaces] = next_RHS[j];
        }
        n_newfaces++;
    }
    return n_newfaces;
}

size_t get_next_level_with_nonzero(\
        uint64_t **faces, const size_t n_faces, uint64_t **maybe_newfaces, \
        uint64_t **newfaces, uint64_t **visited_all, \
        size_t n_visited_all, size_t face_length, int *is_not_newface, uint64_t **LHS, uint64_t **RHS, \
        uint32_t **nonzero_maybe_newfaces, uint32_t **nonzero_newfaces){
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

    uint64_t *current_LHS = NULL;
    uint64_t *current_RHS = NULL;
    uint64_t *next_LHS = NULL;
    uint64_t *next_RHS = NULL;
    if (LHS){
        current_LHS = LHS[1];
        current_RHS = RHS[1];
        next_LHS = LHS[0];
        next_RHS = RHS[0];
    }


    // Step 1:
    for (size_t j = 0; j < n_faces - 1; j++){
        intersection(faces[j], faces[n_faces - 1], maybe_newfaces[j], nonzero_maybe_newfaces[j], face_length);
        if (LHS){
            next_LHS[j] = current_LHS[j] | current_LHS[n_faces -1];
            next_RHS[j] = current_RHS[j] | current_RHS[n_faces -1];
        }
    }

    // For each face we will Step 2 and Step 3.
    for (size_t j = 0; j < n_faces-1; j++){
        is_not_newface[j] = 0;
        // Step 2a:
        for(size_t k = 0; k < j; k++){
            // Testing if maybe_newfaces[j] is contained in different nextface.
            if(is_subset(maybe_newfaces[j], maybe_newfaces[k], nonzero_maybe_newfaces[j])){
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
            if(is_subset(maybe_newfaces[j],maybe_newfaces[k], nonzero_maybe_newfaces[j])){
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
            if(is_subset(maybe_newfaces[j], visited_all[k], nonzero_maybe_newfaces[j])){
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
        nonzero_newfaces[n_newfaces] = nonzero_maybe_newfaces[j];
        if (LHS){
            next_LHS[n_newfaces] = next_LHS[j];
            next_RHS[n_newfaces] = next_RHS[j];
        }
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

size_t bit_repr_to_coatom_repr(uint64_t *face, uint64_t **coatoms, \
                               size_t n_coatoms, uint32_t *nonzero_face, \
                               uint8_t *output){
    /*
    Write the coatom-representation of face in output. Return length.
    ``face_length`` is the length of ``face`` and ``coatoms[i]``
    in terms of uint64_t.
    ``n_coatoms`` length of ``coatoms``.
    */
    uint8_t count_length = 0;
    for (uint8_t i = 0; i < n_coatoms; i++){
        if (is_subset(face, coatoms[i], nonzero_face)){
            // ``face`` is contain in ``coatoms[i]``,
            // then ``i`` is an element in the coatom-represention.
            output[count_length + 2] = i;
            count_length++;
        }
    }
    output[0] = count_length;
    return count_length;
}

inline size_t is_bad_face_cc(uint64_t *face, uint32_t *nonzero_face, int dimension, uint64_t ** coatoms, size_t n_coatoms, size_t face_length, uint64_t *LHS, uint64_t *RHS, uint64_t *current_LHS, uint64_t *current_RHS, size_t *pte, uint8_t *output){
    size_t i;
    for(i=0; i< n_coatoms; i++){
        if ((LHS[i] &~current_LHS[0]) + (RHS[i] & ~current_RHS[0])){
            if(is_subset(face, coatoms[i], nonzero_face)){
                current_LHS[0] = current_LHS[0] | LHS[i];
                current_RHS[0] = current_RHS[0] | RHS[i];
            }
        }
    }
    size_t m = dimension + 2;
    size_t e = m-1-count_atoms(current_RHS, 1) + 1;
    size_t t = m-1-count_atoms(current_LHS, 1);

    if (e > t)
        return 0;

    if (2*e >= m)
        return 0;
    pte[0] = e;
    size_t counter = 0;
    output[1] = e;
    return bit_repr_to_coatom_repr(face, coatoms, n_coatoms, nonzero_face, output);
}
