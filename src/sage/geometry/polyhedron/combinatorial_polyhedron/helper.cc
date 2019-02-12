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

// Copyright: see base.pyx

#include <math.h>
#include <cstdint>
#include <stdlib.h>  // for aligned_alloc in C++11
#include <cstdlib>  // for aligned_alloc in C++17
#include <cstdio>




// as of now, 512bit does not have something like _mm256_testc_si256,
// which is the bottle neck of this function,
// so it does not make sene to implement it


#if __AVX2__
    // 256-bit commands
    #include <immintrin.h>
    #define chunktype __m256i
    const size_t chunksize = 256;
    #define bitwise_intersection(one,two) _mm256_and_si256((one),(two))
    // this is supposed to something as (one) & (two)
    #define bitwise_is_not_subset(one,two) !_mm256_testc_si256((two),(one))
    // this is supposed to something as (one) & ~(two)
    #define store_register(one,two) _mm256_store_si256((__m256i*)&(one),(two))
    // this is supposed to be something as `one = two`,
    // where `two` is a register
    // note that `one` needs to be 32-byte-aligned
    #define load_register(one,two) \
        (one) = _mm256_load_si256((const __m256i*)&(two))
    // this is supposed to be something as `one = two`,
    // where `one` is a register
    // note that `two` needs to be 32-byte-aligned

#elif __SSE4_1__
    // 128-bit commands
    // comments see AVX2
    #include <emmintrin.h>
    #include <smmintrin.h>
    #define chunktype __m128i
    const size_t chunksize = 128;
    #define bitwise_intersection(one,two) _mm_and_si128((one),(two))
    #define bitwise_is_not_subset(one,two) !_mm_testc_si128((two),(one))
    #define store_register(one,two) _mm_store_si128((__m128i*)&(one),(two))
    // note that `one` needs to be 16-byte-aligned
    #define load_register(one,two) \
        (one) = _mm_load_si128((const __m128i*)&(two))
    // note that `two` needs to be 16-byte-aligned

#else
    // 64-bit commands, without intrinsics
    #define chunktype uint64_t
    const size_t chunksize = 64;
    #define bitwise_intersection(one,two) (one) & (two)
    #define bitwise_is_not_subset(one,two) (one) & ~(two)
    #define store_register(one,two) one = two
    #define load_register(one,two) one = two
#endif

#if __POPCNT__
    #include <immintrin.h>
    #if INTPTR_MAX == INT64_MAX // 64-bit
        #define popcount(A) _mm_popcnt_u64(A)
    #else // assuming 32-bit
        #define popcount(A) _mm_popcnt_u32(((uint32_t *) &A)[0]) + _mm_popcnt_u32(((uint32_t *) &A)[1])
    #endif
#else
    #define popcount(A) naive_popcount(A)
#endif

inline unsigned int naive_popcount(uint64_t A){
    // popcount without intrinsics
    unsigned int count = 0;
    while (A){
        count += A & 1;
        A >>= 1;
    }
    return count;
}


inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                         size_t face_length){
    // will set C to be the intersection of A and B, i.e.
    // C = A & B
    // `face_length` is the length of A, B and C in terms of uint64_t
    size_t i;
    chunktype a;
    chunktype b;
    chunktype c;
    for (i = 0; i < face_length; i += chunksize/64){
        load_register(a,A[i]);
        load_register(b,B[i]);
        c = bitwise_intersection(a,b);
        store_register(C[i],c);
    }
}

inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
    //returns 1 if A is a proper subset of B, otherwise returns 0,
    // this is done by checking if there is an element in A, which is not in B
    size_t i;
    chunktype a;
    chunktype b;
    for (i = 0; i < face_length; i += chunksize/64){
        load_register(a,A[i]);
        load_register(b,B[i]);
        if (bitwise_is_not_subset(a,b)){
            return 0;
        }
    }
    return 1;
}


inline size_t CountFaceBits(uint64_t* A, size_t face_length) {
    // counts the number of vertices in a face by counting bits set to one
    size_t i;
    unsigned int count = 0;
    for (i=0;i<face_length;i++){
        count += (size_t) popcount(A[i]);
    }
    return count;
}


inline size_t get_next_level(\
        uint64_t **faces, size_t lenfaces, uint64_t **nextfaces, \
        uint64_t **nextfaces2, uint64_t **forbidden, \
        size_t nr_forbidden, size_t face_length){
    // intersects the first `lenfaces - 1` faces of `faces` with'faces[lenfaces-1]`
    // determines which ones are exactly of one dimension less
    // by considering containment
    // newfaces2 will point at those of exactly one dimension less
    // which are not contained in any of the faces in 'forbidden'
    // returns the number of those faces
    const size_t constlenfaces = lenfaces;
    int addfacearray[constlenfaces - 1] = { };
    size_t j,k, addthisface;
    size_t newfacescounter = 0;
    for (j = 0; j < lenfaces - 1; j++){
        intersection(faces[j],faces[lenfaces - 1],nextfaces[j], face_length);
        addfacearray[j] = 1;
    }
    // we have create all possible intersection with the i_th-face, but some of them might not be of exactly one dimension less
    for (j = 0; j< lenfaces-1; j++){
        for(k = 0; k < j; k++){
            // testing if nextfaces[j] is contained in different nextface
            if(is_subset(nextfaces[j],nextfaces[k], face_length)){
                addfacearray[j] = 0;
                break;
                }
            }
        if (!addfacearray[j]) {
            continue;
        }

        for(k = j+1; k < lenfaces-1; k++){
            // testing if nextfaces[j] is contained in a different nextface
            if(is_subset(nextfaces[j],nextfaces[k], face_length)){
            addfacearray[j] = 0;
            break;
            }
        }
        if (!addfacearray[j]) {
            continue;
        }

        for (k = 0; k < nr_forbidden; k++){
            // we do not want to double count any faces,
            // we have visited all faces in forbidden again, so we do not want to do that again
            if(is_subset(nextfaces[j],forbidden[k], face_length)){
                addfacearray[j] = 0;
                break;
            }
        }
    }

    for (j = 0; j < lenfaces -1; j++){
        // let `newfaces2` point to the newfaces we want to consider
        if (!addfacearray[j]) {
            continue;
        }
        nextfaces2[newfacescounter] = nextfaces[j];
        newfacescounter++;
    }
    return newfacescounter;
}


size_t facet_repr_from_bitrep(uint64_t *face, uint64_t **facets, \
                              size_t *output, size_t nr_facets, \
                              size_t face_length){
    // Writes the facet_repr of the current face in output.
    // Returns the length of the representation.
    size_t counter = 0;
    size_t i;
    for (i = 0; i < nr_facets; i++){
        if (is_subset(face, facets[i], face_length)){
            output[counter] = i;
            counter++;
        }
    }
    return counter;
}

void copy_face(uint64_t *input1, uint64_t *output1, \
               size_t face_length){
    size_t i;
    chunktype input;
    for (i = 0; i < face_length; i += chunksize/64){
        load_register(input, input1[i]);
        store_register(output1[i], input);
    }
}
