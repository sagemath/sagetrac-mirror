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
// which is the bottle neck of this function, so it does not make sene to implement it
// 512 bit commands
// #if __AVX512F__
// #define chunktype __m512i
// #define bitwise_intersection(one,two) _mm512_and_si512((one),(two))


#if __AVX2__
    //256 bit commands
    #include <immintrin.h>
    #define chunktype __m256i
    const unsigned int chunksize = 256;
    #define bitwise_intersection(one,two) _mm256_and_si256((one),(two)) // this is supposed to something as (one) & (two)
    #define bitwise_is_not_subset(one,two) !_mm256_testc_si256((two),(one)) // this is supposed to something as (one) & ~(two)
    #define store_register(one,two) _mm256_store_si256((__m256i*)&(one),(two)) //this is supposed to be something as one = two, where two is a register
    #define load_register(one,two) (one) = _mm256_load_si256((const __m256i*)&(two)) //this is supposed to be somethign as one = two, where one is a register

#elif __SSE4_1__
    //128 bit commands
    #include <emmintrin.h>
    #include <smmintrin.h>
    #define chunktype __m128i
    const unsigned int chunksize = 128;
    #define bitwise_intersection(one,two) _mm_and_si128((one),(two))
    #define bitwise_is_not_subset(one,two) !_mm_testc_si128((two),(one))
    #define store_register(one,two) _mm_store_si128((__m128i*)&(one),(two))
    #define load_register(one,two) (one) = _mm_load_si128((const __m128i*)&(two))

#else
    //64 bit commands
    #define chunktype uint64_t
    const unsigned int chunksize = 64;
    #define bitwise_intersection(one,two) (one) & (two)
    #define bitwise_is_not_subset(one,two) (one) & ~(two)
    #define store_register(one,two) one = two
    #define load_register(one,two) one = two
#endif

#if __POPCNT__
    #include <immintrin.h>
    #if INTPTR_MAX == INT64_MAX //64-bit
        #define popcount(A) _mm_popcnt_u64(A)
    #else //assuming 32-bit
        #define popcount(A) _mm_popcnt_u32(((uint32_t *) &A)[0]) + _mm_popcnt_u32(((uint32_t *) &A)[1])
    #endif
#else
    #define popcount(A) naive_popcount(A)
#endif



static uint64_t vertex_to_bit_dictionary[64];
//this dictionary helps storing a vector of 64 or 32 incidences as uint64_t or uint32_t,
//where each bit represents an incidence


void build_dictionary(){
    unsigned int i = 0;
    uint64_t count = 1;
    for (i=0; i< 64;i++){
        vertex_to_bit_dictionary[64 -i-1] = count;
        count *= 2;
    }
}

inline unsigned int naive_popcount(uint64_t A){
    unsigned int count = 0;
    while (A){
        count += A & 1;
        A >>= 1;
    }
    return count;
}



inline void intersection(uint64_t *A, uint64_t *B, uint64_t *C, \
                         size_t face_length){
    // will set C to be the intersection of A and B
    size_t i;
    chunktype a;
    chunktype b;
    chunktype c;
    for (i = 0; i < face_length; i++){
        load_register(a,A[i*chunksize/64]);
        load_register(b,B[i*chunksize/64]);
        c = bitwise_intersection(a,b);
        store_register(C[i*chunksize/64],c);
    }
}

inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
    //returns 1 if A is a proper subset of B, otherwise returns 0,
    // this is done by checking if there is an element in A, which is not in B
    size_t i;
    chunktype a;
    chunktype b;
    for (i = 0; i < face_length; i++){
        load_register(a,A[i*chunksize/64]);
        load_register(b,B[i*chunksize/64]);
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
    for (i=0;i<face_length*chunksize/64;i++){
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


void char_from_incidence_list(unsigned int *incidence_list, \
                              size_t nr_vertices, \
                              uint64_t *output, size_t face_length){
    size_t entry, position, value, i;
    for (i = 0; i < face_length*chunksize/64; i++){
        output[i] = 0;
    }
    while (nr_vertices--) {
        entry = incidence_list[nr_vertices];
        if (entry){
            value = nr_vertices % 64;
            position = nr_vertices/64;
            output[position] += vertex_to_bit_dictionary[value];
        }
    }
}

void char_from_array(unsigned int* input, unsigned int len, \
                     uint64_t *output, size_t face_length){
    size_t entry, position, value, i;
    for (i = 0; i < face_length*chunksize/64; i++){
        output[i] = 0;
    }
    while (len--) {
        entry = input[len];
        value = entry % 64;
        position = entry/64;
        output[position] += vertex_to_bit_dictionary[value];
    }
}

void make_trivial_face(size_t nr_vertices, uint64_t * output, \
                       size_t face_length){
    build_dictionary();
    size_t entry, position, value, i;
    for (i = 0; i < face_length*chunksize/64; i++){
        output[i] = 0;
    }
    for (entry = 0; entry < nr_vertices; entry++){
        value = entry % 64;
        position = entry/64;
        output[position] += vertex_to_bit_dictionary[value];
    }
}

void get_facets_from_incidence_matrix(\
        unsigned int **incidence_matrix, uint64_t **facets, \
        size_t nr_vertices, size_t nr_facets){
    build_dictionary();
    size_t i;
    size_t length_of_face = ((nr_vertices - 1)/chunksize + 1);
    //this determines the length of the face in terms of chunktype
    for(i = 0; i<nr_facets; i++){
        char_from_incidence_list(incidence_matrix[i], nr_vertices, \
                                 facets[i], length_of_face);
    }
}


void get_vertices_from_incidence_matrix(\
        unsigned int **incidence_matrix, uint64_t **vertices, \
        size_t nr_vertices, size_t nr_facets){
    build_dictionary();
    size_t i, j;
    size_t face_length = ((nr_facets - 1)/chunksize + 1);
    //this determines the length in facet representation in terms of chunktype

    const size_t const_facets = nr_facets;
    unsigned int new_facets_array[const_facets];
    unsigned int length_that_face;
    for(i = 0; i< nr_vertices; i++){
        length_that_face = 0;
        for (j=0; j < nr_facets; j++){
            if (incidence_matrix[j][i]){
                new_facets_array[length_that_face] = j;
                length_that_face++;
            }
        }
        char_from_array(new_facets_array, length_that_face, \
                        vertices[i], face_length);
    }
}


void get_facets_bitrep_from_facets_pointer( \
        unsigned int ** facets_input, unsigned int *len_facets, \
        uint64_t ** facets_output, size_t nr_vertices, \
        size_t nr_facets){
    size_t i;
    size_t face_length = ((nr_vertices - 1)/chunksize + 1);
    //this determines the length of the face in terms of chunktype
    for(i = 0; i<nr_facets; i++){
        char_from_array(facets_input[i], len_facets[i], \
                        facets_output[i], face_length);
    }
}

void get_vertices_bitrep_from_facets_pointer( \
        unsigned int ** facets_input, unsigned int *len_facets, \
        uint64_t ** vertices_output, size_t nr_vertices, \
        size_t nr_facets){
    size_t i,j;
    size_t face_length = ((nr_facets - 1)/chunksize + 1);
    //this determines the length in facet representation in terms of chunktype
    const size_t size_three = nr_facets;
    unsigned int *old_facets_walker = new unsigned int [size_three]();
    unsigned int new_facets_array[size_three];
    unsigned int length_that_face;
    for(i = 0; i<nr_vertices; i++){
        length_that_face = 0;
        for (j=0; j < nr_facets; j++){
            if (i == facets_input[j][old_facets_walker[j]]){
                //testing if vertex i is contained in the j-th facet
                new_facets_array[length_that_face] = j;
                old_facets_walker[j]++;
                length_that_face++;
                if (old_facets_walker[j] >= len_facets[j])
                    old_facets_walker[j]--;
            }
        }
        char_from_array(new_facets_array, length_that_face, \
                        vertices_output[i], face_length);
    }
    delete[] old_facets_walker;
}

size_t facet_repr_from_bitrep(uint64_t *face, uint64_t **facets, \
                              size_t *output, size_t nr_facets, \
                              size_t length_of_face){
    // Writes the facet_repr of the current face in output.
    // Returns the length of the representation.
    size_t counter = 0;
    size_t i;
    for (i = 0; i < nr_facets; i++){
        if (is_subset(face, facets[i], length_of_face)){
            output[counter] = i;
            counter++;
        }
    }
    return counter;
}

size_t vertex_repr_from_bitrep(uint64_t *face, size_t *output, \
                               size_t length_of_face){
    build_dictionary();
    size_t i,j;
    size_t counter = 0;
    uint64_t copy;
    for (i = 0; i < length_of_face*chunksize/64; i++){
        if (face[i]){
            copy = face[i];
            for (j = 0; j < 64; j++){
                if (copy >= vertex_to_bit_dictionary[j]){
                    output[counter] = i*64 + j;
                    counter++;
                    copy -= vertex_to_bit_dictionary[j];
                }
            }
        }
    }
    return counter;
}

void copy_face(uint64_t *input1, uint64_t *output1, \
               size_t length_of_face){
    size_t i;
    chunktype input;
    for (i = 0; i < length_of_face; i++){
        load_register(input, input1[i*chunksize/64]);
        store_register(output1[i*chunksize/64], input);
    }
}


inline int is_smaller(uint64_t * one, uint64_t * two, \
                      size_t length_of_face1){
    // returns 1 if `one1` is smaller than `two1`
    // otherwise returns 0
    size_t i;
    for (i = 0; i < length_of_face1; i++){
        if (one[i] < two[i])
            return 1;
        if (two[i] < one[i])
            return 0;
    }
    return 0;
}

void sort_pointers_loop(uint64_t **input, uint64_t **output1, \
                        uint64_t **output2, size_t nr_faces, \
                        size_t length_of_face1){
    // this is mergesort
    // sorts the faces in input and returns them in output1
    // BEWARE: Input is the same as output1 or output2
    size_t middle = nr_faces/2;
    size_t other = nr_faces - middle;
    size_t i = 0, j = middle, counter = 0;
    if (nr_faces == 1){
        output1[0] = input[0];
        return;
    }
    sort_pointers_loop(input, output2, output1, middle, length_of_face1);
    sort_pointers_loop(&(input[middle]), &(output2[middle]), \
                       &(output1[middle]), other, length_of_face1);
    while ((i < middle) and (j < nr_faces)){
        if (is_smaller(output2[i],output2[j],length_of_face1)){
            output1[counter] = output2[i];
            i++;
            counter++;
        } else {
            output1[counter] = output2[j];
            j++;
            counter++;
        }
    }
    if (i < middle){
        while (i < middle){
            output1[counter] = output2[i];
            i++;
            counter++;
        }
    } else {
        while (j < nr_faces){
            output1[counter] = output2[j];
            j++;
            counter++;
        }
    }
}


void sort_pointers(uint64_t **input, size_t nr_faces, \
                   size_t length_of_face){
    // sorts nr_faces in `input` according to their values
    size_t length_of_face1 = length_of_face*chunksize/64;
    const size_t const_nr_faces = nr_faces;
    uint64_t *input2[const_nr_faces];
    sort_pointers_loop(input, input, input2, nr_faces, length_of_face1);
}


size_t find_face(uint64_t **list, uint64_t *face, size_t nr_faces, \
                 size_t length_of_face){
    // finds the position of `face` in `list` (assumed to be sorted)
    // does not check, if the result is actually correct, i.e. if `face`
    // is even an element of `list`
    // use `is_equal` for that
    size_t length_of_face1 = length_of_face*chunksize/64;
    size_t start = 0;
    size_t middle;
    while (nr_faces > 1){
        middle = nr_faces/2;
        if (is_smaller(face, list[middle + start], length_of_face1)){
            nr_faces = middle;
        } else {
            nr_faces -= middle;
            start += middle;
        }
    }
    return start;
}

int is_equal(uint64_t *one, uint64_t *two, size_t length_of_face){
    size_t i;
    size_t length_of_face1 = length_of_face*chunksize/64;
    for (i = 0; i < length_of_face1; i++){
        if (one[i] != two[i])
            return 0;
    }
    return 1;
}
