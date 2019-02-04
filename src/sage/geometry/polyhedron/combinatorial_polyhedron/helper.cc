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

/*
#if __AVX2__
    //256 bit commands
    #include <immintrin.h>
    #define chunktype __m256i
    const unsigned int chunksize = 256;
    #define bitwise_intersection(one,two) _mm256_and_si256((one),(two)) // this is supposed to something as (one) & (two)
    #define bitwise_is_not_subset(one,two) !_mm256_testc_si256((two),(one)) // this is supposed to something as (one) & ~(two)
    #define store_register(one,two) _mm256_storeu_si256((__m256i*)&(one),(two)) //this is supposed to be something as one = two, where two is a register
    #define load_register(one,two) (one) = _mm256_loadu_si256((const __m256i*)&(two)) //this is supposed to be somethign as one = two, where one is a register
    #define leading_zero_count(one) leading_zero_workaround(one) //the workaround is not extremely fast, but it is not needed often (only for edges)
    #define trailing_zero_count(one) trailing_zero_workaround(one) //the workaround is not extremely fast, but it is not needed often

#elif __SSE4_1__
    //128 bit commands
    #include <emmintrin.h>
    #include <smmintrin.h>
    #define chunktype __m128i
    const unsigned int chunksize = 128;
    #define bitwise_intersection(one,two) _mm_and_si128((one),(two))
    #define bitwise_is_not_subset(one,two) !_mm_testc_si128((two),(one))
    #define store_register(one,two) _mm_storeu_si128((__m128i*)&(one),(two))
    #define load_register(one,two) (one) = _mm_loadu_si128((const __m128i*)&(two))
    #define leading_zero_count(one) leading_zero_workaround(one)
    #define trailing_zero_count(one) trailing_zero_workaround(one)

#elif INTPTR_MAX == INT64_MAX
    //64 bit commands
    #define chunktype uint64_t
    const unsigned int chunksize = 64;
    #define bitwise_intersection(one,two) (one) & (two)
    #define bitwise_is_not_subset(one,two) (one) & ~(two)
    #define store_register(one,two) one = two
    #define load_register(one,two) one = two
    #define leading_zero_count(one) leading_zero_naive3(one)
    #define trailing_zero_count(one) trailing_zero_naive3(one)

#else
    //32 bit commands
    #define chunktype uint32_t
    const unsigned int chunksize = 32;
    #define bitwise_intersection(one,two) (one) & (two)
    #define bitwise_is_not_subset(one,two) (one) & ~(two)
    #define store_register(one,two) one = two
    #define load_register(one,two) one = two
    #define leading_zero_count(one) leading_zero_naive3(one)
    #define trailing_zero_count(one) trailing_zero_naive3(one)
#endif

#if __POPCNT__
    #include <immintrin.h>
    #if INTPTR_MAX == INT64_MAX //64-bit
        #define popcount(A) _mm_popcnt_u64(A)
    #else //assuming 32-bit
        #define popcount(A) _mm_popcnt_u32(A)
    #endif
#else
    #define popcount(A) naive_popcount(A)
#endif


#if (__GNUC__ >= 5)  // checking if GCC has aligned_alloc, this should always be the case for sages build in gcc
#define free_aligned(one) free(one)
#else //otherwise falling back to a manual approach
#define aligned_alloc(one,two) aligned_malloc_workaround(two,one)
#define free_aligned(one) aligned_free_workaround(one)
#endif


const unsigned int maxnumberedges = 16348;//^2
//(the edges will be build as an array of arrays,
//such that we can save up to maxnumberedges*maxnumberedges edges,
//the number should contain a high power of two

const unsigned int maxnumberincidences = 16348;//^2
//the maximal number of incidences between l-faces and k-faces

static uint64_t_or_uint32_t vertex_to_bit_dictionary[bit64or32];
//this dictionary helps storing a vector of 64 or 32 incidences as uint64_t or uint32_t,
//where each bit represents an incidence


void build_dictionary(){
    unsigned int i = 0;
    uint64_t_or_uint32_t count = 1;
    for (i=0; i< bit64or32;i++){
        vertex_to_bit_dictionary[bit64or32 -i-1] = count;
        count *= 2;
    }
}


inline unsigned int leading_zero_naive3(uint64_t_or_uint32_t x){
    //taken from https://codingforspeed.com/counting-the-number-of-leading-zeros-for-a-32-bit-integer-signed-or-unsigned/
    //counts the number of leading zero bits of an uint64_t
    unsigned n = 0;
    if (x == 0) return bit64or32;
    while (1) {
        if (x > vertex_to_bit_dictionary[0]) break;
        n++;
        x <<= 1;
    }
    return n;
}

inline unsigned int leading_zero_workaround(chunktype chunk){
    //counts the number of leading zero bits of a chunktype,
    //where chunktype represents 1,2 or 4 uint64_t or uint32_t depending on the processor
    unsigned int i;
    unsigned int count = 0;
    uint64_t_or_uint32_t A[chunksize/bit64or32];
    store_register(A[0],chunk);
    for (i = 0;i < chunksize/bit64or32;i++){
        count += leading_zero_naive3(A[i]);
        if (count < bit64or32*(i+1)){
            return count;
        }
    }
    return count;
}

inline unsigned int trailing_zero_naive3(uint64_t_or_uint32_t x){
    //counts the number of trailing zero bits of an uint64_t_or_uint32_t
    unsigned n = 0;
    if (x == 0) return bit64or32;
    while (1) {
        if (x % 2) break;
        n ++;
        x >>= 1;
    }
    return n;
}

inline unsigned int trailing_zero_workaround(chunktype chunk){
    //counts the number of trailing zero bits of a chunktype,
    //where chunktype represents 1,2 or 4 uint64_t depending on the processor
    unsigned int i;
    unsigned int count = 0;
    uint64_t_or_uint32_t A[chunksize/bit64or32];
    store_register(A[0],chunk);
    for (i = 0;i < chunksize/bit64or32;i++){
        count += trailing_zero_naive3(A[chunksize/bit64or32-i-1]);
        if (count < bit64or32*(i+1)){
            return count;
        }
    }
    return count;
}

inline unsigned int naive_popcount(uint64_t_or_uint32_t A){
    unsigned int count = 0;
    while (A){
        count += A & 1;
        A >>= 1;
    }
    return count;
}


void * aligned_malloc_workaround(size_t size, int align) {
    //taken from https://github.com/xjw/cpp/blob/master/cpp/memory_alignment.cpp
    //they are a workaround in case that C11 is not available

    // alignment could not be less than 0
    if (size<0) {
        return NULL;
    }
    // allocate necessary memory for
    // alignment +
    // area to store the address of memory returned by malloc
    void *p = malloc(size + align-1 + sizeof(void *));
    if (p == NULL) {
        return NULL;
    }
    // address of the aligned memory according to the align parameter
    void *ptr = (void *) (((unsigned long)p + sizeof(void *) + align-1) & ~(align-1));

    // store th address of mallc() above at the beginning of our total memory area
    *((void **)ptr -1) = p;

    // return the address of aligned memory
    return ptr;
}

void aligned_free_workaround(void *p) {
    //taken from https://github.com/xjw/cpp/blob/master/cpp/memory_alignment.cpp

    // Get address of the memory from start of total memory area
    free ( *( (void **)p - 1) );
}
*/



inline void intersection(void *A1, void *B1, void *C1, size_t face_length){
    // will set C to be the intersection of A and B
    size_t i;
    chunktype *A = (chunktype *) A1;
    chunktype *B = (chunktype *) B1;
    chunktype *C = (chunktype *) C1;
    for (i = 0; i < face_length; i++){
        C[i] = bitwise_intersection(A[i],B[i]);
    }
}

inline int is_subset(void *A1, void *B1, size_t face_length){
    //returns 1 if A is a proper subset of B, otherwise returns 0,
    // this is done by checking if there is an element in A, which is not in B
    size_t i;
    chunktype *A = (chunktype *) A1;
    chunktype *B = (chunktype *) B1;
    for (i = 0; i < face_length; i++){
        if (bitwise_is_not_subset(A[i],B[i])){
            return 0;
        }
    }
    return 1;
}


inline size_t CountFaceBits(void* A2, size_t face_length) {
    // counts the number of vertices in a face by counting bits set to one
    size_t i;
    unsigned int count = 0;
    chunktype *A1 = (chunktype *) A2;
    const unsigned int length_of_conversion_face = face_length*chunksize/64;
    uint64_t_or_uint32_t A[length_of_conversion_face];
    for (i=0;i<face_length;i++){
        store_register(A[i*chunksize/bit64or32],A1[i]);
    }
    for (i=0;i<length_of_conversion_face;i++){
        count += (size_t) popcount(A[i]);
    }
    return count;
}


inline size_t get_next_level(void **faces, size_t lenfaces, void **nextfaces, void **nextfaces2, void **forbidden, size_t nr_forbidden, size_t face_length){
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


unsigned int calculate_dimension(void **faces, unsigned int nr_faces, size_t face_length){
    // before doing pretty much anything, we need to know the dimension of the polyhedron
    // this is done by calculating the dimension of an arbitrary facet
    // this dimension is again calculated by calculating the dimension
    // of an arbitrary facet and so on
    unsigned int i,j,k, newfacescounter, dim;
    if (nr_faces == 0){
        return 0;//this isn't supposed to happen, but maybe the data is malformed
    }
    unsigned int bitcount = CountFaceBits(faces[0], face_length);
    if (bitcount == 0){//if a polyhedron contains only the empty face as facet, then it is of dimension 0
        return 0;
    }
    if (nr_faces == 1){//if there is only one facet and this contains bitcount # of vertices/rays/lines then the polyhedron is of dimension bitcount -1
        return bitcount;
    }
    if (bitcount == 1){
        return 1;
    }
    void *nextfaces_creator[nr_faces-1];
    void  *nextfaces2[nr_faces-1], *nextfaces[nr_faces-1];
    for (i=0; i < (nr_faces-1); i++){
        nextfaces_creator[i] = aligned_alloc(chunksize/8,face_length*chunksize/8);
        nextfaces[i] = nextfaces_creator[i];
    }
    newfacescounter = get_next_level(faces,nr_faces,nextfaces,nextfaces2,nextfaces,0,face_length);//calculates the ridges contained in one facet
    dim =  calculate_dimension(nextfaces2,newfacescounter, face_length) + 1;//calculates the dimension of that facet
    if (dim == 1){
        dim = bitcount;//our face should be a somewhat a vertex, but if the polyhedron is unbounded, than our face will have dimension equal to the number of 'vertices' it contains, where some of the vertices might represent lines
    }
    for (i=0; i < (nr_faces - 1); i++){
        free_aligned(nextfaces_creator[i]);
    }
    return dim;
}

void char_from_incidence_list(unsigned int *incidence_list, size_t nr_vertices, \
                              chunktype *array1, size_t face_length){
    size_t entry, position, value, i;
    const size_t size_array = face_length*chunksize/64;
    uint64_t *array = new uint64_t [size_array]();
    while (nr_vertices--) {
        entry = incidence_list[nr_vertices];
        if (entry){
            value = nr_vertices % 64;
            position = nr_vertices/64;
            array[position] += vertex_to_bit_dictionary[value];
        }
    }
    for (i=0;i<face_length;i++){
        load_register(array1[i],array[i*chunksize/64]);
    }
    delete[] array;
}

void char_from_array(unsigned int* input, unsigned int len, \
                     chunktype *array1, size_t face_length){
    size_t entry, position, value, i;
    const size_t size_array = face_length*chunksize/64;
    uint64_t *array = new uint64_t [size_array]();
    while (len--) {
        entry = input[len];
        value = entry % 64;
        position = entry/64;
        array[position] += vertex_to_bit_dictionary[value];
    }
    for (i=0;i<face_length;i++){
        load_register(array1[i],array[i*chunksize/64]);
    }
    delete[] array;
}

void get_facets_from_incidence_matrix(unsigned int **incidence_matrix, void **facets, \
                                      size_t nr_vertices, size_t nr_facets){
    build_dictionary();
    size_t i;
    size_t length_of_face = ((nr_vertices - 1)/chunksize + 1);
    //this determines the length of the face in terms of chunktype
    for(i = 0; i<nr_facets; i++){
        char_from_incidence_list(incidence_matrix[i], nr_vertices, \
                                (chunktype *) facets[i], length_of_face);
    }
}


void get_vertices_from_incidence_matrix(unsigned int **incidence_matrix, void **vertices, \
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
                        (chunktype *) vertices[i], face_length);
    }
}


void get_facets_bitrep_from_facets_pointer( \
        unsigned int ** facets_input, unsigned int *len_facets, \
        void ** facets_output, size_t nr_vertices, size_t nr_facets){
    size_t i;
    size_t face_length = ((nr_vertices - 1)/chunksize + 1);
    //this determines the length of the face in terms of chunktype
    for(i = 0; i<nr_facets; i++){
        char_from_array(facets_input[i], len_facets[i], \
                        (chunktype *) facets_output[i], face_length);
    }
}

void get_vertices_bitrep_from_facets_pointer( \
        unsigned int ** facets_input, unsigned int *len_facets, \
        void ** vertices_output, size_t nr_vertices, size_t nr_facets){
    size_t i,j;
    size_t face_length = ((nr_facets - 1)/chunksize + 1);
    //this determines the length in facet representation in terms of chunktype
    const size_t size_three = nr_facets;
    unsigned int *old_facets_walker = new unsigned int [size_three]();
    unsigned int new_facets_array[size_three];
    unsigned int length_that_face;
    for(i = 0;i<nr_vertices;i++){
        length_that_face = 0;
        for (j=0;j < nr_facets;j++){
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
                        (chunktype *) vertices_output[i], face_length);
    }
    delete[] old_facets_walker;
}
