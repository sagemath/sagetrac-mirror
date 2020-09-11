#ifndef COMBINATORIAL_BITSETS
#define COMBINATORIAL_BITSETS

#include <cstring>
#include "gmp.h"
#include <cstdint>
#include <cstdio>

// Instead of division by 64.
const size_t index_shift = 6;
// Size of ``uint64_t`` in our case.
const size_t LIMB_BITS = 64;

// Any Bit-representation is assumed to be `chunksize`-Bit aligned.
const size_t chunksize = 64;

// Number of ``uint64_t`` a primary operation can proceed.
const size_t intersection_chunk = 1;
const size_t subset_chunk = 1;
const size_t union_chunk = 1;

inline size_t get_face_length(size_t n_atoms){
    // Obtain the number of ``uint64_t`` needed to fit
    // that many atoms.
    return ((n_atoms - 1)/chunksize + 1) * chunksize/LIMB_BITS;
}

/*
#############################################################################
# Basic functions that can be optimized by intrinsics
#############################################################################
*/

inline int is_subset(uint64_t *A, uint64_t *B){
    return !(A[0] & ~B[0]);
}

inline void intersection(uint64_t* dest, uint64_t *A, uint64_t *B){
    dest[0] = A[0] & B[0];
}

inline void unite(uint64_t* dest, uint64_t *A, uint64_t *B){
    dest[0] = A[0] | B[0];
}

inline size_t count_atoms(uint64_t a){
    a = a - ((a >> 1) & 0x5555555555555555ULL);
    a = (a & 0x3333333333333333ULL) + ((a >> 2) & 0x3333333333333333ULL);
    return ( ((a + (a >> 4)) & 0x0f0f0f0f0f0f0f0fULL) * 0x0101010101010101ULL ) >> 56;
}

/*
#############################################################################
# Creating limb patterns
#############################################################################
*/

inline uint64_t limb_one_set_bit(size_t n){
    /*
    Return a limb with only bit n set.
    */
    return ((uint64_t) 1) << (n % LIMB_BITS);
}

inline uint64_t limb_one_zero_bit(size_t n){
    /*
    Return a limb with all bits set, except for bit n.
    */
    return ~limb_one_set_bit(n);
}

inline uint64_t limb_lower_bits_down(size_t n){
    /*
    Return a limb with the lower n bits set, where n is interpreted
    in [0 .. 64-1].
    */
    return (((uint64_t) 1) << (n % LIMB_BITS)) - 1;
}

/*
#############################################################################
# Bitset Initalization
#############################################################################
*/

void bitset_clear(uint64_t* bits, size_t face_length);

void bitset_copy(uint64_t* dst, uint64_t* src, size_t face_length_dst, size_t face_length_src);

/*
#############################################################################
# Bitset Bit Manipulation
#############################################################################
*/

inline int bitset_in(uint64_t* bits, size_t n){
    /*
    Check if n is in bits.  Return True (i.e., 1) if n is in the
    set, False (i.e., 0) otherwise.
    */
    return (bits[n >> index_shift] & limb_one_set_bit(n)) > 0;
}

inline void bitset_discard(uint64_t* bits, size_t n){
    /*
    Remove n from bits.
    */
    bits[n >> index_shift] &= limb_one_zero_bit(n);
}

inline void bitset_add(uint64_t* bits, size_t n){
    /*
    Add n to bits.
    */
    bits[n >> index_shift] |= limb_one_set_bit(n);
}

void bitset_set_first_n(uint64_t* bits, size_t face_length, size_t n);

/*
#############################################################################
# Bitset Searching
#############################################################################
*/

size_t bitset_next(uint64_t* bits, size_t face_length, size_t n);

inline size_t count_atoms(const uint64_t* A, size_t face_length){
    /*
    Return the number of atoms/vertices in A.
    This is the number of set bits in A.
    ``face_length`` is the length of A in terms of uint64_t.
    */
    size_t i;
    unsigned int count = 0;
    for (i=0; i<face_length; i++){
        count += count_atoms(A[i]);
    }
    return count;
}

/*
#############################################################################
# Miscellaneous
#############################################################################
*/

int test_alignment(uint64_t* a);

/*
#############################################################################
# Bitset Comparison
#############################################################################
*/

int bitset_cmp(uint64_t* a, uint64_t* b, size_t face_length);
    /*
    Compare bitsets a and b.  Return 0 if the two sets are
    identical, and consistently return -1 or 1 for two sets that are
    not equal.
    */

int bitset_isempty(uint64_t* a, size_t face_length);

inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
    /*
    Return ``A & ~B == 0``.
    A is not subset of B, iff there is a vertex in A, which is not in B.
    ``face_length`` is the length of A and B in terms of uint64_t.
    */
    for (size_t i = 0; i < face_length; i+=subset_chunk){
        if (!is_subset(A+i, B+i))
            return 0;
    }
    return 1;
}

/*
#############################################################################
# Bitset Arithmetic
#############################################################################
*/

inline void intersection(uint64_t *dest, uint64_t *A, uint64_t *B, \
                         size_t face_length){
    /*
    Set ``dest = A & B``, i.e. dest is the intersection of A and B.
    ``face_length`` is the length of A, B and dest in terms of uint64_t.
    */
    for (size_t i = 0; i < face_length; i+=intersection_chunk){
        intersection(dest+i, A+i, B+i);
    }
}

inline void unite(uint64_t *dest, uint64_t *A, uint64_t *B, \
                         size_t face_length){
    /*
    Set ``dest = A | B``, i.e. dest is the union of A and B.
    ``face_length`` is the length of A, B and dest in terms of uint64_t.
    */
    for (size_t i = 0; i < face_length; i+=union_chunk){
        unite(dest+i, A+i, B+i);
    }
}

#endif
