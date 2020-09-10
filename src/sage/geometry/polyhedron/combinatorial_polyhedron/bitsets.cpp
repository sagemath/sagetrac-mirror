/*
#############################################################################
# Many parts taken from ``src/sage/data_structures/bitset.pxi``
#############################################################################
*/

#include "bitsets.h"
using namespace std;

/*
#############################################################################
# Bitset Initalization
#############################################################################
*/

void bitset_clear(uint64_t* bits, size_t face_length){
    /*
    Remove all elements from the set.
    */
    memset(bits, 0, face_length*8);
}

inline void bitset_copy(uint64_t* dst, uint64_t* src, size_t face_length){
    /*
    Copy the bitset src over to the bitset dst, overwriting dst.

    We assume ``dst.limbs == src.limbs``.
    */
    memcpy(dst, src, face_length*8);
}

void bitset_copy(uint64_t* dst, uint64_t* src, size_t face_length_dst, size_t face_length_src){
    /*
    Copy the bitset src over to the bitset dst, overwriting dst.

    If ``dst`` is longer, then additional bits are set to zero.

    If ``dst`` is shorter, then only the part that fits will be copied.
    */
    if (face_length_src > face_length_dst)
        face_length_src = face_length_dst;

    bitset_copy(dst, src, face_length_src);
    bitset_clear(dst+face_length_src, face_length_dst-face_length_src);
}

/*
#############################################################################
# Bitset Bit Manipulation
#############################################################################
*/

void bitset_set_first_n(uint64_t* bits, size_t face_length, size_t n){
    /*
    Set exactly the first n bits.
    */
    size_t i;
    size_t index = n >> index_shift;
    for(i = 0; i < index; i++)
        bits[i] = -1;
    if (index < face_length)
        bits[index] = limb_lower_bits_down(n);
    for(i=index+1; i < face_length; i++)
        bits[i] = 0;
}

/*
#############################################################################
# Bitset Searching
#############################################################################
*/

inline size_t _bitset_first_in_limb(uint64_t limb){
    /*
    Given a limb of a bitset, return the index of the first nonzero
    bit. If there are no bits set in the limb, return -1.
    */
    if (limb == 0)
        return -1;
    return mpn_scan1(&limb, 0);
}

inline long _bitset_first_in_limb_nonzero(uint64_t limb){
    /*
    Given a non-zero limb of a bitset, return the index of the first
    nonzero bit.
    */
    return mpn_scan1(&limb, 0);
}

size_t bitset_next(uint64_t* bits, size_t face_length, size_t n){
    /*
    Calculate the index of the next element in the set, starting at
    (and including) n.  Return -1 if there are no elements from n
    onwards.
    */
    if (n >= face_length*LIMB_BITS)
        return -1;

    size_t i = n >> index_shift;
    size_t limb = bits[i] & ~limb_lower_bits_down(n);
    size_t ret = _bitset_first_in_limb(limb);
    if (ret != (size_t) -1)
        return (i << index_shift) | ret;
    for(i++; i < face_length; i++){
        if (bits[i])
            return (i << index_shift) | _bitset_first_in_limb_nonzero(bits[i]);
    }
    return -1;
}

/*
#############################################################################
# Miscellaneous
#############################################################################
*/

int test_alignment(uint64_t* a){
    size_t address = (size_t) a;
    size_t required_alignment = chunksize/8;
    return (address == (address & ~(required_alignment -1)));
}

/*
#############################################################################
# Bitset Comparison
#############################################################################
*/

int bitset_cmp(uint64_t* a, uint64_t* b, size_t face_length){
    /*
    Compare bitsets a and b.  Return 0 if the two sets are
    identical, and consistently return -1 or 1 for two sets that are
    not equal.

    We assume ``a.limbs >= b.limbs``.
    */
    return memcmp(a, b, face_length*8);
}

int bitset_isempty(uint64_t* a, size_t face_length){
    /*
    Test whether bits is empty.  Return True (i.e., 1) if the set is
    empty, False (i.e., 0) otherwise.
    */
    if (a[0])
        return 0;
    if (face_length == 1)
        return 1;
    // Compare bits to itself shifted by 1 limb. If these compare equal,
    // all limbs must be 0.
    return (bitset_cmp(a, a+1, face_length -1) == 0);
}
