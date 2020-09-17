"""
A fast bitset datatype in Cython

Operations between bitsets are only guaranteed to work if the bitsets
have the same size, with the exception of ``bitset_realloc``.  Similarly, you
should not try to access elements of a bitset beyond the size.

AUTHORS:

- Robert Bradshaw (2008)
- Rudi Pendavingh, Stefan van Zwam (2013-06-06): added functions map, lex_cmp,
  pickle, unpickle
- Jeroen Demeyer (2014-09-05): use mpn_* functions from MPIR in the
  implementation (:trac:`13352` and :trac:`16937`)
- Simon King (2014-10-28): ``bitset_rshift`` and ``bitset_lshift`` respecting
  the size of the given bitsets (:trac:`15820`)
"""

#*****************************************************************************
#     Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.string cimport strlen
from cysignals.memory cimport check_calloc, check_reallocarray, sig_malloc, sig_free

from sage.cpython.string cimport char_to_str, str_to_bytes, bytes_to_str
from sage.libs.gmp.mpn cimport *
from sage.data_structures.alternative_bitset cimport *
from cython.operator import preincrement as preinc
from libc.stdint cimport UINT32_MAX

from sage.libs.croaring cimport *

# Doctests for the functions in this file are in sage/data_structures/bitset.pyx

cdef inline mp_limb_t limb_lower_bits_up(mp_bitcnt_t n):
    """
    Return a limb with the lower n bits set, where n is interpreted
    in [1 .. GMP_LIMB_BITS].
    """
    return (<mp_limb_t>(-1)) >> ((<unsigned int>(-n)) % 64)


#############################################################################
# Bitset Initalization
#############################################################################
cdef inline bint bitset_init(bitset_t bits, mp_bitcnt_t size) except -1:
    """
    Allocate an empty bitset of size ``size``.

    Size must be at least 1.
    """
    bits.bits = roaring_bitmap_create()
    bits.size = size

cdef inline int bitset_realloc(bitset_t bits, mp_bitcnt_t size) except -1:
    """
    Reallocate a bitset to size ``size``. If reallocation is larger, new bitset
    does not contain any of the extra bits.
    """
    bits.size = size
    bitset_fix(bits)

cdef inline void bitset_free(bitset_t bits):
    """
    Deallocate the memory in bits.
    """
    roaring_bitmap_free(bits.bits)

cdef inline void bitset_clear(bitset_t bits):
    """
    Remove all elements from the set.
    """
    roaring_bitmap_clear(bits.bits)

cdef inline void bitset_zero(bitset_t bits):
    """
    Remove all elements from the set.

    This function is the same as bitset_clear(bits).
    """
    bitset_clear(bits)

cdef int bitset_copy(bitset_t dst, bitset_t src) except 0:
    """
    Copy the bitset src over to the bitset dst, overwriting dst.

    We assume ``dst.limbs == src.limbs``.
    """
    dst.size = src.size
    return roaring_bitmap_overwrite(dst.bits, src.bits)

cdef inline void bitset_fix(bitset_t bits):
    """
    Clear upper bits which should be zero according to the size.
    """
    return
    cdef long maxi = bitset_last(bits)
    if maxi >= bits.size:
        roaring_bitmap_remove_range(bits.bits, bits.size, maxi+1)

#############################################################################
# Bitset Comparison
#############################################################################

cdef inline bint bitset_equal_bits(bitset_t a, bitset_t b, mp_bitcnt_t n):
    """
    Return ``True`` iff the first n bits of a and b agree.
    """
    bitset_fix(a)
    bitset_fix(b)
    cdef bitset_t c
    bitset_init(c, a.size)
    cdef bitset_t d
    bitset_init(d, a.size)
    bitset_copy(c, a)
    bitset_copy(d, b)
    roaring_bitmap_remove_range(c.bits, n, a.size)
    roaring_bitmap_remove_range(d.bits, n, b.size)
    cdef bint output = roaring_bitmap_equals(c.bits, d.bits)
    bitset_free(c)
    bitset_free(d)
    return output

cdef inline bint bitset_equal_bits_shifted(bitset_t a, bitset_t b, mp_bitcnt_t n, mp_bitcnt_t offset_a, mp_bitcnt_t offset_b):
    """
    Return ``True`` iff the first n bits of *b1 and the bits ranging from
    offset to offset+n of *b2 agree.
    """
    cdef bitset_iterator_t it_a = bitset_iterator_init(a)
    cdef bitset_iterator_t it_b = bitset_iterator_init(b)

    bitset_iterator_move(it_a, offset_a)
    bitset_iterator_move(it_b, offset_b)

    cdef long val_a, val_b, max_val_a
    max_val_a = offset_a + n

    val_a = bitset_iterator_value(it_a)
    val_b = bitset_iterator_value(it_b)

    cdef bint output = True

    while -1 < val_a < max_val_a:
        if not val_a == val_b:
            output = False
            break

    bitset_iterator_free(it_a)
    bitset_iterator_free(it_b)
    return output

cdef inline bint bitset_isempty(bitset_t bits):
    """
    Test whether bits is empty.  Return True (i.e., 1) if the set is
    empty, False (i.e., 0) otherwise.
    """
    return roaring_bitmap_is_empty(bits.bits)

cdef inline bint bitset_is_zero(bitset_t bits):
    """
    Test whether bits is empty (i.e., zero).  Return True (1) if
    the set is empty, False (0) otherwise.

    This function is the same as bitset_is_empty(bits).
    """
    return bitset_isempty(bits)

cdef inline bint bitset_eq(bitset_t a, bitset_t b):
    """
    Compare bitset a and b.  Return True (i.e., 1) if the sets are
    equal, and False (i.e., 0) otherwise.

    We assume ``a.limbs >= b.limbs``.
    """
    return roaring_bitmap_equals(a.bits, b.bits)

cdef inline int bitset_cmp(bitset_t a, bitset_t b):
    """
    Compare bitsets a and b.  Return 0 if the two sets are
    identical, and consistently return -1 or 1 for two sets that are
    not equal.

    We assume ``a.limbs >= b.limbs``.
    """
    return bitset_lex_cmp(a, b)

cdef inline int bitset_lex_cmp(bitset_t a, bitset_t b):
    """
    Compare bitsets ``a`` and ``b`` using lexicographical ordering.

    In this order `a < b` if, for some `k`, the first `k` elements from
    `[0 ... n-1]` are in `a` if and only if they are in `b`, and the
    `(k+1)`st element is in `b` but not ``a``. So `1010 < 1011` and
    `1010111 < 1011000`.

    We assume ``a.limbs == b.limbs``.

    INPUT:

    - ``a`` -- a bitset
    - ``b`` -- a bitset, assumed to have the same size as ``a``.

    OUTPUT:

    Return ``0`` if the two sets are identical, return ``1`` if ``a > b``,
    and return ``-1`` if ``a < b``.
    """
    cdef long i = bitset_first_diff(a, b)
    if i == -1:
        return 0
    if bitset_in(a, i):
        return 1
    else:
        return -1

cdef inline bint bitset_issubset(bitset_t a, bitset_t b):
    """
    Test whether a is a subset of b (i.e., every element in a is also
    in b).

    We assume ``a.limbs <= b.limbs``.
    """
    return roaring_bitmap_is_subset(a.bits, b.bits)

cdef inline bint bitset_issuperset(bitset_t a, bitset_t b):
    """
    Test whether a is a superset of b (i.e., every element in b is also
    in a).

    We assume ``a.limbs >= b.limbs``.
    """
    return bitset_issubset(b, a)


cdef inline bint bitset_are_disjoint(bitset_t a, bitset_t b):
    """
    Tests whether ``a`` and ``b`` have an empty intersection.

    We assume ``a.limbs <= b.limbs``.
    """
    return not roaring_bitmap_intersect(a.bits, b.bits)


#############################################################################
# Bitset Bit Manipulation
#############################################################################

cdef inline bint bitset_in(bitset_t bits, mp_bitcnt_t n):
    """
    Check if n is in bits.  Return True (i.e., 1) if n is in the
    set, False (i.e., 0) otherwise.
    """
    return roaring_bitmap_contains(bits.bits, n)

cdef inline bint bitset_check(bitset_t bits, mp_bitcnt_t n):
    """
    Check if n is in bits.  Return True (i.e., 1) if n is in the
    set, False (i.e., 0) otherwise.

    This function is the same as bitset_in(bits, n).
    """
    return bitset_in(bits, n)

cdef inline bint bitset_not_in(bitset_t bits, mp_bitcnt_t n):
    """
    Check if n is not in bits.  Return True (i.e., 1) if n is not in the
    set, False (i.e., 0) otherwise.
    """
    return not bitset_in(bits, n)

cdef inline bint bitset_remove(bitset_t bits, mp_bitcnt_t n) except -1:
    """
    Remove n from bits.  Raise KeyError if n is not contained in bits.
    """
    if not roaring_bitmap_remove_checked(bits.bits, n):
        raise KeyError(n)

cdef inline void bitset_discard(bitset_t bits, mp_bitcnt_t n):
    """
    Remove n from bits.
    """
    roaring_bitmap_remove_checked(bits.bits, n)

cdef inline void bitset_unset(bitset_t bits, mp_bitcnt_t n):
    """
    Remove n from bits.

    This function is the same as bitset_discard(bits, n).
    """
    bitset_discard(bits, n)


cdef inline void bitset_add(bitset_t bits, mp_bitcnt_t n):
    """
    Add n to bits.
    """
    roaring_bitmap_add(bits.bits, n)

cdef inline void bitset_set(bitset_t bits, mp_bitcnt_t n):
    """
    Add n to bits.

    This function is the same as bitset_add(bits, n).
    """
    bitset_add(bits, n)

cdef inline void bitset_set_to(bitset_t bits, mp_bitcnt_t n, bint b):
    """
    If b is True, add n to bits.  If b is False, remove n from bits.
    """
    if b:
        bitset_add(bits, n)
    else:
        bitset_discard(bits, n)

cdef inline void bitset_flip(bitset_t bits, mp_bitcnt_t n):
    """
    If n is in bits, remove n from bits.  If n is not in bits, add n
    to bits.
    """
    if not roaring_bitmap_remove_checked(bits.bits, n):
        bitset_add(bits, n)

cdef inline void bitset_set_first_n(bitset_t bits, mp_bitcnt_t n):
    """
    Set exactly the first n bits.
    """
    roaring_bitmap_add_range(bits.bits, 0, n)

#############################################################################
# Bitset Searching
#############################################################################

cdef inline long bitset_first(bitset_t a):
    """
    Calculate the index of the first element in the set. If the set
    is empty, returns -1.
    """
    cdef uint32_t mini = roaring_bitmap_minimum(a.bits)
    if mini == UINT32_MAX:
        return -1
    else:
        return mini

cdef inline long bitset_last(bitset_t a):
    """
    Calculate the index of the last element in the set. If the set
    is empty, returns -1.
    """
    cdef uint32_t maxi = roaring_bitmap_maximum(a.bits)
    if maxi == 0:
        return bitset_first(a)
    else:
        return maxi

cdef inline long bitset_first_in_complement(bitset_t a):
    """
    Calculate the index of the first element not in the set.

    Return -1, if it does not exist.
    """
    cdef bitset_t neg
    cdef long first, last, output
    first = bitset_first(a)
    if first > 0:
        return 0
    last = bitset_last(a)
    neg.bits = roaring_bitmap_flip(a.bits, 0, last)
    output = bitset_first(neg)
    bitset_free(neg)
    if output == a.size:
        return -1
    return output

cdef inline long bitset_pop(bitset_t a) except -1:
    """
    Remove and return an arbitrary element from the set. Raise
    KeyError if the set is empty.
    """
    cdef long i = bitset_first(a)
    if i == -1:
        raise KeyError('pop from an empty set')
    bitset_discard(a, i)
    return i

cdef inline long bitset_first_diff(bitset_t a, bitset_t b):
    """
    Calculate the index of the first difference between a and b.  If a
    and b are equal, then return -1.

    We assume ``a.limbs == b.limbs``.
    """
    cdef roaring_bitmap_t* c = roaring_bitmap_xor(a.bits, b.bits)
    cdef mp_bitcnt_t first_diff = roaring_bitmap_minimum(c)
    roaring_bitmap_free(c)
    if first_diff == UINT32_MAX:
        return -1
    else:
        return first_diff

#############################################################################
# Bitset Iterator
#############################################################################

cdef bitset_iterator_t bitset_iterator_init(bitset_t bits):
    """
    Allocate an iterator of self.
    """
    return roaring_create_iterator(bits.bits)

cdef void bitset_iterator_free(bitset_iterator_t it):
    roaring_free_uint32_iterator(it)

cdef void bitset_iterator_move(bitset_iterator_t it, uint32_t val):
    roaring_move_uint32_iterator_equalorlarger(it, val)

cdef inline long bitset_iterator_value(bitset_iterator_t it):
    if it[0].has_value:
        return it[0].current_value
    else:
        -1

cdef inline long bitset_iterator_next(bitset_iterator_t it):
    roaring_advance_uint32_iterator(it)
    return bitset_iterator_value(it)

cdef inline long bitset_next(bitset_t a, mp_bitcnt_t n):
    """
    Calculate the index of the next element in the set, starting at
    (and including) n.  Return -1 if there are no elements from n
    onwards.
    """
    cdef bitset_iterator_t it = bitset_iterator_init(a)
    bitset_iterator_move(it, n)
    cdef long value = bitset_iterator_value(it)
    bitset_iterator_free(it)
    return value

cdef inline long bitset_next_diff(bitset_t a, bitset_t b, mp_bitcnt_t n):
    """
    Calculate the index of the next element that differs between a and
    b, starting at (and including) n.  Return -1 if there are no
    elements differing between a and b from n onwards.

    We assume ``a.limbs == b.limbs``.
    """
    cdef bitset_t diff
    diff.bits = roaring_bitmap_xor(a.bits, b.bits)
    cdef long val = bitset_next(diff, n)
    bitset_free(diff)
    return val

cdef inline long bitset_len(bitset_t bits):
    """
    Calculate the number of items in the set (i.e., the number of nonzero bits).
    """
    return roaring_bitmap_get_cardinality(bits.bits)

cdef inline long bitset_hash(bitset_t bits):
    """
    Calculate a (very naive) hash function.

    This function should not depend on the size of the bitset, only on
    the items in the bitset.
    """
    # Stolen from PyRoaringBitMap
    cdef int64_t h_val = 0
    cdef uint32_t i, count, max_count=256
    cdef bitset_iterator_t iterator  = roaring_create_iterator(bits.bits)
    cdef uint32_t *buff = <uint32_t*>sig_malloc(max_count*4)
    while True:
        count = roaring_read_uint32_iterator(iterator, buff, max_count)
        i = 0
        while i < count:
            h_val = ((h_val << 2) + buff[i])
            # TODO find a good hash formula
            # This one should be better, but is too long:
            # h_val = ((h_val<<16) + buff[i]) % 1748104473534059
            i += 1
        if count != max_count:
            break
    bitset_iterator_free(iterator)
    sig_free(buff)
    if bitset_first(bits) == -1:
        return -1
    return h_val

cdef inline mp_limb_t bitset_get_limb(bitset_t bits, mp_bitcnt_t n):
    r"""
    Get a uint64_t representing the bits 64*n through 64*n + 63.
    """
    cdef bitset_iterator_t it = bitset_iterator_init(bits)
    bitset_iterator_move(it, n*64)
    cdef long val = bitset_iterator_value(it)
    while -1 < val < (n+1)*64:
        output |= (<mp_limb_t> 1 << (val - 64*n))

    return output

cdef inline void bitset_and_limb(bitset_t bits, mp_limb_t limb, mp_bitcnt_t n):
    cdef uint32_t i
    for i in range(64):
        if not (limb & (<mp_limb_t> 1 << (i))):
            bitset_discard(bits, n*64+i)

cdef inline void bitset_or_limb(bitset_t bits, mp_limb_t limb, mp_bitcnt_t n):
    cdef uint32_t i
    for i in range(64):
        if (limb & (<mp_limb_t> 1 << (i))):
            bitset_add(bits, n*64+i)

cdef inline void bitset_xor_limb(bitset_t bits, mp_limb_t limb, mp_bitcnt_t n):
    cdef uint32_t i
    for i in range(64):
        if (limb & (<mp_limb_t> 1 << (i))):
            bitset_flip(bits, n*64+i)

cdef inline void bitset_assign_limb(bitset_t bits, mp_limb_t limb, mp_bitcnt_t n):
    cdef uint32_t i
    for i in range(64):
        if (limb & (<mp_limb_t> 1 << (i))):
            bitset_add(bits, n*64+i)
        else:
            bitset_remove(bits, n*64+i)


#############################################################################
# Bitset Arithmetic
#############################################################################

cdef inline void bitset_complement(bitset_t r, bitset_t a):
    """
    Set r to be the complement of a, overwriting r.

    We assume ``r.limbs == a.limbs``.
    """
    roaring_bitmap_overwrite(r.bits, a.bits)
    roaring_bitmap_flip_inplace(r.bits, 0, a.size)
    r.size = a.size
    bitset_fix(r)

cdef inline void bitset_not(bitset_t r, bitset_t a):
    """
    Set r to be the complement of a, overwriting r.

    We assume ``r.limbs == a.limbs``.

    This function is the same as bitset_complement(r, a).
    """
    bitset_complement(r, a)

cdef inline void bitset_intersection(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the intersection of a and b, overwriting r.

    We assume ``a.limbs >= r.limbs == b.limbs``.
    """
    roaring_bitmap_overwrite(r.bits, a.bits)
    roaring_bitmap_and_inplace(r.bits, b.bits)
    r.size = a.size

cdef inline void bitset_and(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the intersection of a and b, overwriting r.

    We assume ``a.limbs >= r.limbs == b.limbs``.

    This function is the same as bitset_intersection(r, a, b).
    """
    bitset_intersection(r, a, b)

cdef inline void bitset_union(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the union of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    roaring_bitmap_overwrite(r.bits, a.bits)
    roaring_bitmap_or_inplace(r.bits, b.bits)
    r.size = b.size
    if a.size > b.size:
        r.size = a.size

cdef inline void bitset_or(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the union of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.

    This function is the same as bitset_union(r, a, b).
    """
    bitset_union(r, a, b)

cdef inline void bitset_difference(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the difference of a and b (i.e., things in a that are not
    in b), overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    roaring_bitmap_overwrite(r.bits, a.bits)
    roaring_bitmap_or_inplace(r.bits, b.bits)
    r.size = b.size
    if a.size > b.size:
        r.size = a.size

cdef inline void bitset_symmetric_difference(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the symmetric difference of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    roaring_bitmap_overwrite(r.bits, a.bits)
    roaring_bitmap_xor_inplace(r.bits, b.bits)
    r.size = a.size

cdef inline void bitset_xor(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the symmetric difference of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.

    This function is the same as bitset_symmetric_difference(r, a, b).
    """
    bitset_symmetric_difference(r, a, b)


cdef void bitset_rshift(bitset_t r, bitset_t a, mp_bitcnt_t n):
    """
    Shift the bitset ``a`` right by ``n`` bits and store the result in
    ``r``.

    There are no assumptions on the sizes of ``a`` and ``r``.  Bits which are
    shifted outside of the resulting bitset are discarded.
    """
    r.size = a.size
    if n >= a.size:
        return

    cdef bitset_t foo
    bitset_init(foo, a.size)

    cdef bitset_iterator_t it = bitset_iterator_init(a)

    cdef long val = bitset_iterator_value(it)
    while val != -1:
        bitset_add(foo, (<mp_bitcnt_t> val) + n)
        val = bitset_iterator_next(it)

    bitset_iterator_free(it)

    bitset_fix(foo)
    bitset_free(r)
    r.bits = foo.bits

cdef void bitset_lshift(bitset_t r, bitset_t a, mp_bitcnt_t n):
    """
    Shift the bitset ``a`` left by ``n`` bits and store the result in
    ``r``.

    There are no assumptions on the sizes of ``a`` and ``r``.  Bits which are
    shifted outside of the resulting bitset are discarded.
    """
    r.size = a.size
    if n >= r.size:
        return

    cdef bitset_t foo
    bitset_init(foo, a.size)

    cdef bitset_iterator_t it = bitset_iterator_init(a)

    bitset_iterator_move(it, n)
    cdef long val = bitset_iterator_value(it)
    while val != -1:
        bitset_add(foo, (<mp_bitcnt_t> val) - n)
        val = bitset_iterator_next(it)

    bitset_iterator_free(it)

    bitset_fix(foo)
    bitset_free(r)
    r.bits = foo.bits

cdef int bitset_map(bitset_t r, bitset_t a, m) except -1:
    """
    Fill bitset ``r`` so ``r == {m[i] for i in a}``.

    We assume ``m`` has a dictionary interface such that
    ``m[i]`` is an integer in ``[0 ... n-1]`` for all ``i`` in ``a``,
    where ``n`` is the capacity of ``r``.
    """
    bitset_clear(r)
    cdef bitset_iterator_t it = bitset_iterator_init(a)

    cdef long val = bitset_iterator_value(it)
    while val != -1:
        bitset_add(r, (<mp_bitcnt_t> m[val]))
        val = bitset_iterator_next(it)

    bitset_iterator_free(it)

#############################################################################
# Hamming Weights
#############################################################################

cdef inline long bitset_hamming_weight(bitset_t a):
    return bitset_len(a)

#############################################################################
# Bitset Conversion
#############################################################################

cdef char* bitset_chars(char* s, bitset_t bits, char zero=c'0', char one=c'1'):
    """
    Return a string representation of the bitset in s, using zero for
    the character representing the items not in the bitset and one for
    the character representing the items in the bitset.

    The string is both stored in s and returned.  If s is NULL, then a
    new string is allocated.
    """
    cdef mp_bitcnt_t i
    if s == NULL:
        s = <char *>sig_malloc(bits.size + 1)
    for i from 0 <= i < bits.size:
        s[i] = one if bitset_in(bits, i) else zero
    s[bits.size] = 0
    return s


cdef int bitset_from_char(bitset_t bits, char* s, char zero=c'0', char one=c'1') except -1:
    """
    Initialize a bitset with a set derived from the C string s, where one
    represents the character indicating set membership.
    """
    bitset_init(bits, strlen(s))
    cdef mp_bitcnt_t i
    for i from 0 <= i < bits.size:
        bitset_set_to(bits, i, s[i] == one)
    return 0


cdef int bitset_from_str(bitset_t bits, object s, char zero=c'0', char one=c'1') except -1:
    """
    Initialize a bitset with a set derived from the Python str s, where one
    represents the character indicating set membership.
    """
    cdef bytes b = str_to_bytes(s)
    return bitset_from_char(bits, b, zero, one)


cdef bitset_string(bitset_t bits):
    """
    Return a python string representing the bitset.
    """
    return bytes_to_str(bitset_bytes(bits))


cdef bitset_bytes(bitset_t bits):
    """
    Return a python bytes string representing the bitset.

    On Python 2 this is equivalent to bitset_string.
    """
    cdef char* s = bitset_chars(NULL, bits)
    cdef object py_s
    py_s = s
    sig_free(s)
    return py_s


cdef list bitset_list(bitset_t bits):
    """
    Return a list of elements in the bitset.
    """
    cdef long length = bitset_len(bits)
    if length == 0:
        return []
    cdef uint32_t* buf = <uint32_t*> sig_malloc(length*sizeof(uint32_t))
    cdef bitset_iterator_t it = bitset_iterator_init(bits)
    roaring_read_uint32_iterator(it, buf, length)

    cdef list elts = [None]*length
    cdef long i
    for i in range(length):
        elts[i] = buf[i]

    sig_free(buf)
    bitset_iterator_free(it)
    return elts

cdef bitset_pickle(bitset_t bs):
    """
    Convert ``bs`` to a reasonably compact Python structure.

    Useful for pickling objects using bitsets as internal data structure.
    To ensure this works on 32-bit and 64-bit machines, the size of a long
    is stored too.
    """
    raise NotImplementedError

cdef bitset_unpickle(bitset_t bs, tuple input):
    """
    Convert the data into a bitset.

    Companion of ``bitset_pickle()``. Assumption: ``bs`` has been initialized.
    """
    raise NotImplementedError
