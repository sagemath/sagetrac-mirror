#*****************************************************************************
#     Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .bitset_base cimport bitset_t

# Python layer over bitset_t
cdef class FrozenBitset:
    cdef bitset_t _bitset
    cdef FrozenBitset _new(self,long int capacity)
    cdef FrozenBitset _larger_capacity_(self, long size)
    cdef long capacity(self)
    cdef bint isempty(self)
    cdef bint issubset(self, FrozenBitset other) except -1
    cdef bint issuperset(self, FrozenBitset other) except -1
    cdef bint isdisjoint(self, FrozenBitset other) except -1
    cdef _union(self, FrozenBitset other)
    cdef intersection(self, FrozenBitset other)
    cdef difference(self, FrozenBitset other)
    cdef symmetric_difference(self, FrozenBitset other)
    cdef complement(self)
    cdef __copy__(self)

cdef class Bitset(FrozenBitset):
    cdef __copy__(self)
    cdef update(self, FrozenBitset other)
    cdef intersection_update(self, FrozenBitset other)
    cdef difference_update(self, FrozenBitset other)
    cdef symmetric_difference_update(self, FrozenBitset other)
    cdef add(self, unsigned long n)
    cdef remove(self, unsigned long n)
    cdef discard(self, unsigned long n)
    cdef pop(self)
    cdef clear(self)

