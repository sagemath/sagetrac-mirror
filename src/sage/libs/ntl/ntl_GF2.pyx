#*****************************************************************************
#       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"
include "sage/ext/cdefs.pxi"
include 'misc.pxi'
include 'decl.pxi'

from sage.rings.integer cimport Integer
from sage.rings.integer_ring cimport IntegerRing_class

##############################################################################
# GF2: Bits
##############################################################################

cdef class ntl_GF2:
    r"""
    The \class{GF2} represents the field GF(2). Computationally
    speaking, it is not a particularly useful class.  Its main use is
    to make the interfaces to the various finite field classes as
    uniform as possible.
    """
    def __init__(self, v=None):
        r"""
        Initializes a NTL bit.

        EXAMPLES:
            sage: ntl.GF2(1)
            1
            sage: ntl.GF2(int(2))
            0
            sage: ntl.GF2('1')
            1
        """
        if isinstance(v, ntl_GF2):
            self.x = (<ntl_GF2>v).x
        elif PyInt_Check(v) or PyLong_Check(v) or isinstance(v, Integer):
            GF2_conv_long(self.x, int(v) % 2)
        elif v is not None:
            v = str(v)
            sig_on()
            GF2_from_str(&self.x, v)
            sig_off()

    def __cinit__(self):
        GF2_construct(&self.x)

    def __dealloc__(self):
        GF2_destruct(&self.x)

    def __repr__(self):
        """
        Return the string representation of self.

        EXAMPLES:
            sage: str(ntl.GF2(1)) # indirect doctest
            '1'
        """
        return GF2_to_PyString(&self.x)

    def __reduce__(self):
        """
        Serializes self.

        EXAMPLES:
            sage: a = ntl.GF2(1)
            sage: loads(dumps(a))
            1
        """
        return unpickle_class_value, (ntl_GF2, int(self))

    def __richcmp__(self, other, op):
        """
        Compare self to other.

        EXAMPLES:
            sage: a = ntl.GF2(1)
            sage: b = ntl.GF2(0)
            sage: a == b
            False
        """
        if op != 2 and op != 3:
            raise TypeError, "elements in GF(2) are not ordered."

        if not isinstance(other, ntl_GF2):
            other = ntl_GF2(other)

        if not isinstance(self, ntl_GF2):
            self = ntl_GF2(self)

        cdef int t
        t = GF2_equal((<ntl_GF2>self).x, (<ntl_GF2>other).x)
        if op == 2:
            return t == 1
        elif op == 3:
            return t == 0

    def __mul__(self, other):
        """
            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: o*o
            1
            sage: o*z
            0
            sage: z*o
            0
            sage: z*z
            0
        """
        cdef ntl_GF2 r = ntl_GF2.__new__(ntl_GF2)
        if not isinstance(self, ntl_GF2):
            self = ntl_GF2(self)
        if not isinstance(other, ntl_GF2):
            other = ntl_GF2(other)
        GF2_mul(r.x, (<ntl_GF2>self).x, (<ntl_GF2>other).x)
        return r

    def __div__(self, other):
        """
            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: o/o
            1
            sage: o/z
            Traceback (most recent call last):
            ...
            ZeroDivisionError
        """
        cdef ntl_GF2 r
        if not isinstance(self, ntl_GF2):
            self = ntl_GF2(self)
        if not isinstance(other, ntl_GF2):
            other = ntl_GF2(other)
        if GF2_IsZero((<ntl_GF2>other).x):
            raise ZeroDivisionError
        r = ntl_GF2.__new__(ntl_GF2)
        GF2_div(r.x, (<ntl_GF2>self).x, (<ntl_GF2>other).x)
        return r

    def __sub__(self, other):
        """
            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: o-o
            0
            sage: o-z
            1
            sage: z-o
            1
            sage: z-z
            0
        """
        cdef ntl_GF2 r = ntl_GF2.__new__(ntl_GF2)
        if not isinstance(self, ntl_GF2):
            self = ntl_GF2(self)
        if not isinstance(other, ntl_GF2):
            other = ntl_GF2(other)
        GF2_sub(r.x, (<ntl_GF2>self).x, (<ntl_GF2>other).x)
        return r

    def __add__(self, other):
        """
            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: o+o
            0
            sage: o+z
            1
            sage: z+o
            1
            sage: z+z
            0
        """
        cdef ntl_GF2 r = ntl_GF2.__new__(ntl_GF2)
        if not isinstance(self, ntl_GF2):
            self = ntl_GF2(self)
        if not isinstance(other, ntl_GF2):
            other = ntl_GF2(other)
        GF2_add(r.x, (<ntl_GF2>self).x, (<ntl_GF2>other).x)
        return r

    def __neg__(ntl_GF2 self):
        """
            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: -z
            0
            sage: -o
            1
        """
        cdef ntl_GF2 r = ntl_GF2.__new__(ntl_GF2)
        GF2_negate(r.x, self.x)
        return r

    def __pow__(ntl_GF2 self, long e, ignored):
        """
            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: z^2
            0
            sage: o^2
            1
        """
        cdef ntl_GF2 r = ntl_GF2()
        GF2_power(r.x, self.x, e)
        return r

    def __int__(self):
        """
        Return self as an int.

        EXAMPLES:
            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: int(z)
            0
            sage: int(o)
            1
        """
        cdef long l = GF2_conv_to_long(self.x)
        return int(l)

def unpickle_class_value(cls, x):
    """
    Here for unpickling.

    EXAMPLES:
        sage: sage.libs.ntl.ntl_GF2.unpickle_class_value(ntl.GF2,1)
        1
        sage: type(sage.libs.ntl.ntl_GF2.unpickle_class_value(ntl.GF2,1))
        <type 'sage.libs.ntl.ntl_GF2.ntl_GF2'>
    """
    return cls(x)

def unpickle_class_args(cls, x):
    """
    Here for unpickling.

    EXAMPLES:
        sage: sage.libs.ntl.ntl_GF2.unpickle_class_args(ntl.GF2,[1])
        1
        sage: type(sage.libs.ntl.ntl_GF2.unpickle_class_args(ntl.GF2,[1]))
        <type 'sage.libs.ntl.ntl_GF2.ntl_GF2'>
    """
    return cls(*x)

