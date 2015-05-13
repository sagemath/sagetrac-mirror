#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
include 'misc.pxi'
include 'decl.pxi'

from sage.rings.integer import Integer
from sage.rings.integer cimport Integer

zz_pContextDict = {}

cdef class ntl_zz_pContext_class:
    def __init__(self, long v):
        """
        EXAMPLES:
            # You can construct contexts manually.
            sage: c = ntl.zz_pContext(11)
            sage: n1 = ntl.zz_p(12,c)
            sage: n1
            1

            # or You can construct contexts implicitly.
            sage: n2=ntl.zz_p(12, 7)
            sage: n2
            5
            sage: ntl.zz_p(2,3)+ntl.zz_p(1,3)
            0
            sage: n2+n1  # Mismatched moduli:  It will go BOOM!
            Traceback (most recent call last):
            ...
            ValueError: arithmetic operands must have the same modulus.
        """
        pass

    def __cinit__(self, long v):
        if v > NTL_SP_BOUND:
            raise ValueError, "Modulus (=%s) is too big"%v
        zz_pContext_construct_long(&self.x, v)
        zz_pContextDict[repr(v)] = self
        self.p = v

    def __dealloc__(self):
        zz_pContext_destruct(&self.x)

    def __reduce__(self):
        """
        sage: c=ntl.zz_pContext(13)
        sage: loads(dumps(c)) is c
        True
        """
        return ntl_zz_pContext, (self.p,)

    def modulus(self):
        """
        Print the modulus for self.

        EXAMPLES:
            sage: c1 = ntl.zz_pContext(36)
            sage: c1.modulus()
            36
        """
        return self.p

    def restore(self):
        """
        Restore a zz_pContext.

        EXAMPLES:
            sage: c = ntl.zz_pContext(5)
            sage: m = ntl.zz_p(4,7)
            sage: c.restore()
        """
        self.restore_c()

    cdef void restore_c(self):
        """
        Actual code for the above.

        EXAMPLES:
            sage: n = ntl.zz_p(3,5)
            sage: m = ntl.zz_p(4,7)
            sage: n*n ## indirect doctest
            4
        """
        zz_pContext_restore(&self.x)

def ntl_zz_pContext( v ):
    """
    Creation function for a zz_p context.

    EXAMPLES:
        sage: f = ntl.zz_pContext(26)
        sage: f = ntl.zz_pContext(10^100)
        Traceback (most recent call last):
        ...
        ValueError: Modulus (=10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000) is too big
    """
    if v > NTL_SP_BOUND:
        raise ValueError, "Modulus (=%s) is too big"%v
    if isinstance(v, Integer):
        v = mpz_get_si((<Integer>v).value)
    try:
        return zz_pContextDict[repr(v)]
    except KeyError:
        return ntl_zz_pContext_class(v)
