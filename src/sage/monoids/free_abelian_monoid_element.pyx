"""
Abelian Monoid Elements

AUTHORS:

- David Kohel (2005-09)

EXAMPLES:

Recall the example from abelian monoids::

    sage: F = FreeAbelianMonoid(5,names = list("abcde"))
    sage: (a,b,c,d,e) = F.gens()
    sage: a*b^2*e*d
    a*b^2*d*e
    sage: x = b^2*e*d*a^7
    sage: x
    a^7*b^2*d*e
    sage: x.list()
    [7, 2, 0, 1, 1]

The list is a copy, so changing the list does not change the element::

    sage: x.list()[0] = 0
    sage: x
    a^7*b^2*d*e
"""

#*****************************************************************************
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.memory cimport check_allocarray, sig_free
from cysignals.signals cimport sig_on, sig_off
from sage.structure.richcmp cimport rich_to_bool, richcmp
from sage.rings.integer cimport Integer, _Integer_from_mpz
from sage.libs.gmp.mpz cimport *

def is_FreeAbelianMonoidElement(x):
    r"""
    Queries whether ``x`` is an object of type ``FreeAbelianMonoidElement``.

    INPUT:

    - ``x`` -- an object.

    OUTPUT:

    - ``True`` if ``x`` is an object of type ``FreeAbelianMonoidElement``;
      ``False`` otherwise.
    """
    return isinstance(x, FreeAbelianMonoidElement)

cdef class FreeAbelianMonoidElement(MonoidElement):
    def __init__(self, parent, x):
        r"""
        Create the element ``x`` of the FreeAbelianMonoid ``parent``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: F
            Free abelian monoid on 5 generators (a, b, c, d, e)
            sage: F(1)
            1
            sage: a, b, c, d, e = F.gens()
            sage: a^2 * b^3 * a^2 * b^4
            a^4*b^7
            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: a, b, c, d, e = F.gens()
            sage: a in F
            True
            sage: a*b in F
            True
        """
        cdef Py_ssize_t i
        cdef Integer z

        MonoidElement.__init__(self, parent)
        cdef int n = <int> self._parent._FreeAbelianMonoid_class__ngens
        if isinstance(x, (int, Integer)) and x == 1:
            x = tuple([0]*n)
        if isinstance(x, (list, tuple)):
            if len(x) != n:
                raise IndexError("argument length (= %s) must be %s" % (len(x), n))
            self._element_vector = Vector_integer_dense(parent._module, x)
            self._element_vector.set_immutable()
        else:
            raise TypeError("argument x (= %s) is of wrong type"%x)

    def __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: a, b, c, d, e = F.gens()
            sage: x = a^2 * b^3
            sage: copy(x) == x
            True
            sage: copy(x) is x
            False
        """
        cdef type t = type(self)
        cdef FreeAbelianMonoidElement x = <FreeAbelianMonoidElement>(t.__new__(t))
        x._parent = self._parent
        x._element_vector = self._element_vector  # We can share since it is immutable
        return x

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: a, b, c, d, e = F.gens()
            sage: x = a^2 * b^3
            sage: loads(dumps(x)) == x
            True
        """
        return (self._parent, (self.list(),))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: F(1)
            1
            sage: a, b, c, d, e = F.gens()
            sage: a^2 * b^3 * a^2 * b^4
            a^4*b^7
        """
        s = ""
        A = self._parent
        x = A.variable_names()
        for i, val in enumerate(self._element_vector):
            if val == 0:
                continue
            elif val == 1:
                if len(s) > 0: s += "*"
                s += "%s" % x[i]
            else:
                if len(s) > 0: s += "*"
                s += "%s^%s" % (x[i], val)
        if not s:
            s = "1"
        return s

    cpdef _richcmp_(left, right, int op):
        """
        Rich comparison.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: F(1)
            1
            sage: a, b, c, d, e = F.gens()
            sage: x = a^2 * b^3
            sage: F(1) < x
            True
            sage: x > b
            True
            sage: x <= a^4
            True
            sage: x != a*b
            True
            sage: a*b == b*a
            True
            sage: x > a^3*b^2
            False
        """
        return richcmp(left._element_vector,
                       (<FreeAbelianMonoidElement>right)._element_vector,
                       op)
        #cdef Py_ssize_t i
        #cdef int c
        #for i in range(left._n):
        #    c = mpz_cmp(left._element_vector[i], (<FreeAbelianMonoidElement>right)._element_vector[i])
        #    if c < 0:
        #        return rich_to_bool(op, -1)
        #    elif c > 0:
        #        return rich_to_bool(op, 1)
        #return rich_to_bool(op, 0)

    def __mul__(self, y):
        """
        Multiply ``self`` with ``y``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: a, b, c, d, e = F.gens()
            sage: b * a^2 * b^3
            a^2*b^4
        """
        if not isinstance(y, FreeAbelianMonoidElement):
            raise TypeError("argument y (= %s) is of wrong type"%y)
        cdef type t = type(self)
        cdef FreeAbelianMonoidElement x = <FreeAbelianMonoidElement>(t.__new__(t))
        cdef FreeAbelianMonoidElement s, r
        s = self
        r = y
        x._parent = s._parent
        x._element_vector = s._element_vector._add_(r._element_vector)
        return x

    def __pow__(self, n, modulus):
        """
        Raises self to the power of `n`.

        AUTHORS:

        - Tom Boothby (2007-08): Replaced O(log n) time, O(n) space
          algorithm with O(1) time and space"algorithm".

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5,names = list("abcde"))
            sage: (a,b,c,d,e) = F.gens()
            sage: x = a*b^2*e*d; x
            a*b^2*d*e
            sage: x^3
            a^3*b^6*d^3*e^3
            sage: x^0
            1
        """
        if modulus is not None:
            raise NotImplementedError("modulus for exponents not implemented")
        cdef Integer val
        if isinstance(n, int):
            val = Integer(n)
        else:
            if not isinstance(n, Integer):
                raise TypeError("argument n (= %s) must be an integer" % (n,))
            val = <Integer> n
        if val < 0:
            raise IndexError("argument n (= %s) must be positive" % val)
        elif val == 1:
            return self
        cdef type t = type(self)
        cdef FreeAbelianMonoidElement x = <FreeAbelianMonoidElement>(t.__new__(t))
        cdef FreeAbelianMonoidElement s = self
        x._parent = s._parent
        x._element_vector = s._element_vector._rmul_(val)
        return x

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5,names = list("abcde"))
            sage: (a,b,c,d,e) = F.gens()
            sage: x = a*b^2*e*d
            sage: hash(x) == hash(x)
            True
        """
        return hash(tuple(self.list()))

    def list(self):
        """
        Return (a reference to) the underlying list used to represent this
        element. If this is a monoid in an abelian monoid on `n`
        generators, then this is a list of nonnegative integers of length
        `n`.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(5, 'abcde')
            sage: (a, b, c, d, e) = F.gens()
            sage: a.list()
            [1, 0, 0, 0, 0]
        """
        return list(self._element_vector)

