r"""
Other miscellaneous arithmetic function that are implemented in C for speed.
"""
#*****************************************************************************
#       Copyright (C) 2006 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from integer cimport Integer

ZZ_0 = Integer(0)
ZZ_1 = Integer(1)

def two_squares_pyx(unsigned int n):
    r"""
    If ``n`` is a sum of two squares return a 2-tuple ``(i,j)`` of Sage integers
    so that `i^2 + j^2 = n`. Otherwise raise a ``ValueError``.

    .. NOTE::

        The algorithm used here is very naive and only has interest for small
        values of ``n``. For that reason, the input must fit into an ``unsigned
        int`` which might be smaller than `2^{32}=4294967296` or
        `2^{64}=18446744073709551616` depending on your plateform.

    .. SEEALSO::

        :func:`~sage.arith.two_squares` is much more suited for large integers
        of the input.

    EXAMPLES::

        sage: from sage.rings.arith_pyx import two_squares_pyx
        sage: two_squares_pyx(0)
        (0, 0)
        sage: two_squares_pyx(1)
        (0, 1)
        sage: two_squares_pyx(2)
        (1, 1)
        sage: two_squares_pyx(3)
        Traceback (most recent call last):
        ...
        ValueError: 3 is not a sum of 2 squares
        sage: two_squares_pyx(106)
        (5, 9)
    """
    cdef unsigned int i,ii,j,nn

    i = ii = 0
    while ii <= n:
        j = 0
        while j <= i:
            nn = ii + j*j
            if nn >= n:
                break
            j += 1

        if nn == n:
            return (Integer(j),Integer(i))

        i += 1
        ii = i*i

    raise ValueError("%d is not a sum of 2 squares"%n)

def three_squares_pyx(unsigned int n):
    r"""
    If ``n`` is a sum of three squares return a 3-tuple ``(i,j,k)`` of Sage integers
    so that `i^2 + j^2 + k^2 = n`. Otherwise raise a ``ValueError``.

    .. NOTE::

        The algorithm used here is very naive and only has interest for small
        values of ``n``. For that reason, the input must fit into an ``unsigned
        int`` which might be smaller than `2^{32}=4294967296` or
        `2^{64}=18446744073709551616` depending on your plateform.

    .. SEEALSO::

        :func:`~sage.arith.three_squares` is much more suited for large integers
        of the input.

    EXAMPLES::

        sage: from sage.rings.arith_pyx import three_squares_pyx
        sage: three_squares(0)
        (0, 0, 0)
        sage: three_squares(1)
        (0, 0, 1)
        sage: three_squares(2)
        (0, 1, 1)
        sage: three_squares(3)
        (1, 1, 1)
        sage: three_squares(4)
        (0, 0, 2)
        sage: three_squares(5)
        (0, 1, 2)
        sage: three_squares(6)
        (1, 1, 2)
        sage: three_squares(7)
        Traceback (most recent call last):
        ...
        ValueError: 7 is not a sum of 3 squares
        sage: three_squares(107)
        (3, 7, 7)
    """
    cdef unsigned int i,ii,j,jj,k,nn,nnn

    i = ii = 0
    while ii <= n:
        j = jj = 0
        while j <= i:
            k = 0
            nn = ii + jj
            while k <= j:
                nnn = nn + k*k
                if nnn >= n:
                    break
                k += 1
            if nnn == n or nn >= n:
                break
            j += 1
            jj = j*j

        if nnn == n:
            return (Integer(k),Integer(j),Integer(i))

        i += 1
        ii = i*i

    raise ValueError("%d is not a sum of 3 squares"%n)

def four_squares_pyx(unsigned int n):
    r"""
    Return a 4-tuple ``(i,j,k,l)`` of Sage integers such that `i^2 + j^2 + k^2
    +l^2 = n`.

    .. NOTE::

        The algorithm used here is very naive and only has interest for small
        values of ``n``. For that reason, the input must fit into an ``unsigned
        int`` which might be smaller than `2^{32}=4294967296` or
        `2^{64}=18446744073709551616` depending on your plateform.

    .. SEEALSO::

        :func:`~sage.arith.four_squares` is much more suited for large integers
        of the input.

    EXAMPLES::

        sage: from sage.rings.arith_pyx import four_squares_pyx
        sage: all(sum(i**2 for i in four_squares_pyx(n)) == n for n in xrange(500))
        True
    """
    cdef unsigned int i,ii,j,jj,k,kk,l,nn,nnn,nnnn

    i = ii = 0
    while ii <= n:
        j = jj = 0
        while j <= i:
            k = kk = 0
            nn = ii + jj
            while k <= j:
                l = 0
                nnn = nn + kk
                while l <= k:
                    nnnn = nnn + l*l
                    if nnnn >= n:
                        break
                    l += 1
                if nnnn == n or nnn >= n:
                    break
                k += 1
                kk = k*k
            if nnnn == n or nn >= n:
                break
            j += 1
            jj = j*j

        if nnnn == n:
            return (Integer(l),Integer(k),Integer(j),Integer(i))

        i += 1
        ii = i*i
