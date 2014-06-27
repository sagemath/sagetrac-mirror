# -*- coding: utf-8 -*-
"""
Linkage for arithmetic in `\ZZ[x]/<f>` using FLINT's ``fmpz_poly_t``

AUTHOR:

- Martin Albrecht (2014-06)

"""
#*****************************************************************************
#       Copyright (C) 2014 Martin Albrecht <martinralbrecht@googlemail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.flint.fmpz_poly cimport *, fmpz_poly_t

include "sage/ext/stdsage.pxi"

cdef inline celement *celement_new(const fmpz_poly_t h):
    """
    Create a new element in `\\ZZ[x]/<h>`.

    INPUT:

    - ``h`` - the modulus (ignored).

    EXAMPLE::

        sage: P.<x> = ZZ[]; Q.<x> = P.quotient(x^16 + 1)
        sage: f = 2*x + 1
        sage: f<<1
        2*x^2 + x
        sage: f<<16
        -2*x - 1

    """
    cdef celement *f = <celement *>sage_malloc(sizeof(fmpz_poly_t))
    fmpz_poly_init(f)
    return f

cdef inline int celement_delete(fmpz_poly_t f, const fmpz_poly_t h):
    """
    INPUT:

    - ``f`` - an element.
    - ``h`` - the modulus (ignored).

    Delete an element in `\\ZZ[x]/<h>`.

    EXAMPLE::

      sage: P.<x> = ZZ[]; Q.<x> = P.quotient(-x^2 + x - 2)
      sage: for i in range(1024):
      ...      f = Q.random_element()
      ...      del f

    """
    fmpz_poly_clear(f)
    sage_free(f)

cdef inline int celement_construct(fmpz_poly_t f, const fmpz_poly_t h):
    """
    Create a new element in `\\ZZ[x]/<h>`.

    INPUT:

    - ``f`` - memory for an element.
    - ``h`` - the modulus (ignored).

    EXAMPLE::

        sage: P.<x> = ZZ[]; Q.<x> = P.quotient(x^16 + 1)

    """
    fmpz_poly_init(f)

cdef inline int celement_destruct(fmpz_poly_t f, const fmpz_poly_t h):
    """
    Clear the element `f`

    INPUT:

    - ``f`` - an element.
    - ``h`` - the modulus (ignored).

    EXAMPLE::

        sage: P.<x> = ZZ[]; Q.<x> = P.quotient(x^32 + 1)
        sage: del x

    """
    fmpz_poly_clear(f)

cdef inline int celement_red(fmpz_poly_t f, const fmpz_poly_t h) except -2:
    """
    Compute `f %  h`.

    INPUT:

    - ``f`` - an element.
    - ``h`` - the modulus (ignored).

    EXAMPLE::

        sage: P.<x> = ZZ[]; Q.<x> = P.quotient(x^32 + 1)
        sage: x^31 * x
        -1

    """
    fmpz_poly_rem(f, f, h)

cdef inline int celement_gen(fmpz_poly_t f, long i, const fmpz_poly_t h) except -2:
    """
    Return `x ∈ \\ZZ[x]/<h>`

    INPUT:

    - ``f`` - output.
    - ``i`` - index of generator (ignored).
    - ``h`` - the modulus (ignored).

    """
    fmpz_poly_zero(f)
    fmpz_poly_set_coeff_ui(f, 1, 1)

cdef object celement_repr(fmpz_poly_t e, const fmpz_poly_t n):
    raise NotImplementedError

cdef inline int celement_set(fmpz_poly_t f, fmpz_poly_t a, const fmpz_poly_t h) except -2:
    """
    Set `f = a % h`.

    INPUT:

    - ``f`` - output.
    - ``a`` - a polynomial over the integers.
    - ``n`` - the modulus.

    EXAMPLE::

        sage: P.<x> = ZZ[]
        sage: Q.<X> = P.quotient(x^128 + 1)
        sage: Q((x+1)^200) == (x+1)^200 % (x^128 + 1)
        True

    """
    fmpz_poly_set(f, a)
    celement_red(f, h)

cdef inline int celement_set_si(fmpz_poly_t f, long a, const fmpz_poly_t h) except -2:
    """
    Set `f = a % n`.

    INPUT:

    - ``f`` - an element.
    - ``a`` - an integer.
    - ``h`` - the modulus (ignored).

    EXAMPLE::

        sage: P.<x> = ZZ[]
        sage: Q.<X> = P.quotient(x^8 + 1) # indirect doctest

    """
    fmpz_poly_set_coeff_si(f, 0, a)

cdef inline long celement_get_si(fmpz_poly_t res, const fmpz_poly_t n) except -2:
    raise NotImplementedError

cdef inline bint celement_is_zero(fmpz_poly_t f, const fmpz_poly_t h) except -2:
    """
    Return ``True`` if `f == 0`.

    INPUT:

    - ``f`` - an element.
    - ``h`` - the modulus (ignored).

    EXAMPLE::

        sage: P.<x> = ZZ[]
        sage: Q.<X> = P.quotient(x^8 + 1)
        sage: X.is_zero()
        False
        sage: Q(0).is_zero()
        True

    """
    return fmpz_poly_is_zero(f)

cdef inline bint celement_is_one(fmpz_poly_t f, const fmpz_poly_t h) except -2:
    """
    Return ``True`` if `f == 1`.

    INPUT:

    - ``f`` - an element.
    - ``h`` - the modulus (ignored).

    EXAMPLE::

        sage: P.<x> = ZZ[]
        sage: Q.<X> = P.quotient(x^8 + 1)
        sage: X.is_one()
        False
        sage: Q(1).is_one()
        True

    """
    return fmpz_poly_is_one(f)

cdef inline bint celement_equal(fmpz_poly_t f, fmpz_poly_t g, const fmpz_poly_t h) except -2:
    """
    Return ``True`` if `f == g`.

    INPUT:

    - ``f`` - an element.
    - ``g`` - an element.
    - ``h`` - the modulus (ignored).

    EXAMPLE::

        sage: P.<x> = ZZ[]
        sage: Q.<X> = P.quotient(x^8 + 1)
        sage: X == Q(list(X))
        True
        sage: X+1 == X
        False

    """
    return fmpz_poly_equal(f, g)

cdef inline int celement_cmp(fmpz_poly_t f, fmpz_poly_t g, const fmpz_poly_t h) except -2:
    """
    Return an integer < 0  if `f < g`, an integer > 0 if `f > g` and 0 otherwise.

    INPUT:

    - ``f`` - an element.
    - ``g`` - an element.
    - ``h`` - the modulus (ignored).

    EXAMPLE::

        sage: P.<x> = ZZ[]
        sage: Q.<X> = P.quotient(x^16 + 24*x^15 - 9*x^14 + x^13 + x^12 - 5*x^10 - x^7 + 3*x^6 + x^5 - 9*x^4 + x^3 + 49*x^2 + x - 2)
        sage: f = Q.random_element()
        sage: g = -Q(list(-f))
        sage: f < g
        False
        sage: f <= g
        True
        sage: f > X
        True
        sage: f >= X
        True
        sage: f+1 >= f
        True
        sage: f+1 > f
        True

    """
    cdef long deg_f = fmpz_poly_degree(f)
    cdef long delta = deg_f - fmpz_poly_degree(g)
    cdef int i

    if delta > 0:
        return 1
    elif delta < 0:
        return -1

    if fmpz_poly_equal(f, g):
        return 0

    i = deg_f
    cdef fmpz *coeff_f = fmpz_poly_get_coeff_ptr(f, i)
    cdef fmpz *coeff_g = fmpz_poly_get_coeff_ptr(g, i)

    while fmpz_equal(coeff_f, coeff_g) and i > 0:
        i -= 1
        coeff_f = fmpz_poly_get_coeff_ptr(f, i)
        coeff_g = fmpz_poly_get_coeff_ptr(g, i)
    return fmpz_cmp(coeff_f, coeff_g)

cdef long celement_len(fmpz_poly_t f, const fmpz_poly_t h) except -2:
    """
    Return `\deg(f)+1`.

    INPUT:

    - ``f`` - an element.
    - ``h`` - the modulus (ignored).

    EXAMPLE::

        sage: P.<x> = ZZ[]
        sage: Q.<X> = P.quotient(x^16 + 24*x^15 - 9*x^14 + x^13 + x^12 - 5*x^10 - x^7 + 3*x^6 + x^5 - 9*x^4 + x^3 + 49*x^2 + x - 2)
        sage: X.degree()
        1
        sage: (X^15).degree()
        15
        sage: (X^17).degree() == ((x^17) % Q.modulus()).degree() == 15
        True
    """
    return <long>fmpz_poly_length(f)

cdef inline int celement_add(fmpz_poly_t r, fmpz_poly_t f, fmpz_poly_t g, const fmpz_poly_t h) except -2:
    """
    Compute `r = (f + g) % h`.

    INPUT:

    - ``r`` - output.
    - ``f`` - an element.
    - ``g`` - an element.
    - ``h`` - the modulus (ignored).

    EXAMPLE::
      
        sage: P.<x> = ZZ[]
        sage: Q.<Y> = P.quotient(x^8 + 1)
        sage: all([Q((x+1)^i + (2*x+3)^j) == Q((x+1)^i) + Q((2*x+3)^j) for i in range(16) for j in range(16)])
        True

    """
    fmpz_poly_add(r, f, g)

cdef inline int celement_sub(fmpz_poly_t r, fmpz_poly_t f, fmpz_poly_t g, const fmpz_poly_t h) except -2:
    """
    Compute `r = (f - g) % h`.

    INPUT:

    - ``r`` - output.
    - ``f`` - an element.
    - ``g`` - an element.
    - ``h`` - the modulus (ignored).

    EXAMPLE::
      
        sage: P.<x> = ZZ[]
        sage: Q.<Y> = P.quotient(x^8 + 1)
        sage: all([Q((x+1)^i - (2*x+3)^j) == Q((x+1)^i) - Q((2*x+3)^j) for i in range(16) for j in range(16)])
        True

"""
    fmpz_poly_sub(r, f, g)

cdef inline int celement_neg(fmpz_poly_t r, fmpz_poly_t f, const fmpz_poly_t n) except -2:
    """
    Return `r = -f`.

    INPUT:

    - ``r`` - output.
    - ``f`` - an element.
    - ``h`` - the modulus (ignored).

    EXAMPLE::

        sage: P.<x> = ZZ[]
        sage: Q.<Y> = P.quotient(x^8 + 1)
        sage: -Q(Y)
        -Y
    """
    fmpz_poly_neg(r, f)

cdef inline int celement_mul_scalar(fmpz_poly_t r, fmpz_poly_t f, object c, const fmpz_poly_t h) except -2:
    """
    Compute `r = c·f % h`.

    INPUT:

    - ``r`` - output.
    - ``f`` - an element.
    - ``c`` - a scalar.
    - ``h`` - the modulus (ignored).


    EXAMPLE::

        sage: P.<x> = ZZ[]
        sage: Q.<Y> = P.quotient(x^8 + 1)
        sage: all([Q((x+1)^i*2^j) == Q((x+1)^i)*Q(2^j) for i in range(16) for j in range(16)])
        True
      
    """
    cdef Integer c_z = Integer(c)
    fmpz_poly_scalar_mul_mpz(r, f, c_z.value)

cdef inline int celement_mul(fmpz_poly_t r, fmpz_poly_t f, fmpz_poly_t g, const fmpz_poly_t h) except -2:
    """
    Compute `r = (f · g) % h`.

    INPUT:

    - ``r`` - output.
    - ``f`` - an element.
    - ``g`` - an element.
    - ``h`` - the modulus.

    EXAMPLE::

        sage: P.<x> = ZZ[]
        sage: Q.<Y> = P.quotient(x^8 + 1)
        sage: all([Q((x+1)^i*(2*x+3)^j) == Q((x+1)^i)*Q((2*x+3)^j) for i in range(16) for j in range(16)])
        True

    """
    fmpz_poly_mul(r, g, f)
    celement_red(r, h)

cdef inline int celement_mod(fmpz_poly_t r, fmpz_poly_t f, fmpz_poly_t g, const fmpz_poly_t h) except -2:
    """
    Compute `r = (f % g) % h`.

    INPUT:

    - ``r`` - output.
    - ``f`` - an element.
    - ``g`` - an element.
    - ``h`` - the modulus (ignored).

    EXAMPLE::

        sage: P.<x> = ZZ[]
        sage: Q.<Y> = P.quotient(x^32 + 1)
        sage: f = P.random_element(degree=16); f
        x^16 - x^15 + x^14 - 12*x^11 - 2*x^10 - x^9 - 95*x^8 + x^7 + 2*x^6 - x^5 + x^4 + 2*x - 8
        sage: g = P.random_element(degree=8); g
        5*x^6 - 6*x^5 - 4*x^4 + 4*x^3 - x^2 - 2*x - 1
        sage: r = Q(f) % Q(g); r
        Y^16 - Y^15 + Y^14 + 3*Y^11 + 3*Y^9 + 3*Y^8 + 2*Y^7 + 4*Y^6 - 468*Y^5 - 161*Y^4 + 164*Y^3 - 178*Y^2 - 174*Y - 76
        sage: f % g
        x^16 - x^15 + x^14 + 3*x^11 + 3*x^9 + 3*x^8 + 2*x^7 + 4*x^6 - 468*x^5 - 161*x^4 + 164*x^3 - 178*x^2 - 174*x - 76
    """
    fmpz_poly_rem(r, f, g)

cdef inline int celement_pow(fmpz_poly_t r, fmpz_poly_t f, long e, fmpz_poly_t modulus, const fmpz_poly_t h) except -2:
    """
    Compute `r = f^e % h`.

    INPUT:

    - ``r`` - output.
    - ``f`` - an element.
    - ``e`` - an exponent.
    - ``modulus`` - NULL.
    - ``h`` - the modulus (ignored).


    EXAMPLE::

        sage: P.<x> = ZZ[]
        sage: Q.<X> = P.quotient(x^32+ 1)
        sage: f = x^31 + 1
        sage: all([Q(f)^e == Q(f^e) for e in range(128)])
        True

        sage: P.<x> = ZZ[]
        sage: Q.<X> = P.quotient(x^16 + 24*x^15 - 9*x^14 + x^13 + x^12 - 5*x^10 - x^7 + 3*x^6 + x^5 - 9*x^4 + x^3 + 49*x^2 + x - 2)
        sage: f = x^2 + 3*x + 1
        sage: all([Q(f)^e == Q(f^e) for e in range(128)])
        True
    """
    if modulus != NULL:
      raise NotImplementedError

    cdef fmpz_poly_t pow2
    cdef fmpz_poly_t q
    cdef fmpz_poly_t tmp

    fmpz_poly_init(q)
    fmpz_poly_init(tmp)

    if e == 0:
        fmpz_poly_zero(r)
        fmpz_poly_set_coeff_ui(r, 0, 1)
    elif e == 1:
        fmpz_poly_set(r, f)
    elif e == 2:
        fmpz_poly_pow(r, f, 2)
        celement_red(r, h)
    else:
        if r == f:
            fmpz_poly_set(tmp, f)
            f = tmp
        fmpz_poly_init(pow2)
        fmpz_poly_set(pow2, f)
        if e % 2:
            fmpz_poly_set(r, f)
        else:
            fmpz_poly_zero(r)
            fmpz_poly_set_coeff_ui(r, 0, 1)
        e = e >> 1
        while(e != 0):
            fmpz_poly_pow(pow2, pow2, 2)
            if e % 2:
                fmpz_poly_mul(r, r, pow2)
            e = e >> 1
            celement_red(r, h)
        fmpz_poly_clear(pow2)

    celement_red(r, h)
    fmpz_poly_clear(q)
    fmpz_poly_clear(tmp)

cdef inline int celement_div(fmpz_poly_t r, fmpz_poly_t a, fmpz_poly_t b, const fmpz_poly_t n) except -2:
    raise TypeError("Division is undefined.")

cdef inline int celement_floordiv(fmpz_poly_t r, fmpz_poly_t f, fmpz_poly_t g, const fmpz_poly_t h) except -2:
    raise TypeError("Floor division is undefined.")

cdef inline int celement_quorem(fmpz_poly_t q, fmpz_poly_t r, fmpz_poly_t a, fmpz_poly_t b, const fmpz_poly_t n) except -2:
    raise TypeError("Quotients are undefined.")

cdef inline int celement_inv(fmpz_poly_t r, fmpz_poly_t a, const fmpz_poly_t n) except -2:
    raise TypeError("Inversion is undefined.")

cdef inline int celement_gcd(fmpz_poly_t r, fmpz_poly_t a, fmpz_poly_t b, fmpz_poly_t n) except -2:
    raise TypeError("GCD undefined.")

cdef inline int celement_xgcd(fmpz_poly_t r, fmpz_poly_t s, fmpz_poly_t t, fmpz_poly_t a, fmpz_poly_t b, fmpz_poly_t n) except -2:
    raise TypeError("XGCD undefined.")
