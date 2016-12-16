r"""
Cython wrapper for bernmm library

AUTHOR:

    - David Harvey (2008-06): initial version
"""

#*****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#                     2008 David Harvey <dmharvey@math.harvard.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/cdefs.pxi"
include "cysignals/signals.pxi"


cdef extern from "bernmm/bern_rat.h":
    void bern_rat "bernmm::bern_rat" (mpq_t res, long k, int num_threads)

cdef extern from "bernmm/bern_modp.h":
    long bern_modp "bernmm::bern_modp" (long p, long k)



from sage.rings.rational cimport Rational


def bernmm_bern_rat(long k, int num_threads = 1):
    r"""
    Computes k-th Bernoulli number using a multimodular algorithm.
    (Wrapper for bernmm library.)

    INPUT:

    - k -- non-negative integer
    - num_threads -- integer >= 1, number of threads to use

    COMPLEXITY:

        Pretty much quadratic in $k$. See the paper "A multimodular algorithm
        for computing Bernoulli numbers", David Harvey, 2008, for more details.

    EXAMPLES::

        sage: from sage.rings.bernmm import bernmm_bern_rat

        sage: bernmm_bern_rat(0)
        1
        sage: bernmm_bern_rat(1)
        -1/2
        sage: bernmm_bern_rat(2)
        1/6
        sage: bernmm_bern_rat(3)
        0
        sage: bernmm_bern_rat(100)
        -94598037819122125295227433069493721872702841533066936133385696204311395415197247711/33330
        sage: bernmm_bern_rat(100, 3)
        -94598037819122125295227433069493721872702841533066936133385696204311395415197247711/33330

    TESTS::

        sage: lst1 = [ bernoulli(2*k, algorithm='bernmm', num_threads=2) for k in [2932, 2957, 3443, 3962, 3973] ]
        sage: lst2 = [ bernoulli(2*k, algorithm='pari') for k in [2932, 2957, 3443, 3962, 3973] ]
        sage: lst1 == lst2
        True
        sage: [ Zmod(101)(t) for t in lst1 ]
        [77, 72, 89, 98, 86]
        sage: [ Zmod(101)(t) for t in lst2 ]
        [77, 72, 89, 98, 86]
    """
    cdef Rational x

    if k < 0:
        raise ValueError("k must be non-negative")

    x = Rational()
    sig_on()
    bern_rat(x.value, k, num_threads)
    sig_off()

    return x


def bernmm_bern_modp(long p, long k):
    r"""
    Computes $B_k \mod p$, where $B_k$ is the k-th Bernoulli number.

    If $B_k$ is not $p$-integral, returns -1.

    INPUT:

        p -- a prime
        k -- non-negative integer

    COMPLEXITY:

        Pretty much linear in $p$.

    EXAMPLES::

        sage: from sage.rings.bernmm import bernmm_bern_modp

        sage: bernoulli(0) % 5, bernmm_bern_modp(5, 0)
        (1, 1)
        sage: bernoulli(1) % 5, bernmm_bern_modp(5, 1)
        (2, 2)
        sage: bernoulli(2) % 5, bernmm_bern_modp(5, 2)
        (1, 1)
        sage: bernoulli(3) % 5, bernmm_bern_modp(5, 3)
        (0, 0)
        sage: bernoulli(4), bernmm_bern_modp(5, 4)
        (-1/30, -1)
        sage: bernoulli(18) % 5, bernmm_bern_modp(5, 18)
        (4, 4)
        sage: bernoulli(19) % 5, bernmm_bern_modp(5, 19)
        (0, 0)

        sage: p = 10000019; k = 1000
        sage: bernoulli(k) % p
        1972762
        sage: bernmm_bern_modp(p, k)
        1972762

    TESTS:

    Check that bernmm works with the new NTL single precision modular
    arithmetic from :trac:`19874`::

        sage: from sage.rings.bernmm import bernmm_bern_modp
        sage: bernmm_bern_modp(7, 128) == bernoulli(128) % 7
        True
    """
    cdef long x

    if k < 0:
        raise ValueError("k must be non-negative")

    sig_on()
    x = bern_modp(p, k)
    sig_off()

    return x


# ============ end of file
