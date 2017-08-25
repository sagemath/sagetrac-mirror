#cython: boundscheck=False, wraparound=False
"""
Symmetric functions (Cython file)

This provides some low-level Cython implementations of functions
that are useful for symmetric function computations and need to
be optimized.
"""

from cpython.list cimport *
from cpython.dict cimport *

from sage.combinat.permutation_cython cimport next_perm
from sage.combinat.partition import _Partitions

from sage.misc.misc_c import prod
from sage.arith.misc import factorial
from collections import Counter
from sage.rings.integer import Integer
from sage.combinat.permutation import Permutations

cdef inline list add_lists(list la, list mu):
    cdef int i
    return [la[i] + mu[i] for i in xrange(len(la))]

cdef check_and_add(list la, list mu, dict ret):
    cdef int i
    cdef list nu = []
    nu.append(la[0] + mu[0])
    cdef int last = int(nu[0])
    cdef int d
    for i in xrange(1,len(la)):
        d = int(la[i] + mu[i])
        # We do not need to check for < 0 because la and mu are
        #   made from non-negative integers
        if last < d:  # Not a partition
            return
        if d > 0:
            nu.append(d)
        last = d
    # Note that nu[0] will always be > 0 because la and mu are non-empty
    nup = _Partitions(nu)
    if nup in ret:
        PyDict_SetItem(ret, nup, ret[nup] + 1)
    else:
        ret[nup] = 1

cpdef dict product_monomial_monomial(list la, list mu):
    """
    Return the product of monomial symmetric functions ``m[la] * m[mu]``.

    INPUT:

    - ``la``, ``mu`` -- partition given as a list (with no trailing 0's)
      can be mutated
    """
    cdef dict ret = {}
    cdef int ell = len(la) + len(mu)
    cdef list mu2, nu

    if not la:
        return {_Partitions(mu): 1}
    if not mu:
        return {_Partitions(la): 1}

    # Pad with 0's until both have length ell
    la += [0] * (ell - len(la))
    mu += [0] * (ell - len(mu))

    # Flip the order of the partitions because next_perm assumes
    #   lex-ordered is the identity and revlex is the last permutation
    la.reverse()
    mu.reverse()

    check_and_add(la, mu, ret)
    while next_perm(la):
        mu2 = list(mu)
        check_and_add(la, mu2, ret)
        while next_perm(mu2):
            check_and_add(la, mu2, ret)

    return ret

cdef dict c_prime(list mu, list nu):
    cdef dict ret = {}
    cdef int ell = len(mu) + len(nu)

    if not mu:
        return {_Partitions(mu): 1}
    if not nu:
        return {_Partitions(nu): 1}

    # Pad with 0's until both have length ell
    mu += [0] * (ell - len(mu))
    nu += [0] * (ell - len(nu))

    nu.reverse()

    nup = _Partitions(sorted(add_lists(mu, nu), reverse=True))
    ret[nup] = 1

    while next_perm(nu):
        nup = _Partitions(sorted(add_lists(mu, nu), reverse=True))
        if nup in ret:
            PyDict_SetItem(ret, nup, ret[nup] + 1)
        else:
            ret[nup] = 1

    return ret

cdef int c(list la, int n):
    cdef int ret = 1
    cdef int i, m
    cdef dict mult = {}
    for i in xrange(len(la)):
        m = <int> la[i]
        if m in mult:
            PyDict_SetItem(mult, m, mult[m] + 1)
        else:
            PyDict_SetItem(mult, m, 1)
    for m in mult.values():
        ret *= factorial(m)
    return ret * factorial(n - len(la))

cpdef dict pmm(list mu, list nu):
    # Permute the partition with fewer distinct values
    if len(set(mu)) < len(set(nu)):
        mu, nu = nu, mu
    cdef int n = len(mu) + len(nu)
    cdef int denom = c(mu, n)
    cdef dict cp = c_prime(mu, nu)

    for la in cp:
        PyDict_SetItem(cp, la, Integer(c(la._list, n) * cp[la] / denom))

    return cp

