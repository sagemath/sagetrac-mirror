r"""
Poset Discrete Morse Theory

An implementation of Babson and Hersch's [BH2005]_ discrete Morse function
for the order complex of posets.

AUTHORS:

- Jason P. Smith (2018-09-03): initial version

REFFERENCES:

- [BH2005]_
- [SV2006]_
"""

# ****************************************************************************
#       Copyright (C) 2018 Jason P. Smith <jasonsmith.bath@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.misc_c cimport list_diff

def is_PL_ordering(list L):
    r"""
    Check if ``L`` is a PL ordering.

    INPUT:

    - ``L`` -- a list of lists of maximal chains

    OUTPUT:

    Return a pair:

    1. boolean; if ``L`` is a PL ordering;
    2. the index that causes a problem if ``False``, otherwise ``None``.

    EXAMPLES::

        sage: from sage.combinat.posets.discrete_morse_theory import is_PL_ordering
        sage: L = [[1,2,3], [1,2,4], [1,3,5,6], [1,3,4], [7,3,4]]
        sage: is_PL_ordering(L)
        (True, None)
        sage: L = [[1,2,3], [1,2,4], [1,5,6], [7,2,4], [1,8,9]]
        sage: is_PL_ordering(L)
        (False, 4)
    """
    cdef Py_ssize_t i, j, d, len_Li
    cdef list Li, Lj
    for i in range(2, len(L)):
        Li = <list?> L[i]
        len_Li = len(Li)
        d = list_diff(<list?> L[i-1], Li)
        j = i - 2
        while j >= 0:
            Lj = <list?> L[j]
            if d > 0 and (len(Lj) < d or Lj[d-1] != Li[d-1]):
                j = -1
            elif ((len(Lj) < d + 1 and len_Li < d + 1)
                  or (len(Lj) >= d + 1 and len_Li >= d + 1
                      and Lj[:d+1] == Li[:d+1])):
                return (False, i)
            else:
                j -= 1
    return (True, None)

cdef list minimal_skipped_intervals(list L, Py_ssize_t c, bint pure):
    r"""
    Return the minimal skipped intervals (MSI's) of ``L[c]``, where each
    element is a pair ``[s, t]`` which represents the MSI ``L[c][s:t]``.

    INPUT:

    - ``L`` -- a list of lists of maximal chains
    - ``c`` -- an integer
    - ``pure`` -- boolean; if the poset is pure (i.e., all maximal
      chains have the same length)
    """
    cdef list M = []
    cdef list chain = <list> L[c]
    cdef Py_ssize_t rank = len(chain)
    cdef Py_ssize_t s, t
    for t in range(1, rank+1):
        s = 0
        while s + t <= rank:
            if (not any(k[0] >= s and k[1] <= s + t for k in M)
                and is_msi(L, c, s, s+t, pure)):
                M.append([s,s+t])
            s += 1
    return M

cdef bint is_msi(list L, Py_ssize_t c, Py_ssize_t s, Py_ssize_t t, bint pure):
    r"""
    Return if ``L[c][s:t]`` is a minimal skipped interval.

    INPUT:

    - ``L`` -- a list of lists of maximal chains
    - ``c``, ``s``, ``t`` -- integers
    - ``pure`` -- boolean; if the poset is pure
    """
    cdef list chain = <list> L[c]
    cdef list subchain = chain[:s] + chain[t:]
    if pure:
        for j in range(c-1,-1,-1):
            chain = <list> L[j]
            if chain[:s] + chain[t:] == subchain:
                return True
    else:
        for j in range(c-1,-1,-1):
            if interval_removed(subchain, <list> L[j]):
                return True
    return False

cdef bint interval_removed(list A, list B):
    r"""
    Return if ``A`` equals ``B`` after removing a consecutive sublist.
    """
    cdef Py_ssize_t lA = len(A)
    cdef Py_ssize_t lB = len(B)
    cdef Py_ssize_t i
    if lA > lB:
        return False
    if lA == lB:
        return A == B
    return any(A == B[:i] + B[i+lB-lA:] for i in range(lB-lA+1))

def critical_chains(list L):
    r"""
    Return the critical chains, along with the J intervals and minimal
    skipped intervals with respect to ``L``.

    INPUT:

    - ``L`` -- a list of lists of maximal chains

    .. SEEALSO::

        :meth:`sage.combinat.posets.posets.FinitePosets.discrete_morse_theory`

    EXAMPLES::

        sage: from sage.combinat.posets.discrete_morse_theory import critical_chains
        sage: L = [[6, 7, 2, 3], [6, 7, 2, 4], [6, 1, 8, 9, 3], [6, 1, 8, 9, 4],
        ....:      [5, 1, 8, 9, 3], [5, 1, 8, 9, 4]]
        sage: critical_chains(L)
        ([],
         [[], [[2, 4]], [[1, 4]], [[1, 4], [4, 5]], [[0, 1]], [[0, 1], [2, 5]]],
         [[], [[2, 4]], [[1, 4]], [[1, 4], [2, 5]], [[0, 1]], [[0, 1], [2, 5]]])

    An example where the J intervals differ from the minimal
    skipped intervals is given next. Note that in this case there
    are no critical chains even though ``L[4]`` is covered by its
    minimal skipped interval, as we only class a chain as critical
    if it is covered by its J intervals::

        sage: L = [[10, 6, 7], [10, 2, 7], [5, 6, 7], [5, 4, 3], [5, 4, 2, 7],
        ....:      [1, 9, 8, 4, 3], [1, 9, 8, 4, 2, 7]]
        sage: critical_chains(L)
        ([6],
         [[],
          [[1, 2]],
          [[0, 1]],
          [[1, 3]],
          [[0, 2], [2, 3]],
          [[0, 3]],
          [[0, 3], [3, 6]]],
         [[],
          [[1, 2]],
          [[0, 1]],
          [[1, 3]],
          [[0, 2], [1, 3]],
          [[0, 3]],
          [[0, 3], [2, 6]]])
    """
    cdef Py_ssize_t i, j

    # Compute if L corresponds to a pure poset (i.e., all maximal chains
    #   have the same length).
    cdef Py_ssize_t ell = len(L[0])
    cdef bint pure = True
    for i in range(1, len(L)):
        if len(L[i]) != ell:
            pure = False
            break

    cdef list M = [minimal_skipped_intervals(L, i, pure) for i in range(len(L))]
    cdef list J = J_intervals(L, M)

    cdef list crit = []
    cdef set target
    for i in range(len(L)):
        target = set(range(len(<list> L[i])))
        if set(j for k in J[i] for j in range(k[0],k[1])) == target:
            crit.append(i)

    return crit, J, M

cdef list J_intervals(list L, list M):
    r"""
    Return a list of the J intervals of the chains ``L``.

    INPUT:

    - ``L`` -- a list of lists of maximal chains
    - ``M`` -- the minimal skipped intervals of ``L``
    """
    cdef list J = []
    cdef list x, xk
    cdef Py_ssize_t j, k, t
    for i in M:
        x = [list(val) for val in i]
        x.sort(key=lambda y: y[0])  # TODO: Maybe this could be "min"?
        j = 0
        while j < len(x):
            k = j + 1
            while k < len(x):
                xk = <list> x[k]
                if x[j][1] > xk[0]:
                    xk[0] = x[j][1]
                    if xk[0] >= xk[1]:
                        del x[k]
                    else:
                        t = k + 1
                        while t < len(x):
                            if xk[0] >= x[t][0] and xk[1] <= x[t][1]:
                                del x[t]
                            else:
                                t += 1
                k += 1
            j += 1
        J.append(x)
    return J

