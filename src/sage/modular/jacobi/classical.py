r"""
Fourier expansions of classical Jacobi forms.

AUTHOR:

- Martin Raum

EXAMPLES:

To compute a basis of weight `k` and index `m \in \ZZ` Jacobi forms,
we call ``classical_jacobi_forms``.  This compute weight `10`, index
`2` Jacobi forms up to precision `5`.

::

    sage: from sage.modular.jacobi.all import classical_jacobi_forms
    sage: classical_jacobi_forms(10, 2, 5)
    [{(0, 0): 823680,
      (1, 0): -132447744,
      (1, 1): -42172416,
      (1, 2): -329472,
      (2, 0): -47908194048,
      (2, 1): -27622932480,
      (2, 2): -4157936640,
      (3, 0): -1497954991104,
      (3, 1): -1045412020224,
      (3, 2): -318012963840,
      (4, 0): -17346153217536,
      (4, 1): -13216244760576,
      (4, 2): -5574898188288},
     {(1, 0): -22464,
      (1, 1): 9984,
      (1, 2): 1248,
      (2, 0): -19968,
      (2, 1): 149760,
      (2, 2): -149760,
      (3, 0): 913536,
      (3, 1): -1707264,
      (3, 2): 1123200,
      (4, 0): 359424,
      (4, 1): 4253184,
      (4, 2): -2715648}]

For description of Fourier expansions, see the documentation of
classical weak Jacobi forms.
"""

#===============================================================================
#
# Copyright (C) 2010-2014 Martin Raum
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

from sage.combinat.dict_addition import dict_linear_combination
from sage.matrix.all import matrix
from sage.misc.all import isqrt, cached_function
from sage.modular.jacobi.classical_weak import classical_weak_jacobi_fe_indices, classical_weak_jacobi_forms
from sage.rings.all import Integer, ZZ

# We do not implement this separately, because this is the same
# reduction as in the case of weak Jacobi forms.
# from sage.modular.jacobi.classical_weak import classical_jacobi_reduce_fe_index


def classical_jacobi_fe_indices(m, prec, reduced=False):
    r"""
    Indices `(n,r)` of Fourier expansions of Jacobi forms of index `m`.

    INPUT:

    - `m` -- A positive integer.

    - ``prec`` -- A non-negative integer.

    - ``reduce`` -- A boolean (default: ``False``).  If ``True``
      restrict to `0 \le r \le m`.

    OUTPUT:

    A generator of pairs of integers `(n,r)`.

    EXAMPLES::

        sage: from sage.modular.jacobi.all import *
        sage: list(classical_jacobi_fe_indices(2, 3, True))
        [(1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (0, 0)]
        sage: list(classical_jacobi_fe_indices(2, 3, False))
        [(1, 0), (1, 1), (1, -1), (1, 2), (1, -2), (2, 0), (2, 1), (2, -1),
        (2, 2), (2, -2), (2, 3), (2, -3), (0, 0), (2, 4), (2, -4)]
    """
    fm = Integer(4 * m)

    if reduced:
        # positive definite forms
        for n in range(1, prec):
            for r in range(min(m + 1, isqrt(fm * n - 1) + 1)):
                yield (n, r)

        # indefinite forms
        for r in xrange(0, min(m+1, isqrt((prec-1) * fm) + 1) ):
            if fm.divides(r ** 2):
                yield (r ** 2 // fm, r)
    else:
        # positive definite forms
        for n in range(1, prec):
            yield(n, 0)
            for r in range(1, isqrt(fm * n - 1) + 1):
                yield (n, r)
                yield (n, -r)

        # indefinite forms
        if prec > 0:
            yield (0, 0)
        for n in xrange(1, prec):
            if (fm * n).is_square():
                rt_fmm = isqrt(fm * n)
                yield(n, rt_fmm)
                yield(n, -rt_fmm)

    raise StopIteration


@cached_function
def _classical_jacobi_forms_as_weak_jacobi_forms(k, m, algorithm="skoruppa"):
    r"""
    The coordinates of Jacobi forms with respect to a basis of
    classical weak Jacobi forms computed using the given algorithm.

    INPUT:

    - `k` -- An integer.

    - `m` -- A non-negative integer.

    - ``algorithm`` -- Default: ''skoruppa''.  Only ''skoruppa'' is implemented.

    OUTPUT:

    A list of vectors.

    EXAMPLES::

        sage: from sage.modular.jacobi.classical import _classical_jacobi_forms_as_weak_jacobi_forms
        sage: _classical_jacobi_forms_as_weak_jacobi_forms(10, 1)
        [
        (1, 0, 0),
        (0, 0, 1)
        ]
        sage: _classical_jacobi_forms_as_weak_jacobi_forms(12, 1)
        [
        (1, 0, 0),
        (0, 1, 0)
        ]
    """
    prec = m//4 + 1
    weak_forms = classical_weak_jacobi_forms(k, m, prec, algorithm)
    indices = list(classical_weak_jacobi_fe_indices(m, prec, reduced=True))
    weak_index_matrix = \
        matrix(ZZ, [ [ (f[(n,r)] if (n,r) in f else 0) for (n,r) in indices
                       if 4*m*n - r**2 < 0 ] for f in weak_forms] )

    return weak_index_matrix.left_kernel().echelonized_basis()

def classical_jacobi_forms(k, m, prec, algorithm="skoruppa"):
    r"""
    Compute Fourier expansions of classical Jacobi forms of weight `k`
    and index `m`.

    INPUT:

    - `k` -- An integer.

    - `m` -- A non-negative integer.

    - ``prec`` -- A non-negative integer that corresponds to a
      precision of the `q`-expansion.

    - ``algorithm`` -- Default: ''skoruppa''.  Only ''skoruppa'' is
      implemented.

    OUTPUT:

    A list of dictionaries, mapping indices `(n,r)` to rationals.

    EXAMPLES::

        sage: from sage.modular.jacobi.classical import *
        sage: k = 4; m = 2; prec = 5
        sage: classical_jacobi_forms(k, m, prec)
        [{(0, 0): 40320,
          (1, 0): 3386880,
          (1, 1): 2580480,
          (1, 2): 564480,
          (2, 0): 23143680,
          (2, 1): 18063360,
          (2, 2): 11289600,
          (3, 0): 51932160,
          (3, 1): 54190080,
          (3, 2): 33868800,
          (4, 0): 138862080,
          (4, 1): 108380160,
          (4, 2): 95477760}]

    TESTS:

    See ``test_classical.py``.
    """
    weak_forms = classical_weak_jacobi_forms(k, m, prec, algorithm)
    coords = _classical_jacobi_forms_as_weak_jacobi_forms(k, m)

    return [dict_linear_combination(zip(weak_forms, cs)) for cs in coords]
