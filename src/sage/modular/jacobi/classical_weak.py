r"""
Fourier expansions of classical weak Jacobi forms.

AUTHOR:

- Martin Raum

REFERENCES:

.. [Sk84] Skoruppa, Uber den Zusammenhang zwischen Jacobiformen und
   Modulformen halbganzen Gewichts, 1984, University of Bonn.

EXAMPLES:

To compute a basis of weight `k` and index `m \in \ZZ` weak Jacobi
forms, we call ``classical_weak_jacobi_forms``.  This compute weight
`9`, index `2` Jacobi forms up to precision `5`.

::

    sage: from sage.modular.jacobi.all import *
    sage: jforms = classical_weak_jacobi_forms(9, 2, 5)
    sage: jforms
    [{(0, 1): 660,
      (1, 1): -172260,
      (2, 1): -89901900,
      (3, 1): -3699449160,
      (4, 1): -56862841860}]


Fourier expansions are represented as dictionaries.  Fourier
coefficients of Jacobi forms are index by pairs `(n,r)` of integers.
We call `(n,r)` reduced if `0 \le r \le m`.  Fourier expansions are
given in terms of Fourier coefficients of reduced index.  They can be
accessed directly; If a reduced pair does not occur as a key and `n`
does not exceed the prescribed precision, then the corresponding
Fourier coefficient vanishes.

To access Fourier coefficients of non reduced index, we compute the attached reduction by ``classical_jacobi_reduce_fe_index``::

    sage: classical_jacobi_reduce_fe_index((2, 3), 2)
    ((1, 1), -1)

As a result, we obtain `(n', r')` and a sign `s`.  The Fourier
coefficient at `(n,r)` equals `s^k` times the coefficient at `(n',
r')`.  In the specific examples (`k` is odd), the `(2,3)`-th Fourier
coefficient of the first basis element is equal to::

   sage: -jforms[0][(1,1)]
   172260
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

from sage.combinat.all import number_of_partitions
from sage.functions.all import binomial, factorial
from sage.matrix.all import matrix
from sage.misc.all import cached_function
from sage.misc.all import prod, isqrt
from sage.misc.weak_dict import WeakValueDictionary
from sage.modular.all import ModularForms
from sage.rings.all import (PowerSeriesRing, Integer, ZZ, QQ)
from sage.rings import big_oh
import operator


def classical_jacobi_reduce_fe_index((n, r), m):
    r"""
    Reduce an index `(n, r)` of the Fourier expansion of Jacobi forms
    of index `m` suchthat `0 \le r \le m`.

    INPUT:

    - `(n, r)` -- A pair of integers.

    - `m` -- A positive integer.

    OUTPUT:

    A pair `((n', r'), s)`, where `(n', r')` is an index that is
    equivalent to `(n,r)` and `s` is a sign indicating whether `r -
    r'` or `r + r'` is divisible by `2m`.

    EXAMPLES::

        sage: from sage.modular.jacobi.all import *
        sage: classical_jacobi_reduce_fe_index((2,2), 1)
        ((1, 0), 1)
    """
    rred = r % (2*m)

    if rred > m:
        sgn = -1
        rred = 2*m - rred
    else:
        sgn = 1

    nred = n - (r**2 - rred**2)//(4*m)

    return ((nred, rred), sgn)


def classical_weak_jacobi_fe_indices(m, prec, reduced=False):
    r"""
    Indices `(n,r)` of Fourier expansions of weak Jacobi forms of index `m`.

    INPUT:

    - `m` -- A positive integer.

    - ``prec`` -- A non-negative integer.

    - ``reduce`` -- A boolean (default: ``False``).  If ``True``
      restrict to `0 \le r \le m`.

    OUTPUT:

    A generator of pairs `(n,r)`.

    EXAMPLES::

        sage: from sage.modular.jacobi.all import *
        sage: list(classical_weak_jacobi_fe_indices(2, 3, True))
        [(1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (0, 0), (0, 1), (0, 2)]
        sage: list(classical_weak_jacobi_fe_indices(2, 3, False))
        [(1, 0), (1, 1), (1, -1), (1, 2), (1, -2), (2, 0), (2, 1), (2, -1), (2, 2), (2, -2), (2, 3), (2, -3), (0, -1), (1, 3), (2, -4), (0, 0), (2, 4), (1, -3), (0, 1), (0, -2), (0, 2)]
    """
    fm = Integer(4*m)

    if reduced:
        # positive definite forms
        for n in range(1, prec):
            for r in range(min(m + 1, isqrt(fm * n - 1) + 1)):
                yield (n, r)

        # indefinite forms
        for n in xrange(0, prec):
            for r in xrange( isqrt(fm * n - 1) + 1 if n != 0 else 0, m + 1 ):
                yield (n, r)
    else:
        # positive definite forms
        for n in range(1, prec):
            yield(n, 0)
            for r in range(1, isqrt(fm * n - 1) + 1):
                yield (n, r)
                yield (n, -r)

        # indefinite forms
        for n in xrange(0, min(m // 4 + 1, prec)):
            if n == 0:
                r_iteration = range(-m + 1, m + 1)
            else:
                r_iteration = range(-m + 1, -isqrt(fm * n - 1))
                r_iteration += range(isqrt(fm * n - 1) + 1, m + 1)
            for r in r_iteration:
                for l in range((- r - isqrt(r**2 - 4 * m * (n - (prec - 1))) - 1) // (2 * m) + 1,
                                (- r + isqrt(r**2 - 4 * m * (n - (prec - 1)))) // (2 * m) + 1 ):
                    yield (n + l * r + m * l ** 2, r + 2 * m * l)

    raise StopIteration


## TODO: Implement caching by truncation.
def classical_weak_jacobi_forms(k, m, prec, algorithm="skoruppa"):
    r"""
    INPUT:

    - `k` -- An integer.

    - `m` -- A non-negative integer.

    - ``prec`` -- A non-negative integer that corresponds to a precision of
      the q-expansion.

    - ``algorithm`` -- Default: ''skoruppa''.  Only ''skoruppa'' is implemented.

    EXAMPLES::

        sage: from sage.modular.jacobi.all import *
        sage: classical_weak_jacobi_forms(2, 1, 5)
        [{(0, 0): -4,
          (0, 1): 2,
          (1, 0): -984,
          (1, 1): 496,
          (2, 0): -14512,
          (2, 1): 8238,
          (3, 0): -106016,
          (3, 1): 67024,
          (4, 0): -574488,
          (4, 1): 385026}]

    TESTS:

    See ``test_classical_weak.py``.
    """
    if algorithm != "skoruppa":
        raise NotImplementedError("Algorithm {} is not implemented.".format(algorithm))
    factory = ClassicalWeakJacobiFormsFactory(m, prec)

    return [ factory.from_taylor_expansion(fs, k, is_integral=True)
            for fs in _classical_weak_jacobi_taylor_coefficients(k, m)]


@cached_function
def _classical_weak_jacobi_taylor_coefficients(k, m):
    r"""
    A product basis of the echelon bases of

    - `M_k, M_{k + 2}, ..., M_{k + 2 m}` etc. if ``weight`` is even,

    - `M_{k + 1}, ..., M_{k + 2 m - 3}` if ``weight`` is odd.

    INPUT:

    - `k` -- An integer.

    - `m` -- A non-negative integer.

    OUTPUT:

    A list of lists of functions.  When called with with an argument
    ``prec``, each of them returns a `q`-expansion of that precision.

    EXAMPLES::

        sage: from sage.modular.jacobi.classical_weak import _classical_weak_jacobi_taylor_coefficients
        sage: _classical_weak_jacobi_taylor_coefficients(12, 1)
        [[<bound method ModularFormElement.qexp of 1 + 196560*q^2 + 16773120*q^3 + 398034000*q^4 + 4629381120*q^5 + O(q^6)>, <function <lambda> at ...>], [<bound method ModularFormElement.qexp of q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)>, <function <lambda> at ...>], [<function <lambda> at ...>, <bound method ModularFormElement.qexp of 1 - 24*q - 196632*q^2 - 38263776*q^3 - 1610809368*q^4 - 29296875024*q^5 + O(q^6)>]]
    """
    R = PowerSeriesRing(ZZ, 'q')
    q = R.gen()

    if k % 2 == 0:
        nmb_modular_forms = m + 1
        start_weight = k
    else:
        nmb_modular_forms = m - 1
        start_weight = k + 1

    modular_forms = list()
    for (i, k) in enumerate(range(start_weight, start_weight + 2 * nmb_modular_forms, 2)):
        modular_forms += [ [lambda p: big_oh.O(q**p) for _ in range(i)] + [b.qexp] + [lambda p: big_oh.O(q**p) for _ in range(nmb_modular_forms - 1 - i)]
                           for b in ModularForms(1, k).echelon_basis() ]

    return modular_forms

_classical_weak_jacobi_forms_factory_cache = WeakValueDictionary()

def ClassicalWeakJacobiFormsFactory(m, prec):
    r"""
    INPUT:

    - `m` -- A non-negative integer.

    - ``prec`` -- An integer.

    EXAMPLES::

        sage: from sage.modular.jacobi.classical_weak import ClassicalWeakJacobiFormsFactory
        sage: fcty = ClassicalWeakJacobiFormsFactory(2, 5)

    TESTS::

        sage: fcty = ClassicalWeakJacobiFormsFactory(2, 5)
        sage: assert fcty is ClassicalWeakJacobiFormsFactory(2, 5)
    """
    if (m, prec) in _classical_weak_jacobi_forms_factory_cache:
        return _classical_weak_jacobi_forms_factory_cache[(m, prec)]

    for (mc, pc) in _classical_weak_jacobi_forms_factory_cache:
        if m == mc and pc > prec:
            return _classical_weak_jacobi_forms_factory_cache[(m, pc)]._truncate_cache(prec)

    fcty = ClassicalWeakJacobiForms_factory(m, prec)
    _classical_weak_jacobi_forms_factory_cache[(m, prec)] = fcty
    return fcty


class ClassicalWeakJacobiForms_factory:
    r"""
    A factory for Jacobi form of degree 1 and scalar index `m`.
    """
    def __init__(self, m, prec):
        r"""
        INPUT:

        - `m` -- A non-negative integer.

        - ``prec`` -- An integer.

        EXAMPLES::

            sage: from sage.modular.jacobi.classical_weak import ClassicalWeakJacobiForms_factory
            sage: fcty = ClassicalWeakJacobiForms_factory(2, 5)
        """
        self.__prec = prec
        self.__jacobi_index = m

        self.__power_series_ring_ZZ = PowerSeriesRing(ZZ, 'q')
        self.__power_series_ring = PowerSeriesRing(QQ, 'q')

    def jacobi_index(self):
        r"""
        The index of the Jacobi forms that are computed by this factory.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: from sage.modular.jacobi.classical_weak import ClassicalWeakJacobiForms_factory
            sage: fcty = ClassicalWeakJacobiForms_factory(2, 5)
            sage: fcty.jacobi_index()
            2
        """
        return self.__jacobi_index

    def precision(self):
        r"""
        The precision of Fourier expansions that are computed.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: from sage.modular.jacobi.classical_weak import ClassicalWeakJacobiForms_factory
            sage: fcty = ClassicalWeakJacobiForms_factory(2, 5)
            sage: fcty.precision()
            5
        """
        return self.__prec

    def _power_series_ring(self):
        r"""
        An auxiliary power series ring that is cached in the factory.

        OUTPUT:

        An power series ring.

        EXAMPLES::

            sage: from sage.modular.jacobi.classical_weak import ClassicalWeakJacobiForms_factory
            sage: fcty = ClassicalWeakJacobiForms_factory(2, 5)
            sage: fcty._power_series_ring()
            Power Series Ring in q over Rational Field
        """
        return self.__power_series_ring

    def _power_series_ring_ZZ(self):
        r"""
        An auxiliary power series ring over `\ZZ` that is cached in the factory.

        OUTPUT:

        An power series ring.

        EXAMPLES::

            sage: from sage.modular.jacobi.classical_weak import ClassicalWeakJacobiForms_factory
            sage: fcty = ClassicalWeakJacobiForms_factory(2, 5)
            sage: fcty._power_series_ring_ZZ()
            Power Series Ring in q over Integer Ring
        """
        return self.__power_series_ring_ZZ

    def _truncate_cache(self, prec):
        r"""
        A new factory instance whose cache is truncated to a given precision.

        INPUT:

        - ``prec`` -- An integer.

        OUTPUT:

        A factory for classical weak Jacobi forms.

        EXAMPLES::

            sage: from sage.modular.jacobi.classical_weak import ClassicalWeakJacobiForms_factory
            sage: fcty = ClassicalWeakJacobiForms_factory(2, 5)
            sage: fcty_tr = fcty._truncate_cache(3)
            sage: assert isinstance(fcty_tr, ClassicalWeakJacobiForms_factory)
        """
        other = ClassicalWeakJacobiFormsFactory(self.jacobi_index(), prec)

        try:
            other._wronskian_adjoint_even = \
                matrix([[e.truncate(prec) for e in r]
                        for r in self._wronskian_adjoint_even.rows()])
        except AttributeError:
            pass

        try:
            other._wronskian_adjoint_odd = matrix([[e.truncate(prec) for e in r]
                                                   for r in self._wronskian_adjoint_odd.rows()])
        except AttributeError:
            pass

        try:
            other.__wronskian_invdeterminant_even = \
                self.__wronskian_invdeterminant_even.truncate(prec)
        except AttributeError:
            pass

        try:
            other.__wronskian_invdeterminant_odd = \
                self.__wronskian_invdeterminant_odd.truncate(prec)
        except AttributeError:
            pass

        return other

    def _set_wronskian_adjoint(self, wronskian_adjoint, weight_parity=0):
        r"""
        Set the cache for the adjoint of the wronskian. See _wronskian_adjoint.

        INPUT:

        - ``wronskian_adjoint`` -- A list of lists of power series over `\ZZ`.

        - ``weight_parity`` -- An integer (default: `0`).

        TESTS::

            sage: from sage.modular.jacobi.classical_weak import ClassicalWeakJacobiForms_factory
            sage: fcty = ClassicalWeakJacobiForms_factory(2, 5)
            sage: R = fcty._power_series_ring_ZZ()
            sage: adj = [[R(1), R(0)], [R(7), R(-2)]]
            sage: fcty._set_wronskian_adjoint(adj, 0)
            sage: assert all(all(e is f for (e,f) in zip(*rr)) for rr in zip(fcty._wronskian_adjoint(0),adj))

        """
        wronskian_adjoint = [ [ e if e in self.__power_series_ring_ZZ else self.__power_series_ring_ZZ(e)
                                for e in row ]
                              for row in wronskian_adjoint ]

        if weight_parity % 2 == 0:
            self._wronskian_adjoint_even = wronskian_adjoint
        else:
            self._wronskian_adjoint_odd = wronskian_adjoint

    def _wronskian_adjoint(self, weight_parity=0):
        r"""
        The matrix `W^\# \pmod{p}`, mentioned on page 142 of [Sk84]_.
        This matrix is represented by a list of lists of q-expansions.

        The q-expansion is shifted by `-(m + 1) (2*m + 1) / 24` in the
        case of even weights, and it is shifted by `-(m - 1) (2*m - 3)
        / 24` otherwise. This is compensated by the missing q-powers
        returned by _wronskian_invdeterminant.

        INPUT:

        - ``weight_parity`` -- An integer (default: `0`).

        OUTPUT:

        A matrix over a power series ring.

        EXAMPLES::

            sage: from sage.modular.jacobi.classical_weak import ClassicalWeakJacobiForms_factory
            sage: fcty = ClassicalWeakJacobiForms_factory(2, 5)
            sage: fcty._wronskian_adjoint(0)
            [[48 - 720*q - 8400*q^3 + 5040*q^4 + O(q^5),
              -60 + 260*q + 2436*q^3 - 5180*q^4 + O(q^5),
              12 - 20*q - 84*q^3 + 140*q^4 + O(q^5)],
             [3072*q^2 + O(q^5),
              32 - 960*q^2 + 2592*q^4 + O(q^5),
              -8 + 48*q^2 - 72*q^4 + O(q^5)],
             [-960*q^2 - 4032*q^3 + O(q^5),
              -2 - 162*q + 1020*q^2 - 550*q^3 + O(q^5),
              2 + 18*q - 60*q^2 + 22*q^3 + O(q^5)]]
        """
        try:
            if weight_parity % 2 == 0:
                return self._wronskian_adjoint_even
            else:
                return self._wronskian_adjoint_ood

        except AttributeError:
            PS = self.__power_series_ring_ZZ
            m = self.jacobi_index()

            twom = 2 * m

            thetas = dict( ((i, j), dict())
                           for i in xrange(m + 1) for j in xrange(m + 1) )

            ## We want to calculate \hat \theta_{j,l} = sum_r (2 m r + j)**2l q**(m r**2 + j r)
            ## in the case of even weight, and \hat \theta_{j,l} = sum_r (2 m r + j)**(2l + 1) q**(m r**2 + j r),
            ## otherwise.
            for r in xrange(isqrt((self.__prec - 1 + m)//m) + 2):
                for j in (xrange(m + 1) if weight_parity % 2 == 0 else range(1, m)):
                    fact_p = (twom * r + j) ** 2
                    fact_m = (twom * r - j) ** 2
                    if weight_parity % 2 == 0:
                        coeff_p = 2
                        coeff_m = 2
                    else:
                        coeff_p = 2 * (twom * r + j)
                        coeff_m = - 2 * (twom * r - j)

                    for l in (xrange(m + 1) if weight_parity % 2 == 0 else range(1, m)):
                        thetas[(j,l)][m*r**2 + r*j] = coeff_p
                        thetas[(j,l)][m*r**2 - r*j] = coeff_m
                        coeff_p = coeff_p * fact_p
                        coeff_m = coeff_m * fact_m
            if weight_parity % 2 == 0:
                thetas[(0, 0)][0] = 1

            thetas = dict((k, PS(th).add_bigoh(self.__prec))
                           for (k,th) in thetas.iteritems())

            W = matrix( PS, m + 1 if weight_parity % 2 == 0 else (m - 1),
                        [ thetas[(j, l)]
                          for j in (xrange(m + 1) if weight_parity % 2 == 0 else range(1, m))
                          for l in (xrange(m + 1) if weight_parity % 2 == 0 else range(1, m)) ] )


            ## Since the adjoint of matrices with entries in a general ring
            ## is extremely slow for matrices of small size, we hard code the
            ## the cases `m = 2` and `m = 3`.  The expressions are obtained by
            ## computing the adjoint of a matrix with entries `w_{i,j}` in a
            ## polynomial algebra.
            if W.nrows() == 1:
                Wadj = matrix(PS, [[1]])
            elif W.nrows() == 2:
                Wadj = matrix(PS, [[W[1, 1], -W[0, 1]],
                                   [-W[1, 0], W[0, 0]]])

            elif W.nrows() == 3 and self.__prec > 10 ** 5:
                adj00 = W[1, 1] * W[2, 2] - W[2, 1] * W[1, 2]
                adj01 = - W[1,0] * W[2,2] + W[2,0] * W[1,2]
                adj02 = W[1,0] * W[2,1] - W[2,0] * W[1,1]
                adj10 = - W[0,1] * W[2,2] + W[2,1] * W[0,2]
                adj11 = W[0,0] * W[2,2] - W[2,0] * W[0,2]
                adj12 = - W[0,0] * W[2,1] + W[2,0] * W[0,1]
                adj20 = W[0,1] * W[1,2] - W[1,1] * W[0,2]
                adj21 = - W[0,0] * W[1,2] + W[1,0] * W[0,2]
                adj22 = W[0,0] * W[1,1] - W[1,0] * W[0,1]

                Wadj = matrix(PS, [[adj00, adj01, adj02],
                                   [adj10, adj11, adj12],
                                   [adj20, adj21, adj22]])

            elif W.nrows() == 4 and self.__prec > 10 ** 5:
                adj00 = -W[0,2]*W[1,1]*W[2,0] + W[0,1]*W[1,2]*W[2,0] + W[0,2]*W[1,0]*W[2,1] - W[0,0]*W[1,2]*W[2,1] - W[0,1]*W[1,0]*W[2,2] + W[0,0]*W[1,1]*W[2,2]
                adj01 = -W[0,3]*W[1,1]*W[2,0] + W[0,1]*W[1,3]*W[2,0] + W[0,3]*W[1,0]*W[2,1] - W[0,0]*W[1,3]*W[2,1] - W[0,1]*W[1,0]*W[2,3] + W[0,0]*W[1,1]*W[2,3]
                adj02 = -W[0,3]*W[1,2]*W[2,0] + W[0,2]*W[1,3]*W[2,0] + W[0,3]*W[1,0]*W[2,2] - W[0,0]*W[1,3]*W[2,2] - W[0,2]*W[1,0]*W[2,3] + W[0,0]*W[1,2]*W[2,3]
                adj03 = -W[0,3]*W[1,2]*W[2,1] + W[0,2]*W[1,3]*W[2,1] + W[0,3]*W[1,1]*W[2,2] - W[0,1]*W[1,3]*W[2,2] - W[0,2]*W[1,1]*W[2,3] + W[0,1]*W[1,2]*W[2,3]

                adj10 = -W[0,2]*W[1,1]*W[3,0] + W[0,1]*W[1,2]*W[3,0] + W[0,2]*W[1,0]*W[3,1] - W[0,0]*W[1,2]*W[3,1] - W[0,1]*W[1,0]*W[3,2] + W[0,0]*W[1,1]*W[3,2]
                adj11 = -W[0,3]*W[1,1]*W[3,0] + W[0,1]*W[1,3]*W[3,0] + W[0,3]*W[1,0]*W[3,1] - W[0,0]*W[1,3]*W[3,1] - W[0,1]*W[1,0]*W[3,3] + W[0,0]*W[1,1]*W[3,3]
                adj12 = -W[0,3]*W[1,2]*W[3,0] + W[0,2]*W[1,3]*W[3,0] + W[0,3]*W[1,0]*W[3,2] - W[0,0]*W[1,3]*W[3,2] - W[0,2]*W[1,0]*W[3,3] + W[0,0]*W[1,2]*W[3,3]
                adj13 = -W[0,3]*W[1,2]*W[3,1] + W[0,2]*W[1,3]*W[3,1] + W[0,3]*W[1,1]*W[3,2] - W[0,1]*W[1,3]*W[3,2] - W[0,2]*W[1,1]*W[3,3] + W[0,1]*W[1,2]*W[3,3]

                adj20 = -W[0,2]*W[2,1]*W[3,0] + W[0,1]*W[2,2]*W[3,0] + W[0,2]*W[2,0]*W[3,1] - W[0,0]*W[2,2]*W[3,1] - W[0,1]*W[2,0]*W[3,2] + W[0,0]*W[2,1]*W[3,2]
                adj21 = -W[0,3]*W[2,1]*W[3,0] + W[0,1]*W[2,3]*W[3,0] + W[0,3]*W[2,0]*W[3,1] - W[0,0]*W[2,3]*W[3,1] - W[0,1]*W[2,0]*W[3,3] + W[0,0]*W[2,1]*W[3,3]
                adj22 = -W[0,3]*W[2,2]*W[3,0] + W[0,2]*W[2,3]*W[3,0] + W[0,3]*W[2,0]*W[3,2] - W[0,0]*W[2,3]*W[3,2] - W[0,2]*W[2,0]*W[3,3] + W[0,0]*W[2,2]*W[3,3]
                adj23 = -W[0,3]*W[2,2]*W[3,1] + W[0,2]*W[2,3]*W[3,1] + W[0,3]*W[2,1]*W[3,2] - W[0,1]*W[2,3]*W[3,2] - W[0,2]*W[2,1]*W[3,3] + W[0,1]*W[2,2]*W[3,3]

                adj30 = -W[1,2]*W[2,1]*W[3,0] + W[1,1]*W[2,2]*W[3,0] + W[1,2]*W[2,0]*W[3,1] - W[1,0]*W[2,2]*W[3,1] - W[1,1]*W[2,0]*W[3,2] + W[1,0]*W[2,1]*W[3,2]
                adj31 = -W[1,3]*W[2,1]*W[3,0] + W[1,1]*W[2,3]*W[3,0] + W[1,3]*W[2,0]*W[3,1] - W[1,0]*W[2,3]*W[3,1] - W[1,1]*W[2,0]*W[3,3] + W[1,0]*W[2,1]*W[3,3]
                adj32 = -W[1,3]*W[2,2]*W[3,0] + W[1,2]*W[2,3]*W[3,0] + W[1,3]*W[2,0]*W[3,2] - W[1,0]*W[2,3]*W[3,2] - W[1,2]*W[2,0]*W[3,3] + W[1,0]*W[2,2]*W[3,3]
                adj33 = -W[1,3]*W[2,2]*W[3,1] + W[1,2]*W[2,3]*W[3,1] + W[1,3]*W[2,1]*W[3,2] - W[1,1]*W[2,3]*W[3,2] - W[1,2]*W[2,1]*W[3,3] + W[1,1]*W[2,2]*W[3,3]

                Wadj = matrix(PS, [ [adj00, adj01, adj02, adj03],
                                    [adj10, adj11, adj12, adj13],
                                    [adj20, adj21, adj22, adj23],
                                    [adj30, adj31, adj32, adj33] ])
            else:
                Wadj = W.adjoint()

            if weight_parity % 2 == 0:
                wronskian_adjoint = [[Wadj[j, r] for j in xrange(m + 1)]
                                      for r in xrange(m + 1)]
            else:
                wronskian_adjoint = [[Wadj[j, r] for j in xrange(m - 1)]
                                      for r in xrange(m - 1)]

            if weight_parity % 2 == 0:
                self._wronskian_adjoint_even = wronskian_adjoint
            else:
                self._wronskian_adjoint_odd = wronskian_adjoint

            return wronskian_adjoint

    def _set_wronskian_invdeterminant(self, wronskian_invdeterminant,
                                      weight_parity=0):
        r"""
        Set the cache for the inverse determinant of the Wronskian.

        See _wronskian_adjoint.

        INPUT:

        - ``wronskian_invdeterminant`` -- A power series over `\ZZ`.

        - ``weight_parity`` -- An integer (default: `0`).

        TESTS::

            sage: from sage.modular.jacobi.classical_weak import ClassicalWeakJacobiForms_factory
            sage: fcty = ClassicalWeakJacobiForms_factory(2, 5)
            sage: invdet = fcty._power_series_ring_ZZ()(13)
            sage: fcty._set_wronskian_invdeterminant(invdet, 0)
            sage: assert fcty._wronskian_invdeterminant(0) is invdet
        """
        if not wronskian_invdeterminant in self.__power_series_ring_ZZ:
            wronskian_invdeterminant = self.__power_series_ring_ZZ(wronskian_invdeterminant)

        if weight_parity % 2 == 0:
            self._wronskian_invdeterminant_even = wronskian_invdeterminant
        else:
            self._wronskian_invdeterminant_odd = wronskian_invdeterminant

    def _wronskian_invdeterminant(self, weight_parity=0):
        r"""
        The inverse determinant of `W`, which in the considered cases
        is always a negative power of the eta function.

        See [Sk84]_ for details.

        INPUT:

        - ``weight_parity`` -- An integer (default: `0`).

        OUTPUT:

        A power series.

        EXAMPLES::

            sage: from sage.modular.jacobi.classical_weak import ClassicalWeakJacobiForms_factory
            sage: fcty = ClassicalWeakJacobiForms_factory(2, 5)
            sage: fcty._wronskian_invdeterminant(0)
            1 + 15*q + 135*q^2 + 920*q^3 + 5220*q^4 + O(q^5)
        """
        try:
            if weight_parity % 2 == 0:
                return self._wronskian_invdeterminant_even
            else:
                return self._wronskian_invdeterminant_odd

        except AttributeError:
            m = self.jacobi_index()
            if weight_parity % 2 == 0:
                pw = (m + 1) * (2 * m + 1)
            else:
                pw = (m - 1) * (2 * m - 1)

            eta = self.__power_series_ring_ZZ([number_of_partitions(n)
                                               for n in xrange(self.__prec)]).add_bigoh(self.__prec)

            wronskian_invdeterminant = eta ** pw

            if weight_parity % 2 == 0:
                self._wronskian_invdeterminant_even = wronskian_invdeterminant
            else:
                self._wronskian_invdeterminant_odd = wronskian_invdeterminant

            return wronskian_invdeterminant

    def from_taylor_expansion(self, fs, k, is_integral=False):
        r"""
        The Fourier expansion of a weak Jacobi form computed from its
        corrected Taylor coefficients.

        ALGORITHM:

        We combine the theta decomposition and the heat operator as in
        [Sk84]_. This yields a bijections of the space of weak Jacobi
        forms of weight `k` and index `m` with the product of spaces
        of elliptic modular forms `M_k \times S_{k+2} \times .. \times
        S_{k+2m}`.

        INPUT:

        - ``fs`` -- A list of functions that given an integer `p` return the
          q-expansion of a modular form with rational coefficients
          up to precision `p`.  These modular forms correspond to
          the components of the above product.

        - `k` -- An integer. The weight of the weak Jacobi form to be computed.

        - ``is_integral`` -- A boolean. If ``True``, the ``fs`` have integral
          coefficients.

        TESTS::

            sage: from sage.modular.jacobi.classical_weak import *
            sage: factory = ClassicalWeakJacobiFormsFactory(3, 10)
            sage: R.<q> = ZZ[[]]
            sage: expansion = factory.from_taylor_expansion([lambda p: 0 + O(q^p), lambda p: CuspForms(1, 12).gen(0).qexp(p)], 9, True)
            sage: exp_gcd = gcd(expansion.values())
            sage: sorted([ (12 * n - r**2, c/exp_gcd) for ((n, r), c) in expansion.items()])
            [(-4, 0), (-1, 0), (8, 1), (11, -2), (20, -14), (23, 32), (32, 72), (35, -210), (44, -112), (47, 672), (56, -378), (59, -728), (68, 1736), (71, -1856), (80, -1008), (83, 6902), (92, -6400), (95, -5792), (104, 10738), (107, -6564)]
        """
        if is_integral:
            PS = self.__power_series_ring_ZZ
        else:
            PS = self.__power_series_ring

        if (k % 2 == 0 and not len(fs) == self.jacobi_index() + 1) \
          or (k % 2 != 0 and not len(fs) == self.jacobi_index() - 1):
            msg = "fs (which has length {0}) must be a list of {1} Fourier expansions"
            raise ValueError(msg.format(len(fs),
                                        self.jacobi_index() + 1 if k % 2 == 0
                                        else self.jacobi_index() - 1))

        if self.__prec <= 0:  # there are no Fourier indices below the precision
            return dict()

        f_divs = dict()
        for (i, f) in enumerate(fs):
            f_divs[(i, 0)] = PS(f(self.__prec), self.__prec + 10)

        m = self.jacobi_index()

        for i in (xrange(m + 1) if k % 2 == 0 else xrange(m - 1)):
            for j in xrange(1, m - i + 1):
                f_divs[(i,j)] = f_divs[(i, j - 1)].derivative().shift(1)

        phi_divs = list()
        for i in (xrange(m + 1) if k % 2 == 0 else xrange(m - 1)):
            if k % 2 == 0:
                ## This is (13) on page 131 of Skoruppa (change of variables n -> i, r -> j).
                ## The additional factor (2m + k - 1)! is a renormalization to make coefficients integral.
                ## The additional factor (4m)^i stems from the fact that we have used d / d\tau instead of
                ## d^2 / dz^2 when treating the theta series.  Since these are annihilated by the heat operator
                ## the additional factor compensates for this.
                phi_divs.append( sum( f_divs[(j, i - j)]
                                      * ( (4 * self.jacobi_index()) ** i
                                          * binomial(i, j) / 2 ** i  # 2**(self.__jacobi_index - i + 1)
                                          * prod(2 * l + 1 for l in xrange(i))
                                          / factorial(i + k + j - 1)
                                          * factorial(2*self.jacobi_index() + k - 1) )
                                      for j in range(i + 1) ) )
            else:
                phi_divs.append( sum( f_divs[(j, i - j)]
                                      * ( (4 * self.jacobi_index()) ** i
                                          * binomial(i, j) / 2 ** (i + 1)  # 2**(self.__jacobi_index - i + 1)
                                          * prod(2*l + 1 for l in xrange(i + 1))
                                          / factorial(i + k + j)
                                          * factorial(2*self.jacobi_index() + k - 1) )
                                      for j in range(i + 1) ) )

        phi_coeffs = dict()
        for r in (xrange(m + 1) if k % 2 == 0 else xrange(1, m)):
            if k % 2 == 0:
                series = sum(map(operator.mul, self._wronskian_adjoint(k)[r], phi_divs))
            else:
                series = sum(map(operator.mul, self._wronskian_adjoint(k)[r - 1], phi_divs))
            series = self._wronskian_invdeterminant(k) * series

            for n in xrange(self.__prec):
                phi_coeffs[(n, r)] = series[n]

        return phi_coeffs
