r"""
Tests for classical.py.

AUTHOR:

- Martin Raum
"""

#===============================================================================
#
# Copyright (C) 2014 Martin Raum
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

from sage.all import (
    PolynomialRing, LaurentPolynomialRing,
    QQ,
    vector,
    matrix,
    ModularForms
)

from sage.modular.jacobi.classical import (
    classical_jacobi_forms,
    classical_jacobi_fe_indices,
    classical_jacobi_reduce_fe_index
)

def test_classical_jacobi_forms():
    r"""
    Test classical Jacobi forms for correctness.  See individual tests
    for more details.

    .. NOTE:

    This is a test generator to be used by nosetest.

    TESTS::

       sage: from sage.modular.jacobi.test_classical import test_classical_jacobi_forms
       sage: test_classical_jacobi_forms()
       <generator object ...>
    """
    prec = 10
    k_mod = 6

    for k in [8, 7, 12]:
        for m in range(1, 5):
            yield (_test_classical_jacobi_forms,
                   k, m, k_mod, prec)

def _test_classical_jacobi_forms(k, m, k_mod, prec):
    r"""
    Check that the product of Jacobi forms of given weight and
    index with modular forms of given weight is a Jacobi form.

    INPUT:

    - `k` -- An integer.

    - `m -- A positive integer.

    - `k_mod` -- An integer.

    - ``prec`` -- An integer.

    .. NOTE:

    This is a test to be used by nosetest.

    EXAMPLES::

       sage: from sage.modular.jacobi.test_classical import _test_classical_jacobi_forms
       sage: _test_classical_jacobi_forms(4,1,4,3)
    """
    P = PolynomialRing(LaurentPolynomialRing(QQ, 'zeta'), 'q')
    q = P.gen(0)
    zeta = P.base_ring().gen(0)

    indices = list(classical_jacobi_fe_indices(m, prec, reduced=True))
    indices_nonreduced = list(classical_jacobi_fe_indices(m, prec, reduced=False))
    red = lambda nr: classical_jacobi_reduce_fe_index(nr,m)

    psis = [dict([ (nr, (_psi[nr] if nr in _psi else 0)) for nr in indices])
            for _psi in classical_jacobi_forms(k+k_mod, m, prec)]
    psi_span = matrix([[psi[nr] for nr in indices]
                       for psi in psis]).row_module()

    for phi1 in classical_jacobi_forms(k, m, prec):
        phi1_poly = P.zero()
        for nr in indices_nonreduced:
            (n,r) = nr
            (nrred, s) = red(nr)
            if nrred in phi1:
                phi1_poly += s**k * phi1[nrred] * q**n * zeta**r

        for f in ModularForms(1, k_mod).basis():
            f_poly = f.qexp(prec).polynomial()

            phi2_poly = f_poly * phi1_poly
            phi2_vec = vector([phi2_poly[nr[0]][nr[1]] for nr in indices])

            assert not phi2_vec.is_zero()
            assert phi2_vec in psi_span
