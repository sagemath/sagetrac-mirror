r"""
Tests for vector_valued.py.

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

def _test_set__quadratic_forms():
    return [QuadraticForm(matrix([[6]])),
            QuadraticForm(diagonal_matrix([2,2,4])),
            QuadraticForm(matrix(4, [2,0,0,1, 0,2,0,1, 0,0,2,1, 1,1,1,2]))]

def test_vector_valued():
    prec = 20
    k_mod = 4

    for L in _test_set__quadratic_forms():
        for k in [kk - L.dim()/2 for kk in [10,15,21,22]]:

            vforms = vector_valued(k, L, prec)
            yield (_test_vector_valued__multiplication,
                   (prec, k, L, k_mod))

def _test_vector_valued__multiplication(prec, k, L, k_mod):
    L_span = L.matrix().row_module()
    img_forms = vector_valued(k + k_mod, L, prec)

    for f in ModularForms(1, k_mod).basis():
        _multipliy_and_check(prec, f.qexp(prec), vforms, img_forms, L_span)

def test_weakly_holomorphic_vector_valued():
    prec = 20
    k_mod = 4

    for L in _test_set__quadratic_forms():
        for k in [kk - L.dim()/2 for kk in [10,15,21,22]]:
            for order in [-1, -3, 2]:

            vforms = weakly_holomorphic_vector_valued_modular_forms(k, L, order, prec)
            yield (_test_weakly_holomorphic_vector_valued__order,
                   (order, vforms))
            yield (_test_weakly_holomorphic_vector_valued__multiplication,
                   (prec, k, L, order, vforms))

def _test_weakly_holomorphic_vector_valued__order(order, vforms):
    for f in vforms:
        for (mu,fe) in f.items():
            assert all(n >= order for n in fe.keys())

def _test_weakly_holomorphic_vector_valued__multiplication(prec, k, L, order, vforms):
    L_span = L.matrix().row_module()

    if order >= 0:
        img_forms = vector_valued(k, L, prec)
        f = PolynomialRing(QQ, 'q')(1)
        _multipliy_and_check(prec, f, vforms, img_forms, L_span)
    else:
        img_forms = vector_valued(k - 12*order, L, prec)
        f = CuspForms(1,12).gen(0).qexp(prec)
        _multipliy_and_check(prec, f, vforms, img_forms, L_span)
    
def _multipliy_and_check(prec, f, vforms, img_forms, L_span):
    R = PolynomialRing(QQ, 'q')
    q = R.gen(0)

    mu_module = L_span.ambient_module() / L_span
    mu_list = mu_module.list()

    psi_matrix = zero_matrix(QQ, len(mu_list)*prec, len(img_forms))
    for (psi_ix, psi) in enumerate(img_forms):
        for (mu,fe) in psi.items():
            mu_ix = mu_list.index(mu_module()(mu))

            if len(fe) = 0: continue

            n_shift = min(fe.keys())
            n_shift_red = (n_shift.numerator() % n_shift.denominator()) / n_shift.denominator()

            for n in range(n_shift - n_shift_red,
                           prec if n_shift_red == 0 else prec - 1):
                psi_matrix[psi_ix, mu_ix*prec + n] = fe[n + n_shift_red]
    psi_span = psi_matrix.row_module()

    for phi1 in vforms:
        f_poly = f.polynomial()

        phi2_vec = zero_vector(len(mu_list) * prec)
        for (mu,fe) in phi1.items():
            mu_ix = mu_list.index(mu_module()(mu))

            if len(fe) = 0: continue

            n_shift = min(fe.keys())
            n_shift_red = (n_shift.numerator() % n_shift.denominator()) / n_shift.denominator()

            phi1_poly = sum(c * q**(n - n_shift) for (n,c) in fe.items())
            phi2_poly = f_poly * phi1_poly

            for n in range(n_shift - n_shift_red,
                           prec if n_shift_red == 0 else prec - 1):
                phi2_vec[mu_ix*prec + n] = phi2_poly[n + n_shift_red - n_shift]
        assert phi2_vec in psi_span
