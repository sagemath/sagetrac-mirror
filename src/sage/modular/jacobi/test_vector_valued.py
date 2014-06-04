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

from sage.all import (
    matrix, diagonal_matrix, zero_matrix,
    vector, zero_vector,
    ZZ, QQ, PolynomialRing,
    QuadraticForm,
    ModularForms, CuspForms
)

from sage.modular.jacobi.vector_valued import (
    vector_valued_modular_forms,
    vector_valued_modular_forms_weakly_holomorphic,
    vector_valued_modular_forms_weakly_holomorphic_with_principal_part,
    stably_equivalent_positive_definite_quadratic_form,
    _add_vvforms, _mul_scalar_vvform,
    _split_off_hyperbolic, _split_off_E8
)

def _test_set__quadratic_forms():
    return [QuadraticForm(matrix([[6]])),
            QuadraticForm(diagonal_matrix([2,4])),
            QuadraticForm(matrix(3, [2,1,1, 1,2,1, 1,1,2]))]

def _test_set__indefinite_quadratic_forms():
    return [QuadraticForm(matrix([[-2]])),
            QuadraticForm(matrix(2, [-2,-1,-1,-2])),
            QuadraticForm(matrix(2, [-2,0,0,-4])),
            QuadraticForm(matrix(3, [-2,-1,-1, -1,-2,-1, -1,-1,-2])
                          .block_sum(matrix(2, [2,0,0,-4])))]

def test_vector_valued():
    prec = 5

    for L in _test_set__quadratic_forms():
        for k in [k_jac + L.dim()/ZZ(2) for k_jac in [10,15,21,22]]:

            vforms = vector_valued_modular_forms(k, L, prec)
            assert len(vforms) != 0
            yield (_test_vector_valued__multiplication,
                   prec, k, L, 4, vforms)

def _test_vector_valued__multiplication(prec, k, L, k_mod, vforms):
    L_span = L.matrix().row_module()
    img_forms = vector_valued_modular_forms(k + k_mod, L, prec)

    for f in ModularForms(1, k_mod).basis():
        _multipliy_and_check(prec, f.qexp(prec).dict(), vforms, 0, img_forms, L_span)
    
def test_vector_valued_modular_forms_weakly_holomorphic():
    prec = 5
    k_mod = 4

    for L in _test_set__quadratic_forms():
        for k in [k_jac - L.dim()/ZZ(2) for k_jac in [6,15,21,22]]:
            for order in [1,2,3]:
                vforms = vector_valued_modular_forms_weakly_holomorphic(k, L, order, prec)
                if len(vforms) == 0:
                    print k, L, order
                    raise AssertionError()

                yield (_test_vector_valued_modular_forms_weakly_holomorphic__order,
                       order, vforms)
                yield (_test_vector_valued_modular_forms_weakly_holomorphic__multiplication,
                       prec, k, L, order, vforms)

def test_vector_valued_modular_forms_weakly_holomorphic_with_principal_part():
    prec = 5
    for (k, L, pp) in [(11/ZZ(2), QuadraticForm(matrix([[4]])), {(0,): {-1: 2}}),
                       (5, QuadraticForm(matrix(2, [2,1,1,2])), {(0,1): {-1/ZZ(3): 1}, (0,2): {-1/ZZ(3): -1}})]:
        wvvform = vector_valued_modular_forms_weakly_holomorphic_with_principal_part(k, L, pp, prec)

        yield (_test_vector_valued_modular_forms_weakly_holomorphic__principal_part,
               wvvform, pp, L)

        yield (_test_vector_valued_modular_forms_weakly_holomorphic__order,
               1, [wvvform])
        yield (_test_vector_valued_modular_forms_weakly_holomorphic__multiplication,
               prec, k, L, 1, [wvvform])

def _test_vector_valued_modular_forms_weakly_holomorphic__principal_part(wvvform, pp, L):
    from sage.combinat.dict_addition import dict_addition

    pp_computed = {}
    for (mu,fe) in wvvform.items():
        pp_computed[mu] = {}
        for (n,coeff) in fe.items():
            if n < 0: pp_computed[mu][n] = coeff

    L_span = L.matrix().row_module()
    mu_module = L_span.ambient_module() / L_span

    for mu in list(mu_module):
        pp_contribution = {}
        for pp_mu in pp.keys():
            if mu_module(vector(pp_mu)) != mu: continue
            pp_contribution = dict_addition([pp_contribution, pp[pp_mu]])

        pp_computed_contribution = {}
        for pp_mu in pp_computed.keys():
            if mu_module(vector(pp_mu)) != mu: continue
            pp_computed_contribution = dict_addition([pp_computed_contribution,
                                                      pp_computed[pp_mu]])

        assert pp_contribution == pp_computed_contribution

def _test_vector_valued_modular_forms_weakly_holomorphic__order(order, vforms):
    for f in vforms:
        for (mu,fe) in f.items():
            assert all(n >= -order for n in fe.keys())

def _test_vector_valued_modular_forms_weakly_holomorphic__multiplication(prec, k, L, order, vforms):
    L_span = L.matrix().row_module()

    img_forms = vector_valued_modular_forms(k + 12*order, L, prec + order)
    f = (CuspForms(1,12).gen(0).qexp(prec + 2*order)**order).dict()
    _multipliy_and_check(prec, f, vforms, -order, img_forms, L_span)
    
def _multipliy_and_check(prec, f_factor, vforms, valuation, img_forms, L_span):
    ## NOTE: we assume that img_forms consists of expansions that are
    ## regular at infinity
    assert valuation <= 0 and valuation in ZZ

    mu_module = L_span.ambient_module() / L_span
    mu_list = list(mu_module)

    img_prec = prec - valuation
    img_matrix = matrix([_vvform_to_vector(mu_list, 0, img_prec, f)
                         for f in img_forms])
    img_span = img_matrix.row_module()

    for g in vforms:
        fg = _mul_simple_vvform(f_factor, g, L_span)
        fg_vec = _vvform_to_vector(mu_list, 0, img_prec, fg)

        assert not fg_vec.is_zero()
        assert fg_vec in img_span

def _mul_simple_vvform(f, g, L_span):
    res = {}
    for (s,c) in f.items():
        res = _add_vvforms(res, _mul_scalar_vvform(c,_shift_vvform(g, s)), L_span)
    return res

def _shift_vvform(f, s):
    res = {}
    for (mu, fe) in f.items():
        res[mu] = dict((n+s,c) for (n,c) in fe.items())
    return res

def _vvform_to_vector(mu_list, valuation, prec, f):
    mu_module = mu_list[0].parent()
    cmp_len = prec - valuation

    res = zero_vector(QQ, cmp_len*len(mu_list))
    for (mu,fe) in f.items():
        if len(fe) == 0: continue

        mu_ix = mu_list.index(mu_module(vector(mu)))

        n_shift = min(fe.keys())
        n_shift_red = (n_shift.numerator() % n_shift.denominator()) / n_shift.denominator()

        for n in range(cmp_len if n_shift_red == 0 else cmp_len - 1):
            try:
                res[mu_ix*cmp_len + n] = fe[n + n_shift_red + valuation]
            except KeyError:
                pass

    return res

def test_stably_equivalent_positive_definite_quadratic_form():
    for L in _test_set__indefinite_quadratic_forms():
        yield (_test_stably_equivalent_positive_definite_quadratic_form,
               L)

def _test_stably_equivalent_positive_definite_quadratic_form(L):
    M = stably_equivalent_positive_definite_quadratic_form(L, True)

    assert M.is_positive_definite()

    Mrowmodule = M.matrix().row_module()
    Lrowmodule = L.matrix().row_module()
    assert ((Mrowmodule.ambient_module() / Mrowmodule).invariants()
            ==
            (Lrowmodule.ambient_module() / Lrowmodule).invariants())

def test__split_off_hyperbolic():
    for L in _test_set__indefinite_quadratic_forms():
        yield (_test__split_off_hyperbolic,
               L)

def _test__split_off_hyperbolic(L):
    M = _split_off_hyperbolic(L)

    assert M.dim() + 2 == L.dim() + 8

    Mrowmodule = M.matrix().row_module()
    Lrowmodule = L.matrix().row_module()
    assert ((Mrowmodule.ambient_module() / Mrowmodule).invariants()
            ==
            (Lrowmodule.ambient_module() / Lrowmodule).invariants())


def test__split_off_E8():
    E8mat = matrix(ZZ, 8,
            [2, -1, 0, 0,  0, 0, 0, 0,
            -1, 2, -1, 0,  0, 0, 0, 0,
            0, -1, 2, -1,  0, 0, 0, -1,
            0, 0, -1, 2,  -1, 0, 0, 0,
            0, 0, 0, -1,  2, -1, 0, 0,
            0, 0, 0, 0,  -1, 2, -1, 0,
            0, 0, 0, 0,  0, -1, 2, 0,
            0, 0, -1, 0,  0, 0, 0, 2])

    for M in _test_set__quadratic_forms():
        L = QuadraticForm(M.matrix().block_sum(E8mat))
        yield (_test__split_off_E8,
               L, M)

def _test__split_off_E8(L, M):
    M_computed = _split_off_E8(L)

    ## when isomorphism test is available, use it
    assert M.theta_series() == M_computed.theta_series()
