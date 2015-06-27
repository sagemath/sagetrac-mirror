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
    matrix, diagonal_matrix,
    vector, zero_vector,
    ZZ, QQ,
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
    r"""
    A set of quadratic forms used in subsequent tests.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _test_set__quadratic_forms
        sage: _test_set__quadratic_forms()
        [Quadratic form...]
    """
    return [QuadraticForm(matrix([[6]])),
            QuadraticForm(diagonal_matrix([2,4])),
            QuadraticForm(matrix(3, [2,1,1, 1,2,1, 1,1,2]))]

def _test_set__indefinite_quadratic_forms():
    r"""
    A set of indefinite quadratic forms used in subsequent tests.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _test_set__indefinite_quadratic_forms
        sage: _test_set__indefinite_quadratic_forms()
        [Quadratic form...]
    """
    return [QuadraticForm(matrix([[-2]])),
            QuadraticForm(matrix(2, [-2,-1,-1,-2])),
            QuadraticForm(matrix(2, [-2,0,0,-4])),
            QuadraticForm(matrix(3, [-2,-1,-1, -1,-2,-1, -1,-1,-2])
                          .block_sum(matrix(2, [2,0,0,-4])))]

def test_vector_valued_modular_forms():
    r"""
    Test vector valued modular forms.  See individual tests
    for more details.

    .. NOTE:

    This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import test_vector_valued_modular_forms
        sage: test_vector_valued_modular_forms()
        <generator object ...>
    """
    prec = 5

    for L in _test_set__quadratic_forms():
        for k in [k_jac + L.dim()/ZZ(2) for k_jac in [10,15,21,22]]:

            vforms = vector_valued_modular_forms(k, L, prec)
            assert len(vforms) != 0
            yield (_test_vector_valued_modular_forms__multiplication,
                   k, L, 4, prec, vforms)

def _test_vector_valued_modular_forms__multiplication(k, L, k_mod, prec, vforms):
    r"""
    Test vector valued modular forms by taking products with classical
    modular forms.

    INPUT:

    - `k` -- An integer.

    - `m` -- A qudratic form.

    - ``k_mod`` -- An integer.

    - ``prec`` -- An integer.

    - ``vforms`` -- A list of dictionaries, representing vector valued modular forms.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _test_vector_valued_modular_forms__multiplication
        sage: from sage.modular.jacobi.vector_valued import vector_valued_modular_forms
        sage: vforms = vector_valued_modular_forms(7/2,QuadraticForm(ZZ,1,[1]),2)
        sage: _test_vector_valued_modular_forms__multiplication(7/2,QuadraticForm(ZZ,1,[1]),4,2,vforms)
    """
    L_span = L.matrix().row_module()
    img_forms = vector_valued_modular_forms(k + k_mod, L, prec)

    for f in ModularForms(1, k_mod).basis():
        _multipliy_and_check(prec, f.qexp(prec).dict(), vforms, 0, img_forms, L_span)

def test_vector_valued_modular_forms_weakly_holomorphic():
    r"""
    Test weakly holomorphic vector valued modular forms.

    See individual tests for more details.

    .. NOTE::

        This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import test_vector_valued_modular_forms_weakly_holomorphic
        sage: test_vector_valued_modular_forms_weakly_holomorphic()
        <generator object ...>
    """
    prec = 5

    for L in _test_set__quadratic_forms():
        for k in [k_jac - L.dim()/ZZ(2) for k_jac in [6, 15, 21, 22]]:
            for order in [1, 2, 3]:
                vforms = vector_valued_modular_forms_weakly_holomorphic(k, L, order, prec)
                if len(vforms) == 0:
                    print k, L, order
                    raise AssertionError()

                yield (_test_vector_valued_modular_forms_weakly_holomorphic__order,
                       order, vforms)
                yield (_test_vector_valued_modular_forms_weakly_holomorphic__multiplication,
                       k, L, order, prec, vforms)


def test_vector_valued_modular_forms_weakly_holomorphic_with_principal_part():
    r"""
    Test weakly holomorphic vector valued modular forms with given
    principal part.  See individual tests for more details.

    .. NOTE::

        This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import test_vector_valued_modular_forms_weakly_holomorphic_with_principal_part
        sage: test_vector_valued_modular_forms_weakly_holomorphic_with_principal_part()
        <generator object ...>
    """
    prec = 5
    for (k, L, pp) in [(11/ZZ(2), QuadraticForm(matrix([[4]])), {(0,): {-1: 2}}),
                       (5, QuadraticForm(matrix(2, [2,1,1,2])), {(0,1): {-1/ZZ(3): 1}, (0,2): {-1/ZZ(3): -1}})]:
        wvvform = vector_valued_modular_forms_weakly_holomorphic_with_principal_part(k, L, pp, prec)

        yield (_test_vector_valued_modular_forms_weakly_holomorphic__principal_part,
               wvvform, pp, L)

        yield (_test_vector_valued_modular_forms_weakly_holomorphic__order,
               1, [wvvform])
        yield (_test_vector_valued_modular_forms_weakly_holomorphic__multiplication,
               k, L, 1, prec, [wvvform])

def _test_vector_valued_modular_forms_weakly_holomorphic__principal_part(wvvform, pp, L):
    r"""
    Test weakly holomorphic vector valued modular forms with given
    principal part.

    INPUT:

    - ``wvvform`` -- A dictionary representing a weakly holomorphic modular forms.

    - ``pp`` -- A dictionary representing a principal part.

    - `L` -- A quadratic form.
    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _test_vector_valued_modular_forms_weakly_holomorphic__principal_part
        sage: _test_vector_valued_modular_forms_weakly_holomorphic__principal_part({(0,) : {-1: 1}}, {(0,) : {-1: 1}}, QuadraticForm(ZZ,1,[1]))
    """
    from sage.combinat.dict_addition import dict_addition

    pp_computed = {}
    for (mu,fe) in wvvform.items():
        pp_computed[mu] = {}
        for (n,coeff) in fe.items():
            if n < 0:
                pp_computed[mu][n] = coeff

    L_span = L.matrix().row_module()
    mu_module = L_span.ambient_module() / L_span

    for mu in list(mu_module):
        pp_contribution = {}
        for pp_mu in pp.keys():
            if mu_module(vector(pp_mu)) != mu:
                continue
            pp_contribution = dict_addition([pp_contribution, pp[pp_mu]])

        pp_computed_contribution = {}
        for pp_mu in pp_computed.keys():
            if mu_module(vector(pp_mu)) != mu:
                continue
            pp_computed_contribution = dict_addition([pp_computed_contribution,
                                                      pp_computed[pp_mu]])

        assert pp_contribution == pp_computed_contribution


def _test_vector_valued_modular_forms_weakly_holomorphic__order(order, vforms):
    r"""
    Test the order at infinity of weakly holomorphic vector valued modular forms.

    INPUT:

    - ``order`` -- An integer.

    - ``vforms`` -- A dictionary representing a weakly holomorphic modular forms.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _test_vector_valued_modular_forms_weakly_holomorphic__order
        sage: _test_vector_valued_modular_forms_weakly_holomorphic__order(1, [{(0,) : {1: 1}}])
    """
    for f in vforms:
        for (mu,fe) in f.items():
            assert all(n >= -order for n in fe.keys())

def _test_vector_valued_modular_forms_weakly_holomorphic__multiplication(k, L, order, prec, vforms):
    r"""
    Test weakly holomorphic vector valued modular forms by multiplying
    them with a power of the discriminant modular form.

    INPUT:

    - `k` -- An integer.

    - `L` -- A quadratic form.

    - ``order`` -- An integer.

    - ``vforms`` -- A dictionary representing a weakly holomorphic modular forms.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _test_vector_valued_modular_forms_weakly_holomorphic__multiplication
sage: from sage.modular.jacobi.vector_valued import vector_valued_modular_forms
        sage: vforms = vector_valued_modular_forms(7/2,QuadraticForm(ZZ,1,[1]),2)
        sage: _test_vector_valued_modular_forms_weakly_holomorphic__multiplication(7/2, QuadraticForm(ZZ,1,[1]), 1, 2, vforms)
    """

    L_span = L.matrix().row_module()

    img_forms = vector_valued_modular_forms(k + 12*order, L, prec + order)
    f = (CuspForms(1,12).gen(0).qexp(prec + 2*order)**order).dict()
    _multipliy_and_check(prec, f, vforms, -order, img_forms, L_span)

def _multipliy_and_check(prec, f_factor, vforms, valuation, img_forms, L_span):
    r"""
    Multiply vector valued modular forms with a modular form and check
    the image.

    INPUT:

    - ``prec`` -- An integer.

    - ``f_factor`` -- A dictionary representing a modular form.

    - ``vforms`` -- A list of dictionaries, representing vector valued modular forms.

    - ``img_forms`` -- A list of dictionaries, representing vector valued modular forms.

    - ``L_span`` -- A module over the integers.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _test_vector_valued_modular_forms__multiplication
        sage: from sage.modular.jacobi.vector_valued import vector_valued_modular_forms
        sage: vforms = vector_valued_modular_forms(7/2,QuadraticForm(ZZ,1,[1]),2)
        sage: _test_vector_valued_modular_forms__multiplication(7/2,QuadraticForm(ZZ,1,[1]),4,2,vforms) # indirect test
    """
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
    r"""
    Multiply a vector valued modular form `g` by a classical one `f`.

    INPUT:

    - `f` -- A dictionary representing a classical modular form.

    - `g` -- A dictionary representing a vector valued modular form.

    - ``L_span`` -- A module over the integers.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _mul_simple_vvform
        sage: _mul_simple_vvform({1: 2}, {(0,): {0: 1}}, span([vector([2])]))
        {(0,): {1: 2}}
    """
    res = {}
    for (s,c) in f.items():
        res = _add_vvforms(res, _mul_scalar_vvform(c,_shift_vvform(g, s)), L_span)
    return res

def _shift_vvform(f, s):
    r"""
    Shift the exponents of a vector valued modular form.

    INPUT:

    - `f` -- A dictionary representing a classical modular form.

    - `s` -- An integer.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _shift_vvform
        sage: _shift_vvform({(0,): {0: 1}}, 1)
        {(0,): {1: 1}}
    """
    res = {}
    for (mu, fe) in f.items():
        res[mu] = dict((n+s,c) for (n,c) in fe.items())
    return res

def _vvform_to_vector(mu_list, valuation, prec, f):
    r"""
    Convert a  vector valued modular form into a vector of rationals.

    INPUT:

    - ``mu_list`` -- A list of elements of a finite generate abelian group.

    - ``valuation`` -- An integer.

    - ``prec`` -- An integer.

    - `g` -- A dictionary representing a vector valued modular form.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _vvform_to_vector
        sage: mu_module = ZZ**1 / span([vector([2])])
        sage: _vvform_to_vector(list(mu_module), 0, 1, {(0,): {0: 1}, (1,): {1/4: 3}})
        (1, 0)
    """
    mu_module = mu_list[0].parent()
    cmp_len = prec - valuation

    res = zero_vector(QQ, cmp_len * len(mu_list))
    for (mu, fe) in f.items():
        if len(fe) == 0:
            continue

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
    r"""
    Test implementation to find positive definite quadratic forms.
    See individual tests for more details.

    .. NOTE::

        This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import test_stably_equivalent_positive_definite_quadratic_form
        sage: test_stably_equivalent_positive_definite_quadratic_form()
        <generator object ...>
    """
    for L in _test_set__indefinite_quadratic_forms():
        yield (_test_stably_equivalent_positive_definite_quadratic_form,
               L)

def _test_stably_equivalent_positive_definite_quadratic_form(L):
    r"""
    Test implementation to find positive definite quadratic forms.

    INPUT:

    - `L` -- A quadratic form.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _test_stably_equivalent_positive_definite_quadratic_form
        sage: _test_stably_equivalent_positive_definite_quadratic_form(QuadraticForm(ZZ,1,[1]))
    """
    M = stably_equivalent_positive_definite_quadratic_form(L, True)

    assert M.is_positive_definite()

    Mrowmodule = M.matrix().row_module()
    Lrowmodule = L.matrix().row_module()
    assert ((Mrowmodule.ambient_module() / Mrowmodule).invariants()
            ==
            (Lrowmodule.ambient_module() / Lrowmodule).invariants())

def test__split_off_hyperbolic():
    r"""
    Test splitting off of hyperbolic lattices.  See individual tests
    for more details.

    .. NOTE:

    This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import test__split_off_hyperbolic
        sage: test__split_off_hyperbolic()
        <generator object ...>
    """
    for L in _test_set__indefinite_quadratic_forms():
        yield (_test__split_off_hyperbolic,
               L)

def _test__split_off_hyperbolic(L):
    r"""
    Test splitting off of hyperbolic lattices.

    INPUT:

    - `L` -- A quadratic form.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _test__split_off_hyperbolic
        sage: _test__split_off_hyperbolic(QuadraticForm(ZZ,1,[-1]))
    """
    M = _split_off_hyperbolic(L)

    assert M.dim() + 2 == L.dim() + 8

    Mrowmodule = M.matrix().row_module()
    Lrowmodule = L.matrix().row_module()
    assert ((Mrowmodule.ambient_module() / Mrowmodule).invariants()
            ==
            (Lrowmodule.ambient_module() / Lrowmodule).invariants())


def test__split_off_E8():
    r"""
    Test splitting off of `E_8`.  See individual tests for more
    details.

    .. NOTE:

    This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import test__split_off_E8
        sage: test__split_off_E8()
        <generator object ...>
    """
    E8mat = matrix(ZZ, 8,
                   [2, -1, 0, 0, 0, 0, 0, 0,
                    -1, 2, -1, 0, 0, 0, 0, 0,
                    0, -1, 2, -1, 0, 0, 0, -1,
                    0, 0, -1, 2, -1, 0, 0, 0,
                    0, 0, 0, -1, 2, -1, 0, 0,
                    0, 0, 0, 0, -1, 2, -1, 0,
                    0, 0, 0, 0, 0, -1, 2, 0,
                    0, 0, -1, 0, 0, 0, 0, 2])

    for M in _test_set__quadratic_forms():
        L = QuadraticForm(M.matrix().block_sum(E8mat))
        yield (_test__split_off_E8, L, M)


def _test__split_off_E8(L, M):
    r"""
    Test splitting off `E_8`.

    INPUT:

    - `L` -- A quadratic form, containing `E_8`

    - `M` -- A quadratic form.

    TESTS::

        sage: from sage.modular.jacobi.test_vector_valued import _test__split_off_E8
        sage: E8mat = matrix(ZZ, 8, [2, -1, 0, 0,  0, 0, 0, 0,  -1, 2, -1, 0,  0, 0, 0, 0,  0, -1, 2, -1,  0, 0, 0, -1,  0, 0, -1, 2,  -1, 0, 0, 0,  0, 0, 0, -1,  2, -1, 0, 0,  0, 0, 0, 0,  -1, 2, -1, 0,  0, 0, 0, 0,  0, -1, 2, 0,  0, 0, -1, 0,  0, 0, 0, 2])
        sage: _test__split_off_E8(QuadraticForm(E8mat), QuadraticForm(ZZ,0,[]))
    """
    M_computed = _split_off_E8(L)

    ## when isomorphism test is available, use it
    assert M.theta_series() == M_computed.theta_series()
