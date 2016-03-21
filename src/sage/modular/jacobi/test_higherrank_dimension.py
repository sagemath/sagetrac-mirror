r"""
Tests for higherrank_dimension.py.

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

from sage.all import (prod, matrix,
                      ZZ, ComplexField,
                      gcd, lcm,
                      QuadraticForm)
import operator

def _test_set__quadratic_forms():
    r"""
    A set of quadratic forms used in subsequent tests.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank_dimension import _test_set__quadratic_forms
        sage: _test_set__quadratic_forms()
        [Quadratic form...]
    """
    return [QuadraticForm(matrix([[14]])),
            QuadraticForm(-matrix(2, [2, 1, 1, 2])),
            QuadraticForm(-matrix(2, [2, 0, 0, 2])),
            QuadraticForm(-matrix(2, [2, 0, 0, 4])),
            QuadraticForm(matrix(3, [2,1,1, 1,2,1, 1,1,2]))]

def test__discriminant_form():
    r"""
    Test discriminant forms.  See individual tests for more details.

    .. NOTE:

    This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank_dimension import test__discriminant_form
        sage: test__discriminant_form()
        <generator object ...>
    """
    for L in _test_set__quadratic_forms():
        yield (_test__discriminant_form, L)

def _test__discriminant_form(L):
    r"""
    Test discriminant forms.  See individual tests for more details.

    INPUT:

    - `L` -- A quadratic form.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank_dimension import _test__discriminant_form
        sage: _test__discriminant_form(QuadraticForm(ZZ,1,[1]))
    """
    from sage.modular.jacobi.higherrank_dimension import _discriminant_form

    (eds, quad, bil) = _discriminant_form(L)

    assert prod(eds) == L.det()

    l = len(eds)
    for ix in range(l):
        v = l*[0]
        for a in range(eds[ix]):
            v[ix] = a

            # the value of the quadratic form is compatible with the
            # element's order
            assert quad(*v) * (eds[ix] / gcd(a, eds[ix])) ** 2 in ZZ

            # discrimiant forms are non-degenerate
            if not all(e == 0 for e in v):
                w = l*[0]
                for b in range(eds[ix]):
                    w[ix] = b
                    if bil(*(v + w)) not in ZZ:
                        break
                else:
                    raise AssertionError()

            # the descomposition into jordan components is orthogonal
            for ix2 in range(l):
                if ix == ix2:
                    continue
                w = l * [0]
                for b in range(eds[ix2]):
                    w[ix2] = b
                    assert bil(*(v + w)) in ZZ


def test__discriminant_form_pmone():
    r"""
    Test meth:`discriminant_form_pmone`.  See individual tests for more details.

    .. NOTE:

    This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank_dimension import test__discriminant_form_pmone
        sage: test__discriminant_form_pmone()
        <generator object ...>
    """
    for L in _test_set__quadratic_forms():
        yield (_test__discriminant_form_pmone, L)

def _test__discriminant_form_pmone(L):
    r"""
    Test meth:`discriminant_form_pmone`.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank_dimension import _test__discriminant_form_pmone
        sage: _test__discriminant_form_pmone(QuadraticForm(ZZ,1,[1]))
    """
    from sage.modular.jacobi.higherrank_dimension import (_discriminant_form,
                                                          _discriminant_form_pmone)

    (eds, quad, bil) = _discriminant_form(L)
    (singls, pairs) = _discriminant_form_pmone(L, eds)

    is_trivial = lambda x: all(x % e == 0 for (x, e) in zip(x, eds))

    for x in singls:
        assert is_trivial(map(operator.add, x, x))
    for x in pairs:
        assert all(not is_trivial(map(operator.add, x, y))
                   for y in singls+pairs
                   if y != x)

def test__weil_representation():
    r"""
    Test Weil representation.  See individual tests for more details.

    .. NOTE:

    This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank_dimension import test__weil_representation
        sage: test__weil_representation()
        <generator object ...>
    """
    for L in _test_set__quadratic_forms():
        yield (_test__weil_representation, L)

def _test__weil_representation(L):
    r"""
    Test Weil representation.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank_dimension import _test__weil_representation
        sage: _test__weil_representation(QuadraticForm(ZZ,1,[1]))
    """
    from sage.modular.jacobi.higherrank_dimension import (_discriminant_form,
                                                          _discriminant_form_pmone,
                                                          _weil_representation)

    CC = ComplexField(200)
    (eds, quad, bil) = _discriminant_form(L)
    (singls, pairs) = _discriminant_form_pmone(L, eds)
    plus_basis = True

    for plus_basis in [True, False]:
        (S, T) = _weil_representation(CC, L, singls, pairs, plus_basis,
                                      eds, quad, bil)

        if plus_basis:
            assert S.nrows() == S.ncols() == len(singls) + len(pairs)
            assert T.nrows() == T.ncols() == len(singls) + len(pairs)
        else:
            assert S.nrows() == S.ncols() == len(pairs)
            assert T.nrows() == T.ncols() == len(pairs)

        zero_test = lambda m: all(abs(e) < 10 ** -10 for e in m.list())

        assert zero_test(T ** (2 * lcm(eds)) - 1)

        if L.dim() % 2 == 0:
            assert zero_test(S ** 4 - 1)
            assert zero_test((S * T) ** 6 - 1)
        else:
            assert zero_test(S ** 8 - 1)
            assert zero_test((S * T) ** 12 - 1)
