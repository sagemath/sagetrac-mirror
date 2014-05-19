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

def _test_set__quadratic_forms():
    return [QuadraticForm(matrix([[14]])),
            QuadraticForm(-matrix(2, [2, 1, 1, 2])),
            QuadraticForm(-matrix(2, [2, 0, 0, 2])),
            QuadraticForm(-matrix(2, [2, 0, 0, 4])),
            QuadraticForm(matrix(3, [2,1,1, 1,2,1, 1,1,2]))]
def test__discrimant_form():
    for L in _test_set__quadratic_forms():
        yield (_test__discriminant_form, (L,))

def _test__discriminant_form(L):
    from sage.modular.jacobi.higherrank_dimension _discriminant_form

    (eds, quad, bil) = _discriminant_form(L)

    assert prod(eds) == L.det()

    l = len(eds)
    for ix in range(l):
        v = zero_vector(ZZ, l)
        for a in range(eds[ix]):
            v[ix] = a

            # the value of the quadratic form is compatible with the
            # element's order
            assert quad(v) * (eds[ix] / gcd(a, eds[ix]))**2 in ZZ:

            # discrimiant forms are non-degenerate
            w = zero_vector(ZZ, l)
            for b in range(eds[ix]):
                w[ix] = b
                if bil(v, w) not in ZZ: break
            else:
                raise AssertionError()

            # the descomposition into jordan components is orthogonal
            for ix2 in range(l):
                if ix == ix2: continue
                w = zero_vector(ZZ, l)
                for b in range(eds[ix2]):
                    w[ix2] = b
                    assert bil(v, w) in ZZ

def test__discriminant_form_pmone():
    for L in _test_set__quadratic_forms():
        yield (_test__discriminant_form_pmone, (L,))

def _test__discriminant_form_pmone(L):
    from sage.modular.jacobi.higherrank_dimension _discriminant_form_pmone

    (eds, quad, bil) = _discriminant_form(L)
    (singls, pairs) = _discriminant_form_pmone(L, eds)

    is_trivial = lambda x: all(x % e == 0 for (x, e) in zip(x, eds))
     
     for x in singls:
         assert is_trivial(map(operator.add, x, x))
     for x in pairs:
         assert any(is_trivial(map(operator.add, x, y))
                    for y in pairs if y != x)

def test__weil_representation():
    for L in _test_set__quadratic_forms():
        yield (_test__weil_representation, (L,))

def _test__weil_representation(L):
    from sage.modular.jacobi.higherrank_dimension _weil_representation

    CC = ComplexField(100)
    (eds, quad, bil) = _discriminant_form(L)
    (singls, pairs) = _discriminant_form_pmone(L, eds)
    plus_basis = True

    for plus_basis in [True, False]:
        (S, T) = _weil_representation(CC, L, singls, pairs, plus_basis,
                                      discriminant_form_exponents,
                                      disc_quadratic, disc_bilinear)

        if plus_basis:
            assert S.nrows() == S.ncols() == len(singls) + len(pairs)
            assert T.nrows() == T.ncols() == len(singls) + len(pairs)
        else:
            assert S.nrows() == S.ncols() == len(pairs)
            assert T.nrows() == T.ncols() == len(pairs)

        assert (T**lcm(eds)).is_one()
        S4 = S**4
        ST3 = (S*T)**3
        if L.dim() % 2 == 0:
            assert S4.is_one()
            assert ST3.is_one()
        else:
            assert S4.is_one()
            assert ST3.is_one()
            assert (S4**2).is_one()
            assert (ST3**2).is_one()
