r"""
Tests for higherrank.py.

AUTHOR:

- Martin Raum
"""

#===============================================================================
# 
# Copyright (C) 2012-2014 Martin Raum
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

from sage.modules.all import FreeModule
from sage.quadratic_forms.all import QuadraticForm
from sage.rings.all import ZZ

from random import *

from sage.modular.jacobi.higherrank import (
    _higherrank_jacobi_r_classes
)


def _test_set__jacobi_m():
    return [QuadraticForm(matrix(1, [[2]])),
            QuadraticForm(matrix(2, [2, 0, 0, 2])),
            QuadraticForm(matrix(2, [2, 1, 1, 2])),
            QuadraticForm(matrix(3, [2,1,1, 1,2,1, 1,1,2]))]


def test__reduce_higher_rank_jacobi_fe_index():
     for m in _test_set__jacobi_m():
        r_classes = _higherrank_jacobi_r_classes(m)[0]
        m_span = m.matrix().row_module()
        m_adj = QuadraticForm(2 * m.matrix().adjoint())
        fm = FreeModule(ZZ, m.dim())

        for _ in range(10):
            n = random_int()
            for _ in range(10):
                r = fm.random_element()
                yield (_test__reduce_higher_rank_jacobi_fe_index,
                       ((n, r), m, m_adj, r_classes, m_span))

def _test__reduce_higher_rank_jacobi_fe_index((n, r), m, m_adj, r_classes, m_span):
    from sage.modular.jacobi.higherrank import reduce_higher_rank_jacobi_fe_index

    ((nred, rred), s) = reduce_higherrank_jacobi_fe_index((n, r), m_adj, r_classes, m_span)

    if s == 1:
        assert any(vector(r_class[0]) - rred in m_span
                   for r_class in r_classes)
    else:
        assert any(vector(r_class[0]) + rred in m_span
                   for r_class in r_classes)

    assert 2*m.det()*n - m_adj(r) == 2*m.det()*nred - m_adj(rred)

def test__reduce_higher_rank_jacobi_fe_index__r():
    for m in _test_set__jacobi_m():
        r_classes = _higherrank_jacobi_r_classes(m)[0]
        m_span = m.matrix().row_module()
        fm = FreeModule(ZZ, m.dim())

        r = fm.random_element()
        yield (_test__reduce_higher_rank_jacobi_fe_index__r,
               (r, r_classes, m_span))

def test__reduce_higher_rank_jacobi_fe_index__r(r, r_classes, m_span):
    from sage.modular.jacobi.higherrank import _reduce_higher_rank_jacobi_fe_index__r

    (rred, s) = _reduce_higherrank_jacobi_fe_index__r(r, r_classes, m_span)

    if s == 1:
        assert any(vector(r_class[0]) - rred in m_span
                   for r_class in r_classes)
    else:
        assert any(vector(r_class[0]) + rred in m_span
                   for r_class in r_classes)

def test__higherrank_jacobi_r_classes():
    for m in _test_set__jacobi_m():
        yield (_test__higherrank_jacobi_r_classes, (m,))

def _test__higherrank_jacobi_r_classes(m):
    from sage.modular.jacobi.higherrank import _higherrank_jacobi_r_classes

    (r_classes, r_classes_reduction_signs) = _higherrank_jacobi_r_classes(m)
    m_span = m.matrix().span()

    ## cardinality test
    assert (m.det()
            ==
            sum((2 if -1 in signs else 1) for signs in r_classes_reduction_signs)
    )

    ## equivalence in classes
    for (r_class, r_signs) in zip(r_classes, r_classes_reduction_signs):
        r0 = vector(r_class[0])
        for (r,s) in zip(map(vector, r_class), r_signs):
            assert (r - s*r0
                    in
                    m_span
            )

    ## distinction
    for ix in range(len(r_classes)):
        for r_class in r_classes[ix+1:]:
            assert (vector(r_classes[ix][0]) - vector(r_class[0]) not in m_span
                    and
                    vector(r_classes[ix][0]) + vector(r_class[0]) not in m_span
            )

def test_higherrank_jacobi_forms():
    for m in _test_set__jacobi_m():
        for k in [8, 11, 19, 20]:
            prec = 20
            jforms = higherrank_jacobi_forms(k, m, prec)
            yield (_test_higherrank_jacobi_forms__restriction,
                   (k, m, prec, jforms))
            yield (_test_higherrank_jacobi_forms__multiplication,
                   (k, m, prec, jforms))

def _test_higherrank_jacobi_forms__restriction(k, m, prec, jforms):
    rand = Random()
    indices_red = [(nr, reduce_higherrank_jacobi_fe_index(nr, m_adj, r_classes, m_span))
                   for nr in higherrank_jacobi_fe_indices(m, prec, reduced=False)]
    

    max_rst_index = 30
    all_relation_rst_vectors = \
        flatten( m.short_vector_list_up_to_length(max_rst_index+1), max_level = 1 )

    relation_rst_vectors = []
    while len(relation_rst_vectors) < 20
        s = rand.choice(all_relation_rst_vectors)
        if s not in relation_rst_vectors: relation_rst_vectors.append(s)

    for s in relation_rst_vectors:
        jforms_restr = classical_jacobi_forms(k, m(s), prec)
        indices_restr = list(classical_weak_jacobi_fe_indices(m(s), prec, reduced=True))
        restr_span = matrix([[(phi[nr] if nr in phi else 0) for nr in indices_restr]
                             for phi in jforms_restr]).row_module()

        for phi in jforms:
            phi_rst = _restrict_jacobi_form(k, phi, s, indices_red, indices_restr)
            assert (vector([phi_rst[nr] if nr in phi_rst else 0]
                           for nr in indices_restr)
                    in restr_span)

def _restrict_jacobi_form(k, phi, s, indices_red, indices_rest):
    phi_rst = dict()

    for ((n,r), ((nred, rred), s)) in indices_red:
        r_rst = s.dot_product(vector(r))
        if (n, r_rst) in indices_rest:
            try:
                phi_rst[(n, r_rst)] += s**k * phi[nrred]
            except KeyError:
                phi_rst[(n, r_rst)] = s**k * phi[nrred]

    return phi_rst

def _test_higherrank_jacobi_forms__multiplication(k, m, k_mod, prec, jforms):
    l = m.dim()
    (r_classes, _) = _higherrank_jacobi_r_classes(m)
    m_span = m.matrix().row_module()

    P = PolynomialRing(LaurentPolynomialRing(QQ, ['zeta{}'.format(ix) for ix in range(l)]), 'q')
    q = P.gen(0)
    zeta = lambda ix: P.base_ring().gen(ix - 1)

    indices = list(higherrank_jacobi_fe_indices(m, prec, reduced=True))
    indices_nonreduced = list(higherrank_jacobi_fe_indices(m, prec, reduced=False))
    red = lambda nr: reduce_higherrank_jacobi_fe_index(nr, m, r_classes, m_span)

    psis = [dict([ (nr, (_psi[nr] if nr in _psi else 0))
                   for nr in indices]) for _psi in higherrank_jacobi_forms(k+k_mod, m, prec)]
    psi_span = matrix([[psi[nr] for nr in indices]
                       for psi in psis]).row_module()

    for phi1 in jforms:
        phi1_poly = P.zero()
        for nr in indices_nonreduced:
            (nrred, s) = red(nr)
            if nrred in phi1:
                phi1_poly += s**k * phi1[nrred] * q**n * prod(zeta(ix)**r[ix] for ix in range(l)]

        for f in ModularForms(1, k_mod).basis():
            f_poly = f.qexp(prec).polynomial()

            phi2_poly = f_poly * phi1_poly
            phi2_vec = vector([phi2_poly[n][r] for (n,r) in indices])
            assert not phi2_vec.is_zero()
            assert phi2_vec in psi_span
