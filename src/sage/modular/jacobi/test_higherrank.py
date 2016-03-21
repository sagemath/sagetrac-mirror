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

from sage.all import (Set, ZZ, QQ, PolynomialRing, LaurentPolynomialRing,
                      FreeModule, vector, zero_vector,
                      matrix,
                      QuadraticForm,
                      ModularForms,
                      randint, flatten, prod
)

from sage.modular.jacobi.classical import (classical_weak_jacobi_fe_indices,
                                           classical_jacobi_forms)
from sage.modular.jacobi.higherrank import (
    higherrank_jacobi_reduce_fe_index,
    higherrank_jacobi_fe_indices,
    higherrank_jacobi_forms,

    higherrank_jacobi_r_classes,
    _complete_set_of_restriction_vectors
)


def _test_set__jacobi_m():
    r"""
    A set of quadratic forms used in subsequent tests.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import _test_set__jacobi_m
        sage: _test_set__jacobi_m()
        [Quadratic form...]
    """
    return [QuadraticForm(matrix([[2]])),
            QuadraticForm(matrix([[8]])),
            QuadraticForm(matrix(2, [2, 0, 0, 2])),
            QuadraticForm(matrix(2, [2, 1, 1, 2])),
            QuadraticForm(matrix(3, [2,1,1, 1,2,1, 1,1,2]))]


def test_higherrank_jacobi_reduce_fe_index():
    r"""
    Test reduction of Fourier expansion indices in the case of higher
    rank Jacobi indices.  See individual tests for more details.

    .. NOTE:

    This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import test_higherrank_jacobi_reduce_fe_index
        sage: test_higherrank_jacobi_reduce_fe_index()
        <generator object ...>
    """
    for m in _test_set__jacobi_m():
        r_classes = higherrank_jacobi_r_classes(m)[0]
        m_span = m.matrix().row_module()
        m_adj = QuadraticForm(2 * m.matrix().adjoint())
        fm = FreeModule(ZZ, m.dim())

        for _ in range(10):
            n = randint(0, 20)
            for _ in range(10):
                r = tuple(fm.random_element())
                yield (_test_higherrank_jacobi_reduce_fe_index,
                       (n, r), m, r_classes, m_adj, m_span)

def _test_higherrank_jacobi_reduce_fe_index((n, r), m, r_classes, m_adj, m_span):
    r"""
    Test invariants of the of the reduction of Fourier expansion
    indices.

    INPUT:

    - `(n,r)` -- A pair of an integer `n` and a tuple `r` of integers.

    - `m` -- A quadratic form.

    - ``r_classes`` -- See meth:`sage.modular.jacobi.higherrank.higherrank_jacobi_r_classes`.

    - ``m_adj`` -- A quadratic form.

    - ``m_span`` -- A module over the integers.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import _test_higherrank_jacobi_reduce_fe_index
    sage: _test_higherrank_jacobi_reduce_fe_index((1,(0,)), QuadraticForm(ZZ,1,[1]), [[(0,)], [(1,), (-1,)]], QuadraticForm(ZZ,1,[1]), span([vector([2])]))
    """
    ((nred, rred), s) = higherrank_jacobi_reduce_fe_index((n, r), m,
                                                          r_classes, m_adj,
                                                          m_span)

    if s == 1:
        assert vector(r) - vector(rred) in m_span
    else:
        assert vector(r) + vector(rred) in m_span
    assert any(vector(r_class[0]) - vector(rred) in m_span
               for r_class in r_classes)

    assert 2 * m.det() * n - m_adj(r) == 2 * m.det() * nred - m_adj(rred)


def test__higherrank_jacobi_reduce_fe_index__r():
    r"""
    Test reduction of the second component of Fourier expansion
    indices in the case of higher rank Jacobi indices.  See individual
    tests for more details.

    .. NOTE:

    This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import test__higherrank_jacobi_reduce_fe_index__r
        sage: test__higherrank_jacobi_reduce_fe_index__r()
        <generator object ...>
    """
    for m in _test_set__jacobi_m():
        r_classes = higherrank_jacobi_r_classes(m)[0]
        m_span = m.matrix().row_module()
        fm = FreeModule(ZZ, m.dim())

        r = fm.random_element()
        yield (_test__higherrank_jacobi_reduce_fe_index__r,
               r, r_classes, m_span)

def _test__higherrank_jacobi_reduce_fe_index__r(r, r_classes, m_span):
    r"""
    Test invariants of the of the reduction of Fourier expansion
    indices.

    INPUT:

    - `r` -- A tuple `r` of integers.

    - ``r_classes`` -- See meth:`sage.modular.jacobi.higherrank.higherrank_jacobi_r_classes`.

    - ``m_span`` -- A module over the integers.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import _test__higherrank_jacobi_reduce_fe_index__r
    sage: _test__higherrank_jacobi_reduce_fe_index__r((0,), [[(0,)], [(1,), (-1,)]], span([vector([2])]))
    """
    from sage.modular.jacobi.higherrank import _higherrank_jacobi_reduce_fe_index__r

    (rred, s) = _higherrank_jacobi_reduce_fe_index__r(r, r_classes, m_span)

    if s == 1:
        assert vector(r) - vector(rred) in m_span
    else:
        assert vector(r) + vector(rred) in m_span
    assert any(vector(r_class[0]) - vector(rred) in m_span
               for r_class in r_classes)

def test_higherrank_jacobi_r_classes():
    r"""
    Test computation of equivalence classes of second components of
    Fourier expansion indices in the case of higher rank Jacobi
    indices.  See individual tests for more details.

    .. NOTE:

    This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import test_higherrank_jacobi_r_classes
        sage: test_higherrank_jacobi_r_classes()
        <generator object ...>
    """
    for m in _test_set__jacobi_m():
        yield (_test_higherrank_jacobi_r_classes,
               m)

def _test_higherrank_jacobi_r_classes(m):
    r"""
    Test computation of equivalence classes of second components of
    Fourier expansion indices in the case of higher rank Jacobi
    indices.

    INPUT:

    - `m` -- A quadratic form.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import _test_higherrank_jacobi_r_classes
        sage: _test_higherrank_jacobi_r_classes(QuadraticForm(ZZ,1,[1]))
    """
    from sage.modular.jacobi.higherrank import higherrank_jacobi_r_classes

    (r_classes, r_classes_reduction_signs) = higherrank_jacobi_r_classes(m)
    m_span = m.matrix().row_module()

    ## cardinality test
    assert (m.det()
            ==
            sum(2 if -1 in r_signs else 1 for r_signs in r_classes_reduction_signs)
    )
    assert all(len(r_class)
               ==
               len(r_signs)
               for (r_class, r_signs) in zip(r_classes, r_classes_reduction_signs)
    )

    ## normalization of signs
    assert all(r_signs[0]
               ==
               1
               for r_signs in r_classes_reduction_signs)

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

def test__complete_set_of_restriction_vectors():
    r"""
    Test restriction vectors in the case of higher
    rank Jacobi indices.  See individual tests for more details.

    .. NOTE:

    This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import test__complete_set_of_restriction_vectors
        sage: test__complete_set_of_restriction_vectors()
        <generator object ...>
    """
    for m in _test_set__jacobi_m():
        yield (_test__complete_set_of_restriction_vectors,
               m)

def _test__complete_set_of_restriction_vectors(m):
    r"""
    Test restriction vectors in the case of higher
    rank Jacobi indices.

    INPUT:

    - `m` -- A quadratic form.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import _test__complete_set_of_restriction_vectors
        sage: _test__complete_set_of_restriction_vectors(QuadraticForm(ZZ,1,[1]))
    """
    m_span = m.matrix().row_module()
    (r_classes, r_classes_reduction_signs) = higherrank_jacobi_r_classes(m)
    rst_vectors = _complete_set_of_restriction_vectors(m, r_classes, r_classes_reduction_signs, m_span)

    ## minimality of restriction set
    assert len(rst_vectors) == len(r_classes)

    ## nondegenerate restriction vectors
    assert all(m(s) > 0 for (s,_) in rst_vectors)


    ## compute local restrictions and check rank
    rst_dicts_even = []
    rst_dicts_odd = []

    for s_tpl in Set(tuple(s) for (s,_) in rst_vectors):
        s = vector(s_tpl)
        rst_dict_even = {}
        rst_dict_odd = {}

        for (cl_ix, (r_class, r_signs)) in enumerate(zip(r_classes, r_classes_reduction_signs)):
            for (r, r_sign) in zip(r_class, r_signs):
                r_rst = s.dot_product(vector(r))

                if r_rst not in rst_dict_even:
                    rst_dict_even[r_rst] = zero_vector(len(r_classes))
                if r_rst not in rst_dict_odd:
                    rst_dict_odd[r_rst] = zero_vector(len(r_classes))

                rst_dict_even[r_rst][cl_ix] += 1
                ## we have to exclude the zero vector, since it is
                ## automatically zero, which is not dedubtible from
                ## the r_sign list
                if 2*vector(r) not in m_span:
                    rst_dict_odd[r_rst][cl_ix] += r_sign

        rst_dicts_even.append(rst_dict_even)
        rst_dicts_odd.append(rst_dict_odd)

    rst_matrix_even = matrix([v for rst_dict in rst_dicts_even for v in rst_dict.values()])
    rst_matrix_odd = matrix([v for rst_dict in rst_dicts_odd for v in rst_dict.values()])

    pm_fixed_point_size = len([r_class for r_class in r_classes
                               if any(2*vector(r) in m_span for r in r_class)])

    assert rst_matrix_even.rank() == len(r_classes)
    assert rst_matrix_odd.rank() == len(r_classes) - pm_fixed_point_size

def test_higherrank_jacobi_forms():
    r"""
    Test Jacobi forms with higher rank indices.  See individual tests
    for more details.

    .. NOTE:

    This is a test generator to be used by nosetest.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import test_higherrank_jacobi_forms
        sage: test_higherrank_jacobi_forms()
        <generator object ...>
    """
    for m in _test_set__jacobi_m():
        for k in [8, 11, 19, 20]:
            prec = 5
            jforms = higherrank_jacobi_forms(k, m, prec)
            yield (_test_higherrank_jacobi_forms__restriction,
                   k, m, prec, jforms)
            yield (_test_higherrank_jacobi_forms__multiplication,
                   k, m, 6, prec, jforms)

def _test_higherrank_jacobi_forms__restriction(k, m, prec, jforms):
    r"""
    Test Jacobi forms with higher rank indices by restricting them to
    scalar Jacobi indices.

    INPUT:

    - `k` -- An integer.

    - `m` -- A qudratic form.

    - ``prec`` -- An integer.

    - ``jforms`` -- A list of dictionaries, representing Jacobi forms.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import _test_higherrank_jacobi_forms__restriction
        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_forms
        sage: jforms = higherrank_jacobi_forms(4,QuadraticForm(ZZ,1,[1]),3)
        sage: _test_higherrank_jacobi_forms__restriction(4,QuadraticForm(ZZ,1,[1]),3,jforms)
    """
    r_classes = higherrank_jacobi_r_classes(m)[0]
    m_adj = QuadraticForm(2 * m.matrix().adjoint())
    m_span = m.matrix().row_module()

    indices_red = [(nr, higherrank_jacobi_reduce_fe_index(nr, m, r_classes, m_adj, m_span))
                   for nr in higherrank_jacobi_fe_indices(m, prec, r_classes, reduced=False)]


    jforms_rst_dict = dict((m_rst,classical_jacobi_forms(k, m_rst, prec))
                           for m_rst in [1, 2, 3, 4])
    relation_rst_vectors = (
        flatten(m.short_vector_list_up_to_length(5, True)[1:], max_level=1)
    )

    nmb_nonzero_rsts = len(jforms)*[0]
    for s in relation_rst_vectors:
        jforms_rst = jforms_rst_dict[m(s)]
        indices_rst = list(classical_weak_jacobi_fe_indices(m(s), prec, reduced=True))
        rst_span = matrix([[(phi[nr] if nr in phi else 0) for nr in indices_rst]
                             for phi in jforms_rst]).row_module()

        for (phi_ix, phi) in enumerate(jforms):
            phi_rst = _restrict_jacobi_form(k, phi, s, indices_red, indices_rst)
            phi_rst_vec = vector([phi_rst[nr] if nr in phi_rst else 0
                                  for nr in indices_rst])
            if not phi_rst_vec.is_zero():
                assert phi_rst_vec in rst_span
                nmb_nonzero_rsts[phi_ix] += 1

    ## there is a non trivial restriction for each phi
    assert all(nmb != 0 for nmb in nmb_nonzero_rsts)


def _restrict_jacobi_form(k, phi, s, indices_red, indices_rst):
    """
    Restrict a Jacobi forms along `s`.

    INPUT:

    - `k` -- An integer.

    - ``phi`` -- A dictionary, representing the Fourier expansion of a
                 Jacobi form.

    - `s` -- A tuple of integers.

    - ``indices_red`` -- Reductions of Fourier indices of ``phi``.

    - ``indices_rst`` -- A list of Fourier indices of the restriction.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import _test_higherrank_jacobi_forms__restriction
        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_forms
        sage: jforms = higherrank_jacobi_forms(4,QuadraticForm(ZZ,1,[1]),3)
        sage: _test_higherrank_jacobi_forms__restriction(4,QuadraticForm(ZZ,1,[1]),3,jforms) ## indirect test
    """
    phi_rst = dict()

    for ((n,r), (nrred, sign)) in indices_red:
        r_rst = s.dot_product(vector(r))
        if (n, r_rst) in indices_rst:
            try:
                phi_rst[(n, r_rst)] += sign**k * phi[nrred]
            except KeyError:
                phi_rst[(n, r_rst)] = sign**k * phi[nrred]

    return phi_rst

def _test_higherrank_jacobi_forms__multiplication(k, m, k_mod, prec, jforms):
    r"""
    Test Jacobi forms with higher rank indices by taking products with
    modular forms.

    INPUT:

    - `k` -- An integer.

    - `m` -- A qudratic form.

    - ``k_mod`` -- An integer.

    - ``prec`` -- An integer.

    - ``jforms`` -- A list of dictionaries, representing Jacobi forms.

    TESTS::

        sage: from sage.modular.jacobi.test_higherrank import _test_higherrank_jacobi_forms__multiplication
        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_forms
        sage: jforms = higherrank_jacobi_forms(4,QuadraticForm(ZZ,1,[1]),3)
        sage: _test_higherrank_jacobi_forms__multiplication(4,QuadraticForm(ZZ,1,[1]),4,3,jforms)
    """
    l = m.dim()
    (r_classes, _) = higherrank_jacobi_r_classes(m)
    m_adj = QuadraticForm(2 * m.matrix().adjoint())
    m_span = m.matrix().row_module()

    P = PolynomialRing(LaurentPolynomialRing(QQ, ['zeta{}'.format(ix) for ix in range(l)]), 'q')
    q = P.gen(0)
    zeta = lambda ix: P.base_ring().gen(ix)

    indices = list(higherrank_jacobi_fe_indices(m, prec, r_classes, reduced=True))
    indices_nonreduced = list(higherrank_jacobi_fe_indices(m, prec, r_classes, reduced=False))
    red = lambda nr: higherrank_jacobi_reduce_fe_index(nr, m, r_classes, m_adj, m_span)

    psis = [dict([ (nr, (_psi[nr] if nr in _psi else 0))
                   for nr in indices])
            for _psi in higherrank_jacobi_forms(k+k_mod, m, prec)]
    psi_span = matrix([[psi[nr] for nr in indices]
                       for psi in psis]).row_module()

    for phi1 in jforms:
        phi1_poly = P.zero()
        for nr in indices_nonreduced:
            (n,r) = nr
            (nrred, s) = red(nr)
            if nrred in phi1:
                phi1_poly += (s**k * phi1[nrred]
                              * q**n
                              * prod(zeta(ix)**r[ix] for ix in range(l)))

        for f in ModularForms(1, k_mod).basis():
            f_poly = f.qexp(prec).polynomial()

            phi2_poly = f_poly * phi1_poly
            if l == 1:
                phi2_vec = vector([phi2_poly[n][r[0]] for (n,r) in indices])
            else:
                phi2_vec = vector([phi2_poly[n][r] for (n,r) in indices])
            assert not phi2_vec.is_zero()
            assert phi2_vec in psi_span
