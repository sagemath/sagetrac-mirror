r"""
Fourier expansions of (weakly holomorphic) vector valued modular form
for Weil representations.

AUTHOR:

- Martin Raum

REFERENCE:

- [Ra] Martin Raum, Computation of Jacobi forms degree 1 and higher rank index.
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

from sage.functions.all import ceil
from sage.matrix.all import matrix, zero_matrix, identity_matrix
from sage.modular.all import CuspForms
from sage.modular.jacobi.higherrank import (
    higherrank_jacobi_forms,
    higherrank_jacobi_r_classes,
    _higherrank_jacobi_reduce_fe_index__r
)
from sage.modules.all import vector, zero_vector
from sage.quadratic_forms.all import QuadraticForm
from sage.rings.all import ZZ, QQ, PolynomialRing, gcd

from copy import copy


def vector_valued_modular_forms(k, L, prec):
    r"""
    A basis of the space of vector valued modular forms of weight `k`
    and the conjugate of the Weil representation attached to the
    quadratic form `L`.  Compute Fourier expansions up to precision ``prec``.

    INPUT:

    - `k` -- An integer.

    - `L` -- A quadratic form over `\ZZ`.

    - ``prec`` -- A nonnegative integer.

    OUTPUT:

    A list of dictionaries that encode the Fourier coefficients of
    vector valued modular forms. See :meth:`theta_decomposition` for a
    description of these expansions.

    EXAMPLES::

        sage: from sage.modular.jacobi.all import *
        sage: L = QuadraticForm(matrix(2, [2,0,0,2]))
        sage: vector_valued_modular_forms(3, L, 5)
        [{(0, 0): {0: 1, 1: 60, 2: 252, 3: 544, 4: 1020, 5: 1560},
          (0, 1): {3/4: 32, 7/4: 192, 11/4: 480, 15/4: 832, 19/4: 1440},
          (1, 0): {3/4: 32, 7/4: 192, 11/4: 480, 15/4: 832, 19/4: 1440},
          (1, 1): {1/2: 12, 3/2: 160, 5/2: 312, 7/2: 960, 9/2: 876}}]

    ::

        sage: L = QuadraticForm(matrix([-2]))
        sage: vector_valued_modular_forms(9/2, L, 5)
        Traceback (most recent call last):
        ...
        ValueError: Quadratic form must be positive definite. Use L.stabily_equivalent_positive_definite_quadratic_form() inorder to obtain one that is stabily equivalent.

    TESTS:

    See ``test_vector_valued.py``.
    """
    if not L.is_positive_definite():
        raise ValueError("Quadratic form must be positive definite. Use L.stabily_equivalent_positive_definite_quadratic_form() inorder to obtain one that is stabily equivalent.")

    if (2 * k + L.dim()) % 2 != 0:
        return []
    else:
        k_jac = ZZ(k + L.dim() / ZZ(2))

    L_adj = QuadraticForm(2 * L.matrix().adjoint())

    (r_classes,_) = higherrank_jacobi_r_classes(L)
    max_n_shift = max(L_adj(r_class[0]) / ZZ(2*L.det()) for r_class in r_classes)

    return [theta_decomposition(phi, L, r_classes)
            for phi in higherrank_jacobi_forms(k_jac, L, prec + ceil(max_n_shift))]


def vector_valued_modular_forms_weakly_holomorphic(k, L, order, prec):
    r"""
    A basis of weakly holomorphic vector valued modular forms of given
    order, weight `k`, and for the Weil representation associated to
    `L`.

    INPUT:

    - `k` -- An integer.

    - `L` -- A quadratic form over `\ZZ`.

    - ``order`` -- A positive integer, the maximal pole order at infinity.

    - ``prec`` -- A nonnegative integer.

    EXAMPLES::

        sage: from sage.modular.jacobi.all import *
        sage: k = -1
        sage: L = QuadraticForm(matrix(2, [2,1,1,2]))
        sage: order = 1
        sage: B = 5
        sage: vector_valued_modular_forms_weakly_holomorphic(k, L, order, B)
        [{(0, 0): {-1: 1, 0: 24, 1: 98730, 2: 8033552, 3: 268396434, 4: 5498892864},
          (0, 1): {-1/3: 11/3,
           2/3: -49016/3,
           5/3: -6384917/3,
           8/3: -89122696,
           11/3: -2095952744},
          (0, 2): {-1/3: -11/3,
           2/3: 49016/3,
           5/3: 6384917/3,
           8/3: 89122696,
           11/3: 2095952744}},
         {(0, 0): {-1: 0, 0: 1, 1: 21, 2: 198, 3: 1236, 4: 6168},
          (0, 1): {-1/3: -1/18, 2/3: 41/9, 5/3: 446/9, 8/3: 1039/3, 11/3: 5522/3},
          (0, 2): {-1/3: 1/18, 2/3: -41/9, 5/3: -446/9, 8/3: -1039/3, 11/3: -5522/3}}]

    TESTS:

    See ``test_vector_valued.py``.
    """
    if order not in ZZ:
        raise NotImplementedError("Only integral orders allowed.")
    if order <= 0:
        raise ValueError("Order at infinity must be negative.")

    vvforms = vector_valued_modular_forms(k + 12*order, L, prec + order)

    delta = CuspForms(1,12).gen(0).qexp(prec + 1 + order)
    delta_order = (delta.shift(-1)**(-order)).truncate(prec + order)

    R = PolynomialRing(QQ, 'q')
    wvvforms = []
    for f in vvforms:
        wf = {}
        for (mu,fe) in f.items():
            wf[mu] = {}

            fe_shift = min(fe.keys())
            fe_shift_reduced = (fe_shift.numerator() % fe_shift.denominator()) / fe_shift.denominator()

            fe_poly = R([fe[n + fe_shift_reduced]
                         for n in range(prec + order if fe_shift_reduced == 0 else prec + order - 1)])

            fe_poly *= delta_order
            for n in range(-order, prec if fe_shift_reduced == 0 else prec-1):
                wf[mu][n + fe_shift_reduced] = fe_poly[n + order]
        wvvforms.append(wf)

    return wvvforms

def vector_valued_modular_forms_weakly_holomorphic_with_principal_part(k, L, principal_part, prec):
    r"""
    Return a weak holomorphic vector valued modular form with given principal part.

    Raises ``ValueError`` if no such form exists.

    INPUT:

    - `k` -- A half-integer.

    - `L` -- A quadratic form over `\ZZ`.

    - ``principal_part`` -- A dictionary whose keys represent elements
      of the discriminant group and whose values
      are dictionaries corresponding to Fourier
      expansions of a component.  E.g.
      {(0,): {-2: 2, -1: 2}, (1,): {-1/4: 3}}

    - ``prec`` -- A positive integer.

    OUTPUT:

    A dictionary of dictionaries that represents the Fourier expansion
    of a weakly holomorphic modular form.

    EXAMPLES::

        sage: from sage.modular.jacobi.all import *
        sage: k = -1
        sage: L = QuadraticForm(matrix(2, [2, 1, 1, 2]))
        sage: pp = {(0,0): {-1: 1}}
        sage: vector_valued_modular_forms_weakly_holomorphic_with_principal_part(k, L, pp, 5)
        {(0, 0): {-1: 1, 0: 90, 1: 100116, 2: 8046620, 3: 268478010, 4: 5499299952},
         (0, 1): {-1/3: 0,
          2/3: -16038,
          5/3: -2125035,
          8/3: -89099838,
          11/3: -2095831260},
         (0, 2): {-1/3: 0, 2/3: 16038, 5/3: 2125035, 8/3: 89099838, 11/3: 2095831260}}

    TESTS:

    See ``test_vector_valued.py``.
    """
    order = -min(min(fe.keys()) for fe in principal_part.values())
    ## TODO: as soon as better weakly holomorphic forms function is
    ## implemented, we can remove this
    order = ceil(order)
    vvforms = vector_valued_modular_forms_weakly_holomorphic(k, L, order, prec)

    L_span = L.matrix().row_module()
    L_adj = QuadraticForm(2 * L.matrix().adjoint())
    mu_module = L_span.ambient_module() / L_span
    mu_indices = dict([(mu, ix) for (ix, mu) in enumerate(mu_module)])


    pp_matrix = zero_matrix(QQ, len(mu_indices) * order, len(vvforms))
    for (col, vvf) in enumerate(vvforms):
        for (mu, fe) in vvf.items():
            mu_ix = mu_indices[mu_module(vector(mu))]

            n_shift = L_adj(mu) / (2 * L.det())
            n_shift = (n_shift.numerator() % n_shift.denominator()) / n_shift.denominator()
            if n_shift == 0:
                n_shift = 1

            for (n, coeff) in fe.items():
                if n < 0:
                    pp_matrix[mu_ix * order - (n + n_shift), col] = coeff

    pp_vector = vector(QQ, len(mu_indices)*order)
    for (mu, fe) in principal_part.items():
        mu_ix = mu_indices[mu_module(vector(mu))]

        n_shift = L_adj(mu) / (2 * L.det())
        n_shift = (n_shift.numerator() % n_shift.denominator()) / n_shift.denominator()
        if n_shift == 0:
            n_shift = 1

        for (n, coeff) in fe.items():
            if n < 0:
                assert n + n_shift in ZZ
                pp_vector[mu_ix * order - (n + n_shift)] = coeff

    try:
        coords = pp_matrix.solve_right(pp_vector)
    except:
        raise ValueError("Given principal part ({}) can not be constructed for weight {} and index {}".format(principal_part, k, m))

    return _sum_mul_vvforms(coords, vvforms, L_span)


def theta_decomposition(phi, m, r_classes):
    r"""
    Apply the theta decomposition to a Fourier expansion `\phi`.

    INPUT:

    - ``phi`` -- A dictionary representing the Fourier expansion of a
      Jacobi form.

    - `m` -- A quadratic form over `m`.

    - ``r_classes`` -- A list of lists of tuples.

    OUTPUT:

    A dictionary whose keys are lifts of elements of a discrimiant
    group, and whose values are dictionaries whose keys are rationals
    (the exponents of `q`) and whose values are also rationals (the
    corresponding coefficients).


    EXAMPLES::

        sage: from sage.modular.jacobi.all import *
        sage: k = 9
        sage: L = QuadraticForm(matrix(2, [2,1,1,2]))
        sage: jforms = higherrank_jacobi_forms(k, L, 5)
        sage: (r_classes, _) = higherrank_jacobi_r_classes(L)
        sage: theta_decomposition(jforms[0], L, r_classes)
        {(0, 0): {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
         (0, 1): {2/3: -1, 5/3: 16, 8/3: -104, 11/3: 320},
         (0, 2): {2/3: 1, 5/3: -16, 8/3: 104, 11/3: -320}}
    """
    L = m
    L_span = L.matrix().row_module()
    L_adj = QuadraticForm(2 * L.matrix().adjoint())
    mu_list = [tuple(mu.lift()) for mu in L_span.ambient_module() / L_span]

    r_to_mu = {}
    for mu in mu_list:
        try:
            (r, sign) = _higherrank_jacobi_reduce_fe_index__r(mu, r_classes, L_span)
        except:
            raise AssertionError(mu, r_classes)
        try:
            r_to_mu[r].append((mu, sign))
        except KeyError:
            r_to_mu[r] = [(mu, sign)]


    f = dict((mu,{}) for mu in mu_list)
    for ((n,r),coeff) in phi.items():
        disc = n - L_adj(r) / ZZ(2*L.det())
        for (mu, sign) in r_to_mu[r]:
            f[mu][disc] = sign * coeff

    return f


def _mul_scalar_vvform(c, f):
    r"""
    Multiplication of a Fourier expansion by a constant.

    INPUT:


    - `c` -- A constant.

    - `f` -- A dictionary representing the Fourier expansion of a
             weakly holomorphic modular forms.

    OUTPUT:

    A dictionary representing a Fourier expansion.

    EXAMPLES::

        sage: from sage.modular.jacobi.vector_valued import _mul_scalar_vvform
        sage: _mul_scalar_vvform(5, {(0,1): {-1: 4, 0: 2}})
        {(0, 1): {-1: 20, 0: 10}}
    """
    res = {}
    for (mu,fe) in f.items():
        res[mu] = dict((n,c*coeff) for (n,coeff) in fe.items())
    return res


def _add_vvforms(f, g, L_span):
    r"""
    Addition of two Fourier expansions.

    INPUT:

    - `f` -- A dictionary representing the Fourier expansion of a
             weakly holomorphic modular forms.

    - `g` -- A dictionary representing the Fourier expansion of a
             weakly holomorphic modular forms.

    - ``L_span`` -- A module over `\ZZ`.

    OUTPUT:

    A dictionary representing a Fourier expansion.

    EXAMPLES::

        sage: from sage.modular.jacobi.vector_valued import _add_vvforms
        sage: f1 = {(0,1): {-1: 4, 0: 2}}
        sage: f2 = {(0,1): {-1: 2, 2: 3}, (1,1): {2: 2}}
        sage: L_span = (ZZ**2).span([(1,2), (0,3)])
        sage: _add_vvforms(f1, f2, L_span)
        {(0, 1): {-1: 6, 0: 2, 2: 3}, (1, 1): {2: 2}}
    """
    res = copy(g)
    for (mu, fe) in f.items():
        for gmu in g.keys():
            if vector(mu) - vector(gmu) in L_span:
                break
        else:
            res[mu] = fe
            break

        for (n, coeff) in fe.items():
            if n in res[gmu]:
                res[gmu][n] = res[gmu][n] + coeff
            else:
                res[gmu][n] = coeff
    return res


def _sum_mul_vvforms(coefficients, vvforms, L_span):
    r"""
    Linear combination of Fourier expansions.

    INPUT:

    - ``coefficients`` -- A list of constants.

    - `vvforms` -- A list of dictionaries representing the Fourier expansion of a
                   weakly holomorphic modular forms.

    - ``L_span`` -- A module over `\ZZ`.

    OUTPUT:

    A dictionary representing a Fourier expansion.

    EXAMPLES::

        sage: from sage.modular.jacobi.vector_valued import _sum_mul_vvforms
        sage: f1 = {(0,1): {-1: 4, 0: 2}}
        sage: f2 = {(0,1): {-1: 4, 0: -2}}
        sage: f3 = {(0,1): {-1: 2, 2: 3}, (1,1): {2: 2}}
        sage: L_span = (ZZ**2).span([(1,2), (0,3)])
        sage: _sum_mul_vvforms([1,-1,3], [f1,f2,f3], L_span)
        {(0, 1): {-1: 6, 0: 4, 2: 9}, (1, 1): {2: 6}}
    """
    res = {}
    for (c,f) in zip(coefficients, vvforms):
        res = _add_vvforms(res, _mul_scalar_vvform(c,f), L_span)
    return res

################################################################################
## Quadratic forms
################################################################################


def stably_equivalent_positive_definite_quadratic_form(L, split_off_E8=False):
    r"""
    Find the Gram matrix of a positive definite quadratic form
    which is stabily equivalent to `L`.

    INPUT:

    - `L` -- A quadratic form.

    - ``split_off_E8``  -- If ``True``, then split off as many copies
      of E8 as possible.

    EXAMPLES::

        sage: from sage.modular.jacobi.all import *
        sage: m = QuadraticForm(diagonal_matrix([-2]))
        sage: stably_equivalent_positive_definite_quadratic_form(m)
        Quadratic form in 7 variables over Integer Ring with coefficients:
        [ 1 -1 -1 -1 -1 -1 1 ]
        [ * 1 1 1 1 1 0 ]
        [ * * 1 0 1 1 0 ]
        [ * * * 1 1 0 -1 ]
        [ * * * * 1 1 0 ]
        [ * * * * * 1 0 ]
        [ * * * * * * 1 ]

    TESTS:

    See ``test_vector_valued.py``.
    """
    ## This is a workaround since GP crashes
    is_positive_definite = lambda mat: all( mat[:n,:n].det() > 0
                                            for n in range(1, mat.nrows() + 1) )

    while not is_positive_definite(L.matrix()):
        L = _split_off_hyperbolic(L)

    if split_off_E8:
        while True:
            try:
                L = _split_off_E8(L)
            except ValueError, err:
                if err.message == "Lattice does not contain E8":
                    break
                else:
                    raise

    return L.lll()


def _split_off_hyperbolic(L):
    r"""
    Decompose `L \oplus E_8 = U \oplus M` and return `M`.

    INPUT:

    - `L` -- A quadratic form.

    OUTPUT:

    A quadratic form.

    EXAMPLES::

        sage: from sage.modular.jacobi.vector_valued import _split_off_hyperbolic
        sage: m = QuadraticForm(diagonal_matrix([-4, 2]))
        sage: _split_off_hyperbolic(m)
        Quadratic form in 8 variables over Integer Ring with coefficients:
        [ 30 0 -4 4 -8 0 0 0 ]
        [ * 1 0 0 0 0 0 0 ]
        [ * * 1 -1 0 0 0 0 ]
        [ * * * 1 -2 0 0 0 ]
        [ * * * * 2 -1 0 0 ]
        [ * * * * * 1 -1 0 ]
        [ * * * * * * 1 -1 ]
        [ * * * * * * * 1 ]

    TESTS:

    See ``test_vector_valued.py``.
    """
    Lmat = L.matrix()

    ## TODO: This should be implemented in sage.quadratic_forms
    E8mat = matrix(ZZ, 8,
                   [2, -1, 0, 0, 0, 0, 0, 0,
                    -1, 2, -1, 0, 0, 0, 0, 0,
                    0, -1, 2, -1, 0, 0, 0, -1,
                    0, 0, -1, 2, -1, 0, 0, 0,
                    0, 0, 0, -1, 2, -1, 0, 0,
                    0, 0, 0, 0, -1, 2, -1, 0,
                    0, 0, 0, 0, 0, -1, 2, 0,
                    0, 0, -1, 0, 0, 0, 0, 2])
    E8 = QuadraticForm(E8mat)

    ## This is a workaround since GP crashes
    is_positive_definite = lambda L: all( L[:n,:n].det() > 0
                                          for n in range(1, L.nrows() + 1) )

    ## we add corrections to L, so that we obtain a positive definite lattice
    cur_cor = 2
    Lcormat = Lmat + cur_cor * identity_matrix(L.dim())
    while not is_positive_definite(Lcormat):
        cur_cor += 2
        Lcormat = Lmat + cur_cor * identity_matrix(L.dim())
    Lcor = QuadraticForm(Lcormat)

    Lcor_length_inc = 2 * gcd(Lcormat.list() + [e // 2 for e in Lcormat.diagonal()])
    cur_Lcor_length = 0

    ## find a vector that is negative, whose absolute norm is as small
    ## as possible
    n = 0
    for i in range(Lcor.dim()):
        if Lcormat[i, i] - cur_cor < 0:
            if n == 0 or n > cur_cor - Lcormat[i, i]:
                n = (cur_cor - Lcormat[i, i]) // 2
                a_ind = i
    if n != 0:
        a = zero_vector(Lcor.dim())
        a[a_ind] = 1
    else:
        n = 0
        while n == 0:
            cur_Lcor_length += Lcor_length_inc
            short_vectors = Lcor.short_vector_list_up_to_length(cur_Lcor_length+1, True)

            for length in range(cur_Lcor_length - Lcor_length_inc, cur_Lcor_length+1):
                for a in short_vectors[length]:
                    if L(a) < 0:
                        n = -L(a)
                        break
                if n != 0:
                    break

    ## by enumeration of short vectors, find a pair of vectors v, w in E8 of
    ## length n such that their scalar product equals n-1
    short_vectors = E8.short_vector_list_up_to_length(n+1, False)[n]
    for v in short_vectors:
        for w in short_vectors:
            if v * E8mat * w == 2*n-1:
                LE8_mat = Lmat.block_sum(E8mat)
                v_form = vector( list(a) + list(v) ) * LE8_mat
                w_form = vector( list(a) + list(w) ) * LE8_mat
                Lred_basis = matrix(ZZ, [v_form, w_form]).right_kernel().basis_matrix().transpose()
                Lred_basis = matrix(ZZ, Lred_basis)

                return QuadraticForm(Lred_basis.transpose() * LE8_mat * Lred_basis)

    ## There must be a pair of suitable vectors.
    raise RuntimeError( "Unexpectedly reached end of iteration through short vectors in E8." )


def _split_off_E8(L):
    r"""
    Decompose `L = E_8 \oplus M` and return `M`.

    INPUT:

    - `L` -- A positive definite quadratic form.

    OUTPUT:

    A quadratic form.

    EXAMPLES::

        sage: from sage.modular.jacobi.vector_valued import _split_off_E8
        sage: E8mat = matrix(ZZ, 8, [2,-1,0,0,0,0,0,0,  -1,2,-1,0,0,0,0,0,  0,-1,2,-1, 0,0,0,-1,  0,0,-1,2,-1,0,0,0,  0,0,0,-1,2,-1,0,0,  0,0,0,0,-1,2,-1,0,  0,0,0,0,0,-1,2,0,  0,0,-1,0,0,0,0,2])
        sage: L = QuadraticForm(E8mat.block_sum(matrix([[2]])))
        sage: _split_off_E8(L)
        Quadratic form in 1 variables over Integer Ring with coefficients:
        [ 1 ]

    TESTS:

    See ``test_vector_valued.py``.
    """
    assert L.is_positive_definite()

    root_vectors = L.short_vector_list_up_to_length(2, False)[1]
    if len(root_vectors) < 240:
        raise ValueError( "Lattice does not contain E8" )


    Lmat = L.matrix()

    ## Collect pairs of roots with scalar product 0 or -1.
    root_pairs_0 = dict()
    root_pairs_1 = dict()
    for (a_ix, a) in enumerate(root_vectors):
        a_dual = a * Lmat
        for (b_ix, b) in enumerate(root_vectors):
            product = a_dual * b

            if product == 0:
                try:
                    root_pairs_0[a_ix].add(b_ix)
                except KeyError:
                    root_pairs_0[a_ix] = set([b_ix])
            elif product == -1:
                try:
                    root_pairs_1[a_ix].add(b_ix)
                except KeyError:
                    root_pairs_1[a_ix] = set([b_ix])

    ## Build up a chain of roots that form a basis of E8
    root_chain = []
    possible_extensions = dict()
    possible_extensions[-1] = set(range(len(root_vectors)))
    while True:
        i = len(root_chain) - 1

        try:
            root_chain.append(possible_extensions[i].pop())
            extensions = root_pairs_1[root_chain[-1]]
        except KeyError:
            if len(root_chain) == 0:
                raise ValueError( "Lattice does not contain E8" )
            root_chain = root_chain[:-1]
            continue

        for j in range(i + 1):
            extensions = extensions.intersection(root_pairs_0[root_chain[j]])
        possible_extensions[i + 1] = extensions

        if len(root_chain) == 7:
            possible_finals = possible_extensions[2]
            for j in range(3, 7):
                possible_finals = possible_finals.intersection(root_pairs_0[root_chain[j]])

            ## test whether we can extend the chain
            if len(possible_finals) != 0:
                final_vector = possible_finals.pop()

                E8_copy = [root_vectors[root_chain[j]] for j in range(7)]
                E8_copy.append(root_vectors[final_vector])

                E8_copy_dual = matrix(ZZ, E8_copy) * Lmat
                new_basis = E8_copy_dual.right_kernel().basis_matrix()
                new_basis = matrix(ZZ, new_basis)

                return QuadraticForm(new_basis * Lmat * new_basis.transpose())
