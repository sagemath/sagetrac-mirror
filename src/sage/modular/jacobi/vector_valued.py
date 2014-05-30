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


def vector_valued_modular_forms(k, L, prec):
    r"""
    A basis of the space of vector valued modular forms of weight `k`
    and the conjugate of the Weil representation attached to the
    quadratic form `L`.  Compute Fourier expansions up to precision ``prec``.

    INPUT:

    - `k` -- An integer.

    - `L` -- A quadratic form over `\Z`.

    - ``prec`` -- A nonnegative integer.

    OUTPUT:

    A list of dictionaries that encode the Fourier coefficients of
    vector valued modular forms.  See `meth:theta_decomposition` for a
    description of these expansions.

    EXAMPLES::

        sage: L = QuadraticForm(matrix([2]))
        sage: vector_valued_modular_forms(5/2, L, 5)
        ???
        sage: L = QuadraticForm(matrix([-2]))
        sage: vector_valued_modular_forms(5/2, L, 5)
        ???

    TESTS:

    See ``test_vector_valued.py:test_vector_valued_modular_forms``.
    """
    ## TODO: If L is negative definite find equivalent positive
    ## definite quadratic form.
    assert L.is_positive_definite()

    L_adj = QuadraticForm(2 * L.matrix().adjoint())

    (r_classes,_) = _higherrank_jacobi_r_classes(L)
    max_n_shift = max(L_adj(r_class[0]) / ZZ(2*L.det()) for r_class in r_classes)

    return [theta_decomposition(phi, L, r_classes)
            for phi in higher_rank_jacobi_forms(k, L, prec + ceil(max_n_shift))]

def  weakly_holomorphic_vector_valued_modular_forms(k, L, order, prec):
    r"""
    A basis of weakly holomorphic vector valued modular forms of given
    order, weight `k`, and for the Weil representation associated to
    `L`.

    INPUT:

    - `k` -- An integer.

    - `L` -- A quadratic form over `\Z`.

    - ``order`` -- An integer (typically negative).

    - ``prec`` -- A nonnegative integer.

    ..TODO:

    Insert examples.

    TESTS:

    See ``test_vector_valued.py:test_weakly_holomorphic_vector_valued_modular_forms``.
    """
    if order not in ZZ:
        raise NotImplementedError("Only integral orders allowed.")
        
    vvforms = vector_valued_modular_forms(k - 12*order, L, prec - order)
    if order == 0: return vvforms

    delta = CuspForms(1,12).gen(0).qexp(prec+1)
    delta_order = (delta.shift(-1)**order)

    wvvforms = []
    for f in vvforms:
        wf = {}
        for (mu,fe) in f.items():
            wf[mu] = {}

            fe_order = min(fe.keys())
            fe_degree = max(fe.keys())
            fe_poly = R([fe[fe_order + n] for n in range(fe_degree - fe_order + 1)])

            fe_poly *= delta_order
            for n in range(prec - order):
                wf[mu][fe_order + order + n] = fe_poly[n]
        wvvforms.append(wf)

    return wvvforms

def weak_holomorphic_vector_valued_modular_with_principle_part( m, L, principal_part, prec ) :
    r"""
    Return a weak holomorphic vector valued modular form with given principle part.

    Raises ``ValueError`` if no such form exits.

    INPUT:
    
    - `k -- A half-integral.
    
    - `L` -- A quadratic form over `\Z`.

    - ``principal_part`` -- A dictionary whose keys represent elements
                            of the discriminant group and whose values
                            are dictionaries corresponding to Fourier
                            expansions of a component.  E.g.
                            {(0,): {-2: 2, -1: 2}, (1,): {-1/4: 3}}
    
    - ``prec`` -- A positive integer.

    OUTPUT:

    - A dictionary of dictionaries that represents the Fourier
      expansion of a weakly holomorphic modular form.
    
    EXAMPLES::
    
        sage: from sage.modular.jacobi.vector_valued import *
        sage: k = -1; L = QuadraticForm(matrix(2, [2, 1, 1, 2]))
        sage: pp = {(0,0): {-1: 1}}
        sage: weak_holomorphic_vector_valued_modular_with_principle_part(k, L, pps, 10)
        ???
    """
    order = min(min(fe.keys()) for fe in principal_part.items())
    ## TODO: if better weakly holomorphic forms function is
    ## implemented, we can remove this
    order = floor(order)
    vvforms = weakly_holomorphic_vector_valued_modular_forms(k, L, order, prec)
    vvpps = map(_principal_part, vvforms)

    L_span = L.matrix().row_module()
    L_adj = QuadraticForm(2 * L.matrix().adjoint())
    mu_module = L_span.ambient_module() / L_span
    mu_indices = dict(reversed(enumerate(mu_module.list())))
    
    pp_matrix = zero_matrix(QQ, len(r_classes)*abs(order), len(vvforms))
    pp_vector = vector(QQ, len(r_classes)*abs(order))
    for (col,vvf) in enumerate(vvforms):
        for (mu,fe) in vvf.items():
            mu_ix = mu_indices[mu_module(mu)]
            n_shift = L_adj(mu.lift())/(2*L.det())
            n_shift = (n_shift.numerator() % n_shift.denominator()) / n_shift.denominator()
            for (n,coeff) in fe.items():
                pp_matrix[mu_ix*abs(order) + abs(n+n_shift), col] = coeff
    for (mu,fe) in pp.items():
        mu_ix = mu_indices[mu_module(mu)]
        n_shift = L_adj(mu.lift())/(2*L.det())
        n_shift = (n_shift.numerator() % n_shift.denominator()) / n_shift.denominator()
        for (n,coeff) in fe.items():
            pp_vector[mu_ix*abs(order) + abs(n+n_shift), col] = coeff
    
    try :
        coords = pp_matrix.solve_right(pp_vector)
    except :
        raise ValueError( "Given principle part can not be constructed" )

    res = dict()
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

    - A dictionary whose keys are lifts of elements of a discrimiant group, and whose
      values are dictionaries whose keys are rationals (the exponents of `q`) and whose
      values are also rationals (the corresponding coefficients).


    EXAMPLES::

        sage: from sage.modular.jacobi.all import *
        sage: k = 9
        sage: L = QuadraticForm(matrix(2, [2,1,1,2]))
        sage: prec = 10
        sage: jforms = higher_rank_jacobi_forms(k, L, prec)
        sage: r_classes = higherrank_jacobi_r_classes(L)
        sage: theta_decomposition(jforms[0], m, r_classes)
        ???
    """
    L = m
    L_span = L.matrix().row_module()
    L_adj = QuadraticForm(2 * L.matrix().adjoint())
    mu_module = L_span.ambient_module() / L_span

    r_to_mu = {}
    for mu in mu_module:
        (r, sign) = _reduce_higherrank_jacobi_fe_index__r(mu.lift(), r_classes, L_span)
        try:
            r_to_mu[r].append((mu, sign))
        except KeyError:
            r_to_mu[r] = [(mu, sign)]


    f = dict((mu.lift(),{}) for mu in mu_module)
    for ((n,r),coeff) in phi.items():
        disc = n - L_adj(r) / ZZ(2*L.det())
        for (mu, sign) in r_to_mu[r]:
            f[mu][disc] = sign * coeff

    return f

def _principle_part(f):
    r"""
    The principal part of a weakly holomorphic vector valued modular form.

    INPUT:

    - `f` -- A dictionary representing the Fourier expansion of a
             weakly holomorphic modular forms.

    OUTPUT:

    A dictionary representing the principal part.

    ..TODO:

    Insert examples.
    """
    res = {}
    for (mu,fe) in f.items():
        res[mu] = {}
        for (n,coeff) in fe.items():
            if n < 0: res[mu][n] = coeff
    return res

def _mul_scalar_vvform(c, f) :
    r"""
    Multiplication of a Fourier expansion by a constant.

    INPUT:

    - `c` -- A constant.

    - `f` -- A dictionary representing the Fourier expansion of a
             weakly holomorphic modular forms.

    OUTPUT:

    A dictionary representing a Fourier expansion.

    ..TODO:

    Insert examples.
    """
    res = {}
    for (mu,fe) in f.items():
        res[mu] = dict((n,c*coeff) for (n,coeff) in fe.items())
    return res

def _add_vvforms(f, g, L_span) :
    r"""
    Addition of two Fourier expansions.

    INPUT:

    - `f` -- A dictionary representing the Fourier expansion of a
             weakly holomorphic modular forms.

    - `g` -- A dictionary representing the Fourier expansion of a
             weakly holomorphic modular forms.

    - ``L_span`` -- A module over `\Z`.

    OUTPUT:

    A dictionary representing a Fourier expansion.

    ..TODO:

    Insert examples.
    """
    res = copy(g)
    for (mu,fe) in f.items():
        for gmu in g.keys():
            if vector(mu) - vector(gmu) in L_span:
                break
        else:
            res[mu] = fe
            break
        for (n,coeff) in fe:
            if n in res[gmu]:
                res[gmu] = res[gmu] + coeff
            else:
                res[gmu] = coeff
    return res

def _sum_mul_vvforms(coefficients, vvforms, L_span):
    r"""
    Linear combination of Fourier expansions.

    INPUT:

    - ``coefficients`` -- A list of constants.

    - `vvforms` -- A list of dictionaries representing the Fourier expansion of a
                   weakly holomorphic modular forms.

    - ``L_span`` -- A module over `\Z`.

    OUTPUT:

    A dictionary representing a Fourier expansion.

    ..TODO:

    Insert examples.
    """
    res = {}
    for (c,f) in zip(coefficients, vvforms):
        res = _add_vvforms(res, _mul_scalar_vvform(c,f), L_span)
    return res
