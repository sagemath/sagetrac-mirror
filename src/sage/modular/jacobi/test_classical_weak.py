r"""
Tests for classical_weak.py.

AUTHOR:

- Martin Raum
"""

#===============================================================================
# 
# Copyright (C) 2010-2014 Martin Raum
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

from sage.rings.all import (PolynomialRing, LaurentPolynomialRing,
                            PowerSeriesRing,
                            ZZ, QQ, gcd)

from sage.functions.all import factorial, gegenbauer
from sage.matrix.all import matrix
from sage.modular.all import ModularForms, Gamma1
from sage.modules.all import vector, span



from sage.modular.jacobi.classical_weak import (
    classical_weak_jacobi_fe_indices,
    reduce_classical_jacobi_fe_index,
    classical_weak_jacobi_forms,
    classical_weak_jacobi_taylor_coefficients
)

def test_classical_weak_jacobi_forms():
    prec = 20

    for k in [10, 11, 15, 18]:
        for m  in range(1, 4):

            yield (_test_classical_weak_jacobi_forms__taylor_coefficients,
                   k, m, prec)

            if k%2 == 0:
                nu_bound = 20
                yield (_test_classical_weak_jacobi_forms__taylor_coefficient_modularity,
                       10, k, m, prec)

            kmod = 6
            yield (_test_classical_weak_jacobi_forms__multiplication,
                   k, m, k_mod, prec)

            for torsion_point in [0, 1/2]:
                yield (_test_classical_weak_jacobi_forms__torsion_point,
                       torsion_point, k, m, prec)

def _test_classical_weak_jacobi_forms__taylor_coefficients(k, m, prec):
    r"""
    INPUT:
    
    - `k -- An integer.
    
    - `m` -- A non-negative integer.

    - ``prec`` -- A non-negative integer that corresponds to a precision of
                  the q-expansion.
    """
    jacobi_forms = classical_weak_jacobi_forms(k, m, prec)
    taylor_expansions = classical_weak_jacobi_taylor_coefficients(k, m)
    assert all( proj == f
                for (phi, tcs) in zip(jacobi_forms, taylor_expansions)
                for (proj, f) in zip(
                        _taylor_coefficients(phi, k, m, prec),
                        _predicted_taylor_coefficients(tcs, prec) ) )

def _taylor_coefficients(expansion, k, m, prec) :
    r"""
    Normalized Taylor coefficients of a Jacobi form.

    INPUT:

    - ``expansion`` -- A dictionary that corresponds to the Fourier
                       expansion of a classical weak Jacobi form.

    - `k -- An integer.

    - `m -- An integer.

    - ``prec`` -- An integer.

    OUTPUT:

    A list of power series in `q`.
    """
    R = PowerSeriesRing(ZZ, 'q'); q = R.gen(0)

    projs = list()
    for pw in (range(0, 2*m + 1, 2) if k % 2 == 0 else range(1, 2*m - 1, 2)):
        proj = dict( (n, 0) for n in range(prec) )
        for (n, r) in classical_weak_jacobi_fe_indices(m, prec):
            ((nred, rred), sign) = reduce_classical_jacobi_fe_index((n,r), m)
            try :
                proj[n] +=  (sign * r)**pw * expansion[(nred, rred)]
            except (KeyError, ValueError) :
                pass

        projs.append(proj)

    gcd_projs = [gcd(proj.values()) for proj in projs]
    gcd_projs = [g if g != 0 else 1 for g in gcd_projs]
    projs = [sorted(proj.iteritems()) for proj in projs]
    projs = [ R([c for (_, c) in proj]).add_bigoh(prec) / gcd_proj
              for (proj, gcd_proj) in zip(projs, gcd_projs) ]

    return projs

def _predicted_taylor_coefficients(fs, prec) :
    r"""
    Given a list of power series, which are the corrected Taylor coefficients
    of a Jacobi form, return the normalized, uncorrected ones, assuming that
    all but one `f` vanish.

    INPUT:

    - ``fs`` -- A list of power series.

    - ``prec`` -- An integer.

    OUPUT:

    - A list of power series.
    """
    R = PowerSeriesRing(ZZ, 'q'); q = R.gen(0)

    diff = lambda f: f.derivative().shift(1)
    normalize = lambda f: f / gcd(f.list()) if f != 0 else f
    diffnorm = lambda f,l: normalize(reduce(lambda a, g: g(a), l*[diff], f))

    taylor_coefficients = list()
    allf = R(0)
    for f in fs :
        allf = f(prec) + diffnorm(allf, 1)            
        taylor_coefficients.append(allf)

    return taylor_coefficients

def _test_classical_weak_jacobi_forms__taylor_coefficient_modularity(nu_bound, k, m, prec):
    assert k % 2 == 0

    phis = classical_weak_jacobi_forms(k, m, prec)

    fss = [[_jacobi_corrected_taylor_expansions(nu, phi, k, m, prec) for phi in phis]
           for nu in range(0,nu_bound,2)]
    fss_vec = [ [vector(f.padded_list(prec)) for f in fs] for fs in fss ]

    mf_spans = [ span([vector(b.qexp(prec).padded_list(prec)) for b in ModularForms(1, k + 2 * nu).basis()])
                 for nu in range(0,nu_bound,2) ] 

    assert all(f_vec in mf_span
               for (fs_vec, mf_span) in zip(fss_vec, mf_spans)
               for f_vec in fs_vec)


def _corrected_taylor_expansions(nu, phi, k, m, prec):
    r"""
    Return the ``2 nu``-th corrected Taylor coefficient.

    INPUT:

    - ``nu`` -- An integer.  

    - ``phi`` -- A Fourier expansion of a Jacobi form.

    - `k` -- An integer.

    - `m` -- An integer.

    - ``prec`` -- An integer.

    OUTPUT:

    A power series in `q`.

    ..TODO:

    Implement this for odd Taylor coefficients.
    """
    assert nu % 2 == 0

    ## We use EZ85, p.29 (3), the factorial in one of the factors is missing
    factors = [ (-1)**mu * factorial(2*nu) * factorial(k + 2*nu - mu - 2) / ZZ(factorial(mu) * factorial(2*nu - 2*mu) * factorial(k + nu - 2))
                for mu in range(nu + 1) ]
    gegenbauer = lambda n, r: sum( f * r**(2 * nu - 2 * mu) * n**mu 
                                   for (mu,f) in enumerate(factors) )

    coeffs = dict( (n, QQ(0)) for n in range(prec) )
    for (n, r) in classical_weak_jacobi_fe_indices(m, prec):
        (nrred, s) = reduce_classical_jacobi_fe_index((n,r), m)
        coeffs[n] += s**k * gegenbauer(m*n, r) * phi[nrred]

    return PowerSeriesRing(QQ, 'q')(coeffs)

def _test_classical_weak_jacobi_forms__multiplication(k, m, k_mod, prec):
    r"""
    Check that the product of weak Jacobi forms of given weight and
    index with modular forms of given weight is a weak Jacobi form.

    INPUT:

    - `k` -- An integer.

    - `m -- A positive integer.

    - `k_mod` -- An integer.

    - ``prec`` -- An integer.
    """
    P = PolynomialRing(LaurentPolynomialRing(QQ, 'zeta'), 'q')
    q = P.gen(0)
    zeta = P.base_ring().gen(0)

    indices = list(classical_weak_jacobi_fe_indices(m, prec, reduced=True))
    red = lambda nr: reduce_classical_jacobi_fe_index(nr,m)

    psis = [dict([ (nr, (_psi[nr] if nr in _psi else 0)) for nr in indices]) for _psi in classical_weak_jacobi_forms(k+k_mod, m, prec)]
    psi_span = matrix([[psi[nr] for nr in indices] for psi in psis]).row_module()

    for phi1 in classical_weak_jacobi_forms(k, m, prec):
        phi1_poly = P.zero()
        for nr in classical_weak_jacobi_fe_indices(m, prec):
            (nrred, s) = red(nr)
            if nrred in phi1:
                phi1_poly += s**k * phi1[nrred] * q**n * zeta**r

        for f in ModularForms(1, k_mod).basis():
            f_poly = f.qexp(prec).polynomial()

            phi2_poly = f_poly * phi1_poly
            phi2_vec = vector([phi2_poly[nr[0]][nr[1]] for nr in indices])
            assert phi2_vec in psi_span

def _test_classical_weak_jacobi_forms__torsion_point(torsion_point, k, m, prec):
    jforms = classical_weak_jacobi_forms(k, m, prec)
    fs = [_eval_at_torsion_point(torsion_point, phi, k, m, prec) for phi in jforms]
    fs_vec = [vector(f.padded_list(prec)) for f in fs]

    mf_span = span([vector(b.qexp(prec).padded_list(prec))
                    for b in ModularForms(Gamma1(torsion_point.denominator()**3), k).basis()])

    assert all(f_vec in mf_span for f_vec in fs_vec)

def _eval_at_torsion_point(torsion_point, phi, k, m, prec) :
    r"""
    Given a dictonary that represents the Fourier expansion of a
    Jacobi form, return the specialization to ``torsion_point``.

    INPUT:

    - ``torsion_point`` -- A rational.

    - ``phi`` -- A Fourier expansion of a Jacobi form.

    - `k` -- An integer.

    - `m` -- An integer.

    - ``prec`` -- An integer.

    OUPUT:

    - A power series.
    """
    from sage.rings.all import CyclotomicField

    K = CyclotomicField(QQ(torsion_point).denominator()); zeta = K.gen()
    R = PowerSeriesRing(K, 'q'); q = R.gen(0)

    coeffs = dict( (n, QQ(0)) for n in range(prec) )
    for (n, r) in classical_weak_jacobi_fe_indices(m, prec):
        (nrred, s) = reduce_classical_jacobi_fe_index((n,r), m)
        coeffs[n] += s**k * zeta**r * phi[nrred]

    return PowerSeriesRing(K, 'q')(coeffs) 
