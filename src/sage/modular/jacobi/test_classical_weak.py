r"""
Test for classical_weak.py.

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

def test__classical_weak_jacobi_forms(prec, k, m, algorithm="skoruppa") :
    r"""
    INPUT:
    
    - ``prec`` -- A non-negative integer that corresponds to a precision of
                  the q-expansion.
    
    - `k -- An integer.
    
    - `m` -- A non-negative integer.

    - ``algorithm`` -- Default: ''skoruppa''.  Only ''skoruppa'' is implemented.

    TESTS::

        sage: from sage.modular.jacobi.classical_weak import *
        sage: from sage.modular.jacobi.test_classical_weak import *
        sage: [ test__classical_weak_jacobi_forms__skoruppa(40, k, m) for (k, m) in [(10, 1), (12, 2), (9, 3)] ]
        [None, None, None]
        sage: [ test__classical_weak_jacobi_forms__skoruppa(40, k, m) for (k, m) in [(7, 5), (10, 10) ]          # long time
        [None, None]
    """
    jacobi_forms = classical_weak_jacobi_forms(prec, k, m, algorithm)
    taylor_expansions = classical_weak_jacobi_taylor_coefficients(k, m)
    assert all( proj == f
                for (phi, tcs) in zip(jacobi_forms, taylor_expansions)
                for (proj, f) in zip(
                        test__jacobi_taylor_coefficients(phi, k),
                        test__jacobi_predicted_taylor_coefficients(tcs, prec) ) )


def _jacobi_taylor_coefficients(expansion, k, m, prec) :
    r"""
    Compute the renormalized Taylor coefficients of
    a Jacobi form.

    INPUT:

    - ``expansion`` -- A dictionary that corresponds to the Fourier
                       expansion of a classical weak Jacobi form.

    - `k -- An integer.

    - `m -- An integer.

    - ``prec`` -- An integer.

    OUTPUT:

    - A list of power series in `q`.

    TESTS:

    See ``meth:test__classical_weak_jacobi_forms``.
    """
    from sage.rings.arith import gcd

    R = PowerSeriesRing(ZZ, 'q'); q = R.gen(0)

    projs = list()
    for pw in (range(0, 2*m + 1, 2) if k % 2 == 0 else range(1, 2*m - 1, 2)):
        proj = dict( (n, 0) for n in range(prec) )
        for (n, r) in classical_weak_jacobi_fe_indices(m, prec):
            ((nred, rred), sign) = reduce_jacobi_fe_index((n,r), m)
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

def _jacobi_predicted_taylor_coefficients(fs, prec) :
    r"""
    Given a list of power series, which are the corrected Taylor coefficients
    of a Jacobi form, return the renormalized uncorrected ones, assuming that
    all but one `f` vanish.

    INPUT:

    - ``fs`` -- A list of power series.

    - ``prec`` -- An integer.

    OUPUT:

    - A list of power series.

    TESTS:

    See ``meth:test__classical_weak_jacobi_forms``.
    """
    from sage.rings.arith import gcd

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

def _test_jacobi_corrected_taylor_expansions(nu, phi, k) :
    r"""
    Return the ``2 nu``-th corrected Taylor coefficient.

    INPUT:

    - ``nu`` -- An integer.  

    - ``phi`` -- A Fourier expansion of a Jacobi form.

    - `k -- An integer.

    OUTPUT:

    - A power series in `q`.

    ..TODO:

    Implement this for odd Taylor coefficients.

    TESTS::

        sage: from sage.modular.jacobi.classical_weak import *
        sage: from sage.modular.jacobi.test_classical_weak import *
        sage: nu_bound = 10
        sage: precision = 100
        sage: k = 10; m = 7
        sage: phis = classical_weak_jacobi_forms(k, m, prec)
        sage: fss = [ [_test_jacobi_corrected_taylor_expansions(nu, phi, k) for phi in phis] for nu in range(nu_bound) ]
        sage: fss_vec = [ [vector(f.padded_list(prec)) for f in fs] for fs in fss ]
        sage: mf_spans = [ span([vector(b.qexp(prec).padded_list(prec)) for b in ModularForms(1, weight + 2 * nu).basis()]) for nu in range(nu_bound) ] 
        sage: all(f_vec in mf_span for (fs_vec, mf_span) in zip(fss_vec, mf_spans) for f_vec in fs_vec)
        True
    """
    assert nu % 2 = 0

    ## We use EZ85, p.29 (3), the factorial in one of the factors is missing
    factors = [ (-1)**mu * factorial(2 * nu) * factorial(weight + 2 * nu - mu - 2) / ZZ(factorial(mu) * factorial(2 * nu - 2 * mu) * factorial(weight + nu - 2))
                for mu in range(nu + 1) ]
    gegenbauer = lambda n, r: sum( f * r**(2 * nu - 2 * mu) * n**mu 
                                   for (mu,f) in enumerate(factors) )

    coeffs = dict( (n, QQ(0)) for n in range(prec) )
    for (n, r) in classical_weak_jacobi_fe_indices(m, prec):
        (nrred, s) = reduce_fe_index((n,r), m)
        coeffs[n] += s**k * gegenbauer(m*n, r) * phi[nnred]

    return PowerSeriesRing(QQ, 'q')(coeffs)

def test_jacobi_torsion_point(phi, k, torsion_point) :
    r"""
    Given a list of power series, which are the corrected Taylor coefficients
    of a Jacobi form, return the specialization to ``torsion_point``.

    INPUT:

    - ``phi`` -- A Fourier expansion of a Jacobi form.

    - `k -- An integer.

    - ``torsion_point`` -- A rational.

    OUPUT:

    - A power series.

    TESTS:

    See jacobi_form_from_taylor_expansion.

        sage: from sage.modular.jacobi.classical_weak import *
        sage: from sage.modular.jacobi.test_classical_weak import *
        sage: prec = 50
        sage: k = 10; m = 7
        sage: phis = classical_weak_jacobi_forms(k, m, prec)
        sage: fs = [_test_jacobi_torsion_point(phi, k, 2/3) for phi in phis]
        sage: fs_vec = [vector(f.padded_list(prec)) for f in fs]
        sage: mf_span = span([vector(b.qexp(prec).padded_list(prec)) for b in ModularForms(GammaH(9, [4]), k).basis()])
        sage: all(f_vec in mf_span for f_vec in fs_vec)
        True
    """
    from sage.rings.all import CyclotomicField

    K = CyclotomicField(QQ(torsion_point).denominator()); zeta = K.gen()
    R = PowerSeriesRing(K, 'q'); q = R.gen(0)

    ch = JacobiFormD1WeightCharacter(weight)

    coeffs = dict( (n, QQ(0)) for n in range(prec) )
    for (n, r) in classical_weak_jacobi_fe_indices(m, prec):
        (nrred, s) = reduce_fe_index((n,r), m)
        coeffs[n] += s**k * zeta**r * phi[nrred]

    return PowerSeriesRing(K, 'q')(coeffs) 

