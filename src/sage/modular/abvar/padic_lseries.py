# -*- coding: utf-8 -*-
r"""
`p`-adic L-series of modular Jacobians with ordinary reduction at `p`

REFERENCES:

- [MTT]_

- [SW]_

AUTHORS:

- William Stein and Jennifer Balakrishnan (2010-07-01): first version

"""

######################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
######################################################################

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.padics.factory import Qp
from sage.rings.infinity import infinity
from sage.rings.all import PowerSeriesRing

from sage.rings.integer import Integer
from sage.arith.misc import valuation, kronecker_symbol, gcd, prime_divisors
from sage.structure.sage_object import SageObject
from sage.misc.all import verbose, get_verbose
from sage.modular.modsym.modsym import ModularSymbols


class pAdicLseries(SageObject):
    r"""
    The `p`-adic L-series of a modular Jacobian

    EXAMPLES:

    An ordinary example::


    """
    def __init__(self, J, p, normalize='L_ratio'):
        r"""
        INPUT:

        -  ``J`` -- modular abelian variety
        -  ``p`` -- a prime of good reduction
        -  ``normalize`` -- ``'L_ratio'`` (default), ``'period'`` or ``'none'``;
           this describes the way the modular symbols are normalized
        """
        self._J = J
        self._level = J.level()
        self._p = ZZ(p)
        self._normalize = normalize
        if not self._p.is_prime():
            raise ValueError("p (=%s) must be a prime" % p)

    def __eq__(self, other):
        r"""
        Compare self and other.

        TESTS::
            sage: lp1 = J0(23)[0].padic_lseries(5)
            sage: lp2 = J0(23)[0].padic_lseries(7)
            sage: lp3 = J0(29)[0].padic_lseries(5)
            sage: lp1 == lp1
            True
            sage: lp1 == lp2
            False
            sage: lp1 == lp3
            False

        """
        if type(self) != type(other):
            return False
        return (self._J, self._p) == (other._J, other._p)

    def modular_symbols_subspace(self):
        """
        """
        M = ModularSymbols(self._level, 2, sign=1)
        self._modular_symbols_space = M
        S = M.cuspidal_submodule()
        N = S.new_subspace()
        D = N.decomposition()
        e = M([0, infinity])
        if len(D) == 1:
            A = D[0]
        else:
            for x in D:
                if x.rational_period_mapping()(e) == 0:
                    A = x
        self._modular_symbols_subspace = A
        return A

    def real_quadratic_field(self):
        """
        The field of definition of the dual eigenform.
        """
        try:
            A = self._modular_symbols_subspace
        except AttributeError:
            A = self._modular_symbols_subspace = self.modular_symbols_subspace()
        v = self._dual_eigenvector = A.dual_eigenvector()
        K_f = self._real_quadratic_field = v[0].parent()
        return K_f

    def psi(self):
        """
        The embedding `\Q(\alpha) \into \Q_p(a)` sending `\alpha \mapsto a`.
        """
        try:
            K_f = self._real_quadratic_field
        except AttributeError:
            K_f = self._real_quadratic_field = self.real_quadratic_field()
        p = self._p
        Q = Qp(p)
        F = Q.extension(K_f.defining_polynomial(), names='a')
        a = F.gen()
        psi = self._psi = K_f.hom([a])
        return psi

    def modular_symbol(self, r):
        """
        Compute the modular symbol at the cusp `r`.
        """
        v = self._dual_eigenvector
        psi = self._psi
        M = self._modular_symbols_space
        s = M([r, infinity])
        return psi(v.dot_product(s.element()))

    def abelian_variety(self):
        r"""
        Return the abelian variety to which this `p`-adic L-series is associated.

        EXAMPLES::

            sage: L = J0(23)[0].padic_lseries(5)
            sage: L.abelian_variety()
            Simple abelian variety J0(23) of dimension 2
        """
        return self._J

    def prime(self):
        r"""
        Return the prime `p` as in 'p-adic L-function'.

        EXAMPLES::

            sage: L = J0(23)[0].padic_lseries(5)
            sage: L.prime()
            5
        """
        return self._p

    def _repr_(self):
        r"""
        Return print representation.

        EXAMPLES::

            sage: from sage.modular.abvar.padic_lseries import pAdicLseries
            sage: J = J0(23)[0]
            sage: pAdicLseries(J, 5)
            5-adic L-series of Simple abelian variety J0(23) of dimension 2
        """
        s = "%s-adic L-series of %s" % (self._p, self._J)
        if not self._normalize == 'L_ratio':
            s += ' (not normalized)'
        return s

    def ap(self):
        """
        Return the Hecke eigenvalue `a_p`.

        EXAMPLES::

            sage: from sage.modular.abvar.padic_lseries import pAdicLseries
            sage: J = J0(23)[0]
            sage: for p in prime_range(5,30):
            ....:     L = pAdicLseries(J, p)
            ....:     p, L.ap()
            (5, 2*alpha)
            (7, 2*alpha + 2)
            (11, -2*alpha - 4)
            (13, 3)
            (17, -2*alpha + 2)
            (19, -2)
            (23, 1)
            (29, -3)
        """
        try:
            A = self._modular_symbols_subspace
        except AttributeError:
            A = self._modular_symbols_subspace = self.modular_symbols_subspace()
        a_p = self._ap = A.eigenvalue(self._p)
        return a_p

    def is_ordinary(self):
        """
        Check if `p` is an ordinary prime.
        """
        try:
            K_f = self._real_quadratic_field
        except AttributeError:
            K_f = self._real_quadratic_field = self.real_quadratic_field()
        try:
            a_p = self._ap
        except AttributeError:
            a_p = self._ap = self.ap()
        frak_p = [x[0] for x in K_f.factor(self._p)]
        not_in_p = [x for x in frak_p if a_p not in frak_p]
        return bool(not_in_p)

    def measure(self, a, n, prec, quadratic_twist=+1):
        r"""
        Return the measure on `\ZZ_p^{\times}` defined by

           `\mu_{J,\alpha}^+ ( a + p^n \ZZ_p  ) =
           \frac{1}{\alpha^n} \left [\frac{a}{p^n}\right]^{+} -
           \frac{1}{\alpha^{n+1}} \left[\frac{a}{p^{n-1}}\right]^{+}`

        where `[\cdot]^{+}` is the modular symbol. This is used to define
        this `p`-adic L-function (at least when the reduction is good).

        The optional argument ``quadratic_twist`` replaces `J` by the twist in the above formula,
        but the twisted modular symbol is computed using a sum over modular symbols of `J`
        rather then finding the modular symbols for the twist.

        Note that the normalisation is not correct at this
        stage: use  ``_quotient_of periods`` and ``_quotient_of periods_to_twist``
        to correct.

        Note also that this function does not check if the condition
        on the ``quadratic_twist=D`` is satisfied. So the result will only
        be correct if for each prime `\ell` dividing `D`, we have
        `ord_{\ell}(N)<= ord_{\ell}(D)`, where `N` is the level

        INPUT:

        -  ``a`` - an integer

        -  ``n`` - a non-negative integer

        -  ``prec`` - an integer

        -  ``quadratic_twist`` (default = 1) - a fundamental discriminant of a quadratic field,
           should be coprime to the level of `J`

        EXAMPLES::


        """
        if quadratic_twist > 0:
            s = +1
        else:
            s = -1
        try:
            p, alpha, z, w, f = self.__measure_data[(n, prec, s)]
        except (KeyError, AttributeError):
            if not hasattr(self, '__measure_data'):
                self.__measure_data = {}
            p = self._p

            alpha = self.alpha(prec=prec)
            z = 1 / (alpha**n)
            w = p**(n - 1)
            f = self.modular_symbol

            self.__measure_data[(n, prec, s)] = (p, alpha, z, w, f)

        if quadratic_twist == 1:
            return z * f(a / (p * w)) - (z / alpha) * f(a / w)
        else:
            D = quadratic_twist
            chip = kronecker_symbol(D, p)
            if self._E.conductor() % p == 0:
                mu = chip**n * z * sum([kronecker_symbol(D, u) *
                                        f(a / (p * w) + ZZ(u) / D)
                                        for u in range(1, abs(D))])
            else:
                mu = chip**n * sum([kronecker_symbol(D, u) *
                                    (z * f(a / (p * w) + ZZ(u) / D) -
                                     chip * (z / alpha) * f(a / w + ZZ(u) / D))
                                    for u in range(1, abs(D))])
            return s * mu

    def alpha(self, prec=20):
        r"""
        Return a `p`-adic root `\alpha` of the polynomial `x^2 - a_p x
        + p` with `ord_p(\alpha) < 1`.  In the ordinary case this is
        just the unit root.

        INPUT:
        -  ``prec`` - positive integer, the `p`-adic precision of the root.

        EXAMPLES:

        """
        try:
            return self._alpha[prec]
        except AttributeError:
            self._alpha = {}
        except KeyError:
            pass

        p = self._p
        Q = Qp(p)
        try:
            a_p = self._ap
        except AttributeError:
            a_p = self._ap = self.ap()
        try:
            psi = self._psi
        except AttributeError:
            psi = self._psi = self.psi()
        K_f = self._real_quadratic_field
        F = Q.extension(K_f.defining_polynomial(), names='a')
        G = K_f.embeddings(K_f)
        if G[0](K_f.gen()) == K_f.gen():
            conj_map = G[1]
        else:
            conj_map = G[0]
        a_p_conj = conj_map(a_p)
        a_p_conj_padic = psi(a_p_conj)
        a_p_padic = psi(a_p)
        R = F['x']
        x = R.gen()
        f = x**2 - a_p_padic * x + p
        fconj = x**2 - a_p_conj_padic * x + p
        norm_f = f * fconj
        norm_f_basefield = norm_f.change_ring(Q)
        FF = norm_f_basefield().factor()
        root0 = -f.gcd(FF[0][0])[0]
        root1 = -f.gcd(FF[1][0])[0]
        if root0.valuation() < 1:
            padic_lseries_alpha = root0
        else:
            padic_lseries_alpha = root1
        return padic_lseries_alpha

    def order_of_vanishing(self):
        r"""
        Return the order of vanishing of this `p`-adic L-series.

        The output of this function is provably correct, due to a
        theorem of Kato [Ka].  This function will terminate if and only if
        the Mazur-Tate-Teitelbaum analogue [MTT]_ of the BSD conjecture about
        the rank of the curve is true and the subgroup of elements of
        `p`-power order in the Shafarevich-Tate group of this curve is
        finite.  I.e. if this function terminates (with no errors!),
        then you may conclude that the `p`-adic BSD rank conjecture is
        true and that the `p`-part of Sha is finite.

        NOTE: currently `p` must be a prime of good ordinary reduction.

        REFERENCES:

        - [MTT]_

        - [Ka]_

        EXAMPLES::

        """
        try:
            return self.__ord
        except AttributeError:
            pass

        if not self.is_ordinary():
            raise NotImplementedError
        E = self.elliptic_curve()
        if not E.is_good(self.prime()):
            raise ValueError("prime must be of good reduction")
        r = E.rank()
        n = 1
        while True:
            f = self.power_series(n)
            v = f.valuation()
            if v < r:
                raise RuntimeError("while computing p-adic order of vanishing, got a contradiction: the curve is %s, the curve has rank %s, but the p-adic L-series vanishes to order <= %s" % (E, r, v))
            if v == r:
                self.__ord = v
                return v
            n += 1

    def _c_bounds(self, n):
        raise NotImplementedError

    def _prec_bounds(self, n, prec):
        raise NotImplementedError

    def teichmuller(self, prec):
        r"""
        Return Teichmuller lifts to the given precision.

        INPUT:

        - ``prec`` - a positive integer.

        OUTPUT:

        - a list of `p`-adic numbers, the cached Teichmuller lifts

        EXAMPLES::

            sage: L = EllipticCurve('11a').padic_lseries(7)
            sage: L.teichmuller(1)
            [0, 1, 2, 3, 4, 5, 6]
            sage: L.teichmuller(2)
            [0, 1, 30, 31, 18, 19, 48]
        """
        p = self._p
        K = Qp(p, prec, print_mode='series')
        return [Integer(0)] + \
               [a.residue(prec).lift() for a in K.teichmuller_system()]

    def _e_bounds(self, n, prec):
        p = self._p
        prec = max(2, prec)
        R = PowerSeriesRing(ZZ, 'T', prec + 1)
        T = R(R.gen(), prec + 1)
        w = (1 + T)**(p**n) - 1
        return [infinity] + [valuation(w[j], p)
                             for j in range(1, min(w.degree() + 1, prec))]

    def _get_series_from_cache(self, n, prec, D):
        try:
            return self.__series[(n, prec, D)]
        except AttributeError:
            self.__series = {}
        except KeyError:
            for _n, _prec, _D in self.__series.keys():
                if _n == n and _D == D and _prec >= prec:
                    return self.__series[(_n, _prec, _D)].add_bigoh(prec)
        return None

    def _set_series_in_cache(self, n, prec, D, f):
        self.__series[(n, prec, D)] = f

    def _quotient_of_periods_to_twist(self, D):
        r"""
        For a fundamental discriminant `D` of a quadratic number field this computes the constant `\eta` such that
        `\sqrt{D}\cdot\Omega_{E_D}^{+} =\eta\cdot \Omega_E^{sign(D)}`. As in [MTT]_ page 40.
        This is either 1 or 2 unless the condition on the twist is not satisfied, e.g. if we are 'twisting back'
        to a semi-stable curve.

        REFERENCES:

        - [MTT]_

        .. note: No check on precision is made, so this may fail for huge `D`.

        EXAMPLES::

        """
        from sage.functions.all import sqrt
        # This function does not depend on p and could be moved out of
        # this file but it is needed only here

        # Note that the number of real components does not change by twisting.
        if D == 1:
            return 1
        if D > 1:
            Et = self._E.quadratic_twist(D)
            qt = Et.period_lattice().basis()[0] / self._E.period_lattice().basis()[0]
            qt *= sqrt(qt.parent()(D))
        else:
            Et = self._E.quadratic_twist(D)
            qt = Et.period_lattice().basis()[0] / self._E.period_lattice().basis()[1].imag()
            qt *= sqrt(qt.parent()(-D))
        verbose('the real approximation is %s' % qt)
        # we know from MTT that the result has a denominator 1
        return QQ(int(round(8 * qt))) / 8


class pAdicLseriesOrdinary(pAdicLseries):
    def series(self, n=2, quadratic_twist=+1, prec=5):
        r"""
        Return the `n`-th approximation to the `p`-adic L-series as
        a power series in `T` (corresponding to `\gamma-1` with
        `\gamma=1+p` as a generator of `1+p\ZZ_p`).  Each
        coefficient is a `p`-adic number whose precision is provably
        correct.

        Here the normalization of the `p`-adic L-series is chosen
        such that `L_p(J,1) = (1-1/\alpha)^2 L(J,1)/\Omega_J`
        where `\alpha` is the unit root

        INPUT:

        -  ``n`` - (default: 2) a positive integer
        -  ``quadratic_twist`` - (default: +1) a fundamental discriminant
           of a quadratic field, coprime to the
           conductor of the curve
        -  ``prec`` - (default: 5) maximal number of terms of the series
           to compute; to compute as many as possible just
           give a very large number for ``prec``; the result will
           still be correct.

        ALIAS: power_series is identical to series.

        EXAMPLES:

            sage: from sage.modular.abvar.padic_lseries import pAdicLseries
            sage: J = J0(188)[0]
            sage: p = 7
            sage: L = pAdicLseriesOrdinary(J, p)
            sage: L.is_ordinary()
            True
            sage: f = L.power_series(2)
            sage: f[0]
            O(7^20)
            sage: f[1].norm()
            3 + 4*7 + 3*7^2 + 6*7^3 + 5*7^4 + 5*7^5 + 6*7^6 + 4*7^7 + 5*7^8 + 7^10 + 5*7^11 + 4*7^13 + 4*7^14 + 5*7^15 + 2*7^16 + 5*7^17 + 7^18 + 7^19 + O(7^20)

        """
        n = ZZ(n)
        if n < 1:
            raise ValueError("n (=%s) must be a positive integer" % n)
        if not self.is_ordinary():
            raise ValueError("p (=%s) must be an ordinary prime" % self._p)
        # check if the conditions on quadratic_twist are satisfied
        D = ZZ(quadratic_twist)
        if D != 1:
            if D % 4 == 0:
                d = D // 4
                if not d.is_squarefree() or d % 4 == 1:
                    raise ValueError("quadratic_twist (=%s) must be a fundamental discriminant of a quadratic field" % D)
            else:
                if not D.is_squarefree() or D % 4 != 1:
                    raise ValueError("quadratic_twist (=%s) must be a fundamental discriminant of a quadratic field" % D)
            if gcd(D, self._p) != 1:
                raise ValueError("quadratic twist (=%s) must be coprime to p (=%s) " % (D, self._p))
            if gcd(D, self._E.conductor()) != 1:
                for ell in prime_divisors(D):
                    if valuation(self._E.conductor(), ell) > valuation(D, ell):
                        raise ValueError("can not twist a curve of conductor (=%s) by the quadratic twist (=%s)." % (self._E.conductor(), D))
        p = self._p
        if p == 2 and self._normalize:
            print('Warning : For p=2 the normalization might not be correct !')
        # verbose("computing L-series for p=%s, n=%s, and prec=%s" % (p,n,prec))

#        bounds = self._prec_bounds(n, prec)
#        padic_prec = max(bounds[1:]) + 5
        padic_prec = 10
#        verbose("using p-adic precision of %s" % padic_prec)

        res_series_prec = min(p**(n - 1), prec)
        verbose("using series precision of %s" % res_series_prec)

        ans = self._get_series_from_cache(n, res_series_prec, D)
        if ans is not None:
            verbose("found series in cache")
            return ans

        K = QQ
        gamma = K(1 + p)
        R = PowerSeriesRing(K, 'T', res_series_prec)
        T = R(R.gen(), res_series_prec)
        L = R(0)
        one_plus_T_factor = R(1)
        gamma_power = K(1)
        teich = self.teichmuller(padic_prec)
        p_power = p**(n - 1)

        verbose("Now iterating over %s summands" % ((p - 1) * p_power))
        verbose_level = get_verbose()
        count_verb = 0
        for j in range(p_power):
            s = K(0)
            if verbose_level >= 2 and j / p_power * 100 > count_verb + 3:
                verbose("%.2f percent done" % (float(j) / p_power * 100))
                count_verb += 3
            for a in range(1, p):
                b = teich[a] * gamma_power
                s += self.measure(b, n, padic_prec, quadratic_twist=D)
            L += s * one_plus_T_factor
            one_plus_T_factor *= 1 + T
            gamma_power *= gamma

#        verbose("the series before adjusting the precision is %s"%L)
        # Now create series but with each coefficient truncated
        # so it is proven correct:
#        K = Qp(p, padic_prec, print_mode='series')
#        R = PowerSeriesRing(K, 'T', res_series_prec)
#        L = R(L, res_series_prec)
#        aj = L.list()
#        if len(aj):
#            aj = [aj[0].add_bigoh(padic_prec-2)] + [aj[j].add_bigoh(bounds[j]) for j in range(1, len(aj))]
#        L = R(aj, res_series_prec )

#        L /= self._quotient_of_periods_to_twist(D) * self._E.real_components()

#        self._set_series_in_cache(n, res_series_prec, D, L)

        return L

    power_series = series
    
#    def _c_bound(self):
#        try:
#            return self.__c_bound
#        except AttributeError:
#            pass
#        E = self._E
#        p = self._p
#        if E.galois_representation().is_irreducible(p):
#            ans = 0
#        else:
#            m = E.modular_symbol_space(sign=1)
#            b = m.boundary_map().codomain()
#            C = b._known_cusps()  # all known, since computed the boundary map
#            ans = max([valuation(self.modular_symbol(a).denominator(), p) for a in C])
#        self.__c_bound = ans
#        return ans

#    def _prec_bounds(self, n, prec):
#        p = self._p
#        e = self._e_bounds(n-1, prec)
#        c = self._c_bound()
#        return [e[j] - c for j in range(len(e))]
