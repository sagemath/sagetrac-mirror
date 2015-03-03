"""
Periods

Computation of period integrals associated to spaces of modular symbols.
"""

#############################################################################
#       Copyright (C) 2012 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

from sage.modular.periods.periods_cython import extended_period_integral
from sage.modular.arithgroup.all import Gamma0
from sage.rings.all import CDF, ZZ, infinity as oo, next_prime
from sage.matrix.all import matrix
import random


def period_integral(m, k, a, b, c, d, v, eps=None, DEBUG=False):
    """
    Given gamma=[a,b;c,d] in Gamma_0(N), and the coefficients v of the
    q-expansion of a modular forms f in S_k(N,eps), compute the period
    integral <f, X^m*Y^(k-2-m){oo,gamma(oo)}>.

    INPUT:

    - m -- an integer
    - k -- an integer
    - a,b,c,d -- integers (a matrix with determinant 1)
    - v -- a list of complex numbers
    - eps -- ?
    - DEBUG -- (default False) whether to print details

    ALGORITHM: See Algorithm 10.6 of Stein's book "Modular Forms, A
    Computational Approach", where here we use c instead of N*c.

    EXAMPLES::

        sage: from sage.modular.periods.periods import period_integral
        sage: period_integral(1,3,4,-3,15,-11,[complex(0,1)])
        0j

        sage: M = ModularSymbols(11,sign=0)[1]
        sage: coefs = [complex(u) for u in M.q_expansion_basis()[0].list()]
        sage: period_integral(1,3,4,-3,15,-11, coefs)
        (0.7647546611823091-0.034437644328508155j)
    """
    X, Y = ZZ['X,Y'].gens()
    I = complex(0, 1)

    if c == 0:
        # easy special case since gamma fixes oo.
        return 0
    (a, b, c, d) = (int(a), int(b), int(c), int(d))

    if not isinstance(v[0], complex):
        raise ValueError('entries of v should be complex (for efficiency)')
    if not m <= k-2:
        raise ValueError("m must be at most k-2")
    if not a*d-b*c == 1:
        raise ValueError("matrix must have determinant 1")

    alpha = (-d + I)/c
    if alpha.imag != 1.0/c:
        (a, b, c, d) = (-a, -b, -c, -d)
        alpha = (-d + I)/c
        assert alpha.imag == 1.0/c

    gamma_alpha = (a+I)/c

    if DEBUG:
        print "alpha = %s, gamma(alpha) = %s"%(alpha, gamma_alpha)

    # Next, compute eps(d)*(gamma^(-1)(X^m*Y^(k-2-m)))*{alpha,oo} - X^m*Y^(k-2-m)*{gamma(alpha), oo}
    eps_d = complex(1 if eps is None else eps(d))
    ans = 0

    # We have gamma^(-1)P(X,Y) = P(a*X+b*Y,c*X+d*Y) = (a*X+b*Y)^m*(c*X+d*Y)^(k-2-m)
    # Using slow generic multivariate polynomial arithmetic is fine here, because
    # this is not the bottleneck.
    Q = (a*X+b*Y)**m * (c*X + d*Y)**(k-2-m)
    for i in range(k-1):
        coeff = Q.coefficient(X**i * Y**(k-2-i))
        if coeff != 0:
            ans += int(coeff) * extended_period_integral(i, alpha, v)  # this is the bottleneck.

    ans *= eps_d
    ans -= extended_period_integral(m, gamma_alpha, v)

    return ans


class PeriodMapping(object):
    """
    The period mapping associated to a space A of modular symbols
    corresponding to a newform, computed to complex double precision
    using prec terms of the q-expansions.

    EXAMPLES::

        sage: M = ModularSymbols(Gamma1(3),weight=13,sign=0).cuspidal_submodule().decomposition()[0]
        sage: f = M.period_mapping(100)
        sage: A = M.ambient(); [f(A([i,0,oo]))[0] for i in range(12)]
        [27.9280923271, 5.33567468252*I, -1.12953617856, -0.271724173647*I, 0.0784400124002, 0.0301915748497*I, -0.0174311138667, -0.0150957874248*I, 0.0174311138667, 0.0241532598797*I, -0.0380315211637, -0.0663551095597*I]
    """
    def __init__(self, M, prec):
        DEBUG = False
        if not M.dimension() > 0:
            raise ValueError("M must have positive dimension")
        if not prec >= 1:
            raise ValueError("prec must be at least 1")
        if not M.is_cuspidal():
            raise ValueError("M must be cuspidal")
        if not M.is_new():
            raise ValueError("M must be new")
        if not M.sign() == 0:
            raise ValueError("M must have sign 0")

        I = complex(0, 1)
        self._M = M
        self._prec = prec
        f = M.q_eigenform(prec, 'a')
        v0 = f.list()
        embeddings = f.base_ring().embeddings(CDF)
        v = self._v = [[complex(e(a)) for a in v0] for e in embeddings]
        phi = M.rational_period_mapping()

        # Find symbol s = X^i*Y^(k-2-i)*{oo,gamma(oo)}, with gamma in Gamma0(N),
        # such that phi(x) != 0, and with lower left corner as small as possible.
        k = M.weight()
        self._A = A = M.ambient()
        star = A.star_involution()
        N = M.level()
        candidates = []
        done = False
        gens = [g.matrix() for g in Gamma0(N).gens()]
        gamma = gens[0]
        while len(candidates) == 0:
            for m in range(k-2+1):
                gamma *= random.choice(gens)**random.choice([-1, 1])
                if gamma[1, 0] < 0:
                    gamma = -gamma
                c = gamma[1, 0]
                if c:
                    gamma_oo = gamma[0, 0] / c
                    elt = A([m, oo, gamma_oo])
                    t = phi(elt)
                    if t:
                        tstar = phi(star(elt))
                        if t != tstar and t != -tstar:   # not in +1 or -1 subspaces for star
                            candidates.append((abs(c), m, gamma.list(), elt))
                            if abs(c) == N:
                                # can't beat this
                                done = True
                                break
                if done: break
        candidates.sort()
        if DEBUG: print "candidates = ", candidates
        _, self._m, self._gamma, self._s = candidates[0]
        a,b,c,d = self._gamma
        if DEBUG:
            print "Using X^%sY^%s*{oo, [%s,%s;%s,%s](oo)} = %s = %s to compute period integrals"%(
                self._m, k-2-self._m, a,b,c,d, self._s.modular_symbol_rep(), self._s)

        # Evaluate periods of this symbol for each conjugate newform
        periods = [period_integral(self._m, k, a, b, c, d, g, M.character())
                   for g in v]
        if DEBUG:
            print "periods = ", periods

        # Find Hecke operator T that acts with irreducible charpoly on M+.
        p = 2
        Mplus = M.plus_submodule()
        Tplus = Mplus.hecke_matrix(p)
        T = M.hecke_matrix(p)
        if DEBUG:
            print "T_2"
        a = f[2]
        while not Tplus.charpoly().is_irreducible():
            p = next_prime(p)
            assert p < 100, "must be a bug"  # TODO

            a += f[p]
            Tplus += Mplus.hecke_matrix(p)
            T += M.hecke_matrix(p)
            if DEBUG:
                print " + T_%s"%p,
        if DEBUG:
            print "acts irreducibly"

        # Make matrix with rows the periods of ((s+*s)/2)*T^i and ((s-*s)/2)*T^i
        rows = [[z.real for z in periods]]
        z = [emb(a) for emb in embeddings]
        for j in range(len(embeddings)-1):
            rows.append([z[i]*rows[-1][i] for i in range(len(embeddings))])
        P_real = matrix(CDF, rows)
        rows = [[I*z.imag for z in periods]]
        z = [emb(a) for emb in embeddings]
        for j in range(len(embeddings)-1):
            rows.append([z[i]*rows[-1][i] for i in range(len(embeddings))])
        P_imag = matrix(rows)
        P = P_real.stack(P_imag)

        # Map we want is first projection to M, followed by writing element
        # in terms of basis ((s+*s)/2)*T^i, ..., ((s-*s)/2)*T^i) then finally
        # multiplying by period_mat.
        pi = M.projection()

        pimat = pi.matrix()

        Tpows_plus = T.iterates(M.free_module().coordinate_vector(pi(self._s+star(self._s)).element()/2),
                           Mplus.dimension(), rows=True)
        Tpows_minus = T.iterates(M.free_module().coordinate_vector(pi(self._s-star(self._s)).element()/2),
                           Mplus.dimension(), rows=True)
        Tpows = Tpows_plus.stack(Tpows_minus)

        # our matrix acting from right is:    pi * Tpows^(-1)  * P

        self._period_matrix = pimat * Tpows**(-1) * P

    def __call__(self, x):
        """
        EXAMPLES::
        """
        return self._A(x).element() * self._period_matrix

    def __repr__(self):
        """
        Return a string representation of self

        EXAMPLES::

            sage: M = ModularSymbols(Gamma1(11),weight=11,sign=0).cuspidal_submodule().decomposition()[0]
            sage: M.period_mapping(99)  # indirect doctest
            Period mapping associated to Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 100 for Gamma_1(11) of weight 11 with sign 0 and over Rational Field computed to precision 99
        """
        return "Period mapping associated to %s computed to precision %s" % (self._M, self._prec)
