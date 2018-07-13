r"""
Orthogonal groups of torsion quadratic forms.

<Paragraph description>

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Simon Brandhorst (2018-05-15): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap, AbelianGroupElement_gap
from sage.rings.all import ZZ, QQ, Zp
from sage.matrix.all import matrix
from sage.misc.misc_c import prod
from sage.misc.cachefunc import cached_method
from sage.quadratic_forms.genera.normal_form import _min_nonsquare

class AdelicSquareClass(AbelianGroupElement_gap):
    r"""
    Elements of the product over `C_2` x `\QQ_p* / (\QQ_p*)^2` at primes `p`.
    """

    def component(self, p):
        r"""
        Return `(det, spin)` of the `p`-component.

        EXAMPLES::

            sage: from sage.groups.fqf_orthogonal.spinor import GammaA
            sage: gammaA = GammaA([2,5])
            sage: g = gammaA.an_element()
            sage: g
            f1*f2*f3*f4*f5*f6*f7
            sage: g.component(2)
            (-1,2)
            sage: g.component(5)
            (-1, 0)
        """
        p = ZZ(p)
        if not p in self.parent()._S:
            raise ValueError("p is not in S")
        exp = self.exponents()
        det = ZZ(1)
        spin = ZZ(1)
        if p == 2:
            det = (-1)**exp[0]
            spin *= 2**exp[1]
            spin *= 3**(exp[2])
            spin *= 7**(exp[3])
            spin = spin % 8
        else:
            up = _min_nonsquare(p)
            k = 1 + 3*self.parent()._S.index(p)
            det = (-1)**exp[k]
            spin *= p**exp[k+1]
            spin *= up**exp[k+2]
        return (det, spin)

class GammaA(AbelianGroupGap):
    r"""

    - ``S`` -- a finite iterable of primes containing 2 in ascending order

    MATH:

    \Gamma_A = \prod_{p \in S} \Gamma_p

    EXAMPLES::

            sage: from sage.groups.fqf_orthogonal.spinor import GammaA
            sage: gammaA = GammaA([2,5])
    """
    def __init__(self, S):
        r"""
        Constructor.

        TESTS::

            sage: from sage.groups.fqf_orthogonal.spinor import GammaA
            sage: gammaA = GammaA([2,5])
            sage: TestSuite(gammaA).run()
        """
        n = 3*len(S) + 1
        self._S = S
        # first [+-1] x [val mod 2] x [1,3,5,7] or [+-] the kronecker symbol
        # take 3, 7 as generators of (1,3,5,7)
        AbelianGroupGap.__init__(self, n * (2,))

    @staticmethod
    def __classcall_private__(cls, S):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.groups.fqf_orthogonal.spinor import GammaA
            sage: gammaA1 = GammaA([2,5])
            sage: gammaA2 = GammaA([5,2])
            sage: gammaA1 is gammaA2
            True
        """
        if not 2 in S:
            raise ValueError("must contain 2")
        S = [ZZ(s) for s in S]
        S.sort()
        S = tuple(S)
        return super(GammaA, cls).__classcall__(cls, S)

    Element = AdelicSquareClass

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.groups.fqf_orthogonal.spinor import GammaA
            sage: gammaA = GammaA([2,3])
        """
        return "GammaA: " + str(self._S)


    @cached_method
    def gens_by_primes(self):
        r"""
        Return a dictionary which gives the `p` coordinate generators by `p`.
        """
        gamma = []
        for p in self._S:
            if p == 2:
                gammap = self.gens()[:4]
            else:
                k = self._S.index(p)
                gammap = self.gens()[1 + 3*k : 1 + 3*(k+1)]
            gamma.append([p, gammap])
        gamma = dict(gamma)
        return gamma

    def embed(self, det, spin, p):
        r"""
        Return the element of ``self`` corresponding to the input.

        The element is trivial in every coordinate except at `p` where its
        value is given by `(det_p, spin_p)`.

        INPUT:

        - ``det`` -- `1` or `-1`
        - ``spin`` -- a rational number representing a `p`-adic square class
        - ``p`` -- a prime number

        EXAMPLES::

            sage: from sage.groups.fqf_orthogonal.spinor import GammaA
            sage: G = GammaA((2,3))
            sage: G.embed(-1, 1, 3)
            f5
            sage: G.embed(1, 5, 2)
            f3*f4
        """
        det = ZZ(det)
        spin = QQ(spin)
        spin = spin.denominator() * spin.numerator()     # same square class
        p = ZZ(p)
        val = spin.valuation(p) % 2
        u = spin.prime_to_m_part(p)
        dic = self.gens_by_primes()
        g = dic[p][0]**((1 - det) // 2)
        g*= dic[p][1]**val
        if p == 2:
            u = u % 8
            if u == 3:
                g*= dic[p][2]
            if u == 5:
                g*= dic[p][2] * dic[p][3]
            if u == 7:
                g*= dic[p][3]
        else:
            u = ZZ(u).kronecker(p)
            val = val % 2
            g*= dic[p][2]**((1 - u) // 2)
        return g

    def embed_diagonal(self, det, spin):
        r"""
        Return the image of a rational number under the diagonal embedding.

        EXAMPLES::

            sage: from sage.groups.fqf_orthogonal.spinor import GammaA
            sage: G = GammaA((2,3))
            sage: G.embed_diagonal(1, 3)
            f3*f6
            sage: G.embed_diagonal(-1, 3)
            f1*f3*f5*f6
        """
        return self.prod([self.embed(det, spin, p) for p in self._S])

    def gammaS(self):
        r"""
        Return `\Gamma_S`.

        MATH:

        \Gamma_S = \{(d,s) \in \Gamma_\QQ \mid s \in S \}

        EXAMPLES::

            sage: from sage.groups.fqf_orthogonal.spinor import GammaA
            sage: G = GammaA((2,3))
            sage: G.gammaS()
            [f2*f7, f3*f6, f1*f5, f4*f7]
        """
        gens = [self.embed_diagonal(1, s) for s in self._S]
        gens.append(self.embed_diagonal(-1, 1))
        gens.append(self.embed_diagonal(1, -1))
        return gens

    def sigma_sharp(self, rk, det, q):
        r"""
        Return sigma sharp.

        This is the the image of K under `(det, spin)` where `K` is the kernel
        of O(L) --> O(L^v / L)

        INPUT:

        - ``rkL`` -- an non-negative integer

        - ``detL`` -- an integer

        - ``q`` --a torsion quadratic form

        EXAMPLES::

            sage: from sage.groups.fqf_orthogonal.spinor import GammaA
            sage: q = IntegralLattice("U").twist(2).discriminant_group()
            sage: GammaA([2]).sigma_sharp(2, -4, q)
            Subgroup of GammaA: (2,) generated by (f3, f4)
        """
        gens = []
        for p in self._S:
            q_p = q.primary_part(p).normal_form().gram_matrix_quadratic()
            gens_p = sigma_sharp(rk, det, q_p, p)
            gens_p = [self.embed(g[0], g[1], p) for g in gens_p]
            gens += gens_p
        return self.subgroup(gens)

def sigma_sharp(rkL, detL, q, p):
    r"""
    Return `\Sigma^\#(L\otimes \ZZ_p)` of a lattice `L`.

    INPUT:

    - ``rkL`` -- an non-negative integer
    - ``detL`` -- an integer
    - ``q`` -- a rational matrix representing a torsion quadratic form in normal
      form
    - ``p`` -- a prime number

    OUTPUT:

    - a list of tuples `(det_p,spin_p)` with det=+-1 and spin in \QQ.

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.spinor import sigma_sharp
        sage: q = matrix.diagonal([2,2,2])/3
        sage: sigma_sharp(4, -4, q, 3)
        [(-1, -64/27)]
    """
    lq = q.ncols()
    delta = detL * q.det()
    up = _min_nonsquare(p)
    if p != 2:
        if rkL == lq:
            return []
        elif rkL == lq + 1:
            return [(-1, 2*delta)]
        else:
            return [(-1, 1), (1, up), (1, p)]
    if rkL > lq + 1:
        if p == 2:
            return [(1,3), (1,7), (1,2), (-1,1)]
        else:
            return [(1,up), (1, p), (-1,1)]
    from sage.quadratic_forms.genera.normal_form import collect_small_blocks
    blocks = collect_small_blocks(q)
    u = [0, 0, 0, 0]
    v = [0, 0, 0, 0]
    w = [[],[],[],[]]
    for b in blocks[:3]:
        k = b.denominator().valuation(2)
        r = b.ncols()
        if k <= 3:
            if r == 2:
                if b[0,0] == b[1,1] == 0:
                    u[k] += 1
                else:
                    v[k] += 1
            else:
                w[k].append(b[0,0].numerator())
    gamma20 = [(1,3), (1,7), (-1, 1)]
    gamma21 = [(1,3), (1,7)]
    gamma22 = [(1,5)]

    if rkL > lq:
        return gamma20
    if u[1] + v[1] > 0:
        if len(w[1]) > 0:
            return gamma20
        return gamma21
    if len(w[1]) == 1:
        epsilon = w[1][0].numerator()
        delta = (detL*ZZ(2)**(-detL.valuation(2))) / prod([q.det().numerator() for q in blocks if q.denominator().valuation(2) > 1])
        if u[2] + v[2] + len(w[3]) > 0:
            if len(w[2]) > 0:
                return gamma20
            return [(1,5), (-1,delta)]
        if len(w[2]) == 2:
            return gamma20
        if len(w[2]) == 1:
            phi = w[2][0].numerator()
            if epsilon*phi % 4 == 3:
                return [(1,7), (-1,delta)]
            if epsilon*phi % 4 == 1:
                return [(1,3), (-1,delta)]
        return [(-1,delta)]
    if len(w[1]) == 2:
        epsilon = w[1][0].numerator()
        phi = w[1][1].numerator()
        if epsilon*phi % 4 == 3:
            return gamma20
        if len(w[2]) > 0:
            return gamma20
        return [(1,5), (-1,epsilon)]
    assert len(w[1]) == 0
    if u[2] + v[2] > 0:
        return gamma22
    if w[2] == 2:
        return gamma22
    return []
