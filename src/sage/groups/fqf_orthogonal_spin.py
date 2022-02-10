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

    def __repr__(self, sparse=True):
        r"""
        """
        P = self.parent()
        if self.is_one():
            return "()"
        r = ""
        for p in P._S:
            c = self.component(p)
            if not sparse or (p == -1 and c!=1) or (p != -1 and c != (1,1)):
              r += "%s_%s "%(c,p)
        return r

    def component(self, p):
        r"""
        Return `(det, spin)` of the `p`-component.

        EXAMPLES::

            sage: from sage.groups.fqf_orthogonal_spin import GammaA
            sage: gammaA = GammaA([2,5])
            sage: g = gammaA.an_element()
            sage: g
            (-1, 10)_2 (-1, 10)_5
            sage: g.component(2)
            (-1, 10)
            sage: g.component(5)
            (-1, 10)
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
            spin0 = 3**(exp[2])
            spin0 *= 7**(exp[3])
            spin *= (spin0 % 8)
        elif p == -1:
            return (-1)**exp[-1]
        else:
            up = ZZ(_min_nonsquare(p))
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

            sage: from sage.groups.fqf_orthogonal_spin import GammaA
            sage: gammaA = GammaA([2,5])
    """
    def __init__(self, S):
        r"""
        """
        if -1 in S:
            n = 3*(len(S)-1) + 2
            assert S[-1] == -1
        else:
            n = 3*len(S) + 1
        self._S = S
        # first [+-1] x [val mod 2] x [1,3,5,7] or [+-] the kronecker symbol
        # take 3, 7 as generators of (1,3,5,7)
        # the real place -1 is at the last entry
        AbelianGroupGap.__init__(self, n * (2,))

    @staticmethod
    def __classcall_private__(cls, S):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.groups.fqf_orthogonal_spin import GammaA
            sage: gammaA1 = GammaA([2,5])
            sage: gammaA2 = GammaA([5,2])
            sage: gammaA1 is gammaA2
            True
        """
        if not 2 in S:
            raise ValueError("must contain 2")
        S = [ZZ(s) for s in S]
        S.sort()
        if -1 in S:
          # move -1 to the end of S
          S = S[1:]
          S.append(-1)
        S = tuple(S)
        return super(GammaA, cls).__classcall__(cls, S)

    Element = AdelicSquareClass

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.groups.fqf_orthogonal_spin import GammaA
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
            elif p== -1:
                gammap = self.gens()[-1]
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

            sage: from sage.groups.fqf_orthogonal_spin import GammaA
            sage: G = GammaA((2,3))
            sage: G.embed(-1, 1, 3)
            (-1, 1)_3
            sage: G.embed(1, 5, 2)
            (1, 5)_2
        """
        det = ZZ(det)
        spin = QQ(spin)
        spin = spin.denominator() * spin.numerator()     # same square class
        p = ZZ(p)
        dic = self.gens_by_primes()
        if p == -1:
          sign = (det*spin).sign()
          e = (1 - sign)//2
          g = self.gens()[-1]**e
          return g
        val = spin.valuation(p) % 2
        u = spin.prime_to_m_part(p)
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

            sage: from sage.groups.fqf_orthogonal_spin import GammaA
            sage: G = GammaA((2,3))
            sage: G.embed_diagonal(1, 3)
            (1, 3)_2 (1, 3)_3
            sage: G.embed_diagonal(-1, 3)
            (-1, 3)_2 (-1, 3)_3
        """
        return self.prod([self.embed(det, spin, p) for p in self._S])

    def gammaS(self, sign=False):
        r"""
        Return `\Gamma_S`.

        MATH:

        \Gamma_S = \{(d,s) \in \Gamma_\QQ \mid s \in S \}

        MATH:

        \Gamma_S^+ = ker(\Gamma_\QQ \to \{\pm 1\}, (d,s) \mapsto \sign(ds))

        EXAMPLES::

            sage: from sage.groups.fqf_orthogonal_spin import GammaA
            sage: G = GammaA((2,3))
            sage: G.gammaS()
            [(1, 2)_2 (1, 2)_3 ,
             (1, 3)_2 (1, 3)_3 ,
             (-1, 1)_2 (-1, 1)_3 ,
             (1, 7)_2 (1, 2)_3 ]
        """
        #if sign:
        #    G = self.gammaS(False)
        #    A = AbelianGroupGap([2])
        #    imgs = [A([(1-g.component(-1))//2]) for g in G]
        #    return [self(g) for g in self.subgroup(G).hom(imgs).kernel().gens()]
        gens = [self.embed_diagonal(1, s) for s in self._S]
        if sign:
            gens.append(self.embed_diagonal(-1,-1))
        else:
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

            sage: from sage.groups.fqf_orthogonal_spin import GammaA
            sage: q = IntegralLattice("U").twist(2).discriminant_group()
            sage: GammaA([2]).sigma_sharp(2, -4, q)
            Subgroup of GammaA: (2,) generated by (f3, f4)
        """
        gens = []
        for p in [s for s in self._S if s!=-1]:
            q_p = q.primary_part(p).normal_form().gram_matrix_quadratic()
            gens_p = sigma_sharp(rk, det, q_p, p)
            gens_p = [self.embed(g[0], g[1], p) for g in gens_p]
            gens += gens_p
        return self.subgroup(gens)

def sigma_sharp(rkL, detL, q, p):
    r"""
    Return generators for `\Sigma^\#(L\otimes \ZZ_p)` of a lattice `L`.

    INPUT:

    - ``rkL`` -- an non-negative integer
    - ``detL`` -- an integer
    - ``q`` -- a rational matrix representing a torsion quadratic form in normal
      form
    - ``p`` -- a prime number

    OUTPUT:

    - a list of tuples `(det_p,spin_p)` with det=+-1 and spin in \QQ.

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal_spin import sigma_sharp
        sage: q = matrix.diagonal([2,2,2])/3
        sage: sigma_sharp(4, -4, q, 3)
        [(-1, -64/27)]
    """
    lq = q.ncols()
    delta = detL * q.det()
    up = _min_nonsquare(p)
    if p != 2:
        if rkL == lq:
            return [(1, 1)]
        elif rkL == lq + 1:
            return [(-1, 2*delta)]
        else:
            return [(-1, 1), (1, up)]

    # p=2 case
    if rkL > lq + 1:
        if p == 2:
            return [(-1,1), (1,3), (1,7)]
    from sage.quadratic_forms.genera.normal_form import collect_small_blocks
    blocks = collect_small_blocks(q)
    u = [0, 0, 0, 0]
    v = [0, 0, 0, 0]
    w = [[],[],[],[]]
    for b in blocks:
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
    if len(w[2]) == 2:
        return gamma22
    if len(w[1])<= 2:
        return []
    assert False

def reflection(G, v):
    r"""
    Return the matrix represenation of the orthogonal reflection in `v`.

    INPUT:

    - ``v`` -- a vector
    - ``G`` -- a symmetric matrix

    OUTPUT:

    - a matrix

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal_spin import reflection
        sage: G = matrix.identity(2)
        sage: v = vector([1,0])
        sage: reflection(G, v)
        [-1  0]
        [ 0  1]
    """
    n = v.degree()
    E = v.parent().basis()
    vsq = v*G*v
    ref = []
    for k in range(n):
        tk = E[k] - 2/vsq*(E[k]*G*v)*v
        ref.append(tk)
    return matrix(ref)

def spin_exact(G, f):
    r"""
    INPUT:

        - ``G`` - a diagonal matrix
        - ``f`` - a matrix with `f G f^T = G

    OUTPUT:

        - the spinor norm of ``f``
    """
    if not (G.is_square() and f.is_square() and G == f*G*f.T):
        raise ValueError("invalid input")
    n = G.ncols()
    V = QQ**n
    spinor_norm = QQ(1)
    for i in range(n):
        w = V.gen(i)
        v = w*f
        r = v - w
        s = r*G*r
        if not s.is_zero():
            tau = reflection(G, r)
            f = f * tau
            assert w*f == w
            spinor_norm *= s
        else:
            r1 = v + w
            s1 = (r1*G*r1)/2
            r2 = v
            s2 = (r2*G*r2)/2
            assert not s1.is_zero() and not s2.is_zero()
            tau1 = reflection(G, r1)
            tau2 = reflection(G, r2)
            f = f*tau2*tau1
            assert w*f == w
            spinor_norm *= s1*s2
    assert f.is_one()
    return spinor_norm




def det_spin_p(G, T, p, nu):
    r"""
    Return approximations for  `(det_p, spin_p)` of ``T``.

    The algorithm is due to Shimada.
    We follow the conventions of Miranda and Morrison that the quadratic form is defined by
    Q(x) = (x G x.T)/2. Then the spinor norm of the reflection in x is Q(x).

    INPUT:

        - ``G`` -- a diagonal matrix
        - ``T`` -- an isometry up to some precision
        - ``p`` -- a prime number
        - ``nu`` -- an integer giving the valuation of the approximation
        error of ``T``

    EXAMPLES::
    """
    from sage.groups.fqf_orthogonal_lift import _min_val
    from sage.rings.all import Qp
    def mv(A):
        return _min_val(A.change_ring(Qp(p)))
    if p == 2:
        delta = 1
    else:
        delta = 0
    gammaL = [d.valuation(p) for d in G.diagonal()]
    gamma = min(gammaL)
    l = G.ncols()
    E = G.parent()(1)
    reflection_vectors = []

    k = 0
    while k < l:
        g = T.row(k)
        # error estimates
        lambd = mv(g)
        rho = min(delta + nu + gamma, 2*nu + gamma)
        sigma = min( delta + nu + gamma, delta + nu + lambd, 2*nu + gamma)
        kappa = sigma - gammaL[k] - 2*delta
        if (rho <= gammaL[k] + delta) or (kappa < 1 + 2*delta):
            raise ArithmeticError("Recompute with higher precision") #or a ValueError ?
        bm = g - E.row(k)
        qm = bm * G * bm
        if qm.valuation(p) <= gammaL[k] + 2*delta:
            tau1 = reflection(G, bm)
            reflection_vectors.append(bm)
            tau2 = E
        else:
            bp = g + E.row(k)
            qp = bp * G * bp
            assert qp.valuation(p) <= gammaL[k] + 2*delta
            tau1 = reflection(G, bp)
            tau2 = reflection(G, E.row(k))
            reflection_vectors.append(bp)
            reflection_vectors.append(E.row(k))
        lambdaT = mv(T)
        alpha = mv(tau1)
        beta = mv(tau2)
        theta = gamma + min(kappa + 2*min(0,lambd), nu + min(0,lambd), 2*nu)
        nu = min(nu + alpha, lambdaT + theta - gammaL[k] - delta, nu + theta - gammaL[k] - delta) + beta
        T = T * tau1 * tau2
        k += 1
    err = mv(T-E)
    assert err >= nu
    spinor_norm = QQ.prod([v*G*v/2 for v in reflection_vectors])
    determinant = QQ(-1)**(len(reflection_vectors))
    v, u = spinor_norm.val_unit(p)
    if p == 2:
        u = u % 8
    else:
        u = u % p
    spinor_norm = u * p**(v % 2)
    return determinant, spinor_norm

def det_spin_homomorphism(L, sign=False):
    r"""
    Return the det spin homomorphism.

    This has only sense if ``L`` is the orthogonal group of a
    discriminant group of a lattice.

    Input:

    - ``L`` -- an integral lattice
    - ``sign`` -- (Default: ``False``)

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal_spin import det_spin_homomorphism
        sage: L = IntegralLattice("A2").twist(5)
        sage: L = L.direct_sum(L.twist(-1))
        sage: Oq = L.discriminant_group().orthogonal_group()
        sage: det_spin_homomorphism(L).image(Oq).order()
        8
    """
    T = L.discriminant_group()
    rank = L.rank()
    det = L.determinant()
    Oq = T.orthogonal_group()
    S = (2*Oq.domain().order()).prime_divisors()
    L = T.W()
    Gamma = GammaA(S)
    sigma_sharp = Gamma.sigma_sharp(rank, det, T)
    codom = Gamma.quotient(sigma_sharp)

    from sage.groups.fqf_orthogonal_lift import Hensel_qf
    from sage.rings.all import Zp, QQ
    from sage.quadratic_forms.all import QuadraticForm
    det_spin = dict([[f,[]] for f in Oq.gens()])

    for p in S:
        if p == -1:
            continue
        Tp = T.primary_part(p).normal_form()
        if Tp.cardinality() == 1:
            for f in Oq.gens():
                det_spin[f].append((p, (1, 1)))
            continue
        u = matrix([b.vector() for b in Tp.gens()])
        Op = Tp.orthogonal_group()
        q = Tp.gram_matrix_quadratic()
        # Shimada 5.2 Step 2 page 25
        # The discriminant group has full length and a term 1/2 or 3/2
        # to determine the lattice one needs to know the determinant
        val, unit = (q.det()*det).val_unit(ZZ(p))
        assert val==0
        if p == 2 and len(Tp.invariants())== rank and unit % 8 != 1:
            n = q.ncols()
            for i in range(1,n-1):
                if q[i,i+1]==0 and q[i-1,i]==0 and q[i,i].denominator()==2:
                    q[i,i] *= 5
                    break
            else:
                if q[0,0].denominator()==2 and (n==1 or q[0,1]==0):
                    q[0,0] *= 5
                elif q[-1,-1].denominator()==2 and q[-1,-2]==0:
                    q[-1,-1] *= 5
                else:
                    raise AssertionError('bug in det_spin_homomorphism')
        M = q.inverse()
        # diagonalize
        qf = QuadraticForm(QQ, M)
        diag, t = qf.rational_diagonal_form(return_matrix=True)
        diag = diag.Hessian_matrix()
        t = t.T
        q0 = q*q.denominator()
        v1 = t.denominator().valuation(p)
        v2 = t.inverse().denominator().valuation(p)
        v = -v1 -v2 # lower bound for precision loss due to diagonalization
        # compute the spin of lifts of the generators
        for f in Oq.gens():
            # take only the action on the p-part
            fp = Op(f).matrix()
            prec0 = 1
            prec = 25 # initial precision
            # change to the user basis
            g = u * fp * u.inverse()
            while True:
                R = Zp(p, type='fixed-mod', prec=prec+3, print_mode='terse',
                    show_prec=False, print_pos=False)
                g = Hensel_qf(q0.change_ring(R), g.change_ring(R), prec0, prec)
                g = g.change_ring(ZZ)
                try:
                    gg = t*M*g*M.inverse()*t.inverse()
                    det_p, spin_p = det_spin_p(diag, gg, p, prec + v)
                    det_spin[f].append((p, (det_p, spin_p)))
                    break
                except ArithmeticError:
                    # retry with higher precision
                    prec0 = prec
                    prec = 2 * prec
    imgs = [Gamma.prod([Gamma.embed(det=ds[1][0], spin=ds[1][1], p=ds[0])
                        for ds in det_spin[f]])
            for f in Oq.gens()]
    return Oq.hom([codom(g) for g in imgs],codom,check=False)


def genus_group(L):
    r"""
    Compute the genus group of this lattice.

    EXAMPLES::

        sage: L = IntegralLattice(matrix.diagonal([1,-2,64])*2)
        sage: L.genus_group().order()
        2
        sage: L = IntegralLattice(matrix.diagonal([1,-2,64])*2)
        sage: L.genus_group().order()
        2
    """
    sig = L.signature_pair()
    rk = L.rank()
    if False:#L.rank() < 3 or sig[0]*sig[1] == 0:
        raise ValueError("lattice must be indefinite of rank > 2")
    det = L.determinant()
    q = L.discriminant_group()
    from sage.groups.fqf_orthogonal_spin import GammaA
    S = ZZ(2*L.determinant()).prime_factors()
    G = GammaA(S)
    if L.rank() > 2:
        f = det_spin_homomorphism(L)
        codom = f.codomain()
        Gamma = codom._cover
        gammaS = Gamma.gammaS()
        gens = list(f.image(f.domain()).gens()) + [codom(g) for g in gammaS]
        sub = codom.subgroup(gens)
        return codom.quotient(sub)
