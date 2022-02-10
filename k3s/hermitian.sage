# Description of genera
# of hermitian lattices
#my_magma = Magma()
#my_magma.set_server_and_command(command="magma.avx64.exe")

from sage.env import SAGE_ROOT
magma.attach_spec(SAGE_ROOT+"/k3s/unit_quotients/lat.spec")

from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap

# G.ChangeRing(E).sage()
#E0.<a> = CyclotomicField(4)
#K0.<b>, phi0 = E0.subfield(a + a**-1)
#p0 = K0.primes_above(2)[0]
#Erel0.<a0> = E0.relativize(phi0)

my_gram_str = """my_gram:=function(L)
    E:=BaseRing(L);
    P:= PseudoBasis(Module(L));
    PP:= [ L | ];
    for p in P do
        ok, x:= IsPrincipal(p[1]); assert ok;
        Append(~PP, x * p[2]);
    end for;
    return Matrix(E, [[ InnerProduct(x,y): y in PP] : x in PP ] );
    end function;"""

magma.eval(my_gram_str)

# if the trace matrix is integral then the scale of the hermitian lattice
# is divisible by n/Different(E) via Höppner Lemma 1.28
#

def inverseNorm(phi,d):
    r"""
    Return all ideals of norm ``d`` invariant under the involution.

    INPUT:

    - ``phi`` -- an embedding of fields `K` into `E` of degree `2`

    - ``d`` -- an integer

    OUTPUT:

    - a list of ideals of `E`

    EXAMPLES::


    """
    d = ZZ(d)
    E = phi.codomain()
    K = phi.domain()
    Erel = E.relativize(phi,'r')
    D = Erel.relative_different()
    D = E.ideal(D)
    ideals = []
    primes = []
    for p in d.prime_factors():
        v = d.valuation(p)
        for P in K.ideal(p).prime_factors():
            if not E.ideal(P).is_coprime(D):
                P = E.ideal(P).prime_factors()[0]
            else:
                P = E.ideal(P)
            nv = P.norm().valuation(p)
            primes.append([P**e for e in range((v // nv) +1)])
    from sage.misc.mrange import cantor_product
    ideal_group = D.parent()
    for I in cantor_product(*primes):
        I = ideal_group.prod(I)
        if I.norm() == d:
            ideals.append(I)
    return ideals


def all_lattice_with_isometry(n, rk, det, max_scale, min_scale=1, min_norm=1, signatures=None, even=True,return_genera=False):
    r"""
    Return all lattices with an isometry of order `n` with minimal polynomial
    the `n`-th cyclotomic polynomial.

    INPUT:

    - ``n`` -- an integer at least `3`

    - ``rk`` -- an integer

    - ``det`` -- an integer

    - ``max_scale`` -- an integer

    - ``signatures`` -- a list of non-negative integers of length the degree of
      `K` and sum at most the rank. The `k`-th entry is the number of negative
      squares in the completion at the `k`-th real embedding

    EXAMPLES::


    """
    n = ZZ(n)
    rk = ZZ(rk)
    det = ZZ(det)

    if n < 1:
        raise ValueError("n must be at least 1")
    det = ZZ(det)
    max_scale = ZZ(max_scale)

    if rk == 0:
        if return_genera:
            return [Genus(matrix([]))]
        L = IntegralLattice(matrix.zero(0))
        return [LatticeWithIsometry(L,matrix.zero(0),order=n)]

    if n<3:
        if signatures is None:
            pos = 2
            neg = rk - 2
        else:
            neg = signatures[0]
            pos = rk - neg
        genera = all_genera_by_det((pos,neg), det, max_scale,even=even)
        genera = [g for g in genera if min_scale.divides(g.scale()) and min_norm.divides(g.norm())]
        if return_genera:
            # if n==1: raise AssertionError
            # the output does not know that f is trivial ... mhh
            return genera
        f = matrix.identity(ZZ, rk)
        if n == 2:
            f = -f
        rep = [LatticeWithIsometry(IntegralLattice(g),f,order=n) for g in flatten([genus.representatives() for genus in genera])]
        return rep

    E.<a> = CyclotomicField(n)
    K.<b>, phi = E.subfield(a + a**-1)
    Erel.<a0> = E.relativize(phi)
    D = E.different()

    out = []
    max_scaleE = max_scale*n / D
    min_scale = n / D * min_scale
    norm_detE = det*(n/D).norm()**rk

    detE = inverseNorm(phi,norm_detE)
    genera = []
    if signatures is None:
        signatures = [rk-1]+[rk]*(K.degree()-1)
    elif len(signatures) != K.degree():
        raise ValueError("signatures must be a list of length %s" %K.degree())
    elif not all([s in ZZ and rk >= s >= 0 for s in signatures]):
        raise ValueError("signatures must be a list of non-negative integers")
    for d in detE:
        genera += all_hermitian_genera_by_det(rk,signatures,d,phi,max_scaleE)
    genera = [g for g in genera if min_scale.divides(g.scale())]


    if return_genera:
        return genera
    for g in genera:
        for gramE in g.representatives():
            gramZ, iso = trace_lattice(gramE)
            assert gramZ.det() == det
            L = IntegralLattice(gramZ)
            out.append(LatticeWithIsometry(L, iso, order=n,gramE=gramE, magmaRep=g.representative()))
    if even:
        min_norm = lcm(min_norm, 2)
    if min_norm != 1:
        out = [g for g in out if min_norm.divides(gcd(g.L.gram_matrix().diagonal()))]
    return out


def trace_lattice(G, order=2):
    r"""
    Return the trace lattice, the isometry of ``G``

    INPUT:

    An hermitian lattice over a cyclotomic field.

    OUTPUT:

    The trace lattice `1/n Tr \circ h`.
    The isometry corresponding to multiplication by an `n`-th root of unity.
    """
    n = G.ncols()
    E = G.base_ring()
    if E==ZZ or E==QQ:
        if order==1:
            iso = matrix.identity(n)
        elif order==2:
            iso = -matrix.identity(n)
        else:
            raise ValueError()
        return G, iso
    a = E.gen(0)
    iso = matrix.block_diagonal(n * [matrix.companion(E.polynomial())])
    d = E.degree()
    a = E.gen()
    coeff = []
    for i in range(n):
        for ia in range(d):
            for j in range(n):
               for ja in range(d):
                   g = a**(ia-ja)*G[i,j]
                   coeff.append(g.trace())
    gram = matrix(QQ,d*n,coeff)
    iso = iso.T
    assert iso*gram*iso.T == gram
    return gram/a.multiplicative_order(), iso

def eichler(G, u, v, y, mu):
    n = G.ncols()
    d = G.base_ring().degree()
    V = G.base_ring()^n
    uv = (u* G * v.conjugate())
    assert u*G*u.conjugate()==0
    assert u*G*y.conjugate()==0
    assert v*G*y.conjugate()==0
    assert (mu*uv)+(mu*uv).conjugate()==-y*G*y.conjugate()
    assert u*G*u.conjugate()==0
    imgs = [e + (e* G * u.conjugate())/uv.conjugate()*y + (mu*(e*G*u.conjugate())/uv.conjugate() - (e*G*y.conjugate())/uv) * u for e in V.gens()]
    f = matrix(imgs)
    assert G == f*G*f.T.conjugate()
    return restriction_of_scalars(f),G.base_ring().one(),f

def restriction_of_scalars(g):
    n = g.ncols()
    d = g.base_ring().degree()
    F = matrix.zero(QQ,d*n,d*n)
    for i in range(n):
        for j in range(n):
            F[d*i:(i+1)*d,d*j:(j+1)*d] = g[i,j].matrix()
    return F


def symmetry(G, s, sigma):
    # x --> x - h(x,s)sigma^-1 s
    assert s * G * s.conjugate() == sigma + sigma.conjugate()
    n = G.ncols()
    d = G.base_ring().degree()
    V = G.base_ring()^n
    imgs = [e - (e* G * s.conjugate())* sigma^-1 * s for e in V.gens()]
    g = matrix(imgs)
    return restriction_of_scalars(g),g.det(),g

def on_disc(D,g):
    imgs = []
    denom = g.denominator()
    m = D.cardinality()
    for d in D.gens():
        x = d.lift()*g
        x = denom.inverse_mod(m)*denom* x
        imgs.append(D(x))
    return D.hom(imgs)

def to_herm(E, n):
    # transformation matrix which inverts the restriction
    # of scalars E^n/QQ
    d = E.degree()
    k = n//d
    b = matrix([E.gen()^i for i in range(d)]).T
    return matrix.block_diagonal(k*[b])

def skew_element(E):
    if E.degree()>2:
      K,phi = E.maximal_totally_real_subfield()
      EK = E.relativize(phi,"a")
      a = EK.gen()
      to_E,from_E = EK.structure()
    else:
      K = QQ
      EK = E
      a = E.gen()
      to_E = E.hom(E)
    g = matrix(K,2,1,[2,E.gen()+E.gen().conjugate()]).kernel().gen()
    omega = g[0]+a*g[1]
    omega = to_E(omega)
    assert omega+omega.conjugate()==0
    return omega

def vecZtoV(V, x):
  n = x.degree()
  E = V.base_ring()
  d = E.degree()
  assert n % d ==0
  x = x.list()
  return V([E(x[i*d:(i+1)*d]) for i in range(n/d)])
end


def Oq_equiv_herm(L,f, G, Lh):
    """
    INPUT:

    - ``L,f`` -- the trace lattice of ``G`` with ``L`` as an IntegralLattice.
      The order of ``f`` must be at least `3`
    - gram matrix over E
    - hermitian lattice in the genus of Lh


    OUTPUT:

    The image of `O(L,f) \to O(L^\vee/L,\bar f)`.
    """
    ord = f.change_ring(ZZ).multiplicative_order()
    n = G.ncols()
    E = G.base_ring()
    prime = false
    if len(prime_factors(ord))==1:
      prime = true
      A = E.different()
      P = A.prime_factors()[0]
    assert ord>2
    assert G.ncols()>1
    OGL = L.q().O()
    OL = L.orthogonal_group([f])
    GL = OGL.gap().Centralizer(OGL.subgroup(OL.gens()))
    # the ambient group we want to generate
    Oqf = OGL.subgroup(GL.GeneratorsOfGroup())
    K,phi = E.maximal_totally_real_subfield()
    VO = E.maximal_order()^n
    VE = E^n
    D = L.discriminant_group()
    #gens_mat = []
    #gens_matE = []
    #gens_mat.append(f)
    #gens_matE.append(E.gen()*matrix.identity(E,n))
    gens = [Oqf(on_disc(D,f)), Oqf(on_disc(D,-f^0))]
    dets = [E.gen()^n, E(-1)^n]
    S = Oqf.subgroup(gens)
    omega = skew_element(E)
    # create random reflections and eichler isometries until we generate U(L^vee/L)
    # should work most of the time, though not always.
    count = 0
    De = enumerate(D)
    flag = true

    m = pari(L.gram_matrix())
    from sage.env import SAGE_EXTCODE
    gp.read(SAGE_EXTCODE + "/pari/simon/qfsolve.gp")
    m = gp.eval('qflllgram_indefgoon(%s)'%m)
    # convert the output string to sage
    extra = pari(m).sage()[1].columns()
    while S.order()!=Oqf.order():
        count +=1
        if flag:
          try:
            x = next(De)[1]
            s = vecZtoV(VE, x.lift()*x.order())
          except StopIteration:
            print("done discr")
            flag = false
        elif len(extra)>0:
            e = extra.pop()
            s = vecZtoV(VE, e)
        else:
            s = VE(VO.random_element().list())
        if s == 0:
            continue
        sigma = s*G*s.conjugate()/2
        print("Computing spinor norms of Oqf. Total:" +str(Oqf.order())+" Remaining:" +str(Oqf.order()/S.order()) + " number of tries: " +str(count), end="\r")
        sys.stdout.flush()
        if 2.divides(D.cardinality()) and sigma == 0:
            # we might run into trouble if V is locally isotropic at 2 but not globally isotropic, ... so lets see if we are lucky...
            sg = (s*G)
            (m,i) = min((e[1].norm(),e[0]) for e in enumerate(sg) if e[1]!=0)
            v = VE.gen(i)
            ker = (G*matrix([s,v]).T.conjugate()).left_kernel()
            for y in ker.basis():
                y = VE(y)
                mu = -y*G*y.conjugate()/2/(s*G*v.conjugate())
                T,det,TE = eichler(G,s,v,y,mu)
                if gcd(T.denominator(),D.cardinality())==1:
                    Tbar = on_disc(D,T)
                    Tbar = Oqf(Tbar)
                    if not Tbar in S:
                        #gens_mat.append(T)
                        #gens_matE.append(TE)
                        gens.append(Tbar)
                        dets.append(E(1))
                        S = Oqf.subgroup(gens)
        if prime:
            I = A*E.ideal((s*G).list())
            if I.valuation(P)<=0:
                continue
        for i in range(100):
          # we want
          # <s,L> <= (sigma+x*omega) O < E.different()^-1
          sigma1 = sigma + phi(K.random_element())*omega
          if sigma1==0:
              continue
          tau,determ,tauE = symmetry(G,s,sigma1)
          if gcd(tau.denominator(),D.cardinality())!=1:
            continue
          taubar = on_disc(D,tau)
          taubar = Oqf(taubar)
          if not taubar in S:
              #gens_mat.append(tau)
              #gens_matE.append(tauE)
              gens.append(taubar)
              dets.append(determ)
              S = Oqf.subgroup(gens)
    Lhn = magma.HermitianLattice((1/ord)*Lh.GramMatrix())
    Fm,invariants,toF = Lhn.get_the_group(nvals=3)
    Fs = AbelianGroupGap(invariants)
    print("\n")
    imgs = [Fs(toF(g).sage()) for g in dets]
    #return S, imgs, Fs, dets, gens_mat,gens_matE
    delta = S.hom(imgs, codomain=Fs)
    return delta.kernel()

def Ei(phi,p):
    r"""
    p a prime ideal of K
    """
    E = phi.codomain()
    K = phi.domain()
    pE = E.fractional_ideal([phi(e) for e in p.gens()])
    F = pE.factor()
    if E.degree()>2:
        Erel = E.relativize(phi,"a")
        toE,fromE = Erel.structure()
    if len(F)==2:
        # split
        pass
    elif len(F)==1 and F[0][1]==1:
        # inert
        pass
    elif len(F)==1 and F[0][1]>1:
        # ramified
        P = F[0][0]
        if E.degree()>2:
            D = Erel.relative_different()
            D = E.fractional_ideal([toE(e) for e in D.gens()])
        else:
            D = E.different()
        e = D.valuation(P)

    return True



#def to_magma(G):
    #r"""
    #"""
    #Gm = magma(G)
    #Em = Gm.BaseRing()
    #Km = Em.sub(Em.1 + Em.1**-1)
    #E = magma.RelativeField(Km, Em)
    #return Gm.ChangeRing(E)

class GenusHermitian(object):
    r"""

    a list of pairs
    [[real_embeddings],[numbers of negative squares]]
    """
    def __init__(self, local_symbols, phi, rank, signatures):
        r"""
        """
        self._local_symbols = local_symbols
        self._phi = phi
        self._rank = rank
        self._signatures = signatures
        assert [s>= 0 for s in signatures]

        self._Erel = phi.codomain().relativize(phi, 'a')
        self._ramified_primes = self._Erel.relative_discriminant().prime_factors()

        P = len(self.non_norm_primes())
        I = len([s for s in self._signatures if s % 2 == 1])
        if (P + I) % 2 == 1:
            raise ValueError("Invariants violate the product formula.")

    def __eq__(self, other):
        r"""
        """
        s1 = self._signatures
        s2 = other._signatures
        if s1 != s2:
            return False

        if len(self._local_symbols) != len(other._local_symbols):
            return False

        for s in s1._local_symbols:
            if not s in s2._local_symbols:
                return False
        return True

    def __neq__(self, other):
        return not (self==other)

    def non_norm_primes(self):
        r"""
        """
        return [s._prime for s in self._local_symbols if not s._is_local_norm(s.det())]

    def representative(self, rational=False):
        r"""
        Return a representative of this genus.

        Uses Markus Kirschmers Magma programs.

        """
        m = magma
        Em = m(self._phi.codomain())
        d = self._Erel.absolute_degree()
        if d == 2:
          Km = m.RationalsAsNumberField()
          m.IsSubfield(Km,Em)
        else:
          Km = Em.sub(Em.1 + Em.1**-1)
        E = m.RelativeField(Km, Em)
        P = self.non_norm_primes()
        O = Km.Integers()
        #print(Km, O)
        # convert to magma
        Pm = m([self._ideal_to_magma(Km, p) for p in P])
        V = m.HermitianFormWithInvariants(E, self._rank, Pm, self._signatures)
        if rational:
            return V.ChangeRing(Em).sage()
        M = V.MaximalIntegralHermitianLattice()
        #M = 27*M
        #M = M.MaximalIntegralLattice()
        for sym in self._local_symbols:
            p = sym._prime
            g = sym.gram_matrix()
            g = m(g).ChangeRing(E)
            L = m.HermitianLattice(g)
            pm = self._ideal_to_magma(Km, p)
            if d == 2:
                pm = gcd([ZZ(a) for a in p.gens()])
            M = m.FindLattice(M, L, pm)
        return M


    def representatives(self):
        r"""
        Return a list of representatives of this genus.

        Uses Markus Kirschmers Magma programs.

        """
        E = self._phi.codomain()
        if self._rank == 1:
            reps = [self.representative()]
        else:
            reps = self.representative().GenusRepresentatives()
        output = []
        for rep in reps:
            rep = rep.my_gram().ChangeRing(E).sage()
            conv = rep.base_ring().hom(E.gens())
            rep = rep.apply_morphism(conv)
            output.append(rep)
        return output


    def _ideal_to_magma(self,Km, I):
        r"""
        """
        O = Km.Integers()
        P = 0 * O
        for g in I.gens():
            P += Km(g.list())*O
        return P

    def volume(self):
        return product([s._prime**s.gram_matrix().determinant().valuation(s._prime) for s in self._local_symbols])

    def norm(self):
        r"""
        """
        return prod([s.norm() for s in self._local_symbols])

    def scale(self):
        r"""
        """
        return prod([s.scale() for s in self._local_symbols])

    def trace(self):
        r"""
        Return the genus of the trace lattice of self.

        This should be faster than first computing a representative.
        """
        # we are cheesy here
        gram = self.representative().my_gram()
        E = self._phi.codomain()
        gram = gram.ChangeRing(E).sage()
        gram, iso = trace_lattice(gram)
        return Genus(gram)
        # some code that won't work
        from sage.quadratic_forms.genera.genus import LocalGenusSymbol,GenusSymbol_global_ring,Genus_Symbol_p_adic_ring

        rank = self._Erel.absolute_degree() * self._rank
        s_minus = 2*sum(self._signatures)
        signature = (rank - s_minus, s_minus)

        local_symbols = []
        if not any(s._is_dyadic for s in self._local_symbols):
            d = self.volume().norm() % 8
            gp = Genus_Symbol_p_adic_ring(2,[[0,rank,d,0,0]])
            local_symbols.append(gp)
        for sym in self._local_symbols:
            p = ZZ(sym._prime.norm()).prime_factors()[0]
            if not len(self._Erel.primes_above(p))==1:
                print('this may be wrong')
            gp,_ = trace_lattice(sym.gram_matrix())
            gp = LocalGenusSymbol(gp,p)
            local_symbols.append(gp)
        return GenusSymbol_global_ring(signature,local_symbols)


class LocalGenusHermitian(object):
    r"""
    The genus of an hermitian lattice over a complete discrete valuation ring.

    Let `K` be a number field, `p` be a prime ideal in `K`, `E` a degree `2`
    field extension of `K` with involution.

    This class describes the genus of a hermitian lattice over the completion
    of `K` at `p`.

    The format for a block over a non-dyadic prime is
    ``(scale, rank, det)``
    and over a dyadic prime it is
    ``(scale, rank, det, norm)``

    Let  and `P` be the largest ideal of
    `O_E` that contains `p` and is invariant under the involution.
    Then the scale of this genus is `P**scale`.

    A description of the jordan constituents ....

    INPUT:

    - ``p`` - a prime ideal in `K`

    - ``symbol`` -- a list of modular symbols

    - ``phi`` -- the field embedding of `K` into `E`

    """
    def __init__(self, prime, symbol, phi):
        r"""
        """
        self._symbol = tuple(symbol)
        self._prime = prime
        self._phi = phi
        E = phi.codomain()
        K = phi.domain()

#        for g in self._prime.gens():
#            if g.valuation(self._prime) == 1:
#                self._uniformizer = g
#                break


        f = phi(self._prime).factor()
        self._is_ramified = (len(f) == 1 and f[0][1]==2)
        self._is_dyadic = ZZ(2).divides(self._prime.norm())

        if self._is_ramified:
            P = E.primes_above(phi(self._prime))
            assert len(P) == 1
            P = P[0]
            self._pi = E.uniformizer(P)
            for u in K.unit_group().gens_values():
                if not self._is_local_norm(u):
                    self._u = u
                    break
            else:
                raise AssertionError("mhh")
            self._uniformizer = phi.domain()(self._pi * self._pi.conjugate())
            assert self._is_local_norm(self._uniformizer)
        else:
            self._uniformizer = prime.ring().uniformizer(prime)

        if self._is_ramified and self._is_dyadic:
            Erel = E.relativize(phi,'n')
            D = Erel.relative_different()
            self._e = D.valuation(self._pi)

            for u in K.unit_group().gens_values():
                if (not self._is_local_norm(u)
                    and (u - 1).valuation(self._prime) == self._e - 1):
                    self._u = u
                    break
            else:
                raise AssertionError("mhh")


    def __repr__(self):
        r"""
        """
        return str(self._symbol)

    def __eq__(self, other  ):
        r"""
        """
        if self._prime != other._prime:
            return False
        sym1 = self._symbol
        sym2 = other._symbol

        if len(sym1) != len(sym2):
            return False
        n = len(sym1)
        # scales agree
        if any(sym1[k][0]!=sym2[k][0] for k in range(n)):
            return False
        # ranks
        if any(sym1[k][1]!=sym2[k][1] for k in range(n)):
            return False
        if not self._is_dyadic:
            # determinant classes
            return all(sym1[k][1]!=sym2[k][1] for k in range(n))

        # dyadic case
        det1 = prod(sym1[k][2] for k in range(n))
        det2 = prod(sym2[k][2] for k in range(n))
        if det1 != det2:
            return False

        # same norms ???? check this
        if any(sym1[k][3] != sym2[k][3] for k in range(n)):
            return False
        e = self._e
        for i in range(n-1):
            u1 = prod(sym1[k][2] for k in range(i+1))
            u2 = prod(sym2[k][2] for k in range(i+1))
            u = u1*u2

            if u != 1:
                A = sym1[i][3] + sym2[i+1][3] - sym1[i][0]
                if A > e - 1:
                    return False
        return True

    def __neq__(self, other):
        return not (self == other)

    def _is_local_norm(self, x):
        r"""
        Return whether `x` is a norm of `E` i.e. in `N(E)`.

        INPUT:

        an element of `K`
        """
        K = self._phi.domain()
        E = self._phi.codomain()
        x = K(x)
        b = K.gen()
        assert self._phi(b) == E.0 + 1/E.0
        hs = K.hilbert_symbol(x, b**2 - 4, self._prime)
        return hs == ZZ(1)

    def det(self):
        r"""
        Return the norm class of the determinant.

        OUTPUT:

        An element of `K`.
        """
        d = ZZ.prod([b[2] for b in self._symbol])

        v = sum([b[0]*b[1] for b in self._symbol])
        if self._is_ramified:
            v = v // 2
        if d == 1:
            u = 1
        else:
            u = self._u
        return u*self._uniformizer**v

    def rank(self):
        r"""
        Return the rank of this symbol.
        """
        return ZZ.sum(b[1] for b in self._symbol)

    def scale(self):
        r"""
        Return the scale of this symbol.

        This is the minimum of the scales of the jordan blocks.

        OUTPUT:

        an ideal in `O_E`
        """
        if len(self._symbol)==0:
            s = 0
            print(self.rank())
        else:
            s = min([b[0] for b in self._symbol])
        E = self._phi.codomain()
        if self._is_ramified:
            P = E.ideal(self._prime).prime_factors()[0]
            return P**s
        else:
            return E.ideal(self._prime)**s

    def norm(self):
        r"""
        Return the norm of this symbol.

        This is the maximum of the norms of the jordan blocks.

        OUTPUT:

        an ideal in `O_K`
        """
        if self._is_dyadic:
            return min([b[3] for b in self._symbol])
        G = self.gram_matrix() #dirty
        n = G.ncols()
        N = G.diagonal() +[G[i,j]+G[j,i] for i in range(n) for j in range(i)]
        n =  min([x.valuation(self._prime) for x in N])
        return self._prime**n

    def gram_matrix(self):
        r"""
        Return the gram matrix of a representative of this local genus.
        """
        g = [self._gram_from_block(b) for b in self._symbol]
        return matrix.block_diagonal(g).dense_matrix()

    def _gram_from_block(self, block):
        r"""
        Return a gram matrix for this block.

        INPUT:

        - ``block`` -- a tuple of 3 or 4 integers

        OUTPUT:

        - a hermitian matrix over `E` representing this block

        EXAMPLES::

            sage:
        """
        E = self._phi.codomain()
        K = self._phi.domain()
        p = self._phi(self._uniformizer)

        i = block[0]
        m = block[1]
        d = block[2]

        if d == 1:
            u = E(1) # norm class representative of the determinant
        else:
            u = self._u

        # unramified
        if not self._is_ramified:
            return matrix.diagonal(E, m * [p**i])

        pi = self._pi
        H = matrix(E, 2, [0, pi**i, pi.conjugate()**i, 0])
        # ramified non-dyadic
        if not self._is_dyadic:
            if i % 2 == 0:
                return matrix.diagonal(E, (m-1) * [p**(i//2)] + [u*p**(i//2)])
            else:
                assert m % 2 == 0
                return matrix.block_diagonal((m//2) * [H], base_ring=E, subdivide=False)

        # ramified dyadic
        k = block[3]
        e = self._e

        # odd rank
        if m % 2 == 1:
            assert 2*k == i
            r = m // 2
            U = matrix(E, 1, [u*(-1)^r * p**(i//2)])
            return matrix.block_diagonal([U] + r*[H], base_ring=E, subdivide=False)
        # even rank
        r = m // 2 - 1
        if self._is_local_norm((-1)**(m//2)) == (d == 1): # hyperbolic
            assert i+e >= 2*k >= i
            U = matrix(E, 2, [p**k, pi**i, pi.conjugate()**i, 0])
            return matrix.block_diagonal([U] + r*[H], base_ring=E, subdivide=False)
        else:
            assert i+e > 2*k >= i
            u0 = self._u - 1
            # assert K.hilbert_symbol(a,u0) = 1
            U = matrix(E, 2, [p**k, pi**i, pi.conjugate()**i, -p**(i-k)*u0])
            return matrix.block_diagonal([U] + r*[H], base_ring=E, subdivide=False  )

def all_hermitian_genera_by_det(rank, signatures, determinant, phi, max_scale=None):
    r"""
    Return a list of all non-empty global genera with the given conditions.
    """

    from sage.misc.mrange import mrange_iter
    genera = []
    local_symbols = []

    E = phi.codomain()
    K = phi.domain()
    Erel = E.relativize(phi,'a')
    if max_scale == None:
        max_scale = E.ideal(determinant)

    ramified_primes = Erel.relative_discriminant().prime_factors()
    primes = copy(ramified_primes)
    for p in Erel.ideal(determinant).relative_norm().prime_factors():
        if not p in primes:
            primes.append(p)
    ms = max_scale.apply_morphism(Erel.structure()[1]).relative_norm()
    ds = determinant.apply_morphism(Erel.structure()[1]).relative_norm()
    for p in primes:
        det_val = ds.valuation(p)
        mscale_p = ms.valuation(p)
        det_val = det_val / 2
        is_ramified = (p in ramified_primes)
        mscale_p = mscale_p
        if not is_ramified:
            mscale_p /=2
        local_symbol_p = _all_p_adic_genera(p, rank, det_val, mscale_p,
                                            is_ramified, phi)
        local_symbols.append(local_symbol_p)
    # take the cartesian product of the collection of all possible
    # local genus symbols one for each prime
    # and check which combinations produce a global genus
    # TODO:
    # we are overcounting. Find a more
    # clever way to directly match the symbols for different primes.
    for g in mrange_iter(local_symbols):
        # take only non empty genera
        try:
            G = GenusHermitian(g, phi, rank, signatures)
            genera.append(G)
        except ValueError:
            pass
    # for testing
    genera.sort(key=lambda x: [s._symbol for s in x._local_symbols])
    return(genera)

def _all_p_adic_genera(p, rank, det_val, max_scale, is_ramified, phi):
    r"""
    Return all `p`-adic genera with the given conditions.

    This is a helper function for :meth:`all_genera_by_det`.
    No input checks are done.

    Format for a block over a non-dyadic prime is
    ``(scale, rank, det)``
    and over a dyadic prime it is
    ``(scale, rank, det, norm)``.
    Let `p` be a prime ideal in `K` and `P` be the largest ideal of
    `O_E` that contains `p` and is invariant under the involution.
    Then the scale of this genus is `P**scale`.

    INPUT:

    - ``p`` -- a prime ideal in `K`

    - ``rank`` -- the rank of this genus

    - ``det_val`` -- valuation of the determinant at p

    - ``max_scale`` -- an integer the maximal scale of a jordan block

    - ``is_ramified`` -- bool

    - ``phi`` -- the embedding of `K` into `E`

    EXAMPLES::

        sage: E.<a> = CyclotomicField(8)
        sage: K.<b>, phi = E.subfield(a + a**-1)
        sage: p = K.primes_above(3)[0]
        sage: _all_p_adic_genera(p, 2, 4, 2, True, phi)
        [((2, 2, 1),), ((2, 2, -1),)]
        sage: _all_p_adic_genera(p, 2, 4, 4, True, phi)
        [((0, 1, 1), (4, 1, 1)),
        ((0, 1, -1), (4, 1, 1)),
        ((0, 1, 1), (4, 1, -1)),
        ((0, 1, -1), (4, 1, -1)),
        ((2, 2, 1),),
        ((2, 2, -1),)]
    """
    from sage.misc.mrange import cantor_product
    from sage.combinat.integer_lists.invlex import IntegerListsLex
    if is_ramified:
        # the valuation is with respect to p
        # but the scale is with respect to P
        # in the ramified case p = P**2 and thus
        # the determinant of a block is
        #P**(scale*rank) = p**(scale*rank/2)
        det_val *= 2
    scales_rks = [] # contains possibilities for scales and ranks
    for rkseq in IntegerListsLex(rank, length=max_scale+1):   # rank sequences
        # sum(rkseq) = rank
        # len(rkseq) = max_scale + 1
        # now assure that we get the right determinant
        d = 0
        pgensymbol = []
        for i in range(max_scale + 1):
            d += i * rkseq[i]
            # blocks of rank 0 are omitted
            if rkseq[i] != 0:
                pgensymbol.append((ZZ(i), ZZ(rkseq[i])))
        if d == det_val:
            scales_rks.append(pgensymbol)
    # add possible determinant square classes
    if not is_ramified:
        # determinants are all 1
        return [LocalGenusHermitian(p, [b + (ZZ(1),) for b in g], phi)
                                        for g in scales_rks]
    # scale and rank cannot both be odd
    scales_rks = [g for g in scales_rks if all((b[0]*b[1]) % 2 == 0 for b in g)]

    symbols = []
    K = p.ring()
    hyperbolic_det = K.hilbert_symbol(K(-1), K.gen()**2/4 - 1, p)
    if not ZZ(2).divides(p.norm()):   # non-dyadic
        for g in scales_rks:
            n = len(g)
            # possible determinants
            dets = []
            for b in g:
                if b[0] % 2 == 0:
                    dets.append([ZZ(1), ZZ(-1)])
                if b[0] % 2 == 1:
                    dets.append([hyperbolic_det**(b[1]//2)])

            for d in cantor_product(*dets):
                g1 = copy(g)
                for k in range(n):
                    g1[k] += (d[k],)
                symbols.append(LocalGenusHermitian(p, g1, phi))
        return symbols

    e = LocalGenusHermitian(p, (), phi)._e
    for g in scales_rks:
        n = len(g)
        det_norms = []
        for b in g:
            if b[1] % 2 == 1:
                assert b[0] % 2 == 0
                det_norms.append([[ZZ(1), b[0]//2], [ZZ(-1), b[0]//2]])
            if b[1] % 2 == 0:
                dn = []
                i = ZZ(b[0])
                # (i + e) // 2 => k >= i/2
                for k in srange((i/2).ceil(), (i + e)//2):
                    dn.append([ZZ(1), k])
                    dn.append([ZZ(-1), k])
                dn.append([ZZ(hyperbolic_det)**(b[1]//2), (i+e)//2]) #b[1]//2 hyperbolic planes
                if (i+e) % 2 == 1:
                    dn.append([-ZZ(hyperbolic_det)**(b[1]//2), (i+e)//2]) # not hyperbolic
                det_norms.append(dn)

        for dn in cantor_product(*det_norms):
                g1 = copy(g)
                for k in range(n):
                    g1[k] += tuple(dn[k])
                h = LocalGenusHermitian(p, g1, phi)
                if h not in symbols:
                    symbols.append(h)
    return(symbols)

