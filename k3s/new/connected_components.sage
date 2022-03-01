# Helper functions to compute the
# G connected components of some moduli space of K3s
from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap


InducedAut = ("""InducedAut:=function(iso,aut)
local f,g,emb1,emb2,proj1,proj2,aut1,aut2,MyImage;
f:=Range(iso);
g:=Source(iso);
emb1 :=Embedding(g,1);
emb2 :=Embedding(g,2);
proj1:=Projection(g,1);
proj2:=Projection(g,2);
aut1:=proj1*aut[1]*emb1;
aut2:=proj2*aut[2]*emb2;
MyImage:=function(aut,x)
      return Image(aut1,x)*Image(aut2,x);
      end;
aut:= GroupHomomorphismByImagesNC(f,f,
    GeneratorsOfGroup(f),
    List(GeneratorsOfGroup(f), function(i)
          return
            Image(iso,
              MyImage(aut,PreImagesRepresentative(iso,i)
                ) );
        end ));
SetIsInjective( aut, true );
SetIsSurjective( aut, true );
return aut;
end;""")
InducedAut = libgap.function_factory(InducedAut)

def glue_map(L1,L2,L, G1=None,G2=None,alg="sage", notfull=True):
    # check for primitive extension

    D = L/(L1+L2)
    D1 = L1.discriminant_group()
    D2 = L2.discriminant_group()
    self = D1
    other = D2
    Dgap = AbelianGroupGap(D.invariants())
    D1gap = D1.orthogonal_group().domain()
    D2gap = D2.orthogonal_group().domain()

    # the orthogonal_projections with respect to the ambient basis
    A = L.overlattice(L.ambient_module().gens())
    B1 = L1.gens()+A.orthogonal_complement(L1).gens()
    B1 = matrix(B1)
    n = L1.rank()
    k = L1.degree()-L1.rank()
    P1 = matrix.diagonal(n*[1]+k*[0])
    P1 = B1.inverse()*P1*B1

    B2 = L2.gens()+A.orthogonal_complement(L2).gens()
    B2 = matrix(B2)
    n = L2.rank()
    k = L2.degree()-L2.rank()
    P2 = matrix.diagonal(n*[1]+k*[0])
    P2 = B2.inverse()*P2*B2

    imgs1 = [g.lift()*P1 for g in D.gens()]
    imgs1 = [D1(g) for g in imgs1]
    H1g = D1gap.subgroup(imgs1)
    H1 = D1.submodule(imgs1)
    O1 = H1.orthogonal_group([1])
    O1 = H1.orthogonal_group([g.matrix() for g in O1.ambient().gens()])
    H1gap = O1.domain()
    imgs1 = [H1gap(H1(g)) for g in imgs1]
    f1 = Dgap.hom(imgs1,codomain=H1gap)

    imgs2 = [g.lift()*P2 for g in D.gens()]
    imgs2 = [D2(g) for g in imgs2]
    H2 = D2.submodule(imgs2)
    H2g = D2gap.subgroup(imgs2)
    O2 = H2.orthogonal_group([1])
    O2 = H2.orthogonal_group([g.matrix() for g in O2.ambient().gens()])
    H2gap = O2.domain()


    imgs2 = [H2gap(H2(g)) for g in imgs2]
    f2 = Dgap.hom(imgs2,codomain=H2gap)
    f1inv = H1gap.hom([f1.lift(g) for g in f1.codomain().gens()],codomain=f1.domain())

    # now define the glue map on the subgroups.
    phi = H1gap.hom([f2(f1inv(g)) for g in H1gap.gens()],codomain=f2.codomain())
    if G1 is None:
        G1 = D1.O()
    if G2 is None:
        G2 = D2.O()

    def is_contained_in_stab(G,T, D):
        gens = [T(g.lift()) for g in D.gens()]
        return all(all(g*s in D for g in gens) for s in G.gens())

    if alg=="sage":
        S1 = stabilizer_sage(G1, H1g)
        S2 = stabilizer_sage(G2, H2g)
    elif alg=='magma':
        S1 = stabilizer_magma(G1, D1, H1)
        S2 = stabilizer_magma(G2, D2, H2)


    assert is_contained_in_stab(S1,D1,H1)
    assert is_contained_in_stab(S2,D2,H2)


    D = L1.overlattice((L1+L2).gens()).discriminant_group()
    if notfull:
      OD = D.orthogonal_group()
    else:
      OD = D.orthogonal_group([])
    Dg = OD.domain()
    i1 = D1.hom([D(g) for g in D1.gens()])
    i2 = D2.hom([D(g) for g in D2.gens()])

    G12 = G1.gap().DirectProduct(G2.gap())
    embG1 = G12.Embedding(1)
    embG2 = G12.Embedding(2)
    D1 = G1.one().gap().Source()
    D2 = G2.one().gap().Source()
    D12 = D1.DirectProduct(D2)
    embD1 = D12.Embedding(1)
    embD2 = D12.Embedding(2)
    projD1 = D12.Projection(1)
    projD2 = D12.Projection(2)

    act1 = S1.hom([O1(g) for g in S1.gens()],codomain=O1)
    im1 = act1(S1)
    ker1 = [embG1.Image(k.gap()) for k in act1.kernel().gens()]

    act2 = S2.hom([O2(g) for g in S2.gens()],codomain=O2)
    im2 = act2(S2)
    ker2 = [embG2.Image(k.gap()) for k in act2.kernel().gens()]

    imgs = [Dg(
              i1(self.linear_combination_of_smith_form_gens((G1.domain()(projD1.Image(d))).exponents()))
              +i2(other.linear_combination_of_smith_form_gens((G2.domain()(projD2.Image(d))).exponents()))
              ).gap() for d in D12.GeneratorsOfGroup()]
    D12toDg = D12.GroupHomomorphismByImages(Dg.gap(),imgs)


    ###############################
    # compute the stabiliser of L in S1 x S2
    # this is done by hand to avoid an
    # expensive stabilizer computation
    phig = phi.gap()
    im2_phi = O1.subgroup([phig*g1.gap()*phig.Inverse()
                            for g1 in im2.gens()])
    im = im1.subgroup(im1.intersection(im2_phi).gap().SmallGeneratingSet())
    stab = [(act1.lift(x),
              act2.lift(im2(phig.Inverse()*x.gap()*phig)))
            for x in im.gens()]
    stab = [embG1.Image(x[0].gap())*embG2.Image(x[1].gap())
            for x in stab]
    stab += ker1
    stab += ker2   # tends to generate a huge group
    # convert to sage
    stab1 = G1.subgroup([s[0] for s in stab])
    stab2 = G2.subgroup([s[1] for s in stab])
    stab0 = [InducedAut(D12toDg,s) for s in stab]
    stab = G12.Subgroup(stab)

    #stab0 = OD.ambient().subgroup(stab0)
    #stab0 = D.orthogonal_group(stab0.gens())
    if notfull:
      stab0 = OD.subgroup(stab0)
      stab0 = L.discriminant_group().O().subgroup(stab0.gens())
      stabhom = stab.GroupHomomorphismByImagesNC(stab0.gap(),[g.gap() for g in stab0.gens()])
      return H1, H2, phi, stab,[stabhom,G12.Projection(1)],stab0, stab1,stab2
    else:
      return H1, H2, phi, stab,[1,G12.Projection(1)],1, stab1,stab2


def stabilizer_sage(G,D):
    from sage.libs.gap.libgap import libgap
    from sage.env import SAGE_EXTCODE
    gapcode = SAGE_EXTCODE + '/gap/subgroup_orbits.g'
    libgap.Read(gapcode)
    OnSubgroups = libgap.function_factory("OnSubgroups")

    stab = G.gap().Stabilizer(D.gap(),OnSubgroups)
    return G._subgroup_constructor(stab)

def stabilizer_magma(G,T, D):
    r"""
    Compute the stabilizer in G of D<=T.
    """
    from sage.rings.all import GF,ZZ
    if D.cardinality()==1:
        return G
    p = D.invariants()[0]
    field =  GF(p)
    # we continue with magma
    TT = T.submodule((T.cover() & (1/p)*T.relations()).gens())
    G0 = TT.orthogonal_group()
    action_hom = G.hom(G.gens(), codomain=G0,check=False)
    gens = [i.matrix().change_ring(field)
            for i in action_hom.image(G).gens()]
    from sage.interfaces.magma import magma as m
    degree = len(TT.invariants())
    GL = m.GeneralLinearGroup(degree, field)
    G = GL.sub(gens)
    gensD = [TT(g).vector().change_ring(field) for g in D.smith_form_gens()]
    V = m.VectorSpace(field,degree)
    Dmagma = V.sub(gensD)

    stab = G.Stabiliser(Dmagma)
    stab = stab.Generators().SetToSequence()
    stab = [g.Matrix().sage().change_ring(ZZ) for g in stab]
    S = G0.subgroup(stab)
    stab = action_hom.preimage(S)
    return stab



def index_K_Kplus(self):
    if not self.signature_pair()[0]==2:
        raise ValueError("")
    if self.rank()==2:
        OL = self.orthogonal_group()
        Oq = self.discriminant_group().orthogonal_group()
        D = OL.hom([Oq(g) for g in OL.gens()],check=False)
        K = D.kernel()
        A = AbelianGroupGap([2])
        from sage.groups.fqf_orthogonal_spin import spin_exact
        G = self.inner_product_matrix()
        imgs = [g.matrix().det()*spin_exact(G,g.matrix()) for g in K.gens()]
        ds = K.hom([A.gen(0)^((1-g.sign())/2) for g in imgs],codomain=A)
        return K.order()/ds.kernel().order()

    from sage.groups.fqf_orthogonal_spin import GammaA
    T = self.discriminant_group()
    rank = self.rank()
    det = self.determinant()
    S = (2*det).prime_divisors()
    Gamma = GammaA(S)
    sigma_sharp = Gamma.sigma_sharp(rank, det, T)
    gammaS = Gamma.subgroup(Gamma.gammaS(False))
    gammaSplus = Gamma.subgroup(Gamma.gammaS(True))
    c = sigma_sharp.intersection(gammaS).order()
    cplus = sigma_sharp.intersection(gammaSplus).order()
    return c/cplus


def is_real(k3):
    T = k3.transcendental_lattice()
    iK = index_K_Kplus(IntegralLattice(T.gram_matrix()))
    if iK==2:
        return True

    L1 = k3.symplectic_co_invariant_lattice()
    L2 = k3.invariant_lattice()
    L = k3.symplectic_invariant_lattice()
    _,_,_,stab,stabhom, S0, S1,S2= glue_map(T,L2,L,T.image_in_Oq(),L2.image_in_Oq(),alg="magma")

    H = k3.L()
    _,_,_,_,_, _, S1b,_ = glue_map(L,L1,H,S0,L1.image_in_Oq(),notfull=False)
    S = stabhom[0].PreImage(S1b)
    S = stabhom[1].Image(S)
    S = S1._subgroup_constructor(S)
    if T.rank()==2:
        o1 = S.intersection(T.image_in_Oq(False)).order()
        o2 = S.intersection(T.image_in_Oq(True)).order()
        index_S_Splus = o1/o2
        assert index_S_Splus <= 2
        return index_S_Splus == 2

    f = sigma_plus(L)
    C = f.image(S)
    assert C.order()<=2
    return C.order()==2

def sigma_plus(L):
    r"""
    O(D_L) --> Sigma(T) / \sigma^#(L)(Gamma_Q^+ \cap \Sigma(L))
    """
    from sage.groups.fqf_orthogonal_spin import det_spin_homomorphism
    f = det_spin_homomorphism(L)
    SigmaModSigmaSharp = f.codomain()
    gamma = SigmaModSigmaSharp.cover()
    sigma_sharp = SigmaModSigmaSharp.relations()
    gammaSplus = gamma.subgroup(gamma.gammaS(True))
    sigma = gamma.subgroup(tuple(SigmaModSigmaSharp.lift(f(g)) for g in f.domain())+sigma_sharp.gens())
    newrels = sigma.intersection(gammaSplus)
    newrels = SigmaModSigmaSharp.subgroup([SigmaModSigmaSharp(g) for g in newrels.gens()])
    C = cod.quotient(newrels)
    g = SigmaModSigmaSharp.hom([C(g) for g in SigmaModSigmaSharp.gens()],codomain=C)
    return g*f

