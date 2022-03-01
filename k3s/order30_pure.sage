def glue(M1,N1,glue_order,p):
    M = M1.L
    fM = M1.iso
    GM = M1.Oq_equiv()
    N = N1.L
    fN = N1.iso
    GN = N1.Oq_equiv()
    DM = M.discriminant_group(p)
    DN = N.discriminant_group(p)
    ext = extensions(M,fM,N,fN,GM,GN,glue_order,p,DM,DN)
    return [LatticeWithIsometry(e[0],e[1].change_ring(ZZ),Oq_equiv=e[2]) for e in ext]


def split_sig_n(sig_pair, n, q):
    r"""
    """
    assert q.is_prime() and not q.divides(n)
    pos, neg = sig_pair
    rk = pos + neg
    assert pos in [0,1,2,3]
    weights = [euler_phi(n),euler_phi(n*q)]
    m = 2
    ceiling = [rk // d for d in weights]
    floor = [0, 0]
    if pos >= 2:
        # the order must be nq
        floor[1] = 1
    if pos % 2 == 1:
        # there must be an invariant vector of pos. square
        if e != 0:
            return
        floor[0] = 1
    R = IntegerListsLex(max_sum=rk,length=m,floor=floor, ceiling=ceiling)
    R = [[r[k]*weights[k] for k in range(m)] for r in R]
    R = [r for r in R if sum(r)==rk]
    for ranks in R:
        if e == 0:
            assert pos >=1
            sigp = sig_k3(n, ranks[0], 1)
        else:
            assert pos in [0,2]
            sigp = sig_k3(n, ranks[0], 0)
        pospq = (pos//2)*2
        sigpq = sig_k3(n * q, ranks[1], pospq)
        signatures = [sigp,sigpq]
        yield [ranks, signatures]

def splitnq(genus, n, q, k3_unobstructed=True,verbose=5,splitsig=None):
    r"""
    Classify lattices with isometries up to isomorphism
    such that the isometry has minimal polynomial `\Phi_{nq}\Phi_{n}`
    and the lattice is in ``genus``.

    INPUT:

    - ``n`` -- a number
    - ``q`` -- a prime number not dividing n
    - ``k3_unobstructed`` -- (default: ``True``)

    OUTPUT:

    A list consisting of triples ``[A, a, Oa]`` where
    ``A`` is a lattice in ``genus``,
    ``a`` is an isomety in ``a`` and
    ``Oa`` is the image in the discriminant group
    of the stabilizer of ``a``.
    """
    assert not q.divides(n) and q.is_prime()
    if splitsig is None:
        splitsig = split_sig_n(genus.signature_pair(),n,q)
    for rks,sigs in splitsig:
        degsE = euler_phi(n), euler_phi(q*n)
        rksE = [rks[0]//degsE[0],rks[1]//degsE[1]]

        if rks[0] == 0:
            # there is nothing to glue
            for Mh in all_lattice_with_isometry(q*n,rksE[1],genus.det(),
                                                   genus.level(),signatures=sigs[1],
                                                   min_scale=genus.scale()):
                M = Mh.L
                fM = Mh.iso
                if M.genus() == genus and not (k3_unobstructed and is_obstructed(M,fM)):
                    GM =Mh.Oq_equiv()
                    yield M, fM, GM
            continue

        if rks[1] == 0:
            # there is nothing to glue
            for Rh in all_lattice_with_isometry(n,rksE[0],genus.det(),
                                                   genus.level(),signatures=sigs[0],
                                                   min_scale=genus.scale()):
                R = Rh.L
                fR = Rh.iso
                if R.genus() == genus and not (k3_unobstructed and is_obstructed(R,fR)):
                    GR = Rh.Oq_equiv()
                    yield R, fR, GR
            continue

        min_scale = genus.scale()
        min_norm = genus.norm()
        max_level = q * genus.level()

        m = degsE[0]*min(rksE)
        max_det = q^m*genus.det()
        c = cyclotomic_polynomial
        for glue, M_det, R_det in glue_det(genus.det(),
                                           c(n)^rksE[0],
                                           c(n*q)^rksE[1]):
            alwi = all_lattice_with_isometry
            for MG in alwi(q*n, rksE[1], M_det, max_level,
                           signatures=sigs[1],
                           min_scale=min_scale,
                           min_norm=min_norm,
                           return_genera=True):
                try:
                    Mgenus = MG.trace()
                except AttributeError:
                    Mgenus = MG
                for RG in alwi(n, rksE[0], R_det, max_level,
                               signatures=sigs[0], min_scale=min_scale,
                               return_genera=True):
                    try:
                        Rgenus = RG.trace()
                    except AttributeError:
                        Rgenus = RG
                    if not is_admissible(Mgenus, Rgenus, genus, q):
                        continue
                    if verbose>2:
                        print('computing representatives of %s'%MG)
                    tmp = MG.representatives()
                    print("done computing representatives")
                    for ME in tmp:
                        M, fM = trace_lattice(ME, order=q*n)
                        M = IntegralLattice(M)
                        Mh = LatticeWithIsometry(M,fM,order=q*n,gramE=ME,magmaRep=MG.representative())
                        if k3_unobstructed and is_obstructed(M, fM):
                            print("obstructed")
                            continue
                        GM = Mh.Oq_equiv()
                        cm = sage.structure.element.get_coercion_model()
                        cm.reset_cache()
                        if verbose>2:
                            print('computing representatives of %s'%RG)
                        for RE in RG.representatives():
                            R, fR = trace_lattice(RE, order=n)
                            R = IntegralLattice(R)
                            Rh = LatticeWithIsometry(R,fR,order=n,gramE=RE,magmaRep=RG.representative())
                            if k3_unobstructed and is_obstructed(R, fR):
                                continue
                            GR = Rh.Oq_equiv()
                            if verbose >0:
                                print('glueing')
                                print(M.genus())
                                print('and')
                                print(R.genus())
                                print('above %s'%q)
                            print("computing extensions of order %s in"%glue)
                            ext = extensions(M, fM, R, fR, GM,
                                         GR, glue, q,
                                         target_genus=genus)
                            for ex in ext:
                                assert ex[0].genus() == genus
                                yield ex



fi = open("results/purely_ns/order15.txt","r")
k3s = [K3SurfaceAutGrp_from_str(k3) for k3 in fi if k3[:8]!="complete"]
fi.close()

fi = open("results/order15.txt","r")
k3sym = [K3SurfaceAutGrp_from_str(k3) for k3 in fi if k3[:8]!="complete"]
fi.close()

def order30_pure(k3):
  result = []
  g = k3.distinguished_generator()
  C3 = k3.L().kernel_sublattice(g^2+g+1)
  C5 = k3.L().kernel_sublattice(g^4+g^3+g^2+g+1)
  assert C3.rank()*C5.rank()==0
  if C3.rank() >0:
    p=3
  else:
    p=5
  d = k3.transcendental_lattice().det()
  q = d.prime_factors()[0]
  pty = ptype(k3.NS(),g,p)
  print("computing actions on NS")
  onNS = [LatticeWithIsometry(a,b,Oq_equiv=c) for a,b,c in pnq_actions(2,pty)]
  print("done")
  for C30 in splitnq(k3.transcendental_lattice().genus(),15,2,verbose=5):
    C30 = LatticeWithIsometry(C30[0],C30[1],Oq_equiv=C30[2])
    for C2p in onNS:
      print("gluing: %s and %s over %s"%(C30.L.rank(),C2p.L.rank(),d))
      for g in glue(C30,C2p,d,q):
        print("match")
        result.append(K3SurfaceAutGrp(g.L,[],g.iso,30))
  return result

def order30_impure(k3):
  result = []
  g = k3.distinguished_generator()
  N = k3.symplectic_invariant_lattice()
  C3 = N.kernel_sublattice(g^2+g+1)
  C5 = N.kernel_sublattice(g^4+g^3+g^2+g+1)
  assert C3.rank()*C5.rank()==0
  if C3.rank() >0:
    p=3
  else:
    p=5
  d = k3.transcendental_lattice().det()
  pf = [e for e in d.prime_factors() if e!=2]
  assert len(pf)==1
  q = pf[0]
  pty = ptype(N.sublattice((k3.NS() & N).gens()),g,p)
  print("computing actions on NS")
  onNS = [LatticeWithIsometry(a,b,Oq_equiv=c) for a,b,c in pnq_actions(2,pty)]
  print("done")
  for C30 in splitnq(k3.transcendental_lattice().genus(),15,2,verbose=5):
    C30 = LatticeWithIsometry(C30[0],C30[1],Oq_equiv=C30[2])
    for C2p in onNS:
      print("gluing: %s and %s over %s"%(C30.L.rank(),C2p.L.rank(),d))
      for g in glue(C30,C2p,q^d.valuation(q),q):
        print("match")
        result.append(K3SurfaceAut(g.L,cofixed[0].twist(-1),g.iso,g.Oq_equiv()))
  return result
