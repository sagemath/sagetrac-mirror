from sage.interfaces.gap import get_gap_memory_pool_size
memory_gap = get_gap_memory_pool_size()
set_gap_memory_pool_size(4*memory_gap)

libgap.eval("SetRecursionTrapInterval(10000000)");

from sage.quadratic_forms.genera.genus import all_genera_by_det
load("databases/hashimotoslattices.sage")
load("k3_classification/symplectic.sage")
attach("k3_classification/hermitian.sage")

transcendental_values = [n for n in range(2,100) if euler_phi(n) <=20]

def C1part(T, n, c1, Hinv):
    r"""

    INPUT:

    - ``T`` -- the transcendental lattice
    - ``n`` -- an integer - the order of the automorphism
    - ``c1`` -- the rank of the invariant lattice
    - ``Hinv``-- the target lattice

    OUTPUT:

    - a list of possible genera of the invariant part

    """
    assert c1 >= 1
    d = Hinv.determinant()
    t = T.determinant()
    # bounds on |det C1|
    if T.rank() + c1 == Hinv.rank():
        min_glue = (d/t).denominator()
        min_det = ZZ(d/t * min_glue^2)
    max_scale = max(Hinv.discriminant_group().invariants())
    if len(n.prime_factors())==1:
        p = n.prime_factors()[0]
        assert c1 == Hinv.rank()-T.rank()
        m = min(T.rank()//euler_phi(n), c1)
        max_glue = gcd(p^m, t)
        max_det = ZZ(d/t * max_glue^2)
        max_scale *=p
    else:
        # no glue between c1 and cn
        # probably one can improve here by taking a look at the possible
        # polynomials
        # cheap bound
        r = Hinv.rank() - (T.rank() + c1)   # rank of the orthogonal complement
        max_det = d
        min_scale = H.genus().scale()
        min_det = min_scale**c1
        for p in ZZ(n).prime_divisors():
            # bound for numer of cp factors in cofixed part is r//(p-1)
            max_det *= p^(r//(p-1))
            max_scale *= p
    C1possibilities = []
    for det in ZZ(max_det).divisors():
        if min_det.divides(det):
            poss = all_genera_by_det([1, c1-1], det, max_scale=max_scale)
            C1possibilities += poss
    return C1possibilities

def inject_orthogonal_group(iA, iB, C, GA, GB):
    r"""
    Extend the action of ``GA`` on ``A`` trivially to ``B`` inside
    ``C = A + B``.

    INPUT:

    - ``iA``, ``iB`` -- injections of `A`, `B` into `C = A + B`
    - ``GA`` -- something that acts on ``A``

    OUTPUT:

    A subgroup of GL(C).

    EXAMPLES::

        sage:
    """
    assert C.V() == iA.image().V() + iB.image().V()
    Oq = C.orthogonal_group()
    CC = Oq.domain()
    n = len(iA.image().gens())
    ABgen = [(CC(C(g))).gap() for g in iA.image().gens() + iB.image().gens()]
    gensGA_embedded = []
    for f in GA.gens():
        # transport the action to C
        imgs = [iA(a*f) for a in iA.domain().gens()]
        # convert it to gap
        imgs = [(CC(c)).gap() for c in imgs]
        imgs += ABgen[n:]
        f = libgap.GroupHomomorphismByImages(CC.gap(), ABgen, imgs)
        gensGA_embedded.append(f)
    gensGB_embedded = []
    for f in GB.gens():
        # transport the action to C
        imgs = ABgen[:n]
        imgs += [iB(b*f) for b in iB.domain().gens()]
        # convert it to gap
        imgs = [(CC(c)).gap() for c in imgs]
        f = libgap.GroupHomomorphismByImages(CC.gap(), ABgen, imgs)
        gensGB_embedded.append(f)
    return Oq.subgroup(gensGA_embedded+gensGB_embedded)

def extensions_old(T, fT, p, C1, det_target=None):
    r"""
    Return all possible extensions of `(T, fT)` by (C1,id)`
    along the prime `p`.

    INPUT:

    - ``T`` -- transcendental lattice
    - ``fT`` -- an isometry of order ``p^n`` of ``T``
    - ``p`` -- a prime over which to glue
    - ``C1`` -- the invariant lattice
    - ``det_target`` -- an integer

    OUTPUT:

    a list of lattices
    """
    if det_target != None:
        glue_order = ZZ(abs(C1.determinant()*T.determinant()/det_target))
        if glue_order.is_square():
            glue_order = glue_order.sqrt()
        else:
            return []
    DC1 = C1.discriminant_group(p)
    DC1inv = ((1/p)*C1 & C1.dual_lattice()).gens()
    DC1inv = DC1.submodule(DC1inv)
    # compute fixed points of the discriminant_group of T
    Dinv = T.span((fT - fT**0).inverse()) & T.dual_lattice()
    DT = T.discriminant_group()
    DTinv = DT.submodule(Dinv.gens())
    # take the direct sum
    D, iT, iC = DT.direct_sum(DC1)
    if DTinv.cardinality()==1 or DC1inv.cardinality()==1:
        return [D.W()]
    # invariant parts for gluing
    DTinv = D.submodule([iT(g) for g in Dinv.gens()])
    DC1inv = D.submodule([iC(g) for g in DC1inv.gens()])
    # embedd the orthogonal groups
    GC1 = C1.image_in_Oq()
    GC1 = inject_orthogonal_group(iC, iT, D, GC1)
    if p == 2:
        GT = inject_orthogonal_group(iT, iC, D, T.image_in_Oq())
        gens = [g.gap() for g in GC1.gens() + GT.gens()]
        G = D.orthogonal_group().subgroup(gens)
    else:
        # TODO: the stabilizer of the action of T is still missing!
        G = GC1
    # compute the primitive extensions modulo G
    extensions = []
    for glue in D.all_primitive_modulo(DTinv, DC1inv, G, combined=False):
        if det_target != None and glue_order != glue.cardinality():
            continue
        ext = D.W().overlattice([g.lift() for g in glue.gens()])
        extensions.append(ext)
    return extensions


def extensions_old1(T, fT, p, C1, det_target=None):
    r"""
    Return all possible extensions of `(T, fT)` by (C1,id)`
    along the prime `p`.

    INPUT:

    - ``T`` -- transcendental lattice
    - ``fT`` -- an isometry of order ``p^n`` of ``T``
    - ``p`` -- a prime over which to glue
    - ``C1`` -- the invariant lattice
    - ``det_target`` -- an integer

    OUTPUT:

    a list of lattices
    """
    if det_target != None:
        glue_order = ZZ(abs(C1.determinant()*T.determinant()/det_target))
        if glue_order.is_square():
            glue_order = glue_order.sqrt()
        else:
            return []
    DC1 = C1.discriminant_group()
    DC1inv = ((1/p)*C1 & C1.dual_lattice()).gens()
    DC1inv = DC1.submodule(DC1inv)
    # compute fixed points of the discriminant_group of T
    Dinv = T.span((fT - fT**0).inverse()) & T.dual_lattice()
    DT = T.discriminant_group()
    DTinv = DT.submodule(Dinv.gens())

    GC1 = C1.image_in_Oq()
    if p == 2:
        GT = T.image_in_Oq()
    else:
        # TODO: the stabilizer of the action of T is still missing!
        GT = T.orthogonal_group([-fT^0]).gens()
        GDT = DT.orthogonal_group()
        GT = GDT.subgroup(GT)
    # compute the primitive extensions modulo G
    extensions = []
    for glue in DT.all_primitive_prime(DC1,DTinv,DC1inv,GT,GC1):
        if det_target != None and glue_order != glue.cardinality():
            continue
        ext = glue.W().overlattice([g.lift() for g in glue.gens()])
        extensions.append(ext)
    return extensions


#####################################
# prime order case n = p
#####################################

def all_actions_p(Cp, fp, p, H):
    r"""
    Return all possible glues of Cp + C1 that are compatible with the action.

    INPUT:

    - ``Cp`` -- an integral lattice

    - ``fp`` -- action on ``Cp`` of order ``p`` given by a matrix

    - ``p`` -- prime number
    """
    from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring, p_adic_symbol
    assert fp.minimal_polynomial().is_irreducible()
    cand = C1part(Cp, p, H.rank()- Cp.rank(), H)
    #print("testing %s candidates for C1:"%len(cand))
    #print(Cp.gram_matrix())
    actions = []
    genH = H.genus()
    primes = [sym.prime() for sym in genH.local_symbols()]
    genCp = Cp.genus()
    for C1 in cand:
        flag = False
        if not ZZ(C1.determinant()*Cp.determinant()/H.determinant()).is_square():
            continue
        G = C1.direct_sum(Cp.genus())
        # check that G is locally equivalent at all primes except p.
        primesG = set([sym.prime() for sym in G.local_symbols()]+primes)
        if any([genH.local_symbols(p=q)!=G.local_symbols(p=q) for q in primesG if q!=p]):
            continue
        # compute all possible glueings
        for C1l in C1.representatives():
            C1l = IntegralLattice(C1l)
            exts = extensions(Cp, fp, p, C1l, H.determinant())
            for ext in exts:
                if ext.genus() == H.genus():
                    actions.append(ext)
            if len(exts)==0:
                break
            #print("found %s actions"  %len(actions))
    return actions

def all_actions_H(Hinv):
    r"""
    """
    h = Hinv.rank()
    actions = []
    for n in range(6000):
        if euler_phi(n) >= h:
            continue
        n = ZZ(n)
        if n == 2:
            actions += order2_actions(H)
        elif n.is_prime():
            p = n
            actions += orderp_actions(H, p)
        elif len(n.prime_factors())==1:
            continue
        else:
            continue





######## the case of order 2 #########

def order2_actions_old(H):
    r"""
    """
    H = IntegralLattice(H)
    acts = []
    max_scale = H.discriminant_group().invariants()[-1] * 2
    for rkC2 in range(2, H.rank()): # maximal rank is rk H - 1 ... there is still an invariant vector!
        print("rank of anti-invariant: %s"%rkC2)
        actsrk = []
        f = -matrix.identity(rkC2)
        m = min(rkC2, H.rank()-rkC2)
        max_det = ZZ(H.determinant() * 2^m)
        for det in max_det.divisors():
            C2genera = all_genera_by_det([2, rkC2-2], det, max_scale=max_scale)
            #print("Testing %s possible genera for C2" %len(C2genera))
            for C2 in C2genera:
                #print(C2.representative())
                for C2l in C2.representatives():
                    actsrk += all_actions_p(IntegralLattice(C2l), f, 2, H)
        fr =  matrix.block_diagonal([f, matrix.identity(m)])
        acts += [[L, fr] for L in actsrk]
    return acts


############## odd prime order actions ####################

def orderp_actions(H, p, m=1):
    r"""
    INPUT:

    - ``H`` -- an integral lattice; the symplectic invariant lattice

    - ``p`` a prime number
    """
    actions = []
    p = ZZ(p)
    if not p.is_prime():
        raise ValueError("")
    h = H.rank()
    invsH = H.discriminant_group().invariants()
    max_scale = ZZ(p)*H.discriminant_group().invariants()[-1]
    for r in range(m, h//(p-1)+ 1):
        t = r *(p-1)
        if t<2:
            continue
        c1 = h - t
        if c1 < 1:
            continue
        print("Rank of the transcendental lattice is %s"%t)
        max_det = ZZ(H.determinant() * p^(min(r,c1)))
        # there should also be a minimal determinant
        # and a minimal scale
        min_scale = 1
        min_det = 1
        if len(invsH) == h:
            min_scale = invsH[0]
        if len(invsH) > c1:
           min_det = ZZ.prod(invsH[:-c1])

        for det in max_det.divisors():
            if not min_det.divides(det):
                continue
            Tlist = all_lattice_with_isometry(p, r, det, max_scale, min_scale=min_scale)
            for Tf in Tlist:
                T, fp = Tf
                print(T, T.rank(), T.det())
                T = IntegralLattice(T)
                f = matrix.block_diagonal([fp,matrix.identity(c1)])
                acts = all_actions_p(T, fp, p, H)
                acts = [(L,f) for L in acts]
                actions += acts
    return actions



############## composite order ##########

def all_actions_n(Cn, fn, n, H, c1):
    r"""
    """
    assert fn.minimal_polynomial() == cyclotomic_polynomial(n)
    cand = C1part(Cn, n, c1, H)
    #print("testing %s candidates for C1:"%len(cand))
    #print(Cp.gram_matrix())
    actions = []
    for C1 in cand:
        fn1 = matrix.block_diagonal([fn,matrix.identity(C1.rank())])
        if C1.rank() + Cn.rank() == H.rank() and C1.direct_sum(Cn.genus()) != H.genus():
            continue
        for C1l in C1.representatives():
            C1l = IntegralLattice(C1l)
            print(Cn.gram_matrix(), C1l.gram_matrix())
            CnC1 = Cn.direct_sum(C1l)
            # now we need all primitive embeddings
            actions += primitive_embeddings(CnC1,H)
    return actions



def ordern_actions(H, n):
    r"""
    """
    actions = []
    n = ZZ(n)
    if n.is_prime():
        raise ValueError("")
    h = H.rank()
    invsH = H.discriminant_group().invariants()
    max_scale = ZZ(n.radical())*H.discriminant_group().invariants()[-1]
    for rE in range(1, h//euler_phi(n) + 1):
        t = rE * euler_phi(n)  # rank of the transcendental lattice
        for c1 in range(1, h - t + 1):
            r = h - t - c1
            max_det = ZZ(H.determinant() * n.radical()^(min(rE, r)))
            # there should also be a minimal determinant
            # and a minimal scale
            min_scale = H.genus().scale()
            min_det = 1
            if len(invsH) > c1 + r:
                min_det = ZZ.prod(invsH[:-(c1 + r)])
            for det in max_det.divisors():
                if not min_det.divides(det):
                    continue
                Tlist = all_lattice_with_isometry(n, r, det, max_scale, min_scale=min_scale)
                for Tf in Tlist:
                    T, f = Tf
                    T = IntegralLattice(T)
                    tmp = all_actions_n(T, f, n, H, c1)
                    print(tmp)
                    actions += tmp
    return actions







#############################################################
# the K3 side of things
#############################################################


def transport(A, B, a):
    r"""
    Transport a to O(B).

    INPUT:

    - ``A``, ``B`` -- isomorphic torsion quadratic_forms in normal form
    - ``a`` -- an element of the orthogonal group of `A`

    OUTPUT:

    - ``b`` -- an element of the orthogonal group of `B`
    """
    GA = A.gram_matrix_quadratic()
    GB = B.gram_matrix_quadratic()
    assert GA == GB or ((GA+GB).denominator()==1 and 2.divides(gcd((GA+GB).diagonal())))
    a = a.matrix()   # wrt smith form gens
    a = A._to_smith() * a * A._to_gens()      # wrt normal form gens
    a = B._to_gens() * a * B._to_smith() # wrt smith gens
    return B.orthogonal_group()(a)

def MaximalK3surfaceAut(fixed, cofixed, g):
    r"""
    """
    n = fixed.orthogonal_group([g]).order()

    # Glue the K3 lattice from the fixed and cofixed lattice
    H, (ifix, icofix) = fixed.direct_sum(cofixed, return_embeddings=True)
    fixed = H.sublattice(ifix.image().gens())
    cofixed = H.sublattice(icofix.image().gens())
    q_fix = fixed.discriminant_group().normal_form()
    gens = cofixed.twist(-1).discriminant_group().normal_form().gens()
    q_cofix = cofixed.discriminant_group().submodule_with_gens(gens)
    glue = [q_fix.gens()[k].lift() + q_cofix.gens()[k].lift() for k in range(len(q_fix.gens()))]
    H = H.overlattice(glue)
    assert H.determinant() == -1
    assert H.is_even()
    assert H.signature_pair() == (3,19), "%s,%s"%H.signature_pair()

    # create the group G0 and g
    g = matrix.block_diagonal([g, matrix.identity(cofixed.rank())])
    g = fixed.orthogonal_group([g]).gen(0)
    g_q_fix = q_fix.orthogonal_group()(g)
    g_q_cofix = transport(q_fix, q_cofix, g_q_fix)

    O = cofixed.orthogonal_group()
    Oq = q_cofix.orthogonal_group()

    phi = O.hom([Oq(h) for h in O.gens()], check=False)
    G0 = phi.kernel()

    if not g_q_cofix in phi.image(phi.domain()):
        raise ValueError("Action on the discriminant group is not in the image of the orthogonal group")

    g_cofix = phi.lift(g_q_cofix).matrix()
    gn = g.matrix()*g_cofix

    return K3SurfaceAutGrp(H, G0, gn, n)




class K3SurfaceAutGrp(object):
    r"""
    """
    def __init__(self, H, G0, g, n):
        r"""
        """
        self.H = H
        self.G0 = G0
        self.g = g
        self.G = H.orthogonal_group([h.matrix() for h in G0.gens()] + [g])
        self.n = n

    @cached_method
    def transcendental_lattice(self):
        r"""
        """
        R.<x> = QQ[]
        cn = cyclotomic_polynomial(self.n,x)
        g = self.g
        S = self.symplectic_invariant_lattice()
        K = self.H.span(cn(g).integer_kernel().gens())
        return self.H.sublattice((K & S).gens())

    @cached_method
    def picard_lattice(self):
        r"""
        """
        return self.H.orthogonal_complement(self.transcendental_lattice())

    @cached_method
    def symplectic_invariant_lattice(self):
        r"""
        """
        V = self.H
        for g in self.G0.gens():
            g = g.matrix()
            V = V & V.span((g - 1).kernel().gens())
        return self.H.sublattice(V.gens())

    @cached_method
    def invariant_lattice(self):
        r"""
        """
        HG = self.symplectic_invariant_lattice()
        h = HG & HG.span((self.g.matrix() - 1).kernel().gens())
        return self.H.sublattice(h)

    @cached_method
    def symplectic_co_invariant_lattice(self):
        r"""
        """
        return self.H.orthogonal_complement(self.symplectic_invariant_lattice())

    def str(self):
        r"""
        """
        s =  "[" + str(self.n) + ", ["
        s += str(self.H.basis_matrix().list()) + ","
        s += str(self.H.inner_product_matrix().list()) + "],"
        s += str([g.matrix().list() for g in self.G0.gens()]) + ","
        s += str(self.g.list()) + "]"
        return s

def aut_from_str(s):
    r"""

    """
    L = sage_eval(s)
    n = L[0]
    basis = matrix(QQ,22,22,L[1][0])
    gram = matrix(QQ,22,22,L[1][1])
    H = IntegralLattice(gram, basis)

    G0_gens = [matrix(QQ,22,22,h) for h in L[2]]
    G0 = H.orthogonal_group(G0_gens)
    g = matrix(QQ,22,22,L[3])

    return K3SurfaceAutGrp(H, G0, g, n)



def FullGlue(A,B):
    r"""
    """
    AB, (iA, iB) = A.direct_sum(B, return_embeddings=True)
    A = AB.sublattice(iA.image().matrix())
    B = AB.sublattice(iB.image().matrix())

    DA = A.discriminant_group().normal_form()
    tmp = B.discriminant_group().twist(-1).normal_form().gens()
    DB = B.discriminant_group().submodule_with_gens(tmp)

    C = AB.overlattice([DA.gen(k).lift() + DB.gen(k).lift() for k in range(len(DA.gens()))])
    assert C.determinant().abs() == 1
    return C

#######################################

def classify_ord_p(L, p, file_name):
    result = open(file_name,"w")
    classifi = []
    not_realized = []
    for k in L:
        print("Invariant Lattice number %s"%k)
        fix = fixed[k]
        print(fix.gram_matrix())
        print(" ")
        cofix = cofixed[k].twist(-1)
        for a in k3_prime_power(fix.genus(), p, 1): #orderp_actions(fix, p)
            fix_a = a[0]
            g = a[1]
            try:
                aut = MaximalK3surfaceAut(fix_a,cofix, g)
                classifi.append(aut)
                s = aut.str()
                result.write(s+ "\n")
            except ValueError:
                not_realized.append([fix,a])
    result.close()
    return classifi, not_realized

###########################################################
def classify_rank4(L, file_name):
    result = open(file_name,"w")
    classifi = []
    not_realized = []
    for k in L:
        print("Invariant Lattice number %s"%k)
        fix = fixed[k]
        if fix.rank() != 4:
            continue
        print(fix.gram_matrix())
        print(" ")
        cofix = cofixed[k].twist(-1)
        Hacts = k3_prime_power(fix.genus(), 2, 1)
        Hacts += k3_prime_power(fix.genus(), 2, 2)
        Hacts += k3_pq(fix.genus(),2,3)
        for a in Hacts:
            fix_a = a[0]
            g = a[1]
            try:
                aut = MaximalK3surfaceAut(fix_a,cofix, g)
                classifi.append(aut)
                s = aut.str()
                result.write(s+ "\n")
            except ValueError:
                not_realized.append([fix,a])
    result.close()
    return classifi, not_realized, Hacts
###########################################################
p = 3
ord3 = [(fixed[k], cofixed[k], p) for k in range(len(fixed)) if fixed[k].rank()==4]
