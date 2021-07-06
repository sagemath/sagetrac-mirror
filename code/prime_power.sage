from sage.modules.free_quadratic_module_integer_symmetric import FreeQuadraticModule_integer_symmetric
from sage.interfaces.gap import get_gap_memory_pool_size
memory_gap = get_gap_memory_pool_size()
set_gap_memory_pool_size(9048*memory_gap)
libgap.eval("SetRecursionTrapInterval(10000000)");
attach("hermitian.sage")

def k3_prime_power(genus, prime, e, verbose=False,rkT=None):
    r"""
    Returns all conjugacy classes of isometries of order `p^e` in `genus`
    such that the kernel of the `p^e`-part is of signature `(2,n)`.

    INPUT:

    - ``genus`` -- a genus of signature `(3,k)`.
    - ``p`` -- a prime number
    - ``e`` -- a natural number
    """
    assert genus.signature_pair()[0] == 3
    rk = genus.rank()
    weights = [euler_phi(d) for d in (prime^e).divisors()]
    n = len(weights)
    if prime == 2 and e==1:
        floor = [1, 2]     # minimal ranks  due to ample class and 2 form
    else:
        floor = [1] + [0]*(n-2) + [1]
    ceiling = [rk // d for d in weights]
    P = IntegerListsLex(max_sum=rk,length=n,floor=floor, ceiling=ceiling)
    P = [[p[k]*weights[k] for k in range(n)] for p in P]
    P = [p for p in P if sum(p)==rk]
    acts = []
    for ranks in P:
        if rkT is not None and ranks[-1]!=rkT:
            continue
        ranks_E = [ranks[k]//weights[k] for k in range(n)]
        if prime == 2 and e != 1:
            signatures = [[ranks_E[0]-1], [ranks_E[1]]]
            signatures += [[ranks_E[k]]*(weights[k]//2) for k in range(2,n)]
        elif prime == 2 and e == 1:
            signatures = [[ranks_E[0]-1], [ranks_E[1]-1]]
        else:
            signatures = [[ranks_E[0]-1]]
            signatures += [[ranks_E[k]]*(weights[k]//2) for k in range(1,n)]
        signatures[-1][0] -= 1
        for act in prime_power_actions(genus,prime,ranks,signatures,k3_unobstructed=True,verbose=verbose):
            yield act

def enriques_prime_power(genus, prime, e, verbose=False,rkT=None):
    r"""
    Returns all conjugacy classes of isometries of order `p^e` in `genus`
    such that the kernel of the `p^e`-part is of signature `(2,n)`.

    INPUT:

    - ``genus`` -- a genus of signature `(3,k)`.
    - ``p`` -- a prime number
    - ``e`` -- a natural number
    """
    assert genus.signature_pair()[0] == 2
    rk = genus.rank()
    weights = [euler_phi(d) for d in (prime^e).divisors()]
    n = len(weights)
    floor = [0]*n
    ceiling = [rk // d for d in weights]
    P = IntegerListsLex(max_sum=rk,length=n,floor=floor, ceiling=ceiling)
    P = [[p[k]*weights[k] for k in range(n)] for p in P]
    P = [p for p in P if sum(p)==rk]
    acts = []
    for ranks in P:
        if rkT is not None and ranks[-1]!=rkT:
            continue
        ranks_E = [ranks[k]//weights[k] for k in range(n)]
        if prime == 2 and e != 1:
            signatures = [[ranks_E[0]], [ranks_E[1]]]
            signatures += [[ranks_E[k]]*(weights[k]//2) for k in range(2,n)]
        elif prime == 2 and e == 1:
            signatures = [[ranks_E[0]], [ranks_E[1]+2]]
        else:
            signatures = [[ranks_E[0]]]
            signatures += [[ranks_E[k]]*(weights[k]//2) for k in range(1,n)]
        for k in range(0,n):
            if k == 0 or (k == 1 and prime==2):
                eps = 2
            else:
                eps =1
            if signatures[k][0] >= eps:
                signatures[k][0] -= eps
        for act in prime_power_actions(genus,prime,ranks,signatures,k3_unobstructed=True,verbose=verbose):
            yield act


def is_admissible(A, B, C, p):
    r"""
    Return if A+B could glue to C under the condition that
    pC <= A + B.

    INPUT:

    - A,B,C -- Genera

    OUTPUT:

    bool; if False A and B cannot glue to C

    EXAMPLES::

    """
    from sage.quadratic_forms.genera.genus import Genus_Symbol_p_adic_ring
    primes = [sym.prime() for sym in C.local_symbols() if sym.prime() != p]

    for sym in A.local_symbols()+B.local_symbols():
        if sym.prime()!=p and sym.prime() not in primes:
            return False

    glue = ZZ(A.det() * B.det() / C.det()).abs()
    if not glue.is_square():
        return False
    if len(glue.prime_factors()) > 1:
        return False
    g = glue.valuation(p) // 2

    if not C.scale().divides(gcd(A.scale(),B.scale())):
        return False

    if not C.norm().divides(gcd(A.norm(),B.norm())):
        return False

    # check that (A + B)_q = C_q for all q != p
    AB = A.direct_sum(B)
    if not all(AB.local_symbols(q)==C.local_symbols(q)
               for q in primes):
        return False
    # if there is no glue, then A + B = C
    if g==0:
        return AB.local_symbols(p)==C.local_symbols(p)

    if not lcm(A.level(),B.level()).divides(C.level()*p):
        return False

    # if C_p is unimodular we have a necessary and sufficient condition
    if not p.divides(C.det()):
      qA = A.discriminant_form().primary_part(p).normal_form().gram_matrix_quadratic()
      qB = B.discriminant_form().primary_part(p).twist(-1).normal_form().gram_matrix_quadratic()
      return qA == qB

    # check that rationally (A + B)_p = C_p
    # since the determinants match the excess is sufficient
    if AB.local_symbols(p).excess() != C.local_symbols(p).excess():
        return False
    #return True # for bug finding
    Ap = A.local_symbols(p).symbol_tuple_list()
    Bp = B.local_symbols(p).symbol_tuple_list()
    Cp_level = C.level().valuation(p)
    m = Cp_level + 1
    a_max = sum(s[1] for s in Ap if s[0]==m)
    b_max = sum(s[1] for s in Bp if s[0]==m)
    # the ranks of the jordan components
    # of maximum possible scales must match
    if a_max != b_max:
        return False
    a1 = sum(s[1] for s in Ap if s[0] == 1)
    a2 = sum(s[1] for s in Ap if s[0] >= 2)
    b1 = sum(s[1] for s in Bp if s[0] == 1)
    b2 = sum(s[1] for s in Bp if s[0] >= 2)
    if Cp_level != 0:
        ker_min = max(g - a1, g - b1, a_max)
    else:
        ker_min = max(g - a1, g - b1)
    ker_max = min(a2 + a1//2, b2 + b1//2, g)
    # compare kernel dimensions of the glue
    # for the existence of a glue map
    if ker_max < ker_min:
        return False


    ABp = AB.local_symbols(p).symbol_tuple_list()
    Cp = C.local_symbols(p).symbol_tuple_list()


    # if the glue is exactly p^l D_i,
    # then the parts of A and B
    # of scale <= p^(l-1)
    # are basically those of C
    # as they are not glued.
    # we have to be careful at p=2
    # because of sign walking etc.
    if a_max == g:
        if len(Ap)>1:
            Ar = Genus_Symbol_p_adic_ring(p,Ap[:-1])
        else:
            Ar = Genus(matrix([])).local_symbols(p)
        if len(Bp)>1:
            Br = Genus_Symbol_p_adic_ring(p,Bp[:-1])
        else:
            Br = Genus(matrix([])).local_symbols(p)

        ABr = Ar.direct_sum(Br)
        for i in range(Cp_level):
            s1 = ABr.symbol(i)
            s2 = C.local_symbols(p).symbol(i)
            # scales must match
            if s1[1] != s2[1]:
                return False
            # if the jordan component of ABr is odd
            # that of C is odd as well
            if p == 2 and s1[3] > s2[3]:
                return False
            if p != 2 and s1[2] != s2[2]:
                return False

    # p=2
    # if the p^level(C) modular component of C is even
    # then the p^level(C)+1 modular components of A and B
    # have the same pariy
    s3 = C.local_symbols(2).symbol_tuple_list()[-1]
    if p==2 and s3[3]==0:
        s1 = A.local_symbols(2).symbol(s3[0]+1)
        s2 = B.local_symbols(2).symbol(s3[0]+1)
        if s1[3] != s2[3]:
            print('parity')
            return False


    # AB --> C --> 1/p AB cap AB^/vee
    # we use the equivalent
    # pAB --> pC --> AB cap p AB^/vee
    pAB = AB.local_symbols(p).symbol()
    for s in pAB:
        # we rescale AB by p. Thus its gram matrix
        # is rescales by p^2
        s[0] += 2
    pAB = Genus_Symbol_p_adic_ring(p, pAB)


    pC = C.local_symbols(p).symbol()
    for s in pC:
        s[0] += 2
    pC = Genus_Symbol_p_adic_ring(p, pC)

    ab = AB.local_symbols(p).symbol()
    ab0 = [s for s in ab[:1] if s[0] == 0]
    ab1 = [s for s in ab if s[0] >= 1]
    # only the unimodular component of
    # AB cap p AB^/vee changes
    for s in ab0:
        s[0] += 2
    ab0 = Genus_Symbol_p_adic_ring(p, ab0)
    ab1 = Genus_Symbol_p_adic_ring(p, ab1)
    ab = ab1.direct_sum(ab0)
    assert ab.represents(pAB)
    if not (ab.represents(pC) and pC.represents(pAB)):
        return False
    assert C.local_symbols(p).represents(AB.local_symbols(p))

    for i in range(Cp_level+1):
        r1 = sum(s[1] for s in Cp if s[0] <= i-2)
        r2 = sum(s[1] for s in ABp if s[0] <= i)
        r3 = sum(s[1] for s in Cp if s[0] <= i)
        if not r1 <= r2 <= r3:
            print('level obstruction')
            return False

    r1 = sum(s[1] for s in Cp if s[0] <= Cp_level-1)
    r2 = sum(s[1] for s in ABp if s[0] <= Cp_level)
    r3 = sum(s[1] for s in Cp if s[0] <= Cp_level)
    if not r1 <= r2 <= r3:
        print('new level obstruction')
        return False

    # check that the discriminant group
    # of Cp is a subquotient of that of ABp
    # wierd this should be automatic
    # if AB --> C
    if any(r1 > r2 for r1,r2 in
            ((sum(s[1] for s in Cp[k:]),sum(s[1] for s in ABp if s[0]>= Cp[k][0]))
             for k in range(0,len(Cp)))
          ):
        assert False
        print('subquotient')
        return False

    return True


def prime_power_actions(genus, p, ranks, signatures, k3_unobstructed=True,verbose=False):
    r"""
    Returns all conjugacy classes of isometries of order `p^e` in `genus`.

    EXAMPLES::

        sage: target = all_genera_by_det((3,7),2^2)[0]
        sage: p = 3
        sage: ranks = [2, 2, 6]
        sage: signatures = [[0], [1], [0,1,1]]
        sage: prime_power_actions(target, p, ranks, signatures)

    """
    assert p.is_prime()
    G = genus

    i = min(i for i in range(len(ranks)+1) if all(g==0 for g in ranks[i:]))
    if len(ranks)>0 and ranks[-1]==0:
        # the higher orders have rank 0
        # and can be skipped
        ranks = ranks[:i]
        signatures = signatures[:i]
    if len(ranks)==0:
        yield 3*[matrix.identity(0)]
        return
    if len(ranks) == 1:
        # there is nothing to glue
        for M in genus.representatives():
            fM = M ^ 0
            M = LatticeWithIsometry(IntegralLattice(M),fM,order=1)
            GM = M.Oq_equiv()
            yield [M.L, M.iso, GM]
            return
    M_rank = ranks[-1]
    M_signature = signatures[-1]
    R_signatures = signatures[:-1]


    e = len(ranks) - 1
    assert e>=1
    deg_E = euler_phi(p**e)
    M_rk_E = M_rank//deg_E


    R_rank = sum(ranks[:-1])
    if p==2:
        R_neg = sum(flatten(R_signatures[:2])) + 2*sum(flatten(R_signatures[2:]))
    else:
        R_neg = sum(flatten(R_signatures[:1])) + 2*sum(flatten(R_signatures[1:]))
    R_pos = R_rank - R_neg


    max_level = p * G.level()
    min_scale = G.scale()
    min_norm = G.norm()
    M_min_det = min_scale^M_rank
    m = min(M_rk_E*p**(e-1), R_rank) # a bit crude
    M_max_det = p^m*G.det()

    primes = [sym.prime() for sym in G.local_symbols() if sym.prime() != p]
    countM = 0
    countR = 0
    results = []
    for M_det in [d for d in M_max_det.divisors() if M_min_det.divides(d)]:
        for Mh in all_lattice_with_isometry(p^e, M_rk_E, M_det, max_level,
                                           signatures=M_signature,
                                           min_scale=min_scale,min_norm=min_norm):
            cm = sage.structure.element.get_coercion_model()
            cm.reset_cache()
            M = Mh.L
            fM = Mh.iso
            pos, neg = M.signature_pair()
            # In the coinvariant lattice in NS there are no roots
            if k3_unobstructed and all([k3_unobstructed, e>0, pos==0, neg>0,
                   M.maximum()==-2]):
                continue
            if verbose:
                countM += 1
                print("M,fM: %s"%countM)
                print(M.genus())
            GM = Mh.Oq_equiv()
            DM = M.discriminant_group()
            M_max_glue = (M.span((fM^(p^(e-1)) - fM**0).inverse())
                          & G.scale()*M.dual_lattice())
            M_max_glue = DM.submodule(M_max_glue.gens())


            for glue_order in min(p^m,M_max_glue.cardinality()).divisors():
                R_det = G.det() * glue_order^2 / M.determinant()
                if R_det.denominator() != 1:
                    continue
                R_det = ZZ(R_det)
                # the glue group is too small
                if not glue_order.divides(R_det):
                    continue
                for R in all_genera_by_det((R_pos, R_neg),R_det,
                                           max_scale=max_level,
                                           even=G.is_even()):
                    if not is_admissible(M.genus(),R,G,p):
                        continue
                    # now MR is a serious candidate
                    if verbose:
                        countR += 1
                        print("R: %s %s"%(countR,R.local_symbols(p)))
                    # recurse
                    for N, fN, GN in prime_power_actions(R, p, ranks[:-1], R_signatures,
                                                 k3_unobstructed=k3_unobstructed,
                                                 verbose=verbose):
                        ext = extensions(M, fM, N, fN, GM, GN,
                                         glue_order, p,
                                         target_genus=genus)
                        if verbose:
                            print("Found %s matching extensions"%len([ex for ex in ext if ex[0].genus()==G]))
                        for ex in ext:
                            if ex[0].genus()==G:
                                yield ex

                        #if not any(t[0].genus() == G for t in ext):
                        #    print("dead end at level %s"%e)

def ptype(L, f, p):
    r"""

    EXAMPLES::

        sage: gram = matrix(ZZ,8,
        ....:[ (2, 0, -1, 0, 0, 0, 0, 0),
        ....:  (0, 2, 0, -1, 0, 0, 0, 0),
        ....:  (-1, 0, 2, -1, 0, 0, 0, 0),
        ....:  (0, -1, -1, 2, -1, 0, 0, 0),
        ....:  (0, 0, 0, -1, 2, -1, 0, 0),
        ....:  (0, 0, 0, 0, -1, 2, -1, 0),
        ....:  (0, 0, 0, 0, 0, -1, 2, -1),
        ....:  (0, 0, 0, 0, 0, 0, -1, 2)])
        sage: L = IntegralLattice(gram)
        sage:

    """
    g = L.orthogonal_group([f]).gen(0)
    e = g.order().valuation(p)
    # we want the order of f on L!!
    for k in range(1, e+1):
        if all(v==v*g^(p^k) for v in L.gens()):
            break
    else:
        raise AssertionError
    e = k

    R.<x> = ZZ[]
    ptype = []
    for e0 in range(e,0, -1):
        c = cyclotomic_polynomial(p**e0, x)
        C = L.kernel_sublattice(c(f))
        r = x**(p**(e0-1)) - 1
        R = L.kernel_sublattice(r(f))
        glue = (C + R).index_in(L & (C + R).base_extend(QQ))
        # we want the isomorphism class of the glue
        # as a quadratic form
        # We need the projection to C, say
        CR = L.sublattice(C.gens() + R.gens())
        gens = [CR.coordinates(g) for g in L.gens()]
        n = len(C.gens())
        Cv = C.vector_space()
        gens = [Cv.linear_combination_of_basis(g[:n]) for g in gens]
        glueC = C.discriminant_group().submodule(gens)
        glueC = glueC.normal_form()
        q1 = glueC.gram_matrix_quadratic()
        q2 = glueC.twist(-1).normal_form().gram_matrix_quadratic()

        ptype.append([(p,e0),L.genus(),C.genus(),R.genus(),ZZ(glue),[q1,q2]])
        L = R
    ptype.append([(p,0),L.genus(),L.genus(),Genus(matrix([])),ZZ(1),matrix([])])
    ptype.reverse()
    return ptype


def prime_order(p, genus, k3=True,verbose=2,rankCp=None):
    r"""
    """
    G = genus
    pos, neg = G.signature_pair()
    if not 0 <= pos <= 3:
        raise ValueError('number of positive squares must be 0, or 2 or 3')
    rk = genus.rank()
    # order 1 action
    for C1 in G.representatives():
        C1 = IntegralLattice(C1)
        f1 = matrix.identity(rk)
        GC1 = C1.image_in_Oq()
        yield C1, f1, GC1

    countCp = 0
    countC1 = 0

    # fix the rank of Cp
    for k in range(1, rk//(p-1)):
        rkCp = k * (p-1)
        if rankCp is not None and rkCp !=rankCp:
            continue
        if pos > 1 and p != 2:
            sigp = [k-1] + ((p-1)//2 - 1) *[k]
        if pos > 1 and p == 2:
            if k==1:
                continue
            sigp = [k-2]
        if pos <= 1 and p != 2:
            sigp = ((p-1)//2) *[k]
        if pos <= 1 and p == 2:
            sigp = [k]
        rkC1 = rk - rkCp
        detCpmin = G.scale()^rkCp
        m = min(k, rkC1) # a bit crude
        detCpmax = p^m*G.det()
        # fix the determinant of Cp
        for detCp in [d for d in detCpmax.divisors() if detCpmin.divides(d)]:
            cm = sage.structure.element.get_coercion_model()
            cm.reset_cache()
            # fix Cp
            # fix the glue order
            for glue_order in gcd(p^m, detCp).divisors():
                detC1 = G.det() * glue_order^2 / detCp
                if detC1.denominator() != 1:
                    continue
                detC1 = ZZ(detC1)
                # the glue group is too small
                if not glue_order.divides(detC1):
                    continue
                for CpG in all_lattice_with_isometry(p, k, detCp, p*G.level(),
                                                    signatures=sigp,
                                                    min_scale=G.scale(),
                                                    min_norm=G.norm(),
                                                    return_genera=True):
                    try:
                        Cp_genus = CpG.trace()
                    except AttributeError:
                        Cp_genus = CpG
                    posCp, _ = Cp_genus.signature_pair()
                    sigCp = Cp_genus.signature_pair()
                    sigC1 = pos - sigCp[0], neg - sigCp[1]
                    for genusC1 in all_genera_by_det(sigC1, detC1,
                                            max_scale=p*G.level(),
                                            even=G.is_even()):

                        if not is_admissible(Cp_genus, genusC1, G, p):
                            continue
                        # now MR is a serious candidate
                        # recurse
                        f1 = matrix.identity(rkC1)
                        if verbose > 2:
                            print('computing representatives of %s'%CpG)
                        for CpE in CpG.representatives():
                            Cp, fp = trace_lattice(CpE)
                            Cp = IntegralLattice(Cp)
                            # In the coinvariant lattice in NS there are no roots
                            if k3 and posCp==0 and Cp.maximum()==-2:
                                continue
                            if verbose > 0:
                                countCp  += 1
                                print("Cp,fp: %s"%countCp)
                                print(Cp_genus)

                            if verbose > 0:
                                countC1 += 1
                                print("C1: %s %s"%(countC1,genusC1.local_symbols(p)))
                            Cph = LatticeWithIsometry(Cp,fp,order=p,gramE=CpE,magmaRep=CpG.representative())
                            GCp = Cph.Oq_equiv()
                            if verbose > 2:
                                print('computing representatives')
                            for C1 in genusC1.representatives():
                                C1 = IntegralLattice(C1)
                                GC1 = C1.image_in_Oq()
                                ext = extensions(Cp, fp, C1, f1, GCp, GC1,
                                                glue_order, p, target_genus=G)
                                if verbose > 0:
                                    print("Found %s matching extensions"%len(ext))
                                for ex in ext:
                                    yield ex


def isometries(n, genus):
    r"""
    Return all `\Phi_n`-lattices in genus up to isomorphism and powers.

    The signatures are chosen K3 non-symplectic compatible.

    INPUT:

    - ``n`` -- integer `>0`
    - ``genus`` -- a global genus symbol of signatures ``(pos, neg)`` with ``pos <=2``

    OUTPUT:

    a list of pairs (M,fM) where ``M`` is the gram matrix of a lattice and
    ``fM`` the matrix of an isomety of ``M``
    """

    try:
        sig = sig_k3(n,genus.rank(),genus.signature_pair()[0])
    except ValueError:
        return []
    rkE = genus.rank()/euler_phi(n)
    if rkE.denominator() != 1:
        return []
    L = all_lattice_with_isometry(n, rkE, genus.det(), genus.level(),
                                     signatures = sig, min_scale=genus.scale())
    return [M for M in L if M.L.genus()==genus ]

def has_given_invariant_submodule(M, fM, OqfM, qlist, p):
    r"""
    """
    Of = M.orthogonal_group([fM])
    n = Of.order()
    g = OqfM(Of.gen(0))
    DM = M.discriminant_group()
    H = M.dual_lattice() & M.span((fM^(n/p) - 1).inverse())
    H = DM.submodule(H.gens())
    glue_order = p^qlist[0].ncols()
    sg = DM.subgroup_representatives(H, OqfM, g=g, order=glue_order)
    sg = [d.representative().normal_form().gram_matrix_quadratic() for d in sg]
    return any(s in qlist for s in sg)


def next_prime_power(ptype, verbose=0):
    r"""
    """
    i = max(i for i in range(len(ptype)) if ptype[i][2].rank()!=0)
    p, e = ptype[i][0]
    genus = ptype[i][1]
    Cgenus = ptype[i][2]
    Rgenus = ptype[i][3]
    pglue =  ptype[i][4]
    qlist = ptype[i][5]

    if e == 0:
        assert genus == Cgenus
        assert Rgenus.rank() == 0
        for x in prime_order(p, genus, k3=True, verbose=verbose):
            yield x
        return

    for Ch in isometries(p**(e+1), Cgenus):
        fC = Ch.iso
        GC = Ch.Oq_equiv()
        C = Ch.L
        # fight a memory leak
        cm = sage.structure.element.get_coercion_model()
        cm.reset_cache()
        # before we recurse we want to make sure
        # that we could actually glue equivariantly
        if not has_given_invariant_submodule(C, fC, GC, qlist, p):
            continue
        # recurse
        for R, fR, GR in next_prime_power(ptype[:i],verbose=verbose-1):
            ext = extensions(C, fC, R, fR, GC, GR,
                             pglue, p, target_genus=genus,
                             qlist=qlist)
            if verbose > 0:
                print('found %s extensions at level %s'%(len(ext),verbose))
            for ex in ext:
                yield ex

@parallel(ncpus=12)
def pnq_actions_pure(q,ptype, splitsig, k3_unobstructed):
    #magma.attach_spec("sage/code/unit_quotients/lat.spec")
    load("/home/lehrstuhl/ag-brandhorst/brandhorst/sage/code/hermitian.sage")
    acts_iter = pnq_actions(q,ptype,k3_unobstructed=k3_unobstructed,splitsig=splitsig)
    acts = [(g[0],g[1]) for g in acts_iter]
    return acts

def pnq_actions(q, ptype,k3_unobstructed=True,verbose=2,splitsig=None):
    r"""
    Returns all conjugacy classes of isometries of order `p^e` in `genus`.

    EXAMPLES::

        sage: target = all_genera_by_det((3,7),2^2)[0]
        sage: p = 3
        sage: ranks = [2, 2, 6]
        sage: signatures = [[0], [1], [0,1,1]]
        sage: prime_power_actions(target, p, ranks, signatures)

    """  
    # skip those i where Cp^iq + Cp^i is of rank 0
    i = max(i for i in range(len(ptype)) if ptype[i][2].rank()!=0)
    p, e = ptype[i][0]
    genus = ptype[i][1]
    Cgenus = ptype[i][2]
    Rgenus = ptype[i][3]
    pglue =  ptype[i][4]
    qlist = ptype[i][5]

    # find all actions p^e q of types p^e
    # that glue to genus
    if splitsig is not None and len(splitsig)>0:
        splitsig_i = splitsig[i]
        splitsig_new = splitsig[:i]
    else:
        splitsig_i = None
        splitsig_new = None
    for (M, fM, GM) in splitpq(Cgenus,p,e,q,k3_unobstructed=k3_unobstructed,verbose=verbose,splitsig=splitsig_i):
        # base case
        cm = sage.structure.element.get_coercion_model()
        cm.reset_cache()
        if e == 0:
            yield (M, fM, GM)
            continue
        if not has_given_invariant_submodule(M, fM, GM, qlist, p):
            continue
            print('does not have the given invariant submodule')
        if verbose >= 0:
            print('%s^%s * %s -- %s^%s'%(p,e,q,p,e))
            print(M.genus())
        # recurse
        for (N, fN, GN) in pnq_actions(q,ptype[:i],k3_unobstructed=k3_unobstructed,verbose=verbose-1,splitsig=splitsig_new):
            ext = extensions(M, fM, N, fN, GM, GN,
                             pglue, p, target_genus=genus,
                             qlist=qlist)
            for ex in ext:
                assert ex[0].genus() == genus
                if k3_unobstructed:
                    P = PolynomialRing(QQ,"x")
                    x = P.gen()
                    pol = (x^(p^e*q)-1)/(x-1)
                    K = ex[0].kernel_sublattice(pol(ex[1]))
                    if K.signature_pair()[0]!=0: 
                        pol = (x^(p^e*q)-1)/(x-1)/cyclotomic_polynomial(p^e*q,x)
                        K = ex[0].kernel_sublattice(pol(ex[1]))
                    if K.signature_pair()[0]==0 and K.maximum()==-2:
                        print("obstructed")
                        continue
                yield ex

def sig_k3(n, rk, pos):
    r"""
    """
    neg = rk - pos
    degE = euler_phi(n)
    if n == 1:
        assert pos <=1 # invariant ample class
        return [neg]
    if n == 2:
        assert pos in [0,2]
        return [neg]
    else:
        rkE = rk//degE
        assert pos in [0,2]
        if not neg % 2 == 0:
            raise ValueError
        if pos == 0:
            return [rkE]*(degE//2)
        sig = [rkE - 1]
        if degE > 2:
            sig += [rkE]*(degE//2 - 1)
        return sig


def split_sig(sig_pair, p, e, q):
    r"""

        sage: [g for g in split((3,5),2,0,3)]
        [[[6, 2], [[5], [0]]], [[4, 4], [[3], [1]]], [[2, 6], [[1], [2]]]]
    """
    assert p!=q
    pos, neg = sig_pair
    rk = pos + neg
    assert pos in [0,1,2,3]
    weights = [euler_phi(p^e),euler_phi(p^e*q)]
    n = 2
    ceiling = [rk // d for d in weights]
    floor = [0, 0]
    if pos >= 2:
        # the order must be p^eq
        floor[1] = 1
    if pos % 2 == 1:
        # there must be an invariant vector of pos. square
        if e != 0:
            return
        floor[0] = 1
    R = IntegerListsLex(max_sum=rk,length=n,floor=floor, ceiling=ceiling)
    R = [[r[k]*weights[k] for k in range(n)] for r in R]
    R = [r for r in R if sum(r)==rk]
    for ranks in R:
        if e == 0:
            assert pos >=1
            sigp = sig_k3(p**e, ranks[0], 1)
        else:
            assert pos in [0,2]
            sigp = sig_k3(p**e, ranks[0], 0)
        pospq = (pos//2)*2
        sigpq = sig_k3(p**e * q, ranks[1], pospq)
        signatures = [sigp,sigpq]
        yield [ranks, signatures]

def glue_det(det, charpol1, charpol2):
    r"""
    """
    # we assume separable minimal polynomial here
    g1 = charpol1.radical().resultant(charpol2)
    g2 = charpol2.radical().resultant(charpol1)
    glue_max = ZZ(g1).gcd(ZZ(g2))
    for glue in glue_max.divisors():
        # (det1 /glue)* (det2 / glue) == det
        for det1 in (glue*det).divisors():
            if not glue.divides(det1):
                continue
            det2 = det * glue^2 / det1
            try:
                det2 = ZZ(det2)
            except ValueError:
                continue
            if not glue.divides(det2):
                continue
            yield glue, det1, det2



def splitpq(genus, p, e, q, k3_unobstructed=True,verbose=0,splitsig=None):
    r"""
    Classify lattices with isometries up to isomorphism
    such that the isometry has minimal polynomial `\Phi_{p^eq}\Phi_{p^e}`
    and the lattice is in ``genus``.

    INPUT:

    - ``p`` -- a prime number
    - ``e`` -- an exponent
    - ``p`` -- a prime number
    - ``k3_unobstructed`` -- (default: ``True``)

    OUTPUT:

    A list consisting of triples ``[A, a, Oa]`` where
    ``A`` is a lattice in ``genus``,
    ``a`` is an isomety in ``a`` and
    ``Oa`` is the image in the discriminant group
    of the stabilizer of ``a``.
    """
    assert p.is_prime() and q.is_prime()
    primes = [sym.prime() for sym in genus.local_symbols() if sym.prime() != q]

    if splitsig is None:
        splitsig = split_sig(genus.signature_pair(),p,e,q)
    for rks,sigs in splitsig:
        degsE = euler_phi(p**e), euler_phi(q*p**e)
        rksE = [rks[0]//degsE[0],rks[1]//degsE[1]]

        if rks[0] == 0:
            # there is nothing to glue
            for Mh in all_lattice_with_isometry(q*p^e,rksE[1],genus.det(),
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
            for Rh in all_lattice_with_isometry(p^e,rksE[0],genus.det(),
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
                                           c(p^e)^rksE[0],
                                           c(p^e*q)^rksE[1]):
            alwi = all_lattice_with_isometry
            for MG in alwi(q*p^e, rksE[1], M_det, max_level,
                           signatures=sigs[1],
                           min_scale=min_scale,
                           min_norm=min_norm,
                           return_genera=True):
                try:
                    Mgenus = MG.trace()
                except AttributeError:
                    Mgenus = MG
                for RG in alwi(p^e, rksE[0], R_det, max_level,
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
                    for ME in MG.representatives():
                        M, fM = trace_lattice(ME, order=q*p**e)
                        M = IntegralLattice(M)
                        Mh = LatticeWithIsometry(M,fM,order=q*p**e,gramE=ME,magmaRep=MG.representative())
                        if k3_unobstructed and is_obstructed(M, fM):
                            continue
                        GM = Mh.Oq_equiv()
                        cm = sage.structure.element.get_coercion_model()
                        cm.reset_cache()
                        if verbose>2:
                            print('computing representatives of %s'%RG)
                        for RE in RG.representatives():
                            R, fR = trace_lattice(RE, order=p**e)
                            R = IntegralLattice(R)
                            Rh = LatticeWithIsometry(R,fR,order=p**e,gramE=RE,magmaRep=RG.representative())
                            if k3_unobstructed and is_obstructed(R, fR):
                                continue
                            GR = Rh.Oq_equiv()
                            if verbose >0:
                                print('glueing')
                                print(M.genus())
                                print('and')
                                print(R.genus())
                                print('above %s'%q)
                            ext = extensions(M, fM, R, fR, GM,
                                         GR, glue, q,
                                         target_genus=genus)
                            for ex in ext:
                                assert ex[0].genus() == genus
                                yield ex



def is_obstructed(M, fM):
    pos, neg = M.signature_pair()
    assert fM.minpoly().is_irreducible()
    if pos > 0:
        return False
    if neg == 0:
        return False
    if M.maximum() == -2:
        return True



def extensions(M, fM, N , fN, GM, GN, glue_order, p,
               DM_max_glue=None,DN_max_glue=None, target_genus=None,
               qlist=None):
    r"""
    Return all possible extensions of `(M, fM)` by (N,fN)`
    modulo `GM x GN` along the prime `p`.

    OUTPUT:


    EXAMPLES::

        sage: M = IntegralLattice(Matrix([4]))
        sage: N = IntegralLattice(Matrix([52]))
        sage: fM = matrix([-1])
        sage: fN = matrix([1])
        sage: GM = M.image_in_Oq()
        sage: GN = N.image_in_Oq()
        sage: extensions(M,fM,N,fN,GM,GN,2,2)
        [[
        Lattice of degree 2 and rank 2 over Integer Ring
        Basis matrix:
        [1/2 1/2]
        [  0   1]
        Inner product matrix:                             [-1| 0]
        [ 4  0]                                           [--+--]
        [ 0 52]                                         , [ 0| 1],
        Group of isometries of
        Finite quadratic module over Integer Ring with invariants (2, 26)
        Gram matrix of the quadratic form with values in Q/2Z:
        [ 3/2  1/2]
        [ 1/2 1/13]
        generated by 2 elements
        ]]

    """
    from sage.modules.free_quadratic_module_integer_symmetric import FreeQuadraticModule_integer_symmetric
    assert isinstance(M,FreeQuadraticModule_integer_symmetric)
    assert isinstance(N,FreeQuadraticModule_integer_symmetric)
    assert M.degree() == fM.ncols()
    assert N.degree() == fN.ncols()
    results = []
    glue_val = glue_order.valuation(p)

    # TODO: sloppy
    M_ord = MatrixGroup(fM).order()
    N_ord = MatrixGroup(fN).order()
    assert N_ord != M_ord
    #assert N_ord.divides(M_ord)
    assert 1 == gcd(fM.charpoly(),fN.charpoly())
    r = ZZ(fM.charpoly().resultant(fN.charpoly()))
    assert r == 1 or [ZZ(p)] == r.prime_factors()

    DM = M.discriminant_group()
    DN = N.discriminant_group()

    # part of order at most ord(fN)
    if DM_max_glue is None:
        DM_max_glue = M.span((fM**(N_ord) - fM**0).inverse()) & M.dual_lattice()
        DM_max_glue = DM.submodule(DM_max_glue.gens())
    if DN_max_glue is None:
        # any vector is okay as long as it is of norm at most p
        DN_max_glue = DN.submodule([g*(g.order()//p) for g in DN.primary_part(p).gens()])
    H10=None
    H20=None
    if target_genus is not None:
        level = target_genus.level()
        H10 = level*DM
        H20 = level*DN
    fDM = DM.orthogonal_group()(M.orthogonal_group([fM]).gen(0))
    fDN = DN.orthogonal_group()(N.orthogonal_group([fN]).gen(0))
    f = matrix.block_diagonal([fM,fN])
    if min(DM_max_glue.cardinality(),DN_max_glue.cardinality()) < glue_order:
        return []
    for glue, stab in DM.all_primitive_prime_equiv(DN, DM_max_glue, DN_max_glue,
                                                   GM, GN, fDM, fDN, glue_val,
                                                   H10=H10,H20=H20,
                                                   target_genus=target_genus,
                                                   qlist=qlist):
        if hasattr(glue.W(),"overlattice"):
            W = glue.W()
        else:
            W = FreeQuadraticModule_integer_symmetric(
                                    ambient=glue.W().ambient_module(),
                                    basis=glue.W().basis(),
                                    inner_product_matrix=glue.W().inner_product_matrix(),
                                    already_echelonized=False)
        ext = W.overlattice([g.lift() for g in glue.gens()])
        assert ext*f == ext
        ext_Oq = ext.discriminant_group().orthogonal_group()

        stab = ext_Oq.subgroup([ext_Oq(g) for g in stab.gens()])
        results.append([ext, f, stab])
    return results

