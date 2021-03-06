from sage.interfaces.gap import get_gap_memory_pool_size
memory_gap = get_gap_memory_pool_size()
set_gap_memory_pool_size(9048*memory_gap)
libgap.eval("SetRecursionTrapInterval(10000000)");
attach("k3_classification/hermitian.sage")

# TODO:
# Turn things into iterators whenever possible
# tests tests tests tests

def positve_definite_by_orthogonal(genus,prime,e):
    r"""
    Returns all conjugacy classes of isometries of order ``prime^e``
    in ``genus``.

    The computation is using the enumeration of genera and
    computation of orthogonal groups.

    Input:

    - genus -- a positive definite genus
    - prime -- integer
    - e -- integer

    Output:

    - a
    """
    actions = []
    for R in genus.representatives():
        R = IntegralLattice(R)
        G = R.orthogonal_group()
        conj = [g for g in G.conjugacy_classes_representatives() if g.order()==prime^e]
        for f in conj:
            actions.append(LatticeWithIsometry(R,f.matrix()))
    return actions

def positve_definite_prime(genus, prime, e):
    r"""
    Returns all conjugacy classes of isometries of order ``prime^e``
    in ``genus``.

    The computation is using hermitian lattices and equivariant glueings.

    Input:

    - genus -- a positive definite genus
    - prime -- integer
    - e -- integer

    Output:



    """
    assert genus.signature_pair()[1] == 0, 'genus must be positive definite'
    rk = genus.rank()
    weights = [euler_phi(d) for d in (prime^e).divisors()]
    n = len(weights)
    floor = [0] + [0]*(n-2) + [1]
    ceiling = [rk // d for d in weights]
    P = IntegerListsLex(max_sum=rk,length=n,floor=floor, ceiling=ceiling)
    P = [[p[k]*weights[k] for k in range(n)] for p in P]
    P = [p for p in P if sum(p)==rk]
    acts = []
    for ranks in P:
        ranks_E = [ranks[k]//weights[k] for k in range(n)]
        if prime == 2:
            signatures = [[0], [0]]
            signatures += [[0]*(weights[k]//2) for k in range(2,n)]
        else:
            signatures = [[0]]
            signatures += [[0]*(weights[k]//2) for k in range(1,n)]
        acts += prime_power_actions(genus,prime,ranks,signatures)
    results = []
    for act in acts:
        L = act[0]
        f = act[1]
        iso = LatticeWithIsometry(L,f)
        results.append(iso)
    return results


def check_prime_power(genus, prime, e):
    r"""
    Check that
    positive_definite_prime
    and
    positive_definite_by_orthogonal
    return the same number of results
    """
    expected = positve_definite_by_orthogonal(genus,prime,e)
    rep = positve_definite_prime(genus, prime, e)
    if len(rep) != len(expected):
        print(len(rep), expected)
        return False, rep, expected
    return True, rep, expected


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

def k3_pq(genus,p, q):
    r"""
    Returns all conjugacy classes of isometries of order `pq` in `genus`
    such that the kernel of the `pq`-part is of signature `(2,n)`.

    INPUT:

    - ``genus`` -- a genus of signature `(3,k)`.
    - ``p,q`` -- a prime numbers
    """
    n = 4
    rk = genus.rank()
    weights = [1,euler_phi(p),euler_phi(q), euler_phi(p*q)]
    floor = [1,0,0,1]
    ceiling = [rk // d for d in weights]

    P = IntegerListsLex(max_sum=rk,length=4,floor=floor, ceiling=ceiling)
    P = [[i[k]*weights[k] for k in range(n)] for i in P]
    P = [i for i in P if sum(i)==rk]
    acts = []
    for ranks in P:
        ranks_E = [ranks[k]//weights[k] for k in range(n)]
        signatures = [[ranks_E[0]-1]]
        signatures += [[ranks_E[k]]*(max(1,weights[k]//2)) for k in range(1,4)]
        signatures[-1][0] -= 1
        acts += pq_action(genus,p,q, ranks, signatures)
    return acts

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

    # check that rationally (A + B)_p = C_p
    # since the determinants match the excess is sufficient
    if AB.local_symbols(p).excess() != C.local_symbols(p).excess():
        return False

    Ap = A.local_symbols(p).symbol_tuple_list()
    Bp = B.local_symbols(p).symbol_tuple_list()
    a1 = sum(s[1] for s in Ap if s[0] == 1)
    a2 = sum(s[1] for s in Ap if s[0] >= 2)
    b1 = sum(s[1] for s in Bp if s[0] == 1)
    b2 = sum(s[1] for s in Bp if s[0] >= 2)
    Cp_level = C.level().valuation(p)
    a_max = sum(s[1] for s in Ap if s[0]==Cp_level + 1)
    b_max = sum(s[1] for s in Bp if s[0]==Cp_level + 1)

    ker_min = max(g - a1, g - b1, a_max, b_max)
    ker_max = min(a2 + a1//2, b2 + b1//2, g)
    # compare kernel dimensions of the glue
    # for the existence of a glue map
    if ker_max < ker_min:
        return False

    # the ranks of the jordan components
    # of maximum possible scales must match
    if a_max != b_max:
        return False

    ABp = AB.local_symbols(p).symbol_tuple_list()
    Cp = C.local_symbols(p).symbol_tuple_list()
    # check that the discriminant group
    # of Cp is a subquotient of that of ABp
    if any(r1 > r2 for r1,r2 in
            ((sum(s[1] for s in Cp[k:]),sum(s[1] for s in ABp if s[0]>= Cp[k][0]))
             for k in range(0,len(Cp)))
          ):
        print('subquotient')
        False

    # if the glue is exactly p^l D_i,
    # then the lower scale parts of A and B are those of C
    # as they are not glued.
    # we have to be careful at p=2 because of sign walking etc.
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
        for i in range(Cp_level+1):
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

    # the unimodular component of A+B injects
    # into the unimodular component of C
    s1 = A.local_symbols(p).symbol(0)
    s2 = B.local_symbols(p).symbol(0)
    s3 = C.local_symbols(p).symbol(0)
    if s1[1]+s2[1] == s3[1]:
        s1 = Genus_Symbol_p_adic_ring(p,[s1])
        s2 = Genus_Symbol_p_adic_ring(p,[s2])
        s3 = Genus_Symbol_p_adic_ring(p,[s3])
        if s1.direct_sum(s2) != s3:
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
    # we have passed all the checks I could think of
    return True


def prime_power_actions(genus, p, ranks, signatures, k3_unobstructed=True,verbose=False,magm=my_magma):
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
        raise StopIteration
    if len(ranks) == 1:
        # there is nothing to glue
        for M in genus.representatives():
            fM = M ^ 0
            M = IntegralLattice(M)
            GM = Oq_equiv(M, fM, 1)
            yield [M, fM, GM]
            raise StopIteration
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
        for M, fM in all_lattice_with_isometry(p^e, M_rk_E, M_det, max_level,
                                           signatures=M_signature,
                                           min_scale=min_scale,min_norm=min_norm):
            cm = sage.structure.element.get_coercion_model()
            cm.reset_cache()
            M = IntegralLattice(M)
            pos, neg = M.signature_pair()
            # In the coinvariant lattice in NS there are no roots
            if all([k3_unobstructed, e>0, pos==0, neg>0,
                   M.maximum()==-2]):
                continue
            if verbose:
                countM += 1
                print("M,fM: %s"%countM)
                print(M.genus())
            GM = Oq_equiv(M, fM, p^e)
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
                                                 verbose=verbose,magm=magm):
                        N_max_glue = G.scale()*N.dual_lattice() & (1/p * N)
                        N_max_glue = N.discriminant_group().submodule(N_max_glue.gens())
                        ext = extensions(M, fM, N, fN, GM, GN,
                                         glue_order, p,
                                         DM_max_glue=M_max_glue,
                                         DN_max_glue=N_max_glue,
                                         magm=magm, target_genus=genus)
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
        ptype.append([(p,e0),L.genus(),C.genus(),R.genus(),ZZ(glue)])
        L = R
    ptype.append([(p,0),L.genus(),L.genus(),Genus(matrix([])),ZZ(1)])
    ptype.reverse()
    return ptype

def pnq_actions(q, ptype,k3_unobstructed=True):
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

    # find all actions p^e q of types p^e
    # that glue to genus
    for (M, fM, GM) in splitpq(Cgenus,p,e,q):
        # base case
        cm = sage.structure.element.get_coercion_model()
        cm.reset_cache()
        if e == 0:
            yield (M, fM, GM)
            continue
        # recurse
        for (N, fN, GN) in pnq_actions(q,ptype[:-1],k3_unobstructed=k3_unobstructed):
            ext = extensions(M, fM, N, fN, GM, GN, pglue, p)
            for ex in ext:
                if ex[0].genus() == genus:
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
        assert pos % 2 == neg % 2 == 0
        if pos == 0:
            return [rkE]*(degE//2)
        sig = [rkE - 1]
        if degE > 2:
            sig += [rkE]*(degE//2 - 1)
        return sig


def split_sig(sig_pair,p,e,q):
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
        # the order must be p^e q and not just p^e
        floor[1] = 1
    if pos % 2 == 1:
        # there must be an invariant vector of pos. square
        if e != 0:
            raise StopIteration
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

def splitpq(genus,p,e,q, k3_unobstructed=True):
    r"""
    """
    assert p.is_prime() and q.is_prime()
    primes = [sym.prime() for sym in genus.local_symbols() if sym.prime() != q]

    for rks,sigs in split_sig(genus.signature_pair(),p,e,q):
        degsE = euler_phi(p**e), euler_phi(q*p**e)
        rksE = [rks[0]//degsE[0],rks[1]//degsE[1]]

        if rks[1] == 0:
            # there is nothing to glue
            for M, fM in all_lattice_with_isometry(p^e,rksE[0],genus.det(),
                                                   genus.level(),signatures=sigs[0],
                                                   min_scale=genus.scale()):
                M = IntegralLattice(M)
                if M.genus() == genus:
                    pos,neg = M.signature_pair()
                    if k3_unobstructed and e > 0 and pos == 0 and neg > 0:
                        if M.maximum() == -2:
                            continue
                    GM = Oq_equiv(M,fM, p^e)
                    yield M, fM, GM
            continue

        if rks[0] == 0:
            # there is nothing to glue
            for R, fR in all_lattice_with_isometry(q*p^e,rksE[1],genus.det(),
                                                   genus.level(),signatures=sigs[1],
                                                   min_scale=genus.scale()):
                R = IntegralLattice(R)
                if R.genus() == genus:
                    pos,neg = M.signature_pair()
                    if k3_unobstructed and e > 0 and pos == 0 and neg > 0:
                        if M.maximum() == -2:
                            continue
                    GR = Oq_equiv(R,fR, q*p^e)
                    yield R, fR, GR
            continue

        min_scale = genus.scale()
        min_norm = genus.norm()
        max_level = q * genus.level()

        m = degsE[0]*min(rksE)
        max_det = q^m*genus.det()

        for M_det in max_det.divisors():
            for M, fM in all_lattice_with_isometry(q*p^e, rksE[1],
                                           M_det, max_level,
                                           signatures=sigs[1],
                                           min_scale=min_scale,
                                           min_norm=min_norm):
                M = IntegralLattice(M)
                pos,neg = M.signature_pair()
                if k3_unobstructed and e > 0 and pos == 0 and neg > 0:
                    if M.maximum() == -2:
                        continue
                GM = Oq_equiv(M, fM, q*p^e)
                if rks[1]==genus.rank():
                    # there is nothing to glue here
                    if genus==M.genus():
                        yield [M, fM, GM]
                    continue
                DM = M.discriminant_group()
                M_max_glue_group = (M.span((fM^(p^e) - fM**0).inverse())
                                    & M.dual_lattice())
                M_max_glue_group = DM.submodule(M_max_glue_group.gens())
                for glue_order in M_max_glue_group.cardinality().divisors():
                    R_det = genus.det() * glue_order^2 / M.determinant()
                    try:
                        R_det = ZZ(R_det)
                    except TypeError:
                        # the glue group is too small
                        continue
                    if not glue_order.divides(R_det):
                        continue
                    for R,fR in all_lattice_with_isometry(p^e, rksE[0],
                                           R_det, max_level,
                                           signatures=sigs[0],
                                           min_scale=min_scale):
                        R = IntegralLattice(R)
                        if not is_admissible(M.genus(),R.genus(),genus,q):
                            return False
                        # now MR is a serious candidate
                        GR = Oq_equiv(R,fR,p**e)
                        ext = extensions(M, fM, R, fR, GM, GR, glue_order, q,magm=my_magma target_genus=genus)
                        for ex in ext:
                            if ex[0].genus() == genus:
                                yield ex


def Oq_equiv(M, fM, order):
    r"""
    Image of the equivariant part of the orthogonal group.
    """
    if order==1 or order==2:
        return M.image_in_Oq()
    elif prod(M.signature_pair()) == 0:
        #
        G = M.orthogonal_group()
        f = G(fM)
        G = [G(g) for g in G.gap().Centralizer(f.gap()).GeneratorsOfGroup()]
        Oq = M.discriminant_group().orthogonal_group()
        return Oq.subgroup([Oq(g) for g in G])
    else:
        OM = M.orthogonal_group([-fM**0, fM])
        if True:
            # we might be undercounting here
            OGM = M.image_in_Oq()
            G = OGM.gap().Centralizer(OGM.subgroup(OM.gens()))
            return OGM.subgroup(G.GeneratorsOfGroup())
        else:
            # TODO: the stabilizer of the action of fN is still missing!
            # thus we are overcounting!!!
            OGM = M.discriminant_group().orthogonal_group()
            return OGM.subgroup(OM.gens())

def extensions(M, fM, N , fN, GM, GN, glue_order, p,
               DM_max_glue=None,DN_max_glue=None,
               magm=None,target_genus=None):
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
    for glue, stab in DM.all_primitive_prime_equiv(DN, DM_max_glue, DN_max_glue,
                                                   GM, GN, fDM, fDN, glue_val,H10=H10,H20=H20,
                                                   magma=magm,target_genus=target_genus):
        ext = glue.W().overlattice([g.lift() for g in glue.gens()])
        assert ext*f == ext
        ext_Oq = ext.discriminant_group().orthogonal_group()

        stab = ext_Oq.subgroup([ext_Oq(g) for g in stab.gens()])
        results.append([ext, f, stab])
    return results


def pq_action(genus, p, q, ranks, signatures):
    r"""
    Returns all conjugacy classes of isometries of order `pq` in `genus`.

    EXAMPLES::

        sage: target = all_genera_by_det((3,7),2^2)[0]
        sage: p = 3
        sage: ranks = [2, 2, 6]
        sage: signatures = [[0], [1], [0,1,1]]
        sage: prime_power_actions(target, p, ranks, signatures)

    """
    H = genus

    # first we treat the part of order p, q

    rank_1 = ranks[0]
    rank_p = ranks[1]
    rank_q = ranks[2]
    rank_pq = ranks[3]

    # the number of negative squares in the real embeddings of K < E
    sig_1 = signatures[0]
    sig_p = signatures[1]
    sig_q = signatures[2]
    sig_pq = signatures[3]

    rk = H.rank()
    assert sum(ranks) == rk

    e = len(ranks) - 1
    deg_E_pq = euler_phi(p*q)
    deg_E_p = euler_phi(p)
    deg_E_q = euler_phi(q)
    rk_E_pq = rank_pq//deg_E_pq
    rk_E_p = rank_p//deg_E_p
    rk_E_q = rank_q//deg_E_q
    assert rk_E_pq > 0

    max_glue_pq_q = p**(min(rk_E_pq, rk_E_q) * min(deg_E_pq, deg_E_q))
    max_glue_pq_p = q**(min(rk_E_pq, rk_E_p) * min(deg_E_pq, deg_E_p))
    max_glue_pq = max_glue_pq_p * max_glue_pq_q
    max_glue_q_1 = q**(min(rk_E_q, rank_1) * min(deg_E_q, 1))

    primes = [sym.prime() for sym in H.local_symbols()]
    min_scale = H.scale()
    max_level_pq = H.level()*q**min(1,rank_p, rank_pq)*p**min(1,rank_q, rank_pq)
    max_level_q = H.level()*p**min(rank_pq, rank_q, 1)*q**min(rank_1, rank_q, 1)
    min_det_pq = min_scale**rank_pq
    max_det_pq = H.determinant() * max_glue_pq

    results = []
    Cpq_dets = [det_pq for det_pq in max_det_pq.divisors() if min_det_pq.divides(det_pq)]
    Cpq_poss = []
    for det_pq in Cpq_dets:
        Cpq_poss += all_lattice_with_isometry(p*q, rk_E_pq, det_pq, max_level_pq, signatures=sig_pq, min_scale=min_scale)


    for Cpq, fpq in Cpq_poss:
        Cpq = IntegralLattice(Cpq)
        Gpq = Oq_equiv(Cpq, fpq, p*q)
        if rk==rank_pq:
            # there is nothing to do here
            if H==Cpq.genus():
                results.append([Cpq, fpq, Oq_equiv(Cpq,fpq,p*q)])
            continue

        # find the possible Cq
        min_det_q = min_scale**rank_q
        max_det_q = max_glue_pq_q * max_glue_q_1 * H.determinant() / (Cpq.determinant()/(max_glue_pq_p*max_glue_pq_q)).numerator()
        try:
            max_det_q = ZZ(max_det_q)
        except TypeError:
            # the glue group is too small
            continue
        Cq_dets = [det_q for det_q in max_det_q.divisors() if min_det_q.divides(det_q)]
        if rank_q == 0:
            Cq_dets = [1]
        Cq_poss = []
        for det_q in Cq_dets:
            Cq_poss += all_lattice_with_isometry(q, rk_E_q, det_q, max_level_q, signatures=sig_q, min_scale=min_scale)
        for Cq, fq in Cq_poss:
            if Cq.rank() == 0:
                ext = [[Cpq, fpq, Gpq]]
            else:
                Cq = IntegralLattice(Cq)
                Gq = Oq_equiv(Cq,fq,q)
                ext = []
                for glue_order_q in max_glue_pq_q.divisors():
                    ext += extensions(Cpq, fpq, Cq, fq, Gpq, Gq, glue_order_q, p)
            for CpqCq, fpqq, Gpqq in ext:
                # now we need to find the appropriate orthogonal complement
                # CpC1 and glue it along q
                max_glue_R = max_glue_pq_p * max_glue_q_1
                for glue_order_R in max_glue_R.divisors():
                    R_det = H.det() * glue_order_R^2 / CpqCq.determinant()
                    try:
                        R_det = ZZ(R_det)
                    except TypeError:
                        # the glue group is too small
                        continue
                    if not glue_order_R.divides(R_det):
                        continue
                    sig_R = vector(H.signature_pair()) - vector(CpqCq.signature_pair())
                    R_max_level = H.level() * q**max(min([rk_E_pq, rk_E_p, 1]),min([rk_E_q, rank_1, 1]))
                    for R in all_genera_by_det(sig_R, R_det, max_scale=R_max_level):
                        if not min_scale.divides(R.scale()):
                            continue
                        # check that (M + R)_s = G_q for all s != p
                        CpqCqR = CpqCq.genus().direct_sum(R)
                        if not all(CpqCqR.local_symbols(s)==H.local_symbols(s) for s in primes if s != q):
                            continue
                        # double check that the determinants match
                        if not (CpqCqR.det() * H.det()).is_square():
                            continue
                        # check that rationally (M + R)_q = G_q
                        # since the determinants match the excess is sufficient
                        if CpqCqR.local_symbols(q).excess() != H.local_symbols(q).excess():
                            continue
                        # now CpqCqR is a serious candidate
                        # compute the possible order p actions on R
                        for N in prime_power_actions(R, p, [rank_1,rank_p,], [sig_1,sig_p]):
                            N, fN, GN = N
                            ext = extensions(CpqCq, fpqq, N, fN, Gpqq, GN, glue_order_R, q)
                            ext = [ex for ex in ext if ex[0].genus()==H]
                            results += ext
    return results

