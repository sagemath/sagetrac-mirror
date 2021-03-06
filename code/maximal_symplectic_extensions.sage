from sage.interfaces.gap import get_gap_memory_pool_size
load("k3_classification/symplectic.sage")
memory_gap = get_gap_memory_pool_size()
set_gap_memory_pool_size(8*memory_gap)
libgap.eval("SetRecursionTrapInterval(10000000)");


def find_group(fix, cofix, g, ids=False):
    r"""
    """
    gm = g.matrix()
    K = (gm - gm^0).change_ring(ZZ).kernel().basis()
    invariant = fix.sublattice([fix(k) for k in K])
    c = cyclotomic_polynomial(g.order(),x)
    transcendental = fix.sublattice(c(gm).change_ring(ZZ).kernel().basis())
    qfix = fix.discriminant_group().normal_form()
    qcofix = cofix.discriminant_group().normal_form()
    gqfix = qfix.orthogonal_group()(g)
    gqcofix = transport(qfix, qcofix, gqfix)

    O = cofix.orthogonal_group()
    Oq = qcofix.orthogonal_group()
    f = O.hom([Oq(h) for h in O.gens()], check=False)
    glift = f.lift(gqcofix)
    K = f.kernel()
    G = O.subgroup(K.gens() + (glift,))
    try:
        idG = G.gap().IdGroup()
    except ValueError:
        idG = "unknown"
    s = [
        fix,
        K.order(),
        K.gap().IdGroup(),
        transcendental.gram_matrix(),
        invariant.gram_matrix(),
        g.order(),
        G.order(),
        idG]
        #G.structure_description()]
    print(s)
    #pretty_print(s)
    if ids:
        return s
    else:
        return K,G,glift
#    print("n=", g.order(), "polarization:", invariant.gram_matrix())
#    print("transcendental lattice")
#    print(transcendental.gram_matrix())
#    print(G.gap().IdGroup())
#    print(G.structure_description())


def transport(A, B, g):
    r"""
    Transport g to O(B).
    """
    assert A.gram_matrix_quadratic() == B.gram_matrix_quadratic()
    g = g.matrix()   # wrt smith form gens
    g = A._to_smith() * g * A._to_gens()      # wrt normal form gens
    g = B._to_gens() * g * B._to_smith() # wrt smith gens
    return B.orthogonal_group()(g)

R.<x> = QQ[]
def classify_range(r):
    grps = []
    for k in r:
        fix = fixed[k]
        cofix = cofixed[k]
        if fix.rank() == 3:
            # The case that the fixed lattice is positive definite
            # is particularly simple.
            print(fix.gram_matrix())
            Of = fix.orthogonal_group()
            cand = set()
            not_maximal = set()
            for conj in Of.conjugacy_classes():
                g = conj.representative()
                chi = g.matrix().charpoly()
                # there is an invariant ample class
                if not (x-1).divides(chi):
                    continue
                # the transcendental lattice has rank 2
                # and there the minimal polynomal is irreducible and
                # for us (non-symplectic) different from (x-1)
                if ((x-1)^2).divides(chi):
                    continue
                assert (g.matrix().minpoly()//(x-1)).is_irreducible()
                cand.add(conj)
                # we want to get rid of actions from higher orders.
                for k in g.order().divisors():
                    if k == 1:
                        continue
                    not_maximal.add(Of.conjugacy_class(g^k))
            maximal = cand.difference(not_maximal)
            maximal = [g.representative() for g in maximal]
            maximal.sort(key=lambda x: x.order())
            #grps+=[(fix,maximal)]
            for g in maximal:
                Og = Oq_equiv(fix,g.matrix(),g.order())
                grps.append(K3SurfaceAut(fix, cofix.twist(-1), g.matrix(),Og ))
    return grps
