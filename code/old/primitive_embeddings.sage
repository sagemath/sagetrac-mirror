from sage.interfaces.gap import get_gap_memory_pool_size
memory_gap = get_gap_memory_pool_size()
set_gap_memory_pool_size(4*memory_gap)
libgap.eval("SetRecursionTrapInterval(10000000)");



def primitive_embeddings(A, B):
    r"""
    Compute all primitive embeddings of A into the genus of B
    up to equvalence.

    INPUT:

    - ``A, B`` -- Integral lattices

    OUTPUT:

    - a list of lattices
    """
    # notation
    Bm = B.twist(-1)
    DA = A.discriminant_group()
    DBm = Bm.discriminant_group()
    GA = DA.orthogonal_group(tuple())
    GB = DBm.orthogonal_group()
    # compute the primitive extensions modulo OB
    extensions = []
    tmp, i1, i2 = DA.all_primitive(DBm, DA, DBm, GA, GB)
    for glue in tmp:
        ext = glue[0]
        ext = ext.W().overlattice(ext.V().matrix())
        stab = glue[1]
        stab = ext.discriminant_group().orthogonal_group().subgroup(stab)
        extensions.append([ext,stab])
    # compute the orthogonal complements
    sigK = vector(B.signature_pair()) - vector(A.signature_pair())
    Hlist = []
    for ext in extensions:
        E = ext[0]
        stab = ext[1]
        try:
            genusK = E.discriminant_group().twist(-1).genus(sigK)
        except ValueError:
            continue
        for K in genusK.representatives():
            K = IntegralLattice(K)
            GK = K.image_in_Oq()
            GE = stab
            Hlist += primitive_extension_unimodular(E, K, GE, GK)
    # extract B and the embedding of A into B
    results = []
    for H in Hlist:
        Aembed = i1.matrix()*H[1].matrix()
        H = H[0].orthogonal_complement(i2.matrix()*H[1].matrix())
        Aembed = H.sublattice(Aembed)

        Kembed = H.orthogonal_complement(Aembed)
        AK = H.span(Aembed.gens()+Kembed.gens())

        A1 = IntegralLattice(Aembed.gram_matrix())
        K1 = IntegralLattice(Kembed.gram_matrix())
        AK1, (e1,e2) = A1.direct_sum(K1,return_embeddings=True)
        H = AK1.overlattice([AK.coordinate_vector(g) for g in H.gens()])

        assert H.genus() == B.genus()
        assert A1.genus() == Aembed.genus()
        results.append([H,H.sublattice(e1.image().matrix()),H.sublattice(e2.image().matrix())])
    return results


def primitive_extension_unimodular(A, B, GA, GB):
    r"""

    WARNING:

    The size of O(DB) must not be crazy big.

    EXAMPLES::

        sage: A = IntegralLattice("D4").twist(2)
        sage: B = A.twist(-1)
        sage: GA = A.image_in_Oq()
        sage: GB = B.discriminant_group().orthogonal_group().subgroup([])
        sage: primitive_extension_unimodular(A, B, GA, GB)
    """

    DA = A.discriminant_group().normal_form()
    tmp = B.discriminant_group().twist(-1).normal_form().gens()
    DB = B.discriminant_group().submodule_with_gens(tmp)
    n = len(DA.gens())

    domA = GA.domain()
    domB = GB.domain()
    genA = [domA(g).gap() for g in DA.gens()]
    genB = [domB(g).gap() for g in DB.gens()]

    phi = domA.gap().GroupHomomorphismByImages(domB.gap(),genA,genB)

    OqB = DB.orthogonal_group()

    GAphi = [phi.InducedAutomorphism(g) for g in GA.gap().GeneratorsOfGroup()]
    GAphi = OqB.gap().Subgroup(GAphi)
    representatives = OqB.gap().DoubleCosetRepsAndSizes(GB,GAphi)
    extensions = []
    AB, (iA, iB) = A.direct_sum(B, return_embeddings=True)
    for g in representatives:
        g = OqB(g[0], False)
        gens = [DA.gen(k).lift()*iA.matrix().change_ring(QQ) + (DB.gen(k)*g).lift()*iB.matrix().change_ring(QQ) for k in range(n)]
        C = AB.overlattice(gens)
        assert C.determinant().abs() == 1
        extensions.append([C,iA,iB])
    return extensions


