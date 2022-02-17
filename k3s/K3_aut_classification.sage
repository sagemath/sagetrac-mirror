from sage.interfaces.gap import get_gap_memory_pool_size
memory_gap = get_gap_memory_pool_size()
set_gap_memory_pool_size(9048*memory_gap)
libgap.eval("SetRecursionTrapInterval(10000000)");

from sage.env import SAGE_ROOT

dir = SAGE_ROOT + "/k3s/"

from sage.quadratic_forms.genera.genus import genera
load(dir+"K3SurfaceAutGrp.py")
load(dir+"symplectic.sage")
load(dir+"hermitian.sage")
load(dir+"lattice_with_isometry.sage")
load(dir+"prime_power.sage")

transcendental_values = [n for n in range(2,67) if euler_phi(n) <=20]
#############################################################
# the K3 side of things
#############################################################

all_genera_by_det = genera
from sys import stdout
def print_remove(k):
    stdout.write('\r %s'%k)
    stdout.flush()
    #stdout.write('\n')



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

def K3SurfaceAut(A, B, a, Oa):
    r"""
    - A -- lattice; fixed lattice
    - B -- lattice; cofixed lattice
    - a -- matrix; an isometry of A
    - Oa -- image of the stabilizer of a in Oq
    """
    a = A.orthogonal_group([a]).gen(0)
    n = a.order()

    AB = A.direct_sum(B)
    assert AB.signature_pair() == (3,19), [A, B]
    assert A.determinant().abs() == B.determinant().abs()
    i = matrix.identity(A.degree()+B.degree())
    iA = i[:A.degree(),:]
    iB = i[A.degree():,:]

    qA = A.discriminant_group().normal_form()
    OqA = qA.orthogonal_group()

    gens_qB = B.twist(-1).discriminant_group().normal_form().gens()
    qB = B.discriminant_group().submodule_with_gens(gens_qB)
    tmp = qA.gram_matrix_quadratic() +qB.gram_matrix_quadratic()
    assert tmp.denominator() == 1
    assert ZZ(gcd(tmp.diagonal())) % 2 == 0
    OqB = qB.orthogonal_group()
    a_qA = OqA(a)


    A1 = OqA.domain()
    A2 = OqB.domain()
    gens1 = [A1(g).gap() for g in qA.gens()]
    gens2 = [A2(g).gap() for g in qB.gens()]

    # the glue map coming from the normal forms
    phi0 = A1.gap().GroupHomomorphismByImages(A2.gap(), gens1, gens2)

    # notation
    OB = B.orthogonal_group()
    pi = OB.gap().GroupHomomorphismByImagesNC(OqB.gap(),[OqB(h).gap() for h in OB.gens()])
    G0 = pi.Kernel().GeneratorsOfGroup()
    G0 = [matrix.block_diagonal([matrix.identity(A.degree()),OB(h).matrix()]) for h in G0]
    OB_in_OqB = pi.Image(pi.Source())
    Oa_in_OqB = OqB.gap().Subgroup([phi0.InducedAutomorphism(g.gap()) for g in Oa.gens()])

    reps = OqB.gap().DoubleCosetRepsAndSizes(OB_in_OqB, Oa_in_OqB)
    k3s = []
    VA = qA.V().ambient_vector_space()
    VB = qB.V().ambient_vector_space()
    for h in reps:
        h = h[0]
        phih = phi0*h
        bbar = phih.InducedAutomorphism(a_qA.gap())
        h = OqB(h)
        if not OqB(bbar).gap() in OB_in_OqB:
            continue
        gens = [VA(qA.gen(k).lift())*iA + VB((qB.gen(k)*h).lift())*iB for k in range(len(qA.gens()))]

        # create H and g
        H = AB.overlattice(gens)
        b = OB(pi.PreImagesRepresentative(bbar)).matrix()
        g = matrix.block_diagonal([a.matrix(), b])
        try:
            k3s.append(K3SurfaceAutGrp(H, G0, g, n))
        except ValueError:
            pass
    return k3s


#######################################

def classify_ord_p_next(cofix, ptypes, file_name, rw, verbose=0):
    classifi = []
    result = open(file_name, rw)
    result.close()
    not_realized = []
    cofix = cofix.twist(-1)
    k = 0
    for ptype in ptypes:
        print("type number %s"%k)
        k = k+1
        for A, a, Oa in next_prime_power(ptype, verbose=verbose):
            actsg = K3SurfaceAut(A, cofix, a, Oa)
            for aut in actsg:
                classifi.append(aut)
                s = aut.str()
                result = open(file_name,"a")
                result.write(s+ "\n")
                result.close()
            if len(actsg)==0:
                not_realized.append([A,a])
    result = open(file_name,"a")
    result.write("complete \n")
    result.close()
    return classifi, not_realized

def invariant_lattice(L,G):
    F = L
    for g in G.gens():
        F = F.kernel_sublattice(g.matrix()-1)
    return F

def write_order_2class_script():
  fi = open("/home/simon/sage/commands_order2_impure","w")
  for no in range(0,44):
    for rkT in range(2,14):
      s = """./sage -c 'load("k3s/K3_aut_classification.sage"); """
      s += """classify_ord_pe(L=[%s], p=2, e=1,file_name="k3s/results/order2/order2_no%s_rkT%s",rw= "w", verbose=True,rkT=%s)'"""%(no,no,rkT,rkT)
      s+= "\n"
      fi.write(s)
  fi.close()
end

def classify_ord_pe(L, p, e, file_name,rw,verbose=False,rkT=None):
    classifi = []
    result = open(file_name,rw)
    result.close()
    not_realized = []
    for k in L:
        print("Invariant Lattice number %s"%k)
        fix = fixed[k]
        print(fix.gram_matrix())
        print(" ")
        cofix = cofixed[k].twist(-1)
        for Aa in k3_prime_power(fix.genus(), p, e,verbose=verbose,rkT=rkT):
            A, a, Oa = Aa
            actsg = K3SurfaceAut(A, cofix, a, Oa)
            for aut in actsg:
                classifi.append(aut)
                s = aut.str()
                result = open(file_name,"a")
                result.write(s+ "\n")
                result.close()
            if len(actsg)==0:
                not_realized.append([fix,a])
    result = open(file_name,"a")
    result.write("complete \n")
    result.close()
    return classifi, not_realized

def classify_purely_ns_pn(p, e, file_name, rw="w", verbose=True, rkT=None):
    classifi = []
    peactions = open(file_r,'r')
    result = open(file_aw,aw)
    result.close()
    not_realized = []
    k = 0
    types = []
    for str_pe in peactions:
        k+=1
        print(k)
        if str_pe[:8]=="complete":
            continue
        k3 = K3SurfaceAutGrp_from_str(str_pe)
        L = k3.symplectic_invariant_lattice()
        if k3.NS_coinvariant().maximum() == -2:
            # the isometry is obstructed
            continue
        f = k3.distinguished_generator()
        assert p == k3.transcendental_value().prime_factors()[0]
        typ = ptype(L,f,p)
        if typ in types:
            continue
        types.append(typ)
        g = genera((3,19),1,1,even=true)[0]
        for A, a, Oa in pnq_actions(q, typ, verbose=verbose):
            n = a.change_ring(ZZ).multiplicative_order()
            aut = K3SurfaceAutGrp(A, [], a, n)
            classifi.append(aut)
            s = aut.str()
            result = open(file_aw,"a")
            result.write(s+ "\n")
            result.close()
    result = open(file_aw,'a')
    result.write("complete \n")
    result.close()
    return classifi


def classify_purely_ns_nextp(p, file_r, file_aw, aw="w",verbose=2):
    classifi = []
    peactions = open(file_r, 'r')
    result = open(file_aw, aw)
    result.close()
    not_realized = []
    k = 0
    types = []
    for str_pe in peactions:
        k+=1
        print(k)
        if str_pe[:8]=="complete":
            continue
        k3 = K3SurfaceAutGrp_from_str(str_pe)
        L = k3.symplectic_invariant_lattice()
        if k3.NS_coinvariant().maximum() == -2:
            # the isometry is obstructed
            continue
        f = k3.distinguished_generator()
        assert p == k3.transcendental_value().prime_factors()[0]
        typ = ptype(L,f,p)
        if typ in types:
            continue
        types.append(typ)
        g = genera((3,19),1,1,even=true)[0]
        for A, a, Oa in next_prime_power(typ, verbose=verbose):
            n = a.change_ring(ZZ).multiplicative_order()
            aut = K3SurfaceAutGrp(A, [], a, n)
            classifi.append(aut)
            s = aut.str()
            result = open(file_aw,"a")
            result.write(s+ "\n")
            result.close()
    result = open(file_aw,'a')
    result.write("complete \n")
    result.close()
    return classifi

def classify_purely_ns_peq(p, q ,file_r, file_aw, aw="w",verbose=2):
    classifi = []
    peactions = open(file_r,'r')
    result = open(file_aw,aw)
    result.close()
    not_realized = []
    k = 0
    types = []
    for str_pe in peactions:
        k+=1
        print(k)
        if str_pe[:8]=="complete":
            continue
        k3 = K3SurfaceAutGrp_from_str(str_pe)
        L = k3.symplectic_invariant_lattice()
        if k3.NS_coinvariant().maximum() == -2:
            # the isometry is obstructed
            continue
        f = k3.distinguished_generator()
        assert p == k3.transcendental_value().prime_factors()[0]
        typ = ptype(L,f,p)
        if typ in types:
            continue
        types.append(typ)
        g = genera((3,19),1,1,even=true)[0]
        for A, a, Oa in pnq_actions(q, typ, verbose=verbose):
            n = a.change_ring(ZZ).multiplicative_order()
            aut = K3SurfaceAutGrp(A, [], a, n)
            classifi.append(aut)
            s = aut.str()
            result = open(file_aw,"a")
            result.write(s+ "\n")
            result.close()
    result = open(file_aw,'a')
    result.write("complete \n")
    result.close()
    return classifi

def classify(p, q, file_read, file_write, pure, verbose=2):
    assert is_prime(p)
    assert is_prime(q)
    if p == q and not pure:
      classify_mixed_nextp(p, file_read, file_write)
    if p == q and pure:
      classify_purely_ns_nextp(p, file_read, file_write)
    if p != q and not pure:
      classify_ord_peq(p, q, file_read, file_write,"w",verbose=verbose)
    if p != q and pure:
      classify_purely_ns_peq(p, q, file_read, file_write,verbose=verbose)

def classify_ord_peq(p, q, file_r, file_aw,aw='w',verbose=2):
    r"""
    file_r contains the p^e actions
    """
    classifi = []
    peactions = open(file_r,'r')
    result = open(file_aw,aw)
    result.close()
    not_realized = []
    k = 0
    types = []
    for str_pe in peactions:
        if str_pe[:8]=="complete":
            continue
        k+=1
        print(k)
        k3 = K3SurfaceAutGrp_from_str(str_pe)
        L = k3.symplectic_invariant_lattice()
        if k3.NS_coinvariant().maximum() == -2:
            # the isometry is obstructed
            continue
        f = k3.distinguished_generator()
        assert p == k3.transcendental_value().prime_factors()[0]
        typ = ptype(L,f,p)
        if typ in types:
            continue
        types.append(typ)
        cofix = IntegralLattice( k3.symplectic_co_invariant_lattice().gram_matrix())
        for A, a, Oa in pnq_actions(q,typ,verbose=verbose):
            order_a = a.change_ring(ZZ).multiplicative_order()
            print('Found an isometry of order: %s' % order_a)
            assert q.divides(order_a)
            actsg = K3SurfaceAut(A, cofix, a, Oa)
            for aut in actsg:
                classifi.append(aut)
                sres = aut.str()
                result = open(file_aw,'a')
                result.write(sres+ "\n")
                result.close()
            if len(actsg)==0:
                not_realized.append([A,a])
    result = open(file_aw,'a')
    result.write("complete \n")
    result.close()
    return classifi, not_realized

def classify_mixed_nextp(p, file_r, file_aw, aw='w',verbose=2):
    r"""
    file_r contains the p^e actions
    """
    classifi = []
    peactions = open(file_r,'r')
    result = open(file_aw,aw)
    result.close()
    not_realized = []
    k = 0
    types = []
    for str_pe in peactions:
        if str_pe[:8]=="complete":
            continue
        k+=1
        print(k)
        k3 = K3SurfaceAutGrp_from_str(str_pe)
        L = k3.symplectic_invariant_lattice()
        if k3.NS_coinvariant().maximum() == -2:
            # the isometry is obstructed
            continue
        f = k3.distinguished_generator()
        assert p == k3.transcendental_value().prime_factors()[0]
        typ = ptype(L,f,p)
        if typ in types:
            continue
        types.append(typ)
        cofix = IntegralLattice( k3.symplectic_co_invariant_lattice().gram_matrix())
        for A, a, Oa in next_prime_power(typ, verbose=verbose):
            order_a = a.change_ring(ZZ).multiplicative_order()
            print('Found an isometry of order: %s' % order_a)
            actsg = K3SurfaceAut(A, cofix, a, Oa)
            for aut in actsg:
                classifi.append(aut)
                sres = aut.str()
                result = open(file_aw,'a')
                result.write(sres+ "\n")
                result.close()
            if len(actsg)==0:
                not_realized.append([A,a])
    result = open(file_aw,'a')
    result.write("complete \n")
    result.close()
    return classifi, not_realized

#########################
#Utility
############################

def grp_restriction(M,g):
    g = M.orthogonal_group([g]).gen(0)
    g = matrix([M.coordinate_vector(x*g) for x in M.gens()])
    L = M.overlattice(M.ambient_module().gens())
    Orthogonal = L.orthogonal_complement(M)
    B = M.basis_matrix().stack(Orthogonal.basis_matrix())
    identity = matrix.identity(Orthogonal.rank())
    g = matrix.block_diagonal([g, identity])
    g = B.inverse()*g*B
    return g

#acts12 = []
#for k3 in k3s:
    #L = k3.symplectic_invariant_lattice()
    #g = grp_restriction(L,k3.distinguished_generator())
    #t = ptype(L,g,2)
    #ge =pnq_actions(3,t)
    #t
    #for g in ge:
        #acts12.append(g)


def spread_types_to_files(path_read, path, transcendental_value_read):
    os.mkdir(path)
    fi = open(path_read, 'r')
    L = get_ptypes(transcendental_value_read, path_read, returnk3s=True)
    n = len(L)
    for k in range(n):
        s = L[k].str()
        fi = open(path + "/" + str(k), "w")
        fi.write(s)
        fi.close()

def get_ptypes(order, file_path, rkT=None, returnk3s=False):
    assert len(order.prime_factors())==1
    p = order.prime_factors()[0]
    ptypes = []
    k3s = []
    fi = open(file_path,"r")
    for s in fi:
        if s[:8] == "complete":
            continue
        k3 = K3SurfaceAutGrp_from_str(s)
        assert k3.transcendental_value() == order
        if rkT is not None and k3.transcendental_lattice().rank()!=rkT:
            continue
        H = k3.symplectic_invariant_lattice()
        g = k3.distinguished_generator()
        t = ptype(H,g,p)
        if not t in ptypes:
            ptypes.append(t)
            k3s.append(k3)
    fi.close()
    if returnk3s:
      return k3s
    return ptypes


