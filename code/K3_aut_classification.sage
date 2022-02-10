from sage.interfaces.gap import get_gap_memory_pool_size
memory_gap = get_gap_memory_pool_size()
set_gap_memory_pool_size(9048*memory_gap)
libgap.eval("SetRecursionTrapInterval(10000000)");

from sage.quadratic_forms.genera.genus import genera
load("databases/hashimotoslattices.sage")
load("symplectic.sage")
load("hermitian.sage")
load("lattice_with_isometry.sage")
load("prime_power.sage")

transcendental_values = [n for n in range(2,100) if euler_phi(n) <=20]
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
    if B.rank() == 0:
        # purely non-symplectic
        return [K3SurfaceAutGrp(A, A.orthogonal_group([]),a,a.change_ring(ZZ).multiplicative_order())]    
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
    print(B)
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
            k3s.append(K3SurfaceAutGrp(H, H.orthogonal_group(G0), g, n))
        except ValueError:
            pass
    return k3s


class K3SurfaceAutGrp(object):
    r"""

    Input:

    - H an even unimodular lattice of signature (3,19)
    - G0 a group of symplectic automorphisms
    that is H_G0 is negative definite, root free and the action on
    the discriminant group is trivial.
    - g an isometry of H such that g respects the
    orthogonal decomposition H = H_G0 + H^G0
    and g is non-symplectic
    """
    def __init__(self, H, G0, g, n, name=None, check=True):
        r"""
        """
        assert H.det()==-1
        assert H.is_even()
        assert H.signature_pair() == (3,19)
        self.name=name
        self.H = H
        self.G0 = G0
        self.g = g
        self.G = H.orthogonal_group([h.matrix() for h in G0.gens()] + [g])
        self.n = n
        if check and self.NS_coinvariant().max() == -2:
            raise ValueError('This automorphism is obstructed')

    @cached_method
    def transcendental_lattice(self):
        r"""
        """
        R = PolynomialRing(QQ,"x")
        x = R.gen()
        cn = cyclotomic_polynomial(self.n,'x')
        g = self.g
        S = self.symplectic_invariant_lattice()
        K = self.H.kernel_sublattice(cn(g))
        return self.H.sublattice((K & S).gens())

    @cached_method
    def NS(self):
        r"""
        """
        T = self.transcendental_lattice()
        return self.H.orthogonal_complement(T)

    @cached_method
    def symplectic_invariant_lattice(self):
        r"""
        """
        V = self.H
        for g in self.G0.gens():
            g = g.matrix()
            V = V & V.kernel_sublattice(g-1)
        return self.H.sublattice(V.gens())

    @cached_method
    def invariant_lattice(self):
        r"""
        """
        HG = self.symplectic_invariant_lattice()
        return HG.kernel_sublattice(self.g - 1)

    @cached_method
    def NS_coinvariant(self):
        H = self.H
        NS = self.NS()
        inv = self.invariant_lattice()
        return H.sublattice(
            H.orthogonal_complement(inv) & NS)

    @cached_method
    def symplectic_co_invariant_lattice(self):
        r"""
        """
        return self.H.orthogonal_complement(self.symplectic_invariant_lattice())

    def str(self):
        r"""
        A string representation of this K3 group
        for storage.
        """
        s =  "[" + str(self.n) + ", ["
        s += str(self.H.basis_matrix().list()) + ","
        s += str(self.H.inner_product_matrix().list()) + "],"
        s += str([g.matrix().list() for g in self.G0.gens()]) + ","
        s += str(self.g.list()) + "]"
        return s

def aut_from_str(s, check=False):
    r"""

    """
    L = sage_eval(s)
    if type(L[0])==str:
        name = L[0]
        L = L[1:]
    else:
        name=None
    n = L[0]
    deg = ZZ(len(L[1][1]))
    assert deg.is_square()
    deg = deg.sqrt()
    basis = matrix(QQ,22,deg,L[1][0])
    gram = matrix(QQ,deg,deg,L[1][1])
    H = IntegralLattice(gram, basis)

    G0_gens = [matrix(QQ,deg,deg,h) for h in L[2]]
    G0 = H.orthogonal_group(G0_gens)
    g = matrix(QQ,deg,deg,L[3])

    return K3SurfaceAutGrp(H, G0, g, n, check=check,name=name)


#######################################
def get_ptypes(order,fi=None,number=None,rkT=None,orderG0=None):
    assert len(order.prime_factors())==1
    p = order.prime_factors()[0]
    if fi is not None:
        pass
    elif rkT is not None:
        fi = open('results/order%s_no%s_rkT%s.txt'%(order,number,rkT))
    elif number is not None:
        fi = open('results/order%s_no%s.txt'%(order,number))
    else:
        fi = open('results/order%s.txt' % order)
    ptypes = []
    for s in fi:
        if s[:8]=="complete":
            continue
        k3 = aut_from_str(s)
        if rkT is not None and k3.transcendental_lattice().rank()!=rkT:
            continue
        H = k3.symplectic_invariant_lattice()
        if rkT is not None and rkT!=k3.transcendental_lattice().rank():
            continue
        if orderG0 is not None and orderG0 != k3.G0.order():
            continue
        g = k3.g
        t = ptype(H,g,p)
        if not t in ptypes:
            ptypes.append(t)
    return ptypes

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

def classify_ord_pe_parallel(L, p, e, file_name,rw,verbose=False):
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
        for res in k3_prime_power(fix.genus(), p, e,verbose=verbose):
            for Aa in res[1]:
                ipm,basis, a, Oa = Aa
                A = IntegralLattice(ipm,basis)
                Oa = A.discriminant_group().orthogonal_group().subgroup(Oa)
                actsg = K3SurfaceAut(A, cofix, a, Oa)
                for aut in actsg:
                    classifi.append(aut)
                    s = aut.str()
                    result = open(file_name,"a")
                    result.write(s+ "\n")
                    result.close()
                if len(actsg)==0:
                    not_realized.append([fix,a])
    return classifi, not_realized


def classify_purely_ns_pn(p,e,file_name, log_file, rw="w",verbose=True, rkT=None):
    classifi = []
    result = open(file_name,rw)
    result.close()
    not_realized = []
    g = genera((3,19),1,1,even=true)[0]
    for rkT0 in range(2,22):
      if rkT is not None and rkT0!=rkT:
          continue
      for A, a, Oa in k3_prime_power(g, p, e,verbose=verbose, rkT=rkT0):
          aut = K3SurfaceAutGrp(A,A.orthogonal_group([]),a,p^e)
          classifi.append(aut)
          s = aut.str()
          result = open(file_name,"a")
          result.write(s+ "\n")
          result.close()
      log = open(log_file, 'a')
      log.write("completed purely non-symplectic order %s^%s, rkT = %s\n"%(p,e,rkT))
      log.close()
    log = open(log_file, 'a')
    log.write("completed purely non-symplectic order %s^%s \n"%(p,e))
    log.close()
    result = open(file_name,"a")
    result.write("complete \n")
    result.close()
    return classifi

def classify_purely_ns_peq(p, q ,file_r, file_aw, log_file,aw="w",verbose=2,rkT=None):
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
        k3 = aut_from_str(str_pe)
        if rkT is not None and k3.transcendental_lattice().rank()!=rkT:
            continue
        L = k3.symplectic_invariant_lattice()
        if rkT is not None and rkT != k3.transcendental_lattice().rank():
            continue
        if k3.NS_coinvariant().maximum() == -2:
            # the isometry is obstructed
            assert False
        f = k3.g
        assert p == k3.n.prime_factors()[0]
        typ = ptype(L,f,p)
        if typ in types:
            continue
        types.append(typ)
        g = genera((3,19),1,1,even=true)[0]
        for A, a, Oa in pnq_actions(q, typ, verbose=verbose):
            n = a.change_ring(ZZ).multiplicative_order()
            aut = K3SurfaceAutGrp(A, A.orthogonal_group([]), a, n)
            classifi.append(aut)
            s = aut.str()
            result = open(file_aw,"a")
            result.write(s+ "\n")
            result.close()
    result = open(file_aw,'a')
    result.write("complete \n")
    result.close()
    log = open(log_file, 'a')
    log.write("completed purely non-symplectic order %s^%s*%s \n"%(p,e,q))
    log.close()
    return classifi


def classify_purely_ns_6_par(file_r, file_aw, log_file,aw="w",verbose=2,rkT=None, k3_unobstructed=True):
    p = 3
    q = 2
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
        k3 = aut_from_str(str_pe)
        L = k3.symplectic_invariant_lattice()
        if rkT is not None and rkT != k3.transcendental_lattice().rank():
            continue
        if k3.NS_coinvariant().maximum() == -2:
            # the isometry is obstructed
            assert False
        f = k3.g
        assert p == k3.n.prime_factors()[0]
        typ = ptype(L,f,p)
        if typ in types:
            continue
        types.append(typ)
    g = genera((3,19),1,1,even=true)[0]
    inputs = [] 
    for typ in types:
        g1 = typ[0][2]
        g3 = typ[1][2]
        s1 = [a for a in split_sig(g1.signature_pair(),p,0,q)]
        s3 = [a for a in split_sig(g3.signature_pair(),p,1,q)]
        for s11 in s1:
            for s33 in s3:
                inputs.append((q,typ,((s11,),(s33,)),k3_unobstructed))
    for input,output in pnq_actions_pure(inputs):
        log = open(log_file,"a")
        log.write(str(input))
        log.write(str(output))
        log.close()
        for A,a in output: 
            n = a.change_ring(ZZ).multiplicative_order()
            aut = K3SurfaceAutGrp(A, A.orthogonal_group([]), a, n,check=k3_unobstructed)
            classifi.append(aut)
            s = aut.str()
            result = open(file_aw,"a")
            result.write(s+ "\n")
            result.close()
    result = open(file_aw,'a')
    result.write("complete \n")
    result.close()
    return classifi

def classify_ord_peq(p, q, file_r, file_aw,aw='w',verbose=2, rkT=None):
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
        k3 = aut_from_str(str_pe)
        L = k3.symplectic_invariant_lattice()
        if k3.NS_coinvariant().maximum() == -2:
            # the isometry is obstructed
            continue
        if rkT is not None and k3.transcendental_lattice().rank()!=rkT:
            continue
        f = k3.g
        assert p == k3.n.prime_factors()[0]
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


###########################################################
def classify_rank4(L, file_name):
    result = open(file_name,"w")
    result.close()
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
        Hacts += k3_prime_power(fix.genus(), 3, 1)
        Hacts += k3_pq(fix.genus(),2,3)
        for Aa in Hacts:
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
    return classifi, not_realized, Hacts
###########################################################


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
    #g = grp_restriction(L,k3.g)
    #t = ptype(L,g,2)
    #ge =pnq_actions(3,t)
    #t
    #for g in ge:
        #acts12.append(g)
