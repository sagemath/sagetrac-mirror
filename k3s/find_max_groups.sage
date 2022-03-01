from sage.interfaces.gap import get_gap_memory_pool_size
from sage.libs.gap.util import GAPError

memory_gap = get_gap_memory_pool_size()
set_gap_memory_pool_size(16384*memory_gap)
libgap.eval("SetRecursionTrapInterval(10000000)");


def subs(k3):
    G = k3.group().gap()
    G.IsFinite()
    G0 = k3.group().subgroup(k3.symplectic_subgroup()).gap()
    subgroups = []
    for S in G.ConjugacyClassesSubgroups():
        S = S.Representative()
        S0 = G0.Intersection(S)
        n = S.Size()/S0.Size()
        if n > 1:
            subgroups.append((get_id(S0),get_id(S)))
    return subgroups

def get_id(G):
    try:
        idG = G.IdGroup().sage()
    except GAPError:
        idG = (ZZ(G.Size()),"no id")
    return tuple(idG)

def k3_id(k3):
    G = k3.group().gap()
    G0 = k3.symplectic_subgroup().gap()
    G.IsFinite()
    G0.IsFinite()
    k3id = (get_id(G0),get_id(G))
    k3._k3id = k3id
    return k3id

def get_multiplicity(K3s):
    known = []
    for k3 in K3s:
        k3id = k3_id(k3)
        known.append(k3id)
    return known

def get_max_pairs(K3s):
    K3s = copy(K3s)
    K3s.sort(key=lambda x: x.group().order())
    known = []
    maximal = []
    while len(K3s)>0:
        k3 = K3s.pop()
        k3id = k3_id(k3)
        print(len(K3s),k3id)
        if not k3id in known:
            known.append(k3id)
            maximal.append(k3id)
            for s in subs(k3):
                if not s in known:
                    known.append(s)
    return known, maximal

def get_max(K3s):
    K3s = copy(K3s)
    K3s.sort(key=lambda x: x.group().order())
    known = []
    maximal = []
    while len(K3s)>0:
        k3 = K3s.pop()
        k3id = k3_id(k3)[1]
        print(len(K3s),k3id)
        if not k3id in known:
            known.append(k3id)
            maximal.append(k3id)
            for s0,s1 in subs(k3):
                if not s1 in known:
                    known.append(s1)
    return known, maximal


#print('loaded groups')
#known, maximal= get_max(K3s)

def abelian_invs(k3):
    assert k3.group().is_abelian()
    a = k3.group().gap().AbelianInvariants().sage()
    a = AbelianGroup(a).elementary_divisors()
    a0 = k3.symplectic_subgroup().gap().AbelianInvariants().sage()
    a0 = AbelianGroup(a0).elementary_divisors()
    return (a0,a)


def issub(d1,d2):
    if len(d1)>len(d2):
        return false
    return all(ZZ(d1[-k]).divides(d2[-k]) for k in range(len(d1)))



def is_enriques_invol(H,g):
    from sage.quadratic_forms.genera.genus import genera
    G = genera((1,9),2^10,2,True)[0]
    assert G.norm()==4
    assert g!=1 and g^2==1
    g = kernel_sublattice(H,g-1).genus()
    return G == g


def get_enriques_id(k3,g,return_sub=False):
    g = k3.group()(g)
    C = k3.group().gap().Centralizer(g.gap())
    Cg = C.Subgroup([g])
    S = C.FactorGroup(Cg)
    S0 = C.Intersection(k3.symplectic_subgroup())
    S0s = C.Subgroup([g] + [h for h in S0.GeneratorsOfGroup()])
    S0s = S0s.FactorGroup(S0s.Subgroup([g]))
    assert S0.IdGroup()==S0s.IdGroup()
    if return_sub:
        c = S.NaturalHomomorphism()
        return S.Subgroup([c.Image(g) for g in S0.GeneratorsOfGroup()]),S
    return S0.IdGroup().sage(),S.IdGroup().sage()

def get_enriques_ids(k3,return_sub=False):
    ids = []
    C = k3.group().conjugacy_classes_representatives()
    H = k3.H
    for c in C:
        if c.order()==2 and is_enriques_invol(H,c.matrix()):
            ids.append(get_enriques_id(k3,c,return_sub))
    return ids

def get_max_from_ids(ids):
    ids = copy(ids)
    ids.sort()
    ids.reverse()
    known = []
    max = []
    for i in ids:
        G = libgap.SmallGroup(i)
        if not i in known:
            max.append(i)
        S = G.ConjugacyClassesSubgroups()
        for s in S:
            s = s.Representative()
            j = s.IdGroup().sage()
            known.append(tuple(j))
    return max

hashimoto_notation =dict([
 [(1, 1),""],
 [(2, 1),"C_2"],
 [(3, 1),"C_3"],
 [(4, 1),"C_4"],
 [(4, 2),"C_2^2"],
 [(6, 1),"D_6"],
 [(8, 3),"D_8"],
 [(8, 5),"C_2^3"],
 [(10, 1),r"D_{10}"],
 [(12, 3),"A_4"],
 [(12, 4),r"D_{12}"],
 [(16, 11),r"C_2 \times D_8"],
 [(16, 14),"C_2^4"],
 [(18, 4),"A_{3,3}"],
 [(20, 3),"Hol(C_5)"],
 [(21, 1),"C_7:C_3"],
 [(32, 49),"Q_8 * Q_8"],
 [(36, 9), "3^2C_4"],
 [(36, 10),"S_{3,3}"],
 [(48, 29),r"T_{48}"],
 [(48, 48),r"C_2 \times S_4"],
 [(48, 50), "2^4C_3"],
 [(60, 5), "A_5"],
 [(64, 138),r"\Gamma_{25}a_1"],
 [(72, 40),r"N_{72}"],
 [(72, 41),"M_9"],
 [(72, 43),"A_{4,3}"],
 [(96, 227),"2^4D_6"],
 [(120, 34),"S_5"],
 [(168, 42),"L_2(7)"],
 [(192, 955),r"H_{195}"],
 [(192, 1023),r"4^2A_4"],
 [(192, 1493),r"T_{192}"],
 [(288, 1026),r"A_{4,4}"],
 [(360, 118),r"A_6"],
 [(384, 18135),r"F_{384}"],
 [(960, 11357),  r"M_{20}"]
])


def has_same_global_type(k1,k2):
    assert k1.n == k2.n
    if ptype_big(k1) != ptype_big(k2):
        return false
    R.<x> = QQ[]
    for d in [d for d in divisors(k1.n) if d!=1 and d!=k1.n]:
        c = cyclotomic_polynomial(d,x)
        L1 = - k1.H.kernel_sublattice(c(k1.g)).gram_matrix()
        L2 = - k2.H.kernel_sublattice(c(k2.g)).gram_matrix()
        q = QuadraticForm
        if not q(L1).is_globally_equivalent_to(q(L2)):
            return false
    return true

def ptype_small(k):
    n = k.n
    g = k.g
    L = k.H
    typ = []
    for d in divisors(n):
        typ.append(L.kernel_sublattice(g^d-1).genus())
    return typ


def ptype_big(k):
    n = k.n
    g = k.g
    L = k.H
    typ = []
    glues = []
    for d in divisors(n):
        Ld = L.kernel_sublattice(g^d-1)
        Cd = L.kernel_sublattice(cyclotomic_polynomial(d)(g))
        typ.append((Ld.genus(),Cd.genus()))
        if Cd.rank()==0:
            s = 0
        else:
            q = Cd.q()
            s = q.submodule(Ld.proj(Cd).rows()).normal_form().gram_matrix_quadratic()
        glues.append(s)
    return (typ,glues)


def make_names(k3s):
  def sortkey(k):
    no = get_hashimoto_number(k.symplectic_invariant_lattice())
    n = k.transcendental_value()
    T = k.transcendental_lattice()
    F = k.invariant_lattice()
    d = k.transcendental_lattice().rank()/euler_phi(n)
    return (no,n,d,T.det(),T.discriminant_group().invariants(),F.rank(),F.det(),F.discriminant_group().invariants())
  k3s.sort(key=sortkey)
  no0 = None
  n0 = None
  d0 = None
  for k3 in k3s:
    no,n,d = sortkey(k3)[:3]
    if (no,n,d)!=(no0,n0,d0):
      no0,n0,d0 = (no,n,d)
      k = 1
    s = f"{no}.{n}.{d}.{k}"
    k = k+1
    k3._id = s


def make_names_save(folder_write):
  for n in [n for n in range(67,67) if euler_phi(n)<=20 and n!=60]:
    print(n)
    k3s = K3SurfaceAutGrp_from_database(n,True)
    make_names(k3s)
    if n < 10:
      ord = "0"+str(n)
    else:
      ord = str(n)
    fi = open(folder_write+f"/purely_ns/order{ord}.txt","w")
    for k in k3s:
      fi.write(str_new(k))
      fi.write("\n")
    fi.close()

  for n in [n for n in range(14,67) if euler_phi(n)<=12]:
    k3s = K3SurfaceAutGrp_from_database(n,False)
    print(n)
    make_names(k3s)
    if n < 10:
      ord = "0"+str(n)
    else:
      ord = str(n)
    fi = open(folder_write+f"/order{ord}.txt","w")
    for k in k3s:
      fi.write(str_new(k))
      fi.write("\n")
    fi.close()




def str_new(self):
    s =  "["
    s += f"['id', '{self.id()}'],"
    s += f"['transcendental_value', {self.transcendental_value()}],"
    s += f"['connected_components', {self.connected_components()}],"
    s += f"['basis', {self.L().basis_matrix().list()}],"
    s += f"['inner_product_matrix', {self.L().inner_product_matrix().list()}],"
    s += f"['symplectic', {[g.matrix().list() for g in self.symplectic_subgroup().gens()]}],"
    s += f"['distinguished_generator', {self._g.list()}],"
    s += f"['properties',{str([(k,self._prop[k]) for k in []])}]"
    s += "]"
    return s

def adapt_name(k):
  s = k.id().split(".")
  s[2] = eval(s[2]) - 1
  st = f"{s[0]}.{s[1]}.{s[2]}.{s[3]}"
  k._id = st

def adapt_names(folder_write):
  for n in [n for n in range(2,67) if euler_phi(n)<=20 and n!=60]:
    print(n)
    k3s = K3SurfaceAutGrp_from_database(n,True)
    if n < 10:
      ord = "0"+str(n)
    else:
      ord = str(n)
    fi = open(folder_write+f"/purely_ns/order{ord}.txt","w")
    for k in k3s:
      #fpp = tuple(fixed_points_power(k,i) for i in divisors(n)[:-1])
      #k.properties()["fixed points powers"] = fpp
      #k.properties()["fixed points"] = fpp[0]
      fi.write(k.str())
      fi.write("\n")
    fi.close()

  for n in [n for n in range(67,67) if euler_phi(n)<=12]:
    k3s = K3SurfaceAutGrp_from_database(n,False)
    print(n)
    if n < 10:
      ord = "0"+str(n)
    else:
      ord = str(n)
    fi = open(folder_write+f"/order{ord}.txt","w")
    for k in k3s:
      fi.write(k.str())
      fi.write("\n")
    fi.close()
