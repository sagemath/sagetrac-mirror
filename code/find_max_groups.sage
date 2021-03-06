from sage.interfaces.gap import get_gap_memory_pool_size
from sage.libs.gap.util import GAPError

memory_gap = get_gap_memory_pool_size()
set_gap_memory_pool_size(9048*memory_gap)
libgap.eval("SetRecursionTrapInterval(10000000)");

folder = '/home/simon/Dropbox/Math/sage/lattices/results/'
files = ['order%s.txt'%k for k in range(2,17)+[18,20,24,28]]

K3s = [] # the list of all non-symplectic actions
for filename in files:
    fi = open(folder + filename,'r')
    for s in fi:
        if s!='\n':
            K3s.append(aut_from_str(s))
    print(len(K3s, filename))

def subs(k3):
    G = k3.G.gap()
    G.IsFinite()
    G0 = k3.G.subgroup(k3.G0).gap()
    subgroups = []
    for S in G.ConjugacyClassesSubgroups():
        S = S.Representative()
        S0 = G0.Intersection(S)
        n = S.Size()/S0.Size()
        i
        if n > 1 and S0.Size() > 1:
            subgroups.append((get_id(S0),get_id(S)))
    return subgroups

def get_id(G):
    try:
        idG = G.IdGroup().sage()
    except GAPError:
        idG = (ZZ(G.Size()),ZZ(id(G)))
    return tuple(idG)

def k3_id(k3):
    G = k3.G.gap()
    G.IsFinite()
    G0 = k3.G.subgroup(k3.G0).gap()
    return (get_id(G0),get_id(G))

def get_max(K3s):
    K3s = copy(K3s)
    K3s.sort(key=lambda x: x.G.order())
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

#print('loaded groups')
#known, maximal= get_max(K3s)
