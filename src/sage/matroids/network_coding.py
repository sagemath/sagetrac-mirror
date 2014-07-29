from sage.rings.finite_rings.finite_field_prime_modn import *
from sage.matroids.advanced import *
from sage.sets.set import Set
from sage.all import *
msnccons=[  [ [[],[]], [[],[]], [[],[]] ],[] ]
R=None
GF2=GF(2)

class MatroidCodes():
    def __init__(self,M=None,pmaps=None):
        self._M=M
        self._C=[]
        if M is not None and pmaps is not None:
            self.add(pmaps)

    def add(self, pmaps):
        for pmap in pmaps:
            addpcode(self._C, pmap)




def addpcode(list_of_pcodes,pcode):
    """
    adds p-code to a list of p-codes if there is no p-code in ``list_of_pcodes`` that is extension
    of ``pcode``
    """
    for cpcode in list_of_pcodes:
        if is_pmap_deletion(pcode,cpcode):
            return
        elif is_pmap_deletion(cpcode,pcode):
            list_of_pcodes.remove(cpcode)
            list_of_pcodes.append(pcode)
            return
    # brand new pcode
    list_of_pcodes.append(pcode)
    return


def is_pmap_deletion(dict1,dict2):
    """
    checks whether $dict1\subseteq dict2$
    """
    if len(dict1) > len(dict2):
        return False
    idict1=invert(dict1)
    idict2=invert(dict2)
    # check variable set containment
    if Set(idict1.keys()).issubset(Set(idict2.keys())) == True:
        for key in idict1.keys():
            if idict1[key] != idict2[key]: #non-matching variable definitions
                return False
    else:
        return False
    return True


def invert(m):
    """
    Inverts the given network to matroid or matroid to network mapping
    """
    if len(m) > 0:
        if isinstance(m[m.keys()[0]],type(Set([]))):
            # this is network to matroid mapping
            mi={}
            for v in m.keys():
                Ev=m[v]
                for e in Ev:
                    mi[e]=v
        else:
            # this is matroid to network mapping
            mi={v:Set([]) for v in Set(m.values())}
            for e in m.keys():
                for v in Set([m[e]]):
                    mi[v]=mi[v].union(Set([e]))
        return mi
    else:
        return {}

def add2(R,M):
    for N in R:
        if M.is_field_isomorphic(N):
            return
    R.append(M)

def extensions2(M):
    R=[]
    x=len(M)
    for N in M.linear_extensions(len(M)):
        add2(R,N)
    return R

def next2(R):
    R2=[]
    for M in R:
        for N in extensions2(M):
            add2(R2,N)
    return R2

def all2(N, R=10):
    MM={}
    for n in xrange(N+1):
        MM[n,0]=[BinaryMatroid(identity_matrix(GF2,n))]
        for r in xrange(min(n,R+1)):
            MM[r,n-r]=next2(MM[r,n-r-1])
        #print 'n=',n
        #print [len(MM[r,n-r]) for r in xrange(min(n,R)+1)]
    return MM

def ncmatenum(nc,r,N):
    n=r
    MM={}
    MM[n]=[BinaryMatroid(identity_matrix(GF(2),n)),{}]


def init_enum(nc,N):
    """
    returns matroids on 3 elements and all the p-codes they form
    """
    maxrank = N-len(ncinstance_vars(nc))-len(nc[1])
    M3 = all2(3,R=maxrank)
    for key in M3.keys():
        if sum(list(key)) is not 3:
            M3.pop(key)
    p1 = OrderedSetPartitions({0,1,2})[0:6]
    p2 = OrderedSetPartitions({0,1,2})[7:-1]
    p3 = OrderedSetPartitions({0,1,2})[-1]
    pmaps=[]
    for p in p1:
        e1 = list(p[0])[0]
        e2 = list(p[1])[0]
        e3 = list(p[2])[0]
        pmap = {e1:1,e2:2,e3:3}
        pmaps.append(pmap)
    for p in p2:
        e1 = list(p[0])
        e2 = list(p[1])
        for vs in Subsets({1,2,3},2):
            v1=list(vs)[0]
            v2=list(vs)[1]
            pmap1={}
            for e in e1:
                pmap1[e]=v1
            for e in e2:
                pmap1[e]=v2
            pmap2={}
            for e in e1:
                pmap2[e]=v2
            for e in e2:
                pmap2[e]=v1
            pmaps.append(pmap1)
            pmaps.append(pmap2)
    pmaps.extend([{0:1,1:1,2:1},{0:2,1:2,2:2},{0:3,1:3,2:3}])
    size_3_matroidcodes = {}
    for key in M3:
        size_3_matroidcodes[key] = []
        for M in M3[key]:
            Mpcodes=[]
            for pmap in pmaps:
                if is_pcode(M,pmap,nc):
                    Mpmaps.append(pmap)
            size_3_matroidcodes[key].append(MatroidCodes(M,Mpmaps))
    return size_3_matroidcodes

def ncinstance_vars(nc):
    varset=Set([])
    for i in nc:
        for j in nc[0]:
            v=Set(j[1])
            varset=varset.union(v)
    varset.union(Set(nc[1]))
    return varset

def applicable_cons(def_vars,nc):
    cons=[]
    for j in nc[0]:
         if all([k in def_vars for k in j[1]]):
             cons.append(j)
    return cons



def is_pcode(M,pmap,ncinstance):
    ipmap = invert(pmap)
    def_vars = ipmap.keys()
    nodecons=applicable_cons(def_vars,ncinstance)
    for con in nodecons:
        lhs_el = []
        for v in con[0]:
            lhs_el.extend(list(ipmap[v]))
        rhs_el=[]
        for v in con[1]:
            rhs_el.extend(list(ipmap[v]))
        if M.rank(rhs_el) != M.el(lhs_el):
            return False
    def_src = [v for v in def_vars if v in ncinstance[1]]
    if len(def_src) > 0:
        src_sum=0
        for s in def_src:
            src_sum += M.rank(ipmap[s])
        src_joint = M.rank(def_src)
        if src_sum != src_joint:
            return False
    return True


def see(N,mdict):
    """
    Return a dictionary that maps matroids in mdict to their respective extensions
    The keys are tuples ``(r,index)`` where ``r`` is the rank and ``index`` is the
    index of parent in the list it was in ``mdict``
    """


def extendpmaps(M,M_codes,M_child):
    """
    Return a maximal list of p-codes corresponding to a given matroid
    """

    return

def naive_candidates(M_ext,list_of_pmaps):

    cand = []

def butterfly():
 return [[ [[3,4],[1,2,3,4]], [[5,7],[3,5,7]], [[6,7],[4,6,7]] ],[1,2]]




