from sage.rings.finite_rings.finite_field_prime_modn import *
from sage.matroids.advanced import *
from sage.sets.set import Set
from sage.all import *
msnccons=[  [ [[],[]], [[],[]], [[],[]] ],[] ]
#R=None
GF2=GF(2)


def all_1def_pcodes(Mpcode,nc):
    varset = ncinstance_vars(nc)
    for v in varset:
        new_pmap = {0:v,1:v,2:v}
        if is_pcode(Mpcode._M,new_pmap,nc):
             Mpcode._C.append(new_pmap)
    return Mpcode

def all_2def_pcodes(Mpcode,nc):
    varset = ncinstance_vars(nc)
    pmaps2 = [{0: 1, 1: 0, 2: 1}, {0: 0, 1: 1, 2: 0}, {0: 1, 1: 1, 2: 0},
              {0: 0, 1: 0, 2: 1}, {0: 0, 1: 0, 2: 1}, {0: 1, 1: 1, 2: 0},
              {0: 0, 1: 1, 2: 0}, {0: 1, 1: 0, 2: 1}, {0: 1, 1: 0, 2: 0},
              {0: 0, 1: 1, 2: 1}]
    for sub2 in Subsets(varset,2):
        sub2_list = list(sub2)
        for pmap in pmaps2:
            new_pmap = {key:sub2_list[pmap[key]] for key in pmap.keys()}
            if is_pcode(Mpcode._M,new_pmap,nc):
                 Mpcode._C.append(new_pmap)
    return Mpcode

def all_3def_pcodes(Mpcode,nc):
    n3pcode=0
    varset = ncinstance_vars(nc)
    pmaps3=[{0: 0, 1: 1, 2: 2},
            {0: 0, 1: 2, 2: 1},
            {0: 1, 1: 0, 2: 2},
            {0: 1, 1: 2, 2: 0},
            {0: 2, 1: 0, 2: 1},
            {0: 2, 1: 1, 2: 0}]
    for sub3 in Subsets(varset,3):
        sub3_list = list(sub3)
        for pmap in pmaps3:
            new_pmap = {key:sub3_list[pmap[key]] for key in pmap.keys()}
            if is_pcode(Mpcode._M,new_pmap,nc):
                 Mpcode._C.append(new_pmap)
                 n3pcode+=1
    return Mpcode

def pmap_def2():
    # return normalized pmaps defining two variables
    pmaps=[]
    p2 = OrderedSetPartitions({0,1,2})[7:-1]
    for p in p2:
        e1 = list(p[0])
        e2 = list(p[1])
        v1=0
        v2=1
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
    return pmaps

def pmap_def3():
    pmaps=[]
    p1 = OrderedSetPartitions({0,1,2})[0:6]
    for p in p1:
        e1 = list(p[0])[0]
        e2 = list(p[1])[0]
        e3 = list(p[2])[0]
        pmap = {e1:0,e2:1,e3:2}
        pmaps.append(pmap)
    return pmaps

class MatroidCodes():
    def __init__(self,M=None,pmaps=None):
        self._M=M
        self._C=[]
        self._parent=None
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
    maxrank = N-(len(ncinstance_vars(nc))-len(nc[1]))
    M3 = all2(3,maxrank)
    for key in M3.keys():
        if sum(list(key)) is not 3:
            M3.pop(key)
    allvars=ncinstance_vars(nc)
    size_3_matroidcodes = {}
    for key in M3:
        size_3_matroidcodes[key] = []
        for M in M3[key]:
            Mpcode = MatroidCodes(M,[])
            Mpcode = all_1def_pcodes(Mpcode,nc)
            Mpcode = all_2def_pcodes(Mpcode,nc)
            Mpcode = all_3def_pcodes(Mpcode,nc)
            if len(Mpcode._C) > 0:
                size_3_matroidcodes[key].append(Mpcode)
    return size_3_matroidcodes

def init_enum_3subset(nc,N,subset):
    """
    returns matroids on 3 elements and all the p-codes they form
    """
    varset = list(subset)
    p1 = OrderedSetPartitions({0,1,2})[0:6]
    p2 = OrderedSetPartitions({0,1,2})[7:-1]
    p3 = OrderedSetPartitions({0,1,2})[-1]
    pmaps=[]
    for p in p1:
        e1 = list(p[0])[0]
        e2 = list(p[1])[0]
        e3 = list(p[2])[0]
        pmap = {e1:varset[0],e2:varset[1],e3:varset[2]}
        pmaps.append(pmap)
    for p in p2:
        e1 = list(p[0])
        e2 = list(p[1])
        for vs in Subsets(subset,2):
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
    pmaps.extend([{0:varset[0],1:varset[0],2:varset[0]},{0:varset[1],1:varset[1],2:varset[1]},{0:varset[2],1:varset[2],2:varset[2]}])
    return pmaps

def ncinstance_vars(nc):
    varset=Set([])
    for i in nc:
        for j in nc[0]:
            v=Set(j[1])
            varset=varset.union(v)
    varset.union(Set(nc[1]))
    return varset

def extend_mcodes(Mcodes1):
    Mcodes2 = {}
    if len(Mcodes1.keys()) < 1:
        return Mcodes2
    n1 = sum(list(Mcodes1.keys()[0]))
    n2 = n1+1
    Mcodes1.keys[0]
    for key in Mcodes1:
        #ext = next2([Mcode._M for M in ]
        Mcodes2[(list(key)[0],n2-list(key)[0])] = []
        #for Mcode in Mcodes[key]:

def Mcode_candidates(Mcode,maxgnd,nvars):
    candidates = []
    for mcode in Mcode._C:
        # loop over all subsets of defined variables and delete them
        imcode = invert(mcode)
        defvars = set(imcode.keys())
        for delset in Subsets(defvars):
            imcode_copy = imcode.copy()
            for v in delset:
                imcode_copy.pop(v)
            # get candidate pmap
            cand = invert(imcode_copy)
            print cand
            if nvars-len(cand.values()) <=  (maxgnd-len(Mcode._M.groundset())):
                candidates.append(cand)
    return candidates

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
        if M.rank(rhs_el) != M.rank(lhs_el):
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




def extendpmaps(M,M_codes,M_child):
    """
    Return a maximal list of p-codes corresponding to a given matroid
    """

    return

def naive_candidates(M_ext,list_of_pmaps):

    cand = []

def butterfly():
 return [[ [[3,4],[1,3,4]],[[5,6],[1,5,6]],[[7],[4,5,7]],[[8,9],[7,8,9]],[[2],[2,3,8]],[[1],[1,6,9]] ] ,[1,2]]




