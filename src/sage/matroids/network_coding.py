#from sage.rings.finite_rings.finite_field_prime_modn import *
from sage.matroids.advanced import *
from sage.sets.set import Set
R=None
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
        MM[n,0]=[BinaryMatroid(identity_matrix(GF(2),n))]
        for r in xrange(min(n,R+1)):
            MM[r,n-r]=next2(MM[r,n-r-1])
        print [len(MM[r,n-r]) for r in xrange(min(n,R)+1)]
    return MM

msnccons=[  [ [[],[]], [[],[]], [[],[]] ],[] ]

def extendpmaps(M,M_codes,M_e):
    """
    Return a maximal list of p-codes corresponding to a given matroid
    """
    
    return



def addpcode(list_of_pcodes,pcode):
    """
    adds p-code to a list of p-codes if there is no p-code in ``list_of_pcodes`` that is extension
    of ``pcode``
    """
    for cpcode in list_of_pcodes:
        for key in pcode.keys():
            if key in cpcode.keys():
                 if cpcode[key]!=pcode[key]:
                     return
            else: 
                return
    list_of_pcodes.append(pcode)

       
    
def is_pmap_deletion(dict1,dict2):
    if len(dict1) > len(dict2):
        return False 
    idict1=invert(dict1)
    idict2=invert(dict2)
    # check variable set containmentadd
    if Set(idict1.keys()).issubset(Set(idict2.keys())) == True:
        for key in idict1.keys():
            if idict1[key] != idict2[key]: #non-matching variable definitions
                return False  
    else:
        return False
    return True
    
    


def naive_candidates(M_ext,list_of_pmaps):
    cand = []
    
    
def invert(m):
    """
    Inverts the given network to matroid or matroid to network mapping
    """
    if len(m) > 0:
        if isinstance(m[1],type(Set([]))):
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
                

        
                
    
def butterfly():
 return [[ [[3,4],[1,2,3,4]], [[5,7],[3,5,7]], [[6,7],[4,6,7]] ],[1,2]]
 

 

