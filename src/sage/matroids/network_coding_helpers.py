from sage.matroids.advanced import BinaryMatroid, TernaryMatroid, QuaternaryMatroid
from sage.all import * 


def appcons(cons):
    r"""
    Return a dictionary of network constraints applicable to each subset
    of indices appearing in ``cons``.
    
    INPUT:
    
    - A list of lists specifying network constraints specified as 
    ``[list1,list2]`` where ``list1`` and ``list2`` are subsets of 
    random variable indices that are forced to have equal rank (entropy)
    
    OUTPUT:
    
    - A dictionary with subsets of indices appearing in ``cons`` as keys
    and associated set of applicable constraint indices as values
    
    .. NOTE::

            This method does NOT do any checks.
    
    EXAMPLES::
    
        
    """
    allvars=set([])
    for k in cons.keys():
        allvars=allvars|cons[k][0]|cons[k][1]
    apcns={}
    for s in Subsets(allvars):
        c=[]
        for k in cons.keys():
            cvars=cons[k][0]|cons[k][1]
            if cvars.issubset(s):
                c.append(k)
        apcns[tuple(sorted(list(s)))]=c
    return apcns
    
def testcons(M1,pcode,nsrc, appcns,newvar,netcons):
    r"""
    Tests whether a matroid satisfies network constrains under given map
    
    INPUT:
    
    ``M1`` - A representable matroid
    ``pcode`` - A dictionary giving a mapping from ground set of ``M1``
    to the random variable indices 
    ``nsrc`` - Number of sources in a network
    ``appcns`` - A dictionary mapping subsets of random variables (as 
    sorted tuples) to subset of constraint labels in ``netcons.keys()``
    ``newvar`` - newest member of ``pcode.values()`` for which the cons-
    -raints are to be tested
    -``netcons`` - A dictionary specifying network constraints as 
    ``[list1,list2]`` where ``list1`` and ``list2`` are subsets of 
    random variable indices that are forced to have equal rank (entropy)
    
    OUTPUT:
    
    - A boolean indicating whether definition of ``newvar`` satisfies 
    network constraints
    
    .. NOTE::

            This method does NOT do any checks.
    """
    newcons= set(appcns[tuple(sorted(list(set(pcode.values()))))])-set(appcns[tuple(sorted(list(  set(pcode.values())-set([newvar]) )))])
    inv_pcode = {v: k for k, v in pcode.items()}
    if(len(newcons)>0):
        for c in newcons:
            if M1.rank(set([inv_pcode[v] for v in netcons[c][0]]))!=M1.rank(set([inv_pcode[v] for v in netcons[c][1]])):
                return 0
    if newvar<=nsrc:
        if set(range(1,nsrc+1)).issubset(set(pcode.values())):
            srcgnd=[k for k in pcode.keys() if pcode[k]<=nsrc]
            ssum=0
            for i in srcgnd:
                ssum=ssum+M1.rank(set([i]))
            if ssum!=M1.rank(set(srcgnd)):
                return 0
    return 1


    
def extend_matroids_simple(list_of_matroids):
    r"""
    Returns a list of simple single element linear extensions of matroids.

    INPUT:
    
    -``list_of_matroids`` - A list of simple representable matroids
    
    OUTPUT:
    
    A list of pairs ``[Me,p]`` where ``Me`` is a simple linear extension
    of matroid at index ``p`` of ``list of matroids``
    
    .. NOTE::

            This method does NOT do any checks.  
    """
    matroidset=set([])
    extmats=[]
    for p in range(len(list_of_matroids)):
        Mp=list_of_matroids[p]
        for Mc1 in Mp.linear_extensions(element=max(Mp.groundset())+1,simple=True):
            bad=0
            for Mc2 in matroidset:
                if Mc2.is_isomorphic(Mc1):
                    bad=1
                    break
            if bad!=1:
                matroidset=matroidset|set([Mc1])
                extmats.append([Mc1,p])
    return extmats



def certsearch_dfs_withrates(M1,pcode,rate_vector,netcons,appcns,d,nsrc,nvars):
    r"""
    Return the lexicographically smallest matroid-network map that forms
    a partial code 
    
    INPUT:
    
    ``M1`` - A matroid
    ``pcode`` - partial code
    ``netcons`` - A dictionary specifying network constraints as 
    ``[list1,list2]`` where ``list1`` and ``list2`` are subsets of 
    random variable indices that are forced to have equal rank (entropy)
    ``appcns`` - A dictionary mapping subsets of random variables (as 
    sorted tuples) to subset of constraint labels in ``netcons.keys()``
    ``d`` - depth in the DFS tree, used for recursion
    ``nsrc`` - number of sources in the network
    ``nvars`` - number of network random variables
    
    OUTPUT:
    
    A 2-tuple containing a boolean and a dictionary
    
    .. NOTE::

            This method does NOT do any checks. 
    """
    ret=False
    if len(pcode.values())>d:
        ret,pcode1=certsearch_dfs_withrates(M1,pcode,rate_vector,netcons,appcns,d+1,nsrc,nvars)
    if ret==True:
        return ret,pcode1
    else:
        try_g=sorted(list(M1.groundset()))[d]
        if try_g in pcode.keys():
            lex_pivot=pcode[try_g]
        else:
            lex_pivot=0
        
        #clean up pcode
        pcode={g:pcode[g] for g in pcode.keys() if g<try_g}
        lex_vars=sorted([v for v in list(set(range(1,nvars+1))-set(pcode.values())) if v>lex_pivot])
        for v in lex_vars:
            pcode[try_g]=v
            # check if defn works
            if M1.rank([try_g])==rate_vector[v-1] and testcons(M1,pcode,nsrc,appcns,v,netcons)==1:
            #def worked
                # are we finished?
                if len(pcode.values())==len(M1.groundset()):
                    return True,pcode
                else:
                    # go deeper in the tree
                    ret,pcode1=certsearch_dfs_withrates(M1,pcode,rate_vector,netcons,appcns,d+1,nsrc,nvars)
                    if ret==1:
                        return True,pcode1
            else:
                pcode.pop(try_g)
        return False,pcode










def parallel_loopy_extensions(Mlist):
    r"""
    Returns all linear extensions of a list of matroids obtained by 
    adding parallel elements and loops
    
    INPUT:
    
    - ``Mlist`` - A list of representable matroids
    
    OUTPUT:
    
    - A list of pairs ``[Me,p]`` where ``Me`` is a simple linear extension
    of matroid at index ``p`` of ``list of matroids``
    
    .. NOTE::

            This method does NOT do any checks.
    """
    Mset=[]
    for i in range(len(Mlist)):
        M1=Mlist[i]
        newgnd=max(M1.groundset())+1
        extvecs=M1.simplify().representation_vectors()
        extvecs[max(extvecs.keys())+1]=[0]*M1.rank()
        for k in extvecs.keys():
            M2=M1.linear_extension(newgnd,col=extvecs[k])
            goodmat=True
            for j in range(len(Mset)):
                M3=Mset[j][0]
                if M2.is_isomorphic(M3):
                    goodmat=False
                    break
            if goodmat==true:
                Mset.append([M2,i])
    return Mset


def idmatrix_codecert_withrates(rate_vector,nsrc,q):
    r"""
    .. NOTE::

            This method does NOT do any checks.
    """
    nonloopysrc=[i for i in range(1,nsrc+1) if rate_vector[i-1]==1 ]
    loopysrc=[i for i in range(1,nsrc+1) if rate_vector[i-1]==0]
    nloopysrc=len(loopysrc)
    Fq= GF(q,'a')
    I=identity_matrix(Fq,nsrc-nloopysrc)
    C=matrix(Fq,nsrc-nloopysrc,nloopysrc)
    T=block_matrix([[I,C]])
    M=Matroid(T)
    pcode={}
    k=0
    l=len(nonloopysrc)
    for i in range(1,nsrc+1):
        if rate_vector[i-1]==1:
            pcode[k]=i
            k=k+1
        else:
            pcode[l]=i
            l=l+1
    return M,pcode

def parallel_loopy_extensions(Mlist,loopy):
    r"""
    Returns all linear extensions of a list of matroids obtained by 
    adding parallel elements and loops
    
    .. NOTE::

            This method does NOT do any checks.
    """
    Mset=[]
    for i in range(len(Mlist)):
        M1=Mlist[i]
        newgnd=max(M1.groundset())+1
        extvecs=M1.simplify().representation_vectors()
        if loopy==True:
            extvecs[max(extvecs.keys())+1]=[0]*M1.rank()
        for k in extvecs.keys():
            M2=M1.linear_extension(newgnd,col=extvecs[k])
            goodmat=True
            for j in range(len(Mset)):
                M3=Mset[j][0]
                if M2.is_isomorphic(M3):
                    goodmat=False
                    break
            if goodmat==True:
                Mset.append([M2,i])
    return Mset

def matroid2degrees(M):
    r"""
    Return a vector specifying the degree of each unique column in the 
    matroid representation
    
    INPUT:
    
    - ``M`` -- A representable matroid
    
    OUTPUT:
    
    A vector specifying number of times each coulumn occurs in 
    ``M.representation()``, in sorted order.
    
    .. NOTE::

            This method does NOT do any checks.
    """
    repdict=M.representation_vectors()
    d=[repdict.values().count(vector([0]*M.rank()))]
    for v in sorted(M.simplify().representation_vectors().values()):
        d.append(repdict.values().count(v))
    d=vector(d)
    d.set_immutable()
    return d

def findcode_with_seed(M,pcode,rate_vector,netcons,appcns,known_bad,nsrc,nvars):
    r"""
    Recursively test if a matroid has any sequence of parallel and/or loopy 
    linear extensions achieving the specified rate vector
    
    INPUT:
    
    ``M`` - A representable matroid
    ``pcode`` - A dictionary giving a mapping from ground set of ``M1``
    to the random variable indices 
    ``nsrc`` - Number of sources in a network
    ``appcns`` - A dictionary mapping subsets of random variables (as 
    sorted tuples) to a subset of constraint labels in ``netcons.keys()``
    ``newvar`` - newest member of ``pcode.values()`` for which the cons-
    -raints are to be tested
    -``netcons`` - A dictionary specifying network constraints as 
    ``[list1,list2]`` where ``list1`` and ``list2`` are subsets of 
    random variable indices that are forced to have equal rank (entropy)
    
    OUTPUT:
    
    A 3-tuple containing:

    1. A boolean indicating whether ``M`` has a  sequence of parallel and/or 
    loopy linear extensions achieving ``rate_vector``
    2. If 1. is ``True``, pair ``[Mx,pcode]`` where  ``Mx`` is matroid of size
    ``nvars`` obtained from  ``M`` via a sequence of parallel and/or loopy 
    linear extensions 
    3. A list of bad linear extensions specified as vectors of element degrees
    for future reference during recursion 
    
    
    
    .. NOTE::

            This method does NOT do any checks.
    """
    # see if there is any extension of M that is feasible code
    if len(M.groundset())==len(rate_vector):
        return True,[M,pcode],known_bad
    extmats=parallel_loopy_extensions([M],len(M.loops())<len([i for i in rate_vector if i==0]))
    for mat in extmats:
        if matroid2degrees(mat[0]) not in known_bad: 
            ret1,pcode1=certsearch_dfs_withrates(mat[0],pcode,rate_vector,netcons,appcns,nsrc,nsrc,nvars)
            if ret1==True:
                ret2,pcode2,known_bad=findcode_with_seed(mat[0],pcode1,rate_vector,netcons,appcns,known_bad,nsrc,nvars)
                if ret2==True:
                    return True,pcode2,known_bad
            else:
                known_bad.append(matroid2degrees(M))
    return False,[M,pcode],known_bad

def conslist2dict(cons):
    r"""
    Return a dictionary containing constraints
    
    INPUT:
    
    - ``cons`` -- A list of lists specifying network constraints
    
    OUTPUT:
    
    - A dictionary with integers as keys and constraints as values
    
    EXAMPLES::
    
        sage: from sage.matroids.network_coding_helpers import *
        sage: from sage.matroids.networks_catalog import *
        sage: cons=[[[1,2],[1,2,3]],[[1,3],[1,2,3]],[[2,3],[1,2,3]]]
        sage: conslist2dict(cons)
        {0: [{1, 2}, {1, 2, 3}], 1: [{1, 3}, {1, 2, 3}], 2: [{2, 3}, {1, 2, 3}]}
    
    .. NOTE::

            This method does NOT do any checks.
    """
    d={}
    for i in range(len(cons)):
        d[i]=[set(cons[i][0]),set(cons[i][1])]
    return d

def ratecertgen(netcons,rate_vector,nsrc,nvars,q):
    r"""
    Test if a rate vector is achievable with scalar linear network codes
    over a specific field
    
    INPUT:
    
    - ``netcons`` -- A dictionary specifying network constraints
    - ``rate_vector`` -- A 0-1 vector whose achievability is to be tested
    -  ``nsrc`` -- Number of sources in the network
    - ``q`` -- Size of finite field over which achievability is to be tested
    
    OUTPUT:
    
    A 2-tuple containing:
    
    1. A boolean indicating whether ``rate_vector`` is achievable with scalar 
    linear network coding over finite field of size ``q``.
    2. If 1. is ``True``, a pair ``[M,pcode]`` where ``M`` is a matroid of size
    ``nvars`` representable over finite field of size ``q`` and ``pcode`` is a 
    dictionary mapping groundset elements to random variable associated with 
    the network
    
    EXAMPLES::
    
        sage: from sage.matroids.network_coding_helpers import *
        sage: from sage.matroids.networks_catalog import *
        sage: N=FanoNet()
        sage: r,code=ratecertgen(N.constraints,[1,1,1,1,1,1,1],N.nsrc,N.size,2)
        sage: r
        True
        sage: code
        [Binary matroid of rank 3 on 7 elements, type (3, 0),
         {0: 1, 1: 2, 2: 3, 3: 4, 4: 6, 5: 7, 6: 5}]
        sage: code[0].representation()
        [1 0 0 1 1 1 0]
        [0 1 0 1 0 1 1]
        [0 0 1 0 1 1 1]
        sage: code[0].is_isomorphic(matroids.named_matroids.Fano())
        True
        sage: r,code=ratecertgen(N.constraints,[1,1,1,1,1,1,1],N.nsrc,N.size,3)
        sage: r
        False

    
    .. NOTE::

            This method does NOT do any checks.
    """
    max_nsimple=len([r for r in rate_vector if r==1])
    min_nsimple=len([i for i in range(nsrc+1) if rate_vector[i]==1])
    simple_pcodes=[idmatrix_codecert_withrates(rate_vector,nsrc,q)]
    appcns=appcons(netcons)
    for i in range(min_nsimple,max_nsimple+1):
        if len(simple_pcodes)==0:
            return False,[]
        simple_pcodes_new=[]
        extmats=extend_matroids_simple([code[0] for code in simple_pcodes])
        known_bad=[]
        for mtr in extmats:
            pcode=simple_pcodes[mtr[1]][1]
            ret,pcode=certsearch_dfs_withrates(mtr[0],pcode,rate_vector,netcons,appcns,0,nsrc,nvars)
            if ret==True:
                simple_pcodes_new.append([mtr[0],pcode])
                # now search via parallel/loopy extensions of mtr[0]
                ret1,pcode1,known_bad=findcode_with_seed(mtr[0],pcode,rate_vector,netcons,appcns,known_bad,nsrc,nvars)
                if ret1==True:
                    # we are done
                    return True,pcode1
        simple_pcodes=simple_pcodes_new
    return False,[]
