from sage.rings.finite_rings.finite_field_prime_modn import *
from sage.matroids.advanced import *
from sage.sets.set import Set
from sage.all import *
import numpy as np
from itertools import izip

msnccons=[  [ [[],[]], [[],[]], [[],[]] ],[] ]
#R=None
GF2=GF(2)


cpdef all_1def_pcodes(Mpcode,nc):
    varset = ncinstance_vars(nc)
    for v in varset:
        new_pmap = {0:v,1:v,2:v}
        if is_pcode(Mpcode._M,new_pmap,nc):
             Mpcode._C.append(new_pmap)
    return Mpcode

cpdef all_2def_pcodes(Mpcode,nc):
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
                 addpcode(Mpcode._C,new_pmap)
            else:
                for var in set(new_pmap.values()):
                    new_pmap_del = new_pmap.copy()
                    new_pmap_del = {key:new_pmap[key] for key in new_pmap.keys() if new_pmap[key] != var}
                    if len(new_pmap_del.keys()) > 1:
                        if is_pcode(Mpcode._M,new_pmap_del,nc):
                            addpcode(Mpcode._C,new_pmap_del)
    return Mpcode

cpdef all_3def_pcodes(Mpcode,nc):
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
                 #Mpcode._C.append(new_pmap)
                 addpcode(Mpcode._C,new_pmap)
                 n3pcode+=1
            else:
                for var in set(new_pmap.values()):
                    new_pmap_del = {key:new_pmap[key] for key in new_pmap.keys() if new_pmap[key] != var}
                    if is_pcode(Mpcode._M,new_pmap_del,nc):
                        addpcode(Mpcode._C,new_pmap_del)
                    else:
                        for var2 in set(new_pmap_del.values()):
                            new_pmap_del2 = new_pmap_del.copy()
                            for key in new_pmap_del:
                                if new_pmap_del2[key] == var2:
                                    new_pmap_del2.pop(key)
                            if is_pcode(Mpcode._M,new_pmap_del2,nc):
                                addpcode(Mpcode._C,new_pmap_del2)
    return Mpcode

cpdef pmap_def2():
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

cpdef pmap_def3():
    pmaps=[]
    p1 = OrderedSetPartitions({0,1,2})[0:6]
    for p in p1:
        e1 = list(p[0])[0]
        e2 = list(p[1])[0]
        e3 = list(p[2])[0]
        pmap = {e1:0,e2:1,e3:2}
        pmaps.append(pmap)
    return pmaps

class MatroidCodes:
    def __init__(self,M=None,pmaps=None):
        self._M=M
        self._C=[]
        self._parent=None
        if M != None and pmaps != None:
            self.add(pmaps)

    def add(self, pmaps):
        for pmap in pmaps:
            addpcode(self._C, pmap)

cpdef addpcode(list_of_pcodes,pcode):
    """
    adds p-code to a list of p-codes if there is no p-code in ``list_of_pcodes`` that is extension
    of ``pcode``
    """
    #print 'here0'
    for cpcode in list_of_pcodes:
        if is_pmap_deletion(pcode,cpcode):
            #print 'here1'
            return [1,cpcode]
        elif is_pmap_deletion(cpcode,pcode):
            list_of_pcodes.remove(cpcode)
            list_of_pcodes.append(pcode)
            #print 'here2'
            return [2,cpcode]
    # brand new pcode
    #print 'here3'
    list_of_pcodes.append(pcode)
    return [0,None]

cpdef is_pmap_deletion(dict1,dict2):
    """
    checks whether $dict1\subseteq dict2$
    """
    if len(dict1) > len(dict2):
        return False
    if not set(dict1.values()).issubset(dict2.values()):
        return False
    else:
        idict1=invert(dict1)
        idict2=invert(dict2)
        # check variable set containment
        for key in idict1.keys():
            if idict1[key] != idict2[key]: #non-matching variable definitions
                return False
        return True

cpdef invert(m):
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

cpdef add2(R,M):
    for N in R:
        if M.is_field_isomorphic(N):
            return
    R.append(M)

cpdef extensions2(M):
    R=[]
    x=len(M)
    for N in M.linear_extensions(len(M)):
        add2(R,N)
    return R

cpdef next2(R):
    R2=[]
    for M in R:
        for N in extensions2(M):
            add2(R2,N)
    return R2

cpdef add2_p(R,M):
    # add2 with parent tracking
    for N in R:
        if M.is_field_isomorphic(N):
            return False
    R.append(M)
    return True

cpdef extensions2_p(M):
    # extensions2 with parent tracking
    R=[]
    x=len(M)
    for N in M.linear_extensions(len(M)):
        add2(R,N)
    return R

cpdef next2_p(R):
    # next2 with parent tracking, returning Mcodes instead of matroids
    R2=[]
    parents=[]
    i=0
    for M in R:
        for N in extensions2(M):
            if add2_p(R2,N) == True:
                parents.append(i)
        i=i+1
    print parents
    return R2,parents


cpdef all2(N, R=10):
    MM={}
    for n in xrange(N+1):
        MM[n,0]=[BinaryMatroid(identity_matrix(GF2,n))]
        for r in xrange(min(n,R+1)):
            MM[r,n-r]=next2(MM[r,n-r-1])
        #print 'n=',n
        #print [len(MM[r,n-r]) for r in xrange(min(n,R)+1)]
    return MM

cpdef ncmatenum(nc,maxgnd):
    nvars = len(ncinstance_vars(nc))
    mcodes_3 = init_enum(nc,maxgnd)
    mcodes_parent = mcodes_3
    for g in range(4,maxgnd+1):
        mcodes_parent=extend_mcodes(mcodes_parent,nc,maxgnd,nvars)

cpdef add_identity_codes(mcode_n_n,nc,maxgnd):
    """ Return new mcode with to size n+1 rank n+1 matroid
    that corresponds to coloop extension of size n rank n
    matroid. The p-codes for this matroid are computed in p-codes
    of input mcode
    """
    n=len(mcode_n_n._M.groundset())
    child_Mcode = MatroidCodes(BinaryMatroid(identity_matrix(GF2,n+1)),[])
    extendpmaps(child_Mcode,mcode_n_n,nc,maxgnd,len(ncinstance_vars(nc)))
    return child_Mcode

cpdef dfz_matroidal_network(M):
    """
    Create a matroidal network based on Dougherty, Freiling and Zeger's
    Construction
    """
    nodes=[]
    msgs=[]
    dag={}
    f={}
    g={}
    allB=M.bases()
    B1=list(allB[0])
    j=0 # edge count
    # STEP 1: create source nodes in dag and assign them msgs and B1 elements
    for i in xrange(len(B1)):
        dag[i] = []
        msgs.append(i)
        f[i]=B1[i]
        g[B1[i]] = i
    nsrc = i+1
    i+=1 # current node index
    j = i # msg/edge index
    print 'dag',dag
    print 'f',f
    print 'g',g
    used_ckts=set([])
    # STEP 2: find circuits with 1 undefined (in terms of g) element
    while True:
        c = candidate_circuit(M.circuits(),g)
        if c == None:
            break
        else:
            print '2.enforce', c
            used_ckts.add(c)
            # Step 2.(i)
            dag[i]=[]
            msgs.append(i)
            for cx in set(c)-set(set(c)-set(g.keys())):
                f[(g[cx],i)]=cx
                j+=1
                dag[g[cx]].append(i)
            i+=1
            # Step 2.(ii)
            dag[i]=[]
            dag[i-1].append(i)
            f[(i-1,i)] = list(set(c)-set(g.keys()))[0]
            g[list(set(c)-set(g.keys()))[0]] = i
            j+=1
            i+=1
    # STEP 3: for every source msg and containing circuit pair add a receiver
    demands = {}
    used_src=set([])
    for x in M.groundset():
        if g[x] < nsrc:
            for ckt in M.circuits():
                if x in ckt and g[x] not in used_src:
                    used_src.add(g[x])
                    print '3.enforce', ckt,x
                    used_ckts.add(ckt)
                    demands[i] = [g[x]]
                    for y in set(ckt)-set([x]):
                        dag[g[y]].append(i)
                        f[(g[y],i)]=y
                        j+=1
                    i+=1
                if used_src==set(xrange(nsrc)):
                    break
        if used_src==set(xrange(nsrc)):
            break
    # STEP 4:
    print used_ckts
    for ckt in M.circuits():
        if ckt not in used_ckts:
            # encorce dependency
            print '4 WAS USED'
            B = M.basis()
            for b in B:
                dag[g[b]].append(i)
                f[(g[b],i)] = b
                j+=1
            demands[i] = range(nsrc)
            i+=1
    return (dag,f,g,demands)

cpdef dfz_ncinstance(M):
    """ return ncinstance corresponding to dfz matroidal network
    """
    dag,f,g,demands = dfz_matroidal_network(M)
    srcs = [s+1 for s in f.keys() if isinstance(s,type(1))]
    print srcs
    nc = [[],[]]
    nc[1].extend(srcs)
    print 'here'
    for src in srcs:
        in_set = set([src])
        out_set = set([ i+1 for i in dag[src]])
        nc[0].append([list(in_set),list(in_set|out_set)])
    for node in set(dag.keys())-set([s-1 for s in srcs]):
        in_set = set([i+1 for i in dag.keys() if  node in dag[i]])
        out_set = set([i+1 for i in dag[node]])
        nc[0].append([list(in_set),list(in_set|out_set)])
        print 'node',node+1,[list(in_set),list(in_set|out_set)]
    for d in demands.keys():
        in_set = set([i+1 for i in dag.keys() if  d in dag[i]])
        out_set= set([i+1 for i in demands[d]])
        nc[0].append([list(in_set),list(in_set|out_set)])
        print 'd',d,[list(in_set),list(in_set|out_set)]
    return nc

cpdef candidate_circuit(circuits,g):
    for c in circuits:
        if len(set(c)-set(g.keys())):
            return c
    return None

cpdef max_rank(nc,maxgnd):
    """
    Return maximum rank of matroid that must be considered
    """
    maxr = min(maxgnd-ncinstance_vars(nc)+len(nc[1]),floor(maxgnd/2))

cpdef init_enum(nc,N):
    maxrank = N-(len(ncinstance_vars(nc))-len(nc[1]))
    M3 = all2(3,maxrank)
    for key in M3.keys():
        if sum(list(key)) != 3:
            M3.pop(key)
    allvars=ncinstance_vars(nc)
    size_3_matroidcodes = {}
    for key in M3:
        size_3_matroidcodes[key] = []
        for M in M3[key]:
            Mpcode = MatroidCodes(M,[])
            Mpcode = all_3def_pcodes(Mpcode,nc)
            Mpcode = all_2def_pcodes(Mpcode,nc)
            Mpcode = all_1def_pcodes(Mpcode,nc)
            if len(Mpcode._C) > 0:
                size_3_matroidcodes[key].append(Mpcode)
    return size_3_matroidcodes

cpdef init_enum_3subset(nc,N,subset):
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

cpdef ncinstance_vars(nc):
    varset=Set([])
    for i in nc:
        for j in nc[0]:
            v=Set(j[1])
            varset=varset.union(v)
    varset.union(Set(nc[1]))
    return varset

cpdef extend_mcodes(Mcodes1,nc,maxgnd,nvars):
    Mcodes2 = {}
    if len(Mcodes1.keys()) < 1:
        return Mcodes2
    n1 = sum(list(Mcodes1.keys()[0]))
    n2 = n1+1
    for key in Mcodes1:
        print 'extending',key
        Mcodes2[(list(key)[0],n2-list(key)[0])] = []
        # create non-isomorphic extensions of given key matroids
        key_matroids = [Mcode._M for Mcode in Mcodes1[key]]
        key_ext,parents= next2_p(key_matroids)
        for i in xrange(len(key_matroids)):
            children_ind = [k for k in xrange(len(key_ext)) if parents[k]==i]
            children = [key_ext[k] for k in xrange(len(key_ext)) if parents[k]==i]
            print 'parent', Mcodes1[key][parents[k]]._M.representation()
            for k in children_ind:
                child_Mcode = MatroidCodes(key_ext[k],[])
                extendpmaps(child_Mcode,Mcodes1[key][parents[k]],nc,maxgnd,nvars)
                if(len(child_Mcode._C)>0):
                    Mcodes2[(list(key)[0],n2-list(key)[0])].append(child_Mcode)
    return Mcodes2

cpdef mdcs(config):
    """Return mdcs instance corresponding to configuration matrix specified as
    a list of lists"""
    if mdcs_instance_is_valid(config) == False:
        raise ValueError('Not a valid MDCS instance')
    else:
        # create constraints from the matrix
        nc = [[],[]]
        nc[1].extend(range(1,len(config)+1))
        for l in xrange(len(config)):
            level=config[l]
            lsets_all=[set([len(k.binary())-i+len(config) for i in xrange(len(k.binary()))
                       if k.binary()[i]=='1']) for k in level]
            lsets=[s for s in lsets_all if len(s)>0]
            for lset in lsets:
                nc[0].append([list(lset),list(lset|set(range(1,l+2)))])
        return nc

cpdef mdcs_instance_is_valid(config_list):
    """
    Test the validity of a given MDCS instance
    """
    if len(config_list) == 0:
        return False
    max_dec_per_level=len(config_list[0])
    for level in xrange(1,len(config_list)):
        if len(config_list[level]) != max_dec_per_level:
            return False
    # test axioms (C1) and (C5)
    all_levels=[]
    all_enc = set([])
    for level in config_list:
        lsets_all=[set([len(k.binary())-i for i in xrange(len(k.binary())) if k.binary()[i]=='1']) for k in level]
        lsets=[s for s in lsets_all if len(s)>0]
        all_enc=all_enc.union(*lsets)
        # (C5)
        if len(lsets) == 0:
            return False
        all_levels.append(lsets)
        # (C1)
        if len(lsets)>1:
            print lsets
            for x in lsets:
                lsets2=copy(lsets)
                lsets2.remove(x)
                if any([x >= y or x <= y for y in lsets2]):
                    return False
    print 'C1,C5 pass'
    # test (C3)
    print all_enc
    encstr='1'*max(all_enc)
    if set([len(encstr)-i for i in xrange(len(encstr)) if encstr[i]=='1']) != all_enc:
        return False
    print  'C3 pass'
    # test (C4)
    all_enc_fans=[]
    for enc in all_enc:
        e_dec=set([])
        for i in xrange(len(all_levels)):
            l = all_levels[i]
            l_e = set([i+dec for dec in xrange(len(l)) if enc in l[dec]])
            e_dec=e_dec.union(l_e)
        all_enc_fans.append(e_dec)
    print all_enc_fans
    for efan in all_enc_fans:
        others=all_enc_fans[:]
        print efan,others
        others.remove(efan)
        for o in others:
            if efan==o:
                return False
    print 'C4 pass'
    # test axiom (C2)
    print all_levels
    for i in xrange(1,len(all_levels)):
        # test (C2) for all lower levels
        l_i=all_levels[i]
        for j in xrange(i):
            l_j=all_levels[j]
            for lset_i in l_i:
                for lset_j in l_j:
                    if lset_j >= lset_i:
                        print i,j,l_i,l_j,lset_j,lset_i
                        return False
    print 'C2 pass'
    ## print lsets
    return True



cpdef Mcode_candidates(Mcode,maxgnd,nvars):
    candidates = []
    i=0
    for mcode in Mcode._C:
        # loop over all subsets of defined variables and delete them
        imcode = invert(mcode)
        defvars = set(imcode.keys())
        for delset in Subsets(defvars):
            imcode_copy = imcode.copy()
            del_e=[]
            for v in delset:
                del_e.append(imcode_copy.pop(v))
            # get candidate pmap
            cand = invert(imcode_copy)
            #print cand
            if nvars-len(set(cand.values())) <=  (maxgnd-len(Mcode._M.groundset())):
                candidates.append([cand,del_e])
        i+=1
    # keep only the unique candidates
    #unique_ind=list(np.unique(np.array([cand[0] for cand in candidates])))
    return candidates#[candidates[k] for k in unique_ind]

cpdef applicable_cons(def_vars,nc):
    cons=[]
    for j in nc[0]:
         if all([k in def_vars for k in j[1]]):
             cons.append(j)
    return cons



cpdef is_pcode(M,pmap,ncinstance):
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
        src_joint = M.rank(set([]).union(*[ipmap[s] for s in def_src]))
        if src_sum != src_joint:
            return False
    return True




cpdef extendpmaps(child_Mcode,Mcode,nc,maxgnd,nvars):
    """
    Populate child._C with maximal pmap extensions that form p-codes
    """
    child_Mcode._C=[]
    # create a list of unique candidates from all mcodes
    unique_candidates = Mcode_candidates(Mcode,maxgnd,nvars)
    for u in unique_candidates:
        cand=u[0]
        # loop over subsets U of unmapped gndset that have >=1 sized
        # intersection with each subset in u[1]
        for U in Subsets(set(Mcode._M.groundset())-set(cand.keys())):
            if U_is_valid(U,u[1]):
                Ue=U|Set(child_Mcode._M.groundset()-Mcode._M.groundset())
                # loop over all new variable definitions
                for v in set(xrange(1,nvars+1))-set(cand.values()):
                    for e in Ue:
                        cand[e]=v
                    if is_pcode(child_Mcode._M,cand,nc):
                        child_Mcode._C.append(cand.copy())
                        #addpcode(child_Mcode._C,cand.copy())
                    for e in Ue:
                        cand.pop(e)
    return

cpdef U_is_valid(U,delset_maps):
    #test if delset
    for mapset in delset_maps:
        print mapset
        if mapset != None and len(U & mapset) < 1:
            return False
    return True

cpdef naive_candidates(M_ext,list_of_pmaps):

    cand = []

cpdef butterfly():
 return [[ [[3,4],[1,3,4]],[[5,6],[1,5,6]],[[7],[4,5,7]],[[8,9],[7,8,9]],[[2],[2,3,8]],[[1],[1,6,9]] ] ,[1,2]]




