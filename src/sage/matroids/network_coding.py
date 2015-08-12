r"""
A scalar linear network coding achievability prover 

<Paragraph description>

AUTHORS:

- Jayant Apte (2015-07-08): initial version

EXAMPLES::

<Lots and lots of examples>

REFERENCES
==========


..  [Yeung] Raymond W. Yeung. 2008. Information Theory and Network Coding (1 ed.). Springer Publishing Company, Incorporated.
..  [Oxley] James Oxley, "Matroid Theory, Second Edition". Oxford University Press, 2011.

"""

#*****************************************************************************
#       Copyright (C) 2015 Jayant Apte <jayant91089@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.sage_object import SageObject
from sage.matroids.advanced import BinaryMatroid, TernaryMatroid, QuaternaryMatroid
from sage.all import Subsets,GF

class NCinstance(SageObject):
     r"""
     
     The network coding instance class
     
     Network Coding is a paradigm for communication over networks 
     (assumed to be directed acyclic hypergraphs) with error-free links.
     A network code for a network coding instance of size `n` is a 
     collection of `n` discrete random variables, satisfying constraints 
     imposed on entropies of subsets by the network coding instance. 
     """
     def __init__(self, size, nsrc, constraints):
         assert size > 1
         assert nsrc >= 1
         self.size = int(size)  # No. of random variables
         self.nsrc = int(nsrc)  # No. of Sources
         self.constraints = constraints  # constraints
         
     def __repr__(self):
         if len(self.constraints) == 0:
             return "An empty network coding instance of size %s"%(self.size)
         else:
             return "A network coding instance of size %s with %s sources"%(
                     self.size, self.nsrc)
     def empty(self):
         self.constraints = [];



def appcons(cons):
    r"""
	Return a dictionary of network constraints applicable to each
    
    INPUT:
    
    - A list of lists specifying network constraints specified as 
    ``[list1,list2]`` where ``list1`` and ``list2`` are subsets of 
    random variable indices that are forced to have equal rank (entropy)   
	"""
    allvars=set([])
    for k in cons.keys():
        allvars=allvars|cons[k][0]|cons[k][1]
    apcns={}
    for s in Subsets(allvars):
        #print 's',s
        c=[]
        for k in cons.keys():
            cvars=cons[k][0]|cons[k][1]
            if cvars.issubset(s):
                c.append(k)
        #print c
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
	"""
    newcons= set(appcns[tuple(sorted(list(set(pcode.values()))))])-set(appcns[tuple(sorted(list(  set(pcode.values())-set([newvar]) )))])
    #print appcns
    #print 'test',newcons,appcns[tuple(sorted(list(set(pcode.values()))))],appcns[tuple(sorted(list(  set(pcode.values())-set([newvar]) )))]
    inv_pcode = {v: k for k, v in pcode.items()}
    if(len(newcons)>0):
        for c in newcons:
            #print netcons[c][0]
            #print M1.rank(set(inv_pcode[v] for v in netcons[c][0]))
            #print netcons[c][1]
            #print M1.rank(set(inv_pcode[v] for v in netcons[c][1]))
            if M1.rank(set([inv_pcode[v] for v in netcons[c][0]]))!=M1.rank(set([inv_pcode[v] for v in netcons[c][1]])):
                #print 'fail1!!!'
                return 0
    if newvar<=nsrc:
        if set(range(1,nsrc+1)).issubset(set(pcode.values())):
            #print 'here2',set(range(1,nsrc+1)),set(pcode.values())
            srcgnd=[k for k in pcode.keys() if pcode[k]<=nsrc]
            ssum=0
            #print srcgnd,M1.groundset()
            for i in srcgnd:
                ssum=ssum+M1.rank(set([i]))
            if ssum!=M1.rank(set(srcgnd)):
                #print 'fail2!!!',ssum,M1.rank(set(srcgnd)),M1.representation()
                return 0
    return 1

def extend_matroids(list_of_matroids):
    r"""
	Returns a list of single element linear extensions of matroids
    
    INPUT:
    
    -``list_of_matroids`` - A list of representable matroids
	"""
    matroidset=set([])
    extmats=[]
    for p in range(len(list_of_matroids)):
        Mp=list_of_matroids[p]
        for Mc1 in Mp.linear_extensions(element=max(Mp.groundset())+1):
            bad=0
            for Mc2 in matroidset:
                if Mc2.is_isomorphic(Mc1):
                    bad=1
                    break
            if bad!=1:
                matroidset=matroidset|set([Mc1])
                extmats.append([Mc1,p])
    return extmats

def extend_matroids_simple(list_of_matroids):
    r"""
	Returns a list of simple single element linear extensions of matroids
	
    INPUT:
    
    -``list_of_matroids`` - A list of simple representable matroids
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

def idmatrix_allcodes(r,q,netcons,appcns,nsrc,nvars):
    r"""
	Return all partial codes for an identity matrix
    
    INPUT:
    
    ``r`` - rank of identity matrix
    ``q`` - size of finite field
    ``netcons`` - A dictionary specifying network constraints as 
    ``[list1,list2]`` where ``list1`` and ``list2`` are subsets of 
    random variable indices that are forced to have equal rank (entropy)
    ``appcns`` - A dictionary mapping subsets of random variables (as 
    sorted tuples) to subset of constraint labels in ``netcons.keys()``
    ``nsrc`` - number of sources in the network
    ``nvars`` - number of network random variables
	"""
    if q==2:
        mat=identity_matrix(GF2,r)
        M1=BinaryMatroid(mat)
    if q==3:
        mat=identity_matrix(GF3,r)
        M1=TernaryMatroid(mat)
    if q==4:
        mat=identity_matrix(GF4,r)
        M1=QuaternaryMatroid(mat)
    #print M1.representation()
    codelist=[{}]
    for g in M1.groundset():
        codelist1=codelist
        codelist=[]
        for pcode in codelist1:
            for v in set(range(1,nvars+1))-set(pcode.values()):
                pcode[g]=v
                #print pcode
                if testcons(M1,pcode,nsrc,appcns,v,netcons)==1:
                    codelist.append(deepcopy(pcode))
                pcode.pop(g)
    return M1,codelist

def idmatrix_codecert(r,q,netcons,appcns,nsrc,nvars):
    r"""
	Return the lexicographically smallest partial code for a network
    
    INPUT:
    
    ``r`` - rank of identity matrix
    ``q`` - size of finite field
    ``netcons`` - A dictionary specifying network constraints as 
    ``[list1,list2]`` where ``list1`` and ``list2`` are subsets of 
    random variable indices that are forced to have equal rank (entropy)
    ``appcns`` - A dictionary mapping subsets of random variables (as 
    sorted tuples) to subset of constraint labels in ``netcons.keys()``
    ``nsrc`` - number of sources in the network
    ``nvars`` - number of network random variables
	"""
    Fq = GF(q,'a')
    mat=identity_matrix(Fq,r)
    M1=Matroid(mat)
    codelist=[{}]
    lex_gnd=sorted(list(M1.groundset()))
    for g in lex_gnd:
        codelist1=codelist
        codelist=[]
        for pcode in codelist1:
            lex_vars=sorted(list(set(range(1,nvars+1))-set(pcode.values())))
            for v in lex_vars:
                pcode[g]=v
                #print pcode
                if testcons(M1,pcode,nsrc,appcns,v,netcons)==1:
                    codelist.append(deepcopy(pcode))
                    if g==lex_gnd[-1]:
                        return M1,pcode
                pcode.pop(g)
    return M,None

def extendpcode(M1,pcode,newgnd,netcons,appcns,nsrc,nvars):
    r"""
	extend the partial code
	"""
    codelist=[]
    for v in set(range(1,nvars+1))-set(pcode.values()):
        pcode[newgnd]=v
        print 'extend\n',M1.representation(),'\n',pcode
        if testcons(M1,pcode,nsrc,appcns,v,netcons)==1:
            codelist.append(deepcopy(pcode))
        pcode.pop(newgnd)
    return codelist

def certsearch_dfs(M1,pcode,netcons,appcns,d,nsrc,nvars):
    r"""
	Return the lexicographically smallest matroid-network map that forms
    a partial code
    
    INPUT:
    
    ``r`` - rank of identity matrix
    ``q`` - size of finite field
    ``netcons`` - A dictionary specifying network constraints as 
    ``[list1,list2]`` where ``list1`` and ``list2`` are subsets of 
    random variable indices that are forced to have equal rank (entropy)
    ``appcns`` - A dictionary mapping subsets of random variables (as 
    sorted tuples) to subset of constraint labels in ``netcons.keys()``
    ``d`` - depth in the DFS tree, used for recursion
    ``nsrc`` - number of sources in the network
    ``nvars`` - number of network random variables
	"""
    ret=0
    #print 'level d=%d'%d
    if len(pcode.values())>d:
        #print 'first call'
        ret,pcode1=certsearch_dfs(M1,pcode,netcons,appcns,d+1,nsrc,nvars)
    if ret==1:
        return ret,pcode1
    else:
        #print 'explore at level d=%d'%d, sorted(list(M1.groundset()))
        try_g=sorted(list(M1.groundset()))[d]
        #print 'try mapping', try_g
        if try_g in pcode.keys():
            lex_pivot=pcode[try_g]
        else:
            lex_pivot=0
        
        #clean up pcode
        pcode={g:pcode[g] for g in pcode.keys() if g<try_g}
        #print 'clean pcode',pcode
        lex_vars=sorted([v for v in list(set(range(1,nvars+1))-set(pcode.values())) if v>lex_pivot])
        #print 'try vars',lex_vars
        for v in lex_vars:
            #print 'try%d->>%d'%(try_g,v)
            pcode[try_g]=v
            # check if defn works
            if testcons(M1,pcode,nsrc,appcns,v,netcons)==1:
            #def worked
                #print 'pass'
                # are we finished?
                if len(pcode.values())==len(M1.groundset()):
                    return 1,pcode
                else:
                    # go deeper in the tree
                    ret,pcode1=certsearch_dfs(M1,pcode,netcons,appcns,d+1,nsrc,nvars)
                    if ret==1:
                        return 1,pcode1
            else:
                #print 'fail'
                pcode.pop(try_g)
        return 0,pcode


def leftover_certs(M1,pcode,netcons,appcns,d,nsrc,nvars):
    r"""
	returns all matroid network mappings of a given matroid
    
    INPUT
    ``M1`` -- A representabke matroid
    ``pcode`` -- (Default: ``None``) A network-matroid mapping, if 
    specified, only those mappings lexicographically greater than it 
    will be returned
    ``netcons`` -- A dictionary specifying network constraints as 
    ``[list1,list2]`` where ``list1`` and ``list2`` are subsets of 
    random variable indices that are forced to have equal rank (entropy)
    ``appcns`` -- A dictionary mapping subsets of random variables (as 
    sorted tuples) to subset of constraint labels in ``netcons.keys()``
    ``d`` -- depth in the DFS tree, used for recursion
    ``nsrc`` -- number of sources in the network
    ``nvars`` -- number of network random variables
	"""
    if pcode == None:
        pcode = {}
    allpcodelist_ret=[]
    ret=0
    #print 'level d=%d'%d
    if len(pcode.values())>d:
        #print 'first call'
        ret,childpcodelist=certsearch_dfs(M1,pcode,netcons,appcns,d+1,nsrc,nvars)
    if ret==1:
        return ret,pcode1
    else:
        #print 'explore at level d=%d'%d, sorted(list(M1.groundset()))
        try_g=sorted(list(M1.groundset()))[d]
        #print 'try mapping', try_g
        if try_g in pcode.keys():
            lex_pivot=pcode[try_g]
        else:
            lex_pivot=0
        
        #clean up pcode
        pcode={g:pcode[g] for g in pcode.keys() if g<try_g}
        #print 'clean pcode',pcode
        lex_vars=sorted([v for v in list(set(range(1,nvars+1))-set(pcode.values())) if v>lex_pivot])
        #print 'try vars',lex_vars
        for v in lex_vars:
            #print 'try%d->>%d'%(try_g,v)
            pcode[try_g]=v
            # check if defn works
            if testcons(M1,pcode,nsrc,appcns,v,netcons)==1:
            #def worked
                #print 'pass'
                # are we finished?
                if len(pcode.values())==len(M1.groundset()):
                    return 1,pcode
                else:
                    # go deeper in the tree
                    ret,pcode1=certsearch_dfs(M1,pcode,netcons,appcns,d+1,nsrc,nvars)
                    if ret==1:
                        return 1,pcode1
            else:
                #print 'fail'
                pcode.pop(try_g)
        return 0,pcode

def codegen(r,q,netcons,nsrc,nvars):
    r"""
    Returns all valid scalar linear codes obeying ``netcons``
    
    INPUT:
    
    ``r`` -- Rank of matroids to consider
    ``q`` -- Field size
    ``nsrc`` -- Number of sources in the network
    ``nvars`` -- Number of network random variables
    """
    appcns=appcons(netcons)
    M1,codelist=idmatrix_allcodes(r,q,netcons,appcns,nsrc,nvars)
    print 'appcns',appcns
    allpcodes=[[M1,codelist]]
    for i in range(nsrc+1,nvars+1):
        allpcodes_new=[]
        extmats=extend_matroids([code[0] for code in allpcodes])
        for mtr in extmats:
            extcodes=[]
            for pcode in allpcodes[mtr[1]][1]:
                newgnd=list(mtr[0].groundset()-allpcodes[mtr[1]][0].groundset())[0]
                extcodes.extend(extendpcode(mtr[0],pcode,newgnd,netcons,appcns,nsrc,nvars))
            if len(extcodes)>0:
                allpcodes_new.append([mtr[0],extcodes])
        allpcodes=allpcodes_new
    return allpcodes


def codecertgen(r,q,netcons,nsrc,nvars):
    r"""
	Returns a collection of matroids and respective matroid-network 
	mappings
	"""
    appcns=appcons(netcons)
    M1,codecert=idmatrix_codecert(r,q,netcons,appcns,nsrc,nvars)
    allpcodes=[[M1,codecert]]
    print allpcodes
    #print M1
    for i in range(r+1,nvars+1):
        print 'defining element no. %d, with %d matroids'%(i,len(allpcodes))
        allpcodes_new=[]
        extmats=extend_matroids([code[0] for code in allpcodes])
        for mtr in extmats:
            pcode=allpcodes[mtr[1]][1]
            #print 'search for::\n',mtr[0].representation(),pcode
            ret,pcode=certsearch_dfs(mtr[0],pcode,netcons,appcns,0,nsrc,nvars)
            if ret==1:
                allpcodes_new.append([mtr[0],pcode])
        allpcodes=allpcodes_new    
    return allpcodes

def parallel_extensions(Mlist):
    r"""
	Returns all parallel linear extensions of a list of matroids
	"""
    Mset=[]
    for i in range(len(Mlist)):
        M1=Mlist[i]
        newgnd=max(M1.groundset())+1
        extvecs=M1.simplify().representation_vectors()
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


def codecertgen_simple(r,q,netcons,nsrc,nvars,nsimple):
    r"""
	Return a collection of matroids with underlying simple matroid of
	size ``nsimple`` and respective matroid-network mappings
	"""
    appcns=appcons(netcons)
    M1,codecert=idmatrix_codecert(r,q,netcons,appcns,nsrc,nvars)
    allpcodes=[[M1,codecert]]
    print allpcodes
    #print M1
    for i in range(r+1,nsimple+1):
        print 'defining element no. %d, with %d matroids'%(i,len(allpcodes))
        allpcodes_new=[]
        extmats=extend_matroids_simple([code[0] for code in allpcodes])
        for mtr in extmats:
            pcode=allpcodes[mtr[1]][1]
            print 'parent\n',latex(allpcodes[mtr[1]][0].representation()),'search for::\n',latex(mtr[0].representation()),'\n',pcode
            ret,pcode=certsearch_dfs(mtr[0],pcode,netcons,appcns,0,nsrc,nvars)
            if ret==1:
                print 'pass',pcode
                allpcodes_new.append([mtr[0],pcode])
            else:
                print 'fail'
        allpcodes=allpcodes_new
    print 'Done with simple enum'
    simplepcodes=allpcodes
    allpcodes=[]
    print '***********************nonsimple********************'
    for pcode in simplepcodes:
        oldmatlist=[pcode[0]]
        oldpcodelist=[pcode]
        print 'parent',latex(pcode[0].representation()),'\n',pcode[1]
        for i in range(nsimple+1,nvars+1):
            newpcodelist=[]
            parexts=parallel_extensions(oldmatlist)
            for p in parexts:
                print 'parent',latex(oldpcodelist[p[1]][0].representation()),'child','\n',latex(p[0].representation())
                ret,pcode1=certsearch_dfs(p[0],oldpcodelist[p[1]][1],netcons,appcns,0,nsrc,nvars)
                if ret==1:
                    print 'pass',pcode1
                    newpcodelist.append([p[0],pcode1])
                else:
                    print 'fail'
            oldpcodelist=newpcodelist
            oldmatlist=[c[0] for c in newpcodelist]
        allpcodes.extend(newpcodelist)
    return allpcodes
   
def pcodes2rr(allcodes):
    r"""
	Returns a polyhedron corresponding to rate the region
	"""
    ratepts=set([])
    for code in allcodes:
        mtr=code[0]
        pcodes=code[1]
        for pcode in pcodes:
            ipcode={v: k for k, v in pcode.items()}
            print sorted(ipcode.keys())
            ratepts=ratepts|set([tuple([mtr.rank([ipcode[i]]) for i in sorted(ipcode.keys())])])
    return ratepts
    
def pcodecerts2rr(allcodes):
    r"""
	Returns a polyhedron corresponding to rate the region
	"""
    ratepts=set([])
    for code in allcodes:
        mtr=code[0]
        pcode=code[1]
        ipcode={v: k for k, v in pcode.items()}
        print sorted(ipcode.keys())
        ratepts=ratepts|set([tuple([mtr.rank([ipcode[i]]) for i in sorted(ipcode.keys())])])
    return ratepts


def idsc2nsg(mat):
    r"""
	Returns the network symmetry group of an IDSC instance 	 
    
	INPUT:

    - ``mat`` - a matrix specifying IDSC instance.

    OUTPUT:
    
    - The network symmetry group of the IDSC instance
    EXAMPLES::
    """
    G,part=idsc2circgraph(mat)
    Gl=G.line_graph()
    encs=sorted(G.outgoing_edges('S'))
    srcs=sorted(G.incoming_edges('S'))
    dict={}
    j=1
    for i in encs:
        dict[i]=j
        j=j+1
    for i in srcs:
        dict[i]=j
        j=j+1
    for e in set(G.edges())-(set(encs)|set(srcs)):
        dict[e]=e    
    G2=DiGraph()
    for v in dict.keys():
        G2.add_vertex(dict[v])
    for e in Gl.edges():
        G2.add_edge((dict[list(e)[0]],dict[list(e)[1]],None))
    return G2.automorphism_group(partition=[[dict[e] for e  in part],[dict[e] for e in G.edges() if e not in part]])

def idsc2cons(mat):
    r"""
    Return a dictionary containing network constraints.
    INPUT:
    - ``mat`` - A matrix specifyin an IDSC instane.
    
    OUTPUT:
    
    - The network constraints associated with IDSC instance specified by
    ``mat``
    
    EXAMPLES::
    """
    nsrc = ceil(float(log(len(mat))/log(2)))
    x=0
    y=0
    for row in mat:
        for c in row:
            x=x.__or__(Integer(c))
            if Integer(c)!=0:
                y=y+1
    nenc = ceil(float(log(x)/log(2)))
    #print 'nenc',nenc
    ndec = y
    sind = range(1,nsrc+1)
    eind = range(nsrc+1,nsrc+nenc+1)
    dind = range(nsrc+nenc+1,nsrc+nenc+ndec+1)
    dict1={x:[] for x in eind}
    for y in dind:
        dict1[y]=[]
    for s in sind:
        dict1[s]=['S']
    dict1['S']=eind
    dindex=0
    for k in range(len(mat)): 
        row=mat[k]
        #rowwise demands
        rowdem=[i+1 for i in range(len((k+1).binary())) if (k+1).binary()[::-1][i]=='1'] 
        #print rowdem
        for c in row:
            if Integer(c)!=0:
                dict1[dind[dindex]]=rowdem
                ac=Integer(c).binary()[::-1]
                for i in range(len(ac)):
                    if ac[i]=='1':
                        dict1[eind[i]].append(dind[dindex])
                dindex=dindex+1
    consdict={}
    cgraph= DiGraph(dict1)
    j=0
    for e in eind:
        consdict[j]=[set(sind),set(sind)|set([e])]
        j=j+1
    # decoding cons
    for d in dind:
        #print 'dec',d,cgraph.incoming_edges([d]),cgraph.outgoing_edges([d])
        inset=set([list(cgraph.incoming_edges([d])[i])[0] for i in range(len(cgraph.incoming_edges([d])))])
        outset=set([list(cgraph.outgoing_edges([d])[i])[1] for i in range(len(cgraph.outgoing_edges([d])))])
        consdict[j]=[inset,outset|inset]
        #print consdict[j]
        j=j+1
    return consdict

def idsc2circgraph(mat):
    r"""
	Return the circulation graph of an ISDC instance
    
    INPUT:

    - ``mat`` - A matrix specifyin an IDSC instane.
    
    OUTPUT:
    
    A graph encoding the symmetries of the IDSC instance specified by 
    ``mat``
    """
    nsrc = ceil(float(log(len(mat))/log(2)))
    alph='abcdefghijklmnopqrstuvwxyz'
    x=0
    y=0
    for row in mat:
        for c in row:
            x=x.__or__(Integer(c))
            if Integer(c)!=0:
                y=y+1
    nenc = ceil(float(log(x)/log(2)))
    ndec = y
    sind = range(1,nsrc+1)
    eind = [alph[i] for i in range(nenc)]
    dind = [alph[i] for i in range(nenc,ndec+nenc)]
    dict={x:[] for x in eind}
    for y in dind:
        dict[y]=[]
    for s in sind:
        dict[s]=['S']
    dict['S']=eind
    dindex=0
    for k in range(len(mat)):
        row=mat[k]
        rowdem=[i+1 for i in range(len((k+1).binary())) if (k+1).binary()[::-1][i]=='1']
        for c in row:
            if Integer(c)!=0:
                dict[dind[dindex]]=rowdem
                ac=Integer(c).binary()[::-1]
                for i in range(len(ac)):
                    if ac[i]=='1':
                        dict[eind[i]].append(dind[dindex])
                dindex=dindex+1
    return DiGraph(dict)

