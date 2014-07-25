# cython: profile=True

include 'sage/misc/bitset.pxi'
from lean_matrix cimport BinaryMatrix
from sage.rings.all import GF
cdef bint GF2_not_defined = True
cdef GF2, GF2_one, GF2_zero
cdef bint npruned_not_defined = True
cdef npruned


cpdef search_assgn(Ga,Gb,assgn,cand=None):
    """ Implement Ullman's algorithm for exact subgraph isomorphism 
    """
    uGa_S1 = set(Ga['S1']) - set(assgn.keys())
    uGb_S1= set(Gb['S1']) - set(assgn.values())
    uGa_S2 = set(Ga['S2']) - set(assgn.keys())
    uGb_S2= set(Gb['S2']) - set(assgn.values())
    uGa = uGa_S1 | uGa_S2
    uGb = uGb_S1 | uGb_S2
    print 'Ga',Ga
    print 'Gb',Gb
    if cand == None and assgn =={}:
        # first call. populate cand
        cand = {}
        for va in set(Ga['S1']):
            cand[va]=[]
            for vb in set(Gb['S1']):
                print 'Ga:',va,'Gb:',vb,Ga[va][0],'?',Gb[vb][0]
                if Gb[vb][0] >= Ga[va][0]:
                    cand[va].append(vb)
        for va in set(Ga['S2']):
            cand[va]=[]
            for vb in set(Gb['S2']):
                print 'Ga:',va,'Gb:',vb,Ga[va][0],'?',Gb[vb][0]
                if Gb[vb][0] >= Ga[va][0]:
                    cand[va].append(vb)
    print cand
    # loop over uGa
    for v in uGa:
        # check that whenever vertex c of Gb is among candidates for vertex v of 
        # Ga, then every unmapped neighbor of v has at least one candidate among 
        # unmapped neighbors of c.
        u_vn = [vn for vn in Ga[v][1] if vn in uGa] # unmapped neighbours of v
        u_vn_c = [cand[x] for x in u_vn] # candidates of unmapped neighbours of v
        for c in cand[v]:
            badcand = []
            # check if this candidate is valid 
            u_cn =  [cn for cn in uGb if cn in Gb[c][1]] # unmapped neighbours of candidate
            for ucn in u_cn:
                if not any([ucn in uvnc for uvnc in u_vn_c]):
                    # add c to badcand
                    badcand.append(c)
                    
                    
    return True

cpdef findmap(M, N):
    [a,X,Y]=M.has_minor(N)
    N1=M.contract(X).delete(Y)
    NB=N.nonbases()
    return NB._isomorphism(N1.nonbases())
    

cpdef init_iso_matrices(BinaryMatrix M_rmat,BinaryMatrix N_rmat):
    """
    Return a list of two matrices [M1_0, M2_0] s.t. [M1_0 0; 0 M2_0] will be
    the bitwise logical OR of all degree preserving isomorphisms  
    """
    cdef long r1,r2
    cdef BinaryMatrix M1_0,M2_0
    global GF2
    # initialize M1_0 and M2_0
    M1_0 = BinaryMatrix(N_rmat.nrows(),M_rmat.nrows(), ring=GF2)
    M2_0 = BinaryMatrix(N_rmat.ncols(),M_rmat.ncols(),ring=GF2)
    # for lhs vertices in N_rmat(rows), find out rows in M_rmat that have
    # greater than equal degree
    for r1 in xrange(N_rmat.nrows()):
        d1 = N_rmat.row_len(r1)
        for r2 in xrange(M_rmat.nrows()):
            d2 = M_rmat.row_len(r2)
            if d1 <= d2:
                M1_0.set(r1,r2)
    for c1 in xrange(N_rmat.ncols()):
        d1 = col_len(N_rmat,c1)
        for c2 in xrange(M_rmat.ncols()):
            d2 = col_len(M_rmat,c2)
            if d1 <= d2:
                M2_0.set(c1,c2)
    return [M1_0,M2_0]
    
cpdef col_len(BinaryMatrix mat,long c):
    cdef int r,d
    d=0
    for r in xrange(mat.nrows()):
        if mat.is_nonzero(r,c):
            d+=1
    return d

cpdef _check_bin_minor(BinaryMatrix M_rmat,BinaryMatrix N_rmat,indices,M,Npcl, long nloops):
    cdef long i,x
    cdef bitset_t unused_cols
    cdef BinaryMatrix M_rmatT, N_rmatT,M1,M2
    global GF2
    M_rmatT = mat_transpose(M_rmat)
    N_rmatT = mat_transpose(N_rmat)
    #print M_rmat
    #print N_rmat
    bitset_init(unused_cols,M_rmat.nrows()+ M_rmat.ncols())
    [M1_0,M2_0]=init_iso_matrices(M_rmat,N_rmat)
    #print 'M1_0',M1_0,'M2_0',M2_0
    for i in xrange(M_rmat.nrows()+M_rmat.ncols()):
        bitset_add(unused_cols,i)
    M1=copy_mat(M1_0) #BinaryMatrix(M1_0.nrows(),M1_0.ncols(),ring=GF(2))
    M2=BinaryMatrix(M2_0.nrows(),M1_0.ncols(),ring=GF2).augment(copy_mat(M2_0))#BinaryMatrix(M2_0.nrows(),M2_0.ncols(),ring=GF(2))
    if not degrees_are_sane(M1,M2):
        #print 'degrees insane'
        return False
    else:
        #print 'degrees sane'
        #print 'M1',M1,'\n','M2',M2
        x= recurse(unused_cols,0,M_rmat,N_rmat,M_rmatT,N_rmatT,M1,M2,indices, M,Npcl,nloops) 
        #print 'npruned', npruned
        return x== 1

cpdef degrees_are_sane(BinaryMatrix M1,BinaryMatrix M2):
    """ Check if there is an all-zero row in M1 or M2 and return ``False``
    if there is one. 
    """
    cdef long i
    for i in xrange(M1.nrows()):
        if bitset_isempty(M1._M[i]):
            return False
    for i in xrange(M2.nrows()):
        if bitset_isempty(M2._M[i]):
            return False
    return True

cpdef _lp_are_sane(M, N):
    # check if loops of N have consistent mappings
    if len(M.loops()) < len(N.loops()):
        print 'loops insane'
        return False
    # check if there are atleast as many nontrivial parallel class in M as 
    # there are in N
    Mpcl=[k-M.loops() for k in M.flats(1) if len(k-M.loops())>1]
    print 'Mpcl',Mpcl
    Npcl=[k-N.loops() for k in N.flats(1) if len(k-N.loops())>1]
    print 'Npcl',Npcl
    if len(Mpcl) < len(Npcl):
        return False
    for npcl in Npcl:
        if len([mpcl for mpcl in Mpcl if len(mpcl)>=len(npcl)]) == 0:
            return False
    return True

cdef _pcl_have_maps(M, Npcl, iso_map, indices,long nloops):
    # contract the rows that are unmapped
#    print 'indices',indices
#    print 'allrows',indices[4]
#    print 'badrows', set(iso_map.values())
#    print 'contract', set(indices[4])-set(iso_map.values())
    M2=M.contract(set(indices[4])-set(iso_map.values()))
    if len(M2.loops()) < nloops:
        return False
    # for each repmap check if it has big enough parallel class
    for pcl in Npcl:
        rep = list(set(pcl) & set(iso_map.keys()))[0] # representative
        repmap = iso_map[rep]
        repmap_pcl_M = set(M2.closure({repmap}))-M2.loops()
        if len(repmap_pcl_M) < len(pcl):
            # pcl doesn't have consistant map
            return False
    return True
        
cdef iso_mats_to_dict(BinaryMatrix M1, BinaryMatrix M2, indices):
    M_R=indices[0]
    M_C=indices[1]
    N_R=indices[2]
    N_C=indices[3]
    iso_map={}
    for i in xrange(M1.nrows()):
        iso_map[N_R[i]]=M_R[bitset_list(M1._M[i])[0]]
    for i in xrange(M2.nrows()):
        iso_map[N_C[i]]=M_C[bitset_list(M2._M[i])[0]-M1.ncols()]
    return iso_map
                
        
cdef recurse(bitset_t unused_cols,long cur_row, BinaryMatrix M_rmat, BinaryMatrix N_rmat,BinaryMatrix M_rmatT, BinaryMatrix N_rmatT,BinaryMatrix M1,BinaryMatrix M2,indices,M,Npcl,long nloops):
    cdef long i,j
    cpdef bitset_t valid_cols
    cpdef bitset_t curr_cols
    cpdef bitset_t saved_row
    cdef BinaryMatrix iso,saved_M1,saved_M2
    global GF2, GF2_zero, GF2_one, GF2_not_defined
    if GF2_not_defined:
            GF2 = GF(2)
            GF2_zero = GF2(0)
            GF2_one = GF2(1)
            GF2_not_defined = False
    cdef long found = 0
    #print 'cur_row',cur_row
    if cur_row == M1.nrows() + M2.nrows():
        # check if [M1 0; 0 M2] is an induced isomorphism and return True/False accordingly
        #print 'check induced isomorphism'
#        iso= M1.augment(BinaryMatrix(M1.nrows(),M2.ncols()-M1.ncols(),ring=GF2)).stack(M2)
#        print iso
        found = is_weak_induced_isomorphism( M1, M2,   M_rmat,  N_rmat)
        #print 'wii=',found
        # return found
        if found == 1:
            #print 'weak isomorphism'
#            print 'M1',M1
#            print 'M2',M2
            iso_map = iso_mats_to_dict(M1,M2,indices)
            #print iso_map
            
            if _pcl_have_maps(M,Npcl,iso_map,indices,nloops):
                return found
            else:
                return 0
        else:
            return 0
    else:
        # save M1,M2 before pruning
        saved_M1 = copy_mat(M1)
        saved_M2 = copy_mat(M2)
        if prune(cur_row,M1,M2, M_rmat,N_rmat,M_rmatT,N_rmatT) is False:
#            restore_mat(M1,saved_M1)
#            restore_mat(M2,saved_M2)
            if cur_row <M1.nrows():
                restore_mat(M1,saved_M1,cur_row)
                restore_mat(M2,saved_M2,0)
            else:
                restore_mat(M2,saved_M2,cur_row-M1.nrows())
#            if cur_row==2:
#                print M_rmat, N_rmat
#                return 1
#            else:
            return found 
        if cur_row < M1.nrows():
            # cur_row in in M1
            bitset_init(curr_cols,M_rmat.nrows()+M_rmat.ncols())
            bitset_init(valid_cols,M_rmat.nrows()+M_rmat.ncols())
            for i in xrange(M1.ncols()):
                bitset_add(curr_cols,i)
            bitset_intersection(valid_cols, curr_cols, M1._M[cur_row])
            bitset_intersection(valid_cols,valid_cols,unused_cols)
            bitset_init(saved_row,M1.ncols())
            bitset_copy(saved_row,M1._M[cur_row])
        else:
            # cur_row is in M2bitset_pop
            bitset_init(curr_cols,M_rmat.nrows()+M_rmat.ncols())
            bitset_init(valid_cols,M_rmat.nrows()+M_rmat.ncols())
            for i in xrange(M1.ncols(),M1.ncols()+M2.ncols()):
                bitset_add(curr_cols,i)
            bitset_intersection(valid_cols, curr_cols, M2._M[cur_row-M1.nrows()])
            bitset_intersection(valid_cols,valid_cols,unused_cols)
            bitset_init(saved_row,M2.ncols())
            bitset_copy(saved_row,M2._M[cur_row-M1.nrows()])
    iso= M1.augment(BinaryMatrix(M1.nrows(),M2.ncols()-M1.ncols(),ring=GF2)).stack(M2)
    #print iso
    while not bitset_isempty(valid_cols):
        i=bitset_pop(valid_cols)
        bitset_remove(unused_cols,i)
#        if cur_row==0:
#            print cur_row,'->', i
#        if M1.is_nonzero(0,1) :
#            print cur_row,'->', i
#        print iso
#        print cur_row,'->', i
        if cur_row < M1.nrows():
            # set up M1, prune and recurse
            bitset_clear(M1._M[cur_row])
            bitset_add(M1._M[cur_row],i)
            found=recurse(unused_cols,cur_row+1,M_rmat,N_rmat,M_rmatT,N_rmatT,M1,M2,indices,M,Npcl,nloops)
            bitset_add(unused_cols,i)
            if found == 1:
                break
        else:
            # set up M2, prune and recurse
            bitset_clear(M2._M[cur_row-M1.nrows()])
            bitset_add(M2._M[cur_row-M1.nrows()],i)
            found=recurse(unused_cols,cur_row+1,M_rmat,N_rmat,M_rmatT,N_rmatT,M1,M2,indices,M,Npcl,nloops)
            bitset_add(unused_cols,i)
            if found == 1:
                break
    if cur_row <M1.nrows():
        restore_mat(M1,saved_M1,cur_row)
        restore_mat(M2,saved_M2,0)
    else:
        restore_mat(M2,saved_M2,cur_row-M1.nrows())
    return found

cdef restore_mat(BinaryMatrix mat, BinaryMatrix saved_mat,long cur_row):
    cdef long r
    for r in xrange(cur_row,saved_mat.nrows()):
        bitset_copy(mat._M[r],saved_mat._M[r])
    return

cpdef prune(long cur_row,BinaryMatrix M1, BinaryMatrix M2, BinaryMatrix M_rmat, BinaryMatrix N_rmat,BinaryMatrix M_rmatT, BinaryMatrix N_rmatT):
    global npruned,npruned_not_defined
    if npruned_not_defined:
        npruned = 0
        npruned_not_defined=False
    cdef long r,c,flag
    #cdef long npruned = 0
    for r in xrange(max(0,cur_row),M1.nrows()):
        for c in xrange(M1.ncols()):
            if M1.is_nonzero(r,c):
                #print 'check',r,'->>',c
                #print 'n(r)',_neighbours(N_rmat,N_rmatT,r,1)
                for x in _neighbours(N_rmat,N_rmatT,r,1):
                    flag = 0
                    #print 'n(c)',_neighbours(M_rmat,M_rmatT,c,1)
                    for y in _neighbours(M_rmat,M_rmatT,c,1):
                        if M2.is_nonzero(x,y+M1.ncols()):
                            flag = 1
                            break
                    if flag == 0:
                        npruned+=1
                        #print 'bad mapping',r,'->>',c
                        M1.set_unsafe(r, c, 0)
    for r in xrange(max(0,cur_row),M2.nrows()):
        for c in xrange(M1.ncols(),M2.ncols()):
            if M2.is_nonzero(r,c):
                #print 'check',r+M1.nrows(),'->>',c
                #print 'n(r)', _neighbours(N_rmat,N_rmatT,r,0)
                for x in _neighbours(N_rmat,N_rmatT,r,0):
                    flag = 0
                    #print 'n(c)', _neighbours(M_rmat,M_rmatT,c-M1.ncols(),0)
                    for y in _neighbours(M_rmat,M_rmatT,c-M1.ncols(),0):
                        if M1.is_nonzero(x,y):
                            flag = 1
                            break
                    if flag == 0:
                        npruned+=1
                        #print 'bad mapping',r+M1.nrows(),'->>',c
                        M2.set_unsafe(r, c, 0)
    # check degree sanity
    #print 'depth',cur_row,'npruned', npruned
    if degrees_are_sane(M1,M2):
        return True
    else:
        #print 'pruning creates degree insanity'
        return False

cpdef BinaryMatrix copy_mat(BinaryMatrix mat):
    cdef long i
    global GF2
    cdef BinaryMatrix mat_copy
    mat_copy = BinaryMatrix(mat.nrows(),mat.ncols(),ring=GF2)
    for i in xrange(mat.nrows()):
        bitset_copy(mat_copy._M[i], mat._M[i])
    return mat_copy

cpdef _neighbours(BinaryMatrix rmat, BinaryMatrix rmatT, long i, long rc):
    cdef long c
    neighbours=[]
    if rc == 1: #neighbors of lhs , look through ith row
        for c in xrange(rmat.ncols()):
            if rmat.is_nonzero(i, c):
                neighbours.append(c)
    if rc == 0:
        for c in xrange(rmatT.ncols()):
            if rmatT.is_nonzero(i,c):
                neighbours.append(c)
    return neighbours
    
cdef full_set(bitset_t bits):
    "return a set with 0:len-1 elements"
    cdef i
    for  i in xrange(bitset_len(bits)):
        bitset_add(bits,i)

cpdef is_weak_induced_isomorphism(BinaryMatrix M1_1,BinaryMatrix M2_1, BinaryMatrix  M_rmat, BinaryMatrix N_rmat):
    """
    For binary matrices this will in fact answer if M1,M2 specify induced 
    subgraph isomorphism whereas for higher fields, this only verifies a 
    necessary but not sufficient condition 
    """
    cdef long i,k
    cdef int j,*roworders, *colorders
    cdef BinaryMatrix rmat1
    cdef BinaryMatrix M1,M2
    global GF2
    M1 = M1_1.copy()
    M2=M2_1.copy()
    # rmat corresponding to given isomorphism
    rmat1=BinaryMatrix(N_rmat.nrows(),N_rmat.ncols(),ring=GF2)
    roworders = <int* > sage_malloc(N_rmat.nrows() * sizeof(int))
    for i in xrange(M1.nrows()):
        j = bitset_pop(M1._M[i])
        # row i maps to row j
        roworders[i]=j
    colorders = <int* > sage_malloc(N_rmat.ncols() * sizeof(int))
    for i in xrange(M2.nrows()):
        j = bitset_pop(M2._M[i])
        # col i maps to row j
        colorders[i]=j-M1.ncols()
#    for i in range(M1.nrows()):
#        print roworders[i]
#    for i in range(M2.nrows()):
#        print colorders[i]
    for i in xrange(M1.nrows()):
        for k in xrange(M2.nrows()):
            #print roworders[i],',',colorders[k]
            rmat1.set_unsafe(i,k,M_rmat.get_unsafe(roworders[i],colorders[k]))
    sage_free(roworders)
    sage_free(colorders)
    for i in xrange(N_rmat.nrows()):
        if not bitset_eq(N_rmat._M[i], rmat1._M[i]):
            return 0
    return 1
    
cdef mat_transpose(BinaryMatrix mat):
    cdef long i,j
    matT = BinaryMatrix(mat.ncols(),mat.nrows()) 
    for i in xrange(mat.nrows()):
        for j in xrange(mat.ncols()):
            matT.set_unsafe(j,i,mat.get_unsafe(i,j))
    return matT

cpdef bint is_new_rmat(BinaryMatrix cand, used_rmats):
    cdef long i
    for i in xrange(len(used_rmats)):
        if mats_equal(cand,used_rmats[i]):
            return False
    return True

    
cpdef bint mats_equal(BinaryMatrix M1, BinaryMatrix M2):
    cdef long i
    for i in xrange(M1.nrows()):
        if not bitset_eq(M1._M[i],M2._M[i]):
            return False
    return True
    
