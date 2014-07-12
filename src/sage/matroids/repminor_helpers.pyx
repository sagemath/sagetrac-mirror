include 'sage/misc/bitset.pxi'
from lean_matrix cimport LeanMatrix, GenericMatrix, BinaryMatrix, TernaryMatrix, QuaternaryMatrix, IntegerMatrix, generic_identity
from sage.rings.all import ZZ, FiniteField, GF

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

cpdef init_iso_matrices(BinaryMatrix M_rmat,BinaryMatrix N_rmat):
    """
    Return a list of two matrices [M1_0, M2_0] s.t. [M1_0 0; 0 M2_0] will be
    the bitwise logical OR of all degree preserving isomorphisms  
    """
    cdef long r1,r2
    cdef BinaryMatrix M1_0,M2_0
    # initialize M1_0 and M2_0
    M1_0 = BinaryMatrix(N_rmat.nrows(),M_rmat.nrows(), ring=GF(2))
    M2_0 = BinaryMatrix(N_rmat.ncols(),M_rmat.ncols(),ring=GF(2))
    # for lhs vertices in N_rmat(rows), find out rows in M_rmat that have
    # greater than equal degree
    for r1 in range(N_rmat.nrows()):
        d1 = N_rmat.row_len(r1)
        for r2 in range(M_rmat.nrows()):
            d2 = M_rmat.row_len(r2)
            if d1 <= d2:
                M1_0.set(r1,r2)
    for c1 in range(N_rmat.ncols()):
        d1 = col_len(N_rmat,c1)
        for c2 in range(M_rmat.ncols()):
            d2 = col_len(M_rmat,c2)
            if d1 <= d2:
                M2_0.set(c1,c2)
    return [M1_0,M2_0]
    
cpdef col_len(BinaryMatrix mat,long c):
    cdef int r,d
    d=0
    for r in range(mat.nrows()):
        if mat.is_nonzero(r,c):
            d+=1
    return d

cpdef _check_bin_minor(BinaryMatrix M_rmat,BinaryMatrix N_rmat):
    cdef long i
    cdef bitset_t unused_cols
    cdef BinaryMatrix M_rmatT, N_rmatT,M1,M2
    M_rmatT = M_rmat.transpose()
    N_rmatT=N_rmat.transpose()
    bitset_init(unused_cols,M_rmat.nrows()+ M_rmat.ncols())
    [M1_0,M2_0]=init_iso_matrices(M_rmat,N_rmat)
    for i in range(M_rmat.nrows()+M_rmat.ncols()):
        bitset_add(unused_cols,i)
    print bitset_string(unused_cols)
    M1=BinaryMatrix(M1_0.nrows(),M1_0.ncols(),ring=GF(2))
    M2=BinaryMatrix(M2_0.nrows(),M2_0.ncols(),ring=GF(2))
    recurse(unused_cols,0,M_rmat,N_rmat,M_rmatT,N_rmatT,M1,M2,M1_0,M2_0)

cdef recurse(bitset_t unused_cols,long cur_row, BinaryMatrix M_rmat, BinaryMatrix N_rmat,BinaryMatrix M_rmatT, BinaryMatrix N_rmatT,BinaryMatrix M1,BinaryMatrix M2, BinaryMatrix M1_0, BinaryMatrix M2_0):
    cdef long i,j
    cdef bitset_t valid_cols
    cpdef bitset_t curr_cols
    if cur_row == M1.nrows() + M2.nrows():
        # check if [M1 0; 0 M2] is an induced isomorphism and return True/False accordingly
        pass
    elif cur_row < M1.nrows():
        bitset_init(curr_cols,M_rmat.nrows()+M_rmat.ncols())
        bitset_init(valid_cols,M_rmat.nrows()+M_rmat.ncols())
        for i in range(M1.ncols()):
            bitset_add(curr_cols,i)
        bitset_intersection(valid_cols, curr_cols, M1_0._M[cur_row])
    else:
        bitset_init(curr_cols,M_rmat.nrows()+M_rmat.ncols())
        bitset_init(valid_cols,M_rmat.nrows()+M_rmat.ncols())
        for i in range(M1.ncols(),M1.ncols()+M2.ncols()):
            bitset_add(curr_cols,i)
        bitset_intersection(valid_cols, curr_cols, M2_0._M[cur_row-M1_0.nrows()])
    print bitset_string(curr_cols)
    print bitset_string(valid_cols)
    while True:
        try:
            i=pop(valid_cols)
            if cur_row < M1.nrows():
                # set up M1, prune and recurse
                bitset_clear(M1._M[cur_row])
                bitset_add(M1._M[cur_row],i)
                recurse(unused_cols,0,M_rmat,N_rmat,M_rmatT,N_rmatT,M1,M2,M1_0,M2_0)
            else:
                # set up M2, prune and recurse
                recurse(unused_cols,0,M_rmat,N_rmat,M_rmatT,N_rmatT,M1,M2,M1_0,M2_0)
    #M1=prune(M1_0,M2_0,cur_row, M_rmat,N_rmat,M_rmatT,N_rmatT)
    print 'returned M1',M1

cpdef prune(BinaryMatrix M1, BinaryMatrix M2, long cur_row,BinaryMatrix M_rmat, BinaryMatrix N_rmat,BinaryMatrix M_rmatT, BinaryMatrix N_rmatT):
    cdef long r,c,flag
    print 'M1',M1
    print 'M2',M2
    r=cur_row
    for c in range(M1.ncols()):
        if M1.is_nonzero(r,c):
            print r,'mapsto',c,'neighbours1', _neighbours(N_rmat,N_rmatT,r,1)  
            for x in _neighbours(N_rmat,N_rmatT,r,1):
                flag = 0
                print 'neighbours2', _neighbours(M_rmat,M_rmatT,c,1) 
                for y in _neighbours(M_rmat,M_rmatT,c,1):
                    if M2.is_nonzero(x,y):
                        flag = 1
                        break
                if flag == 0:
                    print 'bad neighbour',r,c
                    M1.set_unsafe(r, c, 0)
    return M1

cpdef copy_mat(BinaryMatrix mat):
    cdef long i
    cdef BinaryMatrix mat_copy
    mat_copy = BinaryMatrix(mat.nrows(),mat.ncols(),ring=GF(2))
    for i in range(mat.nrows()):
        bitset_copy(mat_copy._M[i], mat._M[i])
    return mat_copy

cpdef _neighbours(BinaryMatrix rmat, BinaryMatrix rmatT, long i, long rc):
    cdef long c
    neighbours=[]
    if rc == 1: #neighbors of lhs , look through ith row
        for c in range(rmat.ncols()):
            if rmat.is_nonzero(i, c):
                neighbours.append(c)
    if rc == 0:
        for c in range(rmatT.nrows()):
            if rmatT.is_nonzero(i,c):
                neighbours.append(c)
    return neighbours
    
cdef full_set(bitset_t bits):
    "return a set with 0:len-1 elements"
    cdef i
    for  i in range(bitset_len(bits)):
        bitset_add(bits,i)
