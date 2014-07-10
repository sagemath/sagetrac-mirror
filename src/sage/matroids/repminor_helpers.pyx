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
