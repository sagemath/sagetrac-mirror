# cython: profile = True

include 'sage/misc/bitset.pxi'
from lean_matrix cimport BinaryMatrix
from sage.rings.all import GF
cdef bint npruned_not_defined = True
cdef npruned

cpdef init_iso_matrices(BinaryMatrix M_rmat, BinaryMatrix N_rmat):
    """
    Return a list of two matrices [M1_0, M2_0] s.t. [M1_0 0; 0 M2_0] will be
    the bitwise logical OR of all degree preserving isomorphisms
    """
    cdef long r1, r2
    cdef BinaryMatrix M1_0, M2_0
    # initialize M1_0 and M2_0
    M1_0 = BinaryMatrix(N_rmat.nrows(), M_rmat.nrows())
    M2_0 = BinaryMatrix(N_rmat.ncols(), M_rmat.ncols())
    for r1 in xrange(N_rmat.nrows()):
        d1 = N_rmat.row_len(r1)
        for r2 in xrange(M_rmat.nrows()):
            d2 = M_rmat.row_len(r2)
            if d1 <= d2:
                M1_0.set(r1, r2)
    for c1 in xrange(N_rmat.ncols()):
        d1 = col_len(N_rmat, c1)
        for c2 in xrange(M_rmat.ncols()):
            d2 = col_len(M_rmat, c2)
            if d1 <= d2:
                M2_0.set(c1, c2)
    return [M1_0, M2_0]

cpdef col_len(BinaryMatrix mat, long c):
    cdef int r, d
    d = 0
    for r in xrange(mat.nrows()):
        if mat.is_nonzero(r, c):
            d += 1
    return d

cpdef _check_bin_minor(BinaryMatrix M_rmat, BinaryMatrix N_rmat, indices, M,
                       Npcl, long nloops, bint debug=False):
    cdef long i, x
    cdef bitset_t unused_cols
    cdef BinaryMatrix M_rmatT, N_rmatT, M1, M2
    M_rmatT = mat_transpose(M_rmat)
    N_rmatT = mat_transpose(N_rmat)
    bitset_init(unused_cols, M_rmat.nrows() + M_rmat.ncols())
    [M1_0, M2_0] = init_iso_matrices(M_rmat, N_rmat)
    for i in xrange(M_rmat.nrows()+M_rmat.ncols()):
        bitset_add(unused_cols, i)
    M1 = copy_mat(M1_0)
    M2 = BinaryMatrix(M2_0.nrows(), M1_0.ncols()).augment(copy_mat(M2_0))
    if not degrees_are_sane(M1, M2):
        return False
    else:
        x = recurse(unused_cols, 0, M_rmat, N_rmat, M_rmatT, N_rmatT, M1, M2,
                    indices, M, Npcl, nloops)
        return x == 1

cpdef degrees_are_sane(BinaryMatrix M1, BinaryMatrix M2):
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


cdef _pcl_have_maps(M, Npcl, iso_map, indices, long nloops):
    # contract the rows that are unmapped
    M2 = M.contract(set(indices[4]) - set(iso_map.values()))
    if len(M2.loops()) < nloops:
        return False
    # for each repmap check if it has big enough parallel class
    for pcl in Npcl:
        rep = list(set(pcl) & set(iso_map.keys()))[0]  # representative
        repmap = iso_map[rep]
        repmap_pcl_M = set(M2.closure({repmap})) - M2.loops()
        if len(repmap_pcl_M) < len(pcl):
            # pcl doesn't have consistant map
            return False
    return True

cdef iso_mats_to_dict(BinaryMatrix M1, BinaryMatrix M2, indices):
    M_R = indices[0]
    M_C = indices[1]
    N_R = indices[2]
    N_C = indices[3]
    iso_map = {}
    for i in xrange(M1.nrows()):
        iso_map[N_R[i]] = M_R[bitset_list(M1._M[i])[0]]
    for i in xrange(M2.nrows()):
        iso_map[N_C[i]] = M_C[bitset_list(M2._M[i])[0] - M1.ncols()]
    return iso_map


cdef recurse(bitset_t unused_cols, long cur_row, BinaryMatrix M_rmat,
             BinaryMatrix N_rmat, BinaryMatrix M_rmatT, BinaryMatrix N_rmatT,
             BinaryMatrix M1, BinaryMatrix M2, indices, M, Npcl, long nloops):
    cdef long i, j
    cpdef bitset_t valid_cols
    cpdef bitset_t curr_cols
    cpdef bitset_t saved_row
    cdef BinaryMatrix iso, saved_M1, saved_M2
    cdef long found = 0
    if cur_row == M1.nrows() + M2.nrows():
        found = is_weak_induced_isomorphism(M1, M2, M_rmat, N_rmat)
        if found == 1:
            iso_map = iso_mats_to_dict(M1, M2, indices)
            if _pcl_have_maps(M, Npcl, iso_map, indices, nloops):
                return found
            else:
                return 0
        else:
            return 0
    else:
        # save M1, M2 before pruning
        saved_M1 = copy_mat(M1)
        saved_M2 = copy_mat(M2)
        if prune(cur_row, M1, M2, M_rmat, N_rmat, M_rmatT, N_rmatT) is False:
            if cur_row < M1.nrows():
                restore_mat(M1, saved_M1, cur_row)
                restore_mat(M2, saved_M2, 0)
            else:
                restore_mat(M2, saved_M2, cur_row - M1.nrows())
            return found
        if cur_row < M1.nrows():
            # cur_row in in M1
            bitset_init(curr_cols, M_rmat.nrows() + M_rmat.ncols())
            bitset_init(valid_cols, M_rmat.nrows() + M_rmat.ncols())
            for i in xrange(M1.ncols()):
                bitset_add(curr_cols, i)
            bitset_intersection(valid_cols, curr_cols, M1._M[cur_row])
            bitset_intersection(valid_cols, valid_cols, unused_cols)
            bitset_init(saved_row, M1.ncols())
            bitset_copy(saved_row, M1._M[cur_row])
        else:
            # cur_row is in M2
            bitset_init(curr_cols, M_rmat.nrows()+M_rmat.ncols())
            bitset_init(valid_cols, M_rmat.nrows() + M_rmat.ncols())
            for i in xrange(M1.ncols(), M1.ncols() + M2.ncols()):
                bitset_add(curr_cols, i)
            bitset_intersection(valid_cols, curr_cols,
                                M2._M[cur_row - M1.nrows()])
            bitset_intersection(valid_cols, valid_cols, unused_cols)
            bitset_init(saved_row, M2.ncols())
            bitset_copy(saved_row, M2._M[cur_row - M1.nrows()])
    while not bitset_isempty(valid_cols):
        i = bitset_pop(valid_cols)
        bitset_remove(unused_cols, i)
        if cur_row < M1.nrows():
            # set up M1 and recurse
            bitset_clear(M1._M[cur_row])
            bitset_add(M1._M[cur_row], i)
            found = recurse(unused_cols, cur_row + 1, M_rmat, N_rmat, M_rmatT,
                            N_rmatT, M1, M2, indices, M, Npcl, nloops)
            bitset_add(unused_cols, i)
            if found == 1:
                break
        else:
            # set up M2 and recurse
            bitset_clear(M2._M[cur_row - M1.nrows()])
            bitset_add(M2._M[cur_row - M1.nrows()], i)
            found = recurse(unused_cols, cur_row + 1, M_rmat, N_rmat, M_rmatT,
                            N_rmatT, M1, M2, indices, M, Npcl, nloops)
            bitset_add(unused_cols, i)
            if found == 1:
                break
    if cur_row < M1.nrows():
        restore_mat(M1, saved_M1, cur_row)
        restore_mat(M2, saved_M2, 0)
    else:
        restore_mat(M2, saved_M2, cur_row - M1.nrows())
    return found

cdef restore_mat(BinaryMatrix mat, BinaryMatrix saved_mat, long cur_row):
    cdef long r
    for r in xrange(cur_row, saved_mat.nrows()):
        bitset_copy(mat._M[r], saved_mat._M[r])
    return

cpdef prune(long cur_row, BinaryMatrix M1, BinaryMatrix M2,
            BinaryMatrix M_rmat, BinaryMatrix N_rmat, BinaryMatrix M_rmatT,
            BinaryMatrix N_rmatT):
    global npruned, npruned_not_defined
    if npruned_not_defined:
        npruned = 0
        npruned_not_defined = False
    cdef long r, c, flag
    for r in xrange(max(0, cur_row), M1.nrows()):
        for c in xrange(M1.ncols()):
            if M1.is_nonzero(r, c):
                for x in _neighbours(N_rmat, N_rmatT, r, 1):
                    flag = 0
                    for y in _neighbours(M_rmat, M_rmatT, c, 1):
                        if M2.is_nonzero(x, y + M1.ncols()):
                            flag = 1
                            break
                    if flag == 0:
                        npruned += 1
                        M1.set_unsafe(r, c, 0)
    for r in xrange(max(0, cur_row), M2.nrows()):
        for c in xrange(M1.ncols(), M2.ncols()):
            if M2.is_nonzero(r, c):
                for x in _neighbours(N_rmat, N_rmatT, r, 0):
                    flag = 0
                    for y in _neighbours(M_rmat, M_rmatT, c - M1.ncols(), 0):
                        if M1.is_nonzero(x, y):
                            flag = 1
                            break
                    if flag == 0:
                        npruned += 1
                        M2.set_unsafe(r, c, 0)
    # check degree sanity
    if degrees_are_sane(M1, M2):
        return True
    else:
        return False

cpdef BinaryMatrix copy_mat(BinaryMatrix mat):
    cdef long i
    cdef BinaryMatrix mat_copy
    mat_copy = BinaryMatrix(mat.nrows(), mat.ncols())
    for i in xrange(mat.nrows()):
        bitset_copy(mat_copy._M[i], mat._M[i])
    return mat_copy

cpdef _neighbours(BinaryMatrix rmat, BinaryMatrix rmatT, long i, long rc):
    cdef long c
    neighbours = []
    if rc == 1:  # neighbors of lhs , look through ith row
        for c in xrange(rmat.ncols()):
            if rmat.is_nonzero(i, c):
                neighbours.append(c)
    if rc == 0:
        for c in xrange(rmatT.ncols()):
            if rmatT.is_nonzero(i, c):
                neighbours.append(c)
    return neighbours

cdef full_set(bitset_t bits):
    "return a set with 0:len - 1 elements"
    cdef i
    for i in xrange(bitset_len(bits)):
        bitset_add(bits, i)

cpdef is_weak_induced_isomorphism(BinaryMatrix M1_1, BinaryMatrix M2_1,
                                  BinaryMatrix  M_rmat, BinaryMatrix N_rmat):
    """
    For binary matrices this will in fact answer if M1,M2 specify induced
    subgraph isomorphism whereas for higher fields, this only verifies a
    necessary but not sufficient condition
    """
    cdef long i, k
    cdef int j
    cdef int * roworders
    cdef int * colorders
    cdef BinaryMatrix rmat1
    cdef BinaryMatrix M1, M2
    M1 = M1_1.copy()
    M2 = M2_1.copy()
    # rmat corresponding to given isomorphism
    rmat1 = BinaryMatrix(N_rmat.nrows(), N_rmat.ncols())
    roworders = <int * > sage_malloc(N_rmat.nrows() * sizeof(int))
    for i in xrange(M1.nrows()):
        j = bitset_pop(M1._M[i])
        # row i maps to row j
        roworders[i] = j
    colorders = <int * > sage_malloc(N_rmat.ncols() * sizeof(int))
    for i in xrange(M2.nrows()):
        j = bitset_pop(M2._M[i])
        # col i maps to row j
        colorders[i] = j - M1.ncols()
    for i in xrange(M1.nrows()):
        for k in xrange(M2.nrows()):
            rmat1.set_unsafe(i, k, M_rmat.get_unsafe(roworders[i],
                                                     colorders[k]))
    sage_free(roworders)
    sage_free(colorders)
    for i in xrange(N_rmat.nrows()):
        if not bitset_eq(N_rmat._M[i], rmat1._M[i]):
            return 0
    return 1

cdef mat_transpose(BinaryMatrix mat):
    cdef long i, j
    matT = BinaryMatrix(mat.ncols(), mat.nrows())
    for i in xrange(mat.nrows()):
        for j in xrange(mat.ncols()):
            matT.set_unsafe(j, i, mat.get_unsafe(i, j))
    return matT

cpdef bint is_new_rmat(BinaryMatrix cand, used_rmats):
    cdef long i
    for i in xrange(len(used_rmats)):
        if mats_equal(cand, used_rmats[i]):
            return False
    return True

cpdef bint mats_equal(BinaryMatrix M1, BinaryMatrix M2):
    cdef long i
    for i in xrange(M1.nrows()):
        if not bitset_eq(M1._M[i], M2._M[i]):
            return False
    return True
