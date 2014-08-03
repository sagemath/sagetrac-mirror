from lean_matrix cimport  BinaryMatrix
from sage.misc.bitset cimport bitset_t


cdef init_iso_matrices(BinaryMatrix M_rmat, BinaryMatrix N_rmat)
cdef col_len(BinaryMatrix mat, long c)
cdef BinaryMatrix copy_mat(BinaryMatrix mat)
cdef recurse(bitset_t unused_cols, long cur_row, BinaryMatrix M_rmat,
             BinaryMatrix N_rmat, BinaryMatrix M_rmatT, BinaryMatrix N_rmatT,
             BinaryMatrix M1, BinaryMatrix M2, indices, M, Npcl, long nloops)
cdef prune(long cur_row, BinaryMatrix M1, BinaryMatrix M2,
            BinaryMatrix M_rmat, BinaryMatrix N_rmat, BinaryMatrix M_rmatT,
            BinaryMatrix N_rmatT)
cdef _neighbours(BinaryMatrix rmat, BinaryMatrix rmatT, long i, long rc)
cdef degrees_are_sane(BinaryMatrix M1, BinaryMatrix M2)
cdef is_weak_induced_isomorphism(BinaryMatrix M1, BinaryMatrix M2,
                                  BinaryMatrix  M_rmat, BinaryMatrix N_rmat)
cdef mat_transpose(BinaryMatrix mat)
cdef restore_mat(BinaryMatrix mat, BinaryMatrix saved_mat, long cur_row)
cdef iso_mats_to_dict(BinaryMatrix M1, BinaryMatrix M2, indices)
cdef _pcl_have_maps(M, Npcl, iso_map, indices, long nloops)
cpdef _check_bin_minor(BinaryMatrix M_rmat, BinaryMatrix N_rmat, indices, M,
                       Npcl, long nloops, bint debug=*)
cpdef bint is_new_rmat(BinaryMatrix cand, used_rmats)
cdef bint mats_equal(BinaryMatrix M1, BinaryMatrix M2)
