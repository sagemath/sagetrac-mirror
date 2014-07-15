from lean_matrix cimport LeanMatrix, GenericMatrix, BinaryMatrix, TernaryMatrix, QuaternaryMatrix, IntegerMatrix, generic_identity
from sage.misc.bitset cimport bitset_t


cpdef search_assgn(Ga,Gb,assgn,cand=*)
cpdef init_iso_matrices(BinaryMatrix M_rmat, BinaryMatrix N_rmat)
cpdef col_len(BinaryMatrix mat,long c)
cpdef copy_mat(BinaryMatrix mat)
cdef recurse(bitset_t unused_cols,long cur_row, BinaryMatrix M_rmat, BinaryMatrix N_rmat,BinaryMatrix M_rmatT, BinaryMatrix N_rmatT,BinaryMatrix M1,BinaryMatrix M2)
cpdef prune(long cur_row,BinaryMatrix M1, BinaryMatrix M2, BinaryMatrix M_rmat, BinaryMatrix N_rmat,BinaryMatrix M_rmatT, BinaryMatrix N_rmatT)
cpdef _neighbours(BinaryMatrix rmat, BinaryMatrix rmatT, long i, long rc)
cdef full_set(bitset_t bits)
cpdef degrees_are_sane(BinaryMatrix M1,BinaryMatrix M2)
cpdef is_weak_induced_isomorphism(BinaryMatrix M1,BinaryMatrix M2, BinaryMatrix  M_rmat, BinaryMatrix N_rmat)
cdef mat_transpose(BinaryMatrix mat)
cdef restore_mat(BinaryMatrix mat, BinaryMatrix saved_mat)

#cpdef mat_transpose(mat)
