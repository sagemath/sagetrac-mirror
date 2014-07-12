from lean_matrix cimport LeanMatrix, GenericMatrix, BinaryMatrix, TernaryMatrix, QuaternaryMatrix, IntegerMatrix, generic_identity
from sage.misc.bitset cimport bitset_t


cpdef search_assgn(Ga,Gb,assgn,cand=*)
cpdef init_iso_matrices(BinaryMatrix M_rmat, BinaryMatrix N_rmat)
cpdef col_len(BinaryMatrix mat,long c)
cpdef copy_mat(BinaryMatrix mat)
cdef recurse(bitset_t unused_cols,long cur_row, BinaryMatrix M_rmat, BinaryMatrix N_rmat,BinaryMatrix M_rmatT, BinaryMatrix N_rmatT,BinaryMatrix M1,BinaryMatrix M2, BinaryMatrix M1_0, BinaryMatrix M2_0)
cpdef prune(BinaryMatrix M1, BinaryMatrix M2, long cur_row,BinaryMatrix M_rmat, BinaryMatrix N_rmat,BinaryMatrix M_rmatT, BinaryMatrix N_rmatT)
cpdef _neighbours(BinaryMatrix rmat, BinaryMatrix rmatT, long i, long rc)
cdef full_set(bitset_t bits)

#cpdef mat_transpose(mat)
