from lean_matrix cimport LeanMatrix, GenericMatrix, BinaryMatrix, TernaryMatrix, QuaternaryMatrix, IntegerMatrix, generic_identity


cpdef search_assgn(Ga,Gb,assgn,cand=*)
cpdef init_iso_matrices(BinaryMatrix M_rmat, BinaryMatrix N_rmat)
cpdef col_len(BinaryMatrix mat,long c)

#cpdef mat_transpose(mat)
