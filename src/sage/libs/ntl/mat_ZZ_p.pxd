from .types cimport mat_ZZ_p_c, vec_ZZ_p_c, ZZ_p_c

cdef extern from "sage/libs/ntl/ntlwrap.cpp":
    void mat_ZZ_p_add "add"( mat_ZZ_p_c x, mat_ZZ_p_c a, mat_ZZ_p_c b)
    void mat_ZZ_p_sub "sub"( mat_ZZ_p_c x, mat_ZZ_p_c a, mat_ZZ_p_c b)
    void mat_ZZ_p_mul "mul"( mat_ZZ_p_c x, mat_ZZ_p_c a, mat_ZZ_p_c b)
    void mat_ZZ_p_negate "NTL::negate"(mat_ZZ_p_c x, mat_ZZ_p_c a)
    void mat_ZZ_p_power "NTL::power"(mat_ZZ_p_c t, mat_ZZ_p_c x, long e)
    ZZ_p_c mat_ZZ_p_determinant "determinant"(mat_ZZ_p_c m)
    void mat_ZZ_p_transpose "transpose"(mat_ZZ_p_c r, mat_ZZ_p_c m)
    long mat_ZZ_p_IsZero "IsZero"(mat_ZZ_p_c x)
    void mat_ZZ_p_setitem(mat_ZZ_p_c* x, int i, int j, ZZ_p_c* z)

    long mat_ZZ_p_gauss "gauss"(mat_ZZ_p_c A, long w)
    void mat_ZZ_p_solve "solve"(ZZ_p_c d, vec_ZZ_p_c X, mat_ZZ_p_c A, vec_ZZ_p_c b)
    void mat_ZZ_p_inv "inv" (mat_ZZ_p_c X, mat_ZZ_p_c A)

    long mat_ZZ_p_IsIdent "IsIdent"(mat_ZZ_p_c A, long n)
    long mat_ZZ_p_IsDiag "IsDiag"(mat_ZZ_p_c A, long n, ZZ_p_c d)

    void mat_ZZ_p_image "image"(mat_ZZ_p_c X, mat_ZZ_p_c A)
    void mat_ZZ_p_kernel "kernel" (mat_ZZ_p_c X, mat_ZZ_p_c A)

    void vec_ZZ_p_conv_mat_ZZ_p "conv" (vec_ZZ_p_c out, mat_ZZ_p_c inp)
    void mat_ZZ_p_conv_vec_ZZ_p(mat_ZZ_p_c out, vec_ZZ_p_c inp)
