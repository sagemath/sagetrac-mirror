from .types cimport mat_GF2E_c, vec_GF2E_c, GF2E_c

cdef extern from "ccobject.h":
    void mat_GF2E_from_str "_from_str<mat_GF2E>"(mat_GF2E_c* dest, char* s)
    object mat_GF2E_to_PyString "_to_PyString<mat_GF2E>"(mat_GF2E_c *x)

cdef extern from "sage/libs/ntl/ntlwrap.cpp":
    void mat_GF2E_add "add"( mat_GF2E_c x, mat_GF2E_c a, mat_GF2E_c b)
    void mat_GF2E_sub "sub"( mat_GF2E_c x, mat_GF2E_c a, mat_GF2E_c b)
    void mat_GF2E_mul "mul"( mat_GF2E_c x, mat_GF2E_c a, mat_GF2E_c b)
    void mat_GF2E_negate "NTL::negate"(mat_GF2E_c x, mat_GF2E_c a)
    void mat_GF2E_power "NTL::power"(mat_GF2E_c t, mat_GF2E_c x, long e)
    GF2E_c mat_GF2E_determinant "determinant"(mat_GF2E_c m)
    void mat_GF2E_transpose "transpose"(mat_GF2E_c r, mat_GF2E_c m)
    long mat_GF2E_IsZero "IsZero"(mat_GF2E_c x)
    void mat_GF2E_setitem(mat_GF2E_c* x, int i, int j, GF2E_c* z)

    long mat_GF2E_gauss "gauss"(mat_GF2E_c A, long w)
    void mat_GF2E_solve "solve"(GF2E_c d, vec_GF2E_c X, mat_GF2E_c A, vec_GF2E_c b)
    void mat_GF2E_inv "inv" (mat_GF2E_c X, mat_GF2E_c A)

    long mat_GF2E_IsIdent "IsIdent"(mat_GF2E_c A, long n)
    long mat_GF2E_IsDiag "IsDiag"(mat_GF2E_c A, long n, GF2E_c d)

    void mat_GF2E_image "image"(mat_GF2E_c X, mat_GF2E_c A)
    void mat_GF2E_kernel "kernel" (mat_GF2E_c X, mat_GF2E_c A)

    void vec_GF2E_conv_mat_GF2E "conv" (vec_GF2E_c out, mat_GF2E_c inp)
    void mat_GF2E_conv_vec_GF2E(mat_GF2E_c out, vec_GF2E_c inp)
