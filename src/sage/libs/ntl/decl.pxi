cdef extern from "Python.h":
    object PyString_FromString(char *v)
    char* PyString_AsString(object string)

cdef extern from "stdlib.h":
    void free(void *ptr)

cdef extern from "ntl_wrap.h":
    #### ZZ
    struct ZZ
    ZZ* new_ZZ()
    ZZ* str_to_ZZ(char* s)
    void del_ZZ(ZZ* n)
    char* ZZ_to_str(ZZ* x)
    ZZ* ZZ_add(ZZ* x, ZZ* y)
    ZZ* ZZ_sub(ZZ* x, ZZ* y)
    ZZ* ZZ_mul(ZZ* x, ZZ* y)
    ZZ* ZZ_pow(ZZ* x, long e)
    int ZZ_is_one(ZZ* x)
    int ZZ_is_zero(ZZ* x)
    ZZ* ZZ_neg(ZZ* x)
    ZZ* ZZ_copy(ZZ* x)
    # Random-number generation
    void setSeed(ZZ* x)
    ZZ* ZZ_randomBnd(ZZ* x)
    ZZ* ZZ_randomBits(long n)

    #### ZZ_p
    struct ZZ_p
    ZZ_p* new_ZZ_p()
    ZZ_p* str_to_ZZ_p(char* s)
    void del_ZZ_p(ZZ_p* n)
    char* ZZ_p_to_str(ZZ_p* x)
    ZZ_p* ZZ_p_add(ZZ_p* x, ZZ_p* y)
    ZZ_p* ZZ_p_sub(ZZ_p* x, ZZ_p* y)
    ZZ_p* ZZ_p_mul(ZZ_p* x, ZZ_p* y)
    ZZ_p* ZZ_p_pow(ZZ_p* x, long e)
    int ZZ_p_is_one(ZZ_p* x)
    int ZZ_p_is_zero(ZZ_p* x)
    ZZ_p* ZZ_p_neg(ZZ_p* x)
    void ntl_ZZ_set_modulus(ZZ* x)
    int ZZ_p_eq( ZZ_p* x,  ZZ_p* y)
    ZZ_p* ZZ_p_inv(ZZ_p* x)
    ZZ_p* ZZ_p_random()

    #### ZZX
    struct ZZX
    ZZX* ZZX_init()
    ZZX* str_to_ZZX(char* s)
    char* ZZX_repr(ZZX* x)
    void ZZX_dealloc(ZZX* x)
    ZZX* ZZX_copy(ZZX* x)
    void ZZX_setitem(ZZX* x, long i, char* a)
    char* ZZX_getitem(ZZX* x, long i)
    ZZX* ZZX_add(ZZX* x, ZZX* y)
    ZZX* ZZX_sub(ZZX* x, ZZX* y)
    ZZX* ZZX_mul(ZZX* x, ZZX* y)
    ZZX* ZZX_div(ZZX* x, ZZX* y, int* divisible)
    ZZX* ZZX_mod(ZZX* x, ZZX* y)
    void ZZX_quo_rem(ZZX* x, ZZX* other, ZZX** r, ZZX** q)
    ZZX* ZZX_square(ZZX* x)
    int ZZX_equal(ZZX* x, ZZX* y)
    int ZZX_is_zero(ZZX* x)
    int ZZX_is_one(ZZX* x)
    int ZZX_is_monic(ZZX* x)
    ZZX* ZZX_neg(ZZX* x)
    ZZX* ZZX_left_shift(ZZX* x, long n)
    ZZX* ZZX_right_shift(ZZX* x, long n)
    char* ZZX_content(ZZX* x)
    ZZX* ZZX_primitive_part(ZZX* x)
    void ZZX_pseudo_quo_rem(ZZX* x, ZZX* y, ZZX** r, ZZX** q)
    ZZX* ZZX_gcd(ZZX* x, ZZX* y)
    ZZX* ZZX_xgcd(ZZX* x, ZZX* y, ZZ** r, ZZX** s, ZZX** t, int proof)
    long ZZX_degree(ZZX* x)
    ZZ* ZZX_leading_coefficient(ZZX* x)
    char* ZZX_constant_term(ZZX* x)
    void ZZX_set_x(ZZX* x)
    int ZZX_is_x(ZZX* x)
    ZZX* ZZX_derivative(ZZX* x)
    ZZX* ZZX_reverse(ZZX* x)
    ZZX* ZZX_reverse_hi(ZZX* x, long hi)
    ZZX* ZZX_truncate(ZZX* x, long m)
    ZZX* ZZX_multiply_and_truncate(ZZX* x, ZZX* y, long m)
    ZZX* ZZX_square_and_truncate(ZZX* x, long m)
    ZZX* ZZX_invert_and_truncate(ZZX* x, long m)
    ZZX* ZZX_multiply_mod(ZZX* x, ZZX* y,  ZZX* modulus)
    ZZ* ZZX_trace_mod(ZZX* x, ZZX* y)
    char* ZZX_trace_list(ZZX* x)
    ZZ* ZZX_resultant(ZZX* x, ZZX* y, int proof)
    ZZ* ZZX_norm_mod(ZZX* x, ZZX* y, int proof)
    ZZ* ZZX_discriminant(ZZX* x, int proof)
    ZZ* ZZX_polyeval(ZZX* x, ZZ* a)
    ZZX* ZZX_charpoly_mod(ZZX* x, ZZX* y, int proof)
    ZZX* ZZX_minpoly_mod(ZZX* x, ZZX* y)
    void ZZX_clear(ZZX* x)
    void ZZX_preallocate_space(ZZX* x, long n)


    #### ZZ_pX
    struct ZZ_pX
    ZZ_pX* ZZ_pX_init()
    ZZ_pX* str_to_ZZ_pX(char* s)
    char* ZZ_pX_repr(ZZ_pX* x)
    void ZZ_pX_dealloc(ZZ_pX* x)
    ZZ_pX* ZZ_pX_copy(ZZ_pX* x)
    void ZZ_pX_setitem(ZZ_pX* x, long i, char* a)
    char* ZZ_pX_getitem(ZZ_pX* x, long i)
    ZZ_pX* ZZ_pX_add(ZZ_pX* x, ZZ_pX* y)
    ZZ_pX* ZZ_pX_sub(ZZ_pX* x, ZZ_pX* y)
    ZZ_pX* ZZ_pX_mul(ZZ_pX* x, ZZ_pX* y)
    ZZ_pX* ZZ_pX_div(ZZ_pX* x, ZZ_pX* y, int* divisible)
    ZZ_pX* ZZ_pX_mod(ZZ_pX* x, ZZ_pX* y)
    void ZZ_pX_quo_rem(ZZ_pX* x, ZZ_pX* other, ZZ_pX** r, ZZ_pX** q)
    ZZ_pX* ZZ_pX_square(ZZ_pX* x)
    int ZZ_pX_equal(ZZ_pX* x, ZZ_pX* y)
    int ZZ_pX_is_zero(ZZ_pX* x)
    int ZZ_pX_is_one(ZZ_pX* x)
    int ZZ_pX_is_monic(ZZ_pX* x)
    ZZ_pX* ZZ_pX_neg(ZZ_pX* x)
    ZZ_pX* ZZ_pX_left_shift(ZZ_pX* x, long n)
    ZZ_pX* ZZ_pX_right_shift(ZZ_pX* x, long n)
    ZZ_pX* ZZ_pX_gcd(ZZ_pX* x, ZZ_pX* y)
    ZZ_pX* ZZ_pX_xgcd(ZZ_pX** d, ZZ_pX** s, ZZ_pX** t, ZZ_pX* a, ZZ_pX* b)
    ZZ_pX* ZZ_pX_plain_xgcd(ZZ_pX** d, ZZ_pX** s, ZZ_pX** t, ZZ_pX* a, ZZ_pX* b)
    long ZZ_pX_degree(ZZ_pX* x)
    ZZ_p* ZZ_pX_leading_coefficient(ZZ_pX* x)
    char* ZZ_pX_constant_term(ZZ_pX* x)
    void ZZ_pX_set_x(ZZ_pX* x)
    int ZZ_pX_is_x(ZZ_pX* x)
    ZZ_pX* ZZ_pX_derivative(ZZ_pX* x)
    ZZ_pX* ZZ_pX_reverse(ZZ_pX* x)
    ZZ_pX* ZZ_pX_reverse_hi(ZZ_pX* x, long hi)
    ZZ_pX* ZZ_pX_truncate(ZZ_pX* x, long m)
    ZZ_pX* ZZ_pX_multiply_and_truncate(ZZ_pX* x, ZZ_pX* y, long m)
    ZZ_pX* ZZ_pX_square_and_truncate(ZZ_pX* x, long m)
    ZZ_pX* ZZ_pX_invert_and_truncate(ZZ_pX* x, long m)
    ZZ_pX* ZZ_pX_multiply_mod(ZZ_pX* x, ZZ_pX* y,  ZZ_pX* modulus)
    ZZ_p* ZZ_pX_trace_mod(ZZ_pX* x, ZZ_pX* y)
    char* ZZ_pX_trace_list(ZZ_pX* x)
    ZZ_p* ZZ_pX_resultant(ZZ_pX* x, ZZ_pX* y)
    ZZ_p* ZZ_pX_norm_mod(ZZ_pX* x, ZZ_pX* y)
    ZZ_p* ZZ_pX_discriminant(ZZ_pX* x)
    ZZ_pX* ZZ_pX_charpoly_mod(ZZ_pX* x, ZZ_pX* y)
    ZZ_pX* ZZ_pX_minpoly_mod(ZZ_pX* x, ZZ_pX* y)
    void ZZ_pX_clear(ZZ_pX* x)
    void ZZ_pX_preallocate_space(ZZ_pX* x, long n)

    void ZZ_pX_factor(ZZ_pX*** v, long** e, long* n, ZZ_pX* x, long verbose)
    void ZZ_pX_linear_roots(ZZ_p*** v, long* n, ZZ_pX* x)

    #### mat_ZZ
    struct mat_ZZ
    mat_ZZ* new_mat_ZZ(long nrows, long ncols)
    void del_mat_ZZ(mat_ZZ* n)
    char* mat_ZZ_to_str(mat_ZZ* x)
    mat_ZZ* mat_ZZ_add(mat_ZZ* x, mat_ZZ* y)
    mat_ZZ* mat_ZZ_sub(mat_ZZ* x, mat_ZZ* y)
    mat_ZZ* mat_ZZ_mul(mat_ZZ* x, mat_ZZ* y)
    mat_ZZ* mat_ZZ_pow(mat_ZZ* x, long e)
    long mat_ZZ_nrows(mat_ZZ* x)
    long mat_ZZ_ncols(mat_ZZ* x)
    void mat_ZZ_setitem(mat_ZZ* x, int i, int j, ZZ* z)
    ZZ* mat_ZZ_getitem(mat_ZZ* x, int i, int j)
    ZZ* mat_ZZ_determinant(mat_ZZ* x, long deterministic)
    mat_ZZ* mat_ZZ_HNF(mat_ZZ* A, ZZ* D)
    ZZX* mat_ZZ_charpoly(mat_ZZ* A)
    long mat_ZZ_LLL(ZZ **det, mat_ZZ *x, long a, long b, long verbose)
    long mat_ZZ_LLL_U(ZZ **det, mat_ZZ *x, mat_ZZ *U, long a, long b, long verbose)

    #### GF2X
    struct GF2X
    GF2X* new_GF2X()
    GF2X* str_to_GF2X(char* s)
    void del_GF2X(GF2X* n)
    char* GF2X_to_str(GF2X* x)
    GF2X* GF2X_add(GF2X* x, GF2X* y)
    GF2X* GF2X_sub(GF2X* x, GF2X* y)
    GF2X* GF2X_mul(GF2X* x, GF2X* y)
    GF2X* GF2X_pow(GF2X* x, long e)
    int GF2X_eq( GF2X* x,  GF2X* y)
    int GF2X_is_one(GF2X* x)
    int GF2X_is_zero(GF2X* x)
    GF2X* GF2X_neg(GF2X* x)
    GF2X* GF2X_copy(GF2X* x)
    long GF2X_deg(GF2X* x)
    void GF2X_hex(long h)
    char *GF2X_to_bin(GF2X* x)
    char *GF2X_to_hex(GF2X* x)


    #### GF2E
    struct GF2E
    void ntl_GF2E_set_modulus(GF2X *x)
    GF2E* new_GF2E()
    GF2E* str_to_GF2E(char *s)
    void del_GF2E(GF2E *n)
    char *GF2E_to_str(GF2E* x)
    GF2E *GF2E_add(GF2E *x, GF2E *y)
    GF2E *GF2E_sub(GF2E *x, GF2E *y)
    GF2E *GF2E_mul(GF2E *x, GF2E *y)
    GF2E *GF2E_pow(GF2E *x, long e)
    int GF2E_eq( GF2E* x,  GF2E* y)
    int GF2E_is_one(GF2E *x)
    int GF2E_is_zero(GF2E *x)
    GF2E *GF2E_neg(GF2E *x)
    GF2E *GF2E_copy(GF2E *x)
    long GF2E_degree()
    GF2X *GF2E_modulus()
    GF2E *GF2E_random()
    long  GF2E_trace(GF2E *x)
    GF2X *GF2E_ntl_GF2X(GF2E *x)

    #### GF2EX
    struct GF2EX
    GF2EX* new_GF2EX()
    GF2EX* str_to_GF2EX(char* s)
    void del_GF2EX(GF2EX* n)
    char* GF2EX_to_str(GF2EX* x)
    GF2EX* GF2EX_add(GF2EX* x, GF2EX* y)
    GF2EX* GF2EX_sub(GF2EX* x, GF2EX* y)
    GF2EX* GF2EX_mul(GF2EX* x, GF2EX* y)
    GF2EX* GF2EX_pow(GF2EX* x, long e)
    int GF2EX_is_one(GF2EX* x)
    int GF2EX_is_zero(GF2EX* x)
    GF2EX* GF2EX_neg(GF2EX* x)
    GF2EX* GF2EX_copy(GF2EX* x)


    #### mat_GF2E
    struct mat_GF2E
    mat_GF2E* new_mat_GF2E(long nrows, long ncols)
    void del_mat_GF2E(mat_GF2E* n)
    char* mat_GF2E_to_str(mat_GF2E* x)
    mat_GF2E* mat_GF2E_add(mat_GF2E* x, mat_GF2E* y)
    mat_GF2E* mat_GF2E_sub(mat_GF2E* x, mat_GF2E* y)
    mat_GF2E* mat_GF2E_mul(mat_GF2E* x, mat_GF2E* y)
    mat_GF2E* mat_GF2E_pow(mat_GF2E* x, long e)
    long mat_GF2E_nrows(mat_GF2E* x)
    long mat_GF2E_ncols(mat_GF2E* x)
    void mat_GF2E_setitem(mat_GF2E* x, int i, int j, GF2E* z)
    GF2E* mat_GF2E_getitem(mat_GF2E* x, int i, int j)
    GF2E* mat_GF2E_determinant(mat_GF2E* x)
    long mat_GF2E_gauss(mat_GF2E *x, long w)
    long mat_GF2E_is_zero(mat_GF2E *x)
    mat_GF2E* mat_GF2E_transpose(mat_GF2E *x)

