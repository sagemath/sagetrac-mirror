from .types cimport GF2E_c, GF2X_c, GF2_c, GF2XModulus_c, ZZ_c

cdef extern from "ccobject.h":
    void GF2E_from_str "_from_str<GF2E>"(GF2E_c* dest, char* s)
    object GF2E_to_PyString "_to_PyString<GF2E>"(GF2E_c *x)

cdef extern from "sage/libs/ntl/ntlwrap.cpp":
    void GF2E_init "GF2E::init"(GF2X_c x)
    long GF2E_degree "GF2E::degree"()
    GF2XModulus_c GF2E_modulus "GF2E::modulus"()

    int GF2E_IsOne "IsOne"(GF2E_c x)
    int GF2E_IsZero "IsZero"(GF2E_c x)

    void GF2E_add "add"( GF2E_c x, GF2E_c a, GF2E_c b)
    void GF2E_sub "sub"( GF2E_c x, GF2E_c a, GF2E_c b)
    void GF2E_mul "mul"( GF2E_c x, GF2E_c a, GF2E_c b)
    void GF2E_div "div"( GF2E_c x, GF2E_c a, GF2E_c b)
    void GF2E_power "NTL::power"(GF2E_c t, GF2E_c x, long e)
    long GF2E_deg "deg"(GF2E_c x)

    void GF2E_conv_GF2X "conv" (GF2E_c out, GF2X_c inp)
    void GF2E_conv_long "conv" (GF2E_c out, long inp)
    void GF2E_conv_ZZ "conv" (GF2E_c out, ZZ_c inp)
    void GF2E_conv_GF2 "conv" (GF2E_c out, GF2_c inp)
    GF2X_c GF2E_rep "rep"(GF2E_c x)

    GF2E_c GF2E_random "random_GF2E"()

    GF2_c GF2E_trace "trace"(GF2E_c x)
