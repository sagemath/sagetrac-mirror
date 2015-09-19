from libcpp.vector cimport vector

cdef extern from "tdlib/sage_tdlib.hpp":
    int sage_preprocessing_MD(vector[unsigned int] &V_G, vector[unsigned int] &E_G, vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb)
    int sage_preprocessing_FI_TM(vector[unsigned int] &V_G, vector[unsigned int] &E_G, vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb)

    int sage_deltaC_min_d(vector[unsigned int] &V, vector[unsigned int] &E)
    int sage_deltaC_max_d(vector[unsigned int] &V, vector[unsigned int] &E)
    int sage_deltaC_least_c(vector[unsigned int] &V, vector[unsigned int] &E)

    int sage_LBN_deltaC(vector[unsigned int] &V, vector[unsigned int] &E)
    int sage_LBNC_deltaC(vector[unsigned int] &V, vector[unsigned int] &E)
    int sage_LBP_deltaC(vector[unsigned int] &V, vector[unsigned int] &E)
    int sage_LBPC_deltaC(vector[unsigned int] &V, vector[unsigned int] &E)

    int sage_CR_greedy_decomp(vector[unsigned int] &V_G, vector[unsigned int] &E_G, vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb)
    int sage_CR_dynamic_decomp(vector[unsigned int] &V_G, vector[unsigned int] &E_G, vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb)

    int sage_seperator_algorithm(vector[unsigned int] &V_G, vector[unsigned int] &E_G, vector[vector[int]] &V_T, vector[unsigned int] &E_T)

    int sage_MSVS(vector[unsigned int] &V_G, vector[unsigned int] &E_G, vector[vector[int]] &V_T, vector[unsigned int] &E_T)
    void sage_minimalChordal(vector[unsigned int] &V_G, vector[unsigned int] &E_G, vector[unsigned int] &old_elimination_ordering, vector[unsigned int] &new_elimination_ordering)

    int sage_is_valid_decomposition(vector[unsigned int] &V_G, vector[unsigned int] &E_G, vector[vector[int]] &V_T, vector[unsigned int] &E_T)
    int sage_ordering_to_treedec(vector[unsigned int] &V_G, vector[unsigned int] &E_G, vector[vector[int]] &V_T, vector[unsigned int] &E_T, vector[unsigned int] &elim_ordering)
    void sage_treedec_to_ordering(vector[vector[int]] &V, vector[unsigned int] &E, vector[unsigned int] &elim_ordering);
    int sage_get_width(vector[vector[int]] &V_T)
