from libcpp.vector cimport vector

cdef extern from "tdlib/sage_tdlib.hpp":

    int sage_exact_decomposition(vector[unsigned int] &V_G, vector[unsigned int] &E_G, vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb)

