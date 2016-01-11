from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport int64_t

cdef extern from "libopenf4.h":
    cdef vector[string] groebnerBasisF4(int64_t modulus,
                                        int nbVariable,
                                        vector[string] variableName,
                                        vector[string] polynomialList,
                                        int nbThread, int verbose)
    cdef vector[string] groebnerBasisGF2F4(int nbVariable,
                                           vector[string] variableName,
                                           vector[string] polynomialList,
                                           int nbThread, int verbose)
    cdef vector[string] groebnerBasisGF2ExtensionF4(string modulus,
                                                    int nbVariable,
                                                    vector[string] variableName,
                                                    string polyVarName,
                                                    vector[string] polynomialList,
                                                    int nbThread,
                                                    int verbose)
