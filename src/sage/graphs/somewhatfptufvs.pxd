cdef extern from "somewhatfptufvs/externals.hpp":
    struct my_gal:
       pass
       
cdef extern from "somewhatfptufvs/externals.hpp":
    cdef my_gal * gal(int n)
    cdef void gal_an_edge(my_gal *g, int, int, int)
    cdef void free_gal(my_gal *g)
    cdef int UFVS(my_gal *g, int **list)
    cdef int wUFVS(my_gal *g, int **list)
