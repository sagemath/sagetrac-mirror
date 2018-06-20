from libcpp.vector cimport vector
cdef extern from "weisfeiler_lehman/weisfeiler_lehman.cpp":
    pass
cdef extern from "weisfeiler_lehman/weisfeiler_lehman.h" namespace "wl":
    cdef struct GraphNode:
        long long idx, color
        vector[int] adj_list
    
cdef extern from "weisfeiler_lehman/weisfeiler_lehman.h" namespace "wl":
    cdef int prova2(vector[GraphNode] v)
