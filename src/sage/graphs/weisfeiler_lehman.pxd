from libcpp.vector cimport vector
from libcpp cimport bool
cdef extern from "weisfeiler_lehman/weisfeiler_lehman.cpp":
    pass
cdef extern from "weisfeiler_lehman/Tuple.h":
    pass
cdef extern from "weisfeiler_lehman/weisfeiler_lehman.h" namespace "wl":
    cdef struct GraphNode:
        long long idx, color
        vector[int] adj_list

cdef extern from "weisfeiler_lehman/weisfeiler_lehman.h" namespace "wl":
    cdef bool k_WL(const vector[GraphNode]& v, int k)
