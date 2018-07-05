from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
cdef extern from "weisfeiler_lehman/weisfeiler_lehman.cpp":
    pass
cdef extern from "weisfeiler_lehman/Tuple.h":
    pass
cdef extern from "weisfeiler_lehman/weisfeiler_lehman.h" namespace "wl":
    cdef struct GraphNode:
        long long idx, color
        vector[int] adj_list
cdef extern from "weisfeiler_lehman/Tuple.h":
    cdef cppclass Tuple[int]:
        int* content
        Tuple()
        string to_string()
        int* begin()
        int* end()
cdef extern from "weisfeiler_lehman/weisfeiler_lehman.h" namespace "wl":
    cdef unordered_map[Tuple[int], int] k_WL(const vector[GraphNode]& v, int k)
