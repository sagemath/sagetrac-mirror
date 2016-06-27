from matroid cimport Matroid            # We'll need this for later.
from sage.matrix.matrix import Matrix
from sage.matrix.constructor import matrix
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.structure.sage_object import SageObject

cdef class Node(sage.structure.sage_object.SageObject):
    cpdef public __custom_name
    cpdef public _custom_name
    cpdef public _cached_info
    cpdef graph
    cdef int parent_marker
    cdef int f

    # def __init__(self, G, pm, int f):
    cpdef get_graph(self)
    cdef int get_parent_marker(self)
    cpdef get_parent_marker_edge(self)
    cdef int get_f(self)
    cpdef set_f(self, int n)
    cpdef is_polygon(self)
    cpdef is_path(self, P)
    cpdef is_cycle(self, P)
    cpdef typing(self, P)
    cpdef __relink1(self, Z=*, WQ=*)
    cpdef __relink2(self, Z=*, WQ=*)



cdef class CunninghamEdmondsDecomposition(sage.structure.sage_object.SageObject):
    cpdef arborescence  # _T
    cdef list nodes  # _ND
    cdef int root
    cdef int K_1     #These are useful, because several functions should be able to change them, but passing them is inconvinient.
    cdef int K_2     #These are useful, because several functions should be able to change them, but passing them is inconvinient.
    cdef int u_1     #These are useful, because several functions should be able to change them, but passing them is inconvinient.
    cdef int u_2     #These are useful, because several functions should be able to change them, but passing them is inconvinient.
    cdef int next_arborescence_vertex
    cdef int next_vertex
    cdef int next_edge

    cpdef __relink1(Q, Z=*, WQ=*)
    cpdef __typing(self, P, pi)
    cpdef __relink2(Q, Z=*, WQ=*)
    cpdef __hypopath(self, P)
    # cpdef __squeeze(self, L)
    cpdef __update(self, P, C)
    cpdef __is_graphic(self)
    cpdef __merge(self, G)
    cpdef merge(self)
    cpdef __add_cycle(self, cycle)
    cpdef get_arborescence(self)
    cpdef get_nodes(self)
    cpdef get_root(self)
    cpdef __get_pi(self)
    cpdef branch(self, N)
    cdef get_parent(self, N)