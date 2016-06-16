from matroid cimport Matroid            # We'll need this for later.
from sage.matrix.matrix import Matrix
from sage.matrix.constructor import matrix
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.structure.sage_object import SageObject

cdef class Node(sage.structure.sage_object.SageObject):
    cdef G
    cdef pm
    cdef int f

    # def __init__(self, G, pm, int f):
    cdef get_graph(self)
    cdef get_parent_marker(self)
    cdef get_f(self)
    cdef set_f(self, int n)
    cdef is_polygon(self)
    cdef __relink1(self, Z=*, WQ=*, m=*)
    cdef __relink2(self, Z=*, WQ=*)



cdef class CunninghamEdmondsDecomposition(sage.structure.sage_object.SageObject):
    cdef _arborescence  # _T
    cdef dict _nodes  # _ND
    cdef _root
    cdef K_1        #These are useful, because several functions should be able to change them, but passing them is inconvinient.
    cdef K_2        #These are useful, because several functions should be able to change them, but passing them is inconvinient.
    cdef u_1        #These are useful, because several functions should be able to change them, but passing them is inconvinient.
    cdef u_2        #These are useful, because several functions should be able to change them, but passing them is inconvinient.

    cdef __relink1(Q, Z=*, WQ=*, m=*)
    cdef __typing(D_hat, P, pi)
    cdef __relink2(Q, Z=*, WQ=*)
    cdef __hypopath(D, P)
    cdef __update(self, P, C)
    cdef __is_graphic(self)
    cdef __merge(self, G)
    cdef merge(self)
    cdef __add_cycle(self, cycle)