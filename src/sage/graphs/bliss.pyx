include "sage/ext/interrupt.pxi"
include 'sage/ext/stdsage.pxi'
include 'sage/ext/cdefs.pxi'

cdef extern from "/home/azi/bliss-0.72/graph.hh"  namespace "bliss":

    cdef cppclass AbstractGraph:
        pass

    cdef cppclass Graph(AbstractGraph):
        Graph(const unsigned int)

def foo():
    cdef Graph *g = new Graph(10)
