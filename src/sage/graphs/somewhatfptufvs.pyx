r"""
somewhatfptUFVS: a package intended for solving Feedback Vertx Set problems. So far (unweighted?) undirected graphs.
I am sure someone has the time to put together a simple straight of the paper implement for those directed graphs just for completeness.
Aught to be tempted if knew how to get the minimum cut from the implementation in sage for maximum flow.

This module defines the use of this package.

AUTHORS:
		20.01.2014 	Arvid Soldal Sivertsen	-	Pushed a working implementation for (weighted?) undirected graphs. (at least he hope so)

REFERENCE:
		Improved Algorithms for the
		Feedback Vertex Set Problems
	which
		are all fpt algorithms (poor you if you have trubble with the parts of this package that are not good enough to be deemed / pass as fpt)
	by
		Jianer Chen
		Fedor V. Fomin
		Yang Liu
		Songjian Lu
		Yngve Villanger
	
Method
-------
"""

include "sage/ext/interrupt.pxi"
include 'sage/ext/stdsage.pxi'

def minimum_FVS(graph):
    """
    Returns the vertex set of a minimum FVS.

    ......
    
    EXAMPLES::

          sage: C=graphs.PetersenGraph()
          sage: minimum_FVS(C)
          [0, 8, 9]

    Lame TEST::

        sage: g = Graph()
        sage: g.minimum_FVS()
        []
        
.. todo:: assert that all weights are positive

.. todo:: might not be implemented for weighted undirected graphs yet

.. todo:: There is an issue with weighted() only apply for edges, not vertices.

    """
    if graph.order() == 0:
        return []

    cdef my_gal *g
    g=gal(graph.order())
    for e in graph.edge_iterator():
        (u,v,w)=e
        gal_an_edge(g,u,v)

    cdef int* list
    cdef int size
    sig_on()
    if graph.weighted():
      size = wUFVS(g, &list)
    else:
      size = UFVS(g, &list)
    sig_off()
    v = []
    cdef int i
    for i in range(size):
        v.append(list[i])

    sage_free(list)
    free_gal(g)
    return v
    
