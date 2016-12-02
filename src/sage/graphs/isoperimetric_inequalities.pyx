r"""
Cheeger constants and isoperimetric inequalities
"""

from __future__ import print_function

include "sage/data_structures/bitset.pxi"
from libc.stdint cimport uint32_t
from libc.stdlib cimport free
from sage.graphs.base.static_sparse_graph cimport short_digraph, init_short_digraph, free_short_digraph, out_degree

from sage.rings.rational_field import QQ

def cheeger(gg):
    r"""
    Return the cheeger constant of the graph ``gg``

    EXAMPLES::

        sage: from sage.graphs.isoperimetric_inequalities import cheeger

        sage: cheeger(graphs.PetersenGraph())
        1/3
        sage: graphs.PetersenGraph().cheeger_constant()
        1/3

        sage: [cheeger(graphs.CycleGraph(k)) for k in range(2,10)]
        [1, 1, 1/2, 1/2, 1/3, 1/3, 1/4, 1/4]
        sage: [graphs.CycleGraph(k).cheeger_constant() for k in range(2,10)]
        [1, 1, 1/2, 1/2, 1/3, 1/3, 1/4, 1/4]

        sage: [cheeger(graphs.CompleteGraph(k)) for k in range(2,10)]
        [1, 1, 2/3, 3/4, 3/5, 2/3, 4/7, 5/8]
        sage: [graphs.CompleteGraph(k).cheeger_constant() for k in range(2,10)]
        [1, 1, 2/3, 3/4, 3/5, 2/3, 4/7, 5/8]

        sage: [cheeger(graphs.CompleteBipartiteGraph(i,j)) for i in range(2,7) for j in range(2, i)]
        [3/5, 1/2, 3/5, 5/9, 4/7, 5/9, 1/2, 5/9, 1/2, 5/9]
        sage: [graphs.CompleteBipartiteGraph(i,j).cheeger_constant() for i in range(2,7) for j in range(2, i)]
        [3/5, 1/2, 3/5, 5/9, 4/7, 5/9, 1/2, 5/9, 1/2, 5/9]

    More examples::

        sage: G = Graph([(0, 1), (0, 3), (0, 8), (1, 4), (1, 6), (2, 4), (2, 7), (2, 9),
        ....:            (3, 6), (3, 8), (4, 9), (5, 6), (5, 7), (5, 8), (7, 9)])
        sage: cheeger(G)
        1/6
        sage: G.cheeger_constant()
        1/6

        sage: G = Graph([(0, 1), (0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (3, 4), (3, 5)])
        sage: cheeger(G)
        1/2
        sage: G.cheeger_constant()
        1/2
    """
    assert not gg.is_directed()
    gg._scream_if_not_simple()
    if gg.num_verts() == 0:
        raise ValueError("Cheeger constant is not defined for the empty graph")
    elif gg.num_verts() == 1:
        from sage.rings.infinity import infinity
        return infinity

    cdef short_digraph sd           # a copy of the graph gg
    cdef uint32_t * subgraph        # subset of vertices (dense representation)
    cdef bitset_t bitsubgraph       # subset of vertices (sparse representation)
    cdef int k = 0                  # number of vertices in subgraph
    cdef unsigned long vol = 0      # number of edges in the subgraph
    cdef unsigned long boundary = 0 # number of edges in the boundary
    cdef int u = 0                  # current vertex
    cdef unsigned long bmin = 1     # value of boundary for the min
    cdef unsigned long vmin = 1     # value of the volum for the min
    cdef int i,j

    init_short_digraph(sd, gg)

    subgraph = <uint32_t *> check_malloc(sd.n * sizeof(uint32_t))
    subgraph

    bitset_init(bitsubgraph, sd.n)
    bitset_add(bitsubgraph, 0)

    vol = boundary = 0
    while True:
        assert u < sd.n, u
        while True:
            # add vertex u to the subgraph, update the boundary/volume
            # we repeat the operation until we reach the maximum volume
            # or have no more vertex available
            for i in range(sd.neighbors[u+1] - sd.neighbors[u]):
                j = sd.neighbors[u][i]
                if bitset_in(bitsubgraph, j):
                    boundary -= 1
                else:
                    boundary += 1
                vol += 1
            bitset_add(bitsubgraph, u)
            subgraph[k] = u
            u += 1
            k += 1

            if vol > sd.m:
                break

            if boundary * vmin < bmin * vol:
                bmin = boundary
                vmin = vol

            if u == sd.n:
                break

        # backtrack
        k -= 1
        u = subgraph[k]

        bitset_remove(bitsubgraph, u)
        for i in range(sd.neighbors[u+1] - sd.neighbors[u]):
            j = sd.neighbors[u][i]
            if bitset_in(bitsubgraph, j):
                boundary += 1
            else:
                boundary -= 1
            vol -= 1
        u += 1

        if u == sd.n:
            if k == 0:
                # end of the loop
                break
            else:
                # remove one more vertex in order to continue
                k -= 1
                u = subgraph[k]

                bitset_remove(bitsubgraph, u)
                for i in range(sd.neighbors[u+1] - sd.neighbors[u]):
                    j = sd.neighbors[u][i]
                    if bitset_in(bitsubgraph, j):
                        boundary += 1
                    else:
                        boundary -= 1
                    vol -= 1
                u += 1

    free_short_digraph(sd)
    free(subgraph)
    bitset_free(bitsubgraph)

    return QQ((bmin, vmin))
