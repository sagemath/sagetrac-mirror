r"""
Transversal matroids
"""

from __future__ import print_function, absolute_import

from sage.matroids.matroid cimport Matroid
from sage.matroids.basis_exchange_matroid cimport BasisExchangeMatroid

from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph

cdef class TransversalMatroid(BasisExchangeMatroid):
    r"""
    Transversal matroids

    EXAMPLES::

        sage: from sage.matroids.transversal_matroid import TransversalMatroid
        sage: B = BipartiteGraph(graphs.CompleteBipartiteGraph(4,3))
        sage: M = TransversalMatroid(B, [0, 1, 2, 3])
        sage: M.full_rank()
        3

    """

    def __init__(self, B, groundset):
        """
        See class definition for full documentation.
        """
        # Make this work with a bipartite graph as input
        # In a later ticket, make the constructor work with a collection of sets as input
        # Assumptions: `B` is an instance of BipartiteGraph
        #               `groundset` is one side (not optional)
        # Later, the constructor should be able to accept Graph as input
        # And convert to BipartiteGraph if the ground set makes sense
        # TODO: Let it guess the ground set from the largest of B.bipartition()

        self._buckets = [B.neighbors(v) for v in set(B.vertices()) if v not in groundset]

        self._matching = B.matching()

        vertices_in_matching = set([u for u, v, l in self._matching]).union(
            set([v for u, v, l in self._matching]))
        basis = frozenset([v for v in vertices_in_matching if v in groundset])

        BasisExchangeMatroid.__init__(self, groundset, basis)

        # Build a DiGraph for doing basis exchange
        self._D = DiGraph()
        for u, v, l in B.edge_iterator():
            if (u, v, l) in self._matching:
            # For the edges in our matching, orient them as starting from the collections
                if u in self._groundset:
                    self._D.add_edge(v, u, l)
                else:
                    self._D.add_edge(u, v, l)
            else:
            # Otherwise orient them as starting from the ground set
                if u in self._groundset:
                    self._D.add_edge(u, v, l)
                else:
                    self._D.add_edge(v, u, l)


    cdef bint __is_exchange_pair(self, long x, long y) except -1:
        r"""
        Check for `M`-alternating path from `x` to `y`.
        """
        # Question: Do I need to consider exchanges between `x` and `x`?
        if self._D.shortest_path(x, y):
            return True
        else:
            return False

    cdef int __exchange(self, long x, long y) except -1:
        r"""
        Replace ``self.basis() with ``self.basis() - x + y``. Internal method, does no checks.
        """
        shortest_path = self._D.shortest_path(x, y)
        self._D.reverse_edges(shortest_path)

        for u, v, l in shortest_path:
            if (u, v, l) in self._matching:
                self._matching.remove((u, v, l))
            else:
                self._matching.append((v, u, l))

        BasisExchangeMatroid.__exchange(self, x, y)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: B = BipartiteGraph(graphs.CompleteBipartiteGraph(4,3))
            sage: M = TransversalMatroid(B, [0, 1, 2, 3]); M
            Transversal matroid of rank 3 on 4 elements, with 3 subsets.
        """
        S = ("Transversal matroid of rank " + str(self.rank()) + " on "
            + str(self.size()) + " elements, with " + str(len(self._buckets))
            + " subsets.")
        return S
