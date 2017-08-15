r"""
Transversal matroids

A transversal matroid arises from a ground set `E` and a collection `A` of sets
over the ground set. This can be modeled as a bipartite graph `B`, where the
vertices the left are ground set elements, the vertices on the right are the
sets, and edges represent containment. Then a set `X` from the ground set is
independent if and only if `X` has a matching in `B`.

Construction
============

Exposure from the constructor and :mod:`sage.matroids.advanced` will be added later.

To construct a transversal matroid, first import TransversalMatroid from
:mod:`sage.matroids.transversal_matroid`.
The input should be a set system::

    sage: from sage.matroids.transversal_matroid import*
    sage: from sage.matroids.transversal_matroid import *
    sage: sets = [[3,4,5,6,7,8]] * 3
    sage: M = TransversalMatroid(sets); M
    Transversal matroid of rank 3 on 6 elements, with 3 sets.
    sage: M.groundset()
    frozenset({3, 4, 5, 6, 7, 8})
    sage: M.is_isomorphic(matroids.Uniform(3,6))  # currently broken
    True


AUTHORS:

- Zachary Gershkoff (2017-08-07): initial version

REFERENCES
==========

..  [JEB17] \Joseph E. Bonin, Lattices Related to Extensions of Presentations of Transversal Matroids. In The Electronic Journal of Combinatorics (2017), #P1.49

Methods
=======
"""
#*****************************************************************************
#       Copyright (C) 2017 Zachary Gershkoff <zgersh2@lsu.edu>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function, absolute_import

include 'sage/data_structures/bitset.pxi'

from sage.matroids.matroid cimport Matroid
from sage.matroids.basis_exchange_matroid cimport BasisExchangeMatroid
from sage.matroids.minor_matroid import MinorMatroid
from sage.matroids.utilities import newlabel

from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from sage.graphs.bipartite_graph import BipartiteGraph
from sage.rings.integer import Integer

from cpython.object cimport Py_EQ, Py_NE
from copy import copy, deepcopy
from collections import Counter

import networkx as nx

cpdef graph_from_buckets(buckets, groundset):
    r"""
    Construct a bipartite graph from sets over a given ground set.

    INPUT:

    - ``buckets`` -- A frozenset consisting of tuples. Each tuple corresponds
      to a vertex on the right side and its neighbors.
    - ``groundset`` -- A list of vertices on the left side.

    OUTPUT:

    A ``BipartiteGraph``.
    """
    # this will have to be rewritten entirely
    # groundset could be made optional, but we require it for efficiency's sake
    B = BipartiteGraph()
    for e in groundset:
        B.add_vertex(e, left=True)
    for s in buckets:
        v = s[0]
        B.add_vertex(v, right=True)
        for e in s[1]:
            if e not in groundset:
                raise ValueError("buckets do not match ground set")
            B.add_edge(e, v)

    return B

cdef class TransversalMatroid(BasisExchangeMatroid):
    r"""
    The Transversal Matroid class.

    INPUT:

    - ``B`` -- A SageMath graph.
    - ``groundset`` -- (optional) An iterable containing names of ground set
      elements. If the ground set is not specified and ``B`` is an instance of
      ``BipartiteGraph``, the ground set, can be assumed to be the vertices on the
      largest side, or the left side of both sides have the same cardinality.
      If ``B`` is not an instance of ``BipartiteGraph``, the ground set must be
      specified.

    OUTPUT:

    An instance of ``TransversalMatroid``.

    EXAMPLES::

        sage: from sage.matroids.transversal_matroid import TransversalMatroid
        sage: sets = [[0, 1, 2, 3]] * 3
        sage: M = TransversalMatroid(sets)
        sage: M.full_rank()
        3
        sage: M.bases_count()
        4
        sage: sum(1 for b in M.bases())  # some of the code thinks 3 is a coloop
        4

    ::

        sage: from sage.matroids.transversal_matroid import TransversalMatroid
        sage: M = TransversalMatroid(sets=[['a','c']], groundset = ['a', 'c', 'd'])
        sage: M.loops()
        frozenset({'d'})
        sage: M.full_rank()
        1
    """

    def __init__(self, sets, groundset=None, set_labels=None, matching=None):
        """
        See class definition for full documentation.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0,1,2,3], [1,2], [1,3,4]]
            sage: set_labels = [5,6,7]
            sage: M = TransversalMatroid(sets, set_labels=set_labels)
            sage: M.groundset()
            frozenset({0, 1, 2, 3, 4})

        TESTS::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: M = TransversalMatroid([[],[],[]], groundset=range(3)); M
            Transversal matroid of rank 0 on 3 elements, with 0 sets.
        """
        contents = set([e for subset in sets for e in subset])
        if groundset is None:
            groundset = contents
        elif not contents.issubset(groundset):
            raise ValueError("ground set and sets do not match")

        self._sets = Counter([frozenset(s) for s in sets if s])

        element_int_map = {e:i for i, e in enumerate(groundset)}
        int_element_map = {i:e for i, e in enumerate(groundset)}

        # we need a matching and a corresponding graph
        if not set_labels:
            if matching:
                raise ValueError("set labels must be provided if matching is provided")
            set_labels = ['s' + str(i) for i in range(len(sets))]
        else:
            if len(set_labels) != len(sets):
                raise ValueError("set labels do not match sets")
            if not contents.isdisjoint(set_labels):
                raise ValueError("set labels cannot be element labels")
            if matching:
                matching_temp = matching

        if not matching:
            B = BipartiteGraph()
            for e in groundset:
                B.add_vertex(element_int_map[e], left=True)
            for i, s in enumerate(sets):
                new_vertex = set_labels[i]
                for e in s:
                    B.add_edge(new_vertex, element_int_map[e])
            matching_temp = {}
            for u, v, _ in B.matching():
                if u in range(len(groundset)):
                    matching_temp[int_element_map[u]] = v
                else:
                    matching_temp[int_element_map[v]] = u

        self._set_labels = list(set_labels)

        # determine the basis from the matching
        basis = frozenset(matching_temp.keys())

        # This creates self._groundset attribute, among other things
        # It takes the actual ground set labels, not the translated ones
        BasisExchangeMatroid.__init__(self, groundset, basis)

        # matching_temp uses actual ground set labels
        # self._matching will use the translated ones
        self._matching = {element_int_map[e]: matching_temp[e] for e in matching_temp.keys()}

        # Build a DiGraph for doing basis exchange
        self._D = nx.DiGraph()
        # Make sure we get isolated vertices, corresponding to loops
        for v in groundset:
            self._D.add_node(element_int_map[v])

        # For sets in the matching, orient them as starting from the collections
        for u in self._matching.keys():
            self._D.add_edge(self._matching[u], u)

        for i, s in enumerate(sets):
            for e in s:
                if (not (e in matching_temp.keys()) or
                    not (matching_temp[e] == set_labels[i])):
                    self._D.add_edge(element_int_map[e], set_labels[i])

    def digraph(self):
        """
        debugging
        """
        return self._D


    cdef bint __is_exchange_pair(self, long x, long y) except -1:
        r"""
        Check for `M`-alternating path from `x` to `y`.
        """
        if nx.has_path(self._D, y, x):
            return True
        else:
            return False

    cdef int __exchange(self, long x, long y) except -1:
        r"""
        Replace ``self.basis() with ``self.basis() - x + y``. Internal method, does no checks.
        """
        sh = nx.shortest_path(self._D, y, x)
        del self._matching[x]
        for i in range(0, len(sh)-1, 2):
            self._matching[sh[i]] = sh[i+1]

        BasisExchangeMatroid.__exchange(self, x, y)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0, 1, 2, 3]] * 3
            sage: M = TransversalMatroid(sets); M
            Transversal matroid of rank 3 on 4 elements, with 3 sets.
        """
        sets_number = sum(i for i in self._sets.values())
        S = ("Transversal matroid of rank " + str(self.rank()) + " on "
            + str(self.size()) + " elements, with " + str(sets_number)
            + " sets.")
        return S
