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
    sage: M.is_isomorphic(matroids.Uniform(3,6))
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
        sage: sum(1 for b in M.bases())
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
        # update the internal matching
        sh = nx.shortest_path(self._D, y, x)
        del self._matching[x]
        for i in range(0, len(sh)-1, 2):
            self._matching[sh[i]] = sh[i+1]

        # update the graph to reflect this new matching
        sh_edges = []
        sh_edges_r = []
        for i in range(len(sh[:-1])):
            sh_edges.append((sh[i], sh[i+1]))
            sh_edges_r.append((sh[i+1], sh[i]))
        self._D.remove_edges_from(sh_edges)
        self._D.add_edges_from(sh_edges_r)

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

    cpdef sets(self):
        """
        Return the sets of the transversal matroid.

        A transversal matroid can be viewed as a ground set with a collection
        from its powerset. This is represented as a bipartite graph, where
        an edge represents containment.

        OUTPUT:

        A list of lists that correspond to the sets of the transversal matroid.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0,1,2,3], [1,2], [3,4]]
            sage: set_labels = [5,6,7]
            sage: M = TransversalMatroid(sets, set_labels=set_labels)
            sage: sorted(M.sets()) == sorted(sets)
            True
        """
        # Format this in the same way we expect input
        set_list = []
        for s, number in self._sets.iteritems():
            for i in range(number):
                set_list.append(list(s))
        return set_list

    def __richcmp__(left, right, op):
        r"""
        Compare two matroids.

        We take a very restricted view on equality: the objects need to be of
        the exact same type (so no subclassing) and the internal data need to
        be the same. For transversal matroids, in particular, the presentation
        as a bipartite graph must be the same.

        .. WARNING::

            This method is linked to __hash__. If you override one, you MUST override the other!

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [['a','b','c'], ['c','d'], ['a','d'], ['e']]
            sage: M = TransversalMatroid(sets)
            sage: M1 = TransversalMatroid(sets)
            sage: sets2 = [['a','b','c'], ['c','d'], ['a','d','e'], ['e']]
            sage: M2 = TransversalMatroid(sets2)
            sage: M1 == M2
            False
            sage: M1.equals(M2)
            True
        """
        if op not in [Py_EQ, Py_NE]:
            return NotImplemented
        if not isinstance(left, TransversalMatroid) or not isinstance(right, TransversalMatroid):
            return NotImplemented
        if left.__class__ != right.__class__:   # since we have some subclasses, an extra test
            return NotImplemented
        if op == Py_EQ:
            res = True
        if op == Py_NE:
            res = False
        # res gets inverted if matroids are deemed different.
        if (left.groundset() == right.groundset() and
            sorted(left.sets()) == sorted(right.sets())):
            return res
        else:
            return not res

    def __hash__(self):
        r"""
        Return an invariant of the matroid.

        This function is called when matroids are added to a set. It is very
        desirable to override it so it can distinguish matroids on the same
        groundset, which is a very typical use case!

        .. WARNING::

            This method is linked to __richcmp__ (in Cython) and __cmp__ or
            __eq__/__ne__ (in Python). If you override one, you should (and in
            Cython: MUST) override the other!

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets1 = [[0,1,2,3], [1,2], [3,4]]
            sage: M1 = TransversalMatroid(sets1)
            sage: M2 = TransversalMatroid(sets1, set_labels=[5,6,7])
            sage: sets3 = [['a','b','c'], ['c','d'], ['a','d','e']]
            sage: M3 = TransversalMatroid(sets3)
            sage: hash(M1) == hash(M2)
            True
            sage: hash(M1) == hash(M3)
            False
        """
        return hash((self._E, self._sets.iteritems()))

    cpdef __translate_matching(self):
        """
        Return a Python dictionary that can be used as input in __init__().
        """
        matching = {}
        for x in self._matching.keys():
            matching[self._E[x]] = self._matching[x]
        return matching


    def __copy__(self):
        """
        Create a shallow copy.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0,1,2,3], [1,2], [1,3,4]]
            sage: M = TransversalMatroid(sets)
            sage: N = copy(M)  # indirect doctest
            sage: N == M
            True
        """
        cdef TransversalMatroid N
        N = TransversalMatroid(groundset=self._E, sets=self.sets(),
            set_labels=self._set_labels, matching=self.__translate_matching())
        N.rename(getattr(self, '__custom_name'))
        return N

    def __deepcopy__(self, memo):
        """
        Create a deep copy.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: sets = [[0,1,2,3], [1,2], [1,3,4]]
            sage: M = TransversalMatroid(sets)
            sage: N = deepcopy(M)  # indirect doctest
            sage: N == M
            True
        """
        cdef TransversalMatroid N
        N = TransversalMatroid(groundset=deepcopy(self._E, memo), sets=deepcopy(
            self.sets(), memo), set_labels=deepcopy(self._set_labels, memo),
            matching=deepcopy(self.__translate_matching(), memo))
        N.rename(deepcopy(getattr(self, '__custom_name'), memo))
        return N

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        OUTPUT:

        A tuple ``(unpickle, (version, data))``, where ``unpickle`` is the
        name of a function that, when called with ``(version, data)``,
        produces a matroid isomorphic to ``self``. ``version`` is an integer
        (currently 0) and ``data`` is a tuple ``(sets, E, name)`` where
        ``E`` is the groundset of the matroid, ``sets`` is the subsets of the
        transversal, and ``name`` is a custom name.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import *
            sage: sets = [range(6)] * 3
            sage: M = TransversalMatroid(sets)
            sage: M == loads(dumps(M))
            True
            sage: M.rename('U36')
            sage: loads(dumps(M))
            U36
        """
        import sage.matroids.unpickling
        data = (self.sets(), self._E, self._set_labels, self.__translate_matching(),
            getattr(self, '__custom_name'))
        version = 0
        return sage.matroids.unpickling.unpickle_transversal_matroid, (version, data)
