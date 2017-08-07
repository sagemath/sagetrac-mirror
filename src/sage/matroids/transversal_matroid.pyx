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
The input should be a ``BipartiteGraph``. A ``Graph`` may also be used, but in this
case, the ground set must be specified so the method will know the desired
bipartition::

    sage: from sage.matroids.transversal_matroid import *
    sage: G = graphs.CompleteBipartiteGraph(3,6)
    sage: B = BipartiteGraph(G)
    sage: M = TransversalMatroid(B); M
    Transversal matroid of rank 3 on 6 elements, with 3 subsets.
    sage: M.groundset()
    frozenset({3, 4, 5, 6, 7, 8})
    sage: M.is_isomorphic(matroids.Uniform(3,6))
    True
    sage: N = TransversalMatroid(G)
    Traceback (most recent call last):
    ...
    TypeError: ground set must be specified or graph must be BipartiteGraph
    sage: N = TransversalMatroid(G, groundset=[3,4,5,6,7,8])
    sage: N == M
    True

If a ``BipartiteGraph`` is given and both sides have the same cardinality, the left
side is assumed to be the ground set, unless the ground set is specified::

    sage: from sage.matroids.transversal_matroid import *
    sage: B = BipartiteGraph(graphs.CompleteBipartiteGraph(4,4))
    sage: M = TransversalMatroid(B); M.groundset()
    frozenset({0, 1, 2, 3})

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

from sage.matroids.matroid cimport Matroid
from sage.matroids.basis_exchange_matroid cimport BasisExchangeMatroid
from sage.matroids.minor_matroid import MinorMatroid
from sage.matroids.utilities import newlabel

from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from sage.graphs.bipartite_graph import BipartiteGraph

from cpython.object cimport Py_EQ, Py_NE
from copy import copy, deepcopy

cpdef graph_from_buckets(buckets, groundset):
    r"""
    Construct a bipartite graph from sets over a given ground set.

    INPUT:

    - ``buckets`` -- A frozenset consisting of tuples. Each tuple corresponds
      to a vertex on the right side and its neighbors.
    - ``groundset`` -- A list of vertices on the left side.

    OUTPUT:

    A ``BipartiteGraph``.

    EXAMPLES::

        sage: from sage.matroids.transversal_matroid import *
        sage: B = BipartiteGraph(graphs.CompleteBipartiteGraph(4,3))
        sage: M = TransversalMatroid(B)
        sage: buckets = M.sets()
        sage: G = graph_from_buckets(buckets, range(4))
        sage: G == B
        True
    """
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
        sage: B = BipartiteGraph(graphs.CompleteBipartiteGraph(4,3))
        sage: M = TransversalMatroid(B, [0, 1, 2, 3])
        sage: M.full_rank()
        3
        sage: M.bases_count()
        4
        sage: sum(1 for b in M.bases())
        4

    ::

        sage: from sage.matroids.transversal_matroid import TransversalMatroid
        sage: G = Graph([('a', 'b'), ('b', 'c')])
        sage: G.add_vertex('d')
        sage: M = TransversalMatroid(G, groundset=['a','c','d'])
        sage: M.loops()
        frozenset({'d'})
        sage: M.full_rank()
        1
        sage: N = M.delete('d'); N
        Transversal matroid of rank 1 on 2 elements, with 1 subsets.
        sage: N.is_isomorphic(matroids.Uniform(1,2))
        True
    """

    def __init__(self, B, groundset = None):
        """
        See class definition for full documentation.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: edgedict = {5:[0,1,2,3], 6:[1,2], 7:[1,3,4]}
            sage: B = BipartiteGraph(edgedict)
            sage: M = TransversalMatroid(B)
            sage: M.groundset()
            frozenset({0, 1, 2, 3, 4})
            sage: G = Graph(edgedict)
            sage: M1 = TransversalMatroid(G)
            Traceback (most recent call last):
            ...
            TypeError: ground set must be specified or graph must be BipartiteGraph
            sage: N = TransversalMatroid(G, groundset = [0, 5, 1])
            Traceback (most recent call last):
            ...
            ValueError: ground set must specify a bipartition
        """
        # Make this work with a bipartite graph as input
        # In a later ticket, make the constructor work with a collection of sets as input
        if groundset is None:
            if isinstance(B, BipartiteGraph):
                bipartition = B.bipartition()
                if len(bipartition[0]) >= len(bipartition[1]):
                    groundset = bipartition[0]
                else:
                    groundset = bipartition[1]
            else:
                raise TypeError("ground set must be specified or graph must be BipartiteGraph")
        else:
            if not (set(groundset).issubset(set(B.vertices())) and
                len(groundset) == len(set(groundset))):
                raise ValueError("ground set must correspond to vertices")
            if not (B.is_independent_set(groundset) and
                B.is_independent_set(set(B.vertices()).difference(groundset))):
                raise ValueError("ground set must specify a bipartition")

        # put the ground set on the left
        partition = [list(groundset), list(set(B.vertices()).difference(groundset))]

        # put the bucket's name in the bucket, to distinguish between those
        # with the same sets of neighbors
        self._buckets = frozenset([(v, frozenset(B.neighbors(v))) for v in
            B.vertices() if v not in groundset])

        # throw away edge labels
        self._matching = set([(u, v) for (u, v, l) in B.matching()])

        vertices_in_matching = set([u for u, v in self._matching])
        vertices_in_matching.update([v for u, v in self._matching])
        basis = frozenset([v for v in vertices_in_matching if v in groundset])

        # This creates self._groundset attribute, among other things
        BasisExchangeMatroid.__init__(self, groundset, basis)

        # Build a DiGraph for doing basis exchange
        # This will be a simple DiGraph
        self._D = DiGraph()
        # Make sure we get isolated vertices, corresponding to loops
        for v in groundset:
            self._D.add_vertex(v)
        for u, v, l in B.edge_iterator():
            if (u, v) in self._matching:
            # For the edges in our matching, orient them as starting from the collections
                if u in self._groundset:
                    self._D.add_edge(v, u)
                else:
                    self._D.add_edge(u, v)
            else:
            # Otherwise orient them as starting from the ground set
                if u in self._groundset:
                    self._D.add_edge(u, v)
                else:
                    self._D.add_edge(v, u)

    cdef bint __is_exchange_pair(self, long x, long y) except -1:
        r"""
        Check for `M`-alternating path from `x` to `y`.
        """
        e = self._E[x]
        f = self._E[y]
        if self._D.shortest_path(f, e):
            return True
        else:
            return False

    cdef int __exchange(self, long x, long y) except -1:
        r"""
        Replace ``self.basis() with ``self.basis() - x + y``. Internal method, does no checks.
        """
        e = self._E[x]
        f = self._E[y]
        sh = self._D.shortest_path(f, e)
        shortest_path = []
        for i in range(len(sh[:-1])):
            shortest_path.append((sh[i], sh[i+1]))
        self._D.reverse_edges(shortest_path)

        for u, v in shortest_path:
            if (u, v) in self._matching:
                self._matching.remove((u, v))
            else:
                self._matching.add((v, u))

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

    cpdef sets(self):
        """
        Return the sets of the transversal matroid.

        A transversal matroid can be viewed as a ground set with a collection
        from its powerset. This is represented as a bipartite graph, where
        an edge represents containment.

        OUTPUT:

        A frozenset of tuples consisting of a name for the set, and the ground set
        elements it contains.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: edgedict = {5:[0,1,2,3], 6:[1,2], 7:[1,3,4]}
            sage: B = BipartiteGraph(edgedict)
            sage: M = TransversalMatroid(B)
            sage: M.sets()
            frozenset({(5, frozenset({0, 1, 2, 3})),
                       (6, frozenset({1, 2})),
                       (7, frozenset({1, 3, 4}))})
        """
        return self._buckets

    def graph(self):
        """
        Return a bipartite graph representing the transversal matroid.

        The TransversalMatroid object keeps track of a particular correspondence
        between ground set elements and sets as specified by the input. The graph
        returned by this method will reflect this correspondence, as opposed to
        giving a minimal presentation.

        OUTPUT:

        A SageMath graph.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: edgedict = {5:[0,1,2,3], 6:[1,2], 7:[1,3,4]}
            sage: B = BipartiteGraph(edgedict)
            sage: M = TransversalMatroid(B)
            sage: B2 = M.graph()
            sage: B == B2
            True
            sage: B is B2
            False
            sage: edgedict = {5:[0,1,2,3], 6:[1,2], 'a':[1,3,4]}
            sage: B3 = BipartiteGraph(edgedict)
            sage: M = TransversalMatroid(B3)
            sage: B4 = M.graph()
            sage: B4 == B2
            False
            sage: B4.is_isomorphic(B2)
            True
        """
        return graph_from_buckets(self._buckets, self._groundset)

    cpdef _minor(self, contractions, deletions):
        """
        Return a minor.

        Deletions will yield a new transversal matroid. Contractions will have to
        be a MinorMatroid until Gammoid is implemented.

        INPUT:

        - ``contractions`` -- An independent subset of the ground set, as a frozenset.
        - ``deletions`` -- A coindependent subset of the ground set, as a frozenset.

        OUTPUT:

        If ``contractions`` is the empty set, an instance of TransversalMatroid.
        Otherwise, an instance of MinorMatroid.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: edgedict = {5:[0,1,2,3], 6:[1,2], 7:[1,3,4]}
            sage: B = BipartiteGraph(edgedict)
            sage: M1 = TransversalMatroid(B)
            sage: N1 = M1.delete([2,3])
            sage: B.delete_vertices([2,3])
            sage: M2 = TransversalMatroid(B)
            sage: N1 == M2
            True
            sage: M1._minor(deletions=set([3]), contractions=set([4]))
            M / {4}, where M is Transversal matroid of rank 3 on 4 elements, with 3 subsets.

        ::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: G = Graph([('a', 'b'), ('b', 'c')])
            sage: G.add_vertex('d')
            sage: G.add_edge('e','f')
            sage: M = TransversalMatroid(G, groundset=['a','c','d','e'])
            sage: N = M.delete(['d','e']); N
            Transversal matroid of rank 1 on 2 elements, with 1 subsets.
        """
        # if contractions are just coloops, we can just delete them
        if self.corank(contractions) == 0:
            deletions = deletions.union(contractions)
            contractions = set()

        if deletions:
            buckets = set()
            for label, s in self._buckets:
                new_s = frozenset([e for e in s if e not in deletions])
                if new_s:
                # skip over empty buckets, keeping in mind the label is in there
                    buckets.add((label, new_s))
            groundset = self._groundset.difference(deletions)

            N = TransversalMatroid(graph_from_buckets(buckets, groundset))
            # Check if what remains is just coloops
            return N.contract(contractions)
        else:
            N = self

        if contractions:
            return MinorMatroid(N, contractions=contractions, deletions=set())
        else:
            return N

    cpdef transversal_extension(self, element=None, newset=False, sets=[]):
        r"""
        Return a TransversalMatroid extended by an element.

        This will modify the presentation of the transversal matroid by adding
        a new element, and placing this element in the specified sets. It is also
        possible to use this method to create a new set which will have the new
        element as its only member, making it a coloop.

        INPUT:

        - ``element`` -- (optional) The name for the new element.
        - ``newset`` -- (optional) If specified, the element will be
          given its own set. If ``True``, a name will be generated; otherwise
          this value will be used. This will make the element
          a coloop.
        - ``sets`` -- (default: ``None``) An iterable of labels representing the
          sets in the current presentation that the new element will belong to.

        OUTPUT:

        A TransversalMatroid with a ground set element added to specified sets.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: G = Graph([('a', 'b'), ('b', 'c')])
            sage: M = TransversalMatroid(G, groundset=['a','c'])
            sage: M1 = M.transversal_extension(element='d', newset='e')
            sage: M2 = M.transversal_extension(element='d', newset=True)
            sage: M1.coloops()
            frozenset({'d'})
            sage: True in M2.graph().vertices()
            False
            sage: M1.is_isomorphic(M2)
            True
            sage: M3 = M.transversal_extension('d', sets=['b'])
            sage: M3.is_isomorphic(matroids.Uniform(1,3))
            True
            sage: M4 = M.transversal_extension('d', sets=['a'])
            Traceback (most recent call last):
            ...
            ValueError: sets do not match presentation
            sage: M4 = M.transversal_extension('a', sets=['b'])
            Traceback (most recent call last):
            ...
            ValueError: cannot extend by element already in ground set
            sage: M2.transversal_extension(newset='b')
            Traceback (most recent call last):
            ...
            ValueError: newset is already a vertex in the presentation
            sage: M5 = M1.transversal_extension()
            sage: len(M5.loops())
            1
            sage: M2.transversal_extension(element='b')
            Transversal matroid of rank 2 on 4 elements, with 2 subsets.

        ::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: edgedict = {5:[0,1,2,3], 6:[1,2], 7:[1,3,4]}
            sage: B = BipartiteGraph(edgedict)
            sage: M = TransversalMatroid(B)
            sage: N = M.delete(2)
            sage: M1 = N.transversal_extension(element=2, sets=[5,6])
            sage: M1 == M
            True

        ::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: B = BipartiteGraph(graphs.CompleteBipartiteGraph(4,3))
            sage: M = TransversalMatroid(B)
            sage: N = M.transversal_extension(element='a', newset=True, sets=[4])
            sage: N.graph().degree('a')
            2
        """
        sets = set(sets)
        if element is None:
            element = newlabel(self.groundset())
        elif element in self._groundset:
            raise ValueError("cannot extend by element already in ground set")
        labels = [label for label, s in self._buckets]
        if not sets.issubset(labels):
            raise ValueError("sets do not match presentation")

        # If the new element has the same name as a bucket label,
        # ground set labels are more important, so we should change the
        # bucket's label
        labels_map = {l: l for l in labels}
        if element in labels:
            new_label = newlabel(self.groundset().union(labels))
            labels_map[element] = new_label

        # newset should not be a ground set element or existing set
        if newset in self._D:
            # keywords `True` and `False` give us problems here
            if newset is not False and newset is not True:
                raise ValueError("newset is already a vertex in the presentation")

        buckets = set()
        for l, s in self._buckets:
            if l in sets:
                buckets.add((labels_map[l], frozenset(s.union([element]))))
            else:
                buckets.add((labels_map[l], s))

        if newset:
            if newset is True:
                newset = newlabel(self.groundset().union(labels))
            buckets.add((newset, frozenset([element])))

        groundset = self.groundset().union([element])

        return TransversalMatroid(graph_from_buckets(buckets, groundset))

    def transversal_extensions(self, element=None, sets=[]):
        r"""
        Return an iterator of extensions based on the transversal presentation.

        This method will take an extension by adding an element to every possible
        sub-collection of the collection of desired sets. No checking is done
        for equal matroids. It is advised to make sure the presentation has as
        few sets as possible by using
        :meth:`reduce_presentation() <sage.matroids.transversal_matroid.TransversalMatroid.reduce_presentation>`

        INPUT:

        - ``element`` -- (optional) The name of the new element.
        - ``sets`` -- (optional) a list containing names of sets in the matroid's
          presentation. If not specified, every set will be used.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: M = TransversalMatroid(BipartiteGraph(graphs.CompleteBipartiteGraph(3,4)))
            sage: len([N for N in M.transversal_extensions()])
            8
            sage: len([N for N in M.transversal_extensions(sets=[0,1])])
            4
            sage: set(sorted([N.groundset() for N in M.transversal_extensions(element=7)]))
            {frozenset({3, 4, 5, 6, 7})}
        """
        if element is None:
            element = newlabel(self.groundset())
        elif element in self._groundset:
            raise ValueError("cannot extend by element already in ground set")

        labels = [label for label, s in self._buckets]
        if not sets:
            sets = labels
        elif not set(sets).issubset(labels):
            raise ValueError("sets do not match presentation")

        # Adapted from the Python documentation
        from itertools import chain, combinations
        sets_list = list(sets)
        powerset = chain.from_iterable(combinations(sets_list, r) for r in range(len(sets_list)+1))

        for collection in powerset:
            yield self.transversal_extension(element=element, sets=collection)

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
            sage: B = BipartiteGraph(graphs.CubeGraph(3))
            sage: M = TransversalMatroid(B)
            sage: N = TransversalMatroid(M.graph())
            sage: O = TransversalMatroid(BipartiteGraph(graphs.CompleteBipartiteGraph(4,3)))
            sage: M == N
            True
            sage: M == O
            False
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
        if left.groundset() == right.groundset() and left.sets() == right.sets():
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
            sage: B = BipartiteGraph(graphs.CubeGraph(3))
            sage: M = TransversalMatroid(B)
            sage: N = TransversalMatroid(M.graph())
            sage: O = TransversalMatroid(BipartiteGraph(graphs.CompleteBipartiteGraph(4,3)))
            sage: hash(M) == hash(N)
            True
            sage: hash(M) == hash(O)
            False
        """
        return hash((self.groundset(), self.sets()))

    cpdef reduce_presentation(self):
        """
        Return an equal transversal matroid where the number of sets equals the rank.

        Every transversal matroid `M` has a presentation with `r(M)` sets, and if `M`
        has no coloops, then every presentation has `r(M)` sets. This method
        discards extra sets if `M` has coloops.

        INPUT:

        None.

        OUTPUT:

        A ``TransversalMatroid`` instance with a reduced presentation.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: edgelist = [(0, 'a'), (1, 'a'), (2, 'b'), (2, 'c')]
            sage: M = TransversalMatroid(Graph(edgelist), groundset = [0,1,2]); M
            Transversal matroid of rank 2 on 3 elements, with 3 subsets.
            sage: N = M.reduce_presentation(); N
            Transversal matroid of rank 2 on 3 elements, with 2 subsets.
            sage: N.equals(M)
            True
            sage: N == M
            False

        ::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: M = TransversalMatroid(BipartiteGraph(graphs.CompleteBipartiteGraph(3,4)))
            sage: N = M.reduce_presentation()
            sage: M == N
            True
        """
        if len(self.sets()) == self.full_rank():
            return self
        else:
            coloops = self.coloops()
            coloops_to_delete = [e for e in coloops if self._D.degree(e) > 1]
            N = self.delete(coloops_to_delete)
            buckets = set(N.sets())
            for c in coloops:
                buckets.add((newlabel(self._D.vertices()), frozenset([c])))
            return TransversalMatroid(graph_from_buckets(buckets, self.groundset()))

    def __copy__(self):
        """
        Create a shallow copy.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: edgedict = {5:[0,1,2,3], 6:[1,2], 7:[1,3,4]}
            sage: M = TransversalMatroid(BipartiteGraph(edgedict))
            sage: N = copy(M)  # indirect doctest
            sage: N == M
            True
        """
        cdef TransversalMatroid N
        N = TransversalMatroid(groundset=self._E, B=self._D.to_undirected())
        N.rename(getattr(self, '__custom_name'))
        return N

    def __deepcopy__(self, memo):
        """
        Create a deep copy.

        EXAMPLES::

            sage: from sage.matroids.transversal_matroid import TransversalMatroid
            sage: edgedict = {5:[0,1,2,3], 6:[1,2], 7:[1,3,4]}
            sage: M = TransversalMatroid(BipartiteGraph(edgedict))
            sage: N = deepcopy(M)  # indirect doctest
            sage: M == N
            True
        """
        cdef TransversalMatroid N
        N = TransversalMatroid(groundset=deepcopy(self._E, memo), B=deepcopy(
            self._D.to_undirected(), memo))
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
            sage: B = BipartiteGraph(graphs.CompleteBipartiteGraph(6,3))
            sage: M = TransversalMatroid(B)
            sage: M == loads(dumps(M))
            True
            sage: M.rename('U36')
            sage: loads(dumps(M))
            U36
        """
        import sage.matroids.unpickling
        data = (self.sets(), self._E, getattr(self, '__custom_name'))
        version = 0
        return sage.matroids.unpickling.unpickle_transversal_matroid, (version, data)
