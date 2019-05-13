r"""
Spanning trees

This module is a collection of algorithms on spanning trees. Also included in
the collection are algorithms for minimum spanning trees. See the book
[JoynerNguyenCohen2010]_ for descriptions of spanning tree algorithms,
including minimum spanning trees.

.. SEEALSO::

    - :meth:`GenericGraph.min_spanning_tree
      <sage.graphs.generic_graph.GenericGraph.min_spanning_tree>`.

.. TODO::

    - Rewrite :func:`kruskal` to use priority queues.
    - Parallel version of Boruvka's algorithm.
    - Randomized spanning tree construction.

REFERENCES:

.. [Aldous90] \D. Aldous, *The random walk construction of
  uniform spanning trees*, SIAM J Discrete Math 3 (1990),
  450-465.

.. [Broder89] \A. Broder, *Generating random spanning trees*,
  Proceedings of the 30th IEEE Symposium on Foundations of
  Computer Science, 1989, pp. 442-447. :doi:`10.1109/SFCS.1989.63516`,
  <http://www.cs.cmu.edu/~15859n/RelatedWork/Broder-GenRanSpanningTrees.pdf>_

.. [CormenEtAl2001] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest,
  and Clifford Stein. *Introduction to Algorithms*. 2nd edition, The MIT Press,
  2001.

.. [GoodrichTamassia2001] Michael T. Goodrich and Roberto Tamassia.
  *Data Structures and Algorithms in Java*. 2nd edition, John Wiley & Sons,
  2001.

.. [JoynerNguyenCohen2010] David Joyner, Minh Van Nguyen, and Nathann Cohen.
  *Algorithmic Graph Theory*. 2010,
  http://code.google.com/p/graph-theory-algorithms-book/

.. [Sahni2000] Sartaj Sahni. *Data Structures, Algorithms, and Applications
  in Java*. McGraw-Hill, 2000.


Methods
-------
"""

# ****************************************************************************
#       Copyright (c) 2007 Jason Grout <jason-sage@creativetrax.com>
#       Copyright (c) 2009 Mike Hansen <mhansen@gmail.com>
#       Copyright (c) 2010 Gregory McWhirter <gmcwhirt@uci.edu>
#       Copyright (c) 2010 Minh Van Nguyen <nguyenminh2@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import absolute_import

cimport cython
from sage.ext.memory_allocator cimport MemoryAllocator
from sage.sets.disjoint_set cimport DisjointSet_of_hashables


cpdef kruskal(G, wfunction=None, bint check=False):
    r"""
    Minimum spanning tree using Kruskal's algorithm.

    This function assumes that we can only compute minimum spanning trees for
    undirected graphs. Such graphs can be weighted or unweighted, and they can
    have multiple edges (since we are computing the minimum spanning tree, only
    the minimum weight among all `(u,v)`-edges is considered, for each pair
    of vertices `u`, `v`).

    INPUT:

    - ``G`` -- an undirected graph.

    - ``weight_function`` -- function (default: ``None``); a function that
      inputs an edge ``e`` and outputs its weight. An edge has the form
      ``(u,v,l)``, where ``u`` and ``v`` are vertices, ``l`` is a label (that
      can be of any kind).  The ``weight_function`` can be used to transform the
      label into a weight. In particular:

      - if ``weight_function`` is not ``None``, the weight of an edge ``e``
        is ``weight_function(e)``;

      - if ``weight_function`` is ``None`` (default) and ``g`` is weighted
        (that is, ``g.weighted()==True``), the weight of an edge
        ``e=(u,v,l)`` is ``l``, independently on which kind of object ``l``
        is: the ordering of labels relies on Python's operator ``<``;

      - if ``weight_function`` is ``None`` and ``g`` is not weighted, we set
        all weights to 1 (hence, the output can be any spanning tree).

    - ``check`` -- boolean (default: ``False``); whether to first perform sanity
      checks on the input graph ``G``. Default: ``check=False``. If we toggle
      ``check=True``, the following sanity checks are first performed on ``G``
      prior to running Kruskal's algorithm on that input graph:

      - Is ``G`` the null graph?
      - Is ``G`` disconnected?
      - Is ``G`` a tree?
      - Does ``G`` have self-loops?
      - Does ``G`` have multiple edges?

      By default, we turn off the sanity checks for performance reasons. This
      means that by default the function assumes that its input graph is
      connected, and has at least one vertex. Otherwise, you should set
      ``check=True`` to perform some sanity checks and preprocessing on the
      input graph. If ``G`` has multiple edges or self-loops, the algorithm
      still works, but the running-time can be improved if these edges are
      removed. To further improve the runtime of this function, you should call
      it directly instead of using it indirectly via
      :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`.

    OUTPUT:

    The edges of a minimum spanning tree of ``G``, if one exists, otherwise
    returns the empty list.

    .. SEEALSO::

        - :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`
        - :func:`kruskal_iterator`
        - :func:`filter_kruskal` and :func:`filter_kruskal_iterator`

    EXAMPLES:

    An example from pages 727--728 in [Sahni2000]_. ::

        sage: from sage.graphs.spanning_tree import kruskal
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: E = kruskal(G, check=True); E
        [(1, 6, 10), (3, 4, 12), (2, 7, 14), (2, 3, 16), (4, 5, 22), (5, 6, 25)]

    Variants of the previous example. ::

        sage: H = Graph(G.edges(labels=False))
        sage: kruskal(H, check=True)
        [(1, 2, None), (1, 6, None), (2, 3, None), (2, 7, None), (3, 4, None), (4, 5, None)]
        sage: G.allow_loops(True)
        sage: G.allow_multiple_edges(True)
        sage: G
        Looped multi-graph on 7 vertices
        sage: for i in range(20):
        ....:     u = randint(1, 7)
        ....:     v = randint(1, 7)
        ....:     w = randint(0, 20)
        ....:     G.add_edge(u, v, w)
        sage: H = copy(G)
        sage: H
        Looped multi-graph on 7 vertices
        sage: def sanitize(G):
        ....:     G.allow_loops(False)
        ....:     G.allow_multiple_edges(False, keep_label='min')
        sage: sanitize(H)
        sage: H
        Graph on 7 vertices
        sage: sum(e[2] for e in kruskal(G, check=True)) == sum(e[2] for e in kruskal(H, check=True))
        True

    An example from pages 599--601 in [GoodrichTamassia2001]_. ::

        sage: G = Graph({"SFO":{"BOS":2704, "ORD":1846, "DFW":1464, "LAX":337},
        ....: "BOS":{"ORD":867, "JFK":187, "MIA":1258},
        ....: "ORD":{"PVD":849, "JFK":740, "BWI":621, "DFW":802},
        ....: "DFW":{"JFK":1391, "MIA":1121, "LAX":1235},
        ....: "LAX":{"MIA":2342},
        ....: "PVD":{"JFK":144},
        ....: "JFK":{"MIA":1090, "BWI":184},
        ....: "BWI":{"MIA":946}})
        sage: G.weighted(True)
        sage: kruskal(G, check=True)
        [('JFK', 'PVD', 144), ('BWI', 'JFK', 184), ('BOS', 'JFK', 187), ('LAX', 'SFO', 337), ('BWI', 'ORD', 621), ('DFW', 'ORD', 802), ('BWI', 'MIA', 946), ('DFW', 'LAX', 1235)]

    An example from pages 568--569 in [CormenEtAl2001]_. ::

        sage: G = Graph({"a":{"b":4, "h":8}, "b":{"c":8, "h":11},
        ....: "c":{"d":7, "f":4, "i":2}, "d":{"e":9, "f":14},
        ....: "e":{"f":10}, "f":{"g":2}, "g":{"h":1, "i":6}, "h":{"i":7}})
        sage: G.weighted(True)
        sage: T = Graph(kruskal(G, check=True), format='list_of_edges')
        sage: sum(T.edge_labels())
        37
        sage: T.is_tree()
        True

    An example with custom edge labels::

        sage: G = Graph([[0,1,1],[1,2,1],[2,0,10]], weighted=True)
        sage: weight = lambda e:3-e[0]-e[1]
        sage: sorted(kruskal(G, check=True))
        [(0, 1, 1), (1, 2, 1)]
        sage: sorted(kruskal(G, wfunction=weight, check=True))
        [(0, 2, 10), (1, 2, 1)]
        sage: sorted(kruskal(G, wfunction=weight, check=False))
        [(0, 2, 10), (1, 2, 1)]

    TESTS:

    The input graph must not be empty. ::

        sage: from sage.graphs.spanning_tree import kruskal
        sage: kruskal(graphs.EmptyGraph(), check=True)
        []
        sage: kruskal(Graph(), check=True)
        []
        sage: kruskal(Graph(multiedges=True), check=True)
        []
        sage: kruskal(Graph(loops=True), check=True)
        []
        sage: kruskal(Graph(multiedges=True, loops=True), check=True)
        []

    The input graph must be connected. ::

        sage: def my_disconnected_graph(n, ntries, directed=False, multiedges=False, loops=False):
        ....:     G = Graph()
        ....:     k = randint(1, n)
        ....:     G.add_vertices(range(k))
        ....:     if directed:
        ....:         G = G.to_directed()
        ....:     if multiedges:
        ....:         G.allow_multiple_edges(True)
        ....:     if loops:
        ....:         G.allow_loops(True)
        ....:     for i in range(ntries):
        ....:         u = randint(0, k-1)
        ....:         v = randint(0, k-1)
        ....:         if u != v or loops:
        ....:             G.add_edge(u, v)
        ....:     while G.is_connected():
        ....:         u = randint(0, k-1)
        ....:         v = randint(0, k-1)
        ....:         G.delete_edge(u, v)
        ....:     return G
        sage: G = my_disconnected_graph(100, 50, directed=False, multiedges=False, loops=False)  # long time
        sage: kruskal(G, check=True)  # long time
        []
        sage: G = my_disconnected_graph(100, 50, directed=False, multiedges=True, loops=False)  # long time
        sage: kruskal(G, check=True)  # long time
        []
        sage: G = my_disconnected_graph(100, 50, directed=False, multiedges=True, loops=True)  # long time
        sage: kruskal(G, check=True)  # long time
        []

    If the input graph is a tree, then return its edges::

        sage: T = graphs.RandomTree(randint(1, 50))  # long time
        sage: sorted(T.edge_iterator()) == kruskal(T, check=True)  # long time
        True

    If the input is not a Graph::

        sage: kruskal("I am not a graph")
        Traceback (most recent call last):
        ...
        ValueError: The input G must be an undirected graph.
        sage: kruskal(digraphs.Path(10))
        Traceback (most recent call last):
        ...
        ValueError: The input G must be an undirected graph.
    """
    return list(kruskal_iterator(G, wfunction=wfunction, check=check))


def kruskal_iterator(G, wfunction=None, bint check=False):
    """
    Return an iterator implementation of Kruskal algorithm.

    OUTPUT:

    The edges of a minimum spanning tree of ``G``, one by one.

    .. SEEALSO:: :func:`kruskal`

    EXAMPLES::

        sage: from sage.graphs.spanning_tree import kruskal_iterator
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: next(kruskal_iterator(G, check=True))
        (1, 6, 10)
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("The input G must be an undirected graph.")

    # sanity checks
    if check:
        if not G.order():
            return
        if not G.is_connected():
            return
        # G is now assumed to be a nonempty connected graph
        if G.num_verts() == G.num_edges() + 1:
            # G is a tree
            yield from G.edge_iterator()
            return
        g = G.to_simple(to_undirected=False, keep_label='min')
    else:
        g = G

    cdef DisjointSet_of_hashables union_find = DisjointSet_of_hashables(g)
    yield from kruskal_iterator_from_edges(g.edge_iterator(), union_find,
                                           weighted=G.weighted(),
                                           weight_function=wfunction)

def kruskal_iterator_from_edges(edges, union_find, weighted=False, weight_function=None):
    """
    Return an iterator implementation of Kruskal algorithm on list of edges.

    INPUT:

    - ``edges`` -- list of edges

    - ``union_find`` -- a
      :class:`~sage.sets.disjoint_set.DisjointSet_of_hashables` encoding a
      forest

    - ``weighted`` -- boolean (default: ``False``); whether edges are weighted,
      i.e., the label of an edge is a weight

    - ``weight_function`` -- function (default: ``None``); a function that
      inputs an edge ``e`` and outputs its weight. See :func:`kruskal` for more
      details.

    OUTPUT:

    The edges of a minimum spanning tree of ``G``, one by one.

    .. SEEALSO::

        - :func:`kruskal`
        - :func:`filter_kruskal`

    EXAMPLES::

        sage: from sage.graphs.spanning_tree import kruskal_iterator_from_edges
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: union_set=DisjointSet(G.order())
        sage: next(kruskal_iterator_from_edges(G.edges(sort=False), union_set, weighted=G.weighted()))
        (1, 6, 10)
    """
    # We sort edges, as specified.
    if weight_function is None:
        if weighted:
            from operator import itemgetter
            edges = sorted(edges, key=itemgetter(2))
    else:
        edges = sorted(edges, key=weight_function)
    # Kruskal's algorithm
    for e in edges:
         # acyclic test via union-find
         u = union_find.find(e[0])
         v = union_find.find(e[1])
         if u != v:
             yield e
             # merge the trees
             union_find.union(u, v)
             if union_find.number_of_subsets() == 1:
                 return

def filter_kruskal(G, threshold=10000, weight_function=None, bint check=False):
    """
    Minimum spanning tree using Filter Kruskal algorithm.

    This function implements the variant of Kruskal's algorithm proposed in
    [OSS2009]_. Instead of directly sorting the whole set of edges, it
    partitions it in a similar way to quicksort and filter out edges that
    connect vertices of the same tree to reduce the cost of sorting.

    This function assumes that we can only compute minimum spanning trees for
    undirected graphs. Such graphs can be weighted or unweighted, and they can
    have multiple edges (since we are computing the minimum spanning tree, only
    the minimum weight among all `(u,v)`-edges is considered, for each pair of
    vertices `u`, `v`).

    INPUT:

    - ``G`` -- an undirected graph

    - ``weight_function`` -- function (default: ``None``); a function that
      inputs an edge ``e`` and outputs its weight. An edge has the form
      ``(u,v,l)``, where ``u`` and ``v`` are vertices, ``l`` is a label (that
      can be of any kind). The ``weight_function`` can be used to transform the
      label into a weight. In particular:

      - if ``weight_function`` is not ``None``, the weight of an edge ``e``
        is ``weight_function(e)``;

      - if ``weight_function`` is ``None`` (default) and ``g`` is weighted
        (that is, ``g.weighted()==True``), the weight of an edge
        ``e=(u,v,l)`` is ``l``, independently on which kind of object ``l``
        is: the ordering of labels relies on Python's operator ``<``;

      - if ``weight_function`` is ``None`` and ``g`` is not weighted, we set
        all weights to 1 (hence, the output can be any spanning tree).

    - ``threshold`` -- integer (default: 10000); maximum number of edges on
       which to run kruskal algorithm. Above that value, edges are partitioned
       into sets of size at most ``threshold``

    - ``check`` -- boolean (default: ``False``); whether to first perform sanity
      checks on the input graph ``G``. Default: ``check=False``. If we toggle
      ``check=True``, the following sanity checks are first performed on ``G``
      prior to running Kruskal's algorithm on that input graph:

      - Is ``G`` the null graph?
      - Is ``G`` disconnected?
      - Is ``G`` a tree?
      - Does ``G`` have self-loops?
      - Does ``G`` have multiple edges?

    OUTPUT:

    The edges of a minimum spanning tree of ``G``, if one exists, otherwise
    returns the empty list.

    .. SEEALSO::

        - :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`
        - :wikipedia:`Kruskal's_algorithm`
        - :func:`kruskal`
        - :func:`filter_kruskal_iterator`

    EXAMPLES::

        sage: from sage.graphs.spanning_tree import filter_kruskal
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: filter_kruskal(G, check=True)
        [(1, 6, 10), (3, 4, 12), (2, 7, 14), (2, 3, 16), (4, 5, 22), (5, 6, 25)]

        sage: filter_kruskal(Graph(2), check=True)
        []
    """
    return list(filter_kruskal_iterator(G, threshold=threshold, weight_function=weight_function, check=check))

def filter_kruskal_iterator(G, threshold=10000, weight_function=None, bint check=False):
    r"""
    Return an iterator implementation of Filter Kruskal's algorithm.

    OUTPUT:

    The edges of a minimum spanning tree of ``G``, one by one.

    .. SEEALSO::

        - :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`
        - :wikipedia:`Kruskal's_algorithm`
        - :func:`kruskal`
        - :func:`filter_kruskal`

    EXAMPLES:

    The edges of a minimum spanning tree of ``G``, if one exists, otherwise
    returns the empty list.

        sage: from sage.graphs.spanning_tree import filter_kruskal_iterator
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: list(filter_kruskal_iterator(G, threshold=3, check=True))
        [(1, 6, 10), (3, 4, 12), (2, 7, 14), (2, 3, 16), (4, 5, 22), (5, 6, 25)]

    The weights of the spanning trees returned by :func:`kruskal_iterator` and
    :func:`filter_kruskal_iterator` are the same::

        sage: from sage.graphs.spanning_tree import kruskal_iterator
        sage: G = graphs.RandomBarabasiAlbert(50, 2)
        sage: for u, v in G.edge_iterator(labels=False):
        ....:     G.set_edge_label(u, v, randint(1, 10))
        sage: G.weighted(True)
        sage: sum(e[2] for e in kruskal_iterator(G)) == sum(e[2] for e in filter_kruskal_iterator(G, threshold=20))
        True

    TESTS:

    The threshold must be at least 1::

        sage: from sage.graphs.spanning_tree import filter_kruskal_iterator
        sage: next(filter_kruskal_iterator(Graph(), threshold=0))
        Traceback (most recent call last):
        ...
        ValueError: the threshold mut be at least 1

    Check that a threshold of 1 is accepted::

        sage: len(list(filter_kruskal_iterator(graphs.HouseGraph(), threshold=1)))
        4
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the input graph must be undirected")
    if threshold < 1:
        raise ValueError("the threshold mut be at least 1")
    if check:
        if not G.order() or not G.is_connected():
            return
        # G is now assumed to be a nonempty connected graph
        if G.order() == G.size() + 1:
            # G is a tree
            yield from G.edge_iterator()
            return

        g = G.to_simple(to_undirected=False, keep_label='min')
    else:
        g = G

    cdef int m = g.size()
    if m <= threshold:
        yield from kruskal_iterator_from_edges(g.edge_iterator(),
                                               DisjointSet_of_hashables(g),
                                               weighted=G.weighted(),
                                               weight_function=weight_function)
        return

    #
    # Initialize some data structure
    #
    cdef list edges = list(g.edge_iterator())
    # Precompute edge weights to avoid frequent calls to weight_function
    cdef list weight
    if weight_function is None:
        if G.weighted():
            weight = [e[2] for e in edges]
        else:
            weight = [1 for _ in range(m)]
    else:
        weight = [weight_function(e) for e in edges]

    cdef MemoryAllocator mem = MemoryAllocator()
    # Array storing a permutation of the edges.
    # e_index[i] is the position of edge i in list edges
    cdef int* e_index = <int*> mem.allocarray(m, sizeof(int))
    cdef int i, j
    for i in range(m):
        e_index[i] = i
    # Stack of range of edge partitions
    cdef list stack = [(0, m - 1)]
    cdef int begin, end
    # Parameter to  equally divide edges with weight equal the to pivot
    cdef bint ch = True
    # Data structure to record the vertices in each tree of the forest
    cdef DisjointSet_of_hashables union_find = DisjointSet_of_hashables(g)

    #
    # Iteratively partition the list of edges
    #
    while stack:
        begin, end = stack.pop()

        if end - begin < threshold:
            # Filter edges connecting vertices of a same tree
            L = [edges[e_index[i]] for i in range(begin, end + 1)
                 if union_find.find(edges[e_index[i]][0]) != union_find.find(edges[e_index[i]][1])]
            yield from kruskal_iterator_from_edges(L, union_find,
                                                   weighted=G.weighted(),
                                                   weight_function=weight_function)
            if union_find.number_of_subsets() == 1:
                return
            continue

        # Choose a pivot
        pivot = weight[e_index[(begin + end) // 2]]

        # Partition edges with respect to pivot, as in quicksort
        i, j = begin, end
        while i < j:
            while weight[e_index[i]] < pivot and i < j:
                i += 1
            if ch and weight[e_index[i]] == pivot and i < j:
                i += 1
                ch = False
                continue
            while weight[e_index[j]] > pivot and i < j:
                j -= 1
            if not ch and weight[e_index[j]] == pivot and i < j:
                j -= 1
                ch = True
                continue
            if i < j:
                e_index[i], e_index[j] = e_index[j], e_index[i]

        # Record range of edge partitions
        if weight[e_index[i]] <= pivot:
            stack.append((i + 1, end))
            stack.append((begin, i))
        else:
            stack.append((i, end))
            stack.append((begin, i - 1))


cpdef boruvka(G, wfunction=None, bint check=False, bint by_weight=True):
    r"""
    Minimum spanning tree using Boruvka's algorithm.

    This function assumes that we can only compute minimum spanning trees for
    undirected graphs. Such graphs can be weighted or unweighted, and they can
    have multiple edges (since we are computing the minimum spanning tree, only
    the minimum weight among all `(u,v)`-edges is considered, for each pair of
    vertices `u`, `v`).

    INPUT:

    - ``G`` -- an undirected graph.

    - ``wfunction`` -- weight function (default: ``None``); a function that
      inputs an edge ``e`` and outputs its weight. An edge has the form
      ``(u,v,l)``, where ``u`` and ``v`` are vertices, ``l`` is a label (that
      can be of any kind).  The ``wfunction`` can be used to transform the label
      into a weight. In particular:

      - if ``wfunction`` is not ``None``, the weight of an edge ``e`` is
        ``wfunction(e)``;

      - if ``wfunction`` is ``None`` (default) and ``g`` is weighted (that is,
        ``g.weighted()==True``), the weight of an edge ``e=(u,v,l)`` is ``l``,
        independently on which kind of object ``l`` is: the ordering of labels
        relies on Python's operator ``<``;

      - if ``wfunction`` is ``None`` and ``g`` is not weighted, we set all
        weights to 1 (hence, the output can be any spanning tree).

    - ``check`` -- boolean (default: ``False``); whether to first perform sanity
      checks on the input graph ``G``. Default: ``check=False``. If we toggle
      ``check=True``, the following sanity checks are first performed on ``G``
      prior to running Boruvka's algorithm on that input graph:

      - Is ``G`` the null graph or graph on one vertex?
      - Is ``G`` disconnected?
      - Is ``G`` a tree?

      By default, we turn off the sanity checks for performance reasons. This
      means that by default the function assumes that its input graph is
      connected, and has at least one vertex. Otherwise, you should set
      ``check=True`` to perform some sanity checks and preprocessing on the
      input graph.
    
    - ``by_weight`` -- boolean (default: ``False``); whether to find MST by
      using weights of edges provided.  Default: ``by_weight=True``. If
      ``wfunction`` is given, MST is calculated using the weights of edges as
      per the function. If ``wfunction`` is ``None``, the weight of an edge
      ``e=(u,v,l)`` is ``l`` if graph is weighted, or all edge weights are
      considered ``1`` if graph is unweighted. If we toggle ``by_weight=False``,
      all weights are considered as ``1`` and MST is calculated.

    OUTPUT:

    The edges of a minimum spanning tree of ``G``, if one exists, otherwise
    returns the empty list.

    .. SEEALSO::

        - :meth:`~sage.graphs.generic_graph.GenericGraph.min_spanning_tree`

    EXAMPLES:

    An example from pages 727--728 in [Sahni2000]_::

        sage: from sage.graphs.spanning_tree import boruvka
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: E = boruvka(G, check=True); E
        [(1, 6, 10), (2, 7, 14), (3, 4, 12), (4, 5, 22), (5, 6, 25), (2, 3, 16)]
        sage: boruvka(G, by_weight=True)
        [(1, 6, 10), (2, 7, 14), (3, 4, 12), (4, 5, 22), (5, 6, 25), (2, 3, 16)]
        sage: sorted(boruvka(G, by_weight=False))
        [(1, 2, 28), (1, 6, 10), (2, 3, 16), (2, 7, 14), (3, 4, 12), (4, 5, 22)]

    An example with custom edge labels::

        sage: G = Graph([[0,1,1],[1,2,1],[2,0,10]], weighted=True)
        sage: weight = lambda e:3-e[0]-e[1]
        sage: boruvka(G, wfunction=lambda e:3-e[0]-e[1], by_weight=True)
        [(0, 2, 10), (1, 2, 1)]
        sage: boruvka(G, wfunction=lambda e:float(1/e[2]), by_weight=True)
        [(0, 2, 10), (0, 1, 1)]

    An example of disconnected graph with ``check`` disabled::

        sage: from sage.graphs.spanning_tree import boruvka
        sage: G = Graph({1:{2:28}, 3:{4:16}}, weighted=True)
        sage: boruvka(G, check=False)
        []

    TESTS:
    
    If the input graph is a tree, then return its edges::

        sage: T = graphs.RandomTree(randint(1, 10))
        sage: T.edges() == sorted(boruvka(T, check=True))
        True

    Check if the weight of MST returned by Prim's and Boruvka's is the same::

        sage: G = Graph([(u,v,randint(1,5)) for u,v in graphs.CompleteGraph(4).edges(labels=0)], weighted=True)
        sage: G.weighted()
        True
        sage: E1 = G.min_spanning_tree(algorithm='Boruvka')
        sage: E2 = G.min_spanning_tree(algorithm='Prim_Boost')
        sage: sum(e[2] for e in E1) == sum(e[2] for e in E2)
        True

    If the input is not a Graph::

        sage: boruvka("I am not a graph")
        Traceback (most recent call last):
        ...
        ValueError: the input graph must be undirected
        sage: boruvka(digraphs.Path(10))
        Traceback (most recent call last):
        ...
        ValueError: the input graph must be undirected
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the input graph must be undirected")

    if G.order() <= 1:
        return []

    # sanity checks
    if check:
        if not G.is_connected():
            return []
        # G is now assumed to be a nonempty connected graph
        if G.num_verts() == G.num_edges() + 1:
            # G is a tree
            return G.edges(sort=False)

    # Boruvka's algorithm

    # Store the list of active edges as (e, e_weight) in a list
    if by_weight:
        if wfunction is None:
            if G.weighted():
                edge_list = [(e, e[2]) for e in G.edge_iterator()]
            else:
                edge_list = [(e, 1) for e in G.edge_iterator()]
        else:
            edge_list = [(e, wfunction(e)) for e in G.edge_iterator()]
    else:
        edge_list = [(e, 1) for e in G.edge_iterator()]

    # initially, each vertex is a connected component
    cdef DisjointSet_of_hashables partitions = DisjointSet_of_hashables(G.vertex_iterator())
    # a dictionary to store the least weight outgoing edge for each component
    cdef dict cheapest = {}
    cdef list T = []  # stores the edges in minimum spanning tree
    cdef int numConComp = G.order()
    cdef int numConCompPrevIter = numConComp + 1

    # Dictionary to maintain active cheapest edges between pairs of components
    cdef dict components_dict = {}

    while numConComp > 1:
        # Check if number of connected components decreased.
        # Otherwise, the graph is not connected.
        if numConCompPrevIter == numConComp:
            return []
        else:
            numConCompPrevIter = numConComp

        # Iterate over all active edges to identify the cheapest edge between
        # each pair of components (trees of the forest), as well as cheapest
        # active edge incident to a component.
        for e, e_weight in edge_list:
            component1 = partitions.find(e[0])
            component2 = partitions.find(e[1])

            if component1 != component2:
                if component1 in cheapest:
                    if cheapest[component1][1] > e_weight:
                        cheapest[component1] = (e, e_weight)
                else:
                    cheapest[component1] = (e, e_weight)

                if component2 in cheapest:
                    if cheapest[component2][1] > e_weight:
                        cheapest[component2] = (e, e_weight)
                else:
                    cheapest[component2] = (e, e_weight)
                # store the cheapest edge between the two components
                pair = frozenset((component1, component2))
                if pair in components_dict:
                    if components_dict[pair][1] > e_weight:
                        components_dict[pair] = (e, e_weight)
                else:
                    components_dict[pair] = (e, e_weight)

        # Update the list of active edges
        edge_list = components_dict.values()

        # Go through all the current connected components and merge wherever
        # possible
        for v in cheapest:
            e, e_weight = cheapest[v]
            component1 = partitions.find(e[0])
            component2 = partitions.find(e[1])

            if component1 != component2:
                partitions.union(component1, component2)
                T.append(e)
                numConComp = numConComp - 1

        # reset the dictionaries for next iteration
        cheapest = {}
        components_dict = {}

    return T

@cython.binding(True)
def random_spanning_tree(self, output_as_graph=False):
    r"""
    Return a random spanning tree of the graph.

    This uses the Aldous-Broder algorithm ([Broder89]_, [Aldous90]_) to generate
    a random spanning tree with the uniform distribution, as follows.

    Start from any vertex. Perform a random walk by choosing at every step one
    neighbor uniformly at random. Every time a new vertex `j` is met, add the
    edge `(i, j)` to the spanning tree, where `i` is the previous vertex in the
    random walk.

    INPUT:

    - ``output_as_graph`` -- boolean (default: ``False``); whether to return a
      list of edges or a graph

    .. SEEALSO::

        :meth:`~sage.graphs.generic_graph.GenericGraph.spanning_trees_count`
        and :meth:`~sage.graphs.graph.Graph.spanning_trees`

    EXAMPLES::

        sage: G = graphs.TietzeGraph()
        sage: G.random_spanning_tree(output_as_graph=True)
        Graph on 12 vertices
        sage: rg = G.random_spanning_tree(); rg # random
        [(0, 9),
        (9, 11),
        (0, 8),
        (8, 7),
        (7, 6),
        (7, 2),
        (2, 1),
        (1, 5),
        (9, 10),
        (5, 4),
        (2, 3)]
        sage: Graph(rg).is_tree()
        True

    A visual example for the grid graph::

        sage: G = graphs.Grid2dGraph(6, 6)
        sage: pos = G.get_pos()
        sage: T = G.random_spanning_tree(True)
        sage: T.set_pos(pos)
        sage: T.show(vertex_labels=False)

    TESTS::

        sage: G = Graph()
        sage: G.random_spanning_tree()
        Traceback (most recent call last):
        ...
        ValueError: works only for non-empty connected graphs

        sage: G = graphs.CompleteGraph(3).complement()
        sage: G.random_spanning_tree()
        Traceback (most recent call last):
        ...
        ValueError: works only for non-empty connected graphs
    """
    from sage.misc.prandom import randint
    from sage.graphs.graph import Graph

    cdef int N = self.order()

    if not N or not self.is_connected():
        raise ValueError('works only for non-empty connected graphs')

    s = next(self.vertex_iterator())
    cdef set found = set([s])
    cdef int found_nr = 1
    cdef list tree_edges = []
    cdef list neighbors
    while found_nr < N:
        neighbours = self.neighbors(s)
        new_s = neighbours[randint(0, len(neighbours) - 1)]
        if new_s not in found:
            found.add(new_s)
            found_nr += 1
            tree_edges.append((s, new_s))
        s = new_s

    if not output_as_graph:
        return tree_edges
    return Graph(tree_edges)
