r"""
Spanning trees

This module is a collection of algorithms on spanning trees. Also included in
the collection are algorithms for minimum spanning trees. See the book
[JoynerNguyenCohen2010]_ for descriptions of spanning tree algorithms,
including minimum spanning trees.

.. SEEALSO::

   * :meth:`GenericGraph.min_spanning_tree
     <sage.graphs.generic_graph.GenericGraph.min_spanning_tree>`.

**Todo**

* Rewrite :func:`kruskal` to use priority queues. Once Cython has support
  for generators and the ``yield`` statement, rewrite :func:`kruskal` to use
  ``yield``.
* Parallel version of Boruvka's algorithm.
* Randomized spanning tree construction.

REFERENCES:

.. [Aldous90] \D. Aldous, 'The random walk construction of
  uniform spanning trees', SIAM J Discrete Math 3 (1990),
  450-465.

.. [Broder89] \A. Broder, 'Generating random spanning trees',
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

#*****************************************************************************
#       Copyright (c) 2007 Jason Grout <jason-sage@creativetrax.com>
#       Copyright (c) 2009 Mike Hansen <mhansen@gmail.com>
#       Copyright (c) 2010 Gregory McWhirter <gmcwhirt@uci.edu>
#       Copyright (c) 2010 Minh Van Nguyen <nguyenminh2@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

cimport cython


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

    - ``weight_function`` (function) - a function that inputs an edge ``e``
      and outputs its weight. An edge has the form ``(u,v,l)``, where ``u``
      and ``v`` are vertices, ``l`` is a label (that can be of any kind).
      The ``weight_function`` can be used to transform the label into a
      weight. In particular:

      - if ``weight_function`` is not ``None``, the weight of an edge ``e``
        is ``weight_function(e)``;

      - if ``weight_function`` is ``None`` (default) and ``g`` is weighted
        (that is, ``g.weighted()==True``), the weight of an edge
        ``e=(u,v,l)`` is ``l``, independently on which kind of object ``l``
        is: the ordering of labels relies on Python's operator ``<``;

      - if ``weight_function`` is ``None`` and ``g`` is not weighted, we set
        all weights to 1 (hence, the output can be any spanning tree).

    - ``check`` -- Whether to first perform sanity checks on the input
      graph ``G``. Default: ``check=False``. If we toggle ``check=True``, the
      following sanity checks are first performed on ``G`` prior to running
      Kruskal's algorithm on that input graph:

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

    EXAMPLES:

    An example from pages 727--728 in [Sahni2000]_. ::

        sage: from sage.graphs.spanning_tree import kruskal
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: E = kruskal(G, check=True); E
        [(1, 6, 10), (2, 3, 16), (2, 7, 14), (3, 4, 12), (4, 5, 22), (5, 6, 25)]

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
        ....:     E = {}
        ....:     for u, v, _ in G.multiple_edges():
        ....:         E.setdefault(u, v)
        ....:     for u in E:
        ....:         W = sorted(G.edge_label(u, E[u]))
        ....:         for w in W[1:]:
        ....:             G.delete_edge(u, E[u], w)
        ....:     G.allow_multiple_edges(False)
        sage: sanitize(H)
        sage: H
        Graph on 7 vertices
        sage: kruskal(G, check=True) == kruskal(H, check=True)
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
        [('BOS', 'JFK', 187), ('BWI', 'JFK', 184), ('BWI', 'MIA', 946), ('BWI', 'ORD', 621), ('DFW', 'LAX', 1235), ('DFW', 'ORD', 802), ('JFK', 'PVD', 144), ('LAX', 'SFO', 337)]

    An example from pages 568--569 in [CormenEtAl2001]_. ::

        sage: G = Graph({"a":{"b":4, "h":8}, "b":{"c":8, "h":11},
        ....: "c":{"d":7, "f":4, "i":2}, "d":{"e":9, "f":14},
        ....: "e":{"f":10}, "f":{"g":2}, "g":{"h":1, "i":6}, "h":{"i":7}})
        sage: G.weighted(True)
        sage: kruskal(G, check=True)
        [('a', 'b', 4), ('a', 'h', 8), ('c', 'd', 7), ('c', 'f', 4), ('c', 'i', 2), ('d', 'e', 9), ('f', 'g', 2), ('g', 'h', 1)]

    An example with custom edge labels::

        sage: G = Graph([[0,1,1],[1,2,1],[2,0,10]], weighted=True)
        sage: weight = lambda e:3-e[0]-e[1]
        sage: kruskal(G, check=True)
        [(0, 1, 1), (1, 2, 1)]
        sage: kruskal(G, wfunction=weight, check=True)
        [(0, 2, 10), (1, 2, 1)]
        sage: kruskal(G, wfunction=weight, check=False)
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

    If the input graph is a tree, then return its edges. ::

        sage: T = graphs.RandomTree(randint(1, 50))  # long time
        sage: T.edges() == kruskal(T, check=True)  # long time
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
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("The input G must be an undirected graph.")

    sortedE_iter = None
    # sanity checks
    if check:
        if G.order() == 0:
            return []
        if not G.is_connected():
            return []
        # G is now assumed to be a nonempty connected graph
        if G.num_verts() == G.num_edges() + 1:
            # G is a tree
            return G.edges()
        g = G.to_simple(to_undirected=False, keep_label='min')
    else:
        g = G

    # G is assumed to be connected, undirected, and with at least a vertex
    # We sort edges, as specified.
    if wfunction is None:
        if g.weighted():
            from operator import itemgetter
            sortedE_iter = iter(sorted(g.edges(), key=itemgetter(2)))
        else:
            sortedE_iter = iter(sorted(g.edges()))
    else:
        sortedE_iter = iter(sorted(g.edges(), key=wfunction))


    # Kruskal's algorithm
    T = []
    cdef int n = g.order()
    cdef int m = n - 1
    cdef int i = 0  # count the number of edges added so far
    union_find = dict()
    while i < m:
        e = next(sortedE_iter)
        components = []
        # acyclic test via union-find
        for startv in iter(e[0:2]):
            v = startv
            children = []
            # find the component a vertex lives in
            while v in union_find:
                children.append(v)
                v = union_find[v]
            # compress the paths as much as we can for efficiency reasons
            for c in children:
                union_find[c] = v
            components.append(v)
        if components[0] != components[1]:
            i += 1
            # NOTE: Once Cython supports generator and the yield statement,
            # we should replace the following line with a yield statement.
            # That way, we could access the edge of a minimum spanning tree
            # immediately after it is found, instead of waiting for all the
            # edges to be found and return the edges as a list.
            T.append(e)
            # union the components by making one the parent of the other
            union_find[components[0]] = components[1]
    return sorted(T)



from sage.sets.disjoint_set import *

cpdef boruvka(G, wfunction=None, bint check=False):
    r"""
    Minimum spanning tree using Boruvka's algorithm.

    This function assumes that we can only compute minimum spanning trees for
    undirected graphs. Such graphs can be weighted or unweighted, and they can
    have multiple edges (since we are computing the minimum spanning tree, only
    the minimum weight among all `(u,v)`-edges is considered, for each pair
    of vertices `u`, `v`).

    INPUT:

    - ``G`` -- an undirected graph.

    - ``weight_function`` (function) - a function that inputs an edge ``e``
      and outputs its weight. An edge has the form ``(u,v,l)``, where ``u``
      and ``v`` are vertices, ``l`` is a label (that can be of any kind).
      The ``weight_function`` can be used to transform the label into a
      weight. In particular:

      - if ``weight_function`` is not ``None``, the weight of an edge ``e``
        is ``weight_function(e)``;

      - if ``weight_function`` is ``None`` (default) and ``g`` is weighted
        (that is, ``g.weighted()==True``), the weight of an edge
        ``e=(u,v,l)`` is ``l``, independently on which kind of object ``l``
        is: the ordering of labels relies on Python's operator ``<``;

      - if ``weight_function`` is ``None`` and ``g`` is not weighted, we set
        all weights to 1 (hence, the output can be any spanning tree).

    - ``check`` -- Whether to first perform sanity checks on the input
      graph ``G``. Default: ``check=False``. If we toggle ``check=True``, the
      following sanity checks are first performed on ``G`` prior to running
      Boruvka's algorithm on that input graph:

      - Is ``G`` the null graph or graph on one vertex?
      - Is ``G`` disconnected?
      - Is ``G`` a tree?

      By default, we turn off the sanity checks for performance reasons. This
      means that by default the function assumes that its input graph is
      connected, and has at least one vertex. Otherwise, you should set
      ``check=True`` to perform some sanity checks and preprocessing on the
      input graph. 

    OUTPUT:

    The edges of a minimum spanning tree of ``G``, if one exists, otherwise
    returns the empty list.

    .. SEEALSO::

        - :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`

    EXAMPLES:

    An example from pages 727--728 in [Sahni2000]_. ::

        sage: from sage.graphs.spanning_tree import boruvka
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: E = boruvka(G, check=True); E
        [(1, 6, 10), (2, 7, 14), (3, 4, 12), (4, 5, 22), (5, 6, 25), (2, 3, 16)]

    Variants of the previous example. ::

        sage: H = Graph(G.edges(labels=False))
        sage: boruvka(H, check=True)
        [(1, 2, None), (2, 3, None), (3, 4, None), (4, 5, None), (1, 6, None), (2, 7, None)]
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
        ....:     E = {}
        ....:     for u, v, _ in G.multiple_edges():
        ....:         E.setdefault(u, v)
        ....:     for u in E:
        ....:         W = sorted(G.edge_label(u, E[u]))
        ....:         for w in W[1:]:
        ....:             G.delete_edge(u, E[u], w)
        ....:     G.allow_multiple_edges(False)
        sage: sanitize(H)
        sage: H
        Graph on 7 vertices
        sage: boruvka(G, check=True) == boruvka(H, check=True)
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
        sage: boruvka(G, check=True)
        [('JFK', 'PVD', 144),
        ('BWI', 'MIA', 946),
        ('DFW', 'ORD', 802),
        ('BOS', 'JFK', 187),
        ('LAX', 'SFO', 337),
        ('BWI', 'JFK', 184),
        ('BWI', 'ORD', 621),
        ('DFW', 'LAX', 1235)]

    An example from pages 568--569 in [CormenEtAl2001]_. ::

        sage: G = Graph({"a":{"b":4, "h":8}, "b":{"c":8, "h":11},
        ....: "c":{"d":7, "f":4, "i":2}, "d":{"e":9, "f":14},
        ....: "e":{"f":10}, "f":{"g":2}, "g":{"h":1, "i":6}, "h":{"i":7}})
        sage: G.weighted(True)
        sage: boruvka(G, check=True)
        [('a', 'b', 4), ('c', 'i', 2), ('d', 'e', 9), ('c', 'd', 7), ('g', 'h', 1), ('f', 'g', 2), ('b', 'c', 8), ('c', 'f', 4)]

    An example with custom edge labels::

        sage: G = Graph([[0,1,1],[1,2,1],[2,0,10]], weighted=True)
        sage: weight = lambda e:3-e[0]-e[1]
        sage: boruvka(G, check=True)
        [(0, 1, 1), (1, 2, 1)]
        sage: boruvka(G, wfunction=weight, check=True)
        [(0, 2, 10), (1, 2, 1)]
        sage: boruvka(G, wfunction=weight, check=False)
        [(0, 2, 10), (1, 2, 1)]

    TESTS:

    The input graph must not be empty. ::

        sage: from sage.graphs.spanning_tree import boruvka
        sage: boruvka(graphs.EmptyGraph(), check=True)
        []
        sage: boruvka(Graph(), check=True)
        []
        sage: boruvka(Graph(multiedges=True), check=True)
        []
        sage: boruvka(Graph(loops=True), check=True)
        []
        sage: boruvka(Graph(multiedges=True, loops=True), check=True)
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
        sage: boruvka(G, check=True)  # long time
        []
        sage: G = my_disconnected_graph(100, 50, directed=False, multiedges=True, loops=False)  # long time
        sage: boruvka(G, check=True)  # long time
        []
        sage: G = my_disconnected_graph(100, 50, directed=False, multiedges=True, loops=True)  # long time
        sage: boruvka(G, check=True)  # long time
        []

    If the input graph is a tree, then return its edges. ::

        sage: T = graphs.RandomTree(randint(1, 50))  # long time
        sage: T.edges() == boruvka(T, check=True)  # long time
        True

    If the input is not a Graph::

        sage: boruvka("I am not a graph")
        Traceback (most recent call last):
        ...
        ValueError: The input G must be an undirected graph.
        sage: boruvka(digraphs.Path(10))
        Traceback (most recent call last):
        ...
        ValueError: The input G must be an undirected graph.
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("The input G must be an undirected graph.")

    if G.order() <= 1:
            return []

    # sanity checks
    if check:
        if not G.is_connected():
            return []
        # G is now assumed to be a nonempty connected graph
        if G.num_verts() == G.num_edges() + 1:
            # G is a tree
            return G.edges()

    # Boruvka's algorithm

    # initially, each vertex is a connected component
    partitions = DisjointSet(G.vertices())
    cheapest = {} # a dictionary to store the least weight outgoing edge for each component
    T = [] # stores the edges in minimum spanning tree
    numConnectedComponents = len(G.vertices()) 
    numConnectedComponentsPrevIter = len(G.vertices()) + 1

    # Store an edge list to keep track of active edges
    edge_list = G.edges()
    components_dict = {} # to store the pairwise cheapest edges

    while numConnectedComponents > 1:
        # Check if number of connected components decreased.
        # Otherwise, the graph is not connected.
        if (numConnectedComponentsPrevIter == numConnectedComponents):
            return []
        else:
            numConnectedComponentsPrevIter = numConnectedComponents

        # If the two endpoints of current edge belong to
        # same component, ignore the edge. 
        # Else check if current edge has lesser weight than previous
        # cheapest edges of component1 and component2.
        # Before that, check if component1 and 2 are present in 'cheapest' dict.
        # Also, store a dictionary to maintain active cheapest edges between pairs of components
        for e in edge_list:
            component1 = partitions.find(e[0])
            component2 = partitions.find(e[1])

            if component1 != component2 : 
                pair = (component1, component2) if (component1 < component2) else (component2, component1)
                if wfunction is None:
                    if component1 in cheapest:
                        if cheapest[component1][2] > e[2]:
                            cheapest[component1] = e
                    else:
                        cheapest[component1] = e

                    if component2 in cheapest:
                        if cheapest[component2][2] > e[2]:
                            cheapest[component2] = e
                    else:
                        cheapest[component2] = e
                    # store the cheapest edge between the two components
                    if pair in components_dict:
                        if components_dict[pair][2] > e[2]:
                            components_dict[pair] = e
                    else:
                        components_dict[pair] = e
                    
                else: # if a weight function is provided
                    e_weight = wfunction(e)
                    if component1 in cheapest:
                        if wfunction(cheapest[component1]) > e_weight:
                            cheapest[component1] = e
                    else:
                        cheapest[component1] = e

                    if component2 in cheapest:
                        if wfunction(cheapest[component2]) > e_weight:
                            cheapest[component2] = e
                    else:
                        cheapest[component2] = e
                    # store the cheapest edge between the two components
                    if pair in components_dict:
                        if wfunction(components_dict[pair]) > e_weight:
                            components_dict[pair] = e
                    else:
                        components_dict[pair] = e
        
        edge_list = components_dict.values() # active edges

        # Go through all the current connected components
        # and merge wherever possible
        for v in cheapest:
            e = cheapest[v]
            component1 = partitions.find(e[0])
            component2 = partitions.find(e[1])

            if component1 != component2 :
                partitions.union(component1, component2)
                T.append(e)
                numConnectedComponents = numConnectedComponents - 1
         
        # reset the dictionaries for next iteration
        cheapest = {}
        components_dict = {}

    return T

@cython.binding(True)
def random_spanning_tree(self, output_as_graph=False):
    r"""
    Return a random spanning tree of the graph.

    This uses the Aldous-Broder algorithm ([Broder89]_, [Aldous90]_)
    to generate a random spanning tree with the uniform distribution,
    as follows.

    Start from any vertex. Perform a random walk by choosing at every
    step one neighbor uniformly at random. Every time a new vertex `j`
    is met, add the edge `(i, j)` to the spanning tree, where `i` is
    the previous vertex in the random walk.

    INPUT:

    - ``output_as_graph`` -- boolean (default: ``False``) whether to return a
      list of edges or a graph.

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

    if N == 0 or not self.is_connected():
        raise ValueError('works only for non-empty connected graphs')

    s = next(self.vertex_iterator())
    found = set([s])
    cdef int found_nr = 1
    tree_edges = []
    while found_nr < N:
        neighbours = self.neighbors(s)
        new_s = neighbours[randint(0, len(neighbours) - 1)]
        if not(new_s in found):
            found.add(new_s)
            found_nr += 1
            tree_edges += [(s, new_s)]
        s = new_s

    if not output_as_graph:
        return tree_edges
    return Graph(tree_edges)
