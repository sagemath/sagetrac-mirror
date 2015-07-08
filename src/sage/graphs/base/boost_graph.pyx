#*****************************************************************************
#       Copyright (C) 2015 Michele Borassi michele.borassi@imtlucca.it
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
r"""
Interface to run Boost algorithms

Wrapper for a Boost graph. The Boost graphs are Cython C++ variables, and they
cannot be converted to Python objects: as a consequence, only functions defined
with cdef are able to create, read, modify, and delete these graphs.

A very important feature of Boost graph library is that all object are generic:
for instance, adjacency lists can be stored using different data structures,
and (most of) the functions work with all implementations provided. This feature
is implemented in our interface using fused types: however, Cython's support for
fused types is still experimental, and some features are missing. For instance,
there cannot be nested generic function calls, and no variable can have a
generic type, apart from the arguments of a generic function.

All the input functions use pointers, because otherwise we might have problems
with ``delete()``.

**Basic Boost Graph operations:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`clustering_coeff`  | Returns the clustering coefficient of all vertices in the graph.
    :func:`edge_connectivity` | Returns the edge connectivity of the graph.
    :func:`dominator_tree` | Returns a dominator tree of the graph.

Functions
---------
"""

include "sage/ext/interrupt.pxi"

cdef boost_graph_from_sage_graph(BoostGenGraph *g, g_sage):
    r"""
    Initializes the Boost graph ``g`` to be equal to ``g_sage``.

    The Boost graph ``*g`` must represent an empty graph (an exception is raised
    otherwise).
    """

    from sage.graphs.generic_graph import GenericGraph

    if not isinstance(g_sage, GenericGraph):
        raise ValueError("The input parameter must be a Sage graph.")

    if g.num_verts() > 0:
        raise ValueError("The Boost graph in input must be empty")

    N = g_sage.num_verts()
    cdef dict vertex_to_int = {v:i for i,v in enumerate(g_sage.vertices())}

    for i in range(N):
        g.add_vertex()

    for u,v in g_sage.edge_iterator(labels=None):
        g.add_edge(vertex_to_int[u], vertex_to_int[v])

cdef boost_edge_connectivity(BoostVecGenGraph *g):
    r"""
    Computes the edge connectivity of the input Boost graph.

    The output is a pair ``[ec,edges]``, where ``ec`` is the edge connectivity,
    ``edges`` is the list of edges in a minimum cut.
    """
    result = g[0].edge_connectivity()

    cdef int i
    edges = [(result.edges[i], result.edges[i+1])
             for i in range(0, result.edges.size(), 2)]

    return [result.ec, edges]

cpdef edge_connectivity(g):
    r"""
    Computes the edge connectivity of the input graph, using Boost.

    The output is a pair ``[ec,edges]``, where ``ec`` is the edge connectivity,
    ``edges`` is the list of edges in a minimum cut.

    .. SEEALSO::

        :meth:`sage.graphs.generic_graph.GenericGraph.edge_connectivity`

    EXAMPLES:

    Computing the edge connectivity of a clique::

        sage: from sage.graphs.base.boost_graph import edge_connectivity
        sage: g = graphs.CompleteGraph(5)
        sage: edge_connectivity(g)
        [4, [(0, 1), (0, 2), (0, 3), (0, 4)]]

    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph

    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost_und
    cdef BoostVecDiGraph g_boost_dir

    if isinstance(g, Graph):
        boost_graph_from_sage_graph(&g_boost_und, g)
        sig_check()
        return boost_edge_connectivity(&g_boost_und)

    elif isinstance(g, DiGraph):
        from sage.misc.stopgap import stopgap
        stopgap("The edge connectivity of directed graphs is not implemented " +
                "in Boost. The result may be mathematically unreliable.",18753)

        boost_graph_from_sage_graph(&g_boost_dir, g)
        sig_check()
        return boost_edge_connectivity(&g_boost_dir)

    else:
        raise ValueError("The input must be a Sage graph.")

cdef boost_clustering_coeff(BoostGenGraph *g, vertices):
    r"""
    Computes the clustering coefficient of all vertices in the list provided.

    The output is a pair ``[average_clustering_coefficient, clust_of_v]``, where
    ``average_clustering_coefficient`` is the average clustering of the vertices
    in variable ``vertices``, ``clust_of_v`` is a dictionary that associates to
    each vertex (stored as an integer) its clustering coefficient.
    """
    cdef result_cc result
    cdef dict clust_of_v

    if len(vertices) == g.num_verts():
        result = g[0].clustering_coeff_all()
        clust_of_v = {v:result.clust_of_v[v] for v in range(g.num_verts())}
        return [result.average_clustering_coefficient, clust_of_v]

    else:
        clust_of_v = {v:g[0].clustering_coeff(v) for v in vertices}
        return [(sum(clust_of_v.itervalues())/len(clust_of_v)), clust_of_v]


cpdef clustering_coeff(g, vertices = None):
    r"""
    Computes the clustering coefficient of the input graph, using Boost.

    The output is a pair ``[average_clustering_coefficient, clust_of_v]``, where
    ``average_clustering_coefficient`` is the average clustering of the vertices
    in variable ``vertices``, ``clust_of_v`` is a dictionary that associates to
    each vertex its clustering coefficient. If ``vertices`` is ``None``, all
    vertices are considered.

    .. SEEALSO::

        :meth:`sage.graphs.generic_graph.GenericGraph.clustering_coeff`

    INPUT:

    - ``g`` (Graph) - the input graph.

    - ``vertices`` (list) - the list of vertices we need to analyze (if
      ``None``, we will compute the clustering coefficient of all vertices).

    EXAMPLES:

    Computing the clustering coefficient of a clique::

        sage: from sage.graphs.base.boost_graph import clustering_coeff
        sage: g = graphs.CompleteGraph(5)
        sage: clustering_coeff(g)
        [1.0, {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0}]
        sage: clustering_coeff(g, vertices = [0,1,2])
        [1.0, {0: 1.0, 1: 1.0, 2: 1.0}]

    Of a non-clique graph with triangles::

        sage: g = graphs.IcosahedralGraph()
        sage: clustering_coeff(g, vertices=[1,2,3])
        [0.5, {1: 0.5, 2: 0.5, 3: 0.5}]

    With labels::

        sage: g.relabel(list("abcdefghiklm"))
        sage: clustering_coeff(g, vertices="abde")
        [0.5, {'a': 0.5, 'b': 0.5, 'd': 0.5, 'e': 0.5}]
    """
    from sage.graphs.graph import Graph

    sig_on()
    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost
    cdef list g_vertices = g.vertices()
    cdef dict vertex_to_int = {v:i for i,v in enumerate(g_vertices)}

    if not isinstance(g, Graph):
        sig_off()
        raise ValueError("The input must be a Sage graph.")

    boost_graph_from_sage_graph(&g_boost, g)

    if vertices is None:
        vertices = g_vertices

    vertices_boost = [vertex_to_int[v] for v in vertices]
    [average_clustering, clust_v_int] = boost_clustering_coeff(&g_boost, vertices_boost)
    clust_v_sage = {g_vertices[v]: clust_v_int[v] for v in vertices_boost}
    sig_off()
    return [average_clustering, clust_v_sage]


cpdef dominator_tree(g, root, return_dict = False):
    r"""
    Uses Boost to compute the dominator tree of ``g``, rooted at ``root``.

    A node `d` dominates a node `n` if every path from the entry node
    ``root`` to `n` must go through `d`. The immediate dominator of a node
    `n` is the unique node that strictly dominates `n` but does not dominate
    any other node that dominates `n`. A dominator tree is a tree where each
    node's children are those nodes it immediately dominates. For more
    information, see :wikipedia:`Dominator_(graph_theory)`.

    If the graph is connected and undirected, the parent of a vertex `v` is:

     - the root if `v` is in the same biconnected component as the root;

     - the first cut vertex in a path from `v` to the root, otherwise.

    If the graph is not connected, the dominator tree of the whole graph is
    equal to the dominator tree of the connected component of the root.

    If the graph is directed, computing a dominator tree is more complicated,
    and it needs time `O(m\log m)`, where `m` is the number of edges. The
    implementation provided by Boost is the most general one, so it needs time
    `O(m\log m)` even for undirected graphs.

    INPUT:

    - ``g`` (generic_graph) - the input graph.

    - ``root`` (vertex) - the root of the dominator tree.

    - ``return_dict`` (boolean) - if ``True``, the function returns a
      dictionary associating to each vertex its parent in the dominator
      tree. If ``False`` (default), it returns the whole tree, as a ``Graph``
      or a ``DiGraph``.

    OUTPUT:

    The dominator tree, as a graph or as a dictionary, depending on the
    value of ``return_dict``. If the output is a dictionary, it will contain
    ``None`` in correspondence of ``root`` and of vertices that are not
    reachable from ``root``. If the output is a graph, it will not contain
    vertices that are not reachable from ``root``.

    EXAMPLES:

    An undirected grid is biconnected, and its dominator tree is a star
    (everyone's parent is the root)::

        sage: g = graphs.GridGraph([2,2]).dominator_tree((0,0))
        sage: g.to_dictionary()
        {(0, 0): [(0, 1), (1, 0), (1, 1)], (0, 1): [(0, 0)], (1, 0): [(0, 0)], (1, 1): [(0, 0)]}

    If the graph is made by two 3-cycles `C_1,C_2` connected by an edge `(v,w)`,
    with `v \in C_1`, `w \in C_2`, the cut vertices are `v` and `w`, the
    biconnected components are `C_1`, `C_2`, and the edge `(v,w)`. If the root
    is in `C_1`, the parent of each vertex in `C_1` is the root, the parent of
    `w` is `v`, and the parent of each vertex in `C_2` is `w`::

         sage: G = 2 * graphs.CycleGraph(3)
         sage: v = 0
         sage: w = 3
         sage: G.add_edge(v,w)
         sage: G.dominator_tree(1, return_dict=True)
         {0: 1, 1: None, 2: 1, 3: 0, 4: 3, 5: 3}

    An example with a directed graph::

        sage: g = digraphs.Circuit(10).dominator_tree(5)
        sage: g.to_dictionary()
        {0: [1], 1: [2], 2: [3], 3: [4], 4: [], 5: [6], 6: [7], 7: [8], 8: [9], 9: [0]}

    If the output is a dictionary::

        sage: graphs.GridGraph([2,2]).dominator_tree((0,0), return_dict = True)
        {(0, 0): None, (0, 1): (0, 0), (1, 0): (0, 0), (1, 1): (0, 0)}

    TESTS:

    If ``g`` is not a graph, an error is raised::

        sage: from sage.graphs.base.boost_graph import dominator_tree
        sage: dominator_tree('I am not a graph', 0)
        Traceback (most recent call last):
        ...
        ValueError: The input g must be a Sage graph.

    If ``root`` is not a vertex, an error is raised::

        sage: digraphs.TransitiveTournament(10).dominator_tree('Not a vertex!')
        Traceback (most recent call last):
        ...
        ValueError: The input root must be a vertex of g.
        sage: graphs.GridGraph([2,2]).dominator_tree(0)
        Traceback (most recent call last):
        ...
        ValueError: The input root must be a vertex of g.

    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph

    if not isinstance(g, Graph) and not isinstance(g, DiGraph):
        raise ValueError("The input g must be a Sage graph.")
    if not root in g.vertices():
        raise ValueError("The input root must be a vertex of g.")

    sig_on()
    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost_und
    cdef BoostVecDiGraph g_boost_dir
    cdef vector[v_index] result
    cdef dict vertex_to_int = {v:i for i,v in enumerate(g.vertices())}
    cdef dict int_to_vertex = {i:v for i,v in enumerate(g.vertices())}

    if isinstance(g, Graph):
        boost_graph_from_sage_graph(&g_boost_und, g)
        result = <vector[v_index]> g_boost_und.dominator_tree(vertex_to_int[root])

    elif isinstance(g, DiGraph):
        boost_graph_from_sage_graph(&g_boost_dir, g)
        result = <vector[v_index]> g_boost_dir.dominator_tree(vertex_to_int[root])

    sig_off()

    cdef v_index no_parent = -1

    if return_dict:
        return {v:(None if result[vertex_to_int[v]] == no_parent else int_to_vertex[<int> result[vertex_to_int[v]]]) for v in g.vertices()};

    edges = [[int_to_vertex[result[vertex_to_int[v]]], v] for v in g.vertices() if result[vertex_to_int[v]] != no_parent]

    if g.is_directed():
        if len(edges) == 0:
            g = DiGraph()
            g.add_vertex(root)
            return g
        else:
            return DiGraph(edges)
    else:
        if len(edges) == 0:
            g = Graph()
            g.add_vertex(root)
            return g
        else:
            return Graph(edges)
