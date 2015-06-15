#*****************************************************************************
#       Copyright (C) 2015 Michele Borassi michele.borassi@imtlucca.it
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

r"""
Boost Graphs

Wrapper for a Boost graphs. Unfortunately, directed and undirected graphs must
be handled separately.

**Basic Boost Graph operations:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~BoostGraph.add_vertex` | Adds a vertex to the graph.
    :meth:`~BoostGraph.add_edge` | Adds an edge to the graph.
    :meth:`~BoostGraph.add_edge_unsafe` | Adds an edge to the graph without checking that its endpoints exist.
    :meth:`~BoostGraph.num_verts` | Returns the number of vertices.
    :meth:`~BoostGraph.num_edges` | Returns the number of edges.
    :meth:`~BoostGraph.edge_connectivity` | Returns the edge connectivity of the graph.
    :meth:`~BoostDiGraph.add_vertex` | Adds a vertex to the digraph.
    :meth:`~BoostDiGraph.add_edge` | Adds an edge to the digraph.
    :meth:`~BoostDiGraph.add_edge_unsafe` | Adds an edge to the digraph without checking that its endpoints exist.
    :meth:`~BoostDiGraph.num_verts` | Returns the number of vertices.
    :meth:`~BoostDiGraph.num_edges` | Returns the number of edges.
    :meth:`~BoostDiGraph.edge_connectivity` | Returns the edge connectivity of the digraph.


"""
include "sage/ext/interrupt.pxi"

cdef class BoostGraph(object):
    r"""
    Wrapper for an undirected Boost graph.

    The Boost graph is stored in variable ``graph``, which is implemented using
    vectors for storing vertices and neighbors (``vecS``). Vertices are stored
    as numbers from `0` to `n-1`, and a vector of vectors is used to store
    adjacencies. Edge labels are not supported. Multiple edges and self-loops
    are allowed.

    This class lets us create the graph, add vertices, add edges, and print
    number of vertices and number of edges.

    EXAMPLES:

    Creating a `3`-cycle::

        sage: from sage.graphs.base.boost_graph import BoostGraph
        sage: G = BoostGraph()
        sage: G.add_vertex()
        sage: G.add_vertex()
        sage: G.add_vertex()
        sage: G.add_edge(0,1)
        sage: G.add_edge(1,2)
        sage: G.add_edge(2,0)
        sage: print "Vertices:", G.num_verts()
        Vertices: 3
        sage: print "Edges:", G.num_edges()
        Edges: 3

    Multiple edges and self-loops are allowed::

        sage: from sage.graphs.base.boost_graph import BoostGraph
        sage: g = BoostGraph()
        sage: g.add_vertex()
        sage: g.add_vertex()
        sage: g.add_edge(0,1)
        sage: print g.num_edges()
        1
        sage: g.add_edge(0,1)
        sage: print g.num_edges()
        2
        sage: g.add_edge(0,0)
        sage: print g.num_edges()
        3

    """


    def __cinit__(self, G = None):
        """
        Initializes a Boost graph equal to the Sage graph ``G``.

        If ``G`` has edge labels, these labels are ignored. If ``G`` is None,
        the function initializes an empty Boost graph.

        TESTS:

        Generating an empty graph::

            sage: from sage.graphs.base.boost_graph import BoostGraph
            sage: G = BoostGraph()
            sage: print "Vertices:", G.num_verts()
            Vertices: 0
            sage: print "Edges:", G.num_edges()
            Edges: 0

        Converting a random Erdos-Renyi graph::

            sage: from sage.graphs.base.boost_graph import BoostGraph
            sage: G = graphs.RandomGNM(10000,100000)
            sage: H = BoostGraph(G)
            sage: print "Vertices:", H.num_verts()
            Vertices: 10000
            sage: print "Edges:", H.num_edges()
            Edges: 100000

        Given anything else than a graph::

            sage: BoostGraph([])
            Traceback (most recent call last):
            ...
            ValueError: The input parameter must be a Graph.
        """
        cdef int N
        self.graph = new BoostVecGraph()

        if G is None:
            return

        from sage.graphs.graph import Graph

        if not isinstance(G, Graph):
            raise ValueError("The input parameter must be a Graph.")

        N = G.num_verts()
        vertex_to_int = {v:i for i,v in enumerate(G.vertices())}

        for i in range(N):
            self.add_vertex()

        for u,v in G.edge_iterator(labels=None):
            self.add_edge_unsafe(vertex_to_int[u], vertex_to_int[v])

    cpdef add_vertex(self):
        """
        Adds a new vertex to the graph, and stores it in variable vertices to
        access this vertex.

        TESTS:

        Adding one vertex to an empty graph::

            sage: from sage.graphs.base.boost_graph import BoostGraph
            sage: G = BoostGraph()
            sage: G.add_vertex()
            sage: print "Vertices:", G.num_verts()
            Vertices: 1
        """
        self.vertices.push_back(add_vertex(self.graph[0]))

    cpdef add_edge(self, int u, int v):
        """
        Adds a new edge `(u,v)` to the graph, after ensuring that `u` and `v`
        actually exist. If not, an exception is raised.

        TESTS:

        Adding an edge::

            sage: from sage.graphs.base.boost_graph import BoostGraph
            sage: G = BoostGraph()
            sage: G.add_vertex()
            sage: G.add_vertex()
            sage: G.add_edge(0, 1)
            sage: print "Edges:", G.num_edges()
            Edges: 1

        Adding an edge between two non-existing vertices::

            sage: from sage.graphs.base.boost_graph import BoostGraph
            sage: G = BoostGraph()
            sage: G.add_edge(0, 1)
            Traceback (most recent call last):
            ...
            ValueError: Vertex 0 does not exist.
        """
        if num_vertices(self.graph[0]) <= u:
            raise ValueError("Vertex " + str(u) + " does not exist.")
        if num_vertices(self.graph[0]) <= v:
            raise ValueError("Vertex " + str(v) + " does not exist.")
        self.add_edge_unsafe(u, v)

    cpdef add_edge_unsafe(self, int u, int v):
        """
        Adds a new edge `(u,v)` to the graph, without checking the existence of
        the vertices.

        TESTS:

        Adding an edge::

            sage: from sage.graphs.base.boost_graph import BoostGraph
            sage: G = BoostGraph()
            sage: G.add_vertex()
            sage: G.add_vertex()
            sage: G.add_edge_unsafe(0, 1)
            sage: print "Edges:", G.num_edges()
            Edges: 1
        """
        add_edge(self.vertices[u], self.vertices[v], self.graph[0])

    cpdef num_verts(self):
        """
        Returns the number of vertices in the graph.

        TESTS:

        The number of vertices of an empty graph::

            sage: from sage.graphs.base.boost_graph import BoostGraph
            sage: print "Vertices:", BoostGraph().num_verts()
            Vertices: 0
        """

        return num_vertices(self.graph[0])

    cpdef num_edges(self):
        """
        Returns the number of edges in the graph.

        TESTS:

        The number of edges of an empty graph::

            sage: from sage.graphs.base.boost_graph import BoostGraph
            sage: print "Edges:", BoostGraph().num_edges()
            Edges: 0
        """
        return num_edges(self.graph[0])


    cpdef edge_connectivity(self):
        """
        Computes the edge connectivity of this graph.

        See :meth:`generic_graph <sage.graphs.generic_graph.GenericGraph.edge_connectivity>`
        module for more information.

        The output is a pair ``[ec, edges]``, where ``ec`` is the edge
        connectivity and ``edges`` is the list of edges in a minimum cut.

        EXAMPLE::

            sage: from sage.graphs.base.boost_graph import BoostGraph
            sage: G = BoostGraph(graphs.PathGraph(2))
            sage: G.edge_connectivity()
            [1, [[0, 1]]]

        """
        sig_on()
        cdef edge_container_t_und * disconnecting_set = new edge_container_t_und()
        cdef back_insert_iterator[edge_container_t_und] * inserter = new back_insert_iterator[edge_container_t_und](disconnecting_set[0])
        cdef int ec = <int> edge_connectivity(self.graph[0], inserter[0])
        cdef int i
        cdef edge_descriptor_und e
        edges = []

        for i in range(ec):
            e = disconnecting_set.at(i)
            s = <int> source(e, self.graph[0])
            t = <int> target(e, self.graph[0])
            edges.append([s, t])
        sig_off()
        return [ec, edges]


cdef class BoostDiGraph(object):
    r"""
    Wrapper for a directed Boost graph.

    The Boost graph is stored in variable ``graph``, which is implemented using
    vectors for storing vertices and neighbors (``vecS``). Vertices are stored
    as numbers from `0` to `n-1`, and a vector of vectors is used to store
    adjacencies. Edge labels are not supported. Multiple edges and self-loops
    are allowed.

    This class lets us create the graph, add vertices, add edges, and print
    number of vertices and number of edges.

    EXAMPLES:

    Creating a directed `3`-cycle::

        sage: from sage.graphs.base.boost_graph import BoostGraph
        sage: G = BoostGraph()
        sage: G.add_vertex()
        sage: G.add_vertex()
        sage: G.add_vertex()
        sage: G.add_edge(0,1)
        sage: G.add_edge(1,2)
        sage: G.add_edge(2,0)
        sage: print "Vertices:", G.num_verts()
        Vertices: 3
        sage: print "Edges:", G.num_edges()
        Edges: 3

    Multiple edges and self-loops are allowed::

        sage: from sage.graphs.base.boost_graph import BoostDiGraph
        sage: g = BoostDiGraph()
        sage: g.add_vertex()
        sage: g.add_vertex()
        sage: g.add_edge(0,1)
        sage: print g.num_edges()
        1
        sage: g.add_edge(0,1)
        sage: print g.num_edges()
        2
        sage: g.add_edge(0,0)
        sage: print g.num_edges()
        3

    """


    def __cinit__(self, G = None):
        """
        Initializes a Boost digraph equal to the Sage digraph ``G``.

        If ``G`` has edge labels, these labels are ignored. If ``G`` is None,
        the function initializes an empty Boost digraph.

        TESTS:

        Generating an empty graph::

            sage: from sage.graphs.base.boost_graph import BoostDiGraph
            sage: G = BoostDiGraph()
            sage: print "Vertices:", G.num_verts()
            Vertices: 0
            sage: print "Edges:", G.num_edges()
            Edges: 0

        Converting a random circuit graph::

            sage: from sage.graphs.base.boost_graph import BoostDiGraph
            sage: G = digraphs.Circuit(10)
            sage: H = BoostDiGraph(G)
            sage: print "Vertices:", H.num_verts()
            Vertices: 10
            sage: print "Edges:", H.num_edges()
            Edges: 10

        Given anything else than a graph::

            sage: BoostDiGraph([])
            Traceback (most recent call last):
            ...
            ValueError: The input parameter must be a DiGraph.
        """
        cdef int N
        self.graph = new BoostVecDiGraph()

        if G is None:
            return

        from sage.graphs.digraph import DiGraph

        if not isinstance(G, DiGraph):
            raise ValueError("The input parameter must be a DiGraph.")

        N = G.num_verts()
        vertex_to_int = {v:i for i,v in enumerate(G.vertices())}

        for i in range(N):
            self.add_vertex()

        for u,v in G.edge_iterator(labels=None):
            self.add_edge_unsafe(vertex_to_int[u], vertex_to_int[v])

    cpdef add_vertex(self):
        """
        Adds a new vertex to the graph, and stores it in variable vertices
        to access this vertex.

        TESTS:

        Adding one vertex to an empty graph::

            sage: from sage.graphs.base.boost_graph import BoostDiGraph
            sage: G = BoostDiGraph()
            sage: G.add_vertex()
            sage: print "Vertices:", G.num_verts()
            Vertices: 1
        """
        self.vertices.push_back(add_vertex(self.graph[0]))

    cpdef add_edge(self, int u, int v):
        r"""
        Adds a new edge `(u,v)` to the graph, after ensuring that `u` and
        `v` actually exist. If not, an exception is raised.

        TESTS:

        Adding an edge::

            sage: from sage.graphs.base.boost_graph import BoostDiGraph
            sage: G = BoostDiGraph()
            sage: G.add_vertex()
            sage: G.add_vertex()
            sage: G.add_edge(0, 1)
            sage: print "Edges:", G.num_edges()
            Edges: 1

        Adding an edge between two non-existing vertices::

            sage: from sage.graphs.base.boost_graph import BoostDiGraph
            sage: G = BoostDiGraph()
            sage: G.add_edge(0, 1)
            Traceback (most recent call last):
            ...
            ValueError: Vertex 0 does not exist.
        """
        if num_vertices(self.graph[0]) <= u:
            raise ValueError("Vertex " + str(u) + " does not exist.")
        if num_vertices(self.graph[0]) <= v:
            raise ValueError("Vertex " + str(v) + " does not exist.")
        self.add_edge_unsafe(u, v)

    cpdef add_edge_unsafe(self, int u, int v):
        r"""
        Adds a new edge `(u,v)` to the graph, without checking the existence of
        the vertices.

        TESTS:

        Adding two edges::

            sage: from sage.graphs.base.boost_graph import BoostDiGraph
            sage: G = BoostDiGraph()
            sage: G.add_vertex()
            sage: G.add_vertex()
            sage: G.add_edge_unsafe(0, 1)
            sage: G.add_edge_unsafe(1, 0)
            sage: print "Edges:", G.num_edges()
            Edges: 2
        """
        add_edge(self.vertices[u], self.vertices[v], self.graph[0])

    cpdef num_verts(self):
        """
        Returns the number of vertices in the graph.

        TESTS:

        The number of vertices of an empty graph::

            sage: from sage.graphs.base.boost_graph import BoostDiGraph
            sage: print "Vertices:", BoostDiGraph().num_verts()
            Vertices: 0
        """

        return num_vertices(self.graph[0])

    cpdef num_edges(self):
        """
        Returns the number of edges in the graph.

        TESTS:

        The number of edges of an empty graph::

            sage: from sage.graphs.base.boost_graph import BoostDiGraph
            sage: print "Edges:", BoostDiGraph().num_edges()
            Edges: 0
        """
        return num_edges(self.graph[0])

    cpdef edge_connectivity(self):
        """
        Computes the edge connectivity of this graph.

        See :meth:`generic_graph <sage.graphs.generic_graph.GenericGraph.edge_connectivity>`
        module for more information.

        At the moment, the Boost algorithm may contain bugs, and this is
        explicitly outputted when this function is called.

        The output is a pair ``[ec, edges]``, where ``ec`` is the edge
        connectivity and ``edges`` is the list of edges in a minimum cut.

        EXAMPLE::

            sage: from sage.graphs.base.boost_graph import BoostDiGraph
            sage: g = BoostDiGraph(digraphs.Path(3))
            sage: g.edge_connectivity()
            The directed edge connectivity algorithm implemented in the Boost graph library is not reliable. The result could be wrong.
            [1, [[0, 1]]]

        """
        print("The directed edge connectivity algorithm implemented in the " +
              "Boost graph library is not reliable. The result could be wrong.")
        sig_on()
        cdef edge_container_t_dir * disconnecting_set = new edge_container_t_dir()
        cdef back_insert_iterator[edge_container_t_dir] * inserter = new back_insert_iterator[edge_container_t_dir](disconnecting_set[0])
        cdef int ec = <int> edge_connectivity(self.graph[0], inserter[0])
        cdef int i
        cdef edge_descriptor_dir e
        edges = []

        for i in range(ec):
            e = disconnecting_set.at(i)
            s = <int> source(e, self.graph[0])
            t = <int> target(e, self.graph[0])
            edges.append([s, t])
        sig_off()
        return [ec, edges]