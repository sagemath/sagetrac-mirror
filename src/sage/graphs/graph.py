# -*- coding: utf-8 -*-
r"""
Undirected graphs

This module implements functions and operations involving undirected graphs.

{INDEX_OF_METHODS}

AUTHORS:

- Robert L. Miller (2006-10-22): initial version

- William Stein (2006-12-05): Editing

- Robert L. Miller (2007-01-13): refactoring, adjusting for NetworkX-0.33, fixed
   plotting bugs (2007-01-23): basic tutorial, edge labels, loops, multiple
   edges and arcs (2007-02-07): graph6 and sparse6 formats, matrix input

- Emily Kirkmann (2007-02-11): added graph_border option to plot and show

- Robert L. Miller (2007-02-12): vertex color-maps, graph boundaries, graph6
   helper functions in Cython

- Robert L. Miller Sage Days 3 (2007-02-17-21): 3d plotting in Tachyon

- Robert L. Miller (2007-02-25): display a partition

- Robert L. Miller (2007-02-28): associate arbitrary objects to vertices, edge
   and arc label display (in 2d), edge coloring

- Robert L. Miller (2007-03-21): Automorphism group, isomorphism check,
   canonical label

- Robert L. Miller (2007-06-07-09): NetworkX function wrapping

- Michael W. Hansen (2007-06-09): Topological sort generation

- Emily Kirkman, Robert L. Miller Sage Days 4: Finished wrapping NetworkX

- Emily Kirkman (2007-07-21): Genus (including circular planar, all embeddings
   and all planar embeddings), all paths, interior paths

- Bobby Moretti (2007-08-12): fixed up plotting of graphs with edge colors
   differentiated by label

- Jason Grout (2007-09-25): Added functions, bug fixes, and general enhancements

- Robert L. Miller (Sage Days 7): Edge labeled graph isomorphism

- Tom Boothby (Sage Days 7): Miscellaneous awesomeness

- Tom Boothby (2008-01-09): Added graphviz output

- David Joyner (2009-2): Fixed docstring bug related to GAP.

- Stephen Hartke (2009-07-26): Fixed bug in blocks_and_cut_vertices() that
   caused an incorrect result when the vertex 0 was a cut vertex.

- Stephen Hartke (2009-08-22): Fixed bug in blocks_and_cut_vertices() where the
   list of cut_vertices is not treated as a set.

- Anders Jonsson (2009-10-10): Counting of spanning trees and out-trees added.

- Nathann Cohen (2009-09) : Cliquer, Connectivity, Flows and everything that
                             uses Linear Programming and class numerical.MIP

- Nicolas M. Thiery (2010-02): graph layout code refactoring, dot2tex/graphviz
  interface

- David Coudert (2012-04) : Reduction rules in vertex_cover.

- Birk Eisermann (2012-06): added recognition of weakly chordal graphs and
                            long-hole-free / long-antihole-free graphs

- Alexandre P. Zuge (2013-07): added join operation.

- Amritanshu Prasad (2014-08): added clique polynomial

- Julian Rüth (2018-06-21): upgrade to NetworkX 2

- David Coudert (2018-10-07): cleaning

- Amanda Francis, Caitlin Lienkaemper, Kate Collins, Rajat Mittal (2019-03-10):
  methods for computing effective resistance

- Amanda Francis, Caitlin Lienkaemper, Kate Collins, Rajat Mittal (2019-03-19):
  most_common_neighbors and common_neighbors_matrix added.

- Jean-Florent Raymond (2019-04): is_redundant, is_dominating,
   private_neighbors

Graph Format
------------

Supported formats
~~~~~~~~~~~~~~~~~

Sage Graphs can be created from a wide range of inputs. A few examples are
covered here.

- NetworkX dictionary format:

   ::

       sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], \
             5: [7, 8], 6: [8,9], 7: [9]}
       sage: G = Graph(d); G
       Graph on 10 vertices
       sage: G.plot().show()    # or G.show()

- A NetworkX graph:

   ::

       sage: import networkx
       sage: K = networkx.complete_bipartite_graph(12,7)
       sage: G = Graph(K)
       sage: G.degree()
       [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 12, 12, 12, 12, 12, 12, 12]

- graph6 or sparse6 format:

   ::

       sage: s = ':I`AKGsaOs`cI]Gb~'
       sage: G = Graph(s, sparse=True); G
       Looped multi-graph on 10 vertices
       sage: G.plot().show()    # or G.show()

   Note that the ``\`` character is an escape character in Python, and also a
   character used by graph6 strings:

   ::

       sage: G = Graph('Ihe\n@GUA')
       Traceback (most recent call last):
       ...
       RuntimeError: the string (Ihe) seems corrupt: for n = 10, the string is too short

   In Python, the escaped character ``\`` is represented by ``\\``:

   ::

       sage: G = Graph('Ihe\\n@GUA')
       sage: G.plot().show()    # or G.show()

- adjacency matrix: In an adjacency matrix, each column and each row represent a
   vertex. If a 1 shows up in row `i`, column `j`, there is an edge `(i,j)`.

   ::

       sage: M = Matrix([(0,1,0,0,1,1,0,0,0,0),(1,0,1,0,0,0,1,0,0,0), \
       (0,1,0,1,0,0,0,1,0,0), (0,0,1,0,1,0,0,0,1,0),(1,0,0,1,0,0,0,0,0,1), \
       (1,0,0,0,0,0,0,1,1,0), (0,1,0,0,0,0,0,0,1,1),(0,0,1,0,0,1,0,0,0,1), \
       (0,0,0,1,0,1,1,0,0,0), (0,0,0,0,1,0,1,1,0,0)])
       sage: M
       [0 1 0 0 1 1 0 0 0 0]
       [1 0 1 0 0 0 1 0 0 0]
       [0 1 0 1 0 0 0 1 0 0]
       [0 0 1 0 1 0 0 0 1 0]
       [1 0 0 1 0 0 0 0 0 1]
       [1 0 0 0 0 0 0 1 1 0]
       [0 1 0 0 0 0 0 0 1 1]
       [0 0 1 0 0 1 0 0 0 1]
       [0 0 0 1 0 1 1 0 0 0]
       [0 0 0 0 1 0 1 1 0 0]
       sage: G = Graph(M); G
       Graph on 10 vertices
       sage: G.plot().show()    # or G.show()

- incidence matrix: In an incidence matrix, each row represents a vertex and
   each column represents an edge.

   ::

       sage: M = Matrix([(-1, 0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0),
       ....:             ( 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0),
       ....:             ( 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0),
       ....:             ( 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0),
       ....:             ( 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1),
       ....:             ( 0, 0, 0, 0, 0,-1, 0, 0, 0, 1, 1, 0, 0, 0, 0),
       ....:             ( 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 1, 0, 0, 0),
       ....:             ( 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1, 0, 0),
       ....:             ( 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 1, 0),
       ....:             ( 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 1)])
       sage: M
       [-1  0  0  0  1  0  0  0  0  0 -1  0  0  0  0]
       [ 1 -1  0  0  0  0  0  0  0  0  0 -1  0  0  0]
       [ 0  1 -1  0  0  0  0  0  0  0  0  0 -1  0  0]
       [ 0  0  1 -1  0  0  0  0  0  0  0  0  0 -1  0]
       [ 0  0  0  1 -1  0  0  0  0  0  0  0  0  0 -1]
       [ 0  0  0  0  0 -1  0  0  0  1  1  0  0  0  0]
       [ 0  0  0  0  0  0  0  1 -1  0  0  1  0  0  0]
       [ 0  0  0  0  0  1 -1  0  0  0  0  0  1  0  0]
       [ 0  0  0  0  0  0  0  0  1 -1  0  0  0  1  0]
       [ 0  0  0  0  0  0  1 -1  0  0  0  0  0  0  1]
       sage: G = Graph(M); G
       Graph on 10 vertices
       sage: G.plot().show()    # or G.show()
       sage: DiGraph(matrix(2,[0,0,-1,1]), format="incidence_matrix")
       Traceback (most recent call last):
       ...
       ValueError: there must be two nonzero entries (-1 & 1) per column

- a list of edges::

       sage: g = Graph([(1,3),(3,8),(5,2)])
       sage: g
       Graph on 5 vertices

- an igraph Graph::

       sage: import igraph                                # optional - python_igraph
       sage: g = Graph(igraph.Graph([(1,3),(3,2),(0,2)])) # optional - python_igraph
       sage: g                                            # optional - python_igraph
       Graph on 4 vertices

Generators
----------

Use ``graphs(n)`` to iterate through all non-isomorphic graphs of given size::

    sage: for g in graphs(4):
    ....:     print(g.degree_sequence())
    [0, 0, 0, 0]
    [1, 1, 0, 0]
    [2, 1, 1, 0]
    [3, 1, 1, 1]
    [1, 1, 1, 1]
    [2, 2, 1, 1]
    [2, 2, 2, 0]
    [3, 2, 2, 1]
    [2, 2, 2, 2]
    [3, 3, 2, 2]
    [3, 3, 3, 3]

Similarly ``graphs()`` will iterate through all graphs. The complete graph of 4
vertices is of course the smallest graph with chromatic number bigger than
three::

    sage: for g in graphs():
    ....:     if g.chromatic_number() > 3:
    ....:         break
    sage: g.is_isomorphic(graphs.CompleteGraph(4))
    True

For some commonly used graphs to play with, type::

    sage: graphs.[tab]          # not tested

and hit {tab}. Most of these graphs come with their own custom plot, so you can
see how people usually visualize these graphs.

::

    sage: G = graphs.PetersenGraph()
    sage: G.plot().show()    # or G.show()
    sage: G.degree_histogram()
    [0, 0, 0, 10]
    sage: G.adjacency_matrix()
    [0 1 0 0 1 1 0 0 0 0]
    [1 0 1 0 0 0 1 0 0 0]
    [0 1 0 1 0 0 0 1 0 0]
    [0 0 1 0 1 0 0 0 1 0]
    [1 0 0 1 0 0 0 0 0 1]
    [1 0 0 0 0 0 0 1 1 0]
    [0 1 0 0 0 0 0 0 1 1]
    [0 0 1 0 0 1 0 0 0 1]
    [0 0 0 1 0 1 1 0 0 0]
    [0 0 0 0 1 0 1 1 0 0]

::

    sage: S = G.subgraph([0,1,2,3])
    sage: S.plot().show()    # or S.show()
    sage: S.density()
    1/2

::

    sage: G = GraphQuery(display_cols=['graph6'], num_vertices=7, diameter=5)
    sage: L = G.get_graphs_list()
    sage: graphs_list.show_graphs(L)

.. _Graph:labels:

Labels
------

Each vertex can have any hashable object as a label. These are things like
strings, numbers, and tuples. Each edge is given a default label of ``None``,
but if specified, edges can have any label at all. Edges between vertices `u`
and `v` are represented typically as ``(u, v, l)``, where ``l`` is the label for
the edge.

Note that vertex labels themselves cannot be mutable items::

    sage: M = Matrix( [[0,0],[0,0]] )
    sage: G = Graph({ 0 : { M : None } })
    Traceback (most recent call last):
    ...
    TypeError: mutable matrices are unhashable

However, if one wants to define a dictionary, with the same keys and arbitrary
objects for entries, one can make that association::

    sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), \
          2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
    sage: d[2]
    Moebius-Kantor Graph: Graph on 16 vertices
    sage: T = graphs.TetrahedralGraph()
    sage: T.vertices()
    [0, 1, 2, 3]
    sage: T.set_vertices(d)
    sage: T.get_vertex(1)
    Flower Snark: Graph on 20 vertices

Database
--------

There is a database available for searching for graphs that satisfy a certain
set of parameters, including number of vertices and edges, density, maximum and
minimum degree, diameter, radius, and connectivity. To see a list of all search
parameter keywords broken down by their designated table names, type ::

    sage: graph_db_info()
    {...}

For more details on data types or keyword input, enter ::

    sage: GraphQuery?    # not tested

The results of a query can be viewed with the show method, or can be viewed
individually by iterating through the results ::

    sage: Q = GraphQuery(display_cols=['graph6'],num_vertices=7, diameter=5)
    sage: Q.show()
    Graph6
    --------------------
    F?`po
    F?gqg
    F@?]O
    F@OKg
    F@R@o
    FA_pW
    FEOhW
    FGC{o
    FIAHo

Show each graph as you iterate through the results::

    sage: for g in Q:
    ....:     show(g)

Visualization
-------------

To see a graph `G` you are working with, there are three main options. You can
view the graph in two dimensions via matplotlib with ``show()``. ::

    sage: G = graphs.RandomGNP(15,.3)
    sage: G.show()

And you can view it in three dimensions via jmol with ``show3d()``. ::

    sage: G.show3d()

Or it can be rendered with `\LaTeX`.  This requires the right additions to a
standard `\mbox{\rm\TeX}` installation.  Then standard Sage commands, such as
``view(G)`` will display the graph, or ``latex(G)`` will produce a string
suitable for inclusion in a `\LaTeX` document.  More details on this are at the
:mod:`sage.graphs.graph_latex` module. ::

    sage: from sage.graphs.graph_latex import check_tkz_graph
    sage: check_tkz_graph()  # random - depends on TeX installation
    sage: latex(G)
    \begin{tikzpicture}
    ...
    \end{tikzpicture}

Mutability
----------

Graphs are mutable, and thus unusable as dictionary keys, unless
``data_structure="static_sparse"`` is used::

    sage: G = graphs.PetersenGraph()
    sage: {G:1}[G]
    Traceback (most recent call last):
    ...
    TypeError: This graph is mutable, and thus not hashable. Create an immutable copy by `g.copy(immutable=True)`
    sage: G_immutable = Graph(G, immutable=True)
    sage: G_immutable == G
    True
    sage: {G_immutable:1}[G_immutable]
    1

Methods
-------
"""


# ****************************************************************************
#       Copyright (C) 2006-2007 Robert L. Miller <rlmillster@gmail.com>
#                          2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import itertools

from copy import copy
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
import sage.graphs.generic_graph_pyx as generic_graph_pyx
from sage.graphs.generic_graph import GenericGraph
from sage.graphs.digraph import DiGraph
from sage.graphs.independent_sets import IndependentSets
from sage.misc.rest_index_of_methods import doc_index, gen_thematic_rest_table_index
from sage.graphs.views import EdgesView

from sage.misc.lazy_import import lazy_import
from sage.features import PythonModule
lazy_import('sage.graphs.mcqd', ['mcqd'],
            feature=PythonModule('sage.graphs.mcqd', spkg='mcqd'))

from sage.misc.decorators import rename_keyword


class Graph(GenericGraph):
    r"""
    Undirected graph.

    A graph is a set of vertices connected by edges. See the
    :wikipedia:`Graph_(mathematics)` for more information. For a collection of
    pre-defined graphs, see the :mod:`~sage.graphs.graph_generators` module.

    A :class:`Graph` object has many methods whose list can be obtained by
    typing ``g.<tab>`` (i.e. hit the 'tab' key) or by reading the documentation
    of :mod:`~sage.graphs.graph`, :mod:`~sage.graphs.generic_graph`, and
    :mod:`~sage.graphs.digraph`.

    INPUT:

    By default, a :class:`Graph` object is simple (i.e. no *loops* nor *multiple
    edges*) and unweighted. This can be easily tuned with the appropriate flags
    (see below).

    - ``data`` -- can be any of the following (see the ``format`` argument):

      #. ``Graph()`` -- build a graph on 0 vertices.

      #. ``Graph(5)`` -- return an edgeless graph on the 5 vertices 0,...,4.

      #. ``Graph([list_of_vertices, list_of_edges])`` -- returns a graph with
         given vertices/edges.

         To bypass auto-detection, prefer the more explicit
         ``Graph([V, E], format='vertices_and_edges')``.

      #. ``Graph(list_of_edges)`` -- return a graph with a given list of edges
         (see documentation of
         :meth:`~sage.graphs.generic_graph.GenericGraph.add_edges`).

         To bypass auto-detection, prefer the more explicit
         ``Graph(L, format='list_of_edges')``.

      #. ``Graph({1: [2, 3, 4], 3: [4]})`` -- return a graph by associating to
         each vertex the list of its neighbors.

         To bypass auto-detection, prefer the more explicit
         ``Graph(D, format='dict_of_lists')``.

      #. ``Graph({1: {2: 'a', 3:'b'} ,3:{2:'c'}})`` -- return a graph by
         associating a list of neighbors to each vertex and providing its edge
         label.

         To bypass auto-detection, prefer the more explicit
         ``Graph(D, format='dict_of_dicts')``.

         For graphs with multiple edges, you can provide a list of labels
         instead, e.g.: ``Graph({1: {2: ['a1', 'a2'], 3:['b']} ,3:{2:['c']}})``.

      #. ``Graph(a_symmetric_matrix)`` -- return a graph with given (weighted)
         adjacency matrix (see documentation of
         :meth:`~sage.graphs.generic_graph.GenericGraph.adjacency_matrix`).

         To bypass auto-detection, prefer the more explicit ``Graph(M,
         format='adjacency_matrix')``. To take weights into account, use
         ``format='weighted_adjacency_matrix'`` instead.

      #. ``Graph(a_nonsymmetric_matrix)`` -- return a graph with given incidence
         matrix (see documentation of
         :meth:`~sage.graphs.generic_graph.GenericGraph.incidence_matrix`).

         To bypass auto-detection, prefer the more explicit
         ``Graph(M, format='incidence_matrix')``.

      #. ``Graph([V, f])`` -- return a graph from a vertex set ``V`` and a
         *symmetric* function ``f``. The graph contains an edge `u,v` whenever
         ``f(u,v)`` is ``True``.. Example: ``Graph([ [1..10], lambda x,y:
         abs(x-y).is_square()])``

      #. ``Graph(':I`ES@obGkqegW~')`` -- return a graph from a graph6 or sparse6
         string (see documentation of :meth:`graph6_string` or
         :meth:`sparse6_string`).

      #. ``Graph(a_seidel_matrix, format='seidel_adjacency_matrix')`` -- return
         a graph with a given Seidel adjacency matrix (see documentation of
         :meth:`seidel_adjacency_matrix`).

      #. ``Graph(another_graph)`` -- return a graph from a Sage (di)graph,
         `pygraphviz <https://pygraphviz.github.io/>`__ graph, `NetworkX
         <https://networkx.github.io/>`__ graph, or `igraph
         <http://igraph.org/python/>`__ graph.

    - ``pos`` -- a positioning dictionary (cf. documentation of
      :meth:`~sage.graphs.generic_graph.GenericGraph.layout`). For example, to
      draw 4 vertices on a square::

         {0: [-1,-1],
          1: [ 1,-1],
          2: [ 1, 1],
          3: [-1, 1]}

    - ``name`` -- (must be an explicitly named parameter, i.e.,
       ``name="complete")`` gives the graph a name

    - ``loops`` -- boolean (default: ``None``); whether to allow loops (ignored
       if data is an instance of the ``Graph`` class)

    - ``multiedges`` -- boolean (default: ``None``); whether to allow multiple
       edges (ignored if data is an instance of the ``Graph`` class).

    - ``weighted`` -- boolean (default: ``None``); whether graph thinks of
      itself as weighted or not. See
      :meth:`~sage.graphs.generic_graph.GenericGraph.weighted`.

    - ``format`` -- if set to ``None`` (default), :class:`Graph` tries to guess
      input's format. To avoid this possibly time-consuming step, one of the
      following values can be specified (see description above): ``"int"``,
      ``"graph6"``, ``"sparse6"``, ``"rule"``, ``"list_of_edges"``,
      ``"dict_of_lists"``, ``"dict_of_dicts"``, ``"adjacency_matrix"``,
      ``"weighted_adjacency_matrix"``, ``"seidel_adjacency_matrix"``,
      ``"incidence_matrix"``, ``"NX"``, ``"igraph"``.

    - ``sparse`` -- boolean (default: ``True``); ``sparse=True`` is an alias for
      ``data_structure="sparse"``, and ``sparse=False`` is an alias for
      ``data_structure="dense"``.

    - ``data_structure`` -- one of the following (for more information, see
      :mod:`~sage.graphs.base.overview`)

       * ``"dense"`` -- selects the :mod:`~sage.graphs.base.dense_graph`
         backend.

       * ``"sparse"`` -- selects the :mod:`~sage.graphs.base.sparse_graph`
         backend.

       * ``"static_sparse"`` -- selects the
         :mod:`~sage.graphs.base.static_sparse_backend` (this backend is faster
         than the sparse backend and smaller in memory, and it is immutable, so
         that the resulting graphs can be used as dictionary keys).

    - ``immutable`` -- boolean (default: ``False``); whether to create a
      immutable graph. Note that ``immutable=True`` is actually a shortcut for
      ``data_structure='static_sparse'``. Set to ``False`` by default.

    - ``vertex_labels`` -- boolean (default: ``True``); whether to allow any
      object as a vertex (slower), or only the integers `0,...,n-1`, where `n`
      is the number of vertices.

    - ``convert_empty_dict_labels_to_None`` -- this arguments sets the default
       edge labels used by NetworkX (empty dictionaries) to be replaced by
       ``None``, the default Sage edge label. It is set to ``True`` iff a
       NetworkX graph is on the input.

    EXAMPLES:

    We illustrate the first seven input formats (the other two involve packages
    that are currently not standard in Sage):

    #. An integer giving the number of vertices::

        sage: g = Graph(5); g
        Graph on 5 vertices
        sage: g.vertices()
        [0, 1, 2, 3, 4]
        sage: g.edges()
        []

    #. A dictionary of dictionaries::

        sage: g = Graph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}); g
        Graph on 5 vertices

       The labels ('x', 'z', 'a', 'out') are labels for edges. For example,
       'out' is the label for the edge on 2 and 5. Labels can be used as
       weights, if all the labels share some common parent.::

        sage: a, b, c, d, e, f = sorted(SymmetricGroup(3))              # optional - sage.groups
        sage: Graph({b: {d: 'c', e: 'p'}, c: {d: 'p', e: 'c'}})         # optional - sage.groups
        Graph on 4 vertices

    #. A dictionary of lists::

        sage: g = Graph({0:[1,2,3], 2:[4]}); g
        Graph on 5 vertices

    #. A list of vertices and a function describing adjacencies. Note that the
       list of vertices and the function must be enclosed in a list (i.e., [list
       of vertices, function]).

       Construct the Paley graph over GF(13).::

          sage: g=Graph([GF(13), lambda i,j: i!=j and (i-j).is_square()])
          sage: g.vertices()
          [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
          sage: g.adjacency_matrix()
          [0 1 0 1 1 0 0 0 0 1 1 0 1]
          [1 0 1 0 1 1 0 0 0 0 1 1 0]
          [0 1 0 1 0 1 1 0 0 0 0 1 1]
          [1 0 1 0 1 0 1 1 0 0 0 0 1]
          [1 1 0 1 0 1 0 1 1 0 0 0 0]
          [0 1 1 0 1 0 1 0 1 1 0 0 0]
          [0 0 1 1 0 1 0 1 0 1 1 0 0]
          [0 0 0 1 1 0 1 0 1 0 1 1 0]
          [0 0 0 0 1 1 0 1 0 1 0 1 1]
          [1 0 0 0 0 1 1 0 1 0 1 0 1]
          [1 1 0 0 0 0 1 1 0 1 0 1 0]
          [0 1 1 0 0 0 0 1 1 0 1 0 1]
          [1 0 1 1 0 0 0 0 1 1 0 1 0]

       Construct the line graph of a complete graph.::

          sage: g=graphs.CompleteGraph(4)
          sage: line_graph=Graph([g.edges(labels=false), \
                 lambda i,j: len(set(i).intersection(set(j)))>0], \
                 loops=False)
          sage: line_graph.vertices()
          [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
          sage: line_graph.adjacency_matrix()
          [0 1 1 1 1 0]
          [1 0 1 1 0 1]
          [1 1 0 0 1 1]
          [1 1 0 0 1 1]
          [1 0 1 1 0 1]
          [0 1 1 1 1 0]

    #. A graph6 or sparse6 string: Sage automatically recognizes whether a
       string is in graph6 or sparse6 format::

           sage: s = ':I`AKGsaOs`cI]Gb~'
           sage: Graph(s,sparse=True)
           Looped multi-graph on 10 vertices

       ::

           sage: G = Graph('G?????')
           sage: G = Graph("G'?G?C")
           Traceback (most recent call last):
           ...
           RuntimeError: the string seems corrupt: valid characters are
           ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
           sage: G = Graph('G??????')
           Traceback (most recent call last):
           ...
           RuntimeError: the string (G??????) seems corrupt: for n = 8, the string is too long

       ::

          sage: G = Graph(":I'AKGsaOs`cI]Gb~")
          Traceback (most recent call last):
          ...
          RuntimeError: the string seems corrupt: valid characters are
          ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

       There are also list functions to take care of lists of graphs::

           sage: s = ':IgMoqoCUOqeb\n:I`AKGsaOs`cI]Gb~\n:I`EDOAEQ?PccSsge\\N\n'
           sage: graphs_list.from_sparse6(s)
           [Looped multi-graph on 10 vertices, Looped multi-graph on 10 vertices, Looped multi-graph on 10 vertices]

    #. A Sage matrix:
       Note: If format is not specified, then Sage assumes a symmetric square
       matrix is an adjacency matrix, otherwise an incidence matrix.

       - an adjacency matrix::

            sage: M = graphs.PetersenGraph().am(); M
            [0 1 0 0 1 1 0 0 0 0]
            [1 0 1 0 0 0 1 0 0 0]
            [0 1 0 1 0 0 0 1 0 0]
            [0 0 1 0 1 0 0 0 1 0]
            [1 0 0 1 0 0 0 0 0 1]
            [1 0 0 0 0 0 0 1 1 0]
            [0 1 0 0 0 0 0 0 1 1]
            [0 0 1 0 0 1 0 0 0 1]
            [0 0 0 1 0 1 1 0 0 0]
            [0 0 0 0 1 0 1 1 0 0]
            sage: Graph(M)
            Graph on 10 vertices

         ::

            sage: Graph(matrix([[1,2],[2,4]]),loops=True,sparse=True)
            Looped multi-graph on 2 vertices

            sage: M = Matrix([[0,1,-1],[1,0,-1/2],[-1,-1/2,0]]); M
            [   0    1   -1]
            [   1    0 -1/2]
            [  -1 -1/2    0]
            sage: G = Graph(M,sparse=True); G
            Graph on 3 vertices
            sage: G.weighted()
            True

       - an incidence matrix::

            sage: M = Matrix(6, [-1,0,0,0,1, 1,-1,0,0,0, 0,1,-1,0,0, 0,0,1,-1,0, 0,0,0,1,-1, 0,0,0,0,0]); M
            [-1  0  0  0  1]
            [ 1 -1  0  0  0]
            [ 0  1 -1  0  0]
            [ 0  0  1 -1  0]
            [ 0  0  0  1 -1]
            [ 0  0  0  0  0]
            sage: Graph(M)
            Graph on 6 vertices

            sage: Graph(Matrix([[1],[1],[1]]))
            Traceback (most recent call last):
            ...
            ValueError: there must be one or two nonzero entries per column in an incidence matrix, got entries [1, 1, 1] in column 0
            sage: Graph(Matrix([[1],[1],[0]]))
            Graph on 3 vertices

            sage: M = Matrix([[0,1,-1],[1,0,-1],[-1,-1,0]]); M
            [ 0  1 -1]
            [ 1  0 -1]
            [-1 -1  0]
            sage: Graph(M,sparse=True)
            Graph on 3 vertices

            sage: M = Matrix([[0,1,1],[1,0,1],[-1,-1,0]]); M
            [ 0  1  1]
            [ 1  0  1]
            [-1 -1  0]
            sage: Graph(M)
            Traceback (most recent call last):
            ...
            ValueError: there must be one or two nonzero entries per column in an incidence matrix, got entries [1, 1] in column 2

        Check that :trac:`9714` is fixed::

            sage: MA = Matrix([[1,2,0], [0,2,0], [0,0,1]])
            sage: GA = Graph(MA, format='adjacency_matrix')
            sage: MI = GA.incidence_matrix(oriented=False)
            sage: MI
            [2 1 1 0 0 0]
            [0 1 1 2 2 0]
            [0 0 0 0 0 2]
            sage: Graph(MI).edges(labels=None)
            [(0, 0), (0, 1), (0, 1), (1, 1), (1, 1), (2, 2)]

            sage: M = Matrix([[1], [-1]]); M
            [ 1]
            [-1]
            sage: Graph(M).edges()
            [(0, 1, None)]

    #. A Seidel adjacency matrix::

          sage: from sage.combinat.matrices.hadamard_matrix import \
          ....:  regular_symmetric_hadamard_matrix_with_constant_diagonal as rshcd
          sage: m=rshcd(16,1)- matrix.identity(16)
          sage: Graph(m,format="seidel_adjacency_matrix").is_strongly_regular(parameters=True)
          (16, 6, 2, 2)

    #. List of edges, or labelled edges::

          sage: g = Graph([(1,3),(3,8),(5,2)])
          sage: g
          Graph on 5 vertices

          sage: g = Graph([(1,2,"Peace"),(7,-9,"and"),(77,2, "Love")])
          sage: g
          Graph on 5 vertices
          sage: g = Graph([(0, 2, '0'), (0, 2, '1'), (3, 3, '2')], loops=True, multiedges=True)
          sage: g.loops()
          [(3, 3, '2')]

    #. A NetworkX MultiGraph::

          sage: import networkx
          sage: g = networkx.MultiGraph({0:[1,2,3], 2:[4]})
          sage: Graph(g)
          Multi-graph on 5 vertices

    #. A NetworkX graph::

           sage: import networkx
           sage: g = networkx.Graph({0:[1,2,3], 2:[4]})
           sage: DiGraph(g)
           Digraph on 5 vertices

    #. An igraph Graph (see also
       :meth:`~sage.graphs.generic_graph.GenericGraph.igraph_graph`)::

           sage: import igraph                      # optional - python_igraph
           sage: g = igraph.Graph([(0, 1), (0, 2)]) # optional - python_igraph
           sage: Graph(g)                           # optional - python_igraph
           Graph on 3 vertices

       If ``vertex_labels`` is ``True``, the names of the vertices are given by
       the vertex attribute ``'name'``, if available::

           sage: g = igraph.Graph([(0,1),(0,2)], vertex_attrs={'name':['a','b','c']})  # optional - python_igraph
           sage: Graph(g).vertices()                                                   # optional - python_igraph
           ['a', 'b', 'c']
           sage: g = igraph.Graph([(0,1),(0,2)], vertex_attrs={'label':['a','b','c']}) # optional - python_igraph
           sage: Graph(g).vertices()                                                   # optional - python_igraph
           [0, 1, 2]

       If the igraph Graph has edge attributes, they are used as edge labels::

           sage: g = igraph.Graph([(0,1),(0,2)], edge_attrs={'name':['a','b'], 'weight':[1,3]}) # optional - python_igraph
           sage: Graph(g).edges()                                                               # optional - python_igraph
           [(0, 1, {'name': 'a', 'weight': 1}), (0, 2, {'name': 'b', 'weight': 3})]


    When defining an undirected graph from a function ``f``, it is *very*
    important that ``f`` be symmetric. If it is not, anything can happen::

        sage: f_sym = lambda x,y: abs(x-y) == 1
        sage: f_nonsym = lambda x,y: (x-y) == 1
        sage: G_sym = Graph([[4,6,1,5,3,7,2,0], f_sym])
        sage: G_sym.is_isomorphic(graphs.PathGraph(8))
        True
        sage: G_nonsym = Graph([[4,6,1,5,3,7,2,0], f_nonsym])
        sage: G_nonsym.size()
        4
        sage: G_nonsym.is_isomorphic(G_sym)
        False

    By default, graphs are mutable and can thus not be used as a dictionary
    key::

          sage: G = graphs.PetersenGraph()
          sage: {G:1}[G]
          Traceback (most recent call last):
          ...
          TypeError: This graph is mutable, and thus not hashable. Create an immutable copy by `g.copy(immutable=True)`

    When providing the optional arguments ``data_structure="static_sparse"`` or
    ``immutable=True`` (both mean the same), then an immutable graph results. ::

          sage: G_imm = Graph(G, immutable=True)
          sage: H_imm = Graph(G, data_structure='static_sparse')
          sage: G_imm == H_imm == G
          True
          sage: {G_imm:1}[H_imm]
          1

    TESTS::

        sage: Graph(4, format="HeyHeyHey")
        Traceback (most recent call last):
        ...
        ValueError: Unknown input format 'HeyHeyHey'

        sage: Graph(igraph.Graph(directed=True)) # optional - python_igraph
        Traceback (most recent call last):
        ...
        ValueError: An *undirected* igraph graph was expected. To build an directed graph, call the DiGraph constructor.

        sage: m = matrix([[0, -1], [-1, 0]])
        sage: Graph(m, format="seidel_adjacency_matrix")
        Graph on 2 vertices
        sage: m[0,1] = 1
        sage: Graph(m, format="seidel_adjacency_matrix")
        Traceback (most recent call last):
        ...
        ValueError: the adjacency matrix of a Seidel graph must be symmetric

        sage: m[0,1] = -1; m[1,1] = 1
        sage: Graph(m, format="seidel_adjacency_matrix")
        Traceback (most recent call last):
        ...
        ValueError: the adjacency matrix of a Seidel graph must have 0s on the main diagonal

    From a list of vertices and a list of edges::

        sage: G = Graph([[1,2,3], [(1,2)]]); G
        Graph on 3 vertices
        sage: G.edges()
        [(1, 2, None)]

    Check that :trac:`27505` is fixed::

        sage: Graph(Graph().networkx_graph(), weighted=None, format='NX')
        Graph on 0 vertices
    """
    _directed = False

    def __init__(self, data=None, pos=None, loops=None, format=None,
                 weighted=None, data_structure="sparse",
                 vertex_labels=True, name=None,
                 multiedges=None, convert_empty_dict_labels_to_None=None,
                 sparse=True, immutable=False):
        """
        TESTS::

            sage: G = Graph()
            sage: loads(dumps(G)) == G
            True
            sage: a = matrix(2,2,[1,0,0,1])
            sage: Graph(a).adjacency_matrix() == a
            True

            sage: a = matrix(2,2,[2,0,0,1])
            sage: Graph(a,sparse=True).adjacency_matrix() == a
            True

        The positions are copied when the graph is built from another graph ::

            sage: g = graphs.PetersenGraph()
            sage: h = Graph(g)
            sage: g.get_pos() == h.get_pos()
            True

        The position dictionary is not the input one (:trac:`22424`)::

            sage: my_pos = {0:(0,0), 1:(1,1)}
            sage: G = Graph([[0,1], [(0,1)]], pos=my_pos)
            sage: my_pos == G._pos
            True
            sage: my_pos is G._pos
            False

        Or from a DiGraph ::

            sage: d = DiGraph(g)
            sage: h = Graph(d)
            sage: g.get_pos() == h.get_pos()
            True

        Loops are not counted as multiedges (see :trac:`11693`) and edges are
        not counted twice ::

            sage: Graph({1:[1]}).num_edges()
            1
            sage: Graph({1:[2,2]}).num_edges()
            2

        An empty list or dictionary defines a simple graph
        (:trac:`10441` and :trac:`12910`)::

            sage: Graph([])
            Graph on 0 vertices
            sage: Graph({})
            Graph on 0 vertices
            sage: # not "Multi-graph on 0 vertices"

        Verify that the int format works as expected (:trac:`12557`)::

            sage: Graph(2).adjacency_matrix()
            [0 0]
            [0 0]
            sage: Graph(3) == Graph(3,format='int')
            True

        Problem with weighted adjacency matrix (:trac:`13919`)::

            sage: B = {0:{1:2,2:5,3:4},1:{2:2,4:7},2:{3:1,4:4,5:3},3:{5:4},4:{5:1,6:5},5:{6:7}}
            sage: grafo3 = Graph(B,weighted=True)
            sage: matad = grafo3.weighted_adjacency_matrix()
            sage: grafo4 = Graph(matad,format = "adjacency_matrix", weighted=True)
            sage: grafo4.shortest_path(0,6,by_weight=True)
            [0, 1, 2, 5, 4, 6]

        Graphs returned when setting ``immutable=False`` are mutable::

            sage: g = graphs.PetersenGraph()
            sage: g = Graph(g.edges(),immutable=False)
            sage: g.add_edge("Hey", "Heyyyyyyy")

        And their name is set::

            sage: g = graphs.PetersenGraph()
            sage: Graph(g, immutable=True)
            Petersen graph: Graph on 10 vertices

        Check error messages for graphs built from incidence matrices (see
        :trac:`18440`)::

            sage: Graph(matrix([[-1, 1, 0],[1, 0, 0]]))
            Traceback (most recent call last):
            ...
            ValueError: column 1 of the (oriented) incidence matrix contains only one nonzero value
            sage: Graph(matrix([[1,1],[1,1],[1,0]]))
            Traceback (most recent call last):
            ...
            ValueError: there must be one or two nonzero entries per column in an incidence matrix, got entries [1, 1, 1] in column 0
            sage: Graph(matrix([[3,1,1],[0,1,1]]))
            Traceback (most recent call last):
            ...
            ValueError: each column of a non-oriented incidence matrix must sum to 2, but column 0 does not

        Vertex labels are retained in the graph (:trac:`14708`)::

            sage: g = Graph()
            sage: g.add_vertex(0)
            sage: g.set_vertex(0, 'foo')
            sage: g.get_vertices()
            {0: 'foo'}
            sage: Graph(g).get_vertices()
            {0: 'foo'}
        """
        GenericGraph.__init__(self)

        from sage.structure.element import is_Matrix

        if sparse is False:
            if data_structure != "sparse":
                raise ValueError("The 'sparse' argument is an alias for "
                                 "'data_structure'. Please do not define both.")
            data_structure = "dense"

        if multiedges or weighted:
            if data_structure == "dense":
                raise RuntimeError("Multiedge and weighted c_graphs must be sparse.")
        if immutable:
            data_structure = 'static_sparse'

        # If the data structure is static_sparse, we first build a graph
        # using the sparse data structure, then re-encode the resulting graph
        # as a static sparse graph.
        from sage.graphs.base.sparse_graph import SparseGraphBackend
        from sage.graphs.base.dense_graph import DenseGraphBackend
        if data_structure in ["sparse", "static_sparse"]:
            CGB = SparseGraphBackend
        elif data_structure == "dense":
            CGB = DenseGraphBackend
        else:
            raise ValueError("data_structure must be equal to 'sparse', "
                             "'static_sparse' or 'dense'")
        self._backend = CGB(0, directed=False)

        if format is None and isinstance(data, str):
            if data.startswith(">>graph6<<"):
                data = data[10:]
                format = 'graph6'
            elif data.startswith(">>sparse6<<"):
                data = data[11:]
                format = 'sparse6'
            elif data[0] == ':':
                format = 'sparse6'
            else:
                format = 'graph6'
        if format is None and is_Matrix(data):
            if data.is_symmetric():
                format = 'adjacency_matrix'
            else:
                format = 'incidence_matrix'
        if format is None and isinstance(data, Graph):
            format = 'Graph'
        from sage.graphs.digraph import DiGraph
        if format is None and isinstance(data, DiGraph):
            data = data.to_undirected()
            format = 'Graph'
        if (format is None         and
            isinstance(data, list) and
            len(data) >= 2         and
            callable(data[1])):
            format = 'rule'

        if (format is None            and
            isinstance(data, list)    and
            len(data) == 2            and
            isinstance(data[0], list) and    # a list of two lists, the second of
            ((isinstance(data[1], list) and  # which contains iterables (the edges)
              (not data[1] or callable(getattr(data[1][0], "__iter__", None)))) or
             (isinstance(data[1], EdgesView)))):
            format = "vertices_and_edges"

        if format is None and isinstance(data, dict):
            if not data:
                format = 'dict_of_dicts'
            else:
                val = next(iter(data.values()))
                if isinstance(val, (list, EdgesView)):
                    format = 'dict_of_lists'
                elif isinstance(val, dict):
                    format = 'dict_of_dicts'
        if format is None and hasattr(data, 'adj'):
            # the input is a networkx (Multi)(Di)Graph
            format = 'NX'

        if (format is None          and
            hasattr(data, 'vcount') and
            hasattr(data, 'get_edgelist')):
            try:
                import igraph
            except ImportError:
                raise ImportError("The data seems to be a igraph object, but "+
                                  "igraph is not installed in Sage. To install "+
                                  "it, run 'sage -i python_igraph'")
            if format is None and isinstance(data, igraph.Graph):
                format = 'igraph'
        if format is None and isinstance(data, (int, Integer)):
            format = 'int'
        if format is None and data is None:
            format = 'int'
            data = 0

        # Input is a list of edges or an EdgesView
        if format is None and isinstance(data, (list, EdgesView)):
            format = "list_of_edges"
            if weighted is None:
                weighted = False

        if format is None:
            raise ValueError("This input cannot be turned into a graph")

        if format == 'weighted_adjacency_matrix':
            if weighted is False:
                raise ValueError("Format was weighted_adjacency_matrix but weighted was False.")
            if weighted is None:
                weighted = True
            if multiedges is None:
                multiedges = False
            format = 'adjacency_matrix'

        # At this point, 'format' has been set. We build the graph

        if format == 'graph6':
            if weighted is None:
                weighted = False
            self.allow_loops(loops if loops else False, check=False)
            self.allow_multiple_edges(multiedges if multiedges else False, check=False)
            from .graph_input import from_graph6
            from_graph6(self, data)

        elif format == 'sparse6':
            if weighted is None:
                weighted = False
            self.allow_loops(False if loops is False else True, check=False)
            self.allow_multiple_edges(False if multiedges is False else True, check=False)
            from .graph_input import from_sparse6
            from_sparse6(self, data)

        elif format == 'adjacency_matrix':
            from .graph_input import from_adjacency_matrix
            from_adjacency_matrix(self, data, loops=loops, multiedges=multiedges, weighted=weighted)

        elif format == 'incidence_matrix':
            from .graph_input import from_incidence_matrix
            from_incidence_matrix(self, data, loops=loops, multiedges=multiedges, weighted=weighted)

        elif format == 'seidel_adjacency_matrix':
            weighted = False
            self.allow_loops(False)
            self.allow_multiple_edges(False)
            from .graph_input import from_seidel_adjacency_matrix
            from_seidel_adjacency_matrix(self, data)
        elif format == 'Graph':
            if loops is None:
                loops = data.allows_loops()
            if multiedges is None:
                multiedges = data.allows_multiple_edges()
            if weighted is None:
                weighted = data.weighted()
            self.allow_loops(loops, check=False)
            self.allow_multiple_edges(multiedges, check=False)
            if data.get_pos() is not None:
                pos = data.get_pos()
            self.name(data.name())
            self.set_vertices(data.get_vertices())
            data._backend.subgraph_given_vertices(self._backend, data)

        elif format == 'NX':
            from sage.graphs.graph_input import from_networkx_graph
            from_networkx_graph(self, data,
                                weighted=weighted, multiedges=multiedges, loops=loops,
                                convert_empty_dict_labels_to_None=convert_empty_dict_labels_to_None)
            if weighted is None:
                weighted = self.allows_multiple_edges()

        elif format == 'igraph':
            if data.is_directed():
                raise ValueError("An *undirected* igraph graph was expected. "+
                                 "To build an directed graph, call the DiGraph "+
                                 "constructor.")

            self.add_vertices(range(data.vcount()))
            self.add_edges((e.source, e.target, e.attributes()) for e in data.es())

            if vertex_labels and 'name' in data.vertex_attributes():
                vs = data.vs()
                self.relabel({v:vs[v]['name'] for v in self})

        elif format == 'rule':
            f = data[1]
            verts = data[0]
            if loops is None:
                loops = any(f(v,v) for v in verts)
            if weighted is None:
                weighted = False
            self.allow_loops(loops, check=False)
            self.allow_multiple_edges(True if multiedges else False, check=False)
            self.add_vertices(verts)
            self.add_edges(e for e in itertools.combinations(verts,2) if f(*e))
            if loops:
                self.add_edges((v,v) for v in verts if f(v,v))

        elif format == "vertices_and_edges":
            self.allow_multiple_edges(bool(multiedges), check=False)
            self.allow_loops(bool(loops), check=False)
            self.add_vertices(data[0])
            self.add_edges(data[1])

        elif format == 'dict_of_dicts':
            from .graph_input import from_dict_of_dicts
            from_dict_of_dicts(self, data, loops=loops, multiedges=multiedges, weighted=weighted,
                               convert_empty_dict_labels_to_None = False if convert_empty_dict_labels_to_None is None else convert_empty_dict_labels_to_None)

        elif format == 'dict_of_lists':
            from .graph_input import from_dict_of_lists
            from_dict_of_lists(self, data, loops=loops, multiedges=multiedges, weighted=weighted)

        elif format == 'int':
            self.allow_loops(loops if loops else False, check=False)
            self.allow_multiple_edges(multiedges if multiedges else False, check=False)
            if data < 0:
                raise ValueError("The number of vertices cannot be strictly negative!")
            if data:
                self.add_vertices(range(data))

        elif format == 'list_of_edges':
            self.allow_multiple_edges(True if multiedges else False,
                                      check=False)
            self.allow_loops(True if loops else False, check=False)
            self.add_edges(data)
        else:
            raise ValueError("Unknown input format '{}'".format(format))

        if weighted is None:
            weighted = False
        self._weighted = getattr(self, '_weighted', weighted)

        self._pos = copy(pos)

        if format != 'Graph' or name is not None:
            self.name(name)

        if data_structure == "static_sparse":
            from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            ib = StaticSparseBackend(self,
                                     loops = self.allows_loops(),
                                     multiedges = self.allows_multiple_edges())
            self._backend = ib
            self._immutable = True

    ### Formats

    @doc_index("Basic methods")
    def graph6_string(self):
        r"""
        Return the graph6 representation of the graph as an ASCII string.

        This is only valid for simple (no loops, no multiple edges) graphs
        on at most `2^{18}-1=262143` vertices.

        .. NOTE::

            As the graph6 format only handles graphs with vertex set
            `\{0,...,n-1\}`, a :meth:`relabelled copy
            <sage.graphs.generic_graph.GenericGraph.relabel>` will
            be encoded, if necessary.

        .. SEEALSO::

            * :meth:`~sage.graphs.digraph.DiGraph.dig6_string` --
              a similar string format for directed graphs

        EXAMPLES::

            sage: G = graphs.KrackhardtKiteGraph()
            sage: G.graph6_string()
            'IvUqwK@?G'

        TESTS::

            sage: Graph().graph6_string()
            '?'
        """
        n = self.order()
        if n > 262143:
            raise ValueError('graph6 format supports graphs on 0 to 262143 vertices only.')
        elif self.has_loops() or self.has_multiple_edges():
            raise ValueError('graph6 format supports only simple graphs (no loops, no multiple edges)')
        else:
            return generic_graph_pyx.small_integer_to_graph6(n) + generic_graph_pyx.binary_string_to_graph6(self._bit_vector())

    @doc_index("Basic methods")
    def sparse6_string(self):
        r"""
        Return the sparse6 representation of the graph as an ASCII string.

        Only valid for undirected graphs on 0 to 262143 vertices, but loops
        and multiple edges are permitted.

        .. NOTE::

            As the sparse6 format only handles graphs whose vertex set is
            `\{0,...,n-1\}`, a :meth:`relabelled copy
            <sage.graphs.generic_graph.GenericGraph.relabel>` of your graph will
            be encoded if necessary.

        EXAMPLES::

            sage: G = graphs.BullGraph()
            sage: G.sparse6_string()
            ':Da@en'

        ::

            sage: G = Graph(loops=True, multiedges=True, data_structure="sparse")
            sage: Graph(':?', data_structure="sparse") == G
            True

        TESTS::

            sage: G = Graph()
            sage: G.sparse6_string()
            ':?'

        Check that :trac:`18445` is fixed::

            sage: Graph(graphs.KneserGraph(5,2).sparse6_string()).size()
            15

        Graphs with 1 vertex are correctly handled (:trac:`24923`)::

            sage: Graph([(0, 0)], loops=True).sparse6_string()
            ':@^'
            sage: G = Graph(_)
            sage: G.order(), G.size()
            (1, 1)
            sage: Graph([(0, 0), (0, 0)], loops=True, multiedges=True).sparse6_string()
            ':@N'
            sage: H = Graph(_)
            sage: H.order(), H.size()
            (1, 2)

        Sparse6 encoding of canonical graph is unique (:trac:`31026`)::

            sage: G = Graph([(0,1),(1,2),(2,3),(3,0),(0,2)])
            sage: H = Graph([(0,1),(1,2),(2,3),(3,0),(1,3)])
            sage: G == H
            False
            sage: G.is_isomorphic(H)
            True
            sage: G.sparse6_string() == H.sparse6_string()
            False
            sage: G_ = G.canonical_label()
            sage: H_ = H.canonical_label()
            sage: G_ == H_
            True
            sage: G_.sparse6_string() == H_.sparse6_string()
            True

        The method can handle vertices with different types (:trac:`31026`)::

            sage: G = Graph([(1, 'a')])
            sage: H = Graph(G.sparse6_string())
            sage: G.is_isomorphic(H)
            True
            sage: set(G) == set(H)
            False
        """
        n = self.order()
        if not n:
            return ':?'
        if n > 262143:
            raise ValueError('sparse6 format supports graphs on 0 to 262143 vertices only.')
        if n == 1:
            s = '0' * self.size()
        else:
            try:
                V = sorted(self)
            except TypeError:
                V = self
            v_to_int = {v:i for i,v in enumerate(V)}
            edges = [sorted((v_to_int[u], v_to_int[v])) for u,v in self.edge_iterator(labels=False)]
            edges.sort(key=lambda e: (e[1], e[0])) # reverse lexicographic order

            # encode bit vector
            k = int((ZZ(n) - 1).nbits())
            v = 0
            i = 0
            m = 0
            s = ''
            while m < len(edges):
                if edges[m][1] > v + 1:
                    sp = generic_graph_pyx.int_to_binary_string(edges[m][1])
                    sp = '0'*(k-len(sp)) + sp
                    s += '1' + sp
                    v = edges[m][1]
                elif edges[m][1] == v + 1:
                    sp = generic_graph_pyx.int_to_binary_string(edges[m][0])
                    sp = '0'*(k-len(sp)) + sp
                    s += '1' + sp
                    v += 1
                    m += 1
                else:
                    sp = generic_graph_pyx.int_to_binary_string(edges[m][0])
                    sp = '0'*(k-len(sp)) + sp
                    s += '0' + sp
                    m += 1

        # encode s as a 6-string, as in R(x), but padding with 1's
        # pad on the right to make a multiple of 6
        s = s + ( '1' * ((6 - len(s))%6) )

        # split into groups of 6, and convert numbers to decimal, adding 63
        six_bits = ''
        for i in range(0, len(s), 6):
            six_bits += chr( int( s[i:i+6], 2) + 63 )
        return ':' + generic_graph_pyx.small_integer_to_graph6(n) + six_bits

    ### Attributes

    @doc_index("Basic methods")
    def is_directed(self):
        """
        Since graph is undirected, returns False.

        EXAMPLES::

            sage: Graph().is_directed()
            False
        """
        return False

    ### Properties
    @doc_index("Graph properties")
    def is_tree(self, certificate=False, output='vertex'):
        r"""
        Tests if the graph is a tree

        The empty graph is defined to be not a tree.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); whether to return a
          certificate. The method only returns boolean answers when
          ``certificate = False`` (default). When it is set to ``True``, it
          either answers ``(True, None)`` when the graph is a tree or ``(False,
          cycle)`` when it contains a cycle. It returns ``(False, None)`` when
          the graph is empty or not connected.

        - ``output`` -- either ``'vertex'`` (default) or ``'edge'``; whether the
          certificate is given as a list of vertices (``output = 'vertex'``) or
          a list of edges (``output = 'edge'``).

        When the certificate cycle is given as a list of edges, the edges are
        given as `(v_i, v_{i+1}, l)` where `v_1, v_2, \dots, v_n` are the
        vertices of the cycles (in their cyclic order).

        EXAMPLES::

            sage: all(T.is_tree() for T in graphs.trees(15))
            True

        With certificates::

            sage: g = graphs.RandomTree(30)
            sage: g.is_tree(certificate=True)
            (True, None)
            sage: g.add_edge(10,-1)
            sage: g.add_edge(11,-1)
            sage: isit, cycle = g.is_tree(certificate=True)
            sage: isit
            False
            sage: -1 in cycle
            True

        One can also ask for the certificate as a list of edges::

            sage: g = graphs.CycleGraph(4)
            sage: g.is_tree(certificate=True, output='edge')
            (False, [(3, 2, None), (2, 1, None), (1, 0, None), (0, 3, None)])

        This is useful for graphs with multiple edges::

            sage: G = Graph([(1, 2, 'a'), (1, 2, 'b')], multiedges=True)
            sage: G.is_tree(certificate=True)
            (False, [1, 2])
            sage: G.is_tree(certificate=True, output='edge')
            (False, [(1, 2, 'a'), (2, 1, 'b')])

        TESTS:

        :trac:`14434` is fixed::

            sage: g = Graph({0:[1,4,5],3:[4,8,9],4:[9],5:[7,8],7:[9]})
            sage: _,cycle = g.is_tree(certificate=True)
            sage: g.size()
            10
            sage: g.add_cycle(cycle)
            sage: g.size()
            10

        The empty graph::

            sage: graphs.EmptyGraph().is_tree()
            False
            sage: graphs.EmptyGraph().is_tree(certificate=True)
            (False, None)

        :trac:`22912` is fixed::

            sage: G = Graph([(0,0), (0,1)], loops=True)
            sage: G.is_tree(certificate=True)
            (False, [0])
            sage: G.is_tree(certificate=True, output='edge')
            (False, [(0, 0, None)])
        """
        if output not in ['vertex', 'edge']:
            raise ValueError('output must be either vertex or edge')

        if not self.order() or not self.is_connected():
            return (False, None) if certificate else False

        if certificate:
            if self.order() == self.size() + 1:
                return (True, None)

            if self.allows_loops():
                L = self.loop_edges() if output == 'edge' else self.loop_vertices()
                if L:
                    return False, L[:1]

            if self.has_multiple_edges():
                if output == 'vertex':
                    return (False, list(self.multiple_edges(sort=True)[0][:2]))
                edge1, edge2 = self.multiple_edges(sort=True)[:2]
                if edge1[0] != edge2[0]:
                    return (False, [edge1, edge2])
                return (False, [edge1, (edge2[1], edge2[0], edge2[2])])

            if output == 'edge':
                if self.allows_multiple_edges():
                    def vertices_to_edges(x):
                        return [(u[0], u[1], self.edge_label(u[0], u[1])[0])
                                for u in zip(x, x[1:] + [x[0]])]
                else:
                    def vertices_to_edges(x):
                        return [(u[0], u[1], self.edge_label(u[0], u[1]))
                                for u in zip(x, x[1:] + [x[0]])]

            # This code is a depth-first search that looks for a cycle in the
            # graph. We *know* it exists as there are too many edges around.
            seen = {}
            u = next(self.vertex_iterator())
            seen[u] = u
            stack = [(u, v) for v in self.neighbor_iterator(u)]
            while stack:
                u, v = stack.pop()
                if v in seen:
                    continue
                for w in self.neighbor_iterator(v):
                    if u == w:
                        continue
                    elif w in seen:
                        cycle = [w, v]
                        while u != w:
                            cycle.append(u)
                            u = seen[u]
                        cycle.reverse()
                        if output == 'vertex':
                            return (False, cycle)
                        return (False, vertices_to_edges(cycle))
                    else:
                        stack.append((v, w))
                seen[v] = u

        else:
            return self.order() == self.size() + 1

    @doc_index("Graph properties")
    def is_forest(self, certificate=False, output='vertex'):
        """
        Tests if the graph is a forest, i.e. a disjoint union of trees.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); whether to return a
          certificate. The method only returns boolean answers when
          ``certificate = False`` (default). When it is set to ``True``, it
          either answers ``(True, None)`` when the graph is a forest or
          ``(False, cycle)`` when it contains a cycle.

        - ``output`` -- either ``'vertex'`` (default) or ``'edge'``; whether the
          certificate is given as a list of vertices (``output = 'vertex'``) or
          a list of edges (``output = 'edge'``).

        EXAMPLES::

            sage: seven_acre_wood = sum(graphs.trees(7), Graph())
            sage: seven_acre_wood.is_forest()
            True

        With certificates::

            sage: g = graphs.RandomTree(30)
            sage: g.is_forest(certificate=True)
            (True, None)
            sage: (2*g + graphs.PetersenGraph() + g).is_forest(certificate=True)
            (False, [68, 66, 69, 67, 65])
        """
        connected_components = self.connected_components()
        number_of_connected_components = len(connected_components)
        isit = (self.order() ==
                self.size() + number_of_connected_components)

        if not certificate:
            return isit
        else:
            if isit:
                return (True, None)
            # The graph contains a cycle, and the user wants to see it.

            # No need to copy the graph
            if number_of_connected_components == 1:
                return self.is_tree(certificate=True, output=output)

            # We try to find a cycle in each connected component
            for cc in connected_components:
                isit, cycle = self.subgraph(cc).is_tree(certificate=True, output=output)
                if not isit:
                    return (False, cycle)

    @doc_index("Graph properties")
    def is_cactus(self):
        """
        Check whether the graph is cactus graph.

        A graph is called *cactus graph* if it is connected and every pair of
        simple cycles have at most one common vertex.

        There are other definitions, see the :wikipedia:`Cactus_graph`.

        EXAMPLES::

            sage: g = Graph({1: [2], 2: [3, 4], 3: [4, 5, 6, 7], 8: [3, 5], 9: [6, 7]})
            sage: g.is_cactus()
            True

            sage: c6 = graphs.CycleGraph(6)
            sage: naphthalene = c6 + c6
            sage: naphthalene.is_cactus()  # Not connected
            False
            sage: naphthalene.merge_vertices([0, 6])
            sage: naphthalene.is_cactus()
            True
            sage: naphthalene.merge_vertices([1, 7])
            sage: naphthalene.is_cactus()
            False

        TESTS::

            sage: all(graphs.PathGraph(i).is_cactus() for i in range(5))
            True

            sage: Graph('Fli@?').is_cactus()
            False

        Test a graph that is not outerplanar, see :trac:`24480`::

            sage: graphs.Balaban10Cage().is_cactus()
            False
        """
        self._scream_if_not_simple()

        # Special cases
        if self.order() < 4:
            return True

        if self.size() > 3 * (self.order() - 1) / 2:
            return False

        # Every cactus graph is outerplanar
        if not self.is_circular_planar():
            return False

        if not self.is_connected():
            return False

        # the number of faces is 1 plus the number of blocks of order > 2
        B = self.blocks_and_cut_vertices()[0]
        return len(self.faces()) == sum(1 for b in B if len(b) > 2) + 1

    @doc_index("Graph properties")
    def is_biconnected(self):
        """
        Test if the graph is biconnected.

        A biconnected graph is a connected graph on two or more vertices that is
        not broken into disconnected pieces by deleting any single vertex.

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.is_connected`
            - :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cut_vertices`
            - :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cuts_tree`
            - :wikipedia:`Biconnected_graph`

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G.is_biconnected()
            True
            sage: G.add_path([0,'a','b'])
            sage: G.is_biconnected()
            False
            sage: G.add_edge('b', 1)
            sage: G.is_biconnected()
            True

        TESTS::

            sage: Graph().is_biconnected()
            False
            sage: Graph(1).is_biconnected()
            False
            sage: graphs.CompleteGraph(2).is_biconnected()
            True
        """
        if self.order() < 2 or not self.is_connected():
            return False
        if self.blocks_and_cut_vertices()[1]:
            return False
        return True

    @doc_index("Graph properties")
    def is_block_graph(self):
        r"""
        Return whether this graph is a block graph.

        A block graph is a connected graph in which every biconnected component
        (block) is a clique.

        .. SEEALSO::

            - :wikipedia:`Block_graph` for more details on these graphs
            - :meth:`~sage.graphs.graph_generators.GraphGenerators.RandomBlockGraph`
              -- generator of random block graphs
            - :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cut_vertices`
            - :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cuts_tree`


        EXAMPLES::

            sage: G = graphs.RandomBlockGraph(6, 2, kmax=4)
            sage: G.is_block_graph()
            True
            sage: from sage.graphs.isgci import graph_classes
            sage: G in graph_classes.Block
            True
            sage: graphs.CompleteGraph(4).is_block_graph()
            True
            sage: graphs.RandomTree(6).is_block_graph()
            True
            sage: graphs.PetersenGraph().is_block_graph()
            False
            sage: Graph(4).is_block_graph()
            False
        """
        if not self.is_connected():
            return False
        if self.is_clique():
            return True

        B,C = self.blocks_and_cut_vertices()
        return all(self.is_clique(vertices=block) for block in B)

    @doc_index("Graph properties")
    def is_cograph(self):
        """
        Check whether the graph is cograph.

        A cograph is defined recursively: the single-vertex graph is
        cograph, complement of cograph is cograph, and disjoint union
        of two cographs is cograph. There are many other
        characterizations, see the :wikipedia:`Cograph`.

        EXAMPLES::

            sage: graphs.HouseXGraph().is_cograph()
            True
            sage: graphs.HouseGraph().is_cograph()
            False

        .. TODO::

            Implement faster recognition algorithm, as for instance
            the linear time recognition algorithm using LexBFS proposed
            in [Bre2008]_.

        TESTS::

            sage: [graphs.PathGraph(i).is_cograph() for i in range(6)]
            [True, True, True, True, False, False]
            sage: graphs.CycleGraph(5).is_cograph()  # Self-complemented
            False
        """
        # A cograph has no 4-vertex path as an induced subgraph.
        # We will first try to "decompose" graph by complements and
        # split to connected components, and use fairly slow
        # subgraph search if that fails.
        self._scream_if_not_simple()
        if self.order() < 4:
            return True
        if self.density()*2 > 1:
            return self.complement().is_cograph()
        if not self.is_connected():
            return all(part.is_cograph() for part in self.connected_components_subgraphs())
        P4 = Graph({0: [1], 1: [2], 2: [3]})
        return self.subgraph_search(P4, induced=True) is None

    @doc_index("Graph properties")
    def is_apex(self):
        r"""
        Test if the graph is apex.

        A graph is apex if it can be made planar by the removal of a single
        vertex. The deleted vertex is called ``an apex`` of the graph, and a
        graph may have more than one apex. For instance, in the minimal
        nonplanar graphs `K_5` or `K_{3,3}`, every vertex is an apex. The apex
        graphs include graphs that are themselves planar, in which case again
        every vertex is an apex. The null graph is also counted as an apex graph
        even though it has no vertex to remove.  If the graph is not connected,
        we say that it is apex if it has at most one non planar connected
        component and that this component is apex.  See the :wikipedia:`Apex_graph`
        for more information.

        .. SEEALSO::

          - :meth:`~Graph.apex_vertices`
          - :meth:`~sage.graphs.generic_graph.GenericGraph.is_planar`

        EXAMPLES:

        `K_5` and `K_{3,3}` are apex graphs, and each of their vertices is an
        apex::

            sage: G = graphs.CompleteGraph(5)
            sage: G.is_apex()
            True
            sage: G = graphs.CompleteBipartiteGraph(3,3)
            sage: G.is_apex()
            True

        The Petersen graph is not apex::

            sage: G = graphs.PetersenGraph()
            sage: G.is_apex()
            False

        A graph is apex if all its connected components are apex, but at most
        one is not planar::

            sage: M = graphs.Grid2dGraph(3,3)
            sage: K5 = graphs.CompleteGraph(5)
            sage: (M+K5).is_apex()
            True
            sage: (M+K5+K5).is_apex()
            False

        TESTS:

        The null graph is apex::

            sage: G = Graph()
            sage: G.is_apex()
            True

        The graph might be mutable or immutable::

            sage: G = Graph(M+K5, immutable=True)
            sage: G.is_apex()
            True
        """
        # Easy cases: null graph, subgraphs of K_5 and K_3,3
        if self.order() <= 5 or ( self.order() <= 6 and self.is_bipartite() ):
            return True

        return len(self.apex_vertices(k=1)) > 0

    @doc_index("Graph properties")
    def apex_vertices(self, k=None):
        r"""
        Return the list of apex vertices.

        A graph is apex if it can be made planar by the removal of a single
        vertex. The deleted vertex is called ``an apex`` of the graph, and a
        graph may have more than one apex. For instance, in the minimal
        nonplanar graphs `K_5` or `K_{3,3}`, every vertex is an apex. The apex
        graphs include graphs that are themselves planar, in which case again
        every vertex is an apex. The null graph is also counted as an apex graph
        even though it has no vertex to remove.  If the graph is not connected,
        we say that it is apex if it has at most one non planar connected
        component and that this component is apex.  See the
        :wikipedia:`Apex_graph` for more information.

        .. SEEALSO::

          - :meth:`~Graph.is_apex`
          - :meth:`~sage.graphs.generic_graph.GenericGraph.is_planar`

        INPUT:

        - ``k`` -- integer (default: ``None``); when set to ``None``, the method
          returns the list of all apex of the graph, possibly empty if the graph
          is not apex. When set to a positive integer, the method ends as soon
          as `k` apex vertices are found.

        OUTPUT:

        By default, the method returns the list of all apex of the graph. When
        parameter ``k`` is set to a positive integer, the returned list is
        bounded to `k` apex vertices.

        EXAMPLES:

        `K_5` and `K_{3,3}` are apex graphs, and each of their vertices is an
        apex::

            sage: G = graphs.CompleteGraph(5)
            sage: G.apex_vertices()
            [0, 1, 2, 3, 4]
            sage: G = graphs.CompleteBipartiteGraph(3,3)
            sage: G.is_apex()
            True
            sage: G.apex_vertices()
            [0, 1, 2, 3, 4, 5]
            sage: G.apex_vertices(k=3)
            [0, 1, 2]

        A `4\\times 4`-grid is apex and each of its vertices is an apex. When
        adding a universal vertex, the resulting graph is apex and the universal
        vertex is the unique apex vertex ::

            sage: G = graphs.Grid2dGraph(4,4)
            sage: set(G.apex_vertices()) == set(G.vertices())
            True
            sage: G.add_edges([('universal',v) for v in G])
            sage: G.apex_vertices()
            ['universal']

        The Petersen graph is not apex::

            sage: G = graphs.PetersenGraph()
            sage: G.apex_vertices()
            []

        A graph is apex if all its connected components are apex, but at most
        one is not planar::

            sage: M = graphs.Grid2dGraph(3,3)
            sage: K5 = graphs.CompleteGraph(5)
            sage: (M+K5).apex_vertices()
            [9, 10, 11, 12, 13]
            sage: (M+K5+K5).apex_vertices()
            []

        Neighbors of an apex of degree 2 are apex::

            sage: G = graphs.Grid2dGraph(5,5)
            sage: v = (666, 666)
            sage: G.add_path([(1, 1), v, (3, 3)])
            sage: G.is_planar()
            False
            sage: G.degree(v)
            2
            sage: sorted(G.apex_vertices())
            [(1, 1), (2, 2), (3, 3), (666, 666)]


        TESTS:

        The null graph is apex although it has no apex vertex::

            sage: G = Graph()
            sage: G.apex_vertices()
            []

        Parameter ``k`` cannot be a negative integer::

            sage: G.apex_vertices(k=-1)
            Traceback (most recent call last):
            ...
            ValueError: parameter k must be a non negative integer

        The graph might be mutable or immutable::

            sage: G = Graph(M+K5, immutable=True)
            sage: G.apex_vertices()
            [9, 10, 11, 12, 13]
        """
        if k is None:
            k = self.order()
        elif k < 0:
            raise ValueError("parameter k must be a non negative integer")

        # Easy cases: null graph, subgraphs of K_5 and K_3,3
        if self.order() <= 5 or (self.order() <= 6 and self.is_bipartite()):
            it = self.vertex_iterator()
            return [next(it) for _ in range(k)]


        if not self.is_connected():
            # We search for its non planar connected components. If it has more
            # than one such component, the graph is not apex. It is apex if
            # either it has no such component, in which case the graph is
            # planar, or if its unique non planar component is apex.

            P = [H for H in self.connected_components_subgraphs() if not H.is_planar()]
            if not P: # The graph is planar
                it = self.vertex_iterator()
                return [next(it) for _ in range(k)]
            elif len(P) > 1:
                return []
            else:
                # We proceed with the non planar component
                if P[0].is_immutable():
                    H = Graph(P[0].edges(labels=0, sort=False), immutable=False, loops=False, multiedges=False)
                else:
                    H = P[0]

        elif self.is_planar():
            # A planar graph is apex.
            it = self.vertex_iterator()
            return [next(it) for _ in range(k)]

        else:
            # We make a basic copy of the graph since we will modify it
            H = Graph(self.edges(labels=0, sort=False), immutable=False, loops=False, multiedges=False)


        # General case: basic implementation
        #
        # Test for each vertex if its removal makes the graph planar.
        # Obviously, we don't test vertices of degree one. Furthermore, if a
        # vertex of degree 2 is an apex, its neighbors also are. So we start
        # with vertices of degree 2.
        V = {}
        for u in H:
            d = H.degree(u)
            if d > 1:
                if d in V:
                    V[d].append(u)
                else:
                    V[d] = [u]
        apex = set()
        for deg in sorted(V):
            for u in V[deg]:
                if u in apex: # True if neighbor of an apex of degree 2
                    if deg == 2:
                        # We ensure that its neighbors are known apex
                        apex.update(H.neighbor_iterator(u))
                        if len(apex) >= k:
                            return list(apex)[:k]
                    continue

                E = H.edges_incident(u, labels=0)
                H.delete_vertex(u)
                if H.is_planar():
                    apex.add(u)
                    if deg == 2:
                        # The neighbors of an apex of degree 2 also are
                        apex.update(self.neighbor_iterator(u))

                    if len(apex) >= k:
                        return list(apex)[:k]

                H.add_edges(E)

        return list(apex)

    @doc_index("Graph properties")
    def is_overfull(self):
        r"""
        Tests whether the current graph is overfull.

        A graph `G` on `n` vertices and `m` edges is said to be overfull if:

        - `n` is odd

        - It satisfies `2m > (n-1)\Delta(G)`, where `\Delta(G)` denotes the
          maximum degree among all vertices in `G`.

        An overfull graph must have a chromatic index of `\Delta(G)+1`.

        EXAMPLES:

        A complete graph of order `n > 1` is overfull if and only if `n` is
        odd::

            sage: graphs.CompleteGraph(6).is_overfull()
            False
            sage: graphs.CompleteGraph(7).is_overfull()
            True
            sage: graphs.CompleteGraph(1).is_overfull()
            False

        The claw graph is not overfull::

            sage: from sage.graphs.graph_coloring import edge_coloring
            sage: g = graphs.ClawGraph()
            sage: g
            Claw graph: Graph on 4 vertices
            sage: edge_coloring(g, value_only=True)
            3
            sage: g.is_overfull()
            False

        The Holt graph is an example of a overfull graph::

            sage: G = graphs.HoltGraph()
            sage: G.is_overfull()
            True

        Checking that all complete graphs `K_n` for even `0 \leq n \leq 100`
        are not overfull::

            sage: def check_overfull_Kn_even(n):
            ....:     i = 0
            ....:     while i <= n:
            ....:         if graphs.CompleteGraph(i).is_overfull():
            ....:             print("A complete graph of even order cannot be overfull.")
            ....:             return
            ....:         i += 2
            ....:     print("Complete graphs of even order up to %s are not overfull." % n)
            ...
            sage: check_overfull_Kn_even(100)  # long time
            Complete graphs of even order up to 100 are not overfull.

        The null graph, i.e. the graph with no vertices, is not overfull::

            sage: Graph().is_overfull()
            False
            sage: graphs.CompleteGraph(0).is_overfull()
            False

        Checking that all complete graphs `K_n` for odd `1 < n \leq 100`
        are overfull::

            sage: def check_overfull_Kn_odd(n):
            ....:     i = 3
            ....:     while i <= n:
            ....:         if not graphs.CompleteGraph(i).is_overfull():
            ....:             print("A complete graph of odd order > 1 must be overfull.")
            ....:             return
            ....:         i += 2
            ....:     print("Complete graphs of odd order > 1 up to %s are overfull." % n)
            ...
            sage: check_overfull_Kn_odd(100)  # long time
            Complete graphs of odd order > 1 up to 100 are overfull.

        The Petersen Graph, though, is not overfull while
        its chromatic index is `\Delta+1`::

            sage: g = graphs.PetersenGraph()
            sage: g.is_overfull()
            False
            sage: from sage.graphs.graph_coloring import edge_coloring
            sage: max(g.degree()) + 1 ==  edge_coloring(g, value_only=True)
            True
        """
        # # A possible optimized version. But the gain in speed is very little.
        # return bool(self._backend.num_verts() & 1) and (  # odd order n
        #     2 * self._backend.num_edges(self._directed) > #2m > \Delta(G)*(n-1)
        #     max(self.degree()) * (self._backend.num_verts() - 1))
        # unoptimized version
        return (self.order() % 2 == 1) and (
            2 * self.size() > max(self.degree()) * (self.order() - 1))

    @doc_index("Graph properties")
    def is_even_hole_free(self, certificate=False):
        r"""
        Tests whether ``self`` contains an induced even hole.

        A Hole is a cycle of length at least 4 (included). It is said to be even
        (resp. odd) if its length is even (resp. odd).

        Even-hole-free graphs always contain a bisimplicial vertex, which
        ensures that their chromatic number is at most twice their clique number
        [ACHRS2008]_.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); when ``certificate =
          False``, this method only returns ``True`` or ``False``. If
          ``certificate = True``, the subgraph found is returned instead of
          ``False``.

        EXAMPLES:

        Is the Petersen Graph even-hole-free ::

            sage: g = graphs.PetersenGraph()
            sage: g.is_even_hole_free()
            False

        As any chordal graph is hole-free, interval graphs behave the same way::

            sage: g = graphs.RandomIntervalGraph(20)
            sage: g.is_even_hole_free()
            True

        It is clear, though, that a random Bipartite Graph which is not a forest
        has an even hole::

            sage: g = graphs.RandomBipartite(10, 10, .5)
            sage: g.is_even_hole_free() and not g.is_forest()
            False

        We can check the certificate returned is indeed an even cycle::

            sage: if not g.is_forest():
            ....:    cycle = g.is_even_hole_free(certificate=True)
            ....:    if cycle.order() % 2 == 1:
            ....:        print("Error !")
            ....:    if not cycle.is_isomorphic(
            ....:           graphs.CycleGraph(cycle.order())):
            ....:        print("Error !")
            ...
            sage: print("Everything is Fine !")
            Everything is Fine !

        TESTS:

        Bug reported in :trac:`9925`, and fixed by :trac:`9420`::

            sage: g = Graph(':SiBFGaCEF_@CE`DEGH`CEFGaCDGaCDEHaDEF`CEH`ABCDEF', loops=False, multiedges=False)
            sage: g.is_even_hole_free()
            False
            sage: g.is_even_hole_free(certificate=True)
            Subgraph of (): Graph on 4 vertices

        Making sure there are no other counter-examples around ::

            sage: t = lambda x: (Graph(x).is_forest() or
            ....:       isinstance(Graph(x).is_even_hole_free(certificate=True), Graph))
            sage: all( t(graphs.RandomBipartite(10, 10, .5)) for i in range(100) )
            True
        """
        girth = self.girth()

        if girth > self.order():
            start = 4

        elif not girth % 2:
            if not certificate:
                return False
            start = girth

        else:
            start = girth + 1

        from sage.graphs.generators.basic import CycleGraph

        while start <= self.order():

            subgraph = self.subgraph_search(CycleGraph(start), induced=True)

            if subgraph is not None:
                if certificate:
                    return subgraph
                else:
                    return False

            start += 2

        return True

    @doc_index("Graph properties")
    def is_odd_hole_free(self, certificate=False):
        r"""
        Tests whether ``self`` contains an induced odd hole.

        A Hole is a cycle of length at least 4 (included). It is said to be even
        (resp. odd) if its length is even (resp. odd).

        It is interesting to notice that while it is polynomial to check whether
        a graph has an odd hole or an odd antihole [CCLSV2005]_, it is not known
        whether testing for one of these two cases independently is polynomial
        too.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); when ``certificate =
          False``, this method only returns ``True`` or ``False``. If
          ``certificate = True``, the subgraph found is returned instead of
          ``False``.

        EXAMPLES:

        Is the Petersen Graph odd-hole-free ::

            sage: g = graphs.PetersenGraph()
            sage: g.is_odd_hole_free()
            False

        Which was to be expected, as its girth is 5 ::

            sage: g.girth()
            5

        We can check the certificate returned is indeed a 5-cycle::

            sage: cycle = g.is_odd_hole_free(certificate=True)
            sage: cycle.is_isomorphic(graphs.CycleGraph(5))
            True

        As any chordal graph is hole-free, no interval graph has an odd hole::

            sage: g = graphs.RandomIntervalGraph(20)
            sage: g.is_odd_hole_free()
            True
        """
        girth = self.odd_girth()

        if girth > self.order():
            return True
        if girth == 3:
            start = 5
        else:
            if not certificate:
                return False
            start = girth

        from sage.graphs.generators.basic import CycleGraph

        while start <= self.order():

            subgraph = self.subgraph_search(CycleGraph(start), induced=True)

            if subgraph is not None:
                if certificate:
                    return subgraph
                else:
                    return False

            start += 2

        return True

    @doc_index("Graph properties")
    def is_triangle_free(self, algorithm='dense_graph', certificate=False):
        r"""
        Check whether ``self`` is triangle-free

        INPUT:

        - ``algorithm`` -- (default: ``'dense_graph'``) specifies the algorithm
          to use among:

          - ``'matrix'`` -- tests if the trace of the adjacency matrix is
            positive.

          - ``'bitset'`` -- encodes adjacencies into bitsets and uses fast
            bitset operations to test if the input graph contains a
            triangle. This method is generally faster than standard matrix
            multiplication.

          - ``'dense_graph'`` -- use the implementation of
            :mod:`sage.graphs.base.static_dense_graph`

        - ``certificate`` -- boolean (default: ``False``); whether to return a
          triangle if one is found. This parameter is ignored when ``algorithm``
          is ``'matrix'``.

        EXAMPLES:

        The Petersen Graph is triangle-free::

            sage: g = graphs.PetersenGraph()
            sage: g.is_triangle_free()
            True

        or a complete Bipartite Graph::

            sage: G = graphs.CompleteBipartiteGraph(5,6)
            sage: G.is_triangle_free(algorithm='matrix')
            True
            sage: G.is_triangle_free(algorithm='bitset')
            True
            sage: G.is_triangle_free(algorithm='dense_graph')
            True

        a tripartite graph, though, contains many triangles::

            sage: G = (3 * graphs.CompleteGraph(5)).complement()
            sage: G.is_triangle_free(algorithm='matrix')
            False
            sage: G.is_triangle_free(algorithm='bitset')
            False
            sage: G.is_triangle_free(algorithm='dense_graph')
            False

        Asking for a certificate::

            sage: K4 = graphs.CompleteGraph(4)
            sage: K4.is_triangle_free(algorithm='dense_graph', certificate=True)
            (False, [0, 1, 2])
            sage: K4.is_triangle_free(algorithm='bitset', certificate=True)
            (False, [0, 1, 2])

        TESTS:

        Comparison of algorithms::

            sage: for i in range(10): # long time
            ....:     G = graphs.RandomBarabasiAlbert(50,2)
            ....:     bm = G.is_triangle_free(algorithm='matrix')
            ....:     bb = G.is_triangle_free(algorithm='bitset')
            ....:     bd = G.is_triangle_free(algorithm='dense_graph')
            ....:     if bm != bb or bm != bd:
            ....:        print("That's not good!")

        Asking for an unknown algorithm::

            sage: g.is_triangle_free(algorithm='tip top')
            Traceback (most recent call last):
            ...
            ValueError: Algorithm 'tip top' not yet implemented. Please contribute.

        Check the empty graph::

            sage: graphs.EmptyGraph().is_triangle_free()
            True
        """
        if algorithm == 'dense_graph':
            from sage.graphs.base.static_dense_graph import is_triangle_free
            return is_triangle_free(self, certificate=certificate)

        if algorithm == 'bitset':
            if self.order() < 3:
                return (True, []) if certificate else True
            from sage.data_structures.bitset import Bitset
            N = self.order()
            vertex_to_int = {}
            B = {}
            for i, u in enumerate(self):
                vertex_to_int[u] = i
                B[u] = Bitset(capacity=N)
            # map adjacency to bitsets
            for u, v in self.edge_iterator(labels=None):
                if u != v:
                    B[u].add(vertex_to_int[v])
                    B[v].add(vertex_to_int[u])
            # Search for a triangle
            for u, v in self.edge_iterator(labels=None):
                BB = B[u] & B[v]
                if BB:
                    if certificate:
                        for w in self.neighbor_iterator(u):
                            if vertex_to_int[w] in BB:
                                return False, [u, v, w]
                    return False
            return (True, []) if certificate else True

        elif algorithm == 'matrix':
            if self.order() < 3:
                return True
            return (self.adjacency_matrix()**3).trace() == 0

        else:
            raise ValueError("Algorithm '%s' not yet implemented. Please contribute." %(algorithm))

    @doc_index("Graph properties")
    def is_split(self):
        r"""
        Returns ``True`` if the graph is a Split graph, ``False`` otherwise.

        A Graph `G` is said to be a split graph if its vertices `V(G)` can be
        partitioned into two sets `K` and `I` such that the vertices of `K`
        induce a complete graph, and those of `I` are an independent set.

        There is a simple test to check whether a graph is a split graph (see,
        for instance, the book "Graph Classes, a survey" [BLS1999]_ page
        203) :

        Given the degree sequence `d_1 \geq ... \geq d_n` of `G`, a graph is a
        split graph if and only if :

        .. MATH::

            \sum_{i=1}^\omega d_i = \omega (\omega - 1) + \sum_{i=\omega + 1}^nd_i

        where `\omega = max \{i:d_i\geq i-1\}`.


        EXAMPLES:

        Split graphs are, in particular, chordal graphs. Hence, The Petersen
        graph can not be split::

            sage: graphs.PetersenGraph().is_split()
            False

        We can easily build some "random" split graph by creating a complete
        graph, and adding vertices only connected to some random vertices of the
        clique::

            sage: g = graphs.CompleteGraph(10)
            sage: sets = Subsets(Set(range(10)))
            sage: for i in range(10, 25):
            ....:    g.add_edges([(i,k) for k in sets.random_element()])
            sage: g.is_split()
            True

        Another characterisation of split graph states that a graph is a split
        graph if and only if does not contain the 4-cycle, 5-cycle or `2K_2` as
        an induced subgraph. Hence for the above graph we have::

            sage: forbidden_subgraphs = [graphs.CycleGraph(4), graphs.CycleGraph(5), 2 * graphs.CompleteGraph(2)]
            sage: sum(g.subgraph_search_count(H,induced=True) for H in forbidden_subgraphs)
            0
        """
        self._scream_if_not_simple()
        # our degree sequence is numbered from 0 to n-1, so to avoid
        # any mistake, let's fix it :-)
        degree_sequence = [0] + sorted(self.degree(), reverse=True)

        for i, d in enumerate(degree_sequence):
            if d >= i - 1:
                omega = i
            else:
                break

        left = sum(degree_sequence[:omega + 1])
        right = omega * (omega - 1) + sum(degree_sequence[omega + 1:])

        return left == right

    @doc_index("Algorithmically hard stuff")
    def is_perfect(self, certificate=False):
        r"""
        Tests whether the graph is perfect.

        A graph `G` is said to be perfect if `\chi(H)=\omega(H)` hold for any
        induced subgraph `H\subseteq_i G` (and so for `G` itself, too), where
        `\chi(H)` represents the chromatic number of `H`, and `\omega(H)` its
        clique number. The Strong Perfect Graph Theorem [CRST2006]_ gives
        another characterization of perfect graphs:

        A graph is perfect if and only if it contains no odd hole (cycle on an
        odd number `k` of vertices, `k>3`) nor any odd antihole (complement of a
        hole) as an induced subgraph.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); whether to return a
          certificate.

        OUTPUT:

        When ``certificate = False``, this function returns a boolean
        value. When ``certificate = True``, it returns a subgraph of ``self``
        isomorphic to an odd hole or an odd antihole if any, and ``None``
        otherwise.

        EXAMPLES:

        A Bipartite Graph is always perfect ::

            sage: g = graphs.RandomBipartite(8,4,.5)
            sage: g.is_perfect()
            True

        So is the line graph of a bipartite graph::

            sage: g = graphs.RandomBipartite(4,3,0.7)
            sage: g.line_graph().is_perfect() # long time
            True

        As well as the Cartesian product of two complete graphs::

            sage: g = graphs.CompleteGraph(3).cartesian_product(graphs.CompleteGraph(3))
            sage: g.is_perfect()
            True

        Interval Graphs, which are chordal graphs, too ::

            sage: g =  graphs.RandomIntervalGraph(7)
            sage: g.is_perfect()
            True

        The PetersenGraph, which is triangle-free and has chromatic number 3 is
        obviously not perfect::

            sage: g = graphs.PetersenGraph()
            sage: g.is_perfect()
            False

        We can obtain an induced 5-cycle as a certificate::

            sage: g.is_perfect(certificate=True)
            Subgraph of (Petersen graph): Graph on 5 vertices

        TESTS:

        Check that :trac:`13546` has been fixed::

            sage: Graph(':FgGE@I@GxGs', loops=False, multiedges=False).is_perfect()
            False
            sage: g = Graph({0: [2, 3, 4, 5],
            ....:            1: [3, 4, 5, 6],
            ....:            2: [0, 4, 5, 6],
            ....:            3: [0, 1, 5, 6],
            ....:            4: [0, 1, 2, 6],
            ....:            5: [0, 1, 2, 3],
            ....:            6: [1, 2, 3, 4]})
            sage: g.is_perfect()
            False

        TESTS::

            sage: Graph(':Ab').is_perfect()
            Traceback (most recent call last):
            ...
            ValueError: This method is only defined for simple graphs, and yours is not one of them !
            sage: g = Graph()
            sage: g.allow_loops(True)
            sage: g.add_edge(0,0)
            sage: g.edges()
            [(0, 0, None)]
            sage: g.is_perfect()
            Traceback (most recent call last):
            ...
            ValueError: This method is only defined for simple graphs, and yours is not one of them !

        """
        if self.has_multiple_edges() or self.has_loops():
            raise ValueError("This method is only defined for simple graphs,"
                             " and yours is not one of them !")
        if self.is_bipartite():
            return True if not certificate else None

        self_complement = self.complement()
        self_complement.remove_loops()
        self_complement.remove_multiple_edges()

        if self_complement.is_bipartite():
            return True if not certificate else None

        answer = self.is_odd_hole_free(certificate=certificate)
        if not (answer is True):
            return answer

        return self_complement.is_odd_hole_free(certificate=certificate)


    @doc_index("Graph properties")
    def is_edge_transitive(self):
        r"""
        Check if self is an edge transitive graph.

        A graph is edge-transitive if its automorphism group acts transitively
        on its edge set.

        Equivalently, if there exists for any pair of edges `uv,u'v'\in E(G)` an
        automorphism `\phi` of `G` such that `\phi(uv)=u'v'` (note this does not
        necessarily mean that `\phi(u)=u'` and `\phi(v)=v'`).

        .. SEEALSO::

          - :wikipedia:`Edge-transitive_graph`
          - :meth:`~Graph.is_arc_transitive`
          - :meth:`~Graph.is_half_transitive`
          - :meth:`~Graph.is_semi_symmetric`

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.is_edge_transitive()
            True
            sage: C = graphs.CubeGraph(3)
            sage: C.is_edge_transitive()
            True
            sage: G = graphs.GrayGraph()
            sage: G.is_edge_transitive()
            True
            sage: P = graphs.PathGraph(4)
            sage: P.is_edge_transitive()
            False
        """
        from sage.libs.gap.libgap import libgap

        if not self.size():
            return True

        A = self.automorphism_group()
        e = next(self.edge_iterator(labels=False))
        e = [A._domain_to_gap[e[0]], A._domain_to_gap[e[1]]]
        e.sort()
        return libgap(A).OrbitLength(e, libgap.OnSets) == self.size()

    @doc_index("Graph properties")
    def is_arc_transitive(self):
        r"""
        Check if self is an arc-transitive graph

        A graph is arc-transitive if its automorphism group acts transitively on
        its pairs of adjacent vertices.

        Equivalently, if there exists for any pair of edges `uv,u'v'\in E(G)` an
        automorphism `\phi_1` of `G` such that `\phi_1(u)=u'` and
        `\phi_1(v)=v'`, as well as another automorphism `\phi_2` of `G` such
        that `\phi_2(u)=v'` and `\phi_2(v)=u'`

        .. SEEALSO::

          - :wikipedia:`arc-transitive_graph`
          - :meth:`~Graph.is_edge_transitive`
          - :meth:`~Graph.is_half_transitive`
          - :meth:`~Graph.is_semi_symmetric`

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.is_arc_transitive()
            True
            sage: G = graphs.GrayGraph()
            sage: G.is_arc_transitive()
            False
        """
        from sage.libs.gap.libgap import libgap

        if not self.size():
            return True

        A = self.automorphism_group()
        e = next(self.edge_iterator(labels=False))
        e = [A._domain_to_gap[e[0]], A._domain_to_gap[e[1]]]

        return libgap(A).OrbitLength(e,libgap.OnTuples) == 2*self.size()

    @doc_index("Graph properties")
    def is_half_transitive(self):
        """
        Check if self is a half-transitive graph.

        A graph is half-transitive if it is both vertex and edge transitive
        but not arc-transitive.

        .. SEEALSO::

          - :wikipedia:`half-transitive_graph`
          - :meth:`~Graph.is_edge_transitive`
          - :meth:`~Graph.is_arc_transitive`
          - :meth:`~Graph.is_semi_symmetric`

        EXAMPLES:

        The Petersen Graph is not half-transitive::

            sage: P = graphs.PetersenGraph()
            sage: P.is_half_transitive()
            False

        The smallest half-transitive graph is the Holt Graph::

            sage: H = graphs.HoltGraph()
            sage: H.is_half_transitive()
            True
        """
        # A half-transitive graph always has only vertices of even degree
        if any(d % 2 for d in self.degree_iterator()):
            return False

        return (self.is_edge_transitive() and
                self.is_vertex_transitive() and
                not self.is_arc_transitive())

    @doc_index("Graph properties")
    def is_semi_symmetric(self):
        """
        Check if self is semi-symmetric.

        A graph is semi-symmetric if it is regular, edge-transitive but not
        vertex-transitive.

        .. SEEALSO::

          - :wikipedia:`Semi-symmetric_graph`
          - :meth:`~Graph.is_edge_transitive`
          - :meth:`~Graph.is_arc_transitive`
          - :meth:`~Graph.is_half_transitive`

        EXAMPLES:

        The Petersen graph is not semi-symmetric::

            sage: P = graphs.PetersenGraph()
            sage: P.is_semi_symmetric()
            False

        The Gray graph is the smallest possible cubic semi-symmetric graph::

            sage: G = graphs.GrayGraph()
            sage: G.is_semi_symmetric()
            True

        Another well known semi-symmetric graph is the Ljubljana graph::

            sage: L = graphs.LjubljanaGraph()
            sage: L.is_semi_symmetric()
            True
        """
        # A semi-symmetric graph is always bipartite
        if not self.is_bipartite():
            return False

        return (self.is_regular() and
                self.is_edge_transitive() and not
                self.is_vertex_transitive())

    @doc_index("Graph properties")
    def is_path(self):
        r"""
        Check whether ``self`` is a path.

        A connected graph of order `n \geq 2` is a path if it is a tree
        (see :meth:`is_tree`) with `n-2` vertices of degree 2 and two of
        degree 1. By convention, a graph of order 1 without loops is a path,
        but the empty graph is not a path.

        EXAMPLES:

            sage: G = graphs.PathGraph(5)
            sage: G.is_path()
            True
            sage: H = graphs.CycleGraph(5)
            sage: H.is_path()
            False
            sage: D = graphs.PathGraph(5).disjoint_union(graphs.CycleGraph(5))
            sage: D.is_path()
            False
            sage: E = graphs.EmptyGraph()
            sage: E.is_path()
            False
            sage: O = Graph([[1], []])
            sage: O.is_path()
            True
            sage: O.allow_loops(True)
            sage: O.add_edge(1, 1)
            sage: O.is_path()
            False
        """
        order = self.order()
        if order != self.size() + 1:
            return False

        if order <= 1:
            return order == 1

        deg_one_counter = 0
        seen_counter = 0
        for v in self.depth_first_search(next(self.vertex_iterator())):
            seen_counter += 1
            deg = self._backend.degree(v, False)
            if deg == 1:
                deg_one_counter += 1
                if deg_one_counter > 2:
                    return False

            elif deg != 2:
                return False
        return deg_one_counter == 2 and seen_counter == order


    @doc_index("Connectivity, orientations, trees")
    def degree_constrained_subgraph(self, bounds, solver=None, verbose=0,
                                    *, integrality_tolerance=1e-3):
        r"""
        Returns a degree-constrained subgraph.

        Given a graph `G` and two functions `f, g:V(G)\rightarrow \mathbb Z`
        such that `f \leq g`, a degree-constrained subgraph in `G` is
        a subgraph `G' \subseteq G` such that for any vertex `v \in G`,
        `f(v) \leq d_{G'}(v) \leq g(v)`.

        INPUT:

        - ``bounds`` -- (default: ``None``); Two possibilities:

          - A dictionary whose keys are the vertices, and values a pair of
            real values ``(min,max)`` corresponding to the values
            `(f(v),g(v))`.

          - A function associating to each vertex a pair of
            real values ``(min,max)`` corresponding to the values
            `(f(v),g(v))`.

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        OUTPUT:

        - When a solution exists, this method outputs the degree-constrained
          subgraph as a Graph object.

        - When no solution exists, returns ``False``.

        .. NOTE::

            - This algorithm computes the degree-constrained subgraph of minimum
              weight.
            - If the graph's edges are weighted, these are taken into account.
            - This problem can be solved in polynomial time.

        EXAMPLES:

        Is there a perfect matching in an even cycle? ::

            sage: g = graphs.CycleGraph(6)
            sage: bounds = lambda x: [1,1]
            sage: m = g.degree_constrained_subgraph(bounds=bounds)
            sage: m.size()
            3
        """
        self._scream_if_not_simple()
        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

        p = MixedIntegerLinearProgram(maximization=False, solver=solver)
        b = p.new_variable(binary=True)

        if isinstance(bounds,dict):
            f_bounds = lambda x: bounds[x]
        else:
            f_bounds = bounds


        if self.weighted():
            from sage.rings.real_mpfr import RR
            weight = lambda x: x if x in RR else 1
        else:
            weight = lambda x: 1

        for v in self:
            minimum,maximum = f_bounds(v)
            p.add_constraint(p.sum(b[frozenset((x,y))]*weight(l) for x,y,l in self.edges_incident(v)),
                                 min=minimum, max=maximum)

        p.set_objective(p.sum(b[frozenset((x,y))]*weight(l) for x,y,l in self.edge_iterator()))

        try:
            p.solve(log=verbose)
        except MIPSolverException:
            return False

        g = copy(self)
        b = p.get_values(b, convert=bool, tolerance=integrality_tolerance)
        g.delete_edges(e for e in g.edge_iterator(labels=False) if not b[frozenset(e)])
        return g

    ### Orientations

    @doc_index("Connectivity, orientations, trees")
    def strong_orientation(self):
        r"""
        Returns a strongly connected orientation of the current graph.

        An orientation of an undirected graph is a digraph obtained by giving an
        unique direction to each of its edges. An orientation is said to be
        strong if there is a directed path between each pair of vertices.  See
        also the :wikipedia:`Strongly_connected_component`.

        If the graph is 2-edge-connected, a strongly connected orientation
        can be found in linear time. If the given graph is not 2-connected,
        the orientation returned will ensure that each 2-connected component
        has a strongly connected orientation.

        OUTPUT:

        A digraph representing an orientation of the current graph.

        .. NOTE::

            - This method assumes the graph is connected.
            - This algorithm works in O(m).

        EXAMPLES:

        For a 2-regular graph, a strong orientation gives to each vertex an
        out-degree equal to 1::

            sage: g = graphs.CycleGraph(5)
            sage: g.strong_orientation().out_degree()
            [1, 1, 1, 1, 1]

        The Petersen Graph is 2-edge connected. It then has a strongly connected
        orientation::

            sage: g = graphs.PetersenGraph()
            sage: o = g.strong_orientation()
            sage: len(o.strongly_connected_components())
            1

        The same goes for the CubeGraph in any dimension ::

            sage: all(len(graphs.CubeGraph(i).strong_orientation().strongly_connected_components()) == 1 for i in range(2,6))
            True

        A multigraph also has a strong orientation ::

            sage: g = Graph([(1,2),(1,2)], multiedges=True)
            sage: g.strong_orientation()
            Multi-digraph on 2 vertices

        """
        from sage.graphs.digraph import DiGraph
        d = DiGraph(multiedges=self.allows_multiple_edges())
        i = 0

        # The algorithm works through a depth-first search. Any edge
        # used in the depth-first search is oriented in the direction
        # in which it has been used. All the other edges are oriented
        # backward

        v = next(self.vertex_iterator())
        seen = {}
        i = 1

        # Time at which the vertices have been discovered
        seen[v] = i

        # indicates the stack of edges to explore
        next_ = self.edges_incident(v)

        while next_:
            e = next_.pop()

            # Ignore loops
            if e[0] == e[1]:
                continue

            # We assume e[0] to be a `seen` vertex
            e = e if seen.get(e[0], False) is not False else (e[1], e[0], e[2])

            # If we discovered a new vertex
            if seen.get(e[1], False) is False:
                d.add_edge(e)
                next_.extend(ee for ee in self.edges_incident(e[1])
                                 if ((e[0],e[1]) != (ee[0],ee[1])) and ((e[0],e[1]) != (ee[1],ee[0])))
                i += 1
                seen[e[1]] = i

            # Else, we orient the edges backward
            else:
                if seen[e[0]] < seen[e[1]]:
                    d.add_edge(e[1], e[0], e[2])
                else:
                    d.add_edge(e)

        # Case of multiple edges. If another edge has already been inserted, we
        # add the new one in the opposite direction.
        tmp = None
        for e in self.multiple_edges():
            if tmp == (e[0], e[1]):
                if d.has_edge(e[0], e[1]):
                    d.add_edge(e[1], e[0], e[2])
                else:
                    d.add_edge(e)
            tmp = (e[0], e[1])

        return d

    @doc_index("Connectivity, orientations, trees")
    def minimum_outdegree_orientation(self, use_edge_labels=False, solver=None, verbose=0,
                                      *, integrality_tolerance=1e-3):
        r"""
        Returns an orientation of ``self`` with the smallest possible maximum
        outdegree.

        Given a Graph `G`, it is polynomial to compute an orientation `D` of the
        edges of `G` such that the maximum out-degree in `D` is minimized. This
        problem, though, is NP-complete in the weighted case [AMOZ2006]_.

        INPUT:

        - ``use_edge_labels`` -- boolean (default: ``False``)

          - When set to ``True``, uses edge labels as weights to compute the
            orientation and assumes a weight of `1` when there is no value
            available for a given edge.

          - When set to ``False`` (default), gives a weight of 1 to all the
            edges.

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        EXAMPLES:

        Given a complete bipartite graph `K_{n,m}`, the maximum out-degree of an
        optimal orientation is `\left\lceil \frac {nm} {n+m}\right\rceil`::

            sage: g = graphs.CompleteBipartiteGraph(3,4)
            sage: o = g.minimum_outdegree_orientation()
            sage: max(o.out_degree()) == integer_ceil((4*3)/(3+4))
            True
        """
        self._scream_if_not_simple()
        if self.is_directed():
            raise ValueError("Cannot compute an orientation of a DiGraph. "+\
                                 "Please convert it to a Graph if you really mean it.")

        if use_edge_labels:
            from sage.rings.real_mpfr import RR
            def weight(e):
                l = self.edge_label(e)
                return l if l in RR else 1
        else:
            def weight(e):
                return 1

        from sage.numerical.mip import MixedIntegerLinearProgram

        p = MixedIntegerLinearProgram(maximization=False, solver=solver)
        degree = p.new_variable(nonnegative=True)

        # The orientation of an edge is boolean and indicates whether the edge
        # uv goes from u to v ( equal to 0 ) or from v to u ( equal to 1)
        orientation = p.new_variable(binary=True)

        # Whether an edge adjacent to a vertex u counts positively or
        # negatively. To do so, we first fix an arbitrary extremity per edge uv.
        ext = {frozenset(e): e[0] for e in self.edge_iterator(labels=False)}
        def outgoing(u, e, variable):
            if u == ext[frozenset(e)]:
                return variable
            else:
                return 1 - variable

        for u in self:
            p.add_constraint(p.sum(weight(e) * outgoing(u, e, orientation[frozenset(e)])
                                       for e in self.edge_iterator(vertices=[u], labels=False))
                                 - degree['max'], max=0)

        p.set_objective(degree['max'])

        p.solve(log=verbose)

        orientation = p.get_values(orientation, convert=bool, tolerance=integrality_tolerance)

        # All the edges from self are doubled in O
        # ( one in each direction )
        from sage.graphs.digraph import DiGraph
        O = DiGraph(self)

        # Builds the list of edges that should be removed
        edges = []

        for e in self.edge_iterator(labels=None):
            if orientation[frozenset(e)]:
                edges.append(e[::-1])
            else:
                edges.append(e)

        O.delete_edges(edges)

        return O

    @doc_index("Connectivity, orientations, trees")
    def bounded_outdegree_orientation(self, bound, solver=None, verbose=False,
                                      *, integrality_tolerance=1e-3):
        r"""
        Computes an orientation of ``self`` such that every vertex `v` has
        out-degree less than `b(v)`

        INPUT:

        - ``bound`` -- Maximum bound on the out-degree. Can be of three
          different types :

         * An integer `k`. In this case, computes an orientation whose maximum
           out-degree is less than `k`.

         * A dictionary associating to each vertex its associated maximum
           out-degree.

         * A function associating to each vertex its associated maximum
           out-degree.

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        OUTPUT:

        A DiGraph representing the orientation if it exists. A ``ValueError``
        exception is raised otherwise.

        ALGORITHM:

        The problem is solved through a maximum flow :

        Given a graph `G`, we create a ``DiGraph`` `D` defined on `E(G)\cup
        V(G)\cup \{s,t\}`. We then link `s` to all of `V(G)` (these edges having
        a capacity equal to the bound associated to each element of `V(G)`), and
        all the elements of `E(G)` to `t` . We then link each `v \in V(G)` to
        each of its incident edges in `G`. A maximum integer flow of value
        `|E(G)|` corresponds to an admissible orientation of `G`. Otherwise,
        none exists.

        EXAMPLES:

        There is always an orientation of a graph `G` such that a vertex `v` has
        out-degree at most `\lceil \frac {d(v)} 2 \rceil`::

            sage: g = graphs.RandomGNP(40, .4)
            sage: b = lambda v: integer_ceil(g.degree(v)/2)
            sage: D = g.bounded_outdegree_orientation(b)
            sage: all( D.out_degree(v) <= b(v) for v in g )
            True


        Chvatal's graph, being 4-regular, can be oriented in such a way that its
        maximum out-degree is 2::

            sage: g = graphs.ChvatalGraph()
            sage: D = g.bounded_outdegree_orientation(2)
            sage: max(D.out_degree())
            2

        For any graph `G`, it is possible to compute an orientation such that
        the maximum out-degree is at most the maximum average degree of `G`
        divided by 2. Anything less, though, is impossible.

            sage: g = graphs.RandomGNP(40, .4)
            sage: mad = g.maximum_average_degree()

        Hence this is possible ::

            sage: d = g.bounded_outdegree_orientation(integer_ceil(mad/2))

        While this is not::

            sage: try:
            ....:     g.bounded_outdegree_orientation(integer_ceil(mad/2-1))
            ....:     print("Error")
            ....: except ValueError:
            ....:     pass

        TESTS:

        As previously for random graphs, but more intensively::

            sage: for i in range(30):      # long time (up to 6s on sage.math, 2012)
            ....:     g = graphs.RandomGNP(40, .4)
            ....:     b = lambda v: integer_ceil(g.degree(v)/2)
            ....:     D = g.bounded_outdegree_orientation(b)
            ....:     if not (
            ....:          all( D.out_degree(v) <= b(v) for v in g ) or
            ....:          D.size() != g.size()):
            ....:         print("Something wrong happened")

        """
        self._scream_if_not_simple()
        from sage.graphs.digraph import DiGraph
        n = self.order()

        if not n:
            return DiGraph()

        vertices = list(self)
        vertices_id = {y: x for x,y in enumerate(vertices)}

        b = {}

        # Checking the input type. We make a dictionary out of it
        if isinstance(bound, dict):
            b = bound
        else:
            try:
                b = dict(zip(vertices,map(bound, vertices)))

            except TypeError:
                b = dict(zip(vertices, [bound]*n))

        d = DiGraph()

        # Adding the edges (s,v) and ((u,v),t)
        d.add_edges(('s', vertices_id[v], b[v]) for v in vertices)

        d.add_edges(((vertices_id[u], vertices_id[v]), 't', 1)
                     for u,v in self.edges(labels=None) )

        # each v is linked to its incident edges

        for u,v in self.edge_iterator(labels=None):
            u,v = vertices_id[u], vertices_id[v]
            d.add_edge(u, (u,v), 1)
            d.add_edge(v, (u,v), 1)

        # Solving the maximum flow
        value, flow = d.flow('s','t', value_only=False, integer=True,
                             use_edge_labels=True, solver=solver, verbose=verbose,
                             integrality_tolerance=integrality_tolerance)

        if value != self.size():
            raise ValueError("No orientation exists for the given bound")

        D = DiGraph()
        D.add_vertices(vertices)

        # The flow graph may not contain all the vertices, if they are
        # not part of the flow...

        for u in [x for x in range(n) if x in flow]:

            for uu,vv in flow.neighbors_out(u):
                v = vv if vv != u else uu
                D.add_edge(vertices[u], vertices[v])

        # I do not like when a method destroys the embedding ;-)
        D.set_pos(self.get_pos())

        return D

    @doc_index("Connectivity, orientations, trees")
    def orientations(self, data_structure=None, sparse=None):
        r"""
        Return an iterator over orientations of ``self``.

        An *orientation* of an undirected graph is a directed graph such that
        every edge is assigned a direction.  Hence there are `2^s` oriented
        digraphs for a simple graph with `s` edges.

        INPUT:

        - ``data_structure`` -- one of ``"sparse"``, ``"static_sparse"``, or
          ``"dense"``; see the documentation of :class:`Graph` or
          :class:`DiGraph`; default is the data structure of ``self``

        - ``sparse`` -- boolean (default: ``None``); ``sparse=True`` is an alias
          for ``data_structure="sparse"``, and ``sparse=False`` is an alias for
          ``data_structure="dense"``. By default (``None``), guess the most
          suitable data structure.

        .. WARNING::

            This always considers multiple edges of graphs as distinguishable,
            and hence, may have repeated digraphs.

        EXAMPLES::

            sage: G = Graph([[1,2,3], [(1, 2, 'a'), (1, 3, 'b')]], format='vertices_and_edges')
            sage: it = G.orientations()
            sage: D = next(it)
            sage: D.edges()
            [(1, 2, 'a'), (1, 3, 'b')]
            sage: D = next(it)
            sage: D.edges()
            [(1, 2, 'a'), (3, 1, 'b')]

        TESTS::

            sage: G = Graph()
            sage: D = [g for g in G.orientations()]
            sage: len(D)
            1
            sage: D[0]
            Digraph on 0 vertices

            sage: G = Graph(5)
            sage: it = G.orientations()
            sage: D = next(it)
            sage: D.size()
            0

            sage: G = Graph([[1,2,'a'], [1,2,'b']], multiedges=True)
            sage: len(list(G.orientations()))
            4

            sage: G = Graph([[1,2], [1,1]], loops=True)
            sage: len(list(G.orientations()))
            2

            sage: G = Graph([[1,2],[2,3]])
            sage: next(G.orientations())
            Digraph on 3 vertices
            sage: G = graphs.PetersenGraph()
            sage: next(G.orientations())
            An orientation of Petersen graph: Digraph on 10 vertices

        An orientation must have the same ground set of vertices as the original
        graph (:trac:`24366`)::

            sage: G = Graph(1)
            sage: next(G.orientations())
            Digraph on 1 vertex
        """
        if sparse is not None:
            if data_structure is not None:
                raise ValueError("cannot specify both 'sparse' and 'data_structure'")
            data_structure = "sparse" if sparse else "dense"
        if data_structure is None:
            from sage.graphs.base.dense_graph import DenseGraphBackend
            from sage.graphs.base.sparse_graph import SparseGraphBackend
            if isinstance(self._backend, DenseGraphBackend):
                data_structure = "dense"
            elif isinstance(self._backend, SparseGraphBackend):
                data_structure = "sparse"
            else:
                data_structure = "static_sparse"

        name = self.name()
        if name:
            name = 'An orientation of ' + name

        if not self.size():
            D = DiGraph(data=[self.vertices(), []],
                        format='vertices_and_edges',
                        name=name,
                        pos=self._pos,
                        multiedges=self.allows_multiple_edges(),
                        loops=self.allows_loops(),
                        data_structure=data_structure)
            if hasattr(self, '_embedding'):
                D._embedding = copy(self._embedding)
            yield D
            return

        E = [[(u,v,label), (v,u,label)] if u != v else [(u,v,label)]
             for u,v,label in self.edge_iterator()]
        verts = self.vertices()
        for edges in itertools.product(*E):
            D = DiGraph(data=[verts, edges],
                        format='vertices_and_edges',
                        name=name,
                        pos=self._pos,
                        multiedges=self.allows_multiple_edges(),
                        loops=self.allows_loops(),
                        data_structure=data_structure)
            if hasattr(self, '_embedding'):
                D._embedding = copy(self._embedding)
            yield D

    ### Coloring

    @doc_index("Basic methods")
    def bipartite_color(self):
        """
        Return a dictionary with vertices as the keys and the color class
        as the values.

        Fails with an error if the graph is not bipartite.

        EXAMPLES::

            sage: graphs.CycleGraph(4).bipartite_color()
            {0: 1, 1: 0, 2: 1, 3: 0}
            sage: graphs.CycleGraph(5).bipartite_color()
            Traceback (most recent call last):
            ...
            RuntimeError: Graph is not bipartite.

        TESTS::

            sage: Graph().bipartite_color()
            {}
        """
        isit, certificate = self.is_bipartite(certificate=True)

        if isit:
            return certificate
        else:
            raise RuntimeError("Graph is not bipartite.")

    @doc_index("Basic methods")
    def bipartite_sets(self):
        r"""
        Return `(X,Y)` where `X` and `Y` are the nodes in each bipartite set of
        graph `G`.

        Fails with an error if graph is not bipartite.

        EXAMPLES::

            sage: graphs.CycleGraph(4).bipartite_sets()
            ({0, 2}, {1, 3})
            sage: graphs.CycleGraph(5).bipartite_sets()
            Traceback (most recent call last):
            ...
            RuntimeError: Graph is not bipartite.
        """
        color = self.bipartite_color()
        left = set()
        right = set()

        for u,s in color.items():
            if s:
                left.add(u)
            else:
                right.add(u)

        return left, right

    @doc_index("Coloring")
    def chromatic_index(self, solver=None, verbose=0, *, integrality_tolerance=1e-3):
        r"""
        Return the chromatic index of the graph.

        The chromatic index is the minimal number of colors needed to properly
        color the edges of the graph.

        INPUT:

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        This method is a frontend for method
        :meth:`sage.graphs.graph_coloring.edge_coloring` that uses a mixed
        integer-linear programming formulation to compute the chromatic index.

        .. SEEALSO::

            - :wikipedia:`Edge_coloring` for further details on edge coloring
            - :meth:`sage.graphs.graph_coloring.edge_coloring`
            - :meth:`~Graph.fractional_chromatic_index`
            - :meth:`~Graph.chromatic_number`

        EXAMPLES:

        The clique `K_n` has chromatic index `n` when `n` is odd and `n-1` when
        `n` is even::

            sage: graphs.CompleteGraph(4).chromatic_index()
            3
            sage: graphs.CompleteGraph(5).chromatic_index()
            5
            sage: graphs.CompleteGraph(6).chromatic_index()
            5

        The path `P_n` with `n \geq 2` has chromatic index 2::

            sage: graphs.PathGraph(5).chromatic_index()
            2

        The windmill graph with parameters `k,n` has chromatic index `(k-1)n`::

            sage: k,n = 3,4
            sage: G = graphs.WindmillGraph(k,n)
            sage: G.chromatic_index() == (k-1)*n
            True

        TESTS:

        Graphs without vertices or edges::

            sage: Graph().chromatic_index()
            0
            sage: Graph(2).chromatic_index()
            0
        """
        if not self.order() or not self.size():
            return 0

        from sage.graphs.graph_coloring import edge_coloring
        return edge_coloring(self, value_only=True, solver=solver, verbose=verbose,
                             integrality_tolerance=integrality_tolerance)

    @doc_index("Coloring")
    def chromatic_number(self, algorithm="DLX", solver=None, verbose=0,
                         *, integrality_tolerance=1e-3):
        r"""
        Return the minimal number of colors needed to color the vertices of the
        graph.

        INPUT:

        - ``algorithm`` -- Select an algorithm from the following supported
          algorithms:

          - If ``algorithm="DLX"`` (default), the chromatic number is computed
            using the dancing link algorithm. It is inefficient speedwise to
            compute the chromatic number through the dancing link algorithm
            because this algorithm computes *all* the possible colorings to
            check that one exists.

          - If ``algorithm="CP"``, the chromatic number is computed using the
            coefficients of the chromatic polynomial. Again, this method is
            inefficient in terms of speed and it only useful for small graphs.

          - If ``algorithm="MILP"``, the chromatic number is computed using a
            mixed integer linear program. The performance of this implementation
            is affected by whether optional MILP solvers have been installed
            (see the :mod:`MILP module <sage.numerical.mip>`, or Sage's tutorial
            on Linear Programming).

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        .. SEEALSO::

            For more functions related to graph coloring, see the module
            :mod:`sage.graphs.graph_coloring`.

        EXAMPLES::

            sage: G = Graph({0: [1, 2, 3], 1: [2]})
            sage: G.chromatic_number(algorithm="DLX")
            3
            sage: G.chromatic_number(algorithm="MILP")
            3
            sage: G.chromatic_number(algorithm="CP")
            3

        A bipartite graph has (by definition) chromatic number 2::

            sage: graphs.RandomBipartite(50,50,0.7).chromatic_number()
            2

        A complete multipartite graph with k parts has chromatic number `k`::

            sage: all(graphs.CompleteMultipartiteGraph([5]*i).chromatic_number() == i for i in range(2,5))
            True

        The complete graph has the largest chromatic number from all the graphs
        of order `n`. Namely its chromatic number is `n`::

            sage: all(graphs.CompleteGraph(i).chromatic_number() == i for i in range(10))
            True

        The Kneser graph with parameters `(n, 2)` for `n > 3` has chromatic
        number `n-2`::

            sage: all(graphs.KneserGraph(i,2).chromatic_number() == i-2 for i in range(4,6))
            True

        The Flower Snark graph has chromatic index 4 hence its line graph has
        chromatic number 4::

            sage: graphs.FlowerSnark().line_graph().chromatic_number()
            4

        TESTS::

            sage: G = Graph()
            sage: G.chromatic_number(algorithm="DLX")
            0
            sage: G.chromatic_number(algorithm="MILP")
            0
            sage: G.chromatic_number(algorithm="CP")
            0

            sage: G = Graph({0: [1, 2, 3], 1: [2]})
            sage: G.chromatic_number(algorithm="foo")
            Traceback (most recent call last):
            ...
            ValueError: The 'algorithm' keyword must be set to either 'DLX', 'MILP' or 'CP'.
        """
        self._scream_if_not_simple(allow_multiple_edges=True)
        # default built-in algorithm; bad performance
        if algorithm == "DLX":
            from sage.graphs.graph_coloring import chromatic_number
            return chromatic_number(self)
        # Algorithm with good performance, but requires an optional
        # package: choose any of GLPK or CBC.
        elif algorithm == "MILP":
            from sage.graphs.graph_coloring import vertex_coloring
            return vertex_coloring(self, value_only=True, solver=solver, verbose=verbose,
                                   integrality_tolerance=integrality_tolerance)
        # another algorithm with bad performance; only good for small graphs
        elif algorithm == "CP":
            f = self.chromatic_polynomial()
            i = 0
            while not f(i):
                i += 1
            return i
        else:
            raise ValueError("The 'algorithm' keyword must be set to either 'DLX', 'MILP' or 'CP'.")

    @doc_index("Coloring")
    def coloring(self, algorithm="DLX", hex_colors=False, solver=None, verbose=0,
                 *, integrality_tolerance=1e-3):
        r"""
        Return the first (optimal) proper vertex-coloring found.

        INPUT:

        - ``algorithm`` -- Select an algorithm from the following supported
          algorithms:

          - If ``algorithm="DLX"`` (default), the coloring is computed using the
            dancing link algorithm.

          - If ``algorithm="MILP"``, the coloring is computed using a mixed
            integer linear program. The performance of this implementation is
            affected by whether optional MILP solvers have been installed (see
            the :mod:`MILP module <sage.numerical.mip>`).

        - ``hex_colors`` -- boolean (default: ``False``); if ``True``, return a
          dictionary which can easily be used for plotting.

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        .. SEEALSO::

            For more functions related to graph coloring, see the
            module :mod:`sage.graphs.graph_coloring`.

        EXAMPLES::

            sage: G = Graph("Fooba")
            sage: P = G.coloring(algorithm="MILP")
            sage: Q = G.coloring(algorithm="DLX")
            sage: def are_equal_colorings(A, B):
            ....:     return Set(map(Set, A)) == Set(map(Set, B))
            sage: are_equal_colorings(P, [[1, 2, 3], [0, 5, 6], [4]])
            True
            sage: are_equal_colorings(P, Q)
            True
            sage: G.plot(partition=P)
            Graphics object consisting of 16 graphics primitives
            sage: G.coloring(hex_colors=True, algorithm="MILP")
            {'#0000ff': [4], '#00ff00': [0, 6, 5], '#ff0000': [2, 1, 3]}
            sage: H = G.coloring(hex_colors=True, algorithm="DLX")
            sage: H
            {'#0000ff': [4], '#00ff00': [1, 2, 3], '#ff0000': [0, 5, 6]}
            sage: G.plot(vertex_colors=H)
            Graphics object consisting of 16 graphics primitives

        .. PLOT::

            g = Graph("Fooba")
            sphinx_plot(g.plot(partition=g.coloring()))

        TESTS::

            sage: G.coloring(algorithm="foo")
            Traceback (most recent call last):
            ...
            ValueError: The 'algorithm' keyword must be set to either 'DLX' or 'MILP'.
        """
        self._scream_if_not_simple(allow_multiple_edges=True)
        if algorithm == "MILP":
            from sage.graphs.graph_coloring import vertex_coloring
            return vertex_coloring(self, hex_colors=hex_colors, solver=solver, verbose=verbose,
                                   integrality_tolerance=integrality_tolerance)
        elif algorithm == "DLX":
            from sage.graphs.graph_coloring import first_coloring
            return first_coloring(self, hex_colors=hex_colors)
        else:
            raise ValueError("The 'algorithm' keyword must be set to either 'DLX' or 'MILP'.")

    @doc_index("Coloring")
    def chromatic_symmetric_function(self, R=None):
        r"""
        Return the chromatic symmetric function of ``self``.

        Let `G` be a graph. The chromatic symmetric function `X_G` was described
        in [Sta1995]_, specifically Theorem 2.5 states that

        .. MATH::

            X_G = \sum_{F \subseteq E(G)} (-1)^{|F|} p_{\lambda(F)},

        where `\lambda(F)` is the partition of the sizes of the connected
        components of the subgraph induced by the edges `F` and `p_{\mu}` is the
        powersum symmetric function.

        INPUT:

        - ``R`` -- (optional) the base ring for the symmetric functions;
          this uses `\ZZ` by default

        EXAMPLES::

            sage: s = SymmetricFunctions(ZZ).s()
            sage: G = graphs.CycleGraph(5)
            sage: XG = G.chromatic_symmetric_function(); XG
            p[1, 1, 1, 1, 1] - 5*p[2, 1, 1, 1] + 5*p[2, 2, 1]
             + 5*p[3, 1, 1] - 5*p[3, 2] - 5*p[4, 1] + 4*p[5]
            sage: s(XG)
            30*s[1, 1, 1, 1, 1] + 10*s[2, 1, 1, 1] + 10*s[2, 2, 1]

        Not all graphs have a positive Schur expansion::

            sage: G = graphs.ClawGraph()
            sage: XG = G.chromatic_symmetric_function(); XG
            p[1, 1, 1, 1] - 3*p[2, 1, 1] + 3*p[3, 1] - p[4]
            sage: s(XG)
            8*s[1, 1, 1, 1] + 5*s[2, 1, 1] - s[2, 2] + s[3, 1]

        We show that given a triangle `\{e_1, e_2, e_3\}`, we have
        `X_G = X_{G - e_1} + X_{G - e_2} - X_{G - e_1 - e_2}`::

            sage: G = Graph([[1,2],[1,3],[2,3]])
            sage: XG = G.chromatic_symmetric_function()
            sage: G1 = copy(G)
            sage: G1.delete_edge([1,2])
            sage: XG1 = G1.chromatic_symmetric_function()
            sage: G2 = copy(G)
            sage: G2.delete_edge([1,3])
            sage: XG2 = G2.chromatic_symmetric_function()
            sage: G3 = copy(G1)
            sage: G3.delete_edge([1,3])
            sage: XG3 = G3.chromatic_symmetric_function()
            sage: XG == XG1 + XG2 - XG3
            True
        """
        from sage.combinat.sf.sf import SymmetricFunctions
        from sage.combinat.partition import _Partitions
        from sage.misc.misc import powerset
        if R is None:
            R = ZZ
        p = SymmetricFunctions(R).p()
        ret = p.zero()
        for F in powerset(self.edges()):
            la = _Partitions(self.subgraph(edges=F).connected_components_sizes())
            ret += (-1)**len(F) * p[la]
        return ret

    @doc_index("Coloring")
    def chromatic_quasisymmetric_function(self, t=None, R=None):
        r"""
        Return the chromatic quasisymmetric function of ``self``.

        Let `G` be a graph whose vertex set is totally ordered. The chromatic
        quasisymmetric function `X_G(t)` was first described in [SW2012]_. We
        use the equivalent definition given in [BC2018]_:

        .. MATH::

            X_G(t) = \sum_{\sigma=(\sigma_1,\ldots,\sigma_n)}
            t^{\operatorname{asc}(\sigma)}
            M_{|\sigma_1|,\ldots,|\sigma_n|},

        where we sum over all ordered set partitions of the vertex set of `G`
        such that each block `\sigma_i` is an independent (i.e., stable) set of
        `G`, and where `\operatorname{asc}(\sigma)` denotes the number of edges
        `\{u, v\}` of `G` such that `u < v` and `v` appears in a later part of
        `\sigma` than `u`.

        INPUT:

        - ``t`` -- (optional) the parameter `t`; uses the variable `t` in
          `\ZZ[t]` by default

        - ``R`` -- (optional) the base ring for the quasisymmetric functions;
          uses the parent of `t` by default

        EXAMPLES::

            sage: G = Graph([[1,2,3], [[1,3], [2,3]]])
            sage: G.chromatic_quasisymmetric_function()
            (2*t^2+2*t+2)*M[1, 1, 1] + M[1, 2] + t^2*M[2, 1]
            sage: G = graphs.PathGraph(4)
            sage: XG = G.chromatic_quasisymmetric_function(); XG
            (t^3+11*t^2+11*t+1)*M[1, 1, 1, 1] + (3*t^2+3*t)*M[1, 1, 2]
             + (3*t^2+3*t)*M[1, 2, 1] + (3*t^2+3*t)*M[2, 1, 1]
             + (t^2+t)*M[2, 2]
            sage: XG.to_symmetric_function()
            (t^3+11*t^2+11*t+1)*m[1, 1, 1, 1] + (3*t^2+3*t)*m[2, 1, 1]
             + (t^2+t)*m[2, 2]
            sage: G = graphs.CompleteGraph(4)
            sage: G.chromatic_quasisymmetric_function()
            (t^6+3*t^5+5*t^4+6*t^3+5*t^2+3*t+1)*M[1, 1, 1, 1]

        Not all chromatic quasisymmetric functions are symmetric::

            sage: G = Graph([[1,2], [1,5], [3,4], [3,5]])
            sage: G.chromatic_quasisymmetric_function().is_symmetric()
            False

        We check that at `t = 1`, we recover the usual chromatic symmetric
        function::

            sage: p = SymmetricFunctions(QQ).p()
            sage: G = graphs.CycleGraph(5)
            sage: XG = G.chromatic_quasisymmetric_function(t=1); XG
            120*M[1, 1, 1, 1, 1] + 30*M[1, 1, 1, 2] + 30*M[1, 1, 2, 1]
             + 30*M[1, 2, 1, 1] + 10*M[1, 2, 2] + 30*M[2, 1, 1, 1]
             + 10*M[2, 1, 2] + 10*M[2, 2, 1]
            sage: p(XG.to_symmetric_function())
            p[1, 1, 1, 1, 1] - 5*p[2, 1, 1, 1] + 5*p[2, 2, 1]
             + 5*p[3, 1, 1] - 5*p[3, 2] - 5*p[4, 1] + 4*p[5]

            sage: G = graphs.ClawGraph()
            sage: XG = G.chromatic_quasisymmetric_function(t=1); XG
            24*M[1, 1, 1, 1] + 6*M[1, 1, 2] + 6*M[1, 2, 1] + M[1, 3]
             + 6*M[2, 1, 1] + M[3, 1]
            sage: p(XG.to_symmetric_function())
            p[1, 1, 1, 1] - 3*p[2, 1, 1] + 3*p[3, 1] - p[4]
        """
        from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
        from sage.combinat.set_partition_ordered import OrderedSetPartitions
        if t is None:
            t = ZZ['t'].gen()
        if R is None:
            R = t.parent()
        M = QuasiSymmetricFunctions(R).M()
        ret = M.zero()
        V = self.vertices()

        def asc(sigma):
            stat = 0
            for i, s in enumerate(sigma):
                for u in s:
                    stat += sum(1 for p in sigma[i+1:] for v in p
                                if v > u and self.has_edge(u, v))
            return stat

        for sigma in OrderedSetPartitions(V):
            if any(not self.is_independent_set(s) for s in sigma):
                continue
            ret += M.term(sigma.to_composition(), t**asc(sigma))
        return ret

    @doc_index("Leftovers")
    def matching(self, value_only=False, algorithm="Edmonds",
                 use_edge_labels=False, solver=None, verbose=0,
                 *, integrality_tolerance=1e-3):
        r"""
        Return a maximum weighted matching of the graph represented by the list
        of its edges.

        For more information, see the :wikipedia:`Matching_(graph_theory)`.

        Given a graph `G` such that each edge `e` has a weight `w_e`, a maximum
        matching is a subset `S` of the edges of `G` of maximum weight such that
        no two edges of `S` are incident with each other.

        As an optimization problem, it can be expressed as:

        .. MATH::

            \mbox{Maximize : }&\sum_{e\in G.edges()} w_e b_e\\
            \mbox{Such that : }&\forall v \in G,
            \sum_{(u,v)\in G.edges()} b_{(u,v)}\leq 1\\
            &\forall x\in G, b_x\mbox{ is a binary variable}

        INPUT:

        - ``value_only`` -- boolean (default: ``False``); when set to ``True``,
          only the cardinal (or the weight) of the matching is returned

        - ``algorithm`` -- string (default: ``"Edmonds"``)

          - ``"Edmonds"`` selects Edmonds' algorithm as implemented in NetworkX

          - ``"LP"`` uses a Linear Program formulation of the matching problem

        - ``use_edge_labels`` -- boolean (default: ``False``)

          - when set to ``True``, computes a weighted matching where each edge
            is weighted by its label (if an edge has no label, `1` is assumed)

          - when set to ``False``, each edge has weight `1`

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of verbosity:
          set to 0 by default, which means quiet (only useful when ``algorithm
          == "LP"``)

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        OUTPUT:

        - When ``value_only=False`` (default), this method returns the list of
          edges of a maximum matching of `G`.

        - When ``value_only=True``, this method returns the sum of the
          weights (default: ``1``) of the edges of a maximum matching of `G`.
          The type of the output may vary according to the type of the edge
          labels and the algorithm used.

        ALGORITHM:

        The problem is solved using Edmond's algorithm implemented in NetworkX,
        or using Linear Programming depending on the value of ``algorithm``.

        EXAMPLES:

        Maximum matching in a Pappus Graph::

           sage: g = graphs.PappusGraph()
           sage: g.matching(value_only=True)
           9

        Same test with the Linear Program formulation::

           sage: g = graphs.PappusGraph()
           sage: g.matching(algorithm="LP", value_only=True)
           9

        .. PLOT::

            g = graphs.PappusGraph()
            sphinx_plot(g.plot(edge_colors={"red":g.matching()}))

        TESTS:

        When ``use_edge_labels`` is set to ``False``, with Edmonds' algorithm
        and LP formulation::

            sage: g = Graph([(0,1,0), (1,2,999), (2,3,-5)])
            sage: sorted(g.matching())
            [(0, 1, 0), (2, 3, -5)]
            sage: sorted(g.matching(algorithm="LP"))
            [(0, 1, 0), (2, 3, -5)]

        When ``use_edge_labels`` is set to ``True``, with Edmonds' algorithm and
        LP formulation::

            sage: g = Graph([(0,1,0), (1,2,999), (2,3,-5)])
            sage: g.matching(use_edge_labels=True)
            [(1, 2, 999)]
            sage: g.matching(algorithm="LP", use_edge_labels=True)
            [(1, 2, 999)]

        With loops and multiedges::

            sage: edge_list = [(0,0,5), (0,1,1), (0,2,2), (0,3,3), (1,2,6)
            ....: , (1,2,3), (1,3,3), (2,3,3)]
            sage: g = Graph(edge_list, loops=True, multiedges=True)
            sage: g.matching(use_edge_labels=True)
            [(1, 2, 6), (0, 3, 3)]

        TESTS:

        If ``algorithm`` is set to anything different from ``"Edmonds"`` or
        ``"LP"``, an exception is raised::

           sage: g = graphs.PappusGraph()
           sage: g.matching(algorithm="somethingdifferent")
           Traceback (most recent call last):
           ...
           ValueError: algorithm must be set to either "Edmonds" or "LP"
        """
        from sage.rings.real_mpfr import RR
        def weight(x):
            if x in RR:
                return x
            else:
                return 1

        W = {}
        L = {}
        for u,v,l in self.edge_iterator():
            if u is v:
                continue
            fuv = frozenset((u, v))
            if fuv not in L or ( use_edge_labels and W[fuv] < weight(l) ):
                L[fuv] = l
                if use_edge_labels:
                    W[fuv] = weight(l)

        if algorithm == "Edmonds":
            import networkx
            g = networkx.Graph()
            if use_edge_labels:
                for (u, v),w in W.items():
                    g.add_edge(u, v, weight=w)
            else:
                for u, v in L:
                    g.add_edge(u, v)
            d = networkx.max_weight_matching(g)
            if value_only:
                if use_edge_labels:
                    return sum(W[frozenset(e)] for e in d)
                else:
                    return Integer(len(d))
            else:
                return [(u, v, L[frozenset((u, v))]) for u, v in d]

        elif algorithm == "LP":
            g = self
            from sage.numerical.mip import MixedIntegerLinearProgram
            # returns the weight of an edge considering it may not be
            # weighted ...
            p = MixedIntegerLinearProgram(maximization=True, solver=solver)
            b = p.new_variable(binary=True)
            if use_edge_labels:
                p.set_objective(p.sum(w * b[fe] for fe,w in W.items()))
            else:
                p.set_objective(p.sum(b[fe] for fe in L))
            # for any vertex v, there is at most one edge incident to v in
            # the maximum matching
            for v in g:
                p.add_constraint(p.sum(b[frozenset(e)] for e in self.edge_iterator(vertices=[v], labels=False)
                                           if e[0] != e[1]), max=1)

            p.solve(log=verbose)
            b = p.get_values(b, convert=bool, tolerance=integrality_tolerance)
            if value_only:
                if use_edge_labels:
                    return sum(w for fe, w in W.items() if b[fe])
                else:
                    return Integer(sum(1 for fe in L if b[fe]))
            else:
                return [(u, v, L[frozenset((u, v))]) for u, v in L if b[frozenset((u, v))]]

        else:
            raise ValueError('algorithm must be set to either "Edmonds" or "LP"')

    @doc_index("Algorithmically hard stuff")
    def has_homomorphism_to(self, H, core=False, solver=None, verbose=0,
                            *, integrality_tolerance=1e-3):
        r"""
        Checks whether there is a homomorphism between two graphs.

        A homomorphism from a graph `G` to a graph `H` is a function
        `\phi:V(G)\mapsto V(H)` such that for any edge `uv \in E(G)` the pair
        `\phi(u)\phi(v)` is an edge of `H`.

        Saying that a graph can be `k`-colored is equivalent to saying that it
        has a homomorphism to `K_k`, the complete graph on `k` elements.

        For more information, see the :wikipedia:`Graph_homomorphism`.

        INPUT:

        - ``H`` -- the graph to which ``self`` should be sent.

        - ``core`` -- boolean (default: ``False``; whether to minimize the size
          of the mapping's image (see note below). This is set to ``False`` by
          default.

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        .. NOTE::

           One can compute the core of a graph (with respect to homomorphism)
           with this method ::

               sage: g = graphs.CycleGraph(10)
               sage: mapping = g.has_homomorphism_to(g, core = True)
               sage: print("The size of the core is {}".format(len(set(mapping.values()))))
               The size of the core is 2

        OUTPUT:

        This method returns ``False`` when the homomorphism does not exist, and
        returns the homomorphism otherwise as a dictionary associating a vertex
        of `H` to a vertex of `G`.

        EXAMPLES:

        Is Petersen's graph 3-colorable::

            sage: P = graphs.PetersenGraph()
            sage: P.has_homomorphism_to(graphs.CompleteGraph(3)) is not False
            True

        An odd cycle admits a homomorphism to a smaller odd cycle, but not to an
        even cycle::

            sage: g = graphs.CycleGraph(9)
            sage: g.has_homomorphism_to(graphs.CycleGraph(5)) is not False
            True
            sage: g.has_homomorphism_to(graphs.CycleGraph(7)) is not False
            True
            sage: g.has_homomorphism_to(graphs.CycleGraph(4)) is not False
            False
        """
        self._scream_if_not_simple()
        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
        p = MixedIntegerLinearProgram(solver=solver, maximization=False)
        b = p.new_variable(binary=True)

        # Each vertex has an image
        for ug in self:
            p.add_constraint(p.sum(b[ug,uh] for uh in H) == 1)

        nonedges = H.complement().edges(labels=False)
        for ug,vg in self.edges(labels=False):
            # Two adjacent vertices cannot be mapped to the same element
            for uh in H:
                p.add_constraint(b[ug,uh] + b[vg,uh] <= 1)

            # Two adjacent vertices cannot be mapped to no adjacent vertices
            for uh,vh in nonedges:
                p.add_constraint(b[ug,uh] + b[vg,vh] <= 1)
                p.add_constraint(b[ug,vh] + b[vg,uh] <= 1)

        # Minimize the mapping's size
        if core:

            # the value of m is one if the corresponding vertex of h is used.
            m = p.new_variable(nonnegative=True)
            for uh in H:
                for ug in self:
                    p.add_constraint(b[ug,uh] <= m[uh])

            p.set_objective(p.sum(m[vh] for vh in H))

        try:
            p.solve(log=verbose)
        except MIPSolverException:
            return False

        b = p.get_values(b, convert=bool, tolerance=integrality_tolerance)
        mapping = dict(x[0] for x in b.items() if x[1])
        return mapping


    @doc_index("Clique-related methods")
    def fractional_clique_number(self, solver='PPL', verbose=0,
                                 check_components=True, check_bipartite=True):
        r"""
        Return the fractional clique number of the graph.

        A fractional clique is a nonnegative weight function on the vertices of
        a graph such that the sum of the weights over any independent set is at
        most 1. The fractional clique number is the largest total weight of a
        fractional clique, which is equal to the fractional chromatic number by
        LP-duality.

        ALGORITHM:

        The fractional clique number is computed via the Linear Program for
        fractional chromatic number, see :meth:`fractional_chromatic_number
        <sage.graphs.graph_coloring.fractional_chromatic_number>`

        INPUT:

        - ``solver`` -- (default: ``"PPL"``); specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

          .. NOTE::

              The default solver used here is ``"PPL"`` which provides exact
              results, i.e. a rational number, although this may be slower that
              using other solvers.

        - ``verbose`` -- integer (default: `0`); sets the level of verbosity of
          the LP solver

        - ``check_components`` -- boolean (default: ``True``); whether the
          method is called on each biconnected component of `G`

        - ``check_bipartite`` -- boolean (default: ``True``); whether the graph
          is checked for bipartiteness. If the graph is bipartite then we can
          avoid creating and solving the LP.

        EXAMPLES:

        The fractional clique number of a `C_7` is `7/3`::

            sage: g = graphs.CycleGraph(7)
            sage: g.fractional_clique_number()
            7/3
        """
        return self.fractional_chromatic_number(solver=solver, verbose=verbose,
                                                check_components=check_components,
                                                check_bipartite=check_bipartite)

    @doc_index("Leftovers")
    def maximum_average_degree(self, value_only=True, solver=None, verbose=0):
        r"""
        Return the Maximum Average Degree (MAD) of the current graph.

        The Maximum Average Degree (MAD) of a graph is defined as the average
        degree of its densest subgraph. More formally, ``Mad(G) =
        \max_{H\subseteq G} Ad(H)``, where `Ad(G)` denotes the average degree of
        `G`.

        This can be computed in polynomial time.

        INPUT:

        - ``value_only`` -- boolean (default: ``True``);

          - If ``value_only=True``, only the numerical value of the `MAD` is
            returned.

          - Else, the subgraph of `G` realizing the `MAD` is returned.

        - ``solver`` -- (default: ``None``); specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        EXAMPLES:

        In any graph, the `Mad` is always larger than the average degree::

            sage: g = graphs.RandomGNP(20,.3)
            sage: mad_g = g.maximum_average_degree()
            sage: g.average_degree() <= mad_g
            True

        Unlike the average degree, the `Mad` of the disjoint union of two graphs
        is the maximum of the `Mad` of each graphs::

            sage: h = graphs.RandomGNP(20,.3)
            sage: mad_h = h.maximum_average_degree()
            sage: (g+h).maximum_average_degree() == max(mad_g, mad_h)
            True

        The subgraph of a regular graph realizing the maximum average degree is
        always the whole graph ::

            sage: g = graphs.CompleteGraph(5)
            sage: mad_g = g.maximum_average_degree(value_only=False)
            sage: g.is_isomorphic(mad_g)
            True

        This also works for complete bipartite graphs ::

            sage: g = graphs.CompleteBipartiteGraph(3,4)
            sage: mad_g = g.maximum_average_degree(value_only=False)
            sage: g.is_isomorphic(mad_g)
            True
        """
        self._scream_if_not_simple()
        g = self
        from sage.numerical.mip import MixedIntegerLinearProgram

        p = MixedIntegerLinearProgram(maximization=True, solver=solver)

        d = p.new_variable(nonnegative=True)
        one = p.new_variable(nonnegative=True)

        for u,v in g.edge_iterator(labels=False):
            fuv = frozenset((u, v))
            p.add_constraint(one[fuv] - 2 * d[u], max=0)
            p.add_constraint(one[fuv] - 2 * d[v], max=0)

        p.add_constraint(p.sum(d[v] for v in g), max=1)

        p.set_objective(p.sum(one[frozenset(uv)]
                              for uv in g.edge_iterator(labels=False)))

        p.solve(log=verbose)

        # Paying attention to numerical error :
        # The zero values could be something like 0.000000000001
        # so I can not write l > 0
        # And the non-zero, though they should be equal to
        # 1/(order of the optimal subgraph) may be a bit lower

        # setting the minimum to 1/(10 * size of the whole graph )
        # should be safe :-)
        m = 1/(10 *Integer(g.order()))
        d_val = p.get_values(d)
        g_mad = g.subgraph(v for v,l in d_val.items() if l > m)

        if value_only:
            return g_mad.average_degree()
        else:
            return g_mad

    @doc_index("Algorithmically hard stuff")
    def independent_set_of_representatives(self, family, solver=None, verbose=0,
                                           *, integrality_tolerance=1e-3):
        r"""
        Return an independent set of representatives.

        Given a graph `G` and a family `F=\{F_i:i\in [1,...,k]\}` of subsets of
        ``g.vertices()``, an Independent Set of Representatives (ISR) is an
        assignation of a vertex `v_i\in F_i` to each set `F_i` such that `v_i !=
        v_j` if `i<j` (they are representatives) and the set `\cup_{i}v_i` is an
        independent set in `G`.

        It generalizes, for example, graph coloring and graph list coloring.

        (See [ABZ2007]_ for more information.)

        INPUT:

        - ``family`` -- A list of lists defining the family `F` (actually, a
          Family of subsets of ``G.vertices()``).

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        OUTPUT:

        - A list whose `i^{\mbox{th}}` element is the representative of the
          `i^{\mbox{th}}` element of the ``family`` list. If there is no ISR,
          ``None`` is returned.

        EXAMPLES:

        For a bipartite graph missing one edge, the solution is as expected::

           sage: g = graphs.CompleteBipartiteGraph(3,3)
           sage: g.delete_edge(1,4)
           sage: g.independent_set_of_representatives([[0,1,2],[3,4,5]])
           [1, 4]

        The Petersen Graph is 3-colorable, which can be expressed as an
        independent set of representatives problem : take 3 disjoint copies of
        the Petersen Graph, each one representing one color. Then take as a
        partition of the set of vertices the family defined by the three copies
        of each vertex. The ISR of such a family defines a 3-coloring::

            sage: g = 3 * graphs.PetersenGraph()
            sage: n = g.order() / 3
            sage: f = [[i, i + n, i + 2*n] for i in range(n)]
            sage: isr = g.independent_set_of_representatives(f)
            sage: c = [integer_floor(i / n) for i in isr]
            sage: color_classes = [[], [], []]
            sage: for v, i in enumerate(c):
            ....:   color_classes[i].append(v)
            sage: for classs in color_classes:
            ....:   g.subgraph(classs).size() == 0
            True
            True
            True
        """
        from sage.numerical.mip import MixedIntegerLinearProgram
        p = MixedIntegerLinearProgram(solver=solver)

        # Boolean variable indicating whether the vertex is the representative
        # of some set
        vertex_taken = p.new_variable(binary=True)

        # Boolean variable in two dimension whose first element is a vertex and
        # whose second element is one of the sets given as arguments.
        # When true, indicated that the vertex is the representative of the
        # corresponding set
        classss = p.new_variable(binary=True)

        # Associates to the vertices the classes to which they belong
        lists = {v: [] for v in self}
        for i,f in enumerate(family):
            for v in f:
                lists[v].append(i)

            # a classss has exactly one representative
            p.add_constraint(p.sum(classss[v,i] for v in f), max=1, min=1)

        # A vertex represents at most one classss (vertex_taken is binary), and
        # vertex_taken[v]==1 if v is the representative of some classss
        for v in self:
            p.add_constraint(p.sum(classss[v,i] for i in lists[v]) - vertex_taken[v], max=0)

        # Two adjacent vertices can not both be representatives of a set

        for u,v in self.edge_iterator(labels=None):
            p.add_constraint(vertex_taken[u] + vertex_taken[v], max=1)

        p.set_objective(None)

        try:
            p.solve(log=verbose)
        except Exception:
            return None

        classss = p.get_values(classss, convert=bool, tolerance=integrality_tolerance)

        repr = []
        for i,f in enumerate(family):
            for v in f:
                if classss[v,i]:
                    repr.append(v)
                    break

        return repr

    @doc_index("Algorithmically hard stuff")
    def minor(self, H, solver=None, verbose=0, *, integrality_tolerance=1e-3):
        r"""
        Return the vertices of a minor isomorphic to `H` in the current graph.

        We say that a graph `G` has a `H`-minor (or that it has a graph
        isomorphic to `H` as a minor), if for all `h\in H`, there exist disjoint
        sets `S_h \subseteq V(G)` such that once the vertices of each `S_h` have
        been merged to create a new graph `G'`, this new graph contains `H` as a
        subgraph.

        For more information, see the :wikipedia:`Minor_(graph_theory)`.

        INPUT:

        - ``H`` -- The minor to find for in the current graph.

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        OUTPUT:

        A dictionary associating to each vertex of `H` the set of vertices in
        the current graph representing it.

        ALGORITHM:

        Mixed Integer Linear Programming

        COMPLEXITY:

        Theoretically, when `H` is fixed, testing for the existence of a
        `H`-minor is polynomial. The known algorithms are highly exponential in
        `H`, though.

        .. NOTE::

            This function can be expected to be *very* slow, especially where
            the minor does not exist.

        EXAMPLES:

        Trying to find a minor isomorphic to `K_4` in the `4\times 4` grid::

            sage: g = graphs.GridGraph([4,4])
            sage: h = graphs.CompleteGraph(4)
            sage: L = g.minor(h)
            sage: gg = g.subgraph(flatten(L.values(), max_level = 1))
            sage: _ = [gg.merge_vertices(l) for l in L.values() if len(l)>1]
            sage: gg.is_isomorphic(h)
            True

        We can also try to prove this way that the Petersen graph is not planar,
        as it has a `K_5` minor::

            sage: g = graphs.PetersenGraph()
            sage: K5_minor = g.minor(graphs.CompleteGraph(5))                    # long time

        And even a `K_{3,3}` minor::

            sage: K33_minor = g.minor(graphs.CompleteBipartiteGraph(3,3))        # long time

        (It is much faster to use the linear-time test of planarity in this
        situation, though.)

        As there is no cycle in a tree, looking for a `K_3` minor is useless.
        This function will raise an exception in this case::

            sage: g = graphs.RandomGNP(20,.5)
            sage: g = g.subgraph(edges = g.min_spanning_tree())
            sage: g.is_tree()
            True
            sage: L = g.minor(graphs.CompleteGraph(3))
            Traceback (most recent call last):
            ...
            ValueError: This graph has no minor isomorphic to H !
        """
        self._scream_if_not_simple()
        H._scream_if_not_simple()
        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
        p = MixedIntegerLinearProgram(solver=solver)

        # We use frozenset((u, v)) to avoid confusion between (u, v) and (v, u)

        # rs = Representative set of a vertex
        # for h in H, v in G is such that rs[h,v] == 1 if and only if v
        # is a representative of h in self
        rs = p.new_variable(binary=True)

        for v in self:
            p.add_constraint(p.sum(rs[h,v] for h in H), max=1)

        # We ensure that the set of representatives of a
        # vertex h contains a tree, and thus is connected

        # edges represents the edges of the tree
        edges = p.new_variable(binary=True)

        # there can be a edge for h between two vertices
        # only if those vertices represent h
        for u,v in self.edge_iterator(labels=None):
            fuv = frozenset((u, v))
            for h in H:
                p.add_constraint(edges[h,fuv] - rs[h,u], max=0)
                p.add_constraint(edges[h,fuv] - rs[h,v], max=0)

        # The number of edges of the tree in h is exactly the cardinal
        # of its representative set minus 1

        for h in H:
            p.add_constraint(  p.sum(edges[h,frozenset(e)] for e in self.edge_iterator(labels=None))
                             - p.sum(rs[h,v] for v in self), min=-1, max=-1)

        # a tree  has no cycle
        epsilon = 1/(5*Integer(self.order()))
        r_edges = p.new_variable(nonnegative=True)

        for h in H:
            for u,v in self.edge_iterator(labels=None):
                p.add_constraint(r_edges[h,(u,v)] + r_edges[h,(v,u)] - edges[h,frozenset((u,v))], min=0)

            for v in self:
                p.add_constraint(p.sum(r_edges[h,(u,v)] for u in self.neighbor_iterator(v)), max=1-epsilon)

        # Once the representative sets are described, we must ensure
        # there are arcs corresponding to those of H between them
        h_edges = p.new_variable(nonnegative=True)

        for h1, h2 in H.edge_iterator(labels=None):

            for v1, v2 in self.edge_iterator(labels=None):
                fv1v2 = frozenset((v1, v2))
                p.add_constraint(h_edges[(h1,h2),fv1v2] - rs[h2,v2], max=0)
                p.add_constraint(h_edges[(h1,h2),fv1v2] - rs[h1,v1], max=0)

                p.add_constraint(h_edges[(h2,h1),fv1v2] - rs[h1,v2], max=0)
                p.add_constraint(h_edges[(h2,h1),fv1v2] - rs[h2,v1], max=0)

            p.add_constraint(p.sum(h_edges[(h1,h2),frozenset(e)] + h_edges[(h2,h1),frozenset(e)]
                                       for e in self.edge_iterator(labels=None)), min=1)

        p.set_objective(None)

        try:
            p.solve(log=verbose)
        except MIPSolverException:
            raise ValueError("This graph has no minor isomorphic to H !")

        rs = p.get_values(rs, convert=bool, tolerance=integrality_tolerance)

        rs_dict = {}
        for h in H:
            rs_dict[h] = [v for v in self if rs[h,v]]

        return rs_dict

    ### Convexity

    @doc_index("Algorithmically hard stuff")
    def convexity_properties(self):
        r"""
        Return a ``ConvexityProperties`` object corresponding to ``self``.

        This object contains the methods related to convexity in graphs (convex
        hull, hull number) and caches useful information so that it becomes
        comparatively cheaper to compute the convex hull of many different sets
        of the same graph.

        .. SEEALSO::

            In order to know what can be done through this object, please refer
            to module :mod:`sage.graphs.convexity_properties`.

        .. NOTE::

            If you want to compute many convex hulls, keep this object in memory
            ! When it is created, it builds a table of useful information to
            compute convex hulls. As a result ::

                sage: g = graphs.PetersenGraph()
                sage: g.convexity_properties().hull([1, 3])
                [1, 2, 3]
                sage: g.convexity_properties().hull([3, 7])
                [2, 3, 7]

            Is a terrible waste of computations, while ::

                sage: g = graphs.PetersenGraph()
                sage: CP = g.convexity_properties()
                sage: CP.hull([1, 3])
                [1, 2, 3]
                sage: CP.hull([3, 7])
                [2, 3, 7]

            Makes perfect sense.
        """
        from sage.graphs.convexity_properties import ConvexityProperties
        return ConvexityProperties(self)

    # Centrality
    @doc_index("Distances")
    def centrality_degree(self, v=None):
        r"""
        Return the degree centrality of a vertex.

        The degree centrality of a vertex `v` is its degree, divided by
        `|V(G)|-1`. For more information, see the :wikipedia:`Centrality`.

        INPUT:

        - ``v`` -- a vertex (default: ``None``); set to ``None`` (default) to
          get a dictionary associating each vertex with its centrality degree.

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.centrality_closeness`
            - :meth:`~sage.graphs.generic_graph.GenericGraph.centrality_betweenness`

        EXAMPLES::

            sage: (graphs.ChvatalGraph()).centrality_degree()
            {0: 4/11, 1: 4/11, 2: 4/11, 3: 4/11,  4: 4/11,  5: 4/11,
             6: 4/11, 7: 4/11, 8: 4/11, 9: 4/11, 10: 4/11, 11: 4/11}
            sage: D = graphs.DiamondGraph()
            sage: D.centrality_degree()
            {0: 2/3, 1: 1, 2: 1, 3: 2/3}
            sage: D.centrality_degree(v=1)
            1

        TESTS::

            sage: Graph(1).centrality_degree()
            Traceback (most recent call last):
            ...
            ValueError: the centrality degree is not defined on graphs with only one vertex
        """
        from sage.rings.integer import Integer
        n_minus_one = Integer(self.order() - 1)
        if n_minus_one == 0:
            raise ValueError("the centrality degree is not defined "
                             "on graphs with only one vertex")
        if v is None:
            return {v: self.degree(v)/n_minus_one for v in self}
        else:
            return self.degree(v)/n_minus_one

    ### Distances

    @doc_index("Distances")
    def eccentricity(self, v=None, by_weight=False, algorithm=None,
                     weight_function=None, check_weight=True, dist_dict=None,
                     with_labels=False):
        """
        Return the eccentricity of vertex (or vertices) ``v``.

        The eccentricity of a vertex is the maximum distance to any other
        vertex.

        For more information and examples on how to use input variables, see
        :meth:`~GenericGraph.shortest_path_all_pairs`,
        :meth:`~GenericGraph.shortest_path_lengths` and
        :meth:`~GenericGraph.shortest_paths`

        INPUT:

        - ``v`` - either a single vertex or a list of vertices. If it is not
          specified, then it is taken to be all vertices.

        - ``by_weight`` -- boolean (default: ``False``); if ``True``, edge
          weights are taken into account; if False, all edges have weight 1

        - ``algorithm`` -- string (default: ``None``); one of the following
          algorithms:

          - ``'BFS'`` - the computation is done through a BFS centered on each
            vertex successively. Works only if ``by_weight==False``.

          - ``'DHV'`` - the computation is done using the algorithm proposed in
            [Dragan2018]_. Works only if ``self`` has non-negative edge weights
            and ``v is None`` or ``v`` should contain all vertices of ``self``.
            For more information see method
            :func:`sage.graphs.distances_all_pairs.eccentricity` and
            :func:`sage.graphs.base.boost_graph.eccentricity_DHV`.

          - ``'Floyd-Warshall-Cython'`` - a Cython implementation of the
            Floyd-Warshall algorithm. Works only if ``by_weight==False`` and
            ``v is None`` or ``v`` should contain all vertices of ``self``.

          - ``'Floyd-Warshall-Python'`` - a Python implementation of the
            Floyd-Warshall algorithm. Works also with weighted graphs, even with
            negative weights (but no negative cycle is allowed). However, ``v``
            must be ``None`` or ``v`` should contain all vertices of ``self``.

          - ``'Dijkstra_NetworkX'`` - the Dijkstra algorithm, implemented in
            NetworkX. It works with weighted graphs, but no negative weight is
            allowed.

          - ``'Dijkstra_Boost'`` - the Dijkstra algorithm, implemented in Boost
            (works only with positive weights).

          - ``'Johnson_Boost'`` - the Johnson algorithm, implemented in
            Boost (works also with negative weights, if there is no negative
            cycle). Works only if ``v is None`` or ``v`` should contain all
            vertices of ``self``.

          - ``'From_Dictionary'`` - uses the (already computed) distances, that
            are provided by input variable ``dist_dict``.

          - ``None`` (default): Sage chooses the best algorithm:
            ``'From_Dictionary'`` if ``dist_dict`` is not None, ``'BFS'`` for
            unweighted graphs, ``'Dijkstra_Boost'`` if all weights are
            positive, ``'Johnson_Boost'`` otherwise.

        - ``weight_function`` -- function (default: ``None``); a function that
          takes as input an edge ``(u, v, l)`` and outputs its weight. If not
          ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
          and ``by_weight`` is ``True``, we use the edge label ``l`` as a
          weight, if ``l`` is not ``None``, else ``1`` as a weight.

        - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
          that the ``weight_function`` outputs a number for each edge

        - ``dist_dict`` -- a dictionary (default: ``None``); a dict of dicts of
          distances (used only if ``algorithm=='From_Dictionary'``)

        - ``with_labels`` -- boolean (default: ``False``); whether to return a
          list or a dictionary keyed by vertices.

        EXAMPLES::

            sage: G = graphs.KrackhardtKiteGraph()
            sage: G.eccentricity()
            [4, 4, 4, 4, 4, 3, 3, 2, 3, 4]
            sage: G.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: G.eccentricity(7)
            2
            sage: G.eccentricity([7,8,9])
            [2, 3, 4]
            sage: G.eccentricity([7,8,9], with_labels=True) == {8: 3, 9: 4, 7: 2}
            True
            sage: G = Graph( { 0 : [], 1 : [], 2 : [1] } )
            sage: G.eccentricity()
            [+Infinity, +Infinity, +Infinity]
            sage: G = Graph({0:[]})
            sage: G.eccentricity(with_labels=True)
            {0: 0}
            sage: G = Graph({0:[], 1:[]})
            sage: G.eccentricity(with_labels=True)
            {0: +Infinity, 1: +Infinity}
            sage: G = Graph([(0,1,1), (1,2,1), (0,2,3)])
            sage: G.eccentricity(algorithm = 'BFS')
            [1, 1, 1]
            sage: G.eccentricity(algorithm = 'Floyd-Warshall-Cython')
            [1, 1, 1]
            sage: G.eccentricity(by_weight = True, algorithm = 'Dijkstra_NetworkX')
            [2, 1, 2]
            sage: G.eccentricity(by_weight = True, algorithm = 'Dijkstra_Boost')
            [2, 1, 2]
            sage: G.eccentricity(by_weight = True, algorithm = 'Johnson_Boost')
            [2, 1, 2]
            sage: G.eccentricity(by_weight = True, algorithm = 'Floyd-Warshall-Python')
            [2, 1, 2]
            sage: G.eccentricity(dist_dict = G.shortest_path_all_pairs(by_weight = True)[0])
            [2, 1, 2]
            sage: G.eccentricity(by_weight = False, algorithm = 'DHV')
            [1, 1, 1]
            sage: G.eccentricity(by_weight = True, algorithm = 'DHV')
            [2.0, 1.0, 2.0]

        TESTS:

        A non-implemented algorithm::

            sage: G.eccentricity(algorithm = 'boh')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm "boh"

        An algorithm that does not work with edge weights::

            sage: G.eccentricity(by_weight = True, algorithm = 'BFS')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'BFS' does not work with weights
            sage: G.eccentricity(by_weight = True, algorithm = 'Floyd-Warshall-Cython')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'Floyd-Warshall-Cython' does not work with weights

        An algorithm that computes the all-pair-shortest-paths when not all
        vertices are needed::

            sage: G.eccentricity(0, algorithm = 'Floyd-Warshall-Cython')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'Floyd-Warshall-Cython' works only if all eccentricities are needed
            sage: G.eccentricity(0, algorithm = 'Floyd-Warshall-Python')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'Floyd-Warshall-Python' works only if all eccentricities are needed
            sage: G.eccentricity(0, algorithm = 'Johnson_Boost')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'Johnson_Boost' works only if all eccentricities are needed
            sage: G.eccentricity(0, algorithm = 'DHV')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'DHV' works only if all eccentricities are needed
        """
        by_weight, weight_function = self._get_weight_function(by_weight=by_weight,
                                                               weight_function=weight_function,
                                                               check_weight=check_weight)

        if algorithm is None:
            if dist_dict is not None:
                algorithm = 'From_Dictionary'
            elif not by_weight:
                algorithm = 'BFS'
            elif any(float(weight_function(e)) < 0 for e in self.edge_iterator()):
                algorithm = 'Johnson_Boost'
            if algorithm is None:
                algorithm = 'Dijkstra_Boost'
        if algorithm in ['BFS', 'Floyd-Warshall-Cython']:
            if by_weight:
                raise ValueError("algorithm '{}' does not work with weights".format(algorithm))
            # We don't want the default weight function
            weight_function = None

        if v is not None:
            if not isinstance(v, list):
                v = [v]
            v_set = set(v)

        if v is None or all(u in v_set for u in self):
            if v is None:
                v = list(self)
            # If we want to use BFS, we use the Cython routine
            if algorithm == 'BFS':
                from sage.graphs.distances_all_pairs import eccentricity
                algo = 'bounds'
                if with_labels:
                    return dict(zip(v, eccentricity(self, algorithm=algo, vertex_list=v)))
                else:
                    return eccentricity(self, algorithm=algo,vertex_list=v)

            if algorithm == 'DHV':
                if by_weight:
                    from sage.graphs.base.boost_graph import eccentricity_DHV
                    if with_labels:
                        return dict(zip(v, eccentricity_DHV(self, vertex_list=v,
                                                            weight_function=weight_function,
                                                            check_weight=check_weight)))
                    else:
                        return eccentricity_DHV(self, vertex_list=v,
                                                weight_function=weight_function,
                                                check_weight=check_weight)
                else:
                    from sage.graphs.distances_all_pairs import eccentricity
                    if with_labels:
                        return dict(zip(v, eccentricity(self, algorithm=algorithm,
                                                        vertex_list=v)))
                    else:
                        return eccentricity(self, algorithm=algorithm, vertex_list=v)

            if algorithm in ['Floyd-Warshall-Python', 'Floyd-Warshall-Cython', 'Johnson_Boost']:
                dist_dict = self.shortest_path_all_pairs(by_weight, algorithm,
                                                         weight_function,
                                                         check_weight)[0]
                algorithm = 'From_Dictionary'

        elif algorithm in ['Floyd-Warshall-Python', 'Floyd-Warshall-Cython', 'Johnson_Boost','DHV']:
            raise ValueError("algorithm '" + algorithm + "' works only if all" +
                             " eccentricities are needed")

        ecc = {}

        from sage.rings.infinity import Infinity

        for u in v:
            if algorithm == 'From_Dictionary':
                length = dist_dict[u]
            else:
                # If algorithm is wrong, the error is raised by the
                # shortest_path_lengths function
                length = self.shortest_path_lengths(u, by_weight=by_weight,
                                                    algorithm=algorithm,
                                                    weight_function=weight_function,
                                                    check_weight=check_weight)

            if len(length) != self.num_verts():
                ecc[u] = Infinity
            else:
                ecc[u] = max(length.values())

        if with_labels:
            return ecc
        else:
            if len(ecc) == 1:
                # return single value
                v, = ecc.values()
                return v
            return [ecc[u] for u in v]

    @doc_index("Distances")
    def radius(self, by_weight=False, algorithm='DHV', weight_function=None,
               check_weight=True):
        r"""
        Return the radius of the graph.

        The radius is defined to be the minimum eccentricity of any vertex,
        where the eccentricity is the maximum distance to any other
        vertex. For more information and examples on how to use input variables,
        see :meth:`~GenericGraph.shortest_paths` and
        :meth:`~Graph.eccentricity`

        INPUT:

        - ``by_weight`` -- boolean (default: ``False``); if ``True``, edge
          weights are taken into account; if False, all edges have weight 1

        - ``algorithm`` -- string (default: ``'DHV'``).

          - ``'DHV'`` - Radius computation is done using the algorithm proposed
            in [Dragan2018]_. Works for graph with non-negative edge weights.
            For more information see method
            :func:`sage.graphs.distances_all_pairs.radius_DHV` and
            :func:`sage.graphs.base.boost_graph.radius_DHV`.

          - see method :meth:`eccentricity` for the list of remaining algorithms

        - ``weight_function`` -- function (default: ``None``); a function that
          takes as input an edge ``(u, v, l)`` and outputs its weight. If not
          ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
          and ``by_weight`` is ``True``, we use the edge label ``l`` as a
          weight, if ``l`` is not ``None``, else ``1`` as a weight.

        - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
          that the ``weight_function`` outputs a number for each edge

        EXAMPLES:

        The more symmetric a graph is, the smaller (diameter - radius) is::

            sage: G = graphs.BarbellGraph(9, 3)
            sage: G.radius()
            3
            sage: G.diameter()
            6

        ::

            sage: G = graphs.OctahedralGraph()
            sage: G.radius()
            2
            sage: G.diameter()
            2

        TESTS::

            sage: g = Graph()
            sage: g.radius()
            Traceback (most recent call last):
            ...
            ValueError: radius is not defined for the empty graph
        """
        if not self.order():
            raise ValueError("radius is not defined for the empty graph")

        if not algorithm:
            algorithm = 'DHV'

        if algorithm == 'DHV':
            by_weight, weight_function = self._get_weight_function(by_weight=by_weight,
                                                                   weight_function=weight_function,
                                                                   check_weight=check_weight)
            if by_weight:
                from sage.graphs.base.boost_graph import radius_DHV
                return radius_DHV(self, weight_function=weight_function,
                                  check_weight=False)
            else:
                from sage.graphs.distances_all_pairs import radius_DHV
                return radius_DHV(self)

        return min(self.eccentricity(v=None, by_weight=by_weight,
                                     weight_function=weight_function,
                                     check_weight=check_weight,
                                     algorithm=algorithm))

    @doc_index("Distances")
    def diameter(self, by_weight=False, algorithm=None, weight_function=None,
                 check_weight=True):
        r"""
        Return the diameter of the graph.

        The diameter is defined to be the maximum distance between two vertices.
        It is infinite if the graph is not connected.

        For more information and examples on how to use input variables, see
        :meth:`~GenericGraph.shortest_paths` and
        :meth:`~Graph.eccentricity`

        INPUT:

        - ``by_weight`` -- boolean (default: ``False``); if ``True``, edge
          weights are taken into account; if False, all edges have weight 1

        - ``algorithm`` -- string (default: ``None``); one of the following
          algorithms:

          - ``'BFS'``: the computation is done through a BFS centered on each
            vertex successively. Works only if ``by_weight==False``.

          - ``'Floyd-Warshall-Cython'``: a Cython implementation of the
            Floyd-Warshall algorithm. Works only if ``by_weight==False`` and ``v
            is None``.

          - ``'Floyd-Warshall-Python'``: a Python implementation of the
            Floyd-Warshall algorithm. Works also with weighted graphs, even with
            negative weights (but no negative cycle is allowed). However, ``v``
            must be ``None``.

          - ``'Dijkstra_NetworkX'``: the Dijkstra algorithm, implemented in
            NetworkX. It works with weighted graphs, but no negative weight is
            allowed.

          - ``'DHV'`` - diameter computation is done using the algorithm
            proposed in [Dragan2018]_. Works only for non-negative edge weights.
            For more information see method
            :func:`sage.graphs.distances_all_pairs.diameter_DHV` and
            :func:`sage.graphs.base.boost_graph.diameter_DHV`.

          - ``'standard'``, ``'2sweep'``, ``'multi-sweep'``, ``'iFUB'``:
            these algorithms are implemented in
            :func:`sage.graphs.distances_all_pairs.diameter`
            They work only if ``by_weight==False``. See the function
            documentation for more information.

          - ``'Dijkstra_Boost'``: the Dijkstra algorithm, implemented in Boost
            (works only with positive weights).

          - ``'Johnson_Boost'``: the Johnson algorithm, implemented in
            Boost (works also with negative weights, if there is no negative
            cycle).

          - ``None`` (default): Sage chooses the best algorithm: ``'iFUB'`` for
            unweighted graphs, ``'Dijkstra_Boost'`` if all weights are positive,
            ``'Johnson_Boost'`` otherwise.

        - ``weight_function`` -- function (default: ``None``); a function that
          takes as input an edge ``(u, v, l)`` and outputs its weight. If not
          ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
          and ``by_weight`` is ``True``, we use the edge label ``l`` as a
          weight, if ``l`` is not ``None``, else ``1`` as a weight.

        - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
          that the ``weight_function`` outputs a number for each edge

        EXAMPLES:

        The more symmetric a graph is, the smaller (diameter - radius) is::

            sage: G = graphs.BarbellGraph(9, 3)
            sage: G.radius()
            3
            sage: G.diameter()
            6

        ::

            sage: G = graphs.OctahedralGraph()
            sage: G.radius()
            2
            sage: G.diameter()
            2

        TESTS::

            sage: g = Graph()
            sage: g.diameter()
            Traceback (most recent call last):
            ...
            ValueError: diameter is not defined for the empty graph
            sage: g = Graph([(1, 2, {'weight': 1})])
            sage: g.diameter(algorithm='iFUB', weight_function=lambda e: e[2]['weight'])
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'iFUB' does not work on weighted graphs
        """
        if not self.order():
            raise ValueError("diameter is not defined for the empty graph")

        by_weight, weight_function = self._get_weight_function(by_weight=by_weight,
                                                               weight_function=weight_function,
                                                               check_weight=check_weight)
        if not by_weight:
            # We don't want the default weight function
            weight_function = None

        if algorithm is None:
            if by_weight:
                algorithm = 'iFUB'
            else:
                algorithm = 'DHV'
        elif algorithm == 'BFS':
            algorithm = 'standard'

        if algorithm == 'DHV':
            if by_weight:
                from sage.graphs.base.boost_graph import diameter_DHV
                return diameter_DHV(self, weight_function=weight_function,
                                    check_weight=False)
            else:
                from sage.graphs.distances_all_pairs import diameter
                return diameter(self, algorithm=algorithm)

        if algorithm in ['standard', '2sweep', 'multi-sweep', 'iFUB']:
            if by_weight:
                raise ValueError("algorithm '" + algorithm + "' does not work" +
                                 " on weighted graphs")
            from sage.graphs.distances_all_pairs import diameter
            return diameter(self, algorithm=algorithm)

        return max(self.eccentricity(v=list(self), by_weight=by_weight,
                                     weight_function=weight_function,
                                     check_weight=False,
                                     algorithm=algorithm))

    @doc_index("Distances")
    def center(self, by_weight=False, algorithm=None, weight_function=None,
               check_weight=True):
        r"""
        Return the set of vertices in the center of the graph.

        The center is the set of vertices whose eccentricity is equal to the
        radius of the graph, i.e., achieving the minimum eccentricity.

        For more information and examples on how to use input variables,
        see :meth:`~GenericGraph.shortest_paths` and
        :meth:`~Graph.eccentricity`

        INPUT:

        - ``by_weight`` -- boolean (default: ``False``); if ``True``, edge
          weights are taken into account; if False, all edges have weight 1

        - ``algorithm`` -- string (default: ``None``); see method
          :meth:`eccentricity` for the list of available algorithms

        - ``weight_function`` -- function (default: ``None``); a function that
          takes as input an edge ``(u, v, l)`` and outputs its weight. If not
          ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
          and ``by_weight`` is ``True``, we use the edge label ``l`` as a
          weight, if ``l`` is not ``None``, else ``1`` as a weight.

        - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
          that the ``weight_function`` outputs a number for each edge

        EXAMPLES:

        Is Central African Republic in the center of Africa in graph theoretic
        sense? Yes::

            sage: A = graphs.AfricaMap(continental=True)
            sage: sorted(A.center())
            ['Cameroon', 'Central Africa']

        Some other graphs. Center can be the whole graph::

            sage: G = graphs.DiamondGraph()
            sage: G.center()
            [1, 2]
            sage: P = graphs.PetersenGraph()
            sage: P.subgraph(P.center()) == P
            True
            sage: S = graphs.StarGraph(19)
            sage: S.center()
            [0]

        TESTS::

            sage: G = Graph()
            sage: G.center()
            []
            sage: G.add_vertex()
            0
            sage: G.center()
            [0]
        """
        ecc = self.eccentricity(v=list(self), by_weight=by_weight,
                                weight_function=weight_function,
                                algorithm=algorithm,
                                check_weight=check_weight,
                                with_labels=True)
        try:
            r = min(ecc.values())
        except Exception:
            return []
        return [v for v in self if ecc[v] == r]

    @doc_index("Distances")
    def periphery(self, by_weight=False, algorithm=None, weight_function=None,
                  check_weight=True):
        r"""
        Return the set of vertices in the periphery of the graph.

        The periphery is the set of vertices whose eccentricity is equal to the
        diameter of the graph, i.e., achieving the maximum eccentricity.

        For more information and examples on how to use input variables,
        see :meth:`~GenericGraph.shortest_paths` and
        :meth:`~Graph.eccentricity`

        INPUT:

        - ``by_weight`` -- boolean (default: ``False``); if ``True``, edge
          weights are taken into account; if False, all edges have weight 1

        - ``algorithm`` -- string (default: ``None``); see method
          :meth:`eccentricity` for the list of available algorithms

        - ``weight_function`` -- function (default: ``None``); a function that
          takes as input an edge ``(u, v, l)`` and outputs its weight. If not
          ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
          and ``by_weight`` is ``True``, we use the edge label ``l`` as a
          weight, if ``l`` is not ``None``, else ``1`` as a weight.

        - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
          that the ``weight_function`` outputs a number for each edge

        EXAMPLES::

            sage: G = graphs.DiamondGraph()
            sage: G.periphery()
            [0, 3]
            sage: P = graphs.PetersenGraph()
            sage: P.subgraph(P.periphery()) == P
            True
            sage: S = graphs.StarGraph(19)
            sage: S.periphery()
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: G = Graph()
            sage: G.periphery()
            []
            sage: G.add_vertex()
            0
            sage: G.periphery()
            [0]
        """
        ecc = self.eccentricity(v=list(self), by_weight=by_weight,
                                weight_function=weight_function,
                                algorithm=algorithm,
                                check_weight=check_weight,
                                with_labels=True)
        try:
            d = max(ecc.values())
        except Exception:
            return []
        return [v for v in self if ecc[v] == d]

    ### Constructors

    @doc_index("Basic methods")
    def to_directed(self, data_structure=None, sparse=None):
        """
        Return a directed version of the graph.

        A single edge becomes two edges, one in each direction.

        INPUT:

         - ``data_structure`` -- one of ``"sparse"``, ``"static_sparse"``, or
           ``"dense"``. See the documentation of :class:`Graph` or
           :class:`DiGraph`.

         - ``sparse`` -- boolean (default: ``None``); ``sparse=True`` is an
           alias for ``data_structure="sparse"``, and ``sparse=False`` is an
           alias for ``data_structure="dense"``.

        EXAMPLES::

            sage: graphs.PetersenGraph().to_directed()
            Petersen graph: Digraph on 10 vertices

        TESTS:

        Immutable graphs yield immutable graphs::

            sage: Graph([[1, 2]], immutable=True).to_directed()._backend
            <sage.graphs.base.static_sparse_backend.StaticSparseBackend object at ...>

        :trac:`17005`::

            sage: Graph([[1,2]], immutable=True).to_directed()
            Digraph on 2 vertices

        :trac:`22424`::

            sage: G1=graphs.RandomGNP(5,0.5)
            sage: gp1 = G1.graphplot(save_pos=True)
            sage: G2=G1.to_directed()
            sage: G2.delete_vertex(0)
            sage: G2.add_vertex(5)
            sage: gp2 = G2.graphplot()
            sage: gp1 = G1.graphplot()

        Vertex labels will be retained (:trac:`14708`)::

            sage: G = Graph({0: [1, 2], 1: [0]})
            sage: G.set_vertex(0, 'foo')
            sage: D = G.to_directed()
            sage: G.get_vertices()
            {0: 'foo', 1: None, 2: None}
            sage: D.get_vertices()
            {0: 'foo', 1: None, 2: None}
        """
        if sparse is not None:
            if data_structure is not None:
                raise ValueError("The 'sparse' argument is an alias for "
                                 "'data_structure'. Please do not define both.")
            data_structure = "sparse" if sparse else "dense"

        if data_structure is None:
            from sage.graphs.base.dense_graph import DenseGraphBackend
            from sage.graphs.base.sparse_graph import SparseGraphBackend
            if isinstance(self._backend, DenseGraphBackend):
                data_structure = "dense"
            elif isinstance(self._backend, SparseGraphBackend):
                data_structure = "sparse"
            else:
                data_structure = "static_sparse"
        from sage.graphs.digraph import DiGraph
        D = DiGraph(name           = self.name(),
                    pos            = self.get_pos(),
                    multiedges     = self.allows_multiple_edges(),
                    loops          = self.allows_loops(),
                    data_structure = (data_structure if data_structure!="static_sparse"
                                      else "sparse")) # we need a mutable copy

        D.add_vertices(self.vertex_iterator())
        D.set_vertices(self.get_vertices())
        for u,v,l in self.edge_iterator():
            D.add_edge(u,v,l)
            D.add_edge(v,u,l)
        if hasattr(self, '_embedding'):
            D._embedding = copy(self._embedding)
        D._weighted = self._weighted

        if data_structure == "static_sparse":
            D = D.copy(data_structure=data_structure)

        return D

    @doc_index("Basic methods")
    def to_undirected(self):
        """
        Since the graph is already undirected, simply returns a copy of itself.

        EXAMPLES::

            sage: graphs.PetersenGraph().to_undirected()
            Petersen graph: Graph on 10 vertices
        """
        return self.copy()

    @doc_index("Basic methods")
    def join(self, other, labels="pairs", immutable=None):
        r"""
        Return the join of ``self`` and ``other``.

        INPUT:

        - ``labels`` -- (defaults to 'pairs'); if set to 'pairs', each element
          `v` in the first graph will be named `(0, v)` and each element `u` in
          ``other`` will be named `(1, u)` in the result. If set to 'integers',
          the elements of the result will be relabeled with consecutive
          integers.

        - ``immutable`` -- boolean (default: ``None``); whether to create a
          mutable/immutable join. ``immutable=None`` (default) means that the
          graphs and their join will behave the same way.

        .. SEEALSO::

            * :meth:`~sage.graphs.generic_graph.GenericGraph.union`

            * :meth:`~sage.graphs.generic_graph.GenericGraph.disjoint_union`

        EXAMPLES::

            sage: G = graphs.CycleGraph(3)
            sage: H = Graph(2)
            sage: J = G.join(H); J
            Cycle graph join : Graph on 5 vertices
            sage: J.vertices()
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
            sage: J = G.join(H, labels='integers'); J
            Cycle graph join : Graph on 5 vertices
            sage: J.vertices()
            [0, 1, 2, 3, 4]
            sage: J.edges()
            [(0, 1, None), (0, 2, None), (0, 3, None), (0, 4, None), (1, 2, None), (1, 3, None), (1, 4, None), (2, 3, None), (2, 4, None)]

        ::

            sage: G = Graph(3)
            sage: G.name("Graph on 3 vertices")
            sage: H = Graph(2)
            sage: H.name("Graph on 2 vertices")
            sage: J = G.join(H); J
            Graph on 3 vertices join Graph on 2 vertices: Graph on 5 vertices
            sage: J.vertices()
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
            sage: J = G.join(H, labels='integers'); J
            Graph on 3 vertices join Graph on 2 vertices: Graph on 5 vertices
            sage: J.edges()
            [(0, 3, None), (0, 4, None), (1, 3, None), (1, 4, None), (2, 3, None), (2, 4, None)]
        """
        G = self.disjoint_union(other, labels=labels, immutable=False)
        if labels == "integers":
            G.add_edges((u, v) for u in range(self.order())
                        for v in range(self.order(), self.order() + other.order()))
        else:
            G.add_edges(((0, u), (1, v)) for u in self for v in other)

        G.name('%s join %s'%(self.name(), other.name()))

        if immutable is None:
            immutable = self.is_immutable() and other.is_immutable()
        if immutable:
            G = G.copy(immutable=True)

        return G

    @doc_index("Leftovers")
    def seidel_adjacency_matrix(self, vertices=None):
        r"""
        Return the Seidel adjacency matrix of ``self``.

        Returns `J-I-2A`, for `A` the (ordinary) :meth:`adjacency matrix
        <sage.graphs.generic_graph.GenericGraph.adjacency_matrix>` of ``self``,
        `I` the identity matrix, and `J` the all-1 matrix.  It is closely
        related to :meth:`twograph`.

        The matrix returned is over the integers. If a different ring is
        desired, use either the :meth:`sage.matrix.matrix0.Matrix.change_ring`
        method or the :func:`matrix` function.

        INPUT:

        - ``vertices`` -- list of vertices (default: ``None``); the ordering of
          the vertices defining how they should appear in the matrix. By
          default, the ordering given by
          :meth:`~sage.graphs.generic_graph.GenericGraph.vertices` is used.

        EXAMPLES::

            sage: G = graphs.CycleGraph(5)
            sage: G = G.disjoint_union(graphs.CompleteGraph(1))
            sage: G.seidel_adjacency_matrix().minpoly()
            x^2 - 5
        """
        return - self.adjacency_matrix(sparse=False, vertices=vertices) \
               + self.complement().adjacency_matrix(sparse=False, vertices=vertices)

    @doc_index("Leftovers")
    def seidel_switching(self, s, inplace=True):
        r"""
        Return the Seidel switching of ``self`` w.r.t. subset of vertices ``s``.

        Returns the graph obtained by Seidel switching of ``self`` with respect
        to the subset of vertices ``s``. This is the graph given by Seidel
        adjacency matrix `DSD`, for `S` the Seidel adjacency matrix of ``self``,
        and `D` the diagonal matrix with -1s at positions corresponding to
        ``s``, and 1s elsewhere.

        INPUT:

         - ``s`` -- a list of vertices of ``self``.

        - ``inplace`` -- boolean (default: ``True``); whether to do the
          modification inplace, or to return a copy of the graph after
          switching.

        EXAMPLES::

            sage: G = graphs.CycleGraph(5)
            sage: G = G.disjoint_union(graphs.CompleteGraph(1))
            sage: G.seidel_switching([(0,1),(1,0),(0,0)])
            sage: G.seidel_adjacency_matrix().minpoly()
            x^2 - 5
            sage: G.is_connected()
            True

        TESTS::

            sage: H = G.seidel_switching([1,4,5],inplace=False)
            sage: G.seidel_switching([1,4,5])
            sage: G == H
            True
        """
        G = self if inplace else copy(self)
        boundary = self.edge_boundary(s)
        G.add_edges(itertools.product(s, set(self).difference(s)))
        G.delete_edges(boundary)
        if not inplace:
            return G

    @doc_index("Leftovers")
    def twograph(self):
        r"""
        Return the two-graph of ``self``

        Returns the :class:`two-graph <sage.combinat.designs.twographs.TwoGraph>`
        with the triples
        `T=\{t \in \binom {V}{3} : \left| \binom {t}{2} \cap E \right| \text{odd} \}`
        where `V` and `E` are vertices and edges of ``self``, respectively.

        EXAMPLES::

            sage: p=graphs.PetersenGraph()
            sage: p.twograph()
            Incidence structure with 10 points and 60 blocks
            sage: p=graphs.chang_graphs()
            sage: T8 = graphs.CompleteGraph(8).line_graph()
            sage: C = T8.seidel_switching([(0,1,None),(2,3,None),(4,5,None),(6,7,None)],inplace=False)
            sage: T8.twograph() == C.twograph()
            True
            sage: T8.is_isomorphic(C)
            False

        TESTS::

            sage: from sage.combinat.designs.twographs import TwoGraph
            sage: p=graphs.PetersenGraph().twograph()
            sage: TwoGraph(p, check=True)
            Incidence structure with 10 points and 60 blocks

        .. SEEALSO::

            - :meth:`~sage.combinat.designs.twographs.TwoGraph.descendant` --
              computes the descendant graph of the two-graph of self at a vertex

            - :func:`~sage.combinat.designs.twographs.twograph_descendant`
              -- ditto, but much faster.
        """
        from sage.combinat.designs.twographs import TwoGraph
        G = self.relabel(range(self.order()), inplace=False)
        T = []

        # Triangles
        for x,y,z in G.subgraph_search_iterator(Graph({1:[2,3], 2:[3]})):
            if x < y and y < z:
                T.append([x, y, z])

        # Triples with just one edge
        for x,y,z in G.subgraph_search_iterator(Graph({1:[2], 3:[]}), induced=True):
            if x < y:
                T.append([x, y, z])

        T = TwoGraph(T)
        T.relabel({i: v for i,v in enumerate(self.vertices())})

        return T

    ### Visualization

    @doc_index("Basic methods")
    def write_to_eps(self, filename, **options):
        r"""
        Write a plot of the graph to ``filename`` in ``eps`` format.

        INPUT:

         - ``filename`` -- a string
         - ``**options`` -- same layout options as :meth:`.layout`

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.write_to_eps(tmp_filename(ext='.eps'))

        It is relatively simple to include this file in a LaTeX document.
        ``\usepackage{graphics}`` must appear in the preamble, and
        ``\includegraphics{filename}`` will include the file. To compile the
        document to ``pdf`` with ``pdflatex`` or ``xelatex`` the file needs
        first to be converted to ``pdf``, for example with ``ps2pdf filename.eps
        filename.pdf``.
        """
        from sage.graphs.print_graphs import print_graph_eps
        pos = self.layout(**options)
        [xmin, xmax, ymin, ymax] = self._layout_bounding_box(pos)
        for v in pos:
            pos[v] = (1.8*(pos[v][0] - xmin)/(xmax - xmin) - 0.9, 1.8*(pos[v][1] - ymin)/(ymax - ymin) - 0.9)
        if filename[-4:] != '.eps':
            filename += '.eps'
        f = open(filename, 'w')
        f.write( print_graph_eps(self.vertices(), self.edge_iterator(), pos) )
        f.close()

    @doc_index("Algorithmically hard stuff")
    def topological_minor(self, H, vertices=False, paths=False, solver=None, verbose=0,
                          *, integrality_tolerance=1e-3):
        r"""
        Return a topological `H`-minor from ``self`` if one exists.

        We say that a graph `G` has a topological `H`-minor (or that it has a
        graph isomorphic to `H` as a topological minor), if `G` contains a
        subdivision of a graph isomorphic to `H` (i.e.  obtained from `H`
        through arbitrary subdivision of its edges) as a subgraph.

        For more information, see the :wikipedia:`Minor_(graph_theory)`.

        INPUT:

        - ``H`` -- The topological minor to find in the current graph.

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        OUTPUT:

        The topological `H`-minor found is returned as a subgraph `M` of
        ``self``, such that the vertex `v` of `M` that represents a vertex `h\in
        H` has ``h`` as a label (see :meth:`get_vertex
        <sage.graphs.generic_graph.GenericGraph.get_vertex>` and
        :meth:`set_vertex <sage.graphs.generic_graph.GenericGraph.set_vertex>`),
        and such that every edge of `M` has as a label the edge of `H` it
        (partially) represents.

        If no topological minor is found, this method returns ``False``.

        ALGORITHM:

        Mixed Integer Linear Programming.

        COMPLEXITY:

        Theoretically, when `H` is fixed, testing for the existence of a
        topological `H`-minor is polynomial. The known algorithms are highly
        exponential in `H`, though.

        .. NOTE::

            This function can be expected to be *very* slow, especially where
            the topological minor does not exist.

            (CPLEX seems to be *much* more efficient than GLPK on this kind of
            problem)

        EXAMPLES:

        Petersen's graph has a topological `K_4`-minor::

            sage: g = graphs.PetersenGraph()
            sage: g.topological_minor(graphs.CompleteGraph(4))
            Subgraph of (Petersen graph): Graph on ...

        And a topological `K_{3,3}`-minor::

            sage: g.topological_minor(graphs.CompleteBipartiteGraph(3,3))
            Subgraph of (Petersen graph): Graph on ...

        And of course, a tree has no topological `C_3`-minor::

            sage: g = graphs.RandomGNP(15,.3)
            sage: g = g.subgraph(edges = g.min_spanning_tree())
            sage: g.topological_minor(graphs.CycleGraph(3))
            False
        """
        self._scream_if_not_simple()
        H._scream_if_not_simple()
        # Useful alias ...
        G = self

        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
        p = MixedIntegerLinearProgram(solver=solver)

        # This is an existence problem
        p.set_objective(None)

        #######################
        # Vertex representative #
        #######################
        #
        # v_repr[h,g] = 1 if vertex h from H is represented by vertex
        # g from G, 0 otherwise

        v_repr = p.new_variable(binary=True)

        # Exactly one representative per vertex of H
        for h in H:
            p.add_constraint(p.sum(v_repr[h,g] for g in G), min=1, max=1)

        # A vertex of G can only represent one vertex of H
        for g in G:
            p.add_constraint(p.sum(v_repr[h,g] for h in H), max=1)

        ###################
        # Is representent #
        ###################
        #
        # is_repr[v] = 1 if v represents some vertex of H

        is_repr = p.new_variable(binary=True)

        for g in G:
            for h in H:
                p.add_constraint(v_repr[h,g] - is_repr[g], max=0)

        ###################################
        # paths between the representents #
        ###################################
        #
        # For any edge (h1,h2) in H, we have a corresponding path in G
        # between the representatives of h1 and h2. Which means there is
        # a flow of intensity 1 from one to the other.
        # We are then writing a flow problem for each edge of H.
        #
        # The variable flow[(h1,h2),(g1,g2)] indicates the amount of
        # flow on the edge (g1,g2) representing the edge (h1,h2).

        flow = p.new_variable(binary=True)

        # These functions return the balance of flow corresponding to
        # commodity C at vertex v
        def flow_in(C, v):
            return p.sum(flow[C,(v,u)] for u in G.neighbor_iterator(v))

        def flow_out(C, v):
            return p.sum(flow[C,(u,v)] for u in G.neighbor_iterator(v))

        def flow_balance(C, v):
            return flow_in(C,v) - flow_out(C,v)

        for h1,h2 in H.edge_iterator(labels=False):

            for v in G:

                # The flow balance depends on whether the vertex v is a
                # representative of h1 or h2 in G, or a representative of none
                p.add_constraint(flow_balance((h1,h2),v) == v_repr[h1,v] - v_repr[h2,v])

        #############################
        # Internal vertex of a path #
        #############################
        #
        # is_internal[C][g] = 1 if a vertex v from G is located on the
        # path representing the edge (=commodity) C

        is_internal = p.new_variable(binary=True)

        # When is a vertex internal for a commodity ?
        for C in H.edge_iterator(labels=False):
            for g in G:
                p.add_constraint(flow_in(C,g) + flow_out(C,g) - is_internal[C,g], max=1)

        ############################
        # Two paths do not cross ! #
        ############################

        # A vertex can only be internal for one commodity, and zero if
        # the vertex is a representent

        for g in G:
            p.add_constraint(p.sum(is_internal[C,g] for C in H.edge_iterator(labels=False))
                              + is_repr[g], max=1)

        # (The following inequalities are not necessary, but they seem to be of
        # help (the solvers find the answer quicker when they are added)

        # The flow on one edge can go in only one direction. Besides, it can
        # belong to at most one commodity and has a maximum intensity of 1.

        for g1,g2 in G.edge_iterator(labels=None):

            p.add_constraint(   p.sum(flow[C,(g1,g2)] for C in H.edge_iterator(labels=False))
                              + p.sum(flow[C,(g2,g1)] for C in H.edge_iterator(labels=False)),
                                max=1)


        # Now we can solve the problem itself !

        try:
            p.solve(log=verbose)

        except MIPSolverException:
            return False


        minor = G.subgraph(immutable=False)

        is_repr = p.get_values(is_repr, convert=bool, tolerance=integrality_tolerance)
        v_repr = p.get_values(v_repr, convert=bool, tolerance=integrality_tolerance)
        flow = p.get_values(flow, convert=bool, tolerance=integrality_tolerance)

        for u,v in minor.edge_iterator(labels=False):
            used = False
            for C in H.edge_iterator(labels=False):

                if flow[C,(u,v)] or flow[C,(v,u)]:
                    used = True
                    minor.set_edge_label(u, v, C)
                    break
            if not used:
                minor.delete_edge(u, v)

        minor.delete_vertices(v for v in minor if minor.degree(v) == 0)

        for g in minor:
            if is_repr[g]:
                for h in H:
                    if v_repr[h,v]:
                        minor.set_vertex(g, h)
                        break

        return minor

    ### Cliques

    @doc_index("Clique-related methods")
    def cliques_maximal(self, algorithm="native"):
        """
        Return the list of all maximal cliques.

        Each clique is represented by a list of vertices. A clique is an induced
        complete subgraph, and a maximal clique is one not contained in a larger
        one.

        INPUT:

        - ``algorithm`` -- can be set to ``"native"`` (default) to use Sage's
          own implementation, or to ``"NetworkX"`` to use NetworkX'
          implementation of the Bron and Kerbosch Algorithm [BK1973]_.


        .. NOTE::

            This method sorts its output before returning it. If you prefer to
            save the extra time, you can call
            :class:`sage.graphs.independent_sets.IndependentSets` directly.

        .. NOTE::

            Sage's implementation of the enumeration of *maximal* independent
            sets is not much faster than NetworkX' (expect a 2x speedup), which
            is surprising as it is written in Cython. This being said, the
            algorithm from NetworkX appears to be slightly different from this
            one, and that would be a good thing to explore if one wants to
            improve the implementation.

        ALGORITHM:

        This function is based on NetworkX's implementation of the Bron and
        Kerbosch Algorithm [BK1973]_.

        EXAMPLES::

            sage: graphs.ChvatalGraph().cliques_maximal()
            [[0, 1], [0, 4], [0, 6], [0, 9], [1, 2], [1, 5], [1, 7], [2, 3],
             [2, 6], [2, 8], [3, 4], [3, 7], [3, 9], [4, 5], [4, 8], [5, 10],
             [5, 11], [6, 10], [6, 11], [7, 8], [7, 11], [8, 10], [9, 10], [9, 11]]
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2, 2])
            sage: G.cliques_maximal()
            [[0, 1, 2], [0, 1, 3]]
            sage: C = graphs.PetersenGraph()
            sage: C.cliques_maximal()
            [[0, 1], [0, 4], [0, 5], [1, 2], [1, 6], [2, 3], [2, 7], [3, 4],
             [3, 8], [4, 9], [5, 7], [5, 8], [6, 8], [6, 9], [7, 9]]
            sage: C = Graph('DJ{')
            sage: C.cliques_maximal()
            [[0, 4], [1, 2, 3, 4]]

        Comparing the two implementations::

            sage: g = graphs.RandomGNP(20,.7)
            sage: s1 = Set(map(Set, g.cliques_maximal(algorithm="NetworkX")))
            sage: s2 = Set(map(Set, g.cliques_maximal(algorithm="native")))
            sage: s1 == s2
            True
        """
        if algorithm == "native":
            from sage.graphs.independent_sets import IndependentSets
            return list(IndependentSets(self, maximal=True, complement=True))
        elif algorithm == "NetworkX":
            import networkx
            return list(networkx.find_cliques(self.networkx_graph()))
        else:
            raise ValueError("Algorithm must be equal to 'native' or to 'NetworkX'.")

    @doc_index("Clique-related methods")
    def clique_maximum(self,  algorithm="Cliquer", solver=None, verbose=0,
                       *, integrality_tolerance=1e-3):
        """
        Return the vertex set of a maximal order complete subgraph.

        INPUT:

        - ``algorithm`` -- the algorithm to be used :

          - If ``algorithm = "Cliquer"`` (default), wraps the C program
            Cliquer [NO2003]_.

          - If ``algorithm = "MILP"``, the problem is solved through a Mixed
            Integer Linear Program.

            (see :class:`~sage.numerical.mip.MixedIntegerLinearProgram`)

          - If ``algorithm = "mcqd"``, uses the MCQD solver
            (`<http://www.sicmm.org/~konc/maxclique/>`_). Note that the MCQD
            package must be installed.

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        Parameters ``solver`` and ``verbose`` are used only when
        ``algorithm="MILP"``.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        ALGORITHM:

        This function is based on Cliquer [NO2003]_.

        EXAMPLES:

        Using Cliquer (default)::

            sage: C = graphs.PetersenGraph()
            sage: C.clique_maximum()
            [7, 9]
            sage: C = Graph('DJ{')
            sage: C.clique_maximum()
            [1, 2, 3, 4]

        Through a Linear Program::

            sage: len(C.clique_maximum(algorithm="MILP"))
            4

        TESTS:

        Wrong algorithm::

            sage: C.clique_maximum(algorithm="BFS")
            Traceback (most recent call last):
            ...
            NotImplementedError: Only 'MILP', 'Cliquer' and 'mcqd' are supported.

        """
        self._scream_if_not_simple(allow_multiple_edges=True)
        if algorithm == "Cliquer":
            from sage.graphs.cliquer import max_clique
            return max_clique(self)
        elif algorithm == "MILP":
            return self.complement().independent_set(algorithm=algorithm, solver=solver, verbose=verbose,
                                                     integrality_tolerance=integrality_tolerance)
        elif algorithm == "mcqd":
            return mcqd(self)
        else:
            raise NotImplementedError("Only 'MILP', 'Cliquer' and 'mcqd' are supported.")

    @doc_index("Clique-related methods")
    def clique_number(self, algorithm="Cliquer", cliques=None, solver=None, verbose=0,
                      *, integrality_tolerance=1e-3):
        r"""
        Return the order of the largest clique of the graph

        This is also called as the clique number.

        .. NOTE::

            Currently only implemented for undirected graphs. Use ``to_undirected``
            to convert a digraph to an undirected graph.

        INPUT:

        - ``algorithm`` -- the algorithm to be used :

          - If ``algorithm = "Cliquer"``, wraps the C program Cliquer
            [NO2003]_.

          - If ``algorithm = "networkx"``, uses the NetworkX's implementation of
            the Bron and Kerbosch Algorithm [BK1973]_.

          - If ``algorithm = "MILP"``, the problem is solved through a Mixed
            Integer Linear Program.

            (see :class:`~sage.numerical.mip.MixedIntegerLinearProgram`)

          - If ``algorithm = "mcqd"``, uses the MCQD solver
            (`<http://insilab.org/maxclique/>`_). Note that the MCQD
            package must be installed.

        - ``cliques`` -- an optional list of cliques that can be input if
          already computed. Ignored unless ``algorithm=="networkx"``.

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        ALGORITHM:

        This function is based on Cliquer [NO2003]_ and [BK1973]_.

        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.clique_number()
            4
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.clique_number()
            3

        By definition the clique number of a complete graph is its order::

            sage: all(graphs.CompleteGraph(i).clique_number() == i for i in range(1,15))
            True

        A non-empty graph without edges has a clique number of 1::

            sage: all((i*graphs.CompleteGraph(1)).clique_number() == 1 for i in range(1,15))
            True

        A complete multipartite graph with k parts has clique number k::

            sage: all((i*graphs.CompleteMultipartiteGraph(i*[5])).clique_number() == i for i in range(1,6))
            True

        TESTS::

            sage: g = graphs.PetersenGraph()
            sage: g.clique_number(algorithm="MILP")
            2
            sage: for i in range(10):                                            # optional - mcqd
            ....:     g = graphs.RandomGNP(15,.5)                                # optional - mcqd
            ....:     if g.clique_number() != g.clique_number(algorithm="mcqd"): # optional - mcqd
            ....:         print("This is dead wrong !")                          # optional - mcqd
        """
        self._scream_if_not_simple(allow_loops=False)
        if algorithm == "Cliquer":
            from sage.graphs.cliquer import clique_number
            return clique_number(self)
        elif algorithm == "networkx":
            import networkx
            return networkx.graph_clique_number(self.networkx_graph(), cliques)
        elif algorithm == "MILP":
            return len(self.complement().independent_set(algorithm=algorithm, solver=solver, verbose=verbose,
                                                         integrality_tolerance=integrality_tolerance))
        elif algorithm == "mcqd":
            return len(mcqd(self))
        else:
            raise NotImplementedError("Only 'networkx' 'MILP' 'Cliquer' and 'mcqd' are supported.")

    @doc_index("Clique-related methods")
    def cliques_number_of(self, vertices=None, cliques=None):
        """
        Return a dictionary of the number of maximal cliques containing each
        vertex, keyed by vertex.

        This returns a single value if only one input vertex.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        INPUT:

        - ``vertices`` -- the vertices to inspect (default is entire graph)

        - ``cliques`` -- list of cliques (if already computed)


        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.cliques_number_of()
            {0: 1, 1: 1, 2: 1, 3: 1, 4: 2}
            sage: E = C.cliques_maximal()
            sage: E
            [[0, 4], [1, 2, 3, 4]]
            sage: C.cliques_number_of(cliques=E)
            {0: 1, 1: 1, 2: 1, 3: 1, 4: 2}
            sage: F = graphs.Grid2dGraph(2,3)
            sage: F.cliques_number_of()
            {(0, 0): 2, (0, 1): 3, (0, 2): 2, (1, 0): 2, (1, 1): 3, (1, 2): 2}
            sage: F.cliques_number_of(vertices=[(0, 1), (1, 2)])
            {(0, 1): 3, (1, 2): 2}
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_number_of()
            {0: 2, 1: 2, 2: 1, 3: 1}
        """
        import networkx
        return networkx.number_of_cliques(self.networkx_graph(), vertices, cliques)

    @doc_index("Clique-related methods")
    def cliques_get_max_clique_graph(self):
        """
        Return the clique graph.

        Vertices of the result are the maximal cliques of the graph, and edges
        of the result are between maximal cliques with common members in the
        original graph.

        For more information, see the :wikipedia:`Clique_graph`.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        EXAMPLES::

            sage: (graphs.ChvatalGraph()).cliques_get_max_clique_graph()
            Graph on 24 vertices
            sage: ((graphs.ChvatalGraph()).cliques_get_max_clique_graph()).show(figsize=[2,2], vertex_size=20, vertex_labels=False)
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_get_max_clique_graph()
            Graph on 2 vertices
            sage: (G.cliques_get_max_clique_graph()).show(figsize=[2,2])
        """
        import networkx
        return Graph(networkx.make_max_clique_graph(self.networkx_graph(), create_using=networkx.MultiGraph()),
                     multiedges=False)

    @doc_index("Clique-related methods")
    def cliques_get_clique_bipartite(self, **kwds):
        """
        Return a bipartite graph constructed such that maximal cliques are the
        right vertices and the left vertices are retained from the given
        graph. Right and left vertices are connected if the bottom vertex
        belongs to the clique represented by a top vertex.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        EXAMPLES::

            sage: (graphs.ChvatalGraph()).cliques_get_clique_bipartite()
            Bipartite graph on 36 vertices
            sage: ((graphs.ChvatalGraph()).cliques_get_clique_bipartite()).show(figsize=[2,2], vertex_size=20, vertex_labels=False)
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_get_clique_bipartite()
            Bipartite graph on 6 vertices
            sage: (G.cliques_get_clique_bipartite()).show(figsize=[2,2])
        """
        from .bipartite_graph import BipartiteGraph
        import networkx
        return BipartiteGraph(networkx.make_clique_bipartite(self.networkx_graph(), **kwds))

    @doc_index("Algorithmically hard stuff")
    @rename_keyword(deprecation=32238, verbosity='verbose')
    def independent_set(self, algorithm="Cliquer", value_only=False, reduction_rules=True,
                        solver=None, verbose=0, *, integrality_tolerance=1e-3):
        r"""
        Return a maximum independent set.

        An independent set of a graph is a set of pairwise non-adjacent
        vertices. A maximum independent set is an independent set of maximum
        cardinality.  It induces an empty subgraph.

        Equivalently, an independent set is defined as the complement of a
        vertex cover.

        For more information, see the
        :wikipedia:`Independent_set_(graph_theory)` and the
        :wikipedia:`Vertex_cover`.

        INPUT:

        - ``algorithm`` -- the algorithm to be used

          * If ``algorithm = "Cliquer"`` (default), the problem is solved
            using Cliquer [NO2003]_.

            (see the :mod:`Cliquer modules <sage.graphs.cliquer>`)

          * If ``algorithm = "MILP"``, the problem is solved through a Mixed
            Integer Linear Program.

            (see :class:`~sage.numerical.mip.MixedIntegerLinearProgram`)

         * If ``algorithm = "mcqd"``, uses the MCQD solver
           (`<http://www.sicmm.org/~konc/maxclique/>`_). Note that the MCQD
           package must be installed.

        - ``value_only`` -- boolean (default: ``False``); if set to ``True``,
          only the size of a maximum independent set is returned. Otherwise,
          a maximum independent set is returned as a list of vertices.

        - ``reduction_rules`` -- (default: ``True``); specify if the reductions
          rules from kernelization must be applied as pre-processing or not.
          See [ACFLSS04]_ for more details. Note that depending on the instance,
          it might be faster to disable reduction rules.

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        .. NOTE::

            While Cliquer/MCAD are usually (and by far) the most efficient
            implementations, the MILP formulation sometimes proves faster on
            very "symmetrical" graphs.

        EXAMPLES:

        Using Cliquer::

            sage: C = graphs.PetersenGraph()
            sage: C.independent_set()
            [0, 3, 6, 7]

        As a linear program::

            sage: C = graphs.PetersenGraph()
            sage: len(C.independent_set(algorithm="MILP"))
            4

        .. PLOT::

            g = graphs.PetersenGraph()
            sphinx_plot(g.plot(partition=[g.independent_set()]))
        """
        my_cover = self.vertex_cover(algorithm=algorithm, value_only=value_only,
                                         reduction_rules=reduction_rules,
                                         solver=solver, verbose=verbose,
                                         integrality_tolerance=integrality_tolerance)
        if value_only:
            return self.order() - my_cover
        else:
            my_cover = set(my_cover)
            return [u for u in self if u not in my_cover]

    @doc_index("Algorithmically hard stuff")
    @rename_keyword(deprecation=32238, verbosity='verbose')
    def vertex_cover(self, algorithm="Cliquer", value_only=False,
                     reduction_rules=True, solver=None, verbose=0,
                     *, integrality_tolerance=1e-3):
        r"""
        Return a minimum vertex cover of self represented by a set of vertices.

        A minimum vertex cover of a graph is a set `S` of vertices such that
        each edge is incident to at least one element of `S`, and such that `S`
        is of minimum cardinality. For more information, see the
        :wikipedia:`Vertex_cover`.

        Equivalently, a vertex cover is defined as the complement of an
        independent set.

        As an optimization problem, it can be expressed as follows:

        .. MATH::

            \mbox{Minimize : }&\sum_{v\in G} b_v\\
            \mbox{Such that : }&\forall (u,v) \in G.edges(), b_u+b_v\geq 1\\
            &\forall x\in G, b_x\mbox{ is a binary variable}

        INPUT:

        - ``algorithm`` -- string (default: ``"Cliquer"``). Indicating which
          algorithm to use. It can be one of those values.

          - ``"Cliquer"`` will compute a minimum vertex cover using the Cliquer
            package.

          - ``"MILP"`` will compute a minimum vertex cover through a mixed
            integer linear program.

          - ``"mcqd"`` will use the MCQD solver
            (`<http://www.sicmm.org/~konc/maxclique/>`_). Note that the MCQD
            package must be installed.

        - ``value_only`` -- boolean (default: ``False``); if set to ``True``,
          only the size of a minimum vertex cover is returned. Otherwise,
          a minimum vertex cover is returned as a list of vertices.

        - ``reduction_rules`` -- (default: ``True``); specify if the reductions
          rules from kernelization must be applied as pre-processing or not.
          See [ACFLSS04]_ for more details. Note that depending on the instance,
          it might be faster to disable reduction rules.

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        EXAMPLES:

        On the Pappus graph::

           sage: g = graphs.PappusGraph()
           sage: g.vertex_cover(value_only=True)
           9

        .. PLOT::

            g = graphs.PappusGraph()
            sphinx_plot(g.plot(partition=[g.vertex_cover()]))

        TESTS:

        The two algorithms should return the same result::

           sage: g = graphs.RandomGNP(10, .5)
           sage: vc1 = g.vertex_cover(algorithm="MILP")
           sage: vc2 = g.vertex_cover(algorithm="Cliquer")
           sage: len(vc1) == len(vc2)
           True

        The cardinality of the vertex cover is unchanged when reduction rules
        are used. First for trees::

           sage: for i in range(20):
           ....:     g = graphs.RandomTree(20)
           ....:     vc1_set = g.vertex_cover()
           ....:     vc1 = len(vc1_set)
           ....:     vc2 = g.vertex_cover(value_only=True, reduction_rules=False)
           ....:     if vc1 != vc2:
           ....:         print("Error :", vc1, vc2)
           ....:         print("With reduction rules :", vc1)
           ....:         print("Without reduction rules :", vc2)
           ....:         break
           ....:     g.delete_vertices(vc1_set)
           ....:     if g.size():
           ....:         print("This thing is not a vertex cover !")

        Then for random GNP graphs::

           sage: for i in range(20):
           ....:     g = graphs.RandomGNP(50, 0.08)
           ....:     vc1_set = g.vertex_cover()
           ....:     vc1 = len(vc1_set)
           ....:     vc2 = g.vertex_cover(value_only=True, reduction_rules=False)
           ....:     if vc1 != vc2:
           ....:         print("Error :", vc1, vc2)
           ....:         print("With reduction rules :", vc1)
           ....:         print("Without reduction rules :", vc2)
           ....:         break
           ....:     g.delete_vertices(vc1_set)
           ....:     if g.size():
           ....:         print("This thing is not a vertex cover !")

        Testing mcqd::

            sage: graphs.PetersenGraph().vertex_cover(algorithm="mcqd", value_only=True) # optional - mcqd
            6

        Given a wrong algorithm::

            sage: graphs.PetersenGraph().vertex_cover(algorithm="guess")
            Traceback (most recent call last):
            ...
            ValueError: the algorithm must be "Cliquer", "MILP" or "mcqd"

        Ticket :trac:`24287` is fixed::

            sage: G = Graph([(0,1)]*5 + [(1,2)]*2, multiedges=True)
            sage: G.vertex_cover(reduction_rules=True, algorithm='MILP')
            [1]
            sage: G.vertex_cover(reduction_rules=False)
            [1]

        Ticket :trac:`25988` is fixed::

            sage: B = BipartiteGraph(graphs.CycleGraph(6))
            sage: B.vertex_cover(algorithm='Cliquer', reduction_rules=True)
            [1, 3, 5]
        """
        self._scream_if_not_simple(allow_multiple_edges=True)
        g = self

        ppset = []
        folded_vertices = []

        ###################
        # Reduction rules #
        ###################

        if reduction_rules:
            # We apply simple reduction rules allowing to identify vertices that
            # belongs to an optimal vertex cover

            # We first take a copy of the graph without multiple edges, if any.
            g = Graph(data=self.edges(), format='list_of_edges',
                          multiedges=self.allows_multiple_edges())
            g.allow_multiple_edges(False)

            degree_at_most_two = {u for u in g if g.degree(u) <= 2}

            while degree_at_most_two:

                u = degree_at_most_two.pop()
                du = g.degree(u)

                if not du:
                    # RULE 1: isolated vertices are not part of the cover. We
                    # simply remove them from the graph. The degree of such
                    # vertices may have been reduced to 0 while applying other
                    # reduction rules
                    g.delete_vertex(u)

                elif du == 1:
                    # RULE 2: If a vertex u has degree 1, we select its neighbor
                    # v and remove both u and v from g.
                    v = next(g.neighbor_iterator(u))
                    ppset.append(v)
                    g.delete_vertex(u)

                    for w in g.neighbor_iterator(v):
                        if g.degree(w) <= 3:
                            # The degree of w will be at most two after the
                            # deletion of v
                            degree_at_most_two.add(w)

                    g.delete_vertex(v)
                    degree_at_most_two.discard(v)

                elif du == 2:
                    v,w  = g.neighbors(u)

                    if g.has_edge(v, w):
                        # RULE 3: If the neighbors v and w of a degree 2 vertex
                        # u are incident, then we select both v and w and remove
                        # u, v, and w from g.
                        ppset.append(v)
                        ppset.append(w)
                        g.delete_vertex(u)
                        neigh = set(g.neighbors(v) + g.neighbors(w)).difference([v, w])
                        g.delete_vertex(v)
                        g.delete_vertex(w)

                        for z in neigh:
                            if g.degree(z) <= 2:
                                degree_at_most_two.add(z)

                    else:
                        # RULE 4, folded vertices: If the neighbors v and w of a
                        # degree 2 vertex u are not incident, then we contract
                        # edges (u, v), (u, w). Then, if the solution contains u,
                        # we replace it with v and w. Otherwise, we let u in the
                        # solution.
                        neigh = set(g.neighbors(v) + g.neighbors(w)).difference([u, v, w])
                        g.delete_vertex(v)
                        g.delete_vertex(w)
                        for z in neigh:
                            g.add_edge(u,z)

                        folded_vertices.append((u, v, w))

                        if g.degree(u) <= 2:
                            degree_at_most_two.add(u)

                    degree_at_most_two.discard(v)
                    degree_at_most_two.discard(w)


                # RULE 5:
                # TODO: add extra reduction rules


        ##################
        # Main Algorithm #
        ##################

        if not g.order():
            # Reduction rules were sufficients to get the solution
            size_cover_g = 0
            cover_g = set()

        elif algorithm == "Cliquer" or algorithm == "mcqd":
            if g.has_multiple_edges() and not reduction_rules:
                g = copy(g)
                g.allow_multiple_edges(False)

            independent = g.complement().clique_maximum(algorithm=algorithm)
            if value_only:
                size_cover_g = g.order() - len(independent)
            else:
                cover_g = set(uu for uu in g if uu not in independent)

        elif algorithm == "MILP":

            from sage.numerical.mip import MixedIntegerLinearProgram
            p = MixedIntegerLinearProgram(maximization=False, solver=solver)
            b = p.new_variable(binary=True)

            # minimizes the number of vertices in the set
            p.set_objective(p.sum(b[v] for v in g))

            # an edge contains at least one vertex of the minimum vertex cover
            for u,v in g.edge_iterator(labels=None):
                p.add_constraint(b[u] + b[v], min=1)

            p.solve(log=verbose)
            b = p.get_values(b, convert=bool, tolerance=integrality_tolerance)
            if value_only:
                size_cover_g = sum(1 for v in g if b[v])
            else:
                cover_g = set(v for v in g if b[v])
        else:
            raise ValueError('the algorithm must be "Cliquer", "MILP" or "mcqd"')

        #########################
        # Returning the results #
        #########################

        # We finally reconstruct the solution according the reduction rules
        if value_only:
            return len(ppset) + len(folded_vertices) + size_cover_g
        else:
            # RULES 2 and 3:
            cover_g.update(ppset)
            # RULE 4:
            folded_vertices.reverse()
            for u,v,w in folded_vertices:
                if u in cover_g:
                    cover_g.discard(u)
                    cover_g.add(v)
                    cover_g.add(w)
                else:
                    cover_g.add(u)
            return list(cover_g)

    @doc_index("Connectivity, orientations, trees")
    def ear_decomposition(self):
        r"""
        Return an Ear decomposition of the graph.

        An ear of an undirected graph `G` is a path `P` where the two endpoints
        of the path may coincide (i.e., form a cycle), but where otherwise no
        repetition of edges or vertices is allowed, so every internal vertex of
        `P` has degree two in `P`.

        An ear decomposition of an undirected graph `G` is a partition of its
        set of edges into a sequence of ears, such that the one or two endpoints
        of each ear belong to earlier ears in the sequence and such that the
        internal vertices of each ear do not belong to any earlier ear.

        For more information, see the :wikipedia:`Ear_decomposition`.

        This method implements the linear time algorithm presented in
        [Sch2013]_.

        OUTPUT:

        - A nested list representing the cycles and chains of the ear
          decomposition of the graph.

        EXAMPLES:

        Ear decomposition of an outer planar graph of order 13::

            sage: g = Graph('LlCG{O@?GBOMW?')
            sage: g.ear_decomposition()
            [[0, 3, 2, 1, 0],
             [0, 7, 4, 3],
             [0, 11, 9, 8, 7],
             [1, 12, 2],
             [3, 6, 5, 4],
             [4, 6],
             [7, 10, 8],
             [7, 11],
             [8, 11]]

        Ear decomposition of a biconnected graph::

            sage: g = graphs.CycleGraph(4)
            sage: g.ear_decomposition()
            [[0, 3, 2, 1, 0]]

        Ear decomposition of a connected but not biconnected graph::

            sage: G = Graph()
            sage: G.add_cycle([0,1,2])
            sage: G.add_edge(0,3)
            sage: G.add_cycle([3,4,5,6])
            sage: G.ear_decomposition()
            [[0, 2, 1, 0], [3, 6, 5, 4, 3]]

        The ear decomposition of a multigraph with loops is the same as the ear
        decomposition of the underlying simple graph::

            sage: g = graphs.BullGraph()
            sage: g.allow_multiple_edges(True)
            sage: g.add_edges(g.edges())
            sage: g.allow_loops(True)
            sage: u = g.random_vertex()
            sage: g.add_edge(u, u)
            sage: g
            Bull graph: Looped multi-graph on 5 vertices
            sage: h = g.to_simple()
            sage: g.ear_decomposition() == h.ear_decomposition()
            True

        TESTS::

            sage: g=Graph()
            sage: g
            Graph on 0 vertices
            sage: g.ear_decomposition()
            Traceback (most recent call last):
            ...
            ValueError: ear decomposition is defined for graphs of order at least 3

        """
        # Ear decomposition of a graph of order < 3 is [].
        if self.order() < 3:
            raise ValueError("ear decomposition is defined for graphs of order at least 3")

        # List to store the order in which dfs visits vertices.
        dfs_order = []

        # Boolean dict to mark vertices as visited or unvisited during
        # Dfs traversal in graph.
        seen = set()

        # Boolean dict to mark vertices as visited or unvisited in
        # Dfs tree traversal.
        traversed = set()

        # Dictionary to store parent vertex of all the visited vertices.
        # Initialized for the first vertex to be visited.
        parent = {next(self.vertex_iterator()): None}

        # List to store visit_time of vertices in Dfs traversal.
        value = {}

        # List to store all the chains and cycles of the input graph G.
        chains = []

        # DFS() : Function that performs depth first search on input graph G and
        #         stores DFS tree in parent array format.
        def DFS(v):
            """
            Depth first search step from vertex v.
            """
            # make v are visited, update its time of visited and value
            seen.add(v)
            dfs_order.append(v)

            # Traverse though all the neighbor vertices of v
            for u in self.neighbor_iterator(v):
                # if any neighbor is not visited, enter
                if u not in seen:
                    # Set the parent of u in DFS tree as v and continue
                    # exploration
                    parent[u] = v
                    DFS(u)

        # Traverse() : Function that use G-T (non-tree edges) to find cycles
        #              and chains by traversing in DFS tree.
        def traverse(start, pointer):
            # Make the first end of non-tree edge visited
            traversed.add(start)
            chain = [start]

            # Traverse DFS Tree of G and print all the not visited vertices
            # Appending all the vertices in chain
            while True:
                chain.append(pointer)
                if pointer in traversed:
                    break
                traversed.add(pointer)
                pointer = parent[pointer]
            chains.append(chain)

        # Perform ear decomposition on each connected component of input graph.
        for v in self:
            if v not in seen:
              # Start the depth first search from first vertex
                DFS(v)
                value = {u:i for i,u in enumerate(dfs_order)}

                # Traverse all the non Tree edges, according to DFS order
                for u in dfs_order:
                    for neighbor in self.neighbor_iterator(u):
                        if value[u] < value[neighbor] and u != parent[neighbor]:
                            traverse(u, neighbor)

                dfs_order = []

        return chains

    @doc_index("Clique-related methods")
    def cliques_vertex_clique_number(self, algorithm="cliquer", vertices=None,
                                     cliques=None):
        """
        Return a dictionary of sizes of the largest maximal cliques containing
        each vertex, keyed by vertex.

        Returns a single value if only one input vertex.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        INPUT:

         - ``algorithm`` -- either ``cliquer`` or ``networkx``

           - ``cliquer`` -- This wraps the C program Cliquer [NO2003]_.

           - ``networkx`` -- This function is based on NetworkX's implementation
             of the Bron and Kerbosch Algorithm [BK1973]_.

        - ``vertices`` -- the vertices to inspect (default is entire graph).
          Ignored unless ``algorithm=='networkx'``.

        - ``cliques`` -- list of cliques (if already computed).  Ignored unless
          ``algorithm=='networkx'``.

        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.cliques_vertex_clique_number()
            {0: 2, 1: 4, 2: 4, 3: 4, 4: 4}
            sage: E = C.cliques_maximal()
            sage: E
            [[0, 4], [1, 2, 3, 4]]
            sage: C.cliques_vertex_clique_number(cliques=E,algorithm="networkx")
            {0: 2, 1: 4, 2: 4, 3: 4, 4: 4}
            sage: F = graphs.Grid2dGraph(2,3)
            sage: F.cliques_vertex_clique_number(algorithm="networkx")
            {(0, 0): 2, (0, 1): 2, (0, 2): 2, (1, 0): 2, (1, 1): 2, (1, 2): 2}
            sage: F.cliques_vertex_clique_number(vertices=[(0, 1), (1, 2)])
            {(0, 1): 2, (1, 2): 2}
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_vertex_clique_number()
            {0: 3, 1: 3, 2: 3, 3: 3}
        """
        if algorithm == "cliquer":
            from sage.graphs.cliquer import clique_number
            if vertices is None:
                vertices = self
            value = {}
            for v in vertices:
                value[v] = 1 + clique_number(self.subgraph(self.neighbors(v)))
                self.subgraph(self.neighbors(v)).plot()
            return value
        elif algorithm == "networkx":
            import networkx
            return networkx.node_clique_number(self.networkx_graph(), vertices, cliques)
        else:
            raise NotImplementedError("Only 'networkx' and 'cliquer' are supported.")

    @doc_index("Clique-related methods")
    def cliques_containing_vertex(self, vertices=None, cliques=None):
        """
        Return the cliques containing each vertex, represented as a dictionary
        of lists of lists, keyed by vertex.

        Returns a single list if only one input vertex.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        INPUT:

        - ``vertices`` -- the vertices to inspect (default is entire graph)

        - ``cliques`` -- list of cliques (if already computed)

        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.cliques_containing_vertex()
            {0: [[4, 0]], 1: [[4, 1, 2, 3]], 2: [[4, 1, 2, 3]], 3: [[4, 1, 2, 3]], 4: [[4, 0], [4, 1, 2, 3]]}
            sage: E = C.cliques_maximal()
            sage: E
            [[0, 4], [1, 2, 3, 4]]
            sage: C.cliques_containing_vertex(cliques=E)
            {0: [[0, 4]], 1: [[1, 2, 3, 4]], 2: [[1, 2, 3, 4]], 3: [[1, 2, 3, 4]], 4: [[0, 4], [1, 2, 3, 4]]}

            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_containing_vertex()
            {0: [[0, 1, 2], [0, 1, 3]], 1: [[0, 1, 2], [0, 1, 3]], 2: [[0, 1, 2]], 3: [[0, 1, 3]]}

        Since each clique of a 2 dimensional grid corresponds to an edge, the
        number of cliques in which a vertex is involved equals its degree::

            sage: F = graphs.Grid2dGraph(2,3)
            sage: d = F.cliques_containing_vertex()
            sage: all(F.degree(u) == len(cliques) for u,cliques in d.items())
            True
            sage: d = F.cliques_containing_vertex(vertices=[(0, 1)])
            sage: list(d)
            [(0, 1)]
            sage: sorted(sorted(x for x in L) for L in d[(0, 1)])
            [[(0, 0), (0, 1)], [(0, 1), (0, 2)], [(0, 1), (1, 1)]]
        """
        import networkx
        return networkx.cliques_containing_node(self.networkx_graph(), vertices, cliques)

    @doc_index("Clique-related methods")
    def clique_complex(self):
        """
        Return the clique complex of self.

        This is the largest simplicial complex on the vertices of self whose
        1-skeleton is self.

        This is only makes sense for undirected simple graphs.

        EXAMPLES::

            sage: g = Graph({0:[1,2],1:[2],4:[]})
            sage: g.clique_complex()
            Simplicial complex with vertex set (0, 1, 2, 4) and facets {(4,), (0, 1, 2)}

            sage: h = Graph({0:[1,2,3,4],1:[2,3,4],2:[3]})
            sage: x = h.clique_complex()
            sage: x
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {(0, 1, 4), (0, 1, 2, 3)}
            sage: i = x.graph()
            sage: i==h
            True
            sage: x==i.clique_complex()
            True

        """
        if self.is_directed() or self.has_loops() or self.has_multiple_edges():
            raise ValueError("Self must be an undirected simple graph to have a clique complex.")
        import sage.topology.simplicial_complex
        C = sage.topology.simplicial_complex.SimplicialComplex(self.cliques_maximal(), maximality_check=True)
        C._graph = self
        return C

    @doc_index("Clique-related methods")
    def clique_polynomial(self, t=None):
        r"""
        Return the clique polynomial of self.

        This is the polynomial where the coefficient of `t^n` is the number of
        cliques in the graph with `n` vertices. The constant term of the clique
        polynomial is always taken to be one.

        EXAMPLES::

            sage: g = Graph()
            sage: g.clique_polynomial()
            1
            sage: g = Graph({0:[1]})
            sage: g.clique_polynomial()
            t^2 + 2*t + 1
            sage: g = graphs.CycleGraph(4)
            sage: g.clique_polynomial()
            4*t^2 + 4*t + 1

        """
        if t is None:
            R = PolynomialRing(ZZ, 't')
            t = R.gen()
        number_of = [0]*(self.order() + 1)
        for x in IndependentSets(self, complement=True):
            number_of[len(x)] += 1
        return sum(coeff*t**i for i,coeff in enumerate(number_of) if coeff)

    ### Miscellaneous

    @doc_index("Leftovers")
    def cores(self, k=None, with_labels=False):
        r"""
        Return the core number for each vertex in an ordered list.

        (for homomorphisms cores, see the :meth:`Graph.has_homomorphism_to`
        method)

        DEFINITIONS:

        * *K-cores* in graph theory were introduced by Seidman in 1983 and by
          Bollobas in 1984 as a method of (destructively) simplifying graph
          topology to aid in analysis and visualization. They have been more
          recently defined as the following by Batagelj et al:

          *Given a graph `G` with vertices set `V` and edges set `E`, the
          `k`-core of `G` is the graph obtained from `G` by recursively removing
          the vertices with degree less than `k`, for as long as there are any.*

          This operation can be useful to filter or to study some properties of
          the graphs. For instance, when you compute the 2-core of graph G, you
          are cutting all the vertices which are in a tree part of graph.  (A
          tree is a graph with no loops). See the :wikipedia:`K-core`.

          [PSW1996]_ defines a `k`-core of `G` as the largest subgraph (it is
          unique) of `G` with minimum degree at least `k`.

        * Core number of a vertex

          The core number of a vertex `v` is the largest integer `k` such that
          `v` belongs to the `k`-core of `G`.

        * Degeneracy

          The *degeneracy* of a graph `G`, usually denoted `\delta^*(G)`, is the
          smallest integer `k` such that the graph `G` can be reduced to the
          empty graph by iteratively removing vertices of degree `\leq k`.
          Equivalently, `\delta^*(G)=k` if `k` is the smallest integer such that
          the `k`-core of `G` is empty.

        IMPLEMENTATION:

        This implementation is based on the NetworkX implementation of the
        algorithm described in [BZ2003]_.

        INPUT:

        - ``k`` -- integer (default: ``None``);

            * If ``k = None`` (default), returns the core number for each vertex.

            * If ``k`` is an integer, returns a pair ``(ordering, core)``, where
              ``core`` is the list of vertices in the `k`-core of ``self``, and
              ``ordering`` is an elimination order for the other vertices such
              that each vertex is of degree strictly less than `k` when it is to
              be eliminated from the graph.

        - ``with_labels`` -- boolean (default: ``False``); when set to
          ``False``, and ``k = None``, the method returns a list whose `i` th
          element is the core number of the `i` th vertex. When set to ``True``,
          the method returns a dictionary whose keys are vertices, and whose
          values are the corresponding core numbers.

        .. SEEALSO::

           * Graph cores is also a notion related to graph homomorphisms. For
             this second meaning, see :meth:`Graph.has_homomorphism_to`.

        EXAMPLES::

            sage: (graphs.FruchtGraph()).cores()
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
            sage: (graphs.FruchtGraph()).cores(with_labels=True)
            {0: 3, 1: 3, 2: 3, 3: 3, 4: 3, 5: 3, 6: 3, 7: 3, 8: 3, 9: 3, 10: 3, 11: 3}
            sage: set_random_seed(0)
            sage: a = random_matrix(ZZ, 20, x=2, sparse=True, density=.1)
            sage: b = Graph(20)
            sage: b.add_edges(a.nonzero_positions(), loops=False)
            sage: cores = b.cores(with_labels=True); cores
            {0: 3, 1: 3, 2: 3, 3: 3, 4: 2, 5: 2, 6: 3, 7: 1, 8: 3, 9: 3, 10: 3, 11: 3, 12: 3, 13: 3, 14: 2, 15: 3, 16: 3, 17: 3, 18: 3, 19: 3}
            sage: [v for v,c in cores.items() if c >= 2] # the vertices in the 2-core
            [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

        Checking the 2-core of a random lobster is indeed the empty set::

            sage: g = graphs.RandomLobster(20, .5, .5)
            sage: ordering, core = g.cores(2)
            sage: len(core) == 0
            True
        """
        self._scream_if_not_simple()
        # compute the degrees of each vertex
        degrees = self.degree(labels=True)

        # Sort vertices by degree. Store in a list and keep track of where a
        # specific degree starts (effectively, the list is sorted by bins).
        verts = sorted(degrees.keys(), key=lambda x: degrees[x])
        bin_boundaries = [0]
        curr_degree = 0
        for i,v in enumerate(verts):
            if degrees[v] > curr_degree:
                bin_boundaries.extend([i] * (degrees[v] - curr_degree))
                curr_degree = degrees[v]
        vert_pos = {v: pos for pos,v in enumerate(verts)}
        # Set up initial guesses for core and lists of neighbors.
        core = degrees
        nbrs = {v: set(self.neighbors(v)) for v in self}
        # form vertex core building up from smallest
        for v in verts:

            # If all the vertices have a degree larger than k, we can return our
            # answer if k is not None
            if k is not None and core[v] >= k:
                return verts[:vert_pos[v]], verts[vert_pos[v]:]

            for u in nbrs[v]:
                if core[u] > core[v]:
                    nbrs[u].remove(v)

                    # Cleverly move u to the end of the next smallest bin (i.e.,
                    # subtract one from the degree of u). We do this by swapping
                    # u with the first vertex in the bin that contains u, then
                    # incrementing the bin boundary for the bin that contains u.
                    pos = vert_pos[u]
                    bin_start = bin_boundaries[core[u]]
                    vert_pos[u] = bin_start
                    vert_pos[verts[bin_start]] = pos
                    verts[bin_start],verts[pos] = verts[pos],verts[bin_start]
                    bin_boundaries[core[u]] += 1
                    core[u] -= 1

        if k is not None:
            return verts, []

        if with_labels:
            return core
        else:
            return list(core.values())

    @doc_index("Leftovers")
    def modular_decomposition(self, algorithm='habib', style='tuple'):
        r"""
        Return the modular decomposition of the current graph.

        A module of an undirected graph is a subset of vertices such that every
        vertex outside the module is either connected to all members of the
        module or to none of them. Every graph that has a nontrivial module can
        be partitioned into modules, and the increasingly fine partitions into
        modules form a tree. The ``modular_decomposition`` function returns
        that tree.

        INPUT:

        - ``algorithm`` -- string (default: ``'habib'``); specifies the
          algorithm to use among:

          - ``'tedder'`` -- linear time algorithm of [TCHP2008]_

          - ``'habib'`` -- `O(n^3)` algorithm of [HM1979]_. This algorithm is
            much simpler and so possibly less prone to errors.

        - ``style`` -- string (default: ``'tuple'``); specifies the output
          format:

          - ``'tuple'`` -- as nested tuples.

          - ``'tree'`` -- as :class:`~sage.combinat.rooted_tree.LabelledRootedTree`.

        OUTPUT:

        A pair of two values (recursively encoding the decomposition) :

        * The type of the current module :

          * ``"PARALLEL"``
          * ``"PRIME"``
          * ``"SERIES"``

        * The list of submodules (as list of pairs ``(type, list)``,
          recursively...) or the vertex's name if the module is a singleton.

        Crash course on modular decomposition:

        A module `M` of a graph `G` is a proper subset of its vertices such
        that for all `u \in V(G)-M, v,w\in M` the relation `u \sim v
        \Leftrightarrow u \sim w` holds, where `\sim` denotes the adjacency
        relation in `G`. Equivalently, `M \subset V(G)` is a module if all its
        vertices have the same adjacency relations with each vertex outside of
        the module (vertex by vertex).

        Hence, for a set like a module, it is very easy to encode the
        information of the adjacencies between the vertices inside and outside
        the module -- we can actually add a new vertex `v_M` to our graph
        representing our module `M`, and let `v_M` be adjacent to `u\in V(G)-M`
        if and only if some `v\in M` (and hence all the vertices contained in
        the module) is adjacent to `u`. We can now independently (and
        recursively) study the structure of our module `M` and the new graph
        `G-M+\{v_M\}`, without any loss of information.

        Here are two very simple modules :

        * A connected component `C` (or the union of some --but not all-- of
          them) of a disconnected graph `G`, for instance, is a module, as no
          vertex of `C` has a neighbor outside of it.

        * An anticomponent `C` (or the union of some --but not all-- of them) of
          an non-anticonnected graph `G`, for the same reason (it is just the
          complement of the previous graph !).

        These modules being of special interest, the disjoint union of graphs is
        called a Parallel composition, and the complement of a disjoint union is
        called a Series composition. A graph whose only modules are singletons
        is called Prime.

        For more information on modular decomposition, in particular for an
        explanation of the terms "Parallel," "Prime" and "Series," see the
        :wikipedia:`Modular_decomposition`.

        You may also be interested in the survey from Michel Habib and
        Christophe Paul entitled "A survey on Algorithmic aspects of modular
        decomposition" [HP2010]_.

        EXAMPLES:

        The Bull Graph is prime::

            sage: graphs.BullGraph().modular_decomposition()
            (PRIME, [1, 2, 0, 3, 4])

        The Petersen Graph too::

            sage: graphs.PetersenGraph().modular_decomposition()
            (PRIME, [1, 4, 5, 0, 2, 6, 3, 7, 8, 9])

        This a clique on 5 vertices with 2 pendant edges, though, has a more
        interesting decomposition::

            sage: g = graphs.CompleteGraph(5)
            sage: g.add_edge(0,5)
            sage: g.add_edge(0,6)
            sage: g.modular_decomposition(algorithm='habib')
            (SERIES, [(PARALLEL, [(SERIES, [1, 2, 3, 4]), 5, 6]), 0])

        We get an equivalent tree when we use the algorithm of [TCHP2008]_::

            sage: g.modular_decomposition(algorithm='tedder')
            (SERIES, [(PARALLEL, [(SERIES, [4, 3, 2, 1]), 5, 6]), 0])

        We can choose output to be a
        :class:`~sage.combinat.rooted_tree.LabelledRootedTree`::

            sage: g.modular_decomposition(style='tree')
            SERIES[0[], PARALLEL[5[], 6[], SERIES[1[], 2[], 3[], 4[]]]]
            sage: ascii_art(g.modular_decomposition(style='tree'))
              __SERIES
             /      /
            0   ___PARALLEL
               / /     /
              5 6   __SERIES
                   / / / /
                  1 2 3 4

        ALGORITHM:

        When ``algorithm='tedder'`` this function uses python implementation of
        algorithm published by Marc Tedder, Derek Corneil, Michel Habib and
        Christophe Paul [TCHP2008]_. When ``algorithm='habib'`` this function
        uses the algorithm of M. Habib and M. Maurer [HM1979]_.

        .. SEEALSO::

            - :meth:`is_prime` -- Tests whether a graph is prime.

            - :class:`~sage.combinat.rooted_tree.LabelledRootedTree`.

        TESTS:

        Empty graph::

            sage: graphs.EmptyGraph().modular_decomposition(algorithm='habib')
            ()
            sage: graphs.EmptyGraph().modular_decomposition(algorithm='tedder')
            ()
            sage: graphs.EmptyGraph().modular_decomposition(algorithm='habib', style='tree')
            None[]
            sage: graphs.EmptyGraph().modular_decomposition(algorithm='tedder', style='tree')
            None[]

        Singleton Vertex::

            sage: Graph(1).modular_decomposition(algorithm='habib')
            (PRIME, [0])
            sage: Graph(1).modular_decomposition(algorithm='tedder')
            (PRIME, [0])
            sage: Graph(1).modular_decomposition(algorithm='habib', style='tree')
            PRIME[0[]]
            sage: Graph(1).modular_decomposition(algorithm='tedder', style='tree')
            PRIME[0[]]

        Vertices may be arbitrary --- check that :trac:`24898` is fixed::

            sage: md = Graph({(1,2):[(2,3)],(2,3):[(1,2)]}).modular_decomposition()
            sage: md[0]
            SERIES
            sage: sorted(md[1])
            [(1, 2), (2, 3)]

        Unknown algorithm::

            sage: graphs.PathGraph(2).modular_decomposition(algorithm='abc')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be 'habib' or 'tedder'

        Unknown style::

            sage: graphs.PathGraph(2).modular_decomposition(style='xyz')
            Traceback (most recent call last):
            ...
            ValueError: style must be 'tuple' or 'tree'
        """
        from sage.graphs.graph_decompositions.modular_decomposition import (modular_decomposition,
                                                                            NodeType,
                                                                            habib_maurer_algorithm,
                                                                            create_prime_node,
                                                                            create_normal_node)

        self._scream_if_not_simple()

        if not self.order():
            D = None
        elif self.order() == 1:
            D = create_prime_node()
            D.children.append(create_normal_node(self.vertices()[0]))
        else:
            if algorithm == 'habib':
                D = habib_maurer_algorithm(self)
            elif algorithm == 'tedder':
                D = modular_decomposition(self)
            else:
                raise ValueError("algorithm must be 'habib' or 'tedder'")

        if style == 'tuple':
            if D is None:
                return tuple()
            def relabel(x):
                if x.node_type == NodeType.NORMAL:
                    return x.children[0]
                else:
                    return x.node_type, [relabel(y) for y in x.children]
            return relabel(D)
        elif style == 'tree':
            from sage.combinat.rooted_tree import LabelledRootedTree
            if D is None:
                return LabelledRootedTree([])
            def to_tree(x):
                if x.node_type == NodeType.NORMAL:
                    return LabelledRootedTree([], label=x.children[0])
                else:
                    return LabelledRootedTree([to_tree(y) for y in x.children], label=x.node_type)
            return to_tree(D)
        else:
            raise ValueError("style must be 'tuple' or 'tree'")

    @doc_index("Graph properties")
    def is_polyhedral(self):
        """
        Check whether the graph is the graph of the polyhedron.

        By a theorem of Steinitz (Satz 43, p. 77 of [St1922]_), graphs of
        three-dimensional polyhedra are exactly the simple 3-vertex-connected
        planar graphs.

        EXAMPLES::

            sage: C = graphs.CubeGraph(3)
            sage: C.is_polyhedral()
            True
            sage: K33=graphs.CompleteBipartiteGraph(3, 3)
            sage: K33.is_polyhedral()
            False
            sage: graphs.CycleGraph(17).is_polyhedral()
            False
            sage: [i for i in range(9) if graphs.CompleteGraph(i).is_polyhedral()]
            [4]

        .. SEEALSO::

            * :meth:`~sage.graphs.generic_graph.GenericGraph.vertex_connectivity`
            * :meth:`~sage.graphs.generic_graph.GenericGraph.is_planar`
            * :meth:`is_circumscribable`
            * :meth:`is_inscribable`
            * :wikipedia:`Polyhedral_graph`

        TESTS::

            sage: G = Graph([[1, 2, 3, 4], [[1, 2], [1,1]]], loops=True)
            sage: G.is_polyhedral()
            False

            sage: G = Graph([[1, 2, 3], [[1, 2], [3, 1], [1, 2], [2, 3]]], multiedges=True)
            sage: G.is_polyhedral()
            False

        """
        return (not self.has_loops()
                and not self.has_multiple_edges()
                and self.vertex_connectivity(k=3)
                and self.is_planar())

    @doc_index("Graph properties")
    def is_circumscribable(self, solver="ppl", verbose=0):
        """
        Test whether the graph is the graph of a circumscribed polyhedron.

        A polyhedron is circumscribed if all of its facets are tangent to a
        sphere. By a theorem of Rivin ([HRS1993]_), this can be checked by
        solving a linear program that assigns weights between 0 and 1/2 on each
        edge of the polyhedron, so that the weights on any face add to exactly
        one and the weights on any non-facial cycle add to more than one.  If
        and only if this can be done, the polyhedron can be circumscribed.

        INPUT:

        - ``solver`` -- (default: ``"ppl"``); specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        EXAMPLES::

            sage: C = graphs.CubeGraph(3)
            sage: C.is_circumscribable()
            True

            sage: O = graphs.OctahedralGraph()
            sage: O.is_circumscribable()
            True

            sage: TT = polytopes.truncated_tetrahedron().graph()
            sage: TT.is_circumscribable()
            False

        Stellating in a face of the octahedral graph is not circumscribable::

            sage: f = set(flatten(choice(O.faces())))
            sage: O.add_edges([[6, i] for i in f])
            sage: O.is_circumscribable()
            False

        .. SEEALSO::

            * :meth:`is_polyhedral`
            * :meth:`is_inscribable`

        TESTS::

            sage: G = graphs.CompleteGraph(5)
            sage: G.is_circumscribable()
            Traceback (most recent call last):
            ...
            NotImplementedError: this method only works for polyhedral graphs

        .. TODO::

            Allow the use of other, inexact but faster solvers.
        """
        if not self.is_polyhedral():
            raise NotImplementedError('this method only works for polyhedral graphs')

        from sage.numerical.mip import MixedIntegerLinearProgram
        from sage.numerical.mip import MIPSolverException
        # For a description of the algorithm see paper by Rivin and:
        # https://www.ics.uci.edu/~eppstein/junkyard/uninscribable/
        # In order to simulate strict inequalities in the following LP, we
        # introduce a variable c[0] and maximize it. If it is positive, then
        # the LP has a solution, such that all inequalities are strict
        # after removing the auxiliary variable c[0].
        M = MixedIntegerLinearProgram(maximization=True, solver=solver)
        e_var = M.new_variable(nonnegative=True)
        c = M.new_variable()
        M.set_min(c[0], -1)
        M.set_max(c[0], 1)
        M.set_objective(c[0])

        for e in self.edge_iterator(labels=0):
            fe = frozenset(e)
            M.set_max(e_var[fe], ZZ(1)/ZZ(2))
            M.add_constraint(e_var[fe] - c[0], min=0)
            M.add_constraint(e_var[fe] + c[0], max=ZZ(1)/ZZ(2))

        # The faces are completely determined by the graph structure:
        # for polyhedral graph, there is only one way to choose the faces.
        # We add an equality constraint for each face.
        efaces = self.faces()
        vfaces = set(frozenset([e[0] for e in face]) for face in efaces)
        for edges in efaces:
            M.add_constraint(M.sum(e_var[frozenset(e)] for e in edges) == 1)

        # In order to generate all simple cycles of G, which are not faces,
        # we use the "all_simple_cycles" method of directed graphs, generating
        # each cycle twice (in both directions). The set below make sure only
        # one direction gives rise to an (in)equality
        D = self.to_directed()
        inequality_constraints = set()
        for cycle in D.all_simple_cycles():
            if len(cycle) > 3:
                scycle = frozenset(cycle)
                if scycle not in vfaces:
                    edges = (frozenset((cycle[i], cycle[i+1])) for i in range(len(cycle)-1))
                    inequality_constraints.add(frozenset(edges))

        for ieq in inequality_constraints:
            M.add_constraint(M.sum(e_var[fe] for fe in ieq) - c[0] >= 1)

        try:
            solution = M.solve(log=verbose)
        except MIPSolverException as msg:
            if str(msg) == "PPL : There is no feasible solution":
                return False
        return solution > 0

    @doc_index("Graph properties")
    def is_inscribable(self, solver="ppl", verbose=0):
        """
        Test whether the graph is the graph of an inscribed polyhedron.

        A polyhedron is inscribed if all of its vertices are on a sphere.
        This is dual to the notion of circumscribed polyhedron: A Polyhedron is
        inscribed if and only if its polar dual is circumscribed and hence a
        graph is inscribable if and only if its planar dual is circumscribable.

        INPUT:

        - ``solver`` -- (default: ``"ppl"``); specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        EXAMPLES::

            sage: H = graphs.HerschelGraph()
            sage: H.is_inscribable()               # long time (> 1 sec)
            False
            sage: H.planar_dual().is_inscribable() # long time (> 1 sec)
            True

            sage: C = graphs.CubeGraph(3)
            sage: C.is_inscribable()
            True

        Cutting off a vertex from the cube yields an uninscribable graph::

            sage: C = graphs.CubeGraph(3)
            sage: v = next(C.vertex_iterator())
            sage: triangle = [_ + v for _ in C.neighbors(v)]
            sage: C.add_edges(Combinations(triangle, 2))
            sage: C.add_edges(zip(triangle, C.neighbors(v)))
            sage: C.delete_vertex(v)
            sage: C.is_inscribable()
            False

        Breaking a face of the cube yields an uninscribable graph::

            sage: C = graphs.CubeGraph(3)
            sage: face = choice(C.faces())
            sage: C.add_edge([face[0][0], face[2][0]])
            sage: C.is_inscribable()
            False


        .. SEEALSO::

            * :meth:`is_polyhedral`
            * :meth:`is_circumscribable`

        TESTS::

            sage: G = graphs.CompleteBipartiteGraph(3,3)
            sage: G.is_inscribable()
            Traceback (most recent call last):
            ...
            NotImplementedError: this method only works for polyhedral graphs
        """
        if not self.is_polyhedral():
            raise NotImplementedError('this method only works for polyhedral graphs')
        return self.planar_dual().is_circumscribable(solver=solver, verbose=verbose)

    @doc_index("Graph properties")
    def is_prime(self, algorithm='habib'):
        r"""
        Test whether the current graph is prime.

        INPUT:

        - ``algorithm`` -- (default: ``'tedder'``) specifies the algorithm to
          use among:

          - ``'tedder'`` -- Use the linear algorithm of [TCHP2008]_.

          - ``'habib'`` -- Use the $O(n^3)$ algorithm of [HM1979]_. This is
            probably slower, but is much simpler and so possibly less error
            prone.

        A graph is prime if all its modules are trivial (i.e. empty, all of the
        graph or singletons) -- see :meth:`modular_decomposition`.

        EXAMPLES:

        The Petersen Graph and the Bull Graph are both prime::

            sage: graphs.PetersenGraph().is_prime()
            True
            sage: graphs.BullGraph().is_prime()
            True

        Though quite obviously, the disjoint union of them is not::

            sage: (graphs.PetersenGraph() + graphs.BullGraph()).is_prime()
            False

        TESTS::

            sage: graphs.EmptyGraph().is_prime()
            True
        """
        from sage.graphs.graph_decompositions.modular_decomposition import NodeType

        if self.order() <= 1:
            return True

        D = self.modular_decomposition(algorithm=algorithm)

        return D[0] == NodeType.PRIME and len(D[1]) == self.order()

    def _gomory_hu_tree(self, vertices, algorithm=None):
        r"""
        Return a Gomory-Hu tree associated to self.

        This function is the private counterpart of ``gomory_hu_tree()``, with
        the difference that it has an optional argument needed for recursive
        computations, which the user is not interested in defining himself.

        See the documentation of ``gomory_hu_tree()`` for more information.

        INPUT:

        - ``vertices`` -- a set of "real" vertices, as opposed to the fakes one
          introduced during the computations. This variable is useful for the
          algorithm and for recursion purposes.

        - ``algorithm`` -- select the algorithm used by the :meth:`edge_cut`
          method. Refer to its documentation for allowed values and default
          behaviour.

        EXAMPLES:

        This function is actually tested in ``gomory_hu_tree()``, this example
        is only present to have a doctest coverage of 100%::

            sage: g = graphs.PetersenGraph()
            sage: t = g._gomory_hu_tree(frozenset(g.vertices()))
        """
        self._scream_if_not_simple()

        # Small case, not really a problem ;-)
        if len(vertices) == 1:
            g = Graph()
            g.add_vertices(vertices)
            return g

        # Take any two vertices (u,v)
        it = iter(vertices)
        u,v = next(it),next(it)

        # Compute a uv min-edge-cut.
        #
        # The graph is split into U,V with u \in U and v\in V.
        flow,edges,[U,V] = self.edge_cut(u, v, use_edge_labels=True, vertices=True, algorithm=algorithm)

        # One graph for each part of the previous one
        gU,gV = self.subgraph(U, immutable=False), self.subgraph(V, immutable=False)

        # A fake vertex fU (resp. fV) to represent U (resp. V)
        fU = frozenset(U)
        fV = frozenset(V)

        # Each edge (uu,vv) with uu \in U and vv\in V yields:
        # - an edge (uu,fV) in gU
        # - an edge (vv,fU) in gV
        #
        # If the same edge is added several times their capacities add up.

        from sage.rings.real_mpfr import RR
        for uu,vv,capacity in edges:
            capacity = capacity if capacity in RR else 1

            # Assume uu is in gU
            if uu in V:
                uu,vv = vv,uu

            # Create the new edges if necessary
            if not gU.has_edge(uu, fV):
                gU.add_edge(uu, fV, 0)
            if not gV.has_edge(vv, fU):
                gV.add_edge(vv, fU, 0)

            # update the capacities
            gU.set_edge_label(uu, fV, gU.edge_label(uu, fV) + capacity)
            gV.set_edge_label(vv, fU, gV.edge_label(vv, fU) + capacity)

        # Recursion on each side
        gU_tree = gU._gomory_hu_tree(vertices & frozenset(gU), algorithm=algorithm)
        gV_tree = gV._gomory_hu_tree(vertices & frozenset(gV), algorithm=algorithm)

        # Union of the two partial trees
        g = gU_tree.union(gV_tree)

        # An edge to connect them, with the appropriate label
        g.add_edge(u, v, flow)

        return g

    @doc_index("Connectivity, orientations, trees")
    def gomory_hu_tree(self, algorithm=None):
        r"""
        Return a Gomory-Hu tree of self.

        Given a tree `T` with labeled edges representing capacities, it is very
        easy to determine the maximum flow between any pair of vertices :
        it is the minimal label on the edges of the unique path between them.

        Given a graph `G`, a Gomory-Hu tree `T` of `G` is a tree with the same
        set of vertices, and such that the maximum flow between any two vertices
        is the same in `G` as in `T`. See the :wikipedia:`Gomory–Hu_tree`. Note
        that, in general, a graph admits more than one Gomory-Hu tree.

        See also 15.4 (Gomory-Hu trees) from [Sch2003]_.

        INPUT:

        - ``algorithm`` -- select the algorithm used by the :meth:`edge_cut`
          method. Refer to its documentation for allowed values and default
          behaviour.

        OUTPUT:

        A graph with labeled edges

        EXAMPLES:

        Taking the Petersen graph::

            sage: g = graphs.PetersenGraph()
            sage: t = g.gomory_hu_tree()

        Obviously, this graph is a tree::

            sage: t.is_tree()
            True

        Note that if the original graph is not connected, then the Gomory-Hu
        tree is in fact a forest::

            sage: (2*g).gomory_hu_tree().is_forest()
            True
            sage: (2*g).gomory_hu_tree().is_connected()
            False

        On the other hand, such a tree has lost nothing of the initial graph
        connectedness::

            sage: all(t.flow(u,v) == g.flow(u,v) for u,v in Subsets(g.vertices(), 2))
            True

        Just to make sure, we can check that the same is true for two vertices
        in a random graph::

            sage: g = graphs.RandomGNP(20,.3)
            sage: t = g.gomory_hu_tree()
            sage: g.flow(0,1) == t.flow(0,1)
            True

        And also the min cut::

            sage: g.edge_connectivity() == min(t.edge_labels()) or not g.is_connected()
            True

        TESTS:

        :trac:`16475`::

            sage: G = graphs.PetersenGraph()
            sage: for u,v in G.edge_iterator(labels=False):
            ....:     G.set_edge_label(u, v, 1)
            sage: for u, v in [(0, 1), (0, 4), (0, 5), (1, 2), (1, 6), (3, 4), (5, 7), (5, 8)]:
            ....:     G.set_edge_label(u, v, 2)
            sage: T = G.gomory_hu_tree()
            sage: from itertools import combinations
            sage: for u,v in combinations(G,2):
            ....:     assert T.flow(u,v,use_edge_labels=True) == G.flow(u,v,use_edge_labels=True)

            sage: graphs.EmptyGraph().gomory_hu_tree()
            Graph on 0 vertices
        """
        if not self.order():
            return Graph()
        if not self.is_connected():
            g = Graph()
            for cc in self.connected_components_subgraphs():
                g = g.union(cc._gomory_hu_tree(frozenset(cc.vertex_iterator()), algorithm=algorithm))
        else:
            g = self._gomory_hu_tree(frozenset(self.vertex_iterator()), algorithm=algorithm)

        if self.get_pos() is not None:
            g.set_pos(dict(self.get_pos()))
        return g

    @doc_index("Leftovers")
    def two_factor_petersen(self, solver=None, verbose=0, *, integrality_tolerance=1e-3):
        r"""
        Return a decomposition of the graph into 2-factors.

        Petersen's 2-factor decomposition theorem asserts that any `2r`-regular
        graph `G` can be decomposed into 2-factors.  Equivalently, it means that
        the edges of any `2r`-regular graphs can be partitionned in `r` sets
        `C_1,\dots,C_r` such that for all `i`, the set `C_i` is a disjoint union
        of cycles (a 2-regular graph).

        As any graph of maximal degree `\Delta` can be completed into a regular
        graph of degree `2\lceil\frac\Delta 2\rceil`, this result also means
        that the edges of any graph of degree `\Delta` can be partitionned in
        `r=2\lceil\frac\Delta 2\rceil` sets `C_1,\dots,C_r` such that for all
        `i`, the set `C_i` is a graph of maximal degree `2` (a disjoint union of
        paths and cycles).

        INPUT:

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        EXAMPLES:

        The Complete Graph on `7` vertices is a `6`-regular graph, so it can be
        edge-partitionned into `2`-regular graphs::

            sage: g = graphs.CompleteGraph(7)
            sage: classes = g.two_factor_petersen()
            sage: for c in classes:
            ....:     gg = Graph()
            ....:     gg.add_edges(c)
            ....:     print(max(gg.degree())<=2)
            True
            True
            True
            sage: Set(set(classes[0]) | set(classes[1]) | set(classes[2])).cardinality() == g.size()
            True

        ::

            sage: g = graphs.CirculantGraph(24, [7, 11])
            sage: cl = g.two_factor_petersen()
            sage: g.plot(edge_colors={'black':cl[0], 'red':cl[1]})
            Graphics object consisting of 73 graphics primitives

        """
        self._scream_if_not_simple()
        d = self.eulerian_orientation()

        # This new graph is bipartite, and built the following way :
        #
        # To each vertex v of the digraph are associated two vertices,
        # a sink (-1,v) and a source (1,v)
        # Any edge (u,v) in the digraph is then added as ((-1,u),(1,v))

        g = Graph()
        g.add_edges(((-1, u), (1, v)) for u, v in d.edge_iterator(labels=None))

        # This new bipartite graph is now edge_colored
        from sage.graphs.graph_coloring import edge_coloring
        classes = edge_coloring(g, solver=solver, verbose=verbose,
                                integrality_tolerance=integrality_tolerance)

        # The edges in the classes are of the form ((-1,u),(1,v))
        # and have to be translated back to (u,v)
        classes_b = []
        for c in classes:
            classes_b.append([(u,v) for ((uu,u),(vv,v)) in c])

        return classes_b

    @doc_index("Leftovers")
    def kirchhoff_symanzik_polynomial(self, name='t'):
        r"""
        Return the Kirchhoff-Symanzik polynomial of a graph.

        This is a polynomial in variables `t_e` (each of them representing an
        edge of the graph `G`) defined as a sum over all spanning trees:

        .. MATH::

            \Psi_G(t) = \sum_{\substack{T\subseteq V \\ \text{a spanning tree}}} \prod_{e \not\in E(T)} t_e

        This is also called the first Symanzik polynomial or the Kirchhoff
        polynomial.

        INPUT:

        - ``name`` -- name of the variables (default: ``'t'``)

        OUTPUT:

        - a polynomial with integer coefficients

        ALGORITHM:

            This is computed here using a determinant, as explained in Section
            3.1 of [Mar2009a]_.

            As an intermediate step, one computes a cycle basis `\mathcal C` of
            `G` and a rectangular `|\mathcal C| \times |E(G)|` matrix with
            entries in `\{-1,0,1\}`, which describes which edge belong to which
            cycle of `\mathcal C` and their respective orientations.

            More precisely, after fixing an arbitrary orientation for each edge
            `e\in E(G)` and each cycle `C\in\mathcal C`, one gets a sign for
            every incident pair (edge, cycle) which is `1` if the orientation
            coincide and `-1` otherwise.

        EXAMPLES:

        For the cycle of length 5::

            sage: G = graphs.CycleGraph(5)
            sage: G.kirchhoff_symanzik_polynomial()
            t0 + t1 + t2 + t3 + t4

        One can use another letter for variables::

            sage: G.kirchhoff_symanzik_polynomial(name='u')
            u0 + u1 + u2 + u3 + u4

        For the 'coffee bean' graph::

            sage: G = Graph([(0,1,'a'),(0,1,'b'),(0,1,'c')], multiedges=True)
            sage: G.kirchhoff_symanzik_polynomial()
            t0*t1 + t0*t2 + t1*t2

        For the 'parachute' graph::

            sage: G = Graph([(0,2,'a'),(0,2,'b'),(0,1,'c'),(1,2,'d')], multiedges=True)
            sage: G.kirchhoff_symanzik_polynomial()
            t0*t1 + t0*t2 + t1*t2 + t1*t3 + t2*t3

        For the complete graph with 4 vertices::

            sage: G = graphs.CompleteGraph(4)
            sage: G.kirchhoff_symanzik_polynomial()
            t0*t1*t3 + t0*t2*t3 + t1*t2*t3 + t0*t1*t4 + t0*t2*t4 + t1*t2*t4
            + t1*t3*t4 + t2*t3*t4 + t0*t1*t5 + t0*t2*t5 + t1*t2*t5 + t0*t3*t5
            + t2*t3*t5 + t0*t4*t5 + t1*t4*t5 + t3*t4*t5

        REFERENCES:

        [Bro2011]_
        """
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        # The order of the vertices in each tuple matters, so use a list
        edges = list(self.edges(sort=False))
        cycles = self.cycle_basis(output='edge')

        edge2int = {e: j for j, e in enumerate(edges)}
        circuit_mtrx = matrix(ZZ, self.size(), len(cycles))
        for i, cycle in enumerate(cycles):
            for edge in cycle:
                if edge in edges:
                    circuit_mtrx[edge2int[edge], i] = +1
                else:
                    circuit_mtrx[edge2int[(edge[1], edge[0], edge[2])], i] = -1

        D = matrix.diagonal(PolynomialRing(ZZ, name, self.size()).gens())
        return (circuit_mtrx.transpose() * D * circuit_mtrx).determinant()

    @doc_index("Leftovers")
    def magnitude_function(self):
        r"""
        Return the magnitude function of the graph as a rational function.

        This is defined as the sum of all coefficients in the inverse of the
        matrix `Z` whose coefficient `Z_{i,j}` indexed by a pair of vertices
        `(i,j)` is `q^d(i,j)` where `d` is the distance function in the graph.

        By convention, if the distance from `i` to `j` is infinite (for two
        vertices not path connected) then `Z_{i,j}=0`.

        The value of the magnitude function at `q=0` is the cardinality of the
        graph. The magnitude function of a disjoint union is the sum of the
        magnitudes functions of the connected components. The magnitude function
        of a Cartesian product is the product of the magnitudes functions of the
        factors.

        EXAMPLES::

            sage: g = Graph({1:[], 2:[]})
            sage: g.magnitude_function()
            2

            sage: g = graphs.CycleGraph(4)
            sage: g.magnitude_function()
            4/(q^2 + 2*q + 1)

            sage: g = graphs.CycleGraph(5)
            sage: m = g.magnitude_function(); m
            5/(2*q^2 + 2*q + 1)

        One can expand the magnitude as a power series in `q` as follows::

            sage: q = QQ[['q']].gen()
            sage: m(q)
            5 - 10*q + 10*q^2 - 20*q^4 + 40*q^5 - 40*q^6 + ...

        One can also use the substitution `q = exp(-t)` to obtain the magnitude
        function as a function of `t`::

            sage: g = graphs.CycleGraph(6)
            sage: m = g.magnitude_function()
            sage: t = var('t')                                                  # optional - sage.symbolic
            sage: m(exp(-t))                                                    # optional - sage.symbolic
            6/(2*e^(-t) + 2*e^(-2*t) + e^(-3*t) + 1)

        TESTS::

            sage: g = Graph()
            sage: g.magnitude_function()
            0

            sage: g = Graph({1:[]})
            sage: g.magnitude_function()
            1

            sage: g = graphs.PathGraph(4)
            sage: g.magnitude_function()
            (-2*q + 4)/(q + 1)

        REFERENCES:

        .. [Lein] Tom Leinster, *The magnitude of metric spaces*.
           Doc. Math. 18 (2013), 857-905.
        """
        from sage.matrix.constructor import matrix
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.graphs.distances_all_pairs import distances_all_pairs

        ring = PolynomialRing(ZZ, 'q')
        q = ring.gen()
        N = self.order()
        if not N:
            return ring.zero()
        dist = distances_all_pairs(self)
        vertices = list(self)
        Z = matrix(ring, N, N, ring.zero())
        for i in range(N):
            Z[i, i] = ring.one()
        for i in range(N):
            for j in range(i):
                dij = dist[vertices[i]][vertices[j]]
                if dij in ZZ:
                    Z[i, j] = Z[j, i] = q ** dij
                else:
                    Z[i, j] = Z[j, i] = ring.zero()
        return sum(sum(u) for u in ~Z)

    @doc_index("Leftovers")
    def ihara_zeta_function_inverse(self):
        """
        Compute the inverse of the Ihara zeta function of the graph.

        This is a polynomial in one variable with integer coefficients. The
        Ihara zeta function itself is the inverse of this polynomial.

        See the :wikipedia:`Ihara zeta function` for more information.

        ALGORITHM:

        This is computed here as the (reversed) characteristic polynomial of a
        square matrix of size twice the number of edges, related to the
        adjacency matrix of the line graph, see for example Proposition 9 in
        [SS2008]_ and Def. 4.1 in [Ter2011]_.

        The graph is first replaced by its 2-core, as this does not change the
        Ihara zeta function.

        EXAMPLES::

            sage: G = graphs.CompleteGraph(4)
            sage: factor(G.ihara_zeta_function_inverse())
            (2*t - 1) * (t + 1)^2 * (t - 1)^3 * (2*t^2 + t + 1)^3

            sage: G = graphs.CompleteGraph(5)
            sage: factor(G.ihara_zeta_function_inverse())
            (-1) * (3*t - 1) * (t + 1)^5 * (t - 1)^6 * (3*t^2 + t + 1)^4

            sage: G = graphs.PetersenGraph()
            sage: factor(G.ihara_zeta_function_inverse())
            (-1) * (2*t - 1) * (t + 1)^5 * (t - 1)^6 * (2*t^2 + 2*t + 1)^4
            * (2*t^2 - t + 1)^5

            sage: G = graphs.RandomTree(10)
            sage: G.ihara_zeta_function_inverse()
            1

        REFERENCES:

        [HST2001]_
        """
        from sage.matrix.constructor import matrix

        H = self.subgraph(vertices=self.cores(k=2)[1])
        E = list(H.edges(sort=False))
        m = len(E)
        # compute (Hashimoto) edge matrix T
        T = matrix(ZZ, 2 * m, 2 * m, 0)
        for i in range(m):
            for j in range(m):
                if i != j:
                    if E[i][1] == E[j][0]:  # same orientation
                        T[2 * i, 2 * j] = 1
                        T[2 * j + 1, 2 * i + 1] = 1
                    elif E[i][1] == E[j][1]:  # opposite orientation (towards)
                        T[2 * i, 2 * j + 1] = 1
                        T[2 * j, 2 * i + 1] = 1
                    elif E[i][0] == E[j][0]:  # opposite orientation (away)
                        T[2 * i + 1, 2 * j] = 1
                        T[2 * j + 1, 2 * i] = 1
        return T.charpoly('t').reverse()

    @doc_index("Leftovers")
    def perfect_matchings(self, labels=False):
        r"""
        Return an iterator over all perfect matchings of the graph.

        ALGORITHM:

        Choose a vertex `v`, then recurse through all edges incident to `v`,
        removing one edge at a time whenever an edge is added to a matching.

        INPUT:

        - ``labels`` -- boolean (default: ``False``); when ``True``, the edges
          in each perfect matching are triples (containing the label as the
          third element), otherwise the edges are pairs.

        .. SEEALSO::

            :meth:`matching`

        EXAMPLES::

            sage: G=graphs.GridGraph([2,3])
            sage: for m in G.perfect_matchings():
            ....:     print(sorted(m))
            [((0, 0), (0, 1)), ((0, 2), (1, 2)), ((1, 0), (1, 1))]
            [((0, 0), (1, 0)), ((0, 1), (0, 2)), ((1, 1), (1, 2))]
            [((0, 0), (1, 0)), ((0, 1), (1, 1)), ((0, 2), (1, 2))]

            sage: G = graphs.CompleteGraph(4)
            sage: for m in G.perfect_matchings(labels=True):
            ....:     print(sorted(m))
            [(0, 1, None), (2, 3, None)]
            [(0, 2, None), (1, 3, None)]
            [(0, 3, None), (1, 2, None)]

            sage: G = Graph([[1,-1,'a'], [2,-2, 'b'], [1,-2,'x'], [2,-1,'y']])
            sage: sorted(sorted(m) for m in G.perfect_matchings(labels=True))
            [[(-2, 1, 'x'), (-1, 2, 'y')], [(-2, 2, 'b'), (-1, 1, 'a')]]

            sage: G = graphs.CompleteGraph(8)
            sage: mpc = G.matching_polynomial().coefficients(sparse=False)[0]
            sage: len(list(G.perfect_matchings())) == mpc
            True

            sage: G = graphs.PetersenGraph().copy(immutable=True)
            sage: [sorted(m) for m in G.perfect_matchings()]
            [[(0, 1), (2, 3), (4, 9), (5, 7), (6, 8)],
             [(0, 1), (2, 7), (3, 4), (5, 8), (6, 9)],
             [(0, 4), (1, 2), (3, 8), (5, 7), (6, 9)],
             [(0, 4), (1, 6), (2, 3), (5, 8), (7, 9)],
             [(0, 5), (1, 2), (3, 4), (6, 8), (7, 9)],
             [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]]

            sage: list(Graph().perfect_matchings())
            [[]]

            sage: G = graphs.CompleteGraph(5)
            sage: list(G.perfect_matchings())
            []
        """
        if not self:
            yield []
            return
        if self.order() % 2 or any(len(cc) % 2 for cc in self.connected_components()):
            return

        def rec(G):
            """
            Iterator over all perfect matchings of a simple graph `G`.
            """
            if not G:
                yield []
                return
            if G.order() % 2 == 0:
                v = next(G.vertex_iterator())
                Nv = list(G.neighbor_iterator(v))
                G.delete_vertex(v)
                for u in Nv:
                    Nu = list(G.neighbor_iterator(u))
                    G.delete_vertex(u)
                    for partial_matching in rec(G):
                        partial_matching.append((u, v))
                        yield partial_matching
                    G.add_vertex(u)
                    G.add_edges((u, nu) for nu in Nu)
                G.add_vertex(v)
                G.add_edges((v, nv) for nv in Nv)

        # We create a mutable copy of the graph and remove its loops, if any
        G = self.copy(immutable=False)
        G.allow_loops(False)

        # We create a mapping from frozen unlabeled edges to (labeled) edges.
        # This ease for instance the manipulation of multiedges (if any)
        edges = {}
        for e in G.edges(labels=labels):
            f = frozenset(e[:2])
            if f in edges:
                edges[f].append(e)
            else:
                edges[f] = [e]

        # We now get rid of multiple edges, if any
        G.allow_multiple_edges(False)

        # For each unlabeled matching, we yield all its possible labelings
        for m in rec(G):
            for pm in itertools.product(*[edges[frozenset(e)] for e in m]):
                yield pm

    @doc_index("Leftovers")
    def has_perfect_matching(self, algorithm="Edmonds", solver=None, verbose=0,
                             *, integrality_tolerance=1e-3):
        r"""
        Return whether this graph has a perfect matching.
        INPUT:

        - ``algorithm`` -- string (default: ``"Edmonds"``)

          - ``"Edmonds"`` uses Edmonds' algorithm as implemented in NetworkX to
            find a matching of maximal cardinality, then check whether this
            cardinality is half the number of vertices of the graph.

          - ``"LP_matching"`` uses a Linear Program to find a matching of
            maximal cardinality, then check whether this cardinality is half the
            number of vertices of the graph.

          - ``"LP"`` uses a Linear Program formulation of the perfect matching
            problem: put a binary variable ``b[e]`` on each edge `e`, and for
            each vertex `v`, require that the sum of the values of the edges
            incident to `v` is 1.

        - ``solver`` -- string (default: ``None``); specify a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``); sets the level of verbosity:
          set to 0 by default, which means quiet (only useful when
          ``algorithm == "LP_matching"`` or ``algorithm == "LP"``)

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: graphs.PetersenGraph().has_perfect_matching()
            True
            sage: graphs.WheelGraph(6).has_perfect_matching()
            True
            sage: graphs.WheelGraph(5).has_perfect_matching()
            False
            sage: graphs.PetersenGraph().has_perfect_matching(algorithm="LP_matching")
            True
            sage: graphs.WheelGraph(6).has_perfect_matching(algorithm="LP_matching")
            True
            sage: graphs.WheelGraph(5).has_perfect_matching(algorithm="LP_matching")
            False
            sage: graphs.PetersenGraph().has_perfect_matching(algorithm="LP_matching")
            True
            sage: graphs.WheelGraph(6).has_perfect_matching(algorithm="LP_matching")
            True
            sage: graphs.WheelGraph(5).has_perfect_matching(algorithm="LP_matching")
            False

        TESTS::

            sage: G = graphs.EmptyGraph()
            sage: all(G.has_perfect_matching(algorithm=algo) for algo in ['Edmonds', 'LP_matching', 'LP'])
            True

        Be careful with isolated vertices::

            sage: G = graphs.PetersenGraph()
            sage: G.add_vertex(11)
            sage: any(G.has_perfect_matching(algorithm=algo) for algo in ['Edmonds', 'LP_matching', 'LP'])
            False
        """
        if self.order() % 2:
            return False
        if algorithm == "Edmonds":
            return len(self) == 2*self.matching(value_only=True,
                                                use_edge_labels=False,
                                                algorithm="Edmonds")
        elif algorithm == "LP_matching":
            return len(self) == 2*self.matching(value_only=True,
                                                use_edge_labels=False,
                                                algorithm="LP",
                                                solver=solver,
                                                verbose=verbose,
                                                integrality_tolerance=integrality_tolerance)
        elif algorithm == "LP":
            from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
            p = MixedIntegerLinearProgram(solver=solver)
            b = p.new_variable(binary=True)
            for v in self:
                edges = self.edges_incident(v, labels=False)
                if not edges:
                    return False
                p.add_constraint(p.sum(b[frozenset(e)] for e in edges) == 1)
            try:
                p.solve(log=verbose)
                return True
            except MIPSolverException:
                return False
        else:
            raise ValueError('algorithm must be set to "Edmonds", "LP_matching" or "LP"')

    @doc_index("Leftovers")
    def effective_resistance(self, i, j):
        r"""
        Return the effective resistance between nodes `i` and `j`.

        The resistance distance between vertices `i` and `j` of a simple
        connected graph `G` is defined as the effective resistance between the
        two vertices on an electrical network constructed from `G` replacing
        each edge of the graph by a unit (1 ohm) resistor.

        See the :wikipedia:`Resistance_distance` for more information.

        INPUT:

        - ``i``, ``j`` -- vertices of the graph

        OUTPUT: rational number denoting resistance between nodes `i` and `j`

        EXAMPLES:

        Effective resistances in a straight linear 2-tree on 6 vertices ::

            sage: G = Graph([(0,1),(0,2),(1,2),(1,3),(3,5),(2,4),(2,3),(3,4),(4,5)])
            sage: G.effective_resistance(0,1)
            34/55
            sage: G.effective_resistance(0,3)
            49/55
            sage: G.effective_resistance(1,4)
            9/11
            sage: G.effective_resistance(0,5)
            15/11

        Effective resistances in a fan on 6 vertices ::

            sage: H = Graph([(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,2),(2,3),(3,4),(4,5)])
            sage: H.effective_resistance(1,5)
            6/5
            sage: H.effective_resistance(1,3)
            49/55

        .. SEEALSO::

            * :meth:`effective_resistance_matrix` --
              a similar method giving a matrix full of all effective
              resistances between all nodes

            * :meth:`least_effective_resistance` --
              gives node pairs with least effective resistances

            * See :wikipedia:`Resistance_distance` for more details.

        TESTS::

            sage: G = graphs.CompleteGraph(4)
            sage: all(G.effective_resistance(u, v) == 1/2 for u,v in G.edge_iterator(labels=False))
            True
            sage: Graph(1).effective_resistance(0,0)
            0
            sage: G = Graph([(0,1),(1,2)])
            sage: G.effective_resistance(0,2)
            2
            sage: G = Graph([(0,1),(1,2),(2,0)])
            sage: G.effective_resistance(0,2)
            2/3
            sage: G = Graph([(0,1),(0,2),(0,3),(0,4),(0,5),(1,2),(2,3),(3,4),(4,5),(5,1)])
            sage: r = G.effective_resistance(0,3)
            sage: r == fibonacci(2*(5-3)+1)*fibonacci(2*3-1)/fibonacci(2*5)
            True
        """
        from sage.matrix.constructor import matrix
        if i not in self:
            raise ValueError("vertex ({0}) is not a vertex of the graph".format(repr(i)))
        elif j not in self:
            raise ValueError("vertex ({0}) is not a vertex of the graph".format(repr(j)))

        if i == j :
            return 0

        self._scream_if_not_simple()
        if not self.is_connected():
            raise ValueError('the Graph is not a connected graph')

        vert = list(self)
        i1 = vert.index(i)
        i2 = vert.index(j)
        n = self.order()
        L = self.laplacian_matrix(vertices=vert)
        M = L.pseudoinverse()
        Id = matrix.identity(n)
        sigma = matrix(Id[i1] - Id[i2])
        diff = sigma * M * sigma.transpose()

        return diff[0, 0]

    @doc_index("Leftovers")
    def effective_resistance_matrix(self, vertices=None, nonedgesonly=True):
        r"""
        Return a matrix whose (`i` , `j`) entry gives the effective resistance
        between vertices `i` and `j`.

        The resistance distance between vertices `i` and `j` of a simple
        connected graph `G` is defined as the effective resistance between the
        two vertices on an electrical network constructed from `G` replacing
        each edge of the graph by a unit (1 ohm) resistor.

        INPUT:

        - ``nonedgesonly`` -- boolean (default: ``True``); if ``True`` assign
          zero resistance to pairs of adjacent vertices.

        - ``vertices`` -- list (default: ``None``); the ordering of the
          vertices defining how they should appear in the matrix. By default,
          the ordering given by :meth:`GenericGraph.vertices` is used.

        OUTPUT: matrix

        EXAMPLES:

        The effective resistance matrix  for a straight linear 2-tree counting
        only non-adjacent vertex pairs ::

            sage: G = Graph([(0,1),(0,2),(1,2),(1,3),(3,5),(2,4),(2,3),(3,4),(4,5)])
            sage: G.effective_resistance_matrix()
            [    0     0     0 49/55 59/55 15/11]
            [    0     0     0     0  9/11 59/55]
            [    0     0     0     0     0 49/55]
            [49/55     0     0     0     0     0]
            [59/55  9/11     0     0     0     0]
            [15/11 59/55 49/55     0     0     0]

        The same effective resistance matrix, this time including adjacent
        vertices ::

            sage: G.effective_resistance_matrix(nonedgesonly=False)
            [    0 34/55 34/55 49/55 59/55 15/11]
            [34/55     0 26/55 31/55  9/11 59/55]
            [34/55 26/55     0  5/11 31/55 49/55]
            [49/55 31/55  5/11     0 26/55 34/55]
            [59/55  9/11 31/55 26/55     0 34/55]
            [15/11 59/55 49/55 34/55 34/55     0]

        This example illustrates the common neighbors matrix  for a fan on 6
        vertices counting only non-adjacent vertex pairs ::

            sage: H = Graph([(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,2),(2,3),(3,4),(4,5)])
            sage: H.effective_resistance_matrix()
            [    0     0     0     0     0     0     0]
            [    0     0     0 49/55 56/55   6/5 89/55]
            [    0     0     0     0   4/5 56/55 81/55]
            [    0 49/55     0     0     0 49/55 16/11]
            [    0 56/55   4/5     0     0     0 81/55]
            [    0   6/5 56/55 49/55     0     0 89/55]
            [    0 89/55 81/55 16/11 81/55 89/55     0]

        .. SEEALSO::

            * :meth:`least_effective_resistance` --
              gives node pairs with least effective resistances

            * :meth:`effective_resistance` --
              computes effective resistance for a single node pair

            * See :wikipedia:`Resistance_Distance` for more details.

        TESTS::

            sage: graphs.CompleteGraph(4).effective_resistance_matrix()
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]

            sage: G = Graph(multiedges=True, sparse=True)
            sage: G.add_edges([(0, 1)] * 3)
            sage: G.effective_resistance_matrix()
            Traceback (most recent call last):
            ...
            ValueError: This method is not known to work on graphs with
            multiedges. Perhaps this method can be updated to handle them, but
            in the meantime if you want to use it please disallow multiedges
            using allow_multiple_edges().

            sage: graphs.CompleteGraph(4).effective_resistance_matrix(nonedgesonly=False)
            [  0 1/2 1/2 1/2]
            [1/2   0 1/2 1/2]
            [1/2 1/2   0 1/2]
            [1/2 1/2 1/2   0]
            sage: Graph(1).effective_resistance_matrix()
            [0]
            sage: Graph().effective_resistance_matrix()
            Traceback (most recent call last):
            ...
            ValueError: unable to compute effective resistance for an empty Graph object
            sage: G = Graph([(0,1),(1,2),(2,3),(3,0),(0,2)])
            sage: G.effective_resistance_matrix()
            [0 0 0 0]
            [0 0 0 1]
            [0 0 0 0]
            [0 1 0 0]
            sage: G = Graph([(0,1),(0,2),(0,3),(0,4),(0,5),(1,2),(2,3),(3,4),(4,5),(5,1)])
            sage: r = G.effective_resistance_matrix(nonedgesonly=False)[0,3]
            sage: r == fibonacci(2*(5-3)+1)*fibonacci(2*3-1)/fibonacci(2*5)
            True
        """
        from sage.matrix.constructor import matrix
        from sage.rings.rational_field import QQ

        n = self.order()
        if not n:
            raise ValueError('unable to compute effective resistance for an empty Graph object')
        if vertices is None:
            vertices = self.vertices()
        self._scream_if_not_simple()
        if not self.is_connected():
            raise ValueError('the Graph is not a connected graph')

        L = self.laplacian_matrix(vertices=vertices)
        M = L.pseudoinverse()
        d = matrix(M.diagonal()).transpose()
        onesvec = matrix(QQ, n, 1, lambda i, j: 1)
        S = d * onesvec.transpose() + onesvec * d.transpose() - 2 * M
        onesmat = matrix(QQ, n, n, lambda i, j: 1)
        if nonedgesonly:
            B = onesmat - self.adjacency_matrix(vertices=vertices) - matrix.identity(n)
            S = S.elementwise_product(B)

        return S

    @doc_index("Leftovers")
    def least_effective_resistance(self, nonedgesonly=True):
        r"""
        Return a list of pairs of nodes with the least effective resistance.

        The resistance distance between vertices `i` and `j` of a simple
        connected graph `G` is defined as the effective resistance between the
        two vertices on an electrical network constructed from `G` replacing
        each edge of the graph by a unit (1 ohm) resistor.

        INPUT:

        - ``nonedgesonly`` -- Boolean (default: `True`); if true, assign zero
          resistance to pairs of adjacent vertices

        OUTPUT: list

        EXAMPLES:

        Pairs of non-adjacent nodes with least effective resistance in a
        straight linear 2-tree on 6 vertices::

            sage: G = Graph([(0,1),(0,2),(1,2),(1,3),(3,5),(2,4),(2,3),(3,4),(4,5)])
            sage: G.least_effective_resistance()
            [(1, 4)]

        Pairs of (adjacent or non-adjacent) nodes with least effective
        resistance in a straight linear 2-tree on 6 vertices ::

            sage: G.least_effective_resistance(nonedgesonly = False)
            [(2, 3)]

        Pairs of non-adjacent nodes with least effective resistance in a fan on
        6 vertices counting only non-adjacent vertex pairs ::

            sage: H = Graph([(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,2),(2,3),(3,4),(4,5)])
            sage: H.least_effective_resistance()
            [(2, 4)]

        .. SEEALSO::

            * :meth:`effective_resistance_matrix` --
              a similar method giving a matrix full of all effective
              resistances

            * :meth:`effective_resistance` --
              compuetes effective resistance for a single node pair

            * See :wikipedia:`Resistance_distance` for more details.


        TESTS::

            sage: graphs.CompleteGraph(4).least_effective_resistance()
            []
            sage: graphs.CompleteGraph(4).least_effective_resistance(nonedgesonly=False)
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: Graph(1).least_effective_resistance()
            []
            sage: G = Graph([(0,1),(1,2),(2,3),(3,0),(0,2)])
            sage: G.least_effective_resistance()
            [(1, 3)]
        """
        n = self.order()
        if not n:
            raise ValueError('unable to compute least resistance on empty Graph')
        self._scream_if_not_simple()
        if not self.is_connected():
            raise ValueError('the Graph is not a connected graph')
        if nonedgesonly and self.is_clique():
            return []
        verts = list(self)
        verttoidx = {u: i for i, u in enumerate(verts)}
        S = self.effective_resistance_matrix(vertices=verts, nonedgesonly=nonedgesonly)
        if nonedgesonly:
            edges = self.complement().edges(labels=False)
        else:
            edges = [(verts[i], verts[j]) for i in range(n) for j in range(i + 1, n)]

        rmin = min(S[(verttoidx[e[0]], verttoidx[e[1]])] for e in edges)
        return [e for e in edges if S[(verttoidx[e[0]], verttoidx[e[1]])] == rmin]

    @doc_index("Leftovers")
    def common_neighbors_matrix(self, vertices=None, nonedgesonly=True):
        r"""
        Return a matrix of numbers of common neighbors between each pairs.

        The `(i , j)` entry of the matrix gives the number of common
        neighbors between vertices `i` and `j`.

        This method is only valid for simple (no loops, no multiple edges)
        graphs.

        INPUT:

        - ``nonedgesonly``-- boolean (default: ``True``); if ``True``, assigns
          `0` value to adjacent vertices.

        - ``vertices`` -- list (default: ``None``); the ordering of the
          vertices defining how they should appear in the matrix. By default,
          the ordering given by :meth:`GenericGraph.vertices` is used.

        OUTPUT: matrix

        EXAMPLES:

        The common neighbors matrix  for a straight linear 2-tree counting
        only non-adjacent vertex pairs ::

            sage: G1 = Graph()
            sage: G1.add_edges([(0,1),(0,2),(1,2),(1,3),(3,5),(2,4),(2,3),(3,4),(4,5)])
            sage: G1.common_neighbors_matrix(nonedgesonly = True)
            [0 0 0 2 1 0]
            [0 0 0 0 2 1]
            [0 0 0 0 0 2]
            [2 0 0 0 0 0]
            [1 2 0 0 0 0]
            [0 1 2 0 0 0]

        We now show the common neighbors matrix which includes adjacent
        vertices ::

            sage: G1.common_neighbors_matrix(nonedgesonly = False)
            [0 1 1 2 1 0]
            [1 0 2 1 2 1]
            [1 2 0 2 1 2]
            [2 1 2 0 2 1]
            [1 2 1 2 0 1]
            [0 1 2 1 1 0]

        The common neighbors matrix  for a fan on 6 vertices counting only
        non-adjacent vertex pairs ::

            sage: H = Graph([(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,2),(2,3),(3,4),(4,5)])
            sage: H.common_neighbors_matrix()
            [0 0 0 0 0 0 0]
            [0 0 0 2 1 1 1]
            [0 0 0 0 2 1 1]
            [0 2 0 0 0 2 1]
            [0 1 2 0 0 0 1]
            [0 1 1 2 0 0 1]
            [0 1 1 1 1 1 0]

        It is an error to input anything other than a simple graph::

            sage: G = Graph([(0,0)],loops=True)
            sage: G.common_neighbors_matrix()
            Traceback (most recent call last):
            ...
            ValueError: This method is not known to work on graphs with loops.
            Perhaps this method can be updated to handle them, but in the
            meantime if you want to use it please disallow loops using
            allow_loops().

        .. SEEALSO::

            * :meth:`most_common_neighbors` --
              returns node pairs with most shared neighbors

        TESTS::

            sage: G = graphs.CompleteGraph(4)
            sage: M = G.common_neighbors_matrix()
            sage: M.is_zero()
            True
            sage: Graph(1).common_neighbors_matrix()
            [0]
            sage: Graph().common_neighbors_matrix()
            []
            sage: G = Graph([(0,1),(1,2),(2,3),(3,0),(0,2)])
            sage: G.common_neighbors_matrix()
            [0 0 0 0]
            [0 0 0 2]
            [0 0 0 0]
            [0 2 0 0]
        """
        self._scream_if_not_simple()
        if vertices is None:
            vertices = self.vertices()
        A = self.adjacency_matrix(vertices=vertices)
        M = A**2
        for v in range(self.order()):
            M[v, v] = 0
            if nonedgesonly:
                for w in range(v + 1, self.order()):
                    if A[v, w]:
                        M[v, w] = M[w, v] = 0
        return M

    @doc_index("Leftovers")
    def most_common_neighbors(self, nonedgesonly=True):
        r"""
        Return vertex pairs with maximal number of common neighbors.

        This method is only valid for simple (no loops, no multiple edges)
        graphs with order `\geq 2`

        INPUT:

        - ``nonedgesonly``-- boolean (default: ``True``); if ``True``, assigns
          `0` value to adjacent vertices.

        OUTPUT: list of tuples of edge pairs

        EXAMPLES:

        The maximum common neighbor (non-adjacent) pairs for a straight
        linear 2-tree ::

            sage: G1 = Graph([(0,1),(0,2),(1,2),(1,3),(3,5),(2,4),(2,3),(3,4),(4,5)])
            sage: G1.most_common_neighbors()
            [(0, 3), (1, 4), (2, 5)]

        If we include non-adjacent pairs ::

            sage: G1.most_common_neighbors(nonedgesonly = False)
            [(0, 3), (1, 2), (1, 4), (2, 3), (2, 5), (3, 4)]

        The common neighbors matrix  for a fan on 6 vertices counting only
        non-adjacent vertex pairs ::

            sage: H = Graph([(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,2),(2,3),(3,4),(4,5)])
            sage: H.most_common_neighbors()
            [(1, 3), (2, 4), (3, 5)]

        .. SEEALSO::

            * :meth:`common_neighbors_matrix` --
              a similar method giving a matrix of number of common neighbors

        TESTS::

            sage: G=graphs.CompleteGraph(4)
            sage: G.most_common_neighbors()
            []
            sage: G.most_common_neighbors(nonedgesonly=False)
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: Graph(1).most_common_neighbors()
            Traceback (most recent call last):
            ...
            ValueError: this method is defined for graphs with at least 2 vertices
            sage: Graph().most_common_neighbors()
            Traceback (most recent call last):
            ...
            ValueError: this method is defined for graphs with at least 2 vertices
            sage: G = Graph([(0,1),(1,2),(2,3),(3,0),(0,2)])
            sage: G.most_common_neighbors()
            [(1, 3)]
            sage: G.most_common_neighbors(nonedgesonly=False)
            [(0, 2), (1, 3)]
        """
        self._scream_if_not_simple()
        if self.num_verts() < 2:
            raise ValueError('this method is defined for graphs with at least 2 vertices')
        verts = list(self)
        M = self.common_neighbors_matrix(vertices=verts, nonedgesonly=nonedgesonly)
        output = []
        coefficients = M.coefficients()
        if coefficients:
            maximum = max(coefficients)
            for v in range(self.num_verts()):
                for w in range(v + 1, self.num_verts()):
                    if M[v, w] == maximum:
                        output.append((verts[v], verts[w]))
        return output

    @doc_index("Leftovers")
    def arboricity(self, certificate=False):
        r"""
        Return the arboricity of the graph and an optional certificate.

        The arboricity is the minimum number of forests that covers the
        graph.

        See :wikipedia:`Arboricity`

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); whether to return
          a certificate.

        OUTPUT:

        When ``certificate = True``, then the function returns `(a, F)`
        where `a` is the arboricity and `F` is a list of `a` disjoint forests
        that partitions the edge set of `g`. The forests are represented as
        subgraphs of the original graph.

        If ``certificate = False``, the function returns just a integer
        indicating the arboricity.

        ALGORITHM:

        Represent the graph as a graphical matroid, then apply matroid
        :meth:`sage.matroid.partition` algorithm from the matroids module.

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: a,F = G.arboricity(True)
            sage: a
            2
            sage: all([f.is_forest() for f in F])
            True
            sage: len(set.union(*[set(f.edges()) for f in F])) == G.size()
            True

        TESTS::

            sage: g = Graph()
            sage: g.arboricity(True)
            (0, [])
        """
        from sage.matroids.constructor import Matroid
        P = Matroid(self).partition()
        if certificate:
            return (len(P), [self.subgraph(edges=forest) for forest in P])
        else:
            return len(P)

    @doc_index("Graph properties")
    def is_antipodal(self):
        r"""
        Check whether this graph is antipodal.

        A graph `G` of diameter `d` is said to be antipodal if its distance-`d`
        graph is a disjoint union of cliques.

        EXAMPLES::

            sage: G = graphs.JohnsonGraph(10, 5)
            sage: G.is_antipodal()
            True
            sage: H = G.folded_graph()
            sage: H.is_antipodal()
            False

        REFERENCES:

        See [BCN1989]_ p. 438 or [Sam2012]_ for this definition of antipodal
        graphs.

        TESTS::

            sage: G = graphs.PetersenGraph()
            sage: G.is_antipodal()
            False
            sage: G = graphs.HammingGraph(7, 2)
            sage: G.is_antipodal()
            True
            sage: G = Graph([(0,1), (2, 3)])
            sage: G.is_antipodal()
            False
            sage: G = Graph(4)
            sage: G.is_antipodal()
            True
            sage: graphs.CompleteGraph(5).is_antipodal()
            True
            sage: G = Graph()
            sage: G.is_antipodal()
            Traceback (most recent call last):
            ...
            ValueError: diameter is not defined for the empty graph
            sage: G = Graph(1)
            sage: G.is_antipodal()
            True
        """
        G = self.antipodal_graph()

        vertexSet = set(G)
        while vertexSet:
            v = vertexSet.pop()

            # all neighbours of v should be in the same clique as v
            clique = set(G.neighbor_iterator(v, closed=True))
            for u in clique:
                if set(G.neighbor_iterator(u, closed=True)) != clique:
                    return False

            vertexSet.difference_update(clique)

        return True

    @doc_index("Leftovers")
    def folded_graph(self, check=False):
        r"""
        Return the antipodal fold of this graph.

        Given an antipodal graph `G` let `G_d` be its distance-`d` graph.
        Then the folded graph of `G` has a vertex for each maximal clique
        of `G_d` and two cliques are adjacent if there is an edge in `G`
        connecting the two.

        .. SEEALSO::

            :meth:`sage.graphs.graph.is_antipodal`

        INPUT:

        - ``check`` -- boolean (default: ``False``); whether to check if the
          graph is antipodal. If ``check`` is ``True`` and the graph is not
          antipodal, then return ``False``.

        OUTPUT:

        This function returns a new graph and ``self`` is not touched.

        .. NOTE::

            The input is expected to be an antipodal graph.
            You can check that a graph is antipodal using
            :meth:`sage.graphs.graph.is_antipodal`.

        EXAMPLES::

            sage: G = graphs.JohnsonGraph(10, 5)
            sage: H = G.folded_graph(); H
            Folded Johnson graph with parameters 10,5: Graph on 126 vertices
            sage: Gd = G.distance_graph(G.diameter())
            sage: all(i == 1 for i in Gd.degree())
            True
            sage: H.is_distance_regular(True)
            ([25, 16, None], [None, 1, 4])

        This method doesn't check if the graph is antipodal::

            sage: G = graphs.PetersenGraph()
            sage: G.is_antipodal()
            False
            sage: G.folded_graph()  # some garbage
            Folded Petersen graph: Graph on 2 vertices
            sage: G.folded_graph(check=True)
            False

        REFERENCES:

        See [BCN1989]_ p. 438 or [Sam2012]_ for this definition of folded graph.

        TESTS::

            sage: G = Graph(5)
            sage: G.folded_graph()
            Folded Graph: Graph on 1 vertex
            sage: G = graphs.CompleteGraph(5)
            sage: G.folded_graph()
            Folded Complete graph: Graph on 1 vertex
            sage: G = Graph()
            sage: G.folded_graph()
            Traceback (most recent call last):
            ...
            ValueError: diameter is not defined for the empty graph
            sage: G = Graph(1)
            sage: G.folded_graph()
            Folded Graph: Graph on 1 vertex
        """
        G = self.antipodal_graph()

        vertices = set(G)
        newVertices = []
        while vertices:
            v = vertices.pop()
            clique = frozenset(G.neighbor_iterator(v, closed=True))

            if check:
                for u in clique:
                    if frozenset(G.neighbor_iterator(u, closed=True)) != clique:
                        return False

            newVertices.append(clique)
            vertices.difference_update(clique)

        # now newVertices is a map {0, ..., numCliques-1} -> antipodal classes
        numCliques = len(newVertices)
        edges = []
        for i, j in itertools.combinations(range(numCliques), 2):
            if any(self.has_edge(u, v) for u, v in
                   itertools.product(newVertices[i], newVertices[j])):
                edges.append((i, j))

        H = Graph([range(numCliques), edges], format='vertices_and_edges')
        name = self.name() if self.name() != "" else "Graph"
        H.name(f"Folded {name}")
        return H

    @doc_index("Leftovers")
    def antipodal_graph(self):
        r"""
        Return the antipodal graph of ``self``.

        The antipodal graph of a graph `G` has the same vertex set of `G` and
        two vertices are adjacent if their distance in `G` is equal to the
        diameter of `G`.

        OUTPUT:

        A new graph. ``self`` is not touched.

        EXAMPLES::

            sage: G = graphs.JohnsonGraph(10, 5)
            sage: G.antipodal_graph()
            Antipodal graph of Johnson graph with parameters 10,5: Graph on 252 vertices
            sage: G = graphs.HammingGraph(8, 2)
            sage: G.antipodal_graph()
            Antipodal graph of Hamming Graph with parameters 8,2: Graph on 256 vertices

        The antipodal graph of a disconnected graph is its complement::

            sage: G = Graph(5)
            sage: H = G.antipodal_graph()
            sage: H.is_isomorphic(G.complement())
            True

        TESTS::

            sage: G = Graph([(0, 1), (2, 3)])
            sage: H = G.antipodal_graph()
            sage: H.is_isomorphic(Graph([(0, 2), (0, 3), (1, 2), (1, 3)]))
            True
            sage: G = Graph()
            sage: G.antipodal_graph()
            Traceback (most recent call last):
            ...
            ValueError: diameter is not defined for the empty graph
            sage: G = Graph(1)
            sage: G.antipodal_graph()
            Antipodal graph of Graph: Looped graph on 1 vertex
        """
        H = self.distance_graph(self.diameter())

        name = self.name() if self.name() != "" else "Graph"
        H.name(f"Antipodal graph of {name}")
        return H

    @doc_index("Basic methods")
    def bipartite_double(self, extended=False):
        r"""
        Return the (extended) bipartite double of this graph.

        The bipartite double of a graph `G` has vertex set
        `\{ (v,0), (v,1) : v \in G\}` and for any edge `(u, v)` in `G`
        it has edges `((u,0),(v,1))` and `((u,1),(v,0))`.
        Note that this is the tensor product of `G` with `K_2`.

        The extended bipartite double of `G` is the bipartite double of
        `G` after added all edges `((v,0),(v,1))` for all vertices `v`.

        INPUT:

        - ``extended`` -- boolean (default: ``False``); Whether to return the
          extended bipartite double, or only the bipartite double (default)

        OUTPUT:

        A graph; ``self`` is left untouched.

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: H = G.bipartite_double()
            sage: G == graphs.PetersenGraph()  # G is left invariant
            True
            sage: H.order() == 2 * G.order()
            True
            sage: H.size() == 2 * G.size()
            True
            sage: H.is_bipartite()
            True
            sage: H.bipartite_sets() == (set([(v, 0) for v in G]),
            ....: set([(v, 1) for v in G]))
            True
            sage: H.is_isomorphic(G.tensor_product(graphs.CompleteGraph(2)))
            True

        Behaviour with disconnected graphs::

            sage: G1 = graphs.PetersenGraph()
            sage: G2 = graphs.HoffmanGraph()
            sage: G = G1.disjoint_union(G2)
            sage: H = G.bipartite_double()
            sage: H1 = G1.bipartite_double()
            sage: H2 = G2.bipartite_double()
            sage: H.is_isomorphic(H1.disjoint_union(H2))
            True

        .. SEEALSO::

            :wikipedia:`Bipartite_double_cover`,
            `WolframAlpha Bipartite Double
            <https://mathworld.wolfram.com/BipartiteDoubleGraph.html>`_,
            [VDKT2016]_ p. 20 for the extended bipartite double.

        TESTS::

            sage: G = graphs.PetersenGraph()
            sage: H = G.bipartite_double(True)
            sage: G == graphs.PetersenGraph()  # G is left invariant
            True
            sage: H.order() == 2 * G.order()
            True
            sage: H.size() == 2 * G.size() + G.order()
            True
            sage: H.is_bipartite()
            True
            sage: H.bipartite_sets() == (set([(v, 0) for v in G]),
            ....: set([(v, 1) for v in G]))
            True
            sage: H.is_isomorphic(G.tensor_product(graphs.CompleteGraph(2)))
            False

        Test edge cases::

            sage: G = Graph()
            sage: H = G.bipartite_double()
            sage: H.size() + H.order()
            0
            sage: H = G.bipartite_double(True)
            sage: H.size() + H.order()
            0
            sage: G = Graph(1)
            sage: H = G.bipartite_double()
            sage: H.size() == 0 and H.order() == 2
            True
            sage: H = G.bipartite_double(True)
            sage: H.is_isomorphic(Graph([(0, 1)]))
            True
        """
        G = self.tensor_product(Graph([(0, 1)]))

        if extended:
            G.add_edges(((v, 0), (v, 1)) for v in self)

        prefix = "Extended " if extended else ""
        G.name("%sBipartite Double of %s"%(prefix, self.name()))
        return G

    # Aliases to functions defined in other modules
    from sage.graphs.weakly_chordal import is_long_hole_free, is_long_antihole_free, is_weakly_chordal
    from sage.graphs.asteroidal_triples import is_asteroidal_triple_free
    from sage.graphs.chrompoly import chromatic_polynomial
    from sage.graphs.graph_decompositions.rankwidth import rank_decomposition
    from sage.graphs.graph_decompositions.tree_decomposition import treewidth
    from sage.graphs.graph_decompositions.vertex_separation import pathwidth
    from sage.graphs.graph_decompositions.tree_decomposition import treelength
    from sage.graphs.graph_decompositions.clique_separators import atoms_and_clique_separators
    from sage.graphs.matchpoly import matching_polynomial
    from sage.graphs.cliquer import all_max_clique as cliques_maximum
    from sage.graphs.cliquer import all_cliques
    from sage.graphs.spanning_tree import random_spanning_tree
    from sage.graphs.spanning_tree import spanning_trees
    from sage.graphs.graph_decompositions.graph_products import is_cartesian_product
    from sage.graphs.distances_all_pairs import is_distance_regular
    from sage.graphs.base.static_dense_graph import is_strongly_regular
    from sage.graphs.line_graph import is_line_graph
    from sage.graphs.tutte_polynomial import tutte_polynomial
    from sage.graphs.lovasz_theta import lovasz_theta
    from sage.graphs.partial_cube import is_partial_cube
    from sage.graphs.orientations import strong_orientations_iterator, random_orientation
    from sage.graphs.connectivity import bridges, cleave, spqr_tree
    from sage.graphs.connectivity import is_triconnected
    from sage.graphs.comparability import is_comparability
    from sage.graphs.comparability import is_permutation
    from sage.graphs.convexity_properties import geodetic_closure
    from sage.graphs.domination import is_dominating
    from sage.graphs.domination import is_redundant
    from sage.graphs.domination import private_neighbors
    from sage.graphs.domination import minimal_dominating_sets
    from sage.graphs.traversals import (lex_M, maximum_cardinality_search,
                                        maximum_cardinality_search_M)
    from sage.graphs.isoperimetric_inequalities import cheeger_constant, edge_isoperimetric_number, vertex_isoperimetric_number
    from sage.graphs.graph_coloring import fractional_chromatic_number
    from sage.graphs.graph_coloring import fractional_chromatic_index

_additional_categories = {
    "is_long_hole_free"         : "Graph properties",
    "is_long_antihole_free"     : "Graph properties",
    "is_weakly_chordal"         : "Graph properties",
    "is_asteroidal_triple_free" : "Graph properties",
    "chromatic_polynomial"      : "Coloring",
    "rank_decomposition"        : "Algorithmically hard stuff",
    "treewidth"                 : "Algorithmically hard stuff",
    "pathwidth"                 : "Algorithmically hard stuff",
    "treelength"                : "Algorithmically hard stuff",
    "matching_polynomial"       : "Algorithmically hard stuff",
    "all_max_clique"            : "Clique-related methods",
    "cliques_maximum"           : "Clique-related methods",
    "all_cliques"               : "Clique-related methods",
    "atoms_and_clique_separators" : "Clique-related methods",
    "random_spanning_tree"      : "Connectivity, orientations, trees",
    "spanning_trees"            : "Connectivity, orientations, trees",
    "is_cartesian_product"      : "Graph properties",
    "is_distance_regular"       : "Graph properties",
    "is_strongly_regular"       : "Graph properties",
    "is_line_graph"             : "Graph properties",
    "is_partial_cube"           : "Graph properties",
    "is_comparability"          : "Graph properties",
    "is_permutation"            : "Graph properties",
    "tutte_polynomial"          : "Algorithmically hard stuff",
    "lovasz_theta"              : "Leftovers",
    "strong_orientations_iterator" : "Connectivity, orientations, trees",
    "random_orientation"        : "Connectivity, orientations, trees",
    "bridges"                   : "Connectivity, orientations, trees",
    "cleave"                    : "Connectivity, orientations, trees",
    "spqr_tree"                 : "Connectivity, orientations, trees",
    "is_triconnected"           : "Connectivity, orientations, trees",
    "is_dominating"             : "Domination",
    "is_redundant"              : "Domination",
    "private_neighbors"         : "Domination",
    "minimal_dominating_sets"   : "Domination",
    "lex_M"                     : "Traversals",
    "maximum_cardinality_search" : "Traversals",
    "maximum_cardinality_search_M" : "Traversals",
    "cheeger_constant"          : "Expansion properties",
    "edge_isoperimetric_number" : "Expansion properties",
    "vertex_isoperimetric_number" : "Expansion properties",
    "fractional_chromatic_number" : "Coloring",
    "fractional_chromatic_index" : "Coloring",
    "geodetic_closure"          : "Leftovers"
    }

__doc__ = __doc__.replace("{INDEX_OF_METHODS}",gen_thematic_rest_table_index(Graph,_additional_categories))
