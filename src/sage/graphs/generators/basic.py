# -*- coding: utf-8 -*-
r"""
Basic Graphs

The methods defined here appear in :mod:`sage.graphs.graph_generators`.

"""
###########################################################################
#
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###########################################################################
from __future__ import print_function

# import from Sage library
from sage.graphs.graph import Graph
from sage.graphs import graph
from math import sin, cos, pi

def BullGraph():
    r"""
    Returns a bull graph with 5 nodes.

    A bull graph is named for its shape. It's a triangle with horns.
    For more information, see this
    :wikipedia:`Wikipedia article on the bull graph <Bull_graph>`.

    PLOTTING:

    Upon construction, the position dictionary is filled to override the
    spring-layout algorithm. By convention, the bull graph is drawn as a
    triangle with the first node (0) on the bottom. The second and third nodes
    (1 and 2) complete the triangle. Node 3 is the horn connected to 1 and node
    4 is the horn connected to node 2.


    EXAMPLES:

    Construct and show a bull graph::

        sage: g = graphs.BullGraph(); g
        Bull graph: Graph on 5 vertices
        sage: g.show() # long time

    The bull graph has 5 vertices and 5 edges. Its radius is 2, its
    diameter 3, and its girth 3. The bull graph is planar with chromatic
    number 3 and chromatic index also 3. ::

        sage: g.order(); g.size()
        5
        5
        sage: g.radius(); g.diameter(); g.girth()
        2
        3
        3
        sage: g.chromatic_number()
        3

    The bull graph has chromatic polynomial `x(x - 2)(x - 1)^3` and
    Tutte polynomial `x^4 + x^3 + x^2 y`. Its characteristic polynomial
    is `x(x^2 - x - 3)(x^2 + x - 1)`, which follows from the definition of
    characteristic polynomials for graphs, i.e. `\det(xI - A)`, where
    `x` is a variable, `A` the adjacency matrix of the graph, and `I`
    the identity matrix of the same dimensions as `A`. ::

        sage: chrompoly = g.chromatic_polynomial()
        sage: bool(expand(x * (x - 2) * (x - 1)^3) == chrompoly)
        True
        sage: charpoly = g.characteristic_polynomial()
        sage: M = g.adjacency_matrix(); M
        [0 1 1 0 0]
        [1 0 1 1 0]
        [1 1 0 0 1]
        [0 1 0 0 0]
        [0 0 1 0 0]
        sage: Id = identity_matrix(ZZ, M.nrows())
        sage: D = x*Id - M
        sage: bool(D.determinant() == charpoly)
        True
        sage: bool(expand(x * (x^2 - x - 3) * (x^2 + x - 1)) == charpoly)
        True
    """
    edge_list = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 4)]
    pos_dict = {0:(0,0), 1:(-1,1), 2:(1,1), 3:(-2,2), 4:(2,2)}
    return graph.Graph(edge_list, pos=pos_dict, name="Bull graph")

def ButterflyGraph():
    r"""
    Returns the butterfly graph.

    Let `C_3` be the cycle graph on 3 vertices. The butterfly or bowtie
    graph is obtained by joining two copies of `C_3` at a common vertex,
    resulting in a graph that is isomorphic to the friendship graph `F_2`.
    For more information, see this
    `Wikipedia article on the butterfly graph <http://en.wikipedia.org/wiki/Butterfly_graph>`_.

    .. SEEALSO::

        - :meth:`GraphGenerators.FriendshipGraph`

    EXAMPLES:

    The butterfly graph is a planar graph on 5 vertices and having
    6 edges. ::

        sage: G = graphs.ButterflyGraph(); G
        Butterfly graph: Graph on 5 vertices
        sage: G.show()  # long time
        sage: G.is_planar()
        True
        sage: G.order()
        5
        sage: G.size()
        6

    It has diameter 2, girth 3, and radius 1. ::

        sage: G.diameter()
        2
        sage: G.girth()
        3
        sage: G.radius()
        1

    The butterfly graph is Eulerian, with chromatic number 3. ::

        sage: G.is_eulerian()
        True
        sage: G.chromatic_number()
        3
    """
    edge_dict = {
        0: [3,4],
        1: [2,4],
        2: [4],
        3: [4]}
    pos_dict = {
        0: [-1, 1],
        1: [1, 1],
        2: [1, -1],
        3: [-1, -1],
        4: [0, 0]}
    return graph.Graph(edge_dict, pos=pos_dict, name="Butterfly graph")

def CircularLadderGraph(n):
    """
    Returns a circular ladder graph with 2\*n nodes.

    A Circular ladder graph is a ladder graph that is connected at the
    ends, i.e.: a ladder bent around so that top meets bottom. Thus it
    can be described as two parallel cycle graphs connected at each
    corresponding node pair.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the circular
    ladder graph is displayed as an inner and outer cycle pair, with
    the first n nodes drawn on the inner circle. The first (0) node is
    drawn at the top of the inner-circle, moving clockwise after that.
    The outer circle is drawn with the (n+1)th node at the top, then
    counterclockwise as well.

    EXAMPLES: Construct and show a circular ladder graph with 26 nodes

    ::

        sage: g = graphs.CircularLadderGraph(13)
        sage: g.show() # long time

    Create several circular ladder graphs in a Sage graphics array

    ::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:    k = graphs.CircularLadderGraph(i+3)
        ....:    g.append(k)
        sage: for i in range(3):
        ....:    n = []
        ....:    for m in range(3):
        ....:        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:    j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time
    """
    pos_dict = {}
    for i in range(n):
        x = float(cos((pi/2) + ((2*pi)/n)*i))
        y = float(sin((pi/2) + ((2*pi)/n)*i))
        pos_dict[i] = [x,y]
    for i in range(n,2*n):
        x = float(2*(cos((pi/2) + ((2*pi)/n)*(i-n))))
        y = float(2*(sin((pi/2) + ((2*pi)/n)*(i-n))))
        pos_dict[i] = (x,y)

    G = Graph(pos=pos_dict, name="Circular Ladder graph")
    G.add_vertices( range(2*n) )
    G.add_cycle( range(n) )
    G.add_cycle( range(n,2*n) )
    G.add_edges( (i,i+n) for i in range(n) )
    return G

def ClawGraph():
    """
    Returns a claw graph.

    A claw graph is named for its shape. It is actually a complete
    bipartite graph with (n1, n2) = (1, 3).

    PLOTTING: See CompleteBipartiteGraph.

    EXAMPLES: Show a Claw graph

    ::

        sage: (graphs.ClawGraph()).show() # long time

    Inspect a Claw graph

    ::

        sage: G = graphs.ClawGraph()
        sage: G
        Claw graph: Graph on 4 vertices
    """
    edge_list = [(0, 1), (0, 2), (0, 3)]
    pos_dict = {0:(0,1),1:(-1,0),2:(0,0),3:(1,0)}
    return graph.Graph(edge_list, pos=pos_dict, name="Claw graph")

def CycleGraph(n):
    r"""
    Returns a cycle graph with n nodes.

    A cycle graph is a basic structure which is also typically called
    an n-gon.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, each cycle
    graph will be displayed with the first (0) node at the top, with
    the rest following in a counterclockwise manner.

    The cycle graph is a good opportunity to compare efficiency of
    filling a position dictionary vs. using the spring-layout algorithm
    for plotting. Because the cycle graph is very symmetric, the
    resulting plots should be similar (in cases of small n).

    Filling the position dictionary in advance adds O(n) to the
    constructor.

    EXAMPLES: Compare plotting using the predefined layout and
    networkx::

        sage: import networkx
        sage: n = networkx.cycle_graph(23)
        sage: spring23 = Graph(n)
        sage: posdict23 = graphs.CycleGraph(23)
        sage: spring23.show() # long time
        sage: posdict23.show() # long time

    We next view many cycle graphs as a Sage graphics array. First we
    use the ``CycleGraph`` constructor, which fills in the
    position dictionary::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.CycleGraph(i+3)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time

    Compare to plotting with the spring-layout algorithm::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     spr = networkx.cycle_graph(i+3)
        ....:     k = Graph(spr)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time
    """
    pos_dict = {}
    for i in range(n):
        x = float(cos((pi/2) + ((2*pi)/n)*i))
        y = float(sin((pi/2) + ((2*pi)/n)*i))
        pos_dict[i] = (x,y)
    G = graph.Graph(n,pos=pos_dict, name="Cycle graph")
    G.add_cycle(range(n))
    return G

def CompleteGraph(n):
    """
    Returns a complete graph on n nodes.

    A Complete Graph is a graph in which all nodes are connected to all
    other nodes.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, each complete
    graph will be displayed with the first (0) node at the top, with
    the rest following in a counterclockwise manner.

    In the complete graph, there is a big difference visually in using
    the spring-layout algorithm vs. the position dictionary used in
    this constructor. The position dictionary flattens the graph,
    making it clear which nodes an edge is connected to. But the
    complete graph offers a good example of how the spring-layout
    works. The edges push outward (everything is connected), causing
    the graph to appear as a 3-dimensional pointy ball. (See examples
    below).

    EXAMPLES: We view many Complete graphs with a Sage Graphics Array,
    first with this constructor (i.e., the position dictionary
    filled)::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.CompleteGraph(i+3)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time

    We compare to plotting with the spring-layout algorithm::

        sage: import networkx
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     spr = networkx.complete_graph(i+3)
        ....:     k = Graph(spr)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time

    Compare the constructors (results will vary)

    ::

        sage: import networkx
        sage: t = cputime()
        sage: n = networkx.complete_graph(389); spring389 = Graph(n)
        sage: cputime(t)           # random
        0.59203700000000126
        sage: t = cputime()
        sage: posdict389 = graphs.CompleteGraph(389)
        sage: cputime(t)           # random
        0.6680419999999998

    We compare plotting::

        sage: import networkx
        sage: n = networkx.complete_graph(23)
        sage: spring23 = Graph(n)
        sage: posdict23 = graphs.CompleteGraph(23)
        sage: spring23.show() # long time
        sage: posdict23.show() # long time
    """
    pos_dict = {}
    for i in range(n):
        x = float(cos((pi/2) + ((2*pi)/n)*i))
        y = float(sin((pi/2) + ((2*pi)/n)*i))
        pos_dict[i] = (x,y)
    G = graph.Graph(n,pos=pos_dict, name="Complete graph")
    G.add_edges(((i,j) for i in range(n) for j in range(i+1,n)))
    return G

def CompleteBipartiteGraph(n1, n2):
    """
    Returns a Complete Bipartite Graph sized n1+n2, with each of the
    nodes [0,(n1-1)] connected to each of the nodes [n1,(n2-1)] and
    vice versa.

    A Complete Bipartite Graph is a graph with its vertices partitioned
    into two groups, V1 and V2. Each v in V1 is connected to every v in
    V2, and vice versa.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, each complete
    bipartite graph will be displayed with the first n1 nodes on the
    top row (at y=1) from left to right. The remaining n2 nodes appear
    at y=0, also from left to right. The shorter row (partition with
    fewer nodes) is stretched to the same length as the longer row,
    unless the shorter row has 1 node; in which case it is centered.
    The x values in the plot are in domain [0,maxn1,n2].

    In the Complete Bipartite graph, there is a visual difference in
    using the spring-layout algorithm vs. the position dictionary used
    in this constructor. The position dictionary flattens the graph and
    separates the partitioned nodes, making it clear which nodes an
    edge is connected to. The Complete Bipartite graph plotted with the
    spring-layout algorithm tends to center the nodes in n1 (see
    spring_med in examples below), thus overlapping its nodes and
    edges, making it typically hard to decipher.

    Filling the position dictionary in advance adds O(n) to the
    constructor. Feel free to race the constructors below in the
    examples section. The much larger difference is the time added by
    the spring-layout algorithm when plotting. (Also shown in the
    example below). The spring model is typically described as
    `O(n^3)`, as appears to be the case in the NetworkX source
    code.

    EXAMPLES: Two ways of constructing the complete bipartite graph,
    using different layout algorithms::

        sage: import networkx
        sage: n = networkx.complete_bipartite_graph(389,157); spring_big = Graph(n)   # long time
        sage: posdict_big = graphs.CompleteBipartiteGraph(389,157)                    # long time

    Compare the plotting::

        sage: n = networkx.complete_bipartite_graph(11,17)
        sage: spring_med = Graph(n)
        sage: posdict_med = graphs.CompleteBipartiteGraph(11,17)

    Notice here how the spring-layout tends to center the nodes of n1

    ::

        sage: spring_med.show() # long time
        sage: posdict_med.show() # long time

    View many complete bipartite graphs with a Sage Graphics Array,
    with this constructor (i.e., the position dictionary filled)::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.CompleteBipartiteGraph(i+1,4)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time

    We compare to plotting with the spring-layout algorithm::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     spr = networkx.complete_bipartite_graph(i+1,4)
        ....:     k = Graph(spr)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time

    :trac:`12155`::

        sage: graphs.CompleteBipartiteGraph(5,6).complement()
        complement(Complete bipartite graph): Graph on 11 vertices

    TESTS:

    Prevent negative dimensions (:trac:`18530`)::

        sage: graphs.CompleteBipartiteGraph(-1,1)
        Traceback (most recent call last):
        ...
        ValueError: The arguments n1(=-1) and n2(=1) must be positive integers.
        sage: graphs.CompleteBipartiteGraph(1,-1)
        Traceback (most recent call last):
        ...
        ValueError: The arguments n1(=1) and n2(=-1) must be positive integers.
    """
    if n1<0 or n2<0:
        raise ValueError('The arguments n1(={}) and n2(={}) must be positive integers.'.format(n1,n2))

    pos_dict = {}
    c1 = 1 # scaling factor for top row
    c2 = 1 # scaling factor for bottom row
    c3 = 0 # pad to center if top row has 1 node
    c4 = 0 # pad to center if bottom row has 1 node
    if n1 > n2:
        if n2 == 1:
            c4 = (n1-1)/2
        else:
            c2 = ((n1-1)/(n2-1))
    elif n2 > n1:
        if n1 == 1:
            c3 = (n2-1)/2
        else:
            c1 = ((n2-1)/(n1-1))
    for i in range(n1):
        x = c1*i + c3
        y = 1
        pos_dict[i] = (x,y)
    for i in range(n1+n2)[n1:]:
        x = c2*(i-n1) + c4
        y = 0
        pos_dict[i] = (x,y)

    G = Graph(n1+n2, pos=pos_dict, name="Complete bipartite graph")
    G.add_edges((i,j) for i in range(n1) for j in range(n1,n1+n2))
    return G

def CompleteMultipartiteGraph(l):
    r"""
    Returns a complete multipartite graph.

    INPUT:

    - ``l`` -- a list of integers : the respective sizes
      of the components.

    EXAMPLE:

    A complete tripartite graph with sets of sizes
    `5, 6, 8`::

        sage: g = graphs.CompleteMultipartiteGraph([5, 6, 8]); g
        Multipartite Graph with set sizes [5, 6, 8]: Graph on 19 vertices

    It clearly has a chromatic number of 3::

        sage: g.chromatic_number()
        3
    """
    
    n = sum(l) #getting the number of vertices
    r = len(l) #getting the number of partitions
    positions = {}

    if r > 2: #position code gives bad results on bipartite or isolated graphs

        '''
        Produce a layout of the vertices so that vertices in the same
        vertex set are adjacent and clearly separated from vertices in other
        vertex sets.

        This is done by calculating the vertices of an r-gon then
        calculating the slope between adjacent vertices. We then 'walk'
        around the r-gon placing graph vertices in regular intervals between 
        adjacent vertices of the r-gon.

        Makes a nicely organized graph like in this picture: 
        https://commons.wikimedia.org/wiki/File:Turan_13-4.svg
        '''

        points = [[cos(2*pi*i/r),sin(2*pi*i/r)] for i in range(r)]
        slopes = [[(points[(i+1)%r][0]-points[i%r][0]),
                   (points[(i+1)%r][1]-points[i%r][1])] for i in range(r)]

        counter = 0

        for i in range(len(l)):
            vertex_set_size = l[i]+1
            for j in range(1,vertex_set_size):
                x = points[i][0]+slopes[i][0]*j/(vertex_set_size)
                y = points[i][1]+slopes[i][1]*j/(vertex_set_size)
                positions[counter] = (x,y)
                counter += 1

    g = Graph()
    for i in l:
        g = g + CompleteGraph(i)

    g = g.complement()
    g.set_pos(positions)
    g.name("Multipartite Graph with set sizes "+str(l))



    return g

def DiamondGraph():
    """
    Returns a diamond graph with 4 nodes.

    A diamond graph is a square with one pair of diagonal nodes
    connected.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the diamond
    graph is drawn as a diamond, with the first node on top, second on
    the left, third on the right, and fourth on the bottom; with the
    second and third node connected.

    EXAMPLES: Construct and show a diamond graph

    ::

        sage: g = graphs.DiamondGraph()
        sage: g.show() # long time
    """
    pos_dict = {0:(0,1),1:(-1,0),2:(1,0),3:(0,-1)}
    edges = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
    return graph.Graph(edges, pos=pos_dict, name="Diamond Graph")

def EmptyGraph():
    """
    Returns an empty graph (0 nodes and 0 edges).

    This is useful for constructing graphs by adding edges and vertices
    individually or in a loop.

    PLOTTING: When plotting, this graph will use the default
    spring-layout algorithm, unless a position dictionary is
    specified.

    EXAMPLES: Add one vertex to an empty graph and then show::

        sage: empty1 = graphs.EmptyGraph()
        sage: empty1.add_vertex()
        0
        sage: empty1.show() # long time

    Use for loops to build a graph from an empty graph::

        sage: empty2 = graphs.EmptyGraph()
        sage: for i in range(5):
        ....:     empty2.add_vertex() # add 5 nodes, labeled 0-4
        0
        1
        2
        3
        4
        sage: for i in range(3):
        ....:     empty2.add_edge(i,i+1) # add edges {[0:1],[1:2],[2:3]}
        sage: for i in range(4)[1:]:
        ....:     empty2.add_edge(4,i) # add edges {[1:4],[2:4],[3:4]}
        sage: empty2.show() # long time
    """
    return graph.Graph(sparse=True)

def ToroidalGrid2dGraph(n1, n2):
    r"""
    Returns a toroidal 2-dimensional grid graph with `n_1n_2` nodes (`n_1`
    rows and `n_2` columns).

    The toroidal 2-dimensional grid with parameters `n_1,n_2` is the
    2-dimensional grid graph with identical parameters to which are added
    the edges `((i,0),(i,n_2-1))` and `((0,i),(n_1-1,i))`.

    EXAMPLE:

    The toroidal 2-dimensional grid is a regular graph, while the usual
    2-dimensional grid is not ::

        sage: tgrid = graphs.ToroidalGrid2dGraph(8,9)
        sage: print(tgrid)
        Toroidal 2D Grid Graph with parameters 8,9
        sage: grid = graphs.Grid2dGraph(8,9)
        sage: grid.is_regular()
        False
        sage: tgrid.is_regular()
        True
    """

    g = Grid2dGraph(n1,n2, set_positions=False)

    g.add_edges([((i,0),(i,n2-1)) for i in range(n1)] + [((0,i),(n1-1,i)) for i in range(n2)])

    g.name("Toroidal 2D Grid Graph with parameters "+str(n1)+","+str(n2))

    d = g.get_pos()
    n1 += 0.
    n2 += 0.
    uf = (n1/2)*(n1/2)
    vf = (n2/2)*(n2/2)
    for u,v in d:
        x,y = d[(u,v)]
        x +=  0.25*(1.0+u*(u-n1+1)/uf)
        y +=  0.25*(1+v*(v-n2+1)/vf)
        d[(u,v)] = (x,y)

    return g

def Toroidal6RegularGrid2dGraph(n1, n2):
    r"""
    Returns a toroidal 6-regular grid.

    The toroidal 6-regular grid is a 6-regular graph on `n_1\times n_2`
    vertices and its elements have coordinates `(i,j)` for `i \in \{0...i-1\}`
    and `j \in \{0...j-1\}`.

    Its edges are those of the :meth:`ToroidalGrid2dGraph`, to which are
    added the edges between `(i,j)` and `((i+1)\%n_1, (j+1)\%n_2)`.

    INPUT:

    - ``n1, n2`` (integers) -- see above.

    EXAMPLE:

    The toroidal 6-regular grid on `25` elements::

        sage: g = graphs.Toroidal6RegularGrid2dGraph(5,5)
        sage: g.is_regular(k=6)
        True
        sage: g.is_vertex_transitive()
        True
        sage: g.line_graph().is_vertex_transitive()
        True
        sage: g.automorphism_group().cardinality()
        300
        sage: g.is_hamiltonian()
        True

    TESTS:

    Senseless input::

        sage: graphs.Toroidal6RegularGrid2dGraph(5,2)
        Traceback (most recent call last):
        ...
        ValueError: Parameters n1 and n2 must be integers larger than 3 !
        sage: graphs.Toroidal6RegularGrid2dGraph(2,0)
        Traceback (most recent call last):
        ...
        ValueError: Parameters n1 and n2 must be integers larger than 3 !
    """

    if n1 <= 3 or n2 <= 3:
        raise ValueError("Parameters n1 and n2 must be integers larger than 3 !")

    g = ToroidalGrid2dGraph(n1,n2)
    for u,v in g:
        g.add_edge((u,v),((u+1)%n1,(v+1)%n2))

    g.name("Toroidal Hexagonal Grid graph on "+str(n1)+"x"+str(n2)+" elements")
    return g

def Grid2dGraph(n1, n2, set_positions=True):
    r"""
    Returns a `2`-dimensional grid graph with `n_1n_2` nodes (`n_1` rows and
    `n_2` columns).

    A 2d grid graph resembles a `2` dimensional grid. All inner nodes are
    connected to their `4` neighbors. Outer (non-corner) nodes are
    connected to their `3` neighbors. Corner nodes are connected to their
    2 neighbors.

    INPUT:

    - ``n1`` and ``n2`` -- two positive integers

    - ``set_positions`` -- (default: ``True``) boolean use to prevent setting
      the position of the nodes.
    
    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, nodes are
    labelled in (row, column) pairs with `(0, 0)` in the top left corner.
    Edges will always be horizontal and vertical - another advantage of
    filling the position dictionary.

    EXAMPLES: Construct and show a grid 2d graph Rows = `5`, Columns = `7`

    ::

        sage: g = graphs.Grid2dGraph(5,7)
        sage: g.show() # long time

    TESTS:

    Senseless input::

        sage: graphs.Grid2dGraph(5,0)
        Traceback (most recent call last):
        ...
        ValueError: Parameters n1 and n2 must be positive integers !
        sage: graphs.Grid2dGraph(-1,0)
        Traceback (most recent call last):
        ...
        ValueError: Parameters n1 and n2 must be positive integers !

    The graph name contains the dimension::

        sage: g = graphs.Grid2dGraph(5,7)
        sage: g.name()
        '2D Grid Graph for [5, 7]'
    """

    if n1 <= 0 or n2 <= 0:
        raise ValueError("Parameters n1 and n2 must be positive integers !")

    pos_dict = {}
    if set_positions:
        for i in range(n1):
            y = -i
            for j in range(n2):
                x = j
                pos_dict[i, j] = (x, y)

    G = graph.Graph(pos=pos_dict, name="2D Grid Graph for [{}, {}]".format(n1, n2))
    G.add_vertices( (i,j) for i in range(n1) for j in range(n2) )
    G.add_edges( ((i,j),(i+1,j)) for i in range(n1-1) for j in range(n2) )
    G.add_edges( ((i,j),(i,j+1)) for i in range(n1) for j in range(n2-1) )
    return G

def GridGraph(dim_list):
    """
    Returns an n-dimensional grid graph.

    INPUT:


    -  ``dim_list`` - a list of integers representing the
       number of nodes to extend in each dimension.


    PLOTTING: When plotting, this graph will use the default
    spring-layout algorithm, unless a position dictionary is
    specified.

    EXAMPLES::

        sage: G = graphs.GridGraph([2,3,4])
        sage: G.show()  # long time

    ::

        sage: C = graphs.CubeGraph(4)
        sage: G = graphs.GridGraph([2,2,2,2])
        sage: C.show()  # long time
        sage: G.show()  # long time

    TESTS:

    The graph name contains the dimension::

        sage: g = graphs.GridGraph([5, 7])
        sage: g.name()
        'Grid Graph for [5, 7]'
        sage: g = graphs.GridGraph([2, 3, 4])
        sage: g.name()
        'Grid Graph for [2, 3, 4]'
        sage: g = graphs.GridGraph([2, 4, 3])
        sage: g.name()
        'Grid Graph for [2, 4, 3]'

    One dimensional grids (i.e., path) have simple vertex labels::

        sage: g = graphs.GridGraph([5])
        sage: g.vertices()
        [0, 1, 2, 3, 4]

    The graph is correct::

        sage: dim = [randint(1,4) for i in range(4)]
        sage: g = graphs.GridGraph(dim)
        sage: import networkx
        sage: h = Graph( networkx.grid_graph(list(dim)) )
        sage: g.is_isomorphic(h)
        True

    Trivial cases::

        sage: g = graphs.GridGraph([]); g; g.vertices()
        Grid Graph for []: Graph on 0 vertices
        []
        sage: g = graphs.GridGraph([1]); g; g.vertices()
        Grid Graph for [1]: Graph on 1 vertex
        [0]
        sage: g = graphs.GridGraph([2]); g; g.vertices()
        Grid Graph for [2]: Graph on 2 vertices
        [0, 1]
        sage: g = graphs.GridGraph([1,1]); g; g.vertices()
        Grid Graph for [1, 1]: Graph on 1 vertex
        [(0, 0)]
        sage: g = graphs.GridGraph([1, 1, 1]); g; g.vertices()
        Grid Graph for [1, 1, 1]: Graph on 1 vertex
        [(0, 0, 0)]
        sage: g = graphs.GridGraph([1,1,2]); g; g.vertices()
        Grid Graph for [1, 1, 2]: Graph on 2 vertices
        [(0, 0, 0), (0, 0, 1)]

    All dimensions must be positive integers::

        sage: g = graphs.GridGraph([2,-1,3])
        Traceback (most recent call last):
        ...
        ValueError: All dimensions must be positive integers !
    """
    dim = [int(a) for a in dim_list]
    if any(a <= 0 for a in dim):
        raise ValueError("All dimensions must be positive integers !")

    g = Graph()
    n_dim = len(dim)
    if n_dim==1:
        # Vertices are labeled from 0 to dim[0]-1
        g = PathGraph(dim[0])
    elif n_dim==2:
        # We use the Grid2dGraph generator to also get the positions
        g = Grid2dGraph(*dim)
    elif n_dim>2:
        # Vertices are tuples of dimension n_dim, and the graph contains at
        # least vertex (0, 0, ..., 0)
        g.add_vertex(tuple([0]*n_dim))
        import itertools
        for u in itertools.product(*[range(d) for d in dim]):
            for i in range(n_dim):
                if u[i]+1<dim[i]:
                    v = list(u)
                    v[i] = u[i]+1
                    g.add_edge(u, tuple(v))

    g.name("Grid Graph for {}".format(dim))
    return g


def HouseGraph():
    """
    Returns a house graph with 5 nodes.

    A house graph is named for its shape. It is a triangle (roof) over a
    square (walls).

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the house
    graph is drawn with the first node in the lower-left corner of the
    house, the second in the lower-right corner of the house. The third
    node is in the upper-left corner connecting the roof to the wall,
    and the fourth is in the upper-right corner connecting the roof to
    the wall. The fifth node is the top of the roof, connected only to
    the third and fourth.

    EXAMPLES: Construct and show a house graph

    ::

        sage: g = graphs.HouseGraph()
        sage: g.show() # long time
    """
    pos_dict = {0:(-1,0),1:(1,0),2:(-1,1),3:(1,1),4:(0,2)}
    edges = [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 4)]
    return graph.Graph(edges, pos=pos_dict, name="House Graph")

def HouseXGraph():
    """
    Returns a house X graph with 5 nodes.

    A house X graph is a house graph with two additional edges. The
    upper-right corner is connected to the lower-left. And the
    upper-left corner is connected to the lower-right.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the house X
    graph is drawn with the first node in the lower-left corner of the
    house, the second in the lower-right corner of the house. The third
    node is in the upper-left corner connecting the roof to the wall,
    and the fourth is in the upper-right corner connecting the roof to
    the wall. The fifth node is the top of the roof, connected only to
    the third and fourth.

    EXAMPLES: Construct and show a house X graph

    ::

        sage: g = graphs.HouseXGraph()
        sage: g.show() # long time
    """
    pos_dict = {0:(-1,0),1:(1,0),2:(-1,1),3:(1,1),4:(0,2)}
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (2, 4), (3, 4)]
    return graph.Graph(edges, pos=pos_dict, name="House Graph")

def LadderGraph(n):
    """
    Returns a ladder graph with 2\*n nodes.

    A ladder graph is a basic structure that is typically displayed as
    a ladder, i.e.: two parallel path graphs connected at each
    corresponding node pair.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, each ladder
    graph will be displayed horizontally, with the first n nodes
    displayed left to right on the top horizontal line.

    EXAMPLES: Construct and show a ladder graph with 14 nodes

    ::

        sage: g = graphs.LadderGraph(7)
        sage: g.show() # long time

    Create several ladder graphs in a Sage graphics array

    ::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.LadderGraph(i+2)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time
    """
    pos_dict = {}
    for i in range(n):
        pos_dict[i] = (i,1)
    for i in range(n,2*n):
        x = i - n
        pos_dict[i] = (x,0)
    G = Graph(pos=pos_dict, name="Ladder graph")
    G.add_vertices( range(2*n) )
    G.add_path( range(n) )
    G.add_path( range(n,2*n) )
    G.add_edges( (i,i+n) for i in range(n) )
    return G

def LollipopGraph(n1, n2):
    """
    Returns a lollipop graph with n1+n2 nodes.

    A lollipop graph is a path graph (order n2) connected to a complete
    graph (order n1). (A barbell graph minus one of the bells).

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the complete
    graph will be drawn in the lower-left corner with the (n1)th node
    at a 45 degree angle above the right horizontal center of the
    complete graph, leading directly into the path graph.

    EXAMPLES: Construct and show a lollipop graph Candy = 13, Stick =
    4

    ::

        sage: g = graphs.LollipopGraph(13,4)
        sage: g.show() # long time

    Create several lollipop graphs in a Sage graphics array

    ::

        sage: g = []
        sage: j = []
        sage: for i in range(6):
        ....:     k = graphs.LollipopGraph(i+3,4)
        ....:     g.append(k)
        sage: for i in range(2):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time
    """
    pos_dict = {}

    for i in range(n1):
        x = float(cos((pi/4) - ((2*pi)/n1)*i) - n2/2 - 1)
        y = float(sin((pi/4) - ((2*pi)/n1)*i) - n2/2 - 1)
        j = n1-1-i
        pos_dict[j] = (x,y)
    for i in range(n1, n1+n2):
        x = float(i - n1 - n2/2 + 1)
        y = float(i - n1 - n2/2 + 1)
        pos_dict[i] = (x,y)
    G = graph.Graph(dict( (i,range(i+1,n1)) for i in range(n1) ), pos=pos_dict, name="Lollipop Graph")
    G.add_vertices( range(n1+n2) )
    G.add_path( range(n1-1,n1+n2) )
    return G

def PathGraph(n, pos=None):
    """
    Returns a path graph with n nodes. Pos argument takes a string
    which is either 'circle' or 'line', (otherwise the default is
    used). See the plotting section below for more detail.

    A path graph is a graph where all inner nodes are connected to
    their two neighbors and the two end-nodes are connected to their
    one inner neighbors. (i.e.: a cycle graph without the first and
    last node connected).

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the graph may
    be drawn in one of two ways: The 'line' argument will draw the
    graph in a horizontal line (left to right) if there are less than
    11 nodes. Otherwise the 'line' argument will append horizontal
    lines of length 10 nodes below, alternating left to right and right
    to left. The 'circle' argument will cause the graph to be drawn in
    a cycle-shape, with the first node at the top and then about the
    circle in a clockwise manner. By default (without an appropriate
    string argument) the graph will be drawn as a 'circle' if 10 n 41
    and as a 'line' for all other n.

    EXAMPLES: Show default drawing by size: 'line': n 11

    ::

        sage: p = graphs.PathGraph(10)
        sage: p.show() # long time

    'circle': 10 n 41

    ::

        sage: q = graphs.PathGraph(25)
        sage: q.show() # long time

    'line': n 40

    ::

        sage: r = graphs.PathGraph(55)
        sage: r.show() # long time

    Override the default drawing::

        sage: s = graphs.PathGraph(5,'circle')
        sage: s.show() # long time
    """
    G = graph.Graph(n, name="Path graph")

    pos_dict = {}

    # Choose appropriate drawing pattern
    circle = False
    if pos == "circle": circle = True
    elif pos == "line": circle = False
    # Otherwise use default by size of n
    elif 10 < n < 41: circle = True

    # Draw 'circle'
    if circle:
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos_dict[i] = (x,y)
    # Draw 'line'
    else:
        counter = 0 # node index
        rem = n%10 # remainder to appear on last row
        rows = n//10 # number of rows (not counting last row)
        lr = True # left to right

        for i in range(rows): # note that rows doesn't include last row
            y = -i
            for j in range(10):
                if lr:
                    x = j
                else:
                    x = 9 - j
                pos_dict[counter] = (x,y)
                counter += 1
            if lr: lr = False
            else: lr = True
        y = -rows
        for j in range(rem): # last row
            if lr:
                x = j
            else:
                x = 9 - j
            pos_dict[counter] = (x,y)
            counter += 1

    G.add_edges( (i,i+1) for i in range(n-1) )
    G.set_pos(pos_dict)
    return G

def StarGraph(n):
    """
    Returns a star graph with n+1 nodes.

    A Star graph is a basic structure where one node is connected to
    all other nodes.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, each star
    graph will be displayed with the first (0) node in the center, the
    second node (1) at the top, with the rest following in a
    counterclockwise manner. (0) is the node connected to all other
    nodes.

    The star graph is a good opportunity to compare efficiency of
    filling a position dictionary vs. using the spring-layout algorithm
    for plotting. As far as display, the spring-layout should push all
    other nodes away from the (0) node, and thus look very similar to
    this constructor's positioning.

    EXAMPLES::

        sage: import networkx

    Compare the plots::

        sage: n = networkx.star_graph(23)
        sage: spring23 = Graph(n)
        sage: posdict23 = graphs.StarGraph(23)
        sage: spring23.show() # long time
        sage: posdict23.show() # long time

    View many star graphs as a Sage Graphics Array

    With this constructor (i.e., the position dictionary filled)

    ::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.StarGraph(i+3)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time

    Compared to plotting with the spring-layout algorithm

    ::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     spr = networkx.star_graph(i+3)
        ....:     k = Graph(spr)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time
    """
    pos_dict = {}
    pos_dict[0] = (0,0)
    for i in range(1,n+1):
        x = float(cos((pi/2) + ((2*pi)/n)*(i-1)))
        y = float(sin((pi/2) + ((2*pi)/n)*(i-1)))
        pos_dict[i] = (x,y)
    return graph.Graph({0:range(1,n+1)}, pos=pos_dict, name="Star graph")

