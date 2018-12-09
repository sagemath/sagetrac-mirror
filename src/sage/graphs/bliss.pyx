r"""
Interface with bliss: graph (iso/auto)morphism

Implemented functions:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`automorphism_group` | Returns the automorphism group of the given (di)graph
    :meth:`canonical_form` | Computes a canonical certificate for the given (di) graph.

AUTHORS:

    - Jernej Azarija
"""

#*****************************************************************************
#       Copyright (C) 2015 Jernej Azarija
#       Copyright (C) 2015 Nathann Cohen <nathann.cohen@gail.com>
#       Copyright (C) 2018 Christian Stump <christian.stump@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import numpy
from operator import itemgetter

from cpython cimport PyObject
from libc.limits cimport LONG_MAX

cdef extern from "bliss/graph.hh" namespace "bliss":

    cdef cppclass Stats:
        pass

    cdef cppclass AbstractGraph:
        pass

    cdef cppclass Graph(AbstractGraph):
        Graph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void find_automorphisms(Stats&, void (*)(void* , unsigned int,
                    const unsigned int*), void*)
        void change_color(const unsigned int, const unsigned int);
        const unsigned int* canonical_form(Stats&, void (*)(void*,unsigned int,
                    const unsigned int*), void*)

    cdef cppclass Digraph(AbstractGraph):
        Digraph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void find_automorphisms(Stats&, void (*)(void* , unsigned int,
                    const unsigned int*), void*)
        void change_color(const unsigned int, const unsigned int);
        const unsigned int* canonical_form(Stats&, void (*)(void*,unsigned int,
                    const unsigned int*), void*)
        unsigned int get_hash()

cdef void add_gen(void *user_param, unsigned int n, const unsigned int *aut):
    r"""
    This function is called each time a new generator of the automorphism group
    is found.

    This function is used to append the new generators to a Python list. Its
    main job is to translate a permutation into dijoint cycles.

    INPUT:

    - ``user_param`` (``void *``) -- in the current implementation, it points
      toward a Python object which is a pair
      ``(list_of_current_generators,vert_to_integer_labelling)``.

    - ``n`` (int) -- number of points in the graph.

    - ``aut`` (int *) -- an automorphism of the graph.
    """
    cdef int tmp     = 0
    cdef int marker  = 0
    cdef int cur     = 0
    perm        = []
    done        = [False]*n

    gens, int_to_vertex = <object> <PyObject *> user_param

    while True:
        while cur < n and done[cur]:
            cur+=1
        if cur == n:
            break

        marker = tmp = cur
        cycle = [int_to_vertex[cur]]
        done[cur] = True

        while aut[tmp] != marker:
            tmp = aut[tmp]
            done[tmp] = True
            cycle.append(int_to_vertex[tmp])

        perm.append(tuple(cycle))
    gens.append(perm)

cdef void empty_hook(void *user_param , unsigned int n, const unsigned int *aut):
    return

#####################################################
# constructing bliss graphs from edge lists
#####################################################

cdef Graph *bliss_graph_from_labelled_edges(int Vnr, int Lnr, Vout, Vin, labels, partition):
    r"""
    Return a bliss graph from the input data

    For edge labelled graphs, the bliss graph is constructed using Vnr * log(Lnr) many vertices as described in Sect 14 in the `nauty reference manual <http://pallini.di.uniroma1.it/Guide.html>`_.

    .. WARNING::

        the input is not checked for correctness, any wrong input will result in a segfault

    INPUT:

    - ``Vnr`` (number of vertices, such that the vertices are 0 ... Vnr-1)

    - ``Lnr`` (number of labels, such that the labels are 0 ... Lnr-1)

    - ``Vout`` (the list of vertices of outgoing edges)

    - ``Vin`` (the list of vertices of ingoing edges)

    - ``labels`` (the list of edge labels)

    - ``partition`` -- a partition of the vertex set
    """
    cdef Py_ssize_t i, j
    cdef int logLnr
    cdef str binrep
    cdef str ind

    cdef Graph *g
    cdef int x,y, lab

    if Lnr == 1:
        g = new Graph(Vnr)
        if g == NULL:
            raise MemoryError("Allocation Failed")
    else:
        logLnr = len(numpy.binary_repr(Lnr))
        g = new Graph(Vnr*logLnr)
        if g == NULL:
            raise MemoryError("Allocation Failed")
        for i from 0 <= i < Vnr:
            for j from 1 <= j < logLnr:
                g.add_edge((j-1)*Vnr+i,j*Vnr+i)

    cdef int Enr = len(Vout)

    for i from 0 <= i < Enr:
        x   = Vout[i]
        y   = Vin[i]
        if Lnr == 1:
            lab = 0
        else:
            lab = labels[i]

        if lab != 0:
            lab = lab+1
            binrep = numpy.binary_repr(lab, logLnr)

            for j from 0 <= j < logLnr:
                ind = binrep[j]
                if ind == "1":
                    g.add_edge((logLnr-1-j)*Vnr+x,(logLnr-1-j)*Vnr+y)
        else:
            g.add_edge(x,y)

    if not bool(partition):
        partition = [list(range(Vnr))]
    cdef Pnr = len(partition)
    for i from 0 <= i < Pnr:
        for v in partition[i]:
            if Lnr == 1:
                g.change_color(v, i)
            else:
                for j from 0 <= j < logLnr:
                    g.change_color(j*Vnr+v, j*Pnr+i)
    return g

cdef Digraph *bliss_digraph_from_labelled_edges(int Vnr, int Lnr, Vout, Vin, labels, partition):
    r"""
    Return a bliss digraph from the input data

    For edge labelled graphs, the bliss graph is constructed using Vnr * log(Lnr) many vertices as described in Sect 14 in the `nauty reference manual
    <http://pallini.di.uniroma1.it/Guide.html>`_.

    .. WARNING::

        the input is not checked for correctness, any wrong input will result in a segfault

    INPUT:

    - ``Vnr`` (number of vertices, such that the vertices are 0 ... Vnr-1)

    - ``Lnr`` (number of labels, such that the labels are 0 ... Lnr-1)

    - ``Vout`` (the list of vertices of outgoing edges)

    - ``Vin`` (the list of vertices of ingoing edges)

    - ``labels`` (the list of edge labels)

    - ``partition`` -- a partition of the vertex set
    """
    cdef Py_ssize_t i, j
    cdef int logLnr
    cdef str binrep
    cdef str ind

    cdef Digraph *g
    cdef int x,y, lab

    if Lnr == 1:
        g = new Digraph(Vnr)
        if g == NULL:
            raise MemoryError("Allocation Failed")
    else:
        logLnr = len(numpy.binary_repr(Lnr))
        g = new Digraph(Vnr*logLnr)
        if g == NULL:
            raise MemoryError("Allocation Failed")
        for i from 0 <= i < Vnr:
            for j from 1 <= j < logLnr:
                g.add_edge((j-1)*Vnr+i,j*Vnr+i)

    cdef int Enr = len(Vout)

    for i from 0 <= i < Enr:
        x   = Vout[i]
        y   = Vin[i]
        if Lnr == 1:
            lab = 0
        else:
            lab = labels[i]

        if lab != 0:
            lab = lab+1
            binrep = numpy.binary_repr(lab)

            for j from 0 <= j < logLnr:
                ind = binrep[j]
                if ind == "1":
                    g.add_edge((logLnr-1-j)*Vnr+x,(logLnr-1-j)*Vnr+y)
        else:
            g.add_edge(x,y)

    if not bool(partition):
        partition = [list(range(Vnr))]
    cdef Pnr = len(partition)
    for i from 0 <= i < Pnr:
        for v in partition[i]:
            if Lnr == 1:
                g.change_color(v, i)
            else:
                for j from 0 <= j < logLnr:
                    g.change_color(j*Vnr+v, j*Pnr+i)
    return g

#####################################################
# canonical form from graph or edge list
#####################################################

cdef canonical_form_from_edge_list(int Vnr, list Vout, list Vin, int Lnr=1, list labels=[], list partition=None, bint directed=False, bint certificate=False):
    r"""
    Return an unsorted list of labelled edges of a canonical form.

    INPUT:

    - ``Vnr`` -- number of vertices such that the vertices are 0 ... Vnr-1

    - ``Vout`` -- the list of vertices of outgoing edges

    - ``Vin`` -- the list of vertices of ingoing edges

    - ``Lnr`` -- number of labels such that the labels are 0 ... Lnr-1

    - ``labels`` -- the list of edge labels)

    - ``partition`` -- a partition of the vertex set

    - ``directed`` -- boolean flag whether the edges are directed or not

    - ``certificate`` -- boolean flag whether to return the isomorphism to obtain the canonical labelling
    """
    # We need this to convert the numbers from <unsigned int> to
    # <long>. This assertion should be true simply for memory reasons.
    assert <unsigned long>(Vnr) <= <unsigned long>LONG_MAX

    cdef const unsigned int* aut
    cdef Graph* g
    cdef Digraph* d
    cdef Stats s
    cdef dict relabel

    cdef list new_edges = []
    cdef long e, f

    if directed:
        d = bliss_digraph_from_labelled_edges(Vnr, Lnr, Vout, Vin, labels, partition)
        aut = d.canonical_form(s, empty_hook, NULL)
    else:
        g = bliss_graph_from_labelled_edges(Vnr, Lnr, Vout, Vin, labels, partition)
        aut = g.canonical_form(s, empty_hook, NULL)

    for i from 0 <= i < len(Vout):
        x = Vout[i]
        y = Vin[i]
        e = aut[x]
        f = aut[y]
        if Lnr == 1:
            if not bool(labels):
                lab = None
            else:
                lab = labels[0]
            if directed:
                new_edges.append( (e,f,lab) )
            else:
                new_edges.append( (e,f,lab) if e > f else (f,e,lab))
        else:
            lab = labels[i]
            if directed:
                new_edges.append( (e,f,lab) )
            else:
                new_edges.append( (e,f,lab) if e > f else (f,e,lab))

    if certificate:
        relabel = {v: <long>aut[v] for v in range(Vnr)}

    if directed:
        del d
    else:
        del g

    if certificate:
        return new_edges, relabel
    else:
        return new_edges

cpdef canonical_form(G, partition=None, return_graph=False, use_edge_labels=True, certificate=False):
    r"""
    Return the canonical label of ``G``.

    A canonical label ``canonical_form(G)`` of ``G`` is a (di)graph defined on
    `\{0,...,n-1\}` such that ``G`` is isomorphic to ``H`` if and only if
    ``canonical_form(G)`` is equal to ``canonical_form(H)``.

    INPUT:

    - ``G`` -- A graph or digraph.

    - ``partition`` -- A partition of the vertices of ``G`` into color classes.
      Defaults to ``None``.

    - ``return_graph`` -- If set to ``True``, ``canonical_form`` returns the
      canonical graph of ``G``. Otherwise, it returns its set of edges.

    - ``edge_labels`` -- A boolean whether or not to consider edge labels.

    - ``certificate`` -- If set to ``True`` returns the labeling of G into a
      canonical graph.

    TESTS::

        sage: from sage.graphs.bliss import canonical_form                  # optional - bliss
        sage: G = graphs.PetersenGraph()                                    # optional - bliss
        sage: canonical_form(G)                                             # optional - bliss
        [(2, 0, None),
         (2, 1, None),
         (3, 0, None),
         (4, 1, None),
         (5, 3, None),
         (5, 4, None),
         (6, 0, None),
         (6, 4, None),
         (7, 1, None),
         (7, 3, None),
         (8, 2, None),
         (8, 5, None),
         (9, 6, None),
         (9, 7, None),
         (9, 8, None)]

        sage: P = graphs.GeneralizedPetersenGraph(5,2)                      # optional - bliss
        sage: Q = graphs.PetersenGraph()                                    # optional - bliss
        sage: canonical_form(P) == canonical_form(Q)                        # optional - bliss
        True

        sage: canonical_form(Graph(15),return_graph=True)                   # optional - bliss
        Graph on 15 vertices
        sage: g = digraphs.RandomTournament(40)                             # optional - bliss
        sage: g.is_isomorphic(canonical_form(g,return_graph=True))          # optional - bliss
        True

        sage: g1 = graphs.RandomGNP(100,.4)                                 # optional - bliss
        sage: r = Permutations(range(100)).random_element()                 # optional - bliss
        sage: g2 = Graph([(r[u],r[v]) for u,v in g1.edges(labels=False)])   # optional - bliss
        sage: g1 = canonical_form(g1,return_graph=True)                     # optional - bliss
        sage: g2 = canonical_form(g2,return_graph=True)                     # optional - bliss
        sage: g2 == g2                                                      # optional - bliss
        True

        sage: g = Graph({1: [2]})
        sage: g_ = canonical_form(g, return_graph=True, certificate=True)   # optional - bliss
        sage: 0 in g_[0]                                                    # optional - bliss
        True
    """
    # We need this to convert the numbers from <unsigned int> to
    # <long>. This assertion should be true simply for memory reasons.
    cdef unsigned long Vnr = G.order()
    assert Vnr <= <unsigned long>LONG_MAX

    cdef bint directed = G.is_directed()

    cdef int labInd
    cdef list Vout   = []
    cdef list Vin    = []
    cdef list labels = []

    vert2int         = {}
    int2vert         = [None]*Vnr
    edge_labels      = []
    edge_labels_rev  = {}
    cdef int Lnr     = 0

    for i,v in enumerate(G.vertices()):
        vert2int[v] = i
        int2vert[i] = v

    if bool(partition):
        partition = [ [ vert2int[i] for i in part ] for part in partition ]

    for x,y,lab in G.edges(labels=True):
        if use_edge_labels is False:
            lab = None
        try:
            labInd = edge_labels_rev[lab]
        except KeyError:
            labInd = Lnr
            Lnr    = Lnr+1
            edge_labels_rev[lab] = labInd
            edge_labels.append(lab)

        Vout.append(vert2int[x])
        Vin.append(vert2int[y])
        labels.append(labInd)

    lab_relabels = [ lab for _,lab in sorted(edge_labels_rev.iteritems(), key=itemgetter(0)) ]
    labels = [lab_relabels[i] for i in labels]
    new_edges, relabel = canonical_form_from_edge_list(Vnr, Vout, Vin, Lnr, labels, partition, directed, certificate=True)

    new_edges = [ (x,y,edge_labels[lab]) for x,y,lab in new_edges ]
    relabel = { int2vert[i]:j for i,j in relabel.iteritems() }

    if return_graph:
        if directed:
            from sage.graphs.graph import DiGraph
            G = DiGraph(new_edges,loops=G.allows_loops(),multiedges=G.allows_multiple_edges())
        else:
            from sage.graphs.graph import Graph
            G = Graph(new_edges,loops=G.allows_loops(),multiedges=G.allows_multiple_edges())

        G.add_vertices(vert2int.values())
        return (G, relabel) if certificate else G

    return (sorted(new_edges), relabel) if certificate else sorted(new_edges)

#####################################################
# automorphism group from graphs
#####################################################

cdef automorphism_group_gens_from_edge_list(int Vnr, Vout, Vin, int Lnr=1, labels=[], int2vert=[], partition=None, bint directed=False):
    r"""
    Return an unsorted list of labelled edges of a canonical form.

    INPUT:

    - ``Vnr`` -- number of vertices such that the vertices are 0 ... Vnr-1

    - ``Vout`` -- the list of vertices of outgoing edges

    - ``Vin`` -- the list of vertices of ingoing edges

    - ``Lnr`` -- number of labels such that the labels are 0 ... Lnr-1

    - ``labels`` -- the list of edge labels)

    - ``partition`` -- a partition of the vertex set

    - ``directed`` -- boolean flag whether the edges are directed or not
    """
    # We need this to convert the numbers from <unsigned int> to
    # <long>. This assertion should be true simply for memory reasons.
    assert <unsigned long>(Vnr) <= <unsigned long>LONG_MAX

    cdef Graph* g
    cdef Digraph* d
    cdef Stats s

    if int2vert == []:
        int2vert = list(range(Vnr))

    # the following is needed because the the internal graph has
    # size Vnr*logLnr for labelled graphs
    if Lnr != 1:
        logLnr = len(numpy.binary_repr(Lnr))
        int2vert.extend([None]*(Vnr*(logLnr-1)))

    gens = []
    data = (gens, int2vert)

    if directed:
        d = bliss_digraph_from_labelled_edges(Vnr, Lnr, Vout, Vin, labels, partition)
        d.find_automorphisms(s, add_gen, <PyObject *> data)
        del d
    else:
        g = bliss_graph_from_labelled_edges(Vnr, Lnr, Vout, Vin, labels, partition)
        g.find_automorphisms(s, add_gen, <PyObject *> data)
        del g

    return [ [ cyc for cyc in gen if cyc[0] is not None ] for gen in gens ]

cpdef automorphism_group(G, partition=None, use_edge_labels=True):
    """
    Computes the automorphism group of ``G`` subject to the vertex coloring ``partition``, if given.

    The graph ``G`` can be a directed or undirected graph with or without edge labellings.

    Observe the neither the vertex colorings nor the edge colorings are interchangeable.

    INPUT:

    - ``G`` -- A graph

    - ``partition`` -- A partition of the vertices of ``G`` into color classes.
      Defaults to ``None``, which is equivalent to a partition of size 1.

    - ``edge_labels`` -- A boolean whether or not to consider edge labels.

    EXAMPLES::

        sage: from sage.graphs.bliss import automorphism_group                  # optional - bliss

    Computing the automorphism group of a graph or digraph::

        sage: G = graphs.CompleteMultipartiteGraph([1,1,1,2])                   # optional - bliss
        sage: automorphism_group(G).cardinality()                               # optional - bliss
        12
        sage: D = DiGraph(G.edges())                                            # optional - bliss
        sage: automorphism_group(D).cardinality()                               # optional - bliss
        2

    Observe that the order 12 is given by permuting the first three vertices, or the last two
    in the case of a graph, while only the latter two are possible in the case of a directed
    graph.

    Partitioning the vertices into classes::

        sage: G = graphs.CompleteMultipartiteGraph([3,2])                       # optional - bliss
        sage: automorphism_group(G).cardinality()                               # optional - bliss
        12
        sage: automorphism_group(G,partition=[[0],[1],[2],[3,4]]).cardinality() # optional - bliss
        2
        sage: automorphism_group(G,partition=[[0],[1,2],[3,4]]).cardinality()   # optional - bliss
        4

        sage: automorphism_group(G,partition=[[1,2],[0,3],[4]]).cardinality()   # optional - bliss
        2

    Partitioning the edges into classes::

        sage: G = Graph(graphs.CompleteMultipartiteGraph([8,2]), sparse=True)   # optional - bliss
        sage: for i,j,_ in G.edges():                                           # optional - bliss
        ....:     if 0 <= i < 3:                                                # optional - bliss
        ....:         G.set_edge_label(i,j,"A")                                 # optional - bliss
        ....:     if 3 <= i < 6:                                                # optional - bliss
        ....:         G.set_edge_label(i,j,"B")                                 # optional - bliss
        ....:     if 6 <= i < 8:                                                # optional - bliss
        ....:         G.set_edge_label(i,j,"C")                                 # optional - bliss

        sage: factor(automorphism_group(G).cardinality())                       # optional - bliss
        2^4 * 3^2
        sage: automorphism_group(G,[[0],[1],[2,3],[4,5],[6,7],[8],[9]]).cardinality()   # optional - bliss
        4

    TESTS::

        sage: from sage.graphs.bliss import automorphism_group                  # optional - bliss
        sage: G = graphs.PetersenGraph()                                        # optional - bliss
        sage: automorphism_group(G).is_isomorphic(G.automorphism_group())       # optional - bliss
        True

        sage: G = graphs.HeawoodGraph()                                         # optional - bliss
        sage: p = G.bipartite_sets()                                            # optional - bliss
        sage: A = G.automorphism_group(partition=[list(p[0]), list(p[1])])      # optional - bliss
        sage: automorphism_group(G, partition=p).is_isomorphic(A)               # optional - bliss
        True

        sage: G = graphs.CompleteMultipartiteGraph([5,7,11])
        sage: B = automorphism_group(G)                                         # optional - bliss
        sage: B.cardinality() == prod(factorial(n) for n in [5,7,11])           # optional - bliss
        True

        sage: G = Graph(graphs.CompleteMultipartiteGraph([8,8,8,5]),sparse=True)# optional - bliss
        sage: for i,j,_ in G.edges():                                           # optional - bliss
        ....:     if 0 <= i < 3:                                                # optional - bliss
        ....:         G.set_edge_label(i,j,"A")                                 # optional - bliss
        ....:     if 3 <= i < 6:                                                # optional - bliss
        ....:         G.set_edge_label(i,j,"B")                                 # optional - bliss
        ....:     if 6 <= i < 8:                                                # optional - bliss
        ....:         G.set_edge_label(i,j,"C")                                 # optional - bliss
        ....:
        sage: automorphism_group(G).cardinality() == prod( factorial(n) for n in [3,3,2,8,8,5,2] )  # optional - bliss
        True
        sage: automorphism_group(G, use_edge_labels=False).cardinality() == prod( factorial(n) for n in [8,8,8,5,3] )  # optional - bliss
        True
        sage: automorphism_group(G,[[0 .. 7],[8 .. 11],[12 .. 28]]).cardinality() == prod( factorial(n) for n in [3,3,2,4,4,8,5] )  # optional - bliss
        True

        sage: G = Graph()                                                       # optional - bliss
        sage: G.add_edges((i,j,"A") for i in range(0, 2) for j in range(14,20)) # optional - bliss
        sage: G.add_edges((i,j,"B") for i in range(2, 5) for j in range(14,20)) # optional - bliss
        sage: G.add_edges((i,j,"C") for i in range(5, 9) for j in range(14,20)) # optional - bliss
        sage: G.add_edges((i,j,"D") for i in range(9,14) for j in range(14,20)) # optional - bliss
        sage: A = automorphism_group(G)                                         # optional - bliss
        sage: print(A.gens())                                                   # random, optional - bliss
        [(9,13), (18,19), (17,18), (16,17), (15,16), (14,15), (12,9), (11,12), (10,11), (7,8), (6,7), (5,6), (3,4), (2,3), (0,1)]
        sage: A.cardinality() == prod(factorial(n) for n in [2,3,4,5,6])        # optional - bliss
        True

        sage: alpha = "abcdefghijklmnopqrstuvwxyz"

        sage: G = Graph()                                                       # optional - bliss
        sage: G.add_edges((alpha[i],alpha[j],"A") for i in range(0, 2) for j in range(14,20))   # optional - bliss
        sage: G.add_edges((alpha[i],alpha[j],"B") for i in range(2, 5) for j in range(14,20))   # optional - bliss
        sage: G.add_edges((alpha[i],alpha[j],"C") for i in range(5, 9) for j in range(14,20))   # optional - bliss
        sage: G.add_edges((alpha[i],alpha[j],"D") for i in range(9,14) for j in range(14,20))   # optional - bliss
        sage: A = automorphism_group(G)                                         # optional - bliss
        sage: print(A.gens())                                                   # random, optional - bliss
        [('r','t'), ('s','r'), ('p','s'), ('q','p'), ('o','q'), ('l','n'), ('m','l'), ('j','m'), ('k','j'), ('i','h'), ('f','i'), ('g','f'), ('e','d'), ('c','e'), ('a','b')]
        sage: A.cardinality() == prod(factorial(n) for n in [2,3,4,5,6])        # optional - bliss
        True

        sage: gg = graphs.CompleteGraph(5)                                      # optional - bliss
        sage: gg.allow_loops(True)                                              # optional - bliss
        sage: gg.add_edge(0,0)                                                  # optional - bliss
        sage: gg.add_edge(1,1)                                                  # optional - bliss
        sage: automorphism_group(gg).cardinality()                              # optional - bliss
        12
        sage: automorphism_group(gg,[[0],[1,2,3,4]]).cardinality()              # optional - bliss
        6

    Making sure that #25426 is fixed:

        sage: j = matrix([(3, 2, 1, 0, 0),
        ....:  (2, 2, 0, 1, 0),
        ....:  (1, 0, 3, 0, 2),
        ....:  (0, 1, 0, 2, 1),
        ....:  (0, 0, 2, 1, 2)])
        sage: j.automorphisms_of_rows_and_columns()
        [((), ()), ((1,3)(2,5), (1,3)(2,5))]
    """
    # We need this to convert the numbers from <unsigned int> to
    # <long>. This assertion should be true simply for memory reasons.
    cdef unsigned long Vnr = G.order()
    assert Vnr <= <unsigned long>LONG_MAX

    cdef bint directed = G.is_directed()

    cdef int labInd
    cdef list Vout   = []
    cdef list Vin    = []
    cdef list labels = []

    vert2int         = {}
    int2vert         = [None]*Vnr
    edge_labels      = []
    edge_labels_rev  = {}
    cdef int Lnr     = 0

    for i, v in enumerate(G.vertex_iterator()):
        vert2int[v] = i
        int2vert[i] = v

    if bool(partition):
        partition = [ [ vert2int[i] for i in part ] for part in partition ]

    for x,y,lab in G.edge_iterator(labels=True):
        if use_edge_labels is False:
            lab = None
        try:
            labInd = edge_labels_rev[lab]
        except KeyError:
            labInd = Lnr
            Lnr    = Lnr+1
            edge_labels_rev[lab] = labInd
            edge_labels.append(lab)

        Vout.append(vert2int[x])
        Vin.append(vert2int[y])
        labels.append(labInd)

    lab_relabels = [ lab for _,lab in sorted(edge_labels_rev.iteritems(), key=itemgetter(0)) ]
    labels = [lab_relabels[i] for i in labels]

    gens = automorphism_group_gens_from_edge_list(Vnr, Vout, Vin, Lnr, labels, int2vert, partition, directed)

    from sage.groups.perm_gps.permgroup import PermutationGroup
    return PermutationGroup(gens,domain=sorted(G))

#####################################################
# old direct interactions graphs <-> bliss graphs
#####################################################

cdef Graph *bliss_graph(G, partition, vert2int, int2vert):
    r"""
    Return a bliss copy of a graph G

    INPUT:

    - ``G`` (a graph)

    - ``partition`` -- a partition of the vertex set.

    - ``vert2int, int2vert`` -- Two empty dictionaries. The entries of the
      dictionary are later set to record the labeling of our graph. They are
      taken as arguments to avoid technicalities of returning Python objects in
      Cython functions.
    """
    cdef Graph *g = new Graph(G.order())

    if g == NULL:
        raise MemoryError("Allocation Failed")

    for i,v in enumerate(G.vertices()):
        vert2int[v] = i
        int2vert[i] = v

    for x,y in G.edges(labels=False):
       g.add_edge(vert2int[x],vert2int[y])

    if partition:
        for i in xrange(1, len(partition)):
            for v in partition[i]:
                g.change_color(vert2int[v], i)
    return g

cdef Digraph *bliss_digraph(G, partition, vert2int, int2vert):
    r"""
    Return a bliss copy of a digraph G

    INPUT:

    - ``G`` (a digraph)

    - ``partition`` -- a partition of the vertex set.

    - ``vert2int, int2vert`` -- Two empty dictionaries. The entries of the
      dictionary are later set to record the labeling of our graph. They are
      taken as arguments to avoid technicalities of returning Python objects in
      Cython functions.
    """
    cdef Digraph *g = new Digraph(G.order())

    for i,v in enumerate(G.vertices()):
        vert2int[v] = i
        int2vert[i] = v

    if g == NULL:
        raise MemoryError("Allocation Failed")

    for x,y in G.edges(labels=False):
        g.add_edge(vert2int[x],vert2int[y])

    if partition:
        for i in xrange(1, len(partition)):
            for v in partition[i]:
                g.change_color(vert2int[v], i)
    return g
