"""
Wrapper for Boyer's (C) planarity algorithm.
"""

cdef extern from "planarity_c/graph.h":
    ctypedef struct graphNode:
        int v
        int link[2]
    ctypedef graphNode * graphNodeP

    ctypedef struct BM_graph:
        graphNodeP G
        int N
    ctypedef BM_graph * graphP

    cdef int OK, EMBEDFLAGS_PLANAR, NONEMBEDDABLE, NOTOK

    cdef graphP gp_New()
    cdef void gp_Free(graphP *pGraph)
    cdef int gp_InitGraph(graphP theGraph, int N)
    cdef int gp_AddEdge(graphP theGraph, int u, int ulink, int v, int vlink)
    cdef int gp_Embed(graphP theGraph, int embedFlags)
    cdef int gp_SortVertices(graphP theGraph)

def is_planar(g, kuratowski=False, set_pos=False, set_embedding=False, circular=False):
    """
    Calls Boyer's planarity algorithm to determine whether g is
    planar.  If kuratowski is False, returns True if g is planar,
    False otherwise.  If kuratowski is True, returns a tuple, first
    entry is a boolean (whether or not the graph is planar) and second
    entry is a Kuratowski subgraph, i.e. an edge subdivision of
    `K_5` or `K_{3,3}` (if not planar) or ``None`` (if planar).  Also, will set
    an ``_embedding`` attribute for the graph ``g`` if ``set_embedding`` is set
    to True.

    INPUT:
        kuratowski -- If True, return a tuple of a boolean and either None
        or a Kuratowski subgraph (i.e. an edge subdivision of `K_5` or `K_{3,3}`)
        set_pos -- if True, uses Schnyder's algorithm to determine positions
        set_embedding -- if True, records the combinatorial embedding returned
        (see g.get_embedding())
        circular -- if True, test for circular planarity

    EXAMPLES::

        sage: G = graphs.DodecahedralGraph()
        sage: from sage.graphs.planarity import is_planar
        sage: is_planar(G)
        True
        sage: Graph('@').is_planar()
        True

    TESTS:

    We try checking the planarity of all graphs on 7 or fewer
    vertices.  In fact, to try to track down a segfault, we do it
    twice. ::

        sage: import networkx.generators.atlas  # long time
        sage: atlas_graphs = [Graph(i) for i in networkx.generators.atlas.graph_atlas_g()] # long time
        sage: a = [i for i in [1..1252] if atlas_graphs[i].is_planar()] # long time
        sage: b = [i for i in [1..1252] if atlas_graphs[i].is_planar()] # long time
        sage: a == b # long time
        True

    There were some problems with ``set_pos`` stability in the past,
    so let's check if this this runs without exception::

        sage: for i,g in enumerate(atlas_graphs):                         # long time
        ....:     if (not g.is_connected() or i==0):
        ....:         continue
        ....:     _ = g.is_planar(set_embedding=True, set_pos=True)
    """
    if set_pos and not g.is_connected():
        raise ValueError("is_planar() cannot set vertex positions for a disconnected graph")

    # First take care of a trivial cases
    if g.size() == 0: # There are no edges
        if set_embedding:
            g._embedding = dict((v, []) for v in g.vertices())
        return (True, None) if kuratowski else True
    if len(g) == 2 and g.is_connected(): # P_2 is too small to be triangulated
        u,v = g.vertices()
        if set_embedding:
            g._embedding = { u: [v], v: [u] }
        if set_pos:
            g._pos = { u: [0,0], v: [0,1] }
        return (True, None) if kuratowski else True

    # create to and from mappings to relabel vertices to the set {0,...,n-1}
    cdef int i
    listto = g.vertices()
    ffrom = {}
    for vvv in listto:
        ffrom[vvv] = listto.index(vvv)
    to = {}
    for i from 0 <= i < len(listto):
        to[i] = listto[i]
    g.relabel(ffrom)

    cdef graphP theGraph
    theGraph = gp_New()
    cdef int status
    status = gp_InitGraph(theGraph, g.order())
    if status != OK:
        raise RuntimeError("gp_InitGraph status is not ok.")
    for u, v, _ in g.edge_iterator():
        status = gp_AddEdge(theGraph, u, 0, v, 0)
        if status == NOTOK:
            raise RuntimeError("gp_AddEdge status is not ok.")
        elif status == NONEMBEDDABLE:
            # We now know that the graph is nonplanar.
            if not kuratowski:
                return False
            # With just the current edges, we have a nonplanar graph,
            # so to isolate a kuratowski subgraph, just keep going.
            break

    status = gp_Embed(theGraph, EMBEDFLAGS_PLANAR)
    gp_SortVertices(theGraph)

    # use to and from mappings to relabel vertices back from the set {0,...,n-1}
    g.relabel(to)

    if status == NOTOK:
        raise RuntimeError("Status is not ok.")
    elif status == NONEMBEDDABLE:
        # Kuratowski subgraph isolator
        g_dict = {}
        from sage.graphs.graph import Graph
        for i from 0 <= i < theGraph.N:
            linked_list = []
            j = theGraph.G[i].link[1]
            while j >= theGraph.N:
                linked_list.append(to[theGraph.G[j].v])
                j = theGraph.G[j].link[1]
            if len(linked_list) > 0:
                g_dict[to[i]] = linked_list
        G = Graph(g_dict)
        gp_Free(&theGraph)
        if kuratowski:
            return (False, G)
        else:
            return False
    else:
        if not circular:
            if set_embedding:
                emb_dict = {}
                #for i in range(theGraph.N):
                for i from 0 <= i < theGraph.N:
                    linked_list = []
                    j = theGraph.G[i].link[1]
                    while j >= theGraph.N:
                        linked_list.append(to[theGraph.G[j].v])
                        j = theGraph.G[j].link[1]
                    emb_dict[to[i]] = linked_list
                g._embedding = emb_dict
            if set_pos:
                g.set_planar_positions()
        else:
            if set_embedding:
                # Take counter-clockwise embedding if circular planar test
                # Also, pos must be set after removing extra vertex and edges

                # This is separated out here for now because in the circular case,
                # setting positions would have to come into play while the extra
                # "wheel" or "star" is still part of the graph.

                emb_dict = {}
                #for i in range(theGraph.N):
                for i from 0 <= i < theGraph.N:
                    linked_list = []
                    j = theGraph.G[i].link[0]
                    while j >= theGraph.N:
                        linked_list.append(to[theGraph.G[j].v])
                        j = theGraph.G[j].link[0]
                    emb_dict[to[i]] = linked_list
                g._embedding = emb_dict
        gp_Free(&theGraph)
        if kuratowski:
            return (True,None)
        else:
            return True
