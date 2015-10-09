r"""
Functions for reading/building graphs/digraphs.

This module gathers functions needed to build a graph from any other data.

.. NOTE::

    This is a **internal** module of Sage. All features implemented here are
    made available to end-users through the constructors of :class:`Graph` and
    :class:`DiGraph`.

Note that because they are called by the constructors of :class:`Graph` and
:class:`DiGraph`, most of these functions modify a graph inplace.

{INDEX_OF_FUNCTIONS}

Functions
---------

"""

def from_graph6(G, g6_string):
    r"""
    Fill ``G`` with the data of a graph6 string.

    INPUT:

    - ``G`` -- a graph

    - ``g6_string`` -- a graph6 string

    EXAMPLE::

        sage: from sage.graphs.graph_input import from_graph6
        sage: g = Graph()
        sage: from_graph6(g, 'IheA@GUAo')
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True
    """
    from generic_graph_pyx import length_and_string_from_graph6, binary_string_from_graph6

    if not isinstance(g6_string, str):
        raise ValueError('If input format is graph6, then g6_string must be a string.')
    n = g6_string.find('\n')
    if n == -1:
        n = len(g6_string)
    ss = g6_string[:n]
    n, s = length_and_string_from_graph6(ss)
    m = binary_string_from_graph6(s, n)
    expected = n*(n-1)/2 + (6 - n*(n-1)/2)%6
    if len(m) > expected:
        raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too long."%(ss,n))
    elif len(m) < expected:
        raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too short."%(ss,n))
    G.add_vertices(range(n))
    k = 0
    for i in xrange(n):
        for j in xrange(i):
            if m[k] == '1':
                G._backend.add_edge(i, j, None, False)
            k += 1

def from_sparse6(G, g6_string):
    r"""
    Fill ``G`` with the data of a sparse6 string.

    INPUT:

    - ``G`` -- a graph

    - ``g6_string`` -- a sparse6 string

    EXAMPLE::

        sage: from sage.graphs.graph_input import from_sparse6
        sage: g = Graph()
        sage: from_sparse6(g, ':I`ES@obGkqegW~')
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True
    """
    from generic_graph_pyx import length_and_string_from_graph6, int_to_binary_string
    from math import ceil, floor
    from sage.misc.functional import log
    n = g6_string.find('\n')
    if n == -1:
        n = len(g6_string)
    s = g6_string[:n]
    n, s = length_and_string_from_graph6(s[1:])
    if n == 0:
        edges = []
    else:
        k = int(ceil(log(n,2)))
        ords = [ord(i) for i in s]
        if any(o > 126 or o < 63 for o in ords):
            raise RuntimeError("The string seems corrupt: valid characters are \n" + ''.join([chr(i) for i in xrange(63,127)]))
        bits = ''.join([int_to_binary_string(o-63).zfill(6) for o in ords])
        b = []
        x = []
        for i in xrange(int(floor(len(bits)/(k+1)))):
            b.append(int(bits[(k+1)*i:(k+1)*i+1],2))
            x.append(int(bits[(k+1)*i+1:(k+1)*i+k+1],2))
        v = 0
        edges = []
        for i in xrange(len(b)):
            if b[i] == 1:
                v += 1
            if x[i] > v:
                v = x[i]
            else:
                if v < n:
                    edges.append((x[i],v))
    G.add_vertices(range(n))
    G.add_edges(edges)

def from_dig6(G, dig6_string):
    r"""
    Fill ``G`` with the data of a dig6 string.

    INPUT:

    - ``G`` -- a graph

    - ``dig6_string`` -- a dig6 string

    EXAMPLE::

        sage: from sage.graphs.graph_input import from_dig6
        sage: g = DiGraph()
        sage: from_dig6(g, digraphs.Circuit(10).dig6_string())
        sage: g.is_isomorphic(digraphs.Circuit(10))
        True
    """
    from generic_graph_pyx import length_and_string_from_graph6, binary_string_from_dig6
    if not isinstance(dig6_string, str):
        raise ValueError('If input format is dig6, then dig6_string must be a string.')
    n = dig6_string.find('\n')
    if n == -1:
        n = len(dig6_string)
    ss = dig6_string[:n]
    n, s = length_and_string_from_graph6(ss)
    m = binary_string_from_dig6(s, n)
    expected = n**2
    if len(m) > expected:
        raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too long."%(ss,n))
    elif len(m) < expected:
        raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too short."%(ss,n))
    G.add_vertices(range(n))
    k = 0
    for i in xrange(n):
        for j in xrange(n):
            if m[k] == '1':
                G._backend.add_edge(i, j, None, True)
            k += 1

def from_seidel_adjacency_matrix(G, M):
    r"""
    Fill ``G`` with the data of a Seidel adjacency matrix.

    INPUT:

    - ``G`` -- a graph

    - ``M`` -- a Seidel adjacency matrix

    EXAMPLE::

        sage: from sage.graphs.graph_input import from_seidel_adjacency_matrix
        sage: g = Graph()
        sage: from_seidel_adjacency_matrix(g, graphs.PetersenGraph().seidel_adjacency_matrix())
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True
    """
    from sage.matrix.matrix import is_Matrix
    from sage.rings.integer_ring import ZZ
    assert is_Matrix(M)

    if M.base_ring() != ZZ:
        try:
            M = M.change_ring(ZZ)
        except TypeError:
            raise ValueError("Graph's Seidel adjacency matrix must"+
                             " have only 0,1,-1 integer entries")

    if M.is_sparse():
        entries = set(M[i,j] for i,j in M.nonzero_positions())
    else:
        entries = set(M.list())

    if any(e <  -1 or e > 1 for e in entries):
        raise ValueError("Graph's Seidel adjacency matrix must"+
                         " have only 0,1,-1 integer entries")
    if any(i==j for i,j in M.nonzero_positions()):
        raise ValueError("Graph's Seidel adjacency matrix must"+
                         " have 0s on the main diagonal")
    if not M.is_symmetric():
        raise ValueError("Graph's Seidel adjacency matrix must"+
                         " be symmetric")
    G.add_vertices(range(M.nrows()))
    e = []
    for i,j in M.nonzero_positions():
       if i <= j and M[i,j] < 0:
                e.append((i,j))
    G.add_edges(e)

def from_adjacency_matrix(G, M, loops=False, multiedges=False, weighted=False):
    r"""
    Fill ``G`` with the data of an adjacency matrix.

    INPUT:

    - ``G`` -- a graph

    - ``M`` -- an adjacency matrix

    - ``loops``, ``multiedges``, ``weighted`` (booleans) -- whether to consider
      the graph as having loops, multiple edges, or weights. Set to ``False`` by default.

    EXAMPLE::

        sage: from sage.graphs.graph_input import from_adjacency_matrix
        sage: g = Graph()
        sage: from_adjacency_matrix(g, graphs.PetersenGraph().adjacency_matrix())
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True
    """
    from sage.matrix.matrix import is_Matrix
    from sage.rings.integer_ring import ZZ
    assert is_Matrix(M)
    # note: the adjacency matrix might be weighted and hence not
    # necessarily consists of integers
    if not weighted and M.base_ring() != ZZ:
        try:
            M = M.change_ring(ZZ)
        except TypeError:
            if weighted is False:
                raise ValueError("Non-weighted graph's"+
                " adjacency matrix must have only nonnegative"+
                " integer entries")
            weighted = True

    if M.is_sparse():
        entries = set(M[i,j] for i,j in M.nonzero_positions())
    else:
        entries = set(M.list())

    if not weighted and any(e < 0 for e in entries):
        if weighted is False:
            raise ValueError("Non-weighted digraph's"+
            " adjacency matrix must have only nonnegative"+
            " integer entries")
        weighted = True
        if multiedges is None: multiedges = False
    if weighted is None:
        weighted = False

    if multiedges is None:
        multiedges = ((not weighted) and any(e != 0 and e != 1 for e in entries))

    if not loops and any(M[i,i] for i in xrange(M.nrows())):
        if loops is False:
            raise ValueError("Non-looped digraph's adjacency"+
            " matrix must have zeroes on the diagonal.")
        loops = True
    if loops is None:
        loops = False
    G.allow_loops(loops, check=False)
    G.allow_multiple_edges(multiedges, check=False)
    G.add_vertices(range(M.nrows()))
    e = []
    if weighted:
        for i,j in M.nonzero_positions():
            if i <= j:
                e.append((i,j,M[i][j]))
    elif multiedges:
        for i,j in M.nonzero_positions():
            if i <= j:
                e += [(i,j)]*int(M[i][j])
    else:
        for i,j in M.nonzero_positions():
            if i <= j:
                e.append((i,j))
    G.add_edges(e)
    G._weighted = weighted

def from_incidence_matrix(G, M, loops=False, multiedges=False, weighted=False):
    r"""
    Fill ``G`` with the data of an incidence matrix.

    INPUT:

    - ``G`` -- a graph

    - ``M`` -- an incidence matrix

    - ``loops``, ``multiedges``, ``weighted`` (booleans) -- whether to consider
      the graph as having loops, multiple edges, or weights. Set to ``False`` by default.

    EXAMPLE::

        sage: from sage.graphs.graph_input import from_incidence_matrix
        sage: g = Graph()
        sage: from_incidence_matrix(g, graphs.PetersenGraph().incidence_matrix())
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True
    """
    from sage.matrix.matrix import is_Matrix
    assert is_Matrix(M)

    oriented = any(M[pos] < 0 for pos in M.nonzero_positions(copy=False))

    positions = []
    for i in range(M.ncols()):
        NZ = M.nonzero_positions_in_column(i)
        if len(NZ) == 1:
            if oriented:
                raise ValueError("Column {} of the (oriented) incidence "
                                 "matrix contains only one nonzero value".format(i))
            elif M[NZ[0],i] != 2:
                raise ValueError("Each column of a non-oriented incidence "
                                 "matrix must sum to 2, but column {} does not".format(i))
            if loops is None:
                loops = True
            positions.append((NZ[0],NZ[0]))
        elif len(NZ) != 2 or \
             (oriented and not ((M[NZ[0],i] == +1 and M[NZ[1],i] == -1) or \
                                (M[NZ[0],i] == -1 and M[NZ[1],i] == +1))) or \
             (not oriented and (M[NZ[0],i] != 1 or M[NZ[1],i] != 1)):
            msg  = "There must be one or two nonzero entries per column in an incidence matrix. "
            msg += "Got entries {} in column {}".format([M[j,i] for j in NZ], i)
            raise ValueError(msg)
        else:
            positions.append(tuple(NZ))

    if weighted   is None: G._weighted  = False
    if multiedges is None:
        total = len(positions)
        multiedges = (len(set(positions)) < total  )
    G.allow_loops(False if loops is None else loops, check=False)
    G.allow_multiple_edges(multiedges, check=False)
    G.add_vertices(range(M.nrows()))
    G.add_edges(positions)

def from_dict_of_dicts(G, M, loops=False, multiedges=False, weighted=False, convert_empty_dict_labels_to_None=False):
    r"""
    Fill ``G`` with the data of a dictionary of dictionaries.

    INPUT:

    - ``G`` -- a graph

    - ``M`` -- a dictionary of dictionaries.

    - ``loops``, ``multiedges``, ``weighted`` (booleans) -- whether to consider
      the graph as having loops, multiple edges, or weights. Set to ``False`` by default.

    - ``convert_empty_dict_labels_to_None`` (boolean) -- whether to adjust for
      empty dicts instead of None in NetworkX default edge labels.

    EXAMPLE::

        sage: from sage.graphs.graph_input import from_dict_of_dicts
        sage: g = Graph()
        sage: from_dict_of_dicts(g, graphs.PetersenGraph().to_dictionary(edge_labels=True))
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True
    """
    if not all(isinstance(M[u], dict) for u in M):
        raise ValueError("Input dict must be a consistent format.")

    if not loops and any(u in neighb for u,neighb in M.iteritems()):
        if loops is False:
            u = next(u for u,neighb in M.iteritems() if u in neighb)
            raise ValueError("The graph was built with loops=False but input M has a loop at {}.".format(u))
        loops = True
    if loops is None:
        loops = False

    if weighted is None: G._weighted = False
    for u in M:
        for v in M[u]:
            if hash(u) > hash(v):
                if v in M and u in M[v]:
                    if M[u][v] != M[v][u]:
                        raise ValueError("Dict does not agree on edge (%s,%s)"%(u,v))
                    continue
            if multiedges is not False and not isinstance(M[u][v], list):
                if multiedges is None: multiedges = False
                if multiedges:
                    raise ValueError("Dict of dicts for multigraph must be in the format {v : {u : list}}")
    if multiedges is None and len(M) > 0:
        multiedges = True
    G.allow_loops(loops, check=False)
    G.allow_multiple_edges(multiedges, check=False)
    verts = set().union(M.keys(), *M.values())
    G.add_vertices(verts)
    if convert_empty_dict_labels_to_None:
        for u in M:
            for v in M[u]:
                if hash(u) <= hash(v) or v not in M or u not in M[v]:
                    if multiedges:
                        for l in M[u][v]:
                            G._backend.add_edge(u,v,l,False)
                    else:
                        G._backend.add_edge(u,v,M[u][v] if M[u][v] != {} else None,False)
    else:
        for u in M:
            for v in M[u]:
                if hash(u) <= hash(v) or v not in M or u not in M[v]:
                    if multiedges:
                        for l in M[u][v]:
                            G._backend.add_edge(u,v,l,False)
                    else:
                        G._backend.add_edge(u,v,M[u][v],False)

def from_dict_of_lists(G, D, loops=False, multiedges=False, weighted=False):
    r"""
    Fill ``G`` with the data of a dictionary of lists.

    INPUT:

    - ``G`` -- a graph

    - ``D`` -- a dictionary of lists.

    - ``loops``, ``multiedges``, ``weighted`` (booleans) -- whether to consider
      the graph as having loops, multiple edges, or weights. Set to ``False`` by default.

    EXAMPLE::

        sage: from sage.graphs.graph_input import from_dict_of_lists
        sage: g = Graph()
        sage: from_dict_of_lists(g, graphs.PetersenGraph().to_dictionary())
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True
    """
    if not all(isinstance(D[u], list) for u in D):
        raise ValueError("Input dict must be a consistent format.")

    verts = set().union(D.keys(),*D.values())
    if loops is None or loops is False:
        for u in D:
            if u in D[u]:
                if loops is None:
                    loops = True
                elif loops is False:
                    u = next(u for u,neighb in D.iteritems() if u in neighb)
                    raise ValueError("The graph was built with loops=False but input D has a loop at {}.".format(u))
                break
        if loops is None:
            loops = False
    if weighted is None: G._weighted = False
    for u in D:
        if len(set(D[u])) != len(D[u]):
            if multiedges is False:
                v = next((v for v in D[u] if D[u].count(v) > 1))
                raise ValueError("Non-multigraph got several edges (%s,%s)"%(u,v))
            if multiedges is None:
                multiedges = True
    if multiedges is None: multiedges = False
    G.allow_loops(loops, check=False)
    G.allow_multiple_edges(multiedges, check=False)
    G.add_vertices(verts)
    for u in D:
        for v in D[u]:
            if (multiedges or hash(u) <= hash(v) or
                v not in D or u not in D[v]):
                G._backend.add_edge(u,v,None,False)

from sage.misc.rest_index_of_methods import gen_rest_table_index
import sys
__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=gen_rest_table_index(sys.modules[__name__]))
