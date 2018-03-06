r"""
mutation_class

This file contains helper functions for compute the mutation class of a cluster algebra or quiver.

For the compendium on the cluster algebra and quiver package see [MS2011]_

AUTHORS:

- Gregg Musiker
- Christian Stump
"""

# ****************************************************************************
#       Copyright (C) 2011 Gregg Musiker <musiker@math.mit.edu>
#                          Christian Stump <christian.stump@univie.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import print_function
from six.moves import range

import time
from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree, get_orbits
from sage.rings.all import ZZ, infinity
from sage.graphs.all import DiGraph
from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import _edge_list_to_matrix


def _principal_part( mat ):
    """
    Return the principal part of a matrix.

    INPUT:

    - ``mat`` -- a matrix with at least as many rows as columns

    OUTPUT:

    The top square part of the matrix ``mat``.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _principal_part
        sage: M = Matrix([[1,2],[3,4],[5,6]]); M
        [1 2]
        [3 4]
        [5 6]
        sage: _principal_part(M)
        [1 2]
        [3 4]
    """
    n, m = mat.ncols(), mat.nrows() - mat.ncols()
    if m < 0:
        raise ValueError('The input matrix has more columns than rows.')
    elif m == 0:
        return mat
    else:
        return mat.submatrix(0, 0, n, n)

def _matrix_to_digraph( M ):
    """
    Returns the digraph obtained from the matrix ``M``.
    In order to generate a quiver, we assume that ``M`` is skew-symmetrizable.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _matrix_to_digraph
        sage: _matrix_to_digraph(matrix(3,[0,1,0,-1,0,-1,0,1,0]))
        Digraph on 3 vertices
    """
    n = M.ncols()

    dg = DiGraph(sparse=True)
    for i, j in M.nonzero_positions():
        if i >= n:
            a, b = M[i, j], -M[i, j]
        else:
            a, b = M[i, j], M[j, i]
        if a > 0:
            dg._backend.add_edge(i,j,(a,b),True)
        elif i >= n:
            dg._backend.add_edge(j,i,(-a,-b),True)
    if dg.order() < M.nrows():
        for i in [ index for index in range(M.nrows()) if index not in dg ]:
            dg.add_vertex(i)
    return dg

def _dg_canonical_form( dg, n, m ):
    """
    Turns the digraph ``dg`` into its canonical form, and returns the corresponding isomorphism, and the vertex orbits of the automorphism group.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _dg_canonical_form
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['B',4]).digraph(); dg.edges()
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -2))]
        sage: _dg_canonical_form(dg,4,0); dg.edges()
        ({0: 0, 1: 3, 2: 1, 3: 2}, [[0], [3], [1], [2]])
        [(0, 3, (1, -1)), (1, 2, (1, -2)), (1, 3, (1, -1))]
    """
    vertices = [ v for v in dg ]
    if m > 0:
        partition = [ vertices[:n], vertices[n:] ]
    else:
        partition = [ vertices ]
    partition_add, edges = _graph_without_edge_labels(dg,vertices)
    partition += partition_add
    automorphism_group, obsolete, iso = search_tree(dg, partition=partition, lab=True, dig=True, certificate=True)
    orbits = get_orbits( automorphism_group, n+m )
    orbits = [ [ iso[i] for i in orbit] for orbit in orbits ]
    for v in iso.keys():
        if v >= n+m:
            del iso[v]
            v1,v2,label1 = next(dg._backend.iterator_in_edges([v],True))
            w1,w2,label2 = next(dg._backend.iterator_out_edges([v],True))
            dg._backend.del_edge(v1,v2,label1,True)
            dg._backend.del_edge(w1,w2,label2,True)
            dg._backend.del_vertex(v)
            add_index = True
            index = 0
            while add_index:
                l = partition_add[index]
                if v in l:
                    dg._backend.add_edge(v1,w2,edges[index],True)
                    add_index = False
                index += 1
    dg._backend.relabel( iso, True )
    return iso, orbits

def _digraph_to_dig6( dg, hashable=False ):
    """
    Returns the dig6 and edge data of the digraph dg.

    INPUT:

    - ``dg`` -- a digraph
    - ``hashable`` -- (Boolean; optional; default:False) if ``True``, the edge labels are turned into a dict.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_to_dig6
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',4]).digraph()
        sage: _digraph_to_dig6(dg)
        ('COD?', {})
    """
    dig6 = dg.dig6_string()
    D = {}
    for E in dg._backend.iterator_in_edges(dg,True):
        if E[2] != (1,-1):
            D[ (E[0],E[1]) ] = E[2]
    if hashable:
        D = tuple( sorted( D.items() ) )
    return (dig6,D)


def _dig6_to_digraph( dig6 ):
    """
    Returns the digraph obtained from the dig6 and edge data.

    INPUT:

    - ``dig6`` -- a pair ``(dig6, edges)`` where ``dig6`` is a string encoding a digraph and ``edges`` is a dict or tuple encoding edges

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_to_dig6
        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _dig6_to_digraph
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',4]).digraph()
        sage: data = _digraph_to_dig6(dg)
        sage: _dig6_to_digraph(data)
        Digraph on 4 vertices
        sage: _dig6_to_digraph(data).edges()
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -1))]
    """
    dig6, edges = dig6
    dg = DiGraph( dig6 )
    if not isinstance(edges, dict):
        edges = dict( edges )
    for edge in dg._backend.iterator_in_edges(dg,False):
        if edge in edges:
            dg.set_edge_label( edge[0],edge[1],edges[edge] )
        else:
            dg.set_edge_label( edge[0],edge[1], (1,-1) )
    return dg


def _dig6_to_matrix( dig6 ):
    """
    Return the matrix obtained from the dig6 and edge data.

    INPUT:

    - ``dig6`` -- a pair ``(dig6, edges)`` where ``dig6`` is a string
      encoding a digraph and ``edges`` is a dict or tuple encoding edges

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_to_dig6, _dig6_to_matrix
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',4]).digraph()
        sage: data = _digraph_to_dig6(dg)
        sage: _dig6_to_matrix(data)
        [ 0  1  0  0]
        [-1  0 -1  0]
        [ 0  1  0  1]
        [ 0  0 -1  0]
    """
    dg = _dig6_to_digraph(dig6)
    return _edge_list_to_matrix(dg.edges(), list(range(dg.order())), [])


def _dg_is_sink_source( dg, v ):
    """
    Returns True iff the digraph dg has a sink or a source at vertex v.

    INPUT:

    - ``dg`` -- a digraph
    - ``v`` -- a vertex of dg

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _dg_is_sink_source
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',[1,2],1]).digraph()
        sage: _dg_is_sink_source(dg, 0 )
        True
        sage: _dg_is_sink_source(dg, 1 )
        True
        sage: _dg_is_sink_source(dg, 2 )
        False
    """
    in_edges = [ edge for edge in dg._backend.iterator_in_edges([v],True) ]
    out_edges = [ edge for edge in dg._backend.iterator_out_edges([v],True) ]
    return not ( in_edges and out_edges )


def _graph_without_edge_labels(dg,vertices):
    """
    Replaces edge labels in dg other than ``(1,-1)`` by this edge label, and returns the corresponding partition of the edges.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _graph_without_edge_labels
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['B',4]).digraph(); dg.edges()
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -2))]
        sage: _graph_without_edge_labels(dg,range(3)); dg.edges()
        ([[5]], [(1, -2)])
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 5, (1, -1)), (5, 3, (1, -1))]
    """
    edges = [ edge for edge in dg._backend.iterator_in_edges(dg,True) ]
    edge_labels = sorted([ label for v1,v2,label in edges if not label == (1,-1)])
    i = 1
    while i < len(edge_labels):
        if edge_labels[i] == edge_labels[i-1]:
            edge_labels.pop(i)
        else:
            i += 1
    edge_partition = [[] for _ in range(len(edge_labels))]
    i = 0
    new_vertices = []
    for u,v,l in edges:
        while i in vertices or i in new_vertices:
            i += 1
        new_vertices.append(i)
        if not l == (1,-1):
            index = edge_labels.index(l)
            edge_partition[index].append(i)
            dg._backend.add_edge(u,i,(1,-1),True)
            dg._backend.add_edge(i,v,(1,-1),True)
            dg._backend.del_edge(u,v,l,True)
    return [a for a in edge_partition if a != []], edge_labels


def _has_two_cycles( dg ):
    """
    Returns True if the input digraph has a 2-cycle and False otherwise.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _has_two_cycles
        sage: _has_two_cycles( DiGraph([[0,1],[1,0]]))
        True
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: _has_two_cycles( ClusterQuiver(['A',3]).digraph() )
        False
    """
    edge_set = dg.edges(labels=False)
    for (v,w) in edge_set:
        if (w,v) in edge_set:
            return True
    return False


def _is_valid_digraph_edge_set( edges, frozen=0 ):
    """
    Returns True if the input data is the edge set of a digraph for a quiver (no loops, no 2-cycles, edge-labels of the specified format), and returns False otherwise.

    INPUT:

    - ``frozen`` -- (integer; default:0) The number of frozen vertices.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _is_valid_digraph_edge_set
        sage: _is_valid_digraph_edge_set( [[0,1,'a'],[2,3,(1,-1)]] )
        The given digraph has edge labels which are not integral or integral 2-tuples.
        False
        sage: _is_valid_digraph_edge_set( [[0,1,None],[2,3,(1,-1)]] )
        True
        sage: _is_valid_digraph_edge_set( [[0,1,'a'],[2,3,(1,-1)],[3,2,(1,-1)]] )
        The given digraph or edge list contains oriented 2-cycles.
        False
    """
    try:
        dg = DiGraph()
        dg.allow_multiple_edges(True)
        dg.add_edges( edges )

        # checks if the digraph contains loops
        if dg.has_loops():
            print("The given digraph or edge list contains loops.")
            return False

        # checks if the digraph contains oriented 2-cycles
        if _has_two_cycles( dg ):
            print("The given digraph or edge list contains oriented 2-cycles.")
            return False

        # checks if all edge labels are 'None', positive integers or tuples of positive integers
        if not all( i is None or ( i in ZZ and i > 0 ) or ( isinstance(i, tuple) and len(i) == 2 and i[0] in ZZ and i[1] in ZZ ) for i in dg.edge_labels() ):
            print("The given digraph has edge labels which are not integral or integral 2-tuples.")
            return False

        # checks if all edge labels for multiple edges are 'None' or positive integers
        if dg.has_multiple_edges():
            for e in set( dg.multiple_edges(labels=False) ):
                if not all( i is None or ( i in ZZ and i > 0 ) for i in dg.edge_label( e[0], e[1] ) ):
                    print("The given digraph or edge list contains multiple edges with non-integral labels.")
                    return False

        n = dg.order() - frozen
        if n < 0:
            print("The number of frozen variables is larger than the number of vertices.")
            return False

        if [ e for e in dg.edges(labels=False) if e[0] >= n] != []:
            print("The given digraph or edge list contains edges within the frozen vertices.")
            return False

        return True
    except Exception:
        print("Could not even build a digraph from the input data.")
        return False
