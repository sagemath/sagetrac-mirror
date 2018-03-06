r"""
mutation_class

This file contains helper functions for compute the mutation class of a cluster algebra or quiver.

For the compendium on the cluster algebra and quiver package see [MS2011]_

AUTHORS:

- Christian Stump
"""

# ****************************************************************************
#       Copyright (C) 2018 Gregg Musiker <musiker@math.mit.edu>
#                          Christian Stump <christian.stump@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ****************************************************************************
import time
from copy import copy
from sage.matrix.special import zero_matrix
from sage.matrix.constructor import Matrix
from sage.graphs.digraph import DiGraph
from sage.rings.infinity import Infinity
from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree, get_orbits

def _matrix_mutation_class_iter(M, depth=Infinity, show_depth=False, return_paths=False, up_to_equivalence=True, sink_source=False):
    """
    Returns an iterator for mutation class of matrix with respect to several parameters.

    INPUT:

    - ``M`` -- a skew-symmetrizable matrix representing a valued quiver, possibly with coefficients
    - ``depth`` -- a positive integer or infinity specifying (roughly) how many steps away from the initial seed to mutate
    - ``show_depth`` -- if True, indicates that a running count of the depth is to be displayed
    - ``up_to_equivalence``  -- if True, only one digraph for each graph-isomorphism class is recorded
    - ``sink_source`` -- if True, only mutations at sinks or sources are applied

    EXAMPLES::

        sage: tba
    """
    timer = time.time()
    depth_counter = 0
    n = M.ncols()
    m = M.nrows() - n

    if up_to_equivalence:
        Mhash = matrix_canonical_hash(M, n, m)
    else:
        M.set_immutable()
        Mhash = hash(M)

    have_seen = set([Mhash])
    to_check  = [ (M, []) ]

    if return_paths:
        yield (M, [])
    else:
        yield M

    if show_depth:
        timer2 = time.time()
        dc = str(depth_counter)
        dc += ' ' * (5-len(dc))
        nr = str(len(have_seen))
        nr += ' ' * (10-len(nr))
        print("Depth: %s found: %s time: %.2f s" % (dc, nr, timer2 - timer))

    while to_check and depth_counter < depth:
        to_check_new = []
        for (M, history) in to_check:
            for i in range(n):
                if not sink_source or _matrix_is_sink_source(M, i):
                    M_new = _fast_copy(M, n, m)
                    M_new.mutate(i)

                    if up_to_equivalence:
                        Mhash_new = matrix_canonical_hash(M_new, n, m)
                    else:
                        M_new.set_immutable()
                        Mhash_new = hash(M_new)

                    if Mhash_new not in have_seen:
                        have_seen.add(Mhash_new)
                        to_check_new.append((M_new, history+[i]))
                        if return_paths:
                            yield to_check_new[-1][:2]
                        else:
                            yield to_check_new[-1][0]
        to_check = to_check_new

        depth_counter += 1
        if show_depth and bool(to_check):
            timer2 = time.time()
            dc = str(depth_counter)
            dc += ' ' * (5-len(dc))
            nr = str(len(have_seen))
            nr += ' ' * (10-len(nr))
            print("Depth: %s found: %s time: %.2f s" % (dc, nr, timer2 - timer))

def matrix_canonical_hash(M, n, m):
    dg,partition = _matrix_to_unlabelled_digraph(M, n, m)
    dg_canon = dg.canonical_label(partition=partition, algorithm="bliss", return_graph=False)
    return hash(tuple(sorted(dg_canon)))

def _matrix_to_unlabelled_digraph(M, n, m):
    dg = DiGraph(sparse=True)
    dg.add_vertices(range(n+m))

    edge_labels = {}
    new_vertex  = n+m

    new_partition = []

    for i, j in M.nonzero_positions():
        if i < n:
            a, b = M[i, j], M[j, i]
        else:
            a, b = M[i, j], -M[i, j]
        if a > 0:
            if a == 1 and b == -1:
                dg._backend.add_edge(i,j,None,True)
            else:
                try:
                    x = edge_labels[(a,b)]
                except KeyError:
                    x = len(new_partition)
                    edge_labels[(a,b)] = x
                    new_partition.append([])
                finally:
                    dg.add_vertex(new_vertex)
                    dg._backend.add_edge(i,new_vertex,None,True)
                    dg._backend.add_edge(new_vertex,j,None,True)
                    new_partition[x].append(new_vertex)
                    new_vertex += 1
        elif i >= n:
            if a == -1 and b == 1:
                dg._backend.add_edge(j,i,None,True)
            else:
                a = -a
                b = -b
                try:
                    x = edge_labels[(a,b)]
                except KeyError:
                    x = len(new_partition)
                    edge_labels[(a,b)] = x
                    new_partition.append([])
                finally:
                    dg.add_vertex(new_vertex)
                    dg._backend.add_edge(j,new_vertex,None,True)
                    dg._backend.add_edge(new_vertex,i,None,True)
                    new_partition[x].append(new_vertex)
                    new_vertex += 1
    partition = [list(range(n))]
    if m > 0:
        partition.append(list(range(n,n+m)))
    if new_partition:
        partition.extend(new_partition)
    return dg, partition

def _fast_copy(M, n, m):
    Mnew = zero_matrix(n+m,n)
    for (i,j),val in M.dict().iteritems():
        Mnew[i,j] = val
    return Mnew

def _matrix_is_sink_source(M, i):
    """
    Returns True iff the quiver of M has a sink or a source at vertex v.

    INPUT:

    - ``M`` -- a matrix
    - ``i`` -- an index

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _matrix_is_sink_source
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: tba
    """
    inedge  = False
    outedge = False
    for (a,b), val in M.dict().iteritems():
        if a == i:
            if val > 0:
                inedge = True
            else:
                outedge = True
        elif b == i:
            if val > 0:
                outedge = True
            else:
                inedge = True
        if inedge and outedge:
            return False
    return True
