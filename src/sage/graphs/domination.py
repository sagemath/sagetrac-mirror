r"""
Functions related to domination.

This file contains varied functions related to domination as well as an
implementation of the algorithm described in _Enumerating minimal
dominating sets in $K_t$-free graphs and variants_ by Marthe Bonamy,
Oscar Defrain, Marc Heinrich, Michał Pilipczuk, and Jean-Florent Raymond
:arxiv:`1810.00789` to enumerate the minimal dominating sets of a graph.

EXAMPLES::

to be added

AUTHORS:

- Jean-Florent Raymond (2019-03-04): initial version
"""

# ****************************************************************************
#       Copyright (C) 2019 Jean-Florent Raymond  <raymond@tu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy

def _parent(G, dom, V_prev):
    r'''
    Return an subset of dom that is irredundant in V_prev.

    For internal use.

    INPUT:
    
    - ``dom`` -- an iterable of vertices of ``self``
    - ``V_prev`` -- an iterable of vertices of ``self``

    OUTPUT:

    Return the list obtained from ``dom`` by iteratively removing those
    vertices of mininum index that have no private neighbor in V_prev.

    EXAMPLES:

        sage: from sage.graphs.domination import _parent
        sage: G = graphs.PathGraph(4)
        sage: G.add_vertices([4, 5])
        sage: G.add_edges([(4, 1), (5, 2)])
        sage: _parent(G, [0, 2, 4, 5], [1, 2])
        [4, 5]
        sage: _parent(G, [0, 2, 4, 5], [1, 3])
        [2]

    .. WARNING:
    
    We assume that vertices are sortable (i.e. they can be compared).
    '''

    # The list where we search vertices
    D_start = sorted(dom, reverse=True)

    # The list to be output at the end, that we construct:
    D_end = []

    while D_start:
        v = D_start.pop() # element of min index
        priv = set(G.neighbor_iterator(v, closed=True))
        # We remove the vertices already dominated
        # by other vertices of (D_end union D_start)
        priv.difference_update(*(G.neighbor_iterator(u, closed=True)
                                 for u in D_start if u != v))
        priv.difference_update(*(G.neighbor_iterator(u, closed=True)
                                 for u in D_end if u != v))
        # Now priv is the private neighborhood of v
        # in G wrt D_start + D_end
        if priv.intersection(V_prev) != set():
            # if v has a private in V_prev, we keep it
            D_end.append(v)

    return D_end


def _peel(G, A):
    r'''
    Return a peeling of a vertex iterable of a graph.

    For internal use.
    Given a graph `G` and a subset `A` of its vertices, a peeling
    of `(G,A)` is a list $[(u_0, V_0), \dots, (u_{p+1}, V_{p+1})]$ such
    that $u_0$ and $u_{p+1}$ are `None`, $V_0$ is the empty set,
    $V_{p+1} = V(G)$, $V_p = A$ and for every $i \in $\{1, \dots, p\}$,
    $V_{i-1} = V_i \setminus N[v_i]$, for some $u_i\in V_i$.

    INPUT:
    
    - ``G`` -- a graph
    - ``A`` -- a set of vertices of ``G``

    OUTPUT:

    A peeling of `(G,A)`.

    EXAMPLES:

        sage: from sage.graphs.domination import _peel
        sage: G = Graph(10); _peel(G, {0, 1, 2, 3, 4})
        [(None, set()),
        (4, {4}),
        (3, {3, 4}),
        (2, {2, 3, 4}),
        (1, {1, 2, 3, 4}),
        (0, {0, 1, 2, 3, 4}),
        (None, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9})]


        sage: from sage.graphs.domination import _peel
        sage: G = graphs.PathGraph(10); _peel(G, set((i for i in range(10) if i%2==0)))
        [(None, set()),
        (8, {8}),
        (6, {6, 8}),
        (4, {4, 6, 8}),
        (2, {2, 4, 6, 8}),
        (0, {0, 2, 4, 6, 8}),
        (None, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9})]

    '''
    Acomp = set(G)
    Acomp.difference_update(A) # Acomp  = V - A

    peeling = [(None, set(G))]
    H = copy(G)
    H.delete_vertices(list(Acomp))
    del Acomp
    
    while H:
        ui = next(H.vertex_iterator())  # pick some vertex of H
        Vi = set(H.vertex_iterator())
        peeling.append((ui, Vi))
        H.delete_vertices(H.neighbor_iterator(ui, closed=True))
    peeling.append((None, set()))
    peeling.reverse()
    return peeling


def _cand_ext_enum(G, dom, u_next, V_next):
    r'''
    Return an iterator over the candidate extensions of ``dom`` to ``V_next``.

    For internal use. (See description below.)
    
    INPUT:

    - `G` -- a graph
    - `dom` -- an iterable over some vertices of `G`
    - `u_next` -- a vertex of `G`
    - `V_next` -- an iterable over some vertices of `G`
    
    OUTPUT:

    An iterator over all sets `X` that dominate
    $N(u_next) \cap V_next \setminus N[dom]$ and are inclusion-wise
    minimal (hereafter called candidate extensions of ``dom`` to ``V_next``).

    Intuitively, we have a partial (irredundant) dominating set `dom` and we
    want to extend it to `V_next`. For this we enumerate all the minimal
    dominating sets of $N(u_next) \cap V_next \setminus N[dom]$,
    which will contain all the extensions of `dom` into a minimal dominating
    set of `V_next`. (Note that by doing so we also enumerate sets `X` such
    that $dom \cup X$ is dominating `V_next` but is not minimal.)
    '''

    def _aux_with_rep(G, dom, u_next, V_next):
        # Auxilliary routine. Does the same as _cand_ext_enum,
        # except that the same set may be output several times
        # (up to |G| times).
        #
        # In order to later remove duplicates, we here output pairs
        # (dom,i) where dom is the output candidate extension and i counts
        # how many elements have already been output.
        #

        dominated_by_dom = set().union(
            (G.neighbor_iterator(u, closed=True) for u in dom))
        # S = neighbors of u_next in V_next that are not yet dominated:
        S = set().intersection(G.neighbor_iterator(u_next),
                             V_next - dominated_by_dom)
        del dominated_by_dom

        # True iff u_next is dominated by dom:
        u_next_dom_by_dom = any(
            v in dom for v in G.neighbor_iterator(u_next, closed=True))

        if u_next_dom_by_dom:
            # In this case, u_next is already dominated by dom,
            # so only S has to be dominated.
            # We enumerate the minimal DSs of the bicolored graph G(S):

            cand_ext_index = 0
            for X in minimal_dominating_sets(G, S):
                yield (X, cand_ext_index)
                cand_ext_index += 1

        elif not S:
            # In this case, only u_next has to be dominated
            cand_ext_index = 0
            for w in G.neighbor_iterator(u_next, closed=True):
                # Notice that the case w = u_next is included
                yield (set(w), cand_ext_index)
                cand_ext_index += 1

        else:
            # In this case, both u_next and S have to be dominated

            S_plus = copy(S)
            S_plus.add(u_next)  # S_plus = S + {u_next}

            yield ({u_next}, 0)  # The trivial extension
            # Start from 1 because we already output the 0-th elt:
            cand_ext_index = 1

            for w in G.neighbor_iterator(u_next):
                # Note that w never belongs to dom,
                # as we are not in the first case of the if statement

                S_minus = S
                S_minus.difference_update(G.neighbor_iterator(w, closed=True))
                # S_minus: vertices of S that still need to be
                # dominated, assuming w is included in the DS

                for Q in minimal_dominating_sets(G, S_minus):
                    # Closed neighborhood of Q:
                    NQ = set().union(
                        (self.neighbor_iterator(u, closed=True) for u in Q))
                    Nw_minus = set.intersection(
                        set(G.neighbor_iterator(w, closed=True)), S_plus)
                    if not NQ >= Nw_minus:
                        # If Nw_minus is not included in i.e. if w has
                        # a private neighbor in V_next wrt Q + {w}:
                        sQ = set(Q)
                        sQ.add(w)
                        yield (sQ, cand_ext_index)
                        cand_ext_index += 1
    #
    # End of aux_with_rep routine

    # Here we use aux_with_rep twice to enumerate the minimal
    # dominating sets while avoiding repeated outputs
    for (X, i) in _aux_with_rep(G, dom, u_next, V_next):
        for (Y, j) in _aux_with_rep(G, dom, u_next, V_next):
            if j >= i:
                # This is the first time we meet X: we output it
                yield X
                break
            elif set(Y) == set(X):
                # X has already been output in the past: we ignore it
                break


def minimal_dominating_sets(G, to_dominate=None):
    r'''
    Return an iterator over the minimal dominating sets of the graph.

    INPUT:

    - `G` -- a graph
    - `to_dominate` -- vertex iterable or None (default: `None`)

    OUTPUT:

    An iterator over the inclusion-minimal sets of vertices of `G`
    that dominate `to_dominate`.

    ALGORITHM:

    The algorithm is described in :arxiv:`1810.00789`.

    EXAMPLES:

        sage: G = graphs.ButterflyGraph()
        sage: sorted(list(minimal_dominating_sets(G)))
        [{0, 1}, {1, 3}, {0, 2}, {2, 3}, {4}]
        sage: sorted(list(minimal_dominating_sets(G, [0,3])))
        [{0}, {3}, {4}]
        sage: sorted(list(minimal_dominating_sets(G, [4])))
        [{4}, {0}, {1}, {2}, {3}]

        sage: sorted(list(minimal_dominating_sets(graphs.PetersenGraph())))
        [{0, 2, 6},
        {0, 9, 3},
        {0, 8, 7},
        {1, 3, 7},
        {1, 4, 5},
        {8, 1, 9},
        {8, 2, 4},
        {9, 2, 5},
        {3, 5, 6},
        {4, 6, 7},
        {0, 8, 2, 9},
        {0, 3, 6, 7},
        {1, 3, 5, 9},
        {8, 1, 4, 7},
        {2, 4, 5, 6},
        {0, 1, 2, 3, 4},
        {0, 1, 2, 5, 7},
        {0, 1, 4, 6, 9},
        {0, 1, 5, 6, 8},
        {0, 3, 4, 5, 8},
        {0, 4, 5, 7, 9},
        {1, 2, 3, 6, 8},
        {1, 2, 6, 7, 9},
        {2, 3, 4, 7, 9},
        {2, 3, 5, 7, 8},
        {3, 4, 6, 8, 9},
        {5, 6, 7, 8, 9}]

    TESTS:

        sage: list(minimal_dominating_sets(Graph())):
        [set()]

    ::

        sage: from sage.combinat.subset import Subsets
        sage: def minimal_dominating_sets_naive(G):
        ....:     return (S for S in Subsets(G.vertices())
        ....:             if not(G.is_redundant(S)) and G.is_dominating(S))
        sage: G = graphs.RandomGNP(5, 0.5)
        sage: sorted(list(minimal_dominating_sets(G))) == sorted(list(minimal_dominating_sets_naive(G))) # long time
        True

    ::

        sage: findstat([(G, sum(1 for _ in minimal_dominating_sets(G))) for n in range(6) for G in graphs(n)], depth=0) # optional - internet
        0: (St001302: The number of minimally dominating sets of vertices of a graph., [], 52)

    .. WARNING:
    
    We assume that vertices are sortable (i.e. they can be compared).
    '''

    def tree_search(H, plng, dom, i):
        # Internal routine
        # Recursively generates the leaves descendant of (dom,i)
        # H: graph; plng: peeling

        if i == len(plng) - 2:
            # we reached a leaf, i.e. dom is a minimal DS of vertices_to_dominate
            # '-2' because the last cell of plng is used
            # for the vertices of H - vertices_to_dominate
            yield dom
            return

        u_next, V_next = plng[i+1]

        if H.is_dominating(dom, V_next):  # if dom dominates V_next
            # then dom is its unique extension: we recurse on it
            for Di in tree_search(H, plng, dom, i + 1):
                yield Di
            return

        # For every candidate extension
        for can_ext in _cand_ext_enum(H, dom, u_next, V_next):

            # We complete dom with can_ext -> canD
            canD = set.union(set(can_ext), set(dom))

            if (not H.is_redundant(canD, V_next)) and set(dom) == set(_parent(H, canD, plng[i][1])):
                # If canD is a legitimate child of dom and is not
                # redundant, we recurse on it:
                for Di in tree_search(H, plng, canD, i + 1):
                    yield Di
    ##
    # end of tree-search routine

    if to_dominate is None:
        vertices_to_dominate = set(G)
    else:
        vertices_to_dominate = set(to_dominate)

    if not vertices_to_dominate:  # base case: vertices_to_dominate is empty
        yield set()  # the empty set/list is the only minimal DS of the empty set of vertex
        return

    peeling = _peel(G, vertices_to_dominate)

    for dom in tree_search(G, peeling, set(), 0):
        # we generate the leaves of the search tree that are descendant of (empty set, 0)
        yield dom

