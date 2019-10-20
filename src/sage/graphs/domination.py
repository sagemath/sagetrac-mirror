# -*- coding: utf-8 -*-
r"""
Domination

This module implements methods related to the notion of domination in graphs,
and more precisely:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~minimal_dominating_sets` | Return an iterator over the minimal dominating sets of a graph.
    :meth:`~is_dominating` | Check whether ``dom`` is a dominating set of ``G``.
    :meth:`~is_redundant` | Check whether a ``dom`` has redundant vertices.
    :meth:`~private_neighbors` | Return the private neighbors of a vertex with repect to other vertices.


EXAMPLES:

We enumerate the minimal dominating sets of the 5-star graph::

    sage: g = graphs.StarGraph(5)
    sage: list(g.minimal_dominating_sets())
    [{0}, {1, 2, 3, 4, 5}]

Now only those that dominate the middle vertex::

    sage: list(g.minimal_dominating_sets([0]))
    [{0}, {1}, {2}, {3}, {4}, {5}]

Now the minimal dominating sets of the 5-path graph::

    sage: g = graphs.PathGraph(5)
    sage: list(g.minimal_dominating_sets())
    [{0, 2, 4}, {1, 4}, {0, 3}, {1, 3}]

We count the minimal dominating sets of the Petersen graph::

    sage: sum(1 for _ in graphs.PetersenGraph().minimal_dominating_sets())
    27

AUTHORS:

- Jean-Florent Raymond (2019-03-04): initial version


Methods:
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


def is_dominating(G, dom, focus=None):
    r"""
    Check whether ``dom`` is a dominating set of ``G``.

    We say that a set `D` of vertices of a graph `G` dominates a set `S` if
    every vertex of `S` either belongs to `D` or is adjacent to a vertex of `D`.
    Also, `D` is a dominating set of `G` if it dominates `V(G)`.

    INPUT:

    - ``dom`` -- iterable of vertices of ``G``; the vertices of the supposed
      dominating set.

    - ``focus`` -- iterable of vertices of ``G`` (default: ``None``); if
      specified, this method checks instead if ``dom`` dominates the vertices in
      ``focus``.

    EXAMPLES::

        sage: g = graphs.CycleGraph(5)
        sage: g.is_dominating([0,1], [4, 2])
        True

        sage: g.is_dominating([0,1])
        False
    """
    to_dom = set(G) if focus is None else set(focus)

    for v in dom:
        if not to_dom:
            return True
        to_dom.difference_update(G.neighbor_iterator(v, closed=True))

    return not to_dom

def is_redundant(G, dom, focus=None):
    r"""
    Check whether a ``dom`` has redundant vertices.

    For a graph `G` and sets `D` and `S` of vertices, we say that a vertex `v
    \in D` is *redundant* in `S` if `v` has no private neighbor with repect to
    `D` in `S`.  In other words, there is no vertex in `S` that is dominated by
    `v` but not by `D \setminus \{v\}`.

    INPUT:

    - ``dom`` -- iterable of vertices of ``G``; where we look for redundant
      vertices.

    - ``focus`` -- iterable of vertices of ``G`` (default: ``None``); if
      specified, this method checks instead whether ``dom`` has a redundant
      vertex in ``focus``.

    .. WARNING::

        The assumption is made that ``focus`` (if provided) does not contain
        repeated vertices.

    EXAMPLES::

        sage: G = graphs.CubeGraph(3)
        sage: G.is_redundant(['000', '101'], ['011'])
        True
        sage: G.is_redundant(['000', '101'])
        False
    """
    dom = list(dom)
    focus = list(G) if focus is None else list(focus)

    # dominator[v] (for v in focus) will be equal to:
    #  - (0, None) if v has no neighbor in dom
    #  - (1, u) if v has a unique neighbor in dom, u
    #  - (2, None) if v has >= 2 neighbors in dom

    # Initialization
    dominator = {v: (0, None) for v in focus}

    for x in dom:
        for v in G.neighbor_iterator(x, closed=True):
            if v in focus:
                # remember about x only if we never encountered
                # neighbors of v so far
                if dominator[v][0] == 0:
                    dominator[v] = (1, x)
                elif dominator[v][0] == 1:
                    dominator[v] = (2, None)

    # We now compute the subset of vertices of dom that have a private neighbor
    # in focus. A vertex v in dom has a private neighbor in focus if it is
    # dominated by a unique vertex, that is if dominator[v][0] == 1
    with_private = set()
    for v in focus:
        if dominator[v][0] == 1:
            with_private.add(dominator[v][1])

    # By construction with_private is a subset of dom and we assume the elements
    # of dom to be unique, so the following is equivalent to checking
    # with_private != set(dom)
    return len(with_private) != len(dom)

def private_neighbors(G, vertex, dom):
    r"""
    Return the private neighbors of a vertex with repect to other vertices.

    A private neighbor of a vertex `v` with respect to a vertex subset `D`
    is a closed neighbor of `v` that is not dominated by a vertex of `D
    \setminus \{v\}`.

    INPUT:

    - ``vertex`` -- a vertex of ``G``.

    - ``dom`` -- iterable of vertices of ``G``; the vertices possibly stealing
      private neighbors from ``vertex``.

    OUTPUT:

    Return the closed neighbors of ``vertex`` that are not closed neighbors
    of any other vertex of ``dom``.

    EXAMPLES::

        sage: g = graphs.PathGraph(5)
        sage: list(g.private_neighbors(1, [1, 3, 4]))
        [1, 0]

        sage: list(g.private_neighbors(1, [3, 4]))
        [1, 0]

        sage: list(g.private_neighbors(1, [3, 4, 0]))
        []
    """
    # The set of all vertices that are dominated by vertex_subset - vertex:
    closed_neighborhood_vs = set()
    for u in dom:
        if u != vertex:
            closed_neighborhood_vs.update(
                G.neighbor_iterator(u, closed=True))

    return (neighbor
            for neighbor in G.neighbor_iterator(vertex, closed=True)
            if neighbor not in closed_neighborhood_vs)



# ==============================================================================
# Enumeration of minimal dominating set as described in [BDHPR2019]_
# ==============================================================================

def _parent(G, dom, V_prev):
    r"""
    Return a subset of dom that is irredundant in ``V_prev``.

    For internal use.

    INPUT:

    - ``G`` -- a graph

    - ``dom`` -- an iterable of vertices of ``G``

    - ``V_prev`` -- an iterable of vertices of ``G``

    OUTPUT:

    Return the list obtained from ``dom`` by iteratively removing those
    vertices of mininum index that have no private neighbor in ``V_prev``.

    TESTS::

        sage: from sage.graphs.domination import _parent
        sage: G = graphs.PathGraph(4)
        sage: G.add_vertices([4, 5])
        sage: G.add_edges([(4, 1), (5, 2)])
        sage: _parent(G, [0, 2, 4, 5], [1, 2])
        [4, 5]
        sage: _parent(G, [0, 2, 4, 5], [1, 3])
        [2]
    """
    # The list where we search vertices
    D_start = sorted(dom, reverse=True)

    # The list to be output at the end, that we construct:
    D_end = []

    while D_start:
        v = D_start.pop()  # element of min index
        priv = set(G.neighbor_iterator(v, closed=True))
        # We remove the vertices already dominated
        # by other vertices of (D_end union D_start)
        priv.difference_update(*(G.neighbor_iterator(u, closed=True)
                                 for u in D_start if u != v))
        priv.difference_update(*(G.neighbor_iterator(u, closed=True)
                                 for u in D_end if u != v))
        # Now priv is the private neighborhood of v in G wrt D_start + D_end
        if priv.intersection(V_prev) != set():
            # if v has a private in V_prev, we keep it
            D_end.append(v)

    return D_end


def _peel(G, A):
    r"""
    Return a peeling of a vertex iterable of a graph.

    For internal use.
    Given a graph `G` and a subset `A` of its vertices, a peeling of `(G,A)` is
    a list `[(u_0, V_0), \dots, (u_p, V_p)]` such that `u_0` is ``None``,
    `V_0` is the empty set, `V_p = A` and for every `i \in \{1, \dots, p\}`,
    `V_{i-1} = V_i \setminus N[v_i]`, for some vertex `u_i` of `V_i`.

    INPUT:

    - ``G`` -- a graph

    - ``A`` -- a set of vertices of `G`

    OUTPUT:

    A peeling of `(G, A)`.

    TESTS::

        sage: from sage.graphs.domination import _peel
        sage: G = Graph(10); _peel(G, {0, 1, 2, 3, 4})
        [(None, set()),
        (4, {4}),
        (3, {3, 4}),
        (2, {2, 3, 4}),
        (1, {1, 2, 3, 4}),
        (0, {0, 1, 2, 3, 4})]


        sage: from sage.graphs.domination import _peel
        sage: G = graphs.PathGraph(10); _peel(G, set((i for i in range(10) if i%2==0)))
        [(None, set()),
        (8, {8}),
        (6, {6, 8}),
        (4, {4, 6, 8}),
        (2, {2, 4, 6, 8}),
        (0, {0, 2, 4, 6, 8})]

    """
    Acomp = set(G)
    Acomp.difference_update(A)  # Acomp  = V - A

    peeling = []
    H = copy(G)
    H.delete_vertices(list(Acomp))
    del Acomp

    while H:
        ui = next(H.vertex_iterator())  # pick some vertex of H
        Vi = set(H)
        peeling.append((ui, Vi))
        H.delete_vertices(H.neighbor_iterator(ui, closed=True))
    peeling.append((None, set()))
    peeling.reverse()
    return peeling


def _cand_ext_enum(G, to_dom, u_next):
    r"""
    Return the minimal dominating sets of ``to_dom``.

    For internal use.
    Assumption: ``u_next`` dominates ``to_dom``.

    INPUT:

    - ``G`` -- a graph

    - ``to_dom`` -- a ``set()`` of vertices of ``G``

    - ``u_next`` -- a vertex of ``G`` that dominates ``to_dom``

    OUTPUT:

    An iterator over the minimal dominating sets of ``to_dom``.

    TESTS::

        sage: from sage.graphs.domination import _cand_ext_enum
        sage: g = graphs.DiamondGraph()
        sage: l = list(_cand_ext_enum(g, {0, 2, 3}, 1,))
        sage: len(l) == 3 and {1} in l and {2} in l and {0, 3} in l
        True
    """

    def _aux_with_rep(H, to_dom, u_next):
        """
        Return the minimal dominating sets of ``to_dom``, with the
        assumption that ``u_next`` dominates ``to_som``.

        .. WARNING::

            The same output may be output several times (up to |H| times).

        In order to later remove duplicates, we here output pairs ``(ext, i)``
        where ``ext`` is the output candidate extension and ``i`` counts how
        many elements have already been output.
        """
        if u_next not in to_dom:
            # In this case, enumerating the minimal DSs of the subset
            # to_dom is a smaller instance as it excludes u_next:

            cand_ext_index = 0

            for ext in H.minimal_dominating_sets(to_dom):
                yield (ext, cand_ext_index)
                cand_ext_index += 1

        elif to_dom == {u_next}:
            # In this case, only u_next has to be dominated
            cand_ext_index = 0
            for w in H.neighbor_iterator(u_next, closed=True):
                # Notice that the case w = u_next is included
                yield ({w}, cand_ext_index)
                cand_ext_index += 1

        else:
            # In this case, both u_next and to_dom-u_next have to be dominated

            # We first output the trivial output
            # (as to_dom is subset of N(u_next)):
            yield ({u_next}, 0)
            # Start from 1 because we already output the 0-th elt:
            cand_ext_index = 1

            # When u_next is not in the DS, one of its neighbors w should be:
            for w in H.neighbor_iterator(u_next):

                remains_to_dom = set(to_dom)
                remains_to_dom.difference_update(
                    H.neighbor_iterator(w, closed=True))
                # Here again we recurse on a smaller instance at it
                # excludes u_next (and w)
                for Q in H.minimal_dominating_sets(remains_to_dom):
                    ext = set(Q)
                    ext.add(w)
                    # By construction w dominates u_next and Q dominates
                    # to_dom - N[w], so ext dominates to_dom: it is a
                    # valid output iff it is not redundant
                    if not H.is_redundant(ext):
                        yield (ext, cand_ext_index)
                        cand_ext_index += 1
    #
    # End of aux_with_rep routine

    # Here we use aux_with_rep twice to enumerate the minimal
    # dominating sets while avoiding repeated outputs
    for (X, i) in _aux_with_rep(G, to_dom, u_next):
        for (Y, j) in _aux_with_rep(G, to_dom, u_next):
            if j >= i:
                # This is the first time we meet X: we output it
                yield X
                break
            elif Y == X:  # These are sets
                # X has already been output in the past: we ignore it
                break


def minimal_dominating_sets(G, to_dominate=None, work_on_copy=True):
    r"""
    Return an iterator over the minimal dominating sets of a graph.

    INPUT:

    - ``G`` -- a graph.

    - ``to_dominate`` -- vertex iterable or ``None`` (default: ``None``);
      the set of vertices to be dominated.

    - ``work_on_copy`` -- boolean (default: ``True``); whether or not to work on
      a copy of the input graph; if set to ``False``, the input graph will be
      modified (relabeled).
    
    OUTPUT:

    An iterator over the inclusion-minimal sets of vertices of ``G``.
    If ``to_dominate`` is provided, return an iterator over the
    inclusion-minimal sets of vertices that dominate the vertices of
    ``to_dominate``.

    ALGORITHM: The algorithm described in [BDHPR2019]_.

    EXAMPLES::

        sage: G = graphs.ButterflyGraph()
        sage: ll = list(G.minimal_dominating_sets())
        sage: pp = [{0, 1}, {1, 3}, {0, 2}, {2, 3}, {4}]
        sage: len(pp) == len(pp) and all(x in pp for x in ll) and all(x in ll for x in pp)
        True

        sage: ll = list(G.minimal_dominating_sets([0,3]))
        sage: pp = [{0}, {3}, {4}]
        sage: len(pp) == len(pp) and all(x in pp for x in ll) and all(x in ll for x in pp)
        True

        sage: ll = list(G.minimal_dominating_sets([4]))
        sage: pp = [{4}, {0}, {1}, {2}, {3}]
        sage: len(pp) == len(pp) and all(x in pp for x in ll) and all(x in ll for x in pp)
        True

    ::

        sage: ll = list(graphs.PetersenGraph().minimal_dominating_sets())
        sage: pp = [{0, 2, 6},
        ....: {0, 9, 3},
        ....: {0, 8, 7},
        ....: {1, 3, 7},
        ....: {1, 4, 5},
        ....: {8, 1, 9},
        ....: {8, 2, 4},
        ....: {9, 2, 5},
        ....: {3, 5, 6},
        ....: {4, 6, 7},
        ....: {0, 8, 2, 9},
        ....: {0, 3, 6, 7},
        ....: {1, 3, 5, 9},
        ....: {8, 1, 4, 7},
        ....: {2, 4, 5, 6},
        ....: {0, 1, 2, 3, 4},
        ....: {0, 1, 2, 5, 7},
        ....: {0, 1, 4, 6, 9},
        ....: {0, 1, 5, 6, 8},
        ....: {0, 8, 3, 4, 5},
        ....: {0, 9, 4, 5, 7},
        ....: {8, 1, 2, 3, 6},
        ....: {1, 2, 9, 6, 7},
        ....: {9, 2, 3, 4, 7},
        ....: {8, 2, 3, 5, 7},
        ....: {8, 9, 3, 4, 6},
        ....: {8, 9, 5, 6, 7}]
        sage: len(pp) == len(pp) and all(x in pp for x in ll) and all(x in ll for x in pp)
        True

    TESTS:

    The empty graph is handled correctly::

        sage: list(Graph().minimal_dominating_sets())
        [set()]

    Test on all graphs on 6 vertices::

        sage: from sage.combinat.subset import Subsets
        sage: def minimal_dominating_sets_naive(G):
        ....:     return (S for S in Subsets(G.vertices())
        ....:             if not(G.is_redundant(S)) and G.is_dominating(S))
        sage: def big_check(n):
        ....:     for G in graphs(n):
        ....:         ll = list(G.minimal_dominating_sets())
        ....:         pp = list(minimal_dominating_sets_naive(G))
        ....:         if (len(pp) != len(pp)
        ....:             or any(x not in pp for x in ll)
        ....:             or any(x not in ll for x in pp)):
        ....:             return False
        ....:     return True
        sage: big_check(6)  # long time
        True

    Outputs are unique::

        sage: def check_uniqueness(g):
        ....:     counter_1 = 0
        ....:     for dom_1 in g.minimal_dominating_sets():
        ....:         counter_1 += 1
        ....:         counter_2 = 0
        ....:         for dom_2 in g.minimal_dominating_sets():
        ....:             counter_2 += 1
        ....:             if counter_2 >= counter_1:
        ....:                 break
        ....:             if dom_1 == dom_2:
        ....:                 return False
        ....:     return True
        sage: check_uniqueness(graphs.RandomGNP(9, 0.5))
        True
    """

    def tree_search(H, plng, dom, i):
        r"""
        Enumerate minimal dominating sets recursively.

        INPUT:

        - ``H`` -- a graph

        - ``plng`` -- a peeling of H (result of :func:`_peel`)

        - ``dom`` -- a minimal dominating set of ``plng[i][1]``

        - ``i`` -- an integer, the current position in ``plng``

        OUTPUT:

        An iterator over those minimal dominating sets (in ``H``) of
        ``plng[-1][1]`` that are children of ``dom`` (with respect to the
        :func:`parent` function).

        ALGORITHM:

        We iterate over those minimal dominating sets of ``plng[i + 1][1]`` that
        are children of dom and call recursively on each. The fact that we
        iterate over children (with respect to the `parent` function) ensures
        that we do not have repeated outputs.
        """
        if i == len(plng) - 1:
            # we reached a leaf, i.e. dom is a minimal dominating set
            # of plng[i][1] = plng[-1][1]
            yield dom
            return

        u_next, V_next = plng[i + 1]

        if H.is_dominating(dom, V_next):
            # if dom dominates V_next
            # then dom is its unique extension: we recurse on it
            for Di in tree_search(H, plng, dom, i + 1):
                yield Di
            return

        # Otherwise, V_next - <what dom dominates> is what we have to dominate
        to_dom = V_next - set().union(
            *(G.neighbor_iterator(vert, closed=True)
              for vert in dom))

        for can_ext in _cand_ext_enum(H, to_dom, u_next):

            # We complete dom with can_ext -> canD
            canD = set().union(can_ext, dom)

            if (not H.is_redundant(canD, V_next)) and set(dom) == set(_parent(H, canD, plng[i][1])):
                # By construction, can_ext is a dominating set of
                # `V_next - N[dom]`, so canD dominates V_next.
                # If canD is a legitimate child of dom and is not redundant, we
                # recurse on it:
                for Di in tree_search(H, plng, canD, i + 1):
                    yield Di
    ##
    # end of tree-search routine

    int_to_vertex = list(G)
    vertex_to_int = {u: i for i, u in enumerate(int_to_vertex)}

    if work_on_copy:
        G.relabel(perm=vertex_to_int)
    else:
        G = G.relabel(perm=vertex_to_int, inplace=False)

    if to_dominate is None:
        vertices_to_dominate = set(G)
    else:
        vertices_to_dominate = set(to_dominate)

    if not vertices_to_dominate:
        # base case: vertices_to_dominate is empty
        # the empty set/list is the only minimal DS of the empty set
        yield set()
        return

    peeling = _peel(G, vertices_to_dominate)

    for dom in tree_search(G, peeling, set(), 0):
        yield {int_to_vertex[v] for v in dom}
