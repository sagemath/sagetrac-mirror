"""
Given a list of edges `edges` defining a simple graph `G`, `edges[:i]`
defines a subgraph `G_i`. Define as `active nodes` of `G_i`
the nodes of `G_i` which have degree less than in `G`.
Define as number of active nodes for given `edges` the maximum among the 
number of active nodes of all `G_i`; the active node number of `G`
is the minimum number of active nodes for all possible orderings
of the edges of `G`.
To compute efficiently the matching polynomial using the
algorithm by Butera and Pernici, one must find an ordering of edges
with low number of active nodes.
Here a greedy algorithm is used to compute a generally not optimal
ordering of edges to compute the matching polynomial.
A similar algorithm is used to find an order of vertices to
compute the independence polynomial.

"""

from collections import defaultdict, deque

def _add_links_no_incr(links0, d1, d2):
    """
    Return links added without increasing the number of active nodes

    links0 : links forming `d2`
    d1 : adjacency dict of the graph
    d2 : adjacency dict of the part of the graph covered

    The active nodes are the nodes in `d2` with degree less than in `d1`.
    The added links are appended to links0

    Examples
    ========

    >>> from collections import defaultdict
    >>> from active_nodes import _add_links_no_incr
    >>> links0=[(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)]
    >>> d1 = {0:[1,4,5], 1:[0,2,6], 2:[1,3,7], 3:[2,4,8], 4:[0,3,9], 5:[0,7,8], 6:[1,8,9], 7:[2,5,9], 8:[3,5,6], 9:[4,6,7]}
    >>> d2 = defaultdict(list, {0:[1,4], 1:[0,2], 2:[1,3], 3:[2,4], 4:[3,0]})
    [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9), (5, 7), (5, 8), (6, 8), (6, 9), (7, 9)]
    """
    links = []
    hit = True
    while hit:
        hit = False
        for k, a in d2.items():
            # if there is only one edge missing to complete a node, add it
            if len(a) == len(d1[k]) - 1:
                a1 = d1[k]
                for x in a1: 
                    if x not in a:
                        k1 = x
                        break
                if k > k1:
                    k, k1 = k1, k
                d2[k].append(k1)
                d2[k1].append(k)
                hit = True
                links.append((k, k1))
        # it two nodes adjacent in d1 are present in d2, and there
        # is no edge between them in d2, add it
        for k1, a in d1.items():
            if not d2[k1]:
                continue
            for k2 in a:
                if d2[k2] and k1 not in d2[k2]:
                    d2[k1].append(k2)
                    d2[k2].append(k1)
                    hit = True
                    if k1 < k2:
                        links.append((k1, k2))
                    else:
                        links.append((k2, k1))

    links0.extend(links)
    return links

def _append_link(dx, links, k1, k2):
    """
    add the edge `(k1, k2)` to links and update the dict `dx`
    """
    assert not (k1,k2) in links
    assert not (k2,k1) in links
    dx[k1].append(k2)
    dx[k2].append(k1)
    if k1 < k2:
        links.append((k1, k2))
    else:
        links.append((k2, k1))

def _short_path_active_nodes(adj_list, s, act):
    """
    Return a short path starting in `s` and ending in a node in `act`

    adj_list : adjacency list of the graph
    s : an active node
    act : list of active nodes

    Perform a breadth-first search.
    Return None if there is not such a path.

    """
    prev = {s: None}
    q = deque([s])
    while q:
        u = q.popleft()
        for v in adj_list[u]:
            if v in prev:
                continue
            prev[v] = u
            if v in act:
                # found a path from s to another active node
                a1 = [v]
                while 1:
                    v = prev[v]
                    a1.append(v)
                    if v in act:
                        return a1
                assert 0
            q.append(v)
    return None

def _short_path_active_nodes_all(d, dx, active):
    """
    return a short path between two active nodes

    d : dictionary of the graph G
    Gx : dict of a subgraph of G
    active : active nodes

    The returned path has the endpoints in `active` and the rest of the
    nodes in `G - Gx`

    Examples
    ========

    >>> from collections import defaultdict
    >>> from active_nodes import _short_path_active_nodes_all
    >>> d={0: [1, 2, 5], 1: [0, 3, 5], 2: [0, 3, 4, 5], 3: [1, 2], 4: [2, 5], 5: [0, 1, 2, 4]}
    >>> dx=defaultdict(list, {0: [1], 1: [0]})
    >>> active=[0, 1]
    >>> _short_path_active_nodes_all(d, dx, active)
    [1, 5, 0]
    """
    ax = {}
    # adjacency list of complement of G and Gx
    for k, v in d.items():
        v1 = dx[k]
        ak = [y for y in v if y not in v1]
        ax[k] = ak
    min_length = 10000
    min_path = None
    for i in active:
        a1 = _short_path_active_nodes(ax, i, active)
        if a1 is not None and len(a1) < min_length:
            min_length = len(a1)
            min_path = a1
    return min_path

def _add_paths(d, dx, links, path):
    """
    Return the list of added edges obtained adding `path`

    d : dict of the graph G
    dx : dict of the subgraph Gx (updated)
    links : list of edges forming Gx (updated)
    path : path whose edges must be added

    Examples
    ========

    >>> from active_nodes import _add_paths
    >>> d={0: [1,2,3,5], 1: [0,3,4,5], 2: [0,4], 3: [0,1], 4: [1,2], 5: [0,1]}
    >>> dx = defaultdict(list,{0:[1,3], 1:[0,3], 3:[1,0]})
    >>> links = [(0,1)]
    >>> path = [1,3,0]
    >>> _add_paths(d, dx, links, path)
    [(1, 5), (0, 5), (0, 2), (1, 4), (2, 4)]
    """
    added = []
    for i in range(len(path) - 1):
        k1, k2 = path[i], path[i+1]
        if k2 < k1:
            k2, k1 = k1, k2
        _append_link(dx, links, k1, k2)
        added.append((k1, k2))
    r = _add_links_no_incr(links, d, dx)
    added.extend(r)
    return added

def _d_relabel(d):
    """
    relabel `keys` to be in `range(len(keys))`
    return the relabeled graph and the mapping from old to new keys

    Examples
    ========

    >>> from hobj import _d_relabel
    >>> _d_relabel({1:[2,3], 2:[1,3], 3:[1,2]})
    ({0: [1, 2], 1: [0, 2], 2: [0, 1]}, {1: 0, 2: 1, 3: 2})
    """
    dt = {}
    keys = list(d.keys())
    dt = dict(zip(keys, range(len(d))))
    d1 = {}
    for k, v in d.items():
        d1[dt[k]] = [dt[i] for i in v]
    return d1, dt

def ordered_links(d, k0, k1):
    """
    Return ordered links starting from the link (k0, k1)

    Parameters
    ==========

    d : dict for the graph
    k0, k1: two adjacent nodes of the graphs

    Notes
    =====

    The vertices in `d` must be numbered in `0,..., len(d)-1`

    The algorithm attempts to order the edges in such a way that
    the number of active nodes is kept small.

    Examples
    ========

    >>> from active_nodes import ordered_links
    >>> d = {0:[1,4], 1:[0,2], 2:[1,3], 3:[2,4], 4:[0,3]}
    >>> ordered_links(d, 0, 1)
    [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]
    """
    d0 = d
    assert k0 in d
    assert k1 in d[k0]
    links_all = []
    # dt relabels indices in range(len(d))
    while 1:
        # di maps to the original labelling
        dx = defaultdict(list)
        links = []
        _append_link(dx, links, k0, k1)
        if len(d[k0]) == 1:
            d1 = {}
            for k, v in d.items():
                d1[k] = v[:]
            d = d1
            del d[k0]
            d[k1].remove(k0)
            if k0 < k1:
                links_all.append((k0, k1))
            else:
                links_all.append((k1, k0))
            if len(d[k1]):
                k0 = k1
                k1 = d[k0][0]
            else:
                del d[k1]
                if not d:
                    break
                for kx in d.keys():
                    if len(d[kx]) > 0:
                        k0 = kx
                        k1 = d[k0][0]
                        break
            continue

        r = _add_links_no_incr(links, d, dx)
        while 1:
            active = [k for k in dx if 0 < len(dx[k]) < len(d[k])]
            if not active:
                break
            short_path = _short_path_active_nodes_all(d, dx, active)
            if short_path is None:
                break
            # add the edges in short_path to links
            for i in range(len(short_path) - 1):
                k1x, k2x = short_path[i], short_path[i+1]
                if k2x < k1x:
                    k2x, k1x = k1x, k2x
                _append_link(dx, links, k1x, k2x)
            _add_links_no_incr(links, d, dx)

        for i, j in links:
            if i < j:
                links_all.append((i, j))
            else:
                links_all.append((j, i))
        # edges in `d`
        edges = []
        for k1 in d:
            for k2 in d[k1]:
                if k1 < k2:
                    edges.append((k1,k2))
        # edges1 edges in `d` but not in links
        edges1 = []
        for edge in edges:
            if edge not in links:
                edges1.append(edge)
        if not edges1:
            break
        # initialize the next loop with new d, k0, k1 made from edge1
        d = defaultdict(list)
        for i, j in edges1:
            d[i].append(j)
            d[j].append(i)
        k0 = d.keys()[0]
        k1 = d[k0][0]

    return links_all


def num_active_nodes(d, links):
    """
    number of the active nodes for the graph defined by ``links``

    Parameters
    ==========

    d : dict for the graph
    links : list of edges of the graph

    Examples
    ========

    >>> from active_nodes import ordered_links
    >>> from active_nodes import ordered_links, num_active_nodes
    >>> d = {0:[1,4], 1:[0,2], 2:[1,3], 3:[2,4], 4:[0,3]}
    >>> links = ordered_links(d, 0, 1)
    >>> num_active_nodes(d, links)
    2
    """
    dx = defaultdict(list)
    max_active = 0
    for edge in links:
        k1, k2 = edge
        dx[k1].append(k2)
        dx[k2].append(k1)
        nactivex = len([k for k in dx if 0 < len(dx[k]) < len(d[k])])
        if max_active < nactivex:
            max_active = nactivex
    return max_active

