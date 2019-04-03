def tree_diameter(G, endpoints=False, path=False):
    """
    Returns diameter of unweighted tree.
    Input:

    - ``G`` -- graph

    - ``endpoints`` -- boolean(default: ``False``); whether to return
    enpoints of diameter or not

    - ``path`` -- boolean(default: ``False``); whether to return the
    diameter path

    EXAMPLES::

        sage: from sage.graphs.tree_algorithms import tree_diameter
        sage: G = Graph([(1,2),(2,3),(3,4)])
        sage: tree_diameter(G)
        3

        sage: G = Graph([('a','b'),('b','c'),('c','d')])
        sage: tree_diameter(G)
        3

        sage: G = Graph([(1,2),(2,3),(2,4)])
        sage: tree_diameter(G)
        2

        sage: G = Graph([(1,2),(2,3),(2,4)])
        sage: tree_diameter(G, endpoints=True, path=True)
        (2, (3, 1), [3, 2, 1])

        sage: G = Graph([('a','b'),('b','c'),('c','d')])
        sage: tree_diameter(G, endpoints=True, path=True)
        (3, ('d', 'a'), ['d', 'c', 'b', 'a'])

    """
    u = next(G.vertex_iterator())
    path_lengths = G.shortest_path_lengths(u, algorithm='BFS')
    from operator import itemgetter
    u = max(path_lengths.iteritems(), key=itemgetter(1))[0]
    path_lengths = G.shortest_path_lengths(u, algorithm='BFS')
    (v, dist) = max(path_lengths.iteritems(), key=itemgetter(1))
    if not endpoints and not path:
        return dist
    if endpoints and not path:
        return (dist, (u, v))
    dpath = G.shortest_path(u,v,algorithm='BFS')
    if not endpoints and path:
        return (dist, dpath)
    return (dist, (u, v), dpath)

def tree_center(G):
    """
    Return center of the given tree(G).

    Input:

    - ``G`` -- graph

    EXAMPLES::

        sage: from sage.graphs.tree_algorithms import tree_center
        sage: G = Graph([('a','b'),('b','c'),('c','d')])
        sage: tree_center(G)
        ('c', 'b')

        sage: G = Graph([(1,2),(2,3),(2,4)])
        sage: tree_center(G)
        (2,)

    """
    (dist, dpath) = tree_diameter(G, path = True)
    if dist%2 == 1:
        return (dpath[dist//2], dpath[dist//2+1])
    else:
        return (dpath[dist//2],)

def generateCertificate(g, h, children1, children2, subLabel1, subLabel2):
    """
    Returns the certificate for isomorphic rooted trees.

    Input:

    - ``g`` -- vertex; root of the first tree

    - ``h`` -- vertex; root of the second tree

    - ``children1`` -- dictionary; stores the children of vertices in
    first tree

    - ``children2`` -- dictionary; stores the children of vertices in
    second tree

    -- ``subLabel1`` -- dictionary; stores the label assigned to the
    vertex of first tree

    -- ``subLabel2`` -- dictionary; stores the label assigned to the
    vertex of second tree

    """
    st1 = []
    st2 = []
    cert = dict()
    st1.append(g)
    st2.append(h)
    while st1 and st2:
        x = st1[-1]
        y = st2[-1]
        cert[x] = y
        st1.pop()
        st2.pop()
        children1[x].sort(key=lambda z: subLabel1[z])
        children2[y].sort(key=lambda z: subLabel2[z])
        for i in range(len(children1[x])):
            st1.append(children1[x][i])
            st2.append(children2[y][i])
    return cert

def tree_dfs(u, G):
    """
    Returns a tuple containing number of levels, vertices in each level,
    parent of each vertex and children of each vertex after doing dfs on
    the given tree ``G``.

    Input:

    - ``u`` -- vertex; vertex from which dfs is performed(or root of tree)

    - ``G`` -- graph

    Output:

    - ``level`` -- integer; number of levels in tree

    - ``L`` -- list of list; contains list of vertices contained in each
    level

    - ``par`` -- dictionary; contains parent of each vertex in dfs
    traversal

    - ``children`` -- dictionary; contains list of children corresponding
    to each vertex.

    """
    neighbours = G.neighbor_iterator
    stack = [u]
    diameter = tree_diameter(G)
    L = [[] for i in range(diameter+1)]
    L[0].append(u)
    par = dict()
    par[u] = -1
    children = dict()
    for ver in G.vertex_iterator():
        children[ver] = []
    dis = dict()
    dis[u] = 0
    level = 0
    while stack:
        v = stack[-1]
        stack.pop()
        for x in neighbours(v):
            if par[v] != x:
                stack.append(x)
                children[v].append(x)
                par[x] = v
                dis[x] = dis[v]+1
                level = max(level,dis[x])
                L[dis[x]].append(x)
    return (level, L, par, children)

def rooted_tree_isomorphism(G, H, gc, hc, certificate=False):
    """
    Returns if rooted tree G(root gc) is isomorphic to rooted
    tree(root hc).

    Input:

    - ``G`` -- graph

    - ``H`` -- graph

    - ``gc`` -- vertex; root of G

    - ``hc`` -- vertex; root of H

    """
    (h1, L1, par1, children1) = tree_dfs(gc, G)
    (h2, L2, par2, children2) = tree_dfs(hc, H)
    if h1 != h2:
        return (False, None)
    h = h1 = h2
    label1 = dict()
    subLabel1 = dict()
    label2 = dict()
    subLabel2 = dict()
    for i in range(0,h):
        if len(L1[i]) != len(L2[i]):
            return (False, None)
    for x in G.vertex_iterator():
        subLabel1[x] = []
    for x in H.vertex_iterator():
        subLabel2[x] = []
    for i in range(len(L1[h])):
        label1[L1[h][i]] = label2[L2[h][i]] = 0
    for i in range(h-1,0,-1):
        for y in L1[i+1]:
            subLabel1[par1[y]].append(label1[y])
        for y in L2[i+1]:
            subLabel2[par2[y]].append(label2[y])
        L1[i].sort(key=lambda x: subLabel1[x])
        L2[i].sort(key=lambda x: subLabel2[x])
        for j in range(len(L1[i])):
            if subLabel1[L1[i][j]] != subLabel2[L2[i][j]]:
                return (False, None)
        c = 0
        label1[L1[i][0]] = label2[L2[i][0]] = 0
        for j in range(1, len(L1[i])):
            if subLabel1[L1[i][j]] != subLabel1[L1[i][j-1]]:
                c = c+1
            label1[L1[i][j]] = c
            label2[L2[i][j]] = c
    if subLabel1[gc] != subLabel2[hc]:
        return (False, None)
    if certificate:
        return (True, generateCertificate(gc, hc, children1, children2, subLabel1, subLabel2))
    return (True, None)

def tree_isomorphism(G, H, certificate=False):
    """
    Returns if trees ``G`` and ``H`` are isomorphic with certificate in
    O(n log n)

    Input:

    - ``G`` -- graph

    - ``H`` -- graph

    - ``certificate`` -- boolean(default: False); True if certificate has
    to be returned else False

    EXAMPLES::

        sage: from sage.graphs.tree_algorithms import tree_isomorphism
        sage: G = Graph([(1,2),(2,3),(2,4)])
        sage: H = Graph([(2,3),(3,4),(3,1)])
        sage: tree_isomorphism(G,H,certificate=True)
        (True, {1: 1, 2: 3, 3: 2, 4: 4})

        sage: G = Graph([(1,2),(2,3),(2,4)])
        sage: H = Graph([(2,3),(3,4),(4,1)])
        sage: tree_isomorphism(G,H,certificate=True)
        (False, None)

    """
    gc = tree_center(G)
    hc = tree_center(H)
    if len(gc) != len(hc):
        return (False, None) if certificate else False
    for p in gc:
        for q in hc:
            f = rooted_tree_isomorphism(G, H, p, q, certificate)
            if f[0]:
                return f if certificate else True
    return (False, None) if certificate else False
