def tree_diameter(G, endpoints=False, path=False):
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
    (dist, dpath) = tree_diameter(G, path = True)
    if dist%2 == 1:
        return (dpath[dist//2], dpath[dist//2+1])
    else:
        return (dpath[dist//2],)

def generateCertificate(g, h, children1, children2, subLabel1, subLabel2):
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
        cnt = 0
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

def tree_isomorphism(G, H, certificate):
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
