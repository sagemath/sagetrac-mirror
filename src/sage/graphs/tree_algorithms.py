def tree_diameter(G, endpoints=false, path=false):
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

def tree_centre(G):
    (dist, dpath) = tree_diameter(G, path = 'true')
    if dist%2 == 1:
        return (dpath[dist//2], dpath[dist//2+1])
    else:
        return (dpath[dist//2],)

def dfs(u, G):
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

def rooted_tree_isomorphism(G, H, r1, r2):
    (h1, L1, par1, children1) = dfs(r1, G)
    (h2, L2, par2, children2) = dfs(r2, H)
    if h1 != h2:
        return false
    h = h1 = h2
    label1 = dict()
    subLabel1 = dict()
    label2 = dict()
    subLabel2 = dict()
    for i in range(0,h):
        if len(L1[i]) != len(L2[i]):
            return false
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
        L1[i] = sorted(L1[i], key=lambda x: subLabel1[x])
        L2[i] = sorted(L2[i], key=lambda x: subLabel2[x])
        cnt = 0
        for j in range(len(L1[i])):
            if subLabel1[L1[i][j]] != subLabel2[L2[i][j]]:
                return false
        c = 0
        label1[L1[i][0]] = label2[L2[i][0]] = 0
        for j in range(1, len(L1[i])):
            if subLabel1[L1[i][j]] != subLabel1[L1[i][j-1]]:
                c = c+1
            label1[L1[i][j]] = c
            label2[L2[i][j]] = c
    if subLabel1[r1] != subLabel2[r2]:
        return false
    return true

def tree_isomorphism(G, H):
    gc = tree_centre(G)
    hc = tree_centre(H)
    if len(gc) != len(hc):
        return false
    for p in gc:
        for q in hc:
            f = rooted_tree_isomorphism(G, H, p, q)
            if f:
                return f
    return false
