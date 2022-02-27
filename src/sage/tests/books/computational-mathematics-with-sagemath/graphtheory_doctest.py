## -*- encoding: utf-8 -*-
"""
This file (./graphtheory_doctest.sage) was *autogenerated* from ./graphtheory.tex,
with sagetex.sty version 2011/05/27 v2.3.1.
It contains the contents of all the sageexample environments from this file.
You should be able to doctest this file with:
sage -t ./graphtheory_doctest.sage
It is always safe to delete this file; it is not used in typesetting your
document.

Sage example in ./graphtheory.tex, line 12::

  sage: set_random_seed(1588)
  sage: sage.numerical.backends.generic_backend.default_solver = "Glpk"

Sage example in ./graphtheory.tex, line 59::

  sage: g = Graph()

Sage example in ./graphtheory.tex, line 87::

  sage: g.order(), g.size()
  (0, 0)
  sage: g.add_vertex(0)
  sage: g.order(), g.size()
  (1, 0)
  sage: g.add_vertices([1, 2, 5, 9])
  sage: g.order(), g.size()
  (5, 0)
  sage: g.add_edges([(1,5), (9,2), (2,5), (1,9)])
  sage: g.order(), g.size()
  (5, 4)
  sage: g.add_edge("Madrid", "Edinburgh")
  sage: g.order(), g.size()
  (7, 5)

Sage example in ./graphtheory.tex, line 128::

  sage: g.delete_vertex(0)
  sage: g.delete_edges([(1,5), (2,5)])
  sage: g.order(), g.size()
  (6, 3)
  sage: sorted(g.vertices(sort=False), key=str)
  [1, 2, 5, 9, 'Edinburgh', 'Madrid']
  sage: sorted(g.edges(sort=False, labels=False), key=str)
  [('Edinburgh', 'Madrid'), (1, 9), (2, 9)]

Sage example in ./graphtheory.tex, line 159::

  sage: g = Graph({
  ....:      0: [],
  ....:      1: [5, 9],
  ....:      2: [1, 5, 9],
  ....:      'Edinburgh': ['Madrid']})

Sage example in ./graphtheory.tex, line 212::

  sage: P = graphs.PetersenGraph()
  sage: C = graphs.ChvatalGraph()

Sage example in ./graphtheory.tex, line 252::

  sage: P = graphs.PetersenGraph()
  sage: P.is_planar()
  False
  sage: P.minor(graphs.CompleteBipartiteGraph(3,3))  # random
  {0: [1], 1: [8], 2: [4], 3: [6, 7, 9], 4: [2, 3], 5: [0, 5]}
  sage: P.minor(graphs.CompleteGraph(5))  # random
  {0: [1, 6], 1: [0, 5], 2: [2, 7], 3: [4, 9], 4: [3, 8]}
  sage: P.girth()
  5
  sage: P.is_regular(3)
  True
  sage: P.chromatic_number()
  3
  sage: P.is_vertex_transitive()
  True
  sage: P.show()

Sage example in ./graphtheory.tex, line 332::

  sage: K = graphs.KneserGraph(5, 2); P = graphs.PetersenGraph()
  sage: K.is_isomorphic(P)
  True

Sage example in ./graphtheory.tex, line 340::

  sage: all( graphs.KneserGraph(n,k).chromatic_number() == n - 2*k + 2
  ....:      for n in range(5, 9) for k in range(2, n // 2) )
  True

Sage example in ./graphtheory.tex, line 459::

  sage: H = graphs.ClawGraph()
  sage: def test():
  ....:    g = graphs.RandomGNP(20,2/5)
  ....:    return not g.subgraph_search(H, induced=True) is None
  sage: sum( test() for i in range(100) ) >= 80
  True

Sage example in ./graphtheory.tex, line 485::

  sage: P = graphs.PetersenGraph()
  sage: H = graphs.HoffmanSingletonGraph()
  sage: U = P + H; U2 = P.disjoint_union(H)
  sage: U.is_isomorphic(U2)
  True

Sage example in ./graphtheory.tex, line 496::

  sage: C = graphs.ChvatalGraph()
  sage: U = 3 * C; U2 = C.disjoint_union(C.disjoint_union(C))
  sage: U2.is_isomorphic(U)
  True

Sage example in ./graphtheory.tex, line 506::

  sage: U = 3*P + 2*C

Sage example in ./graphtheory.tex, line 513::

  sage: all( (CC.is_isomorphic(P) or CC.is_isomorphic(C))
  ....:       for CC in U.connected_components_subgraphs() )
  True

Sage example in ./graphtheory.tex, line 522::

  sage: sum( CC.is_isomorphic(P)
  ....:      for CC in U.connected_components_subgraphs() )
  3
  sage: sum( CC.is_isomorphic(C)
  ....:      for CC in U.connected_components_subgraphs() )
  2

Sage example in ./graphtheory.tex, line 568::

  sage: C = graphs.ChvatalGraph(); C.show()

Sage example in ./graphtheory.tex, line 583::

  sage: C.show(partition = [C.independent_set()])

Sage example in ./graphtheory.tex, line 597::

  sage: C.show(vertex_colors = {
  ....:   "red" : [0, 1, 2],    "blue" : [3, 4, 5],
  ....:   "yellow" : [6, 7, 8], "purple" : [9, 10, 11]})

Sage example in ./graphtheory.tex, line 608::

  sage: C.coloring(hex_colors = True)
  {'#00ffff': [3, 8, 5],
     '#7f00ff': [11],
     '#7fff00': [1, 4, 6, 9],
     '#ff0000': [0, 2, 7, 10]}
  sage: C.show(vertex_colors = C.coloring(hex_colors = True))

Sage example in ./graphtheory.tex, line 644::

  sage: L = [graphs.CompleteGraph(i) for i in range(3,3+10)]
  sage: for number, G in enumerate(L):
  ....:     G.plot().save(tmp_filename(ext=".png"))

Sage example in ./graphtheory.tex, line 782::

  sage: P5 = graphs.PathGraph(5); House = graphs.HouseGraph()
  sage: P5.complement().is_isomorphic(House)
  True
  sage: P4 = graphs.PathGraph(4); P4.complement().is_isomorphic(P4)
  True
  sage: C5 = graphs.CycleGraph(5); C5.complement().is_isomorphic(C5)
  True

Sage example in ./graphtheory.tex, line 852::

  sage: n = 5; Path = graphs.PathGraph(n)
  sage: Grid = Path.cartesian_product(Path)
  sage: Grid.is_isomorphic(graphs.GridGraph([n,n]))
  True

Sage example in ./graphtheory.tex, line 1144::

  sage: n = 30; p = 0.3; trials = 50
  sage: def equality(G):
  ....:     return G.edge_connectivity() == min(G.degree())
  sage: sum(equality(graphs.RandomGNP(n,p)) for i in range(trials))/trials
  1

Sage example in ./graphtheory.tex, line 1319::

  sage: g = graphs.ChvatalGraph(); cycle = g.hamiltonian_cycle()
  sage: g.show(vertex_labels = False); cycle.show(vertex_labels = False)

Sage example in ./graphtheory.tex, line 1586::

  sage: set_random_seed(3291)

Sage example in ./graphtheory.tex, line 1589::

  sage: n = 100; p = 5/n; g = graphs.RandomGNP(n, p)

Sage example in ./graphtheory.tex, line 1597::

  sage: # Set of available colors.
  sage: # In the worst-case scenario up to n colors suffice
  sage: available_colors = Set(range(n))

Sage example in ./graphtheory.tex, line 1611::

  sage: # This dictionary contains the color associated
  sage: # with each vertex of the graph
  sage: color = {}
  sage: for u in g:
  ....:    forbidden = Set([color[v] for v in g.neighbors(u)
  ....:                              if v in color])
  ....:    color[u] = min(available_colors - forbidden)

Sage example in ./graphtheory.tex, line 1625::

  sage: # Number of colors used
  sage: max(color.values()) + 1
  6

Sage example in ./graphtheory.tex, line 1635::

  sage: P = Permutations([0,1,2,3]); P.random_element()
  [2, 0, 1, 3]

Sage example in ./graphtheory.tex, line 1646::

  sage: available_colors = Set(range(n))

Sage example in ./graphtheory.tex, line 1656::

  sage: n_tests = 30
  sage: vertices = g.vertices()
  sage: P = Permutations(range(n))
  sage: best_coloring = {}
  sage: best_chromatic_number = +oo

Sage example in ./graphtheory.tex, line 1678::

  sage: for t in range(n_tests):
  ....:    # Random ordering of vertices
  ....:    p = P.random_element()
  ....:    color = {}
  ....:    for i in range(g.order()):
  ....:          u = vertices[p[i]]
  ....:          forbidden = Set([color[v] for v in g.neighbors(u)
  ....:                      if v in color])
  ....:          color[u] = min(available_colors - forbidden)
  ....:    # Update the best coloring
  ....:    if max(color.values()) + 1 < best_chromatic_number:
  ....:          best_coloring = color
  ....:          best_chromatic_number = 1 + max(color.values())

Sage example in ./graphtheory.tex, line 1697::

  sage: best_chromatic_number # Number of colors used
  4

Sage example in ./graphtheory.tex, line 1718::

  sage: def greedy_coloring(g, permutation):
  ....:     n = g.order()
  ....:     available_colors = Set(range(n))
  ....:     vertices = g.vertices()
  ....:     color = {}
  ....:     for i in range(n):
  ....:         u = vertices[permutation[i]]
  ....:         forbidden = Set([color[v] for v in g.neighbors(u)
  ....:                         if v in color])
  ....:         color[u] = min(available_colors - forbidden)
  ....:     return max(color.values()) + 1, color

Sage example in ./graphtheory.tex, line 1736::

  sage: set_random_seed(0)

Sage example in ./graphtheory.tex, line 1746::

  sage: P = Permutations(range(g.order()))
  sage: n_colors, coloration = min([greedy_coloring(g, P.random_element())
  ....:     for i in range(50)], key=lambda c: c[0])
  sage: n_colors
  4

Sage example in ./graphtheory.tex, line 1782::

  sage: n = 20; k = 4; g = graphs.RandomGNP(n, 0.5)
  sage: g = g.subgraph(edges = g.min_spanning_tree())

Sage example in ./graphtheory.tex, line 1786::

  sage: while True:
  ....:    _, edges, [S,Sb] = g.edge_connectivity(vertices = True)
  ....:    cardinality = len(edges)
  ....:    if cardinality < k:
  ....:         CP = cartesian_product([S, Sb])
  ....:         g.add_edges([CP.random_element()
  ....:                      for i in range(k - len(edges))])
  ....:    else:
  ....:         break

Sage example in ./graphtheory.tex, line 1844::

  sage: g = graphs.RandomGNP(40, 0.4)
  sage: P = Permutations(range(g.order()))
  sage: mean = sum( 1/(g.degree(v)+1) for v in g )

Sage example in ./graphtheory.tex, line 1849::

  sage: while True:
  ....:    n = P.random_element()
  ....:    S = [v for v in g if all( n[v] < n[u] for u in g.neighbors(v))]
  ....:    if len(S) >= mean:
  ....:        break

Sage example in ./graphtheory.tex, line 1989::

  sage: def find_induced(H, G):
  ....:  # the function from V(H) to V(G) we aim to define:
  ....:  f = {}
  ....:  # set of vertices of G not yet used by f:
  ....:  G_remain = G.vertices()
  ....:  # set of vertices having no representative yet:
  ....:  H_remain = H.vertices()
  ....:  # while the function is not complete:
  ....:  while H_remain:
  ....:      v = H_remain.pop(0) # look for the next vertex of H
  ....:      # and its potential images in G
  ....:      candidates = [u for u in G_remain
  ....:                    if all(H.has_edge(h,v) == G.has_edge(f_h,u)
  ....:                           for h, f_h in f.items())]
  ....:      # if no candidate is found, we abort immediately
  ....:      if not candidates:
  ....:          raise ValueError("No copy of H has been found in G")
  ....:      # otherwise we select the first candidate
  ....:      f[v] = candidates[0]
  ....:      G_remain.remove(f[v])
  ....:  return f

Sage example in ./graphtheory.tex, line 2012::

  sage: set_random_seed(3)

Sage example in ./graphtheory.tex, line 2021::

  sage: H = graphs.PetersenGraph()
  sage: G = graphs.RandomGNP(500,0.5)
  sage: find_induced(H,G)
  {0: 0, 1: 4, 2: 3, 3: 7, 4: 35, 5: 10, 6: 67, 7: 108, 8: 240, 9: 39}

Sage example in ./graphtheory.tex, line 2070::

  sage: n = 100; V = range(n+1)
  sage: G = Graph()
  sage: G.add_edges([
  ....:    (i,j) for i,j in Subsets(V,2) if is_square(abs(i-j)) ])

Sage example in ./graphtheory.tex, line 2084::

  sage: X = G.independent_set(); X  # random with python3
  [4, 6, 9, 11, 16, 21, 23, 26, 28, 33, 38, 43, 50,
  56, 61, 71, 76, 78, 83, 88, 93, 95, 98, 100]
  sage: G.is_independent_set(X)
  True

Sage example in ./graphtheory.tex, line 2139::

  sage: tasks = {0: [2, 5, 3, 7], 1: [0, 1, 4],
  ....:   2: [5, 0, 4],           3: [0, 1],
  ....:   4: [8],                 5: [2],
  ....:   6: [8, 9, 7],           7: [5, 8, 7],
  ....:   8: [2, 5, 3, 6, 4],     9: [2, 5, 8, 6, 1]}
  sage: G = Graph()
  sage: for i in tasks:
  ....:   G.add_edges(("w" + str(i), "t" + str(j)) for j in tasks[i])

Sage example in ./graphtheory.tex, line 2182::

  sage: M = Graph(G.matching())
  sage: txt = "t{} assigned to {}"
  sage: for i in tasks:                                        # random
  ....:     print(txt.format(i, M.neighbors('t' + str(i))[0])) # random
  t0 assigned to w2
  t1 assigned to w3
  t2 assigned to w5
  t3 assigned to w8
  t4 assigned to w1
  t5 assigned to w7
  t6 assigned to w9
  t7 assigned to w0
  t8 assigned to w4
  t9 assigned to w6

Sage example in ./graphtheory.tex, line 2232::

  sage: n = 10
  sage: G = graphs.CompleteGraph(n)
  sage: from sage.graphs.graph_coloring import edge_coloring
  sage: for day, matches in enumerate(edge_coloring(G)):
  ....:    print("Matches of day {}: {}".format(day, matches))
  Matches of day 0: [(0, 9), (1, 8), (2, 7), (3, 6), (4, 5)]
  Matches of day 1: [(0, 2), (1, 9), (3, 8), (4, 7), (5, 6)]
  Matches of day 2: [(0, 4), (1, 3), (2, 9), (5, 8), (6, 7)]
  Matches of day 3: [(0, 6), (1, 5), (2, 4), (3, 9), (7, 8)]
  Matches of day 4: [(0, 8), (1, 7), (2, 6), (3, 5), (4, 9)]
  Matches of day 5: [(0, 1), (2, 8), (3, 7), (4, 6), (5, 9)]
  Matches of day 6: [(0, 3), (1, 2), (4, 8), (5, 7), (6, 9)]
  Matches of day 7: [(0, 5), (1, 4), (2, 3), (6, 8), (7, 9)]
  Matches of day 8: [(0, 7), (1, 6), (2, 5), (3, 4), (8, 9)]

Sage example in ./graphtheory.tex, line 2260::

  sage: g = graphs.CompleteGraph(10)
  sage: g.show(edge_colors=edge_coloring(g, hex_colors=True))

"""

