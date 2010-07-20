r"""
Graphs associated to the triangle group (2,3,Infinity)

A subgroup of a free product of cyclic group is also a free product of cyclic
group. This module implements various method relative to those groups:
    - enumerate subgroups of finite index
    - structure of subgroups from their cosets (~free generators of subgroup)
    - find spanning trees with optimal property (~fundamental domain)

This method can be applied to Fuchsian groups for which the underlying surface
is a free group as PSL(2,Z) = C_2 * C_3 and more generally the triangle groups
D(p,q,infty) = C_p * C_q.

TODO:

- see how we can feet the general graph scheme (ask Rob)
- more general triangle graphs
- how do we build the coset graph for Gamma0, Gamma1 and others...
"""
from sage.structure.sage_object import SageObject
from sage.plot.misc import options
from sage.rings.integer import Integer
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.matrix.constructor import matrix

class TriangleGraph_2_3_infinity(SageObject):
    r"""
    A (connected) triangle graph (2,3,infinity) corresponds to a (transitive)
    group action of the free product C2*C3 on a finite set. It is an oriented
    graph such that

    - there are two type of edges (type 2 and type 3)
    - each vertex has exactly one incoming and one outgoing edges of each type
      (hence four adjacent edges)
    - the induced graph with edges of type 2 is made of loops and cycles of
      length 2
    - the induced graph with edges of type 3 is made of loops and cycles of
      length 3

    It can alterantively be defined by a pair of permutations (s2,s3) on [0,n-1]
    such that

    - s2 is of order 2
    - s3 is of order 3


    The data structure is the following:

     * _n is the number ov vertices
     * _incoming_edges is a list of two lists of n integers
     * _outgoing_edges is a list of two lists of n integers

    EXAMPLES::

        sage: t = TriangleGraph_2_3_infinity()
        Triangle graph (2,3,infinty) with 1 vertices
        sage: t.create_cycle(0,type=3)
        sage: t
        Triangle graph(2,3,infinity) with 3 vertices
    """
    def __init__(self, data=None, check=True):
        if data is None:
            self._n = 1
            self._incoming_edges = [0, 0]
            self._outgoing_edges = [0, 0]

        elif isinstance(data, (int,Integer)):
            self._n = data
            self._incoming_edges = [range(data), range(data)]
            self._outgoing_edges = [range(data), range(data)]

        elif isinstance(data, (list,tuple)):
            if len(data) != 2: 
                raise ValueError, "must have 2 list"
            if len(data[0]) != len(data[1]):
                raise ValueError, "the two list have different length"
            self._n = len(data[0])
            self._incoming_edges = [[0]*self._n, [0]*self._n]
            self._outgoing_edges = [[0]*self._n, [0]*self._n]

            for i,j in enumerate(data[0]):
                self._outgoing_edges[0][i] = j
                self._incoming_edges[0][j] = i
            for i,j in enumerate(data[1]):
                self._outgoing_edges[1][i] = j
                self._incoming_edges[1][j] = i

            if check: self._check()
        
        else:
            raise ValueError, "data must be None, a number or two lists"

    def _check(self):
        r"""
        perform a check to verify that I am a triangle(2,3,infinty) graph.
        """
        # check the lengths of the lists
        if (len(self._incoming_edges[0]) != self._n or
            len(self._incoming_edges[1]) != self._n or
            len(self._outgoing_edges[0]) != self._n or
            len(self._outgoing_edges[1]) != self._n):
            raise TypeError, "wrong length"

        # check that incoming are the reversed of outgoing 
        for i in xrange(self._n):
            if self._outgoing_edges[0][self._incoming_edges[0][i]] != i:
                raise TypeError, "error for type 2 edges on %d" %i
            if self._incoming_edges[0][self._outgoing_edges[0][i]] != i:
                raise TypeError, "error for type 2 edges on %d" %i
            if self._outgoing_edges[1][self._incoming_edges[1][i]] != i:
                raise TypeError, "error for type 3 edges on %d" %i
            if self._incoming_edges[1][self._outgoing_edges[1][i]] != i:
                raise TypeError, "error for type 3 edges on %d" %i

        # check length of the cycles
        for i in xrange(self._n):
            i1 = self._outgoing_edges[0][i]
            if self._outgoing_edges[0][i1] != i:
                raise TypeError, "not a 2-loop from %d" %i

            i1 = self._outgoing_edges[1][i]
            i2 = self._outgoing_edges[1][i1]
            if self._outgoing_edges[1][i2] != i:
                raise TypeError, "not a 3-loop from %d" %i

    def _repr_(self):
        r"""
        Return a string representation of this graph.
        """
        return "Triangle graph (2,3,infinty) with %d vertices" %self._n

    # numerical property

    def num_edges(self):
        r"""
        Number of edges

        (very unuseful)
        """
        return self._n * 4

    def num_verts(self):
        r"""
        Number of vertices
        """
        return self._n

    def num_loops(self, type=None):
        if type is None:
            i2 = 0; i3 = 0
            for i in xrange(self._n):
                if self._outgoing_edges[0][i] == i: i2 += 1
                if self._outgoing_edges[1][i] == i: i3 += 1
            return (i2,i3)

        elif type == 2:
            i2 = 0
            for i in xrange(self._n):
                if self._outgoing_edges[0][i] == i: i2 += 1
            return i2

        elif type == 3:
            i3 = 0
            for i in xrange(self._n):
                if self._outgoing_edges[1][i] == i: i3 += 1

        else:
            raise TypeError, "type must be 2 or 3"


    def graph(self):
        r"""
        (temporary) returns the underlying graph
        """
        from sage.graphs.digraph import DiGraph
        res = DiGraph(multiedges=True,loops=True)
        res.add_vertices(range(self._n))
        res.add_edges((i,self._outgoing_edges[0][i],2) for i in xrange(self._n))
        res.add_edges((i,self._outgoing_edges[1][i],3) for i in xrange(self._n))

        return res

    @options(graph_border=True,color_by_label=True)
    def plot(self,**options):
        r"""
        plot this graph
        """
        g = self.graph()
        return g.plot(**options)

    def show(self, **options):
        return self.plot().show(**options)


    # access to edges

    def __getitem__(self, i):
        r"""
        either i
        either i,type
        """
        if isinstance(i, (int,Integer)):
            return (self._outgoing_edges[0][i],self._outgoing_edges[1][i])
        elif isinstance(i, (tuple,list)):
            if i[1] == 2:
                return self._outgoing_edges[0][i[0]]
            elif i[1] == 3:
                return self._outgoing_edges[1][i[0]]

    def outgoing_edges(self, i, type=None):
        if type is None:
            return [(i,self._outgoing_edges[0][i],2),
                    (i,self._outgoing_edges[1][i],3)]
        elif type == 2:
            return (i, self._outgoing_edges[0][i],2)
        elif type == 3:
            return (i, self._outgoing_edges[1][i],3)
        else:
            return []
            
    def incoming_edges(self,i):
        if type is None:
            return [(i, self._incoming_edges[0][i],2),
                    (i, self._incoming_edges[1][i],3)]
        elif type == 2:
            return (i, self._incoming_edges[0][i],2)
        elif type == 3:
            return (i, self._incoming_edges[0][i],3)
        else:
            return []

    def get_permutations(self, type=None):
        r"""
        Return the permutations associated to the edges

        An alternative representation of a triangle graph is by a pair of
        permutations of order 2 and 3. This function returns the corresponding
        permutations.

        INPUT:

        - ``type`` - 2 or 3 - the type of the edge
        """
        if type is None:
            return self._outgoing_edges[0][:], self._outgoing_edges[1][:]
        elif type == 2:
            return self._outgoing_edges[0][:]
        elif type == 3:
            return self._outgoing_edges[1][:]

    # graph modification
    #   - create cycle: add vertices and one cycle
    #   - fusion cycle: do not touch vertices and add one cycle
    #   - delete cycle: do not touch vertices and delete one cycle
    # the latter operation could result with a non connected graph

    def create_cycle(self, vertex, type):
        r"""
        Create a cycle with self loop on it from the vertex i

        This method inserts a non trivial cycle (i,n) (or (i,n,n+1)) and trivial
        for the other type (n) (or (n)(n+1)).

        The vertex on which the operation is performed must be a fixed point of
        the corresponding permutation (or in other word it must have a loop on
        it).

        INPUT::

        - ``vertex`` - a vertex

        - ``type`` - a type of edge
        """
        if type == 2:
            if self._outgoing_edges[0][vertex] != vertex:
                raise ValueError, "no loop on vertex %d" %vertex

            self._outgoing_edges[0][vertex] = self._n
            self._outgoing_edges[0].append(vertex)
            self._incoming_edges[0][vertex] = self._n
            self._incoming_edges[0].append(vertex)

            self._outgoing_edges[1].append(self._n)
            self._incoming_edges[1].append(self._n)

            self._n += 1

        elif type == 3:
            if self._outgoing_edges[1][vertex] != vertex:
                raise ValueError, "no loop on vertes %d" %vertex

            self._outgoing_edges[1][vertex] = self._n
            self._outgoing_edges[1].append(self._n+1)
            self._outgoing_edges[1].append(vertex)

            self._incoming_edges[1][vertex] = self._n+1
            self._incoming_edges[1].append(vertex)
            self._incoming_edges[1].append(self._n)

            self._outgoing_edges[0].append(self._n)
            self._outgoing_edges[0].append(self._n+1)
            self._incoming_edges[0].append(self._n)
            self._incoming_edges[0].append(self._n+1)

            self._n += 2

    def fusion_cycle(self, vertices):
        r"""
        Create a cycle from 2 on a list of vertices with loops.

        (i)(j) -> (i,j) or (i)(j)(k) -> (i,j,k)

        INPUT:

        - ``vertices`` - list of vertices
        """
        if not isinstance(vertices, (tuple,list)):
            raise TypeError, "need a list or a tuple"
        if len(vertices) == 2:
            if any(self._outgoing_edges[0][v] != v for v in vertices):
                raise ValueError, "need loops on the two vertices"

            v0,v1 = vertices
            self._outgoing_edges[0][v0] = v1
            self._outgoing_edges[0][v1] = v0

            self._incoming_edges[0][v0] = v1
            self._incoming_edges[0][v1] = v0

        elif len(l) == 3:
            if any(self._outgoing_edges[1][v] != v for v in vertices):
                raise ValueError, "need loops on the three vertices"

            v0,v1,v2 = vertices
            self._outgoing_edges[1][v0] = v1
            self._outgoing_edges[1][v1] = v2
            self._outgoing_edges[1][v2] = v0

            self._incoming_edges[1][v0] = v2
            self._incoming_edges[1][v1] = v0
            self._incoming_edges[1][v2] = v1

        else:
            raise ValueError, "the length of the list of vertices must be 2 or 3"

    def delete_cycle(self, vertex, type):
        r"""
        Delete the cycle of type ``type`` that contains the vertex ``vertex``
        """
        if not isinstance(vertex, (int,Integer)):
            raise TypeError, "wait for an integer"
        if type == 2:
            if self._outgoing_edges[0][vertex] == vertex:
                raise ValueError, "not a bigon"
            v1 = self._outgoing_edges[0][vertex]
            self._outgoing_edges[0][vertex] = vertex
            self._incoming_edges[0][vertex] = vertex
            self._outgoing_edges[0][v1] = v1
            self._incoming_edges[0][v1] = v1
        
        elif type == 3:
            if self._outgoing_edges[1][vertex] == vertex:
                raise ValueError, "not a triangle"
            v1 = self._outgoing_edges[1][vertex]
            v2 = self._incoming_edges[1][vertex]

            self._outgoing_edges[vertex] = vertex
            self._incoming_edges[vertex] = vertex
            self._outgoing_edges[v1] = v1
            self._incoming_edges[v1] = v1
            self._outgoing_edges[v2] = v2
            self._incoming_edges[v2] = v2

        else:
            raise ValueError, "type must be 2 or 3"


    def spanning_tree(self, root=0, id=None, g2=None, g3=None, verbose=False):
        r"""
        Build a spanning tree respecting Kulkarni's method

        INPUT:

        - ``root`` - the right coset that corresponds to the identity

        - ``g2`` - an element of order 2

        - ``g3`` - an element of order 3


        OUTPUT:

        - ``tree`` - a spanning tree

        - ``gens`` - a list of free generators (C2 gens, C3 gens, Cinfinity gens)

        - ``lifts`` - lifts of the element in SL(2,ZZ)

        TODO:
        - add some random stuff to be able to get all possible spanning trees
        - do not calcul the distance at the begining
        """
        if id is None:
            id = matrix([[1,0],[0,1]])
        if g2 is None:
            g2 = matrix([[0,-1],[1,0]])
        if g3 is None:
            g3 = matrix([[1,-1],[1,0]])
        g3_inv = g3**-1

        distances = self.graph().to_undirected().shortest_path_lengths(root)

        if verbose:
            print "distances to the root\n"
            print distances

        # the tree and the lift
        from sage.graphs.digraph import DiGraph
        tree = DiGraph(multiedges=False,loops=False)
        lifts = [None] * self._n
        gens = [[],[],[]]

        x0 = root
        lifts[x0] = id
        tree.add_vertex(x0)
        l = set([x0])

        while True:
            # complete the triangle in the tree
            if self._outgoing_edges[1][x0] != x0:
                if verbose:
                    print "complete triangle from %d" %x0

                x1 = self._outgoing_edges[1][x0]
                x2 = self._outgoing_edges[1][x1]
                
                if verbose: 
                    print "  find %d %d" %(x1,x2)
                    
                tree.add_vertices([x1,x2])
                l.add(x1); tree.add_edge(x0,x1,'3')
                l.add(x2); tree.add_edge(x0,x2,'3inv')
                lifts[x1] = lifts[x0] * g3
                lifts[x2] = lifts[x0] * g3_inv
               
            else:
                if verbose:
                    print "elliptic 3 point at %d" %x0

                gens[1].append(lifts[x0] * g3 * lifts[x0]**-1)
            
            # at this step all triangles are completed
            # perform S edges link and gens
            if verbose: print "perform g2 links"
            while l:
                x1 = l.pop()
                x0 = self._outgoing_edges[0][x1]

                if x1 != x0:
                    if verbose:
                        print "g2 link %d %d" %(x0,x1)
                    
                    if tree.has_vertex(x0):
                        gens[2].append(lifts[x1] * g2 * lifts[x0]**-1)
                        if x0 in l: l.remove(x0)  # x0 must be in l
                    else:
                        tree.add_edge(x1,x0,'2')
                        lifts[x0] = lifts[x1]*g2
                        x1 = self._outgoing_edges[1][x0]
                        x2 = self._incoming_edges[1][x0]

                        if verbose: 
                            print "new vertex", x0
                            print "3-neighbors are %d %d"  %(x1,x2)
                        if (distances[x1] > distances[x0] or distances[x2] > distances[x0]):
                            if verbose:
                                print "the neighbors are far from the root...break g2 links"
                            break
                        
            else:
                if verbose: print "quit"
                break
            
        return tree,gens,lifts
