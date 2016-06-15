r"""
Yang-Baxter Graphs
"""
#*****************************************************************************
#       Copyright (C) 2009 Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.graphs.digraph import DiGraph
from sage.structure.sage_object import SageObject
from sage.misc.lazy_attribute import lazy_attribute
from sage.combinat.partition import Partition
from sage.combinat.permutation import Permutation

def YangBaxterGraph(partition=None, root=None, operators=None):
    r"""
    Construct the Yang-Baxter graph from ``root`` by repeated application of
    ``operators``, or the Yang-Baxter graph associated to ``partition``.

    INPUT:

    The user needs to provide either ``partition`` or both ``root`` and
    ``operators``, where

    - ``partition`` -- a partition of a positive integer

    - ``root`` -- the root vertex

    - ``operator`` - a function that maps vertices `u` to a list of
      tuples of the form `(v, l)` where `v` is a successor of `u` and `l` is
      the label of the edge from `u` to `v`.

    OUTPUT:

    - Either:

      - :class:`YangBaxterGraph_partition` - if partition is defined
      - :class:`YangBaxterGraph_generic` - if partition is ``None``

    EXAMPLES:

    The Yang-Baxter graph defined by a partition `[p_1,\dots,p_k]` is
    the labelled directed graph with vertex set obtained by
    bubble-sorting `(p_k-1,p_k-2,\dots,0,\dots,p_1-1,p_1-2,\dots,0)`;
    there is an arrow from `u` to `v` labelled by `i` if `v` is
    obtained by swapping the `i`-th and `(i+1)`-th elements of `u`.
    For example, if the partition is `[3,1]`, then we begin with
    `(0,2,1,0)` and generate all tuples obtained from it by swapping
    two adjacent entries if they are increasing::

        sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
        sage: bubbleswaps = [SwapIncreasingOperator(i) for i in range(3)]
        sage: Y = YangBaxterGraph(root=(0,2,1,0), operators=bubbleswaps); Y
        Yang-Baxter graph with root vertex (0, 2, 1, 0)
        sage: Y.vertices()
        [(2, 0, 1, 0), (2, 1, 0, 0), (0, 2, 1, 0)]

    The ``partition`` keyword is a shorthand for the above construction.

    ::

        sage: Y = YangBaxterGraph(partition=[3,1]); Y
        Yang-Baxter graph of [3, 1], with top vertex (0, 2, 1, 0)
        sage: Y.vertices()
        [(0, 2, 1, 0), (2, 0, 1, 0), (2, 1, 0, 0)]

    The permutahedron can be realized as a Yang-Baxter graph.

    ::

        sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
        sage: swappers = [SwapIncreasingOperator(i) for i in range(3)]
        sage: Y = YangBaxterGraph(root=(1,2,3,4), operators=swappers); Y
        Yang-Baxter graph with root vertex (1, 2, 3, 4)
        sage: Y.plot()
        Graphics object consisting of 97 graphics primitives

    The Cayley graph of a finite group can be realized as a Yang-Baxter graph.

    ::

        sage: def left_multiplication_by(g):
        ...       return lambda h : h*g
        sage: G = CyclicPermutationGroup(4)
        sage: operators = [ left_multiplication_by(gen) for gen in G.gens() ]
        sage: Y = YangBaxterGraph(root=G.identity(), operators=operators); Y
        Yang-Baxter graph with root vertex ()
        sage: Y.plot(edge_labels=False)
        Graphics object consisting of 9 graphics primitives

        sage: G = SymmetricGroup(4)
        sage: operators = [left_multiplication_by(gen) for gen in G.gens()]
        sage: Y = YangBaxterGraph(root=G.identity(), operators=operators); Y
        Yang-Baxter graph with root vertex ()
        sage: Y.plot(edge_labels=False)
        Graphics object consisting of 96 graphics primitives

    AUTHORS:

    - Franco Saliola (2009-04-23)
    """
    if partition is None:
        return YangBaxterGraph_generic(root=root, operators=operators)
    else:
        return YangBaxterGraph_partition(partition=Partition(partition))

##### General class for Yang-Baxter Graphs ################################

class YangBaxterGraph_generic(SageObject):
    def __init__(self, root, operators):
        r"""
        A class to model the Yang-Baxter graph defined by ``root`` and
        ``operators``.

        INPUT:

        - ``root`` -- the root vertex of the graph

        - ``operators`` -- a list of callables that map vertices to (new)
          vertices.


        .. NOTE::

           This is a lazy implementation: the digraph is only computed
           when it is needed.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(4)]
            sage: Y = YangBaxterGraph(root=(1,0,2,1,0), operators=ops); Y
            Yang-Baxter graph with root vertex (1, 0, 2, 1, 0)
            sage: loads(dumps(Y)) == Y
            True

        AUTHORS:

        - Franco Saliola (2009-04-23)
        """
        self._root = root
        self._operators = operators

    def _successors(self, u):
        r"""
        Return a list of tuples for the form `(op(u), op)`, where op
        is one of the operators defining ``self``.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(4)]
            sage: Y = YangBaxterGraph(root=(1,0,2,1,0), operators=ops)
            sage: Y._successors((1,0,2,1,0))
            [((1, 2, 0, 1, 0), Swap-if-increasing at position 1)]
        """
        successors = set()
        for op in self._operators:
            v = op(u)
            if v != u:
                successors.add((v, op))
        return list(successors)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(2)]
            sage: Y = YangBaxterGraph(root=(1,2,3), operators=ops)
            sage: Y.__repr__()
            'Yang-Baxter graph with root vertex (1, 2, 3)'
        """
        return "Yang-Baxter graph with root vertex %s" % (self._root,)

    @lazy_attribute
    def _digraph(self):
        r"""
        Constructs the underlying digraph and stores the result as an
        attribute.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(2)]
            sage: Y = YangBaxterGraph(root=(1,2,3), operators=ops)
            sage: Y._digraph
            Digraph on 6 vertices
        """
        digraph = DiGraph()
        digraph.add_vertex(self._root)
        queue = [self._root]
        while queue:
            u = queue.pop()
            for (v, l) in self._successors(u):
                if v not in digraph:
                    queue.append(v)
                digraph.add_edge(u, v, l)
        return digraph

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(2)]
            sage: Y = YangBaxterGraph(root=(1,2,3), operators=ops)
            sage: hash(Y)
            1028420699          # 32-bit
            7656306018247013467 # 64-bit
        """
        # TODO: this is ugly but unavoidable: the Yang Baxter graphs are being
        # used in containers but are mutable.
        return hash(self._digraph.copy(immutable=True))

    def __eq__(self, other):
        r"""
        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(4)]
            sage: Y1 = YangBaxterGraph(root=(1,0,2,1,0), operators=ops)
            sage: Y2 = YangBaxterGraph(root=(2,0,2,1,0), operators=ops)
            sage: Y3 = YangBaxterGraph(root=(1,0,2,1,0), operators=ops)
            sage: Y1.__eq__(Y2)
            False
            sage: Y2.__eq__(Y2)
            True
            sage: Y1.__eq__(Y1)
            True
            sage: Y3.__eq__(Y1)
            True
            sage: Y3.__eq__(Y2)
            False
        """
        return type(self) is type(other) and self._digraph == other._digraph

    def __ne__(self, other):
        r"""
        Test non-equality.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(4)]
            sage: Y1 = YangBaxterGraph(root=(1,0,2,1,0), operators=ops)
            sage: Y2 = YangBaxterGraph(root=(2,0,2,1,0), operators=ops)
            sage: Y3 = YangBaxterGraph(root=(1,0,2,1,0), operators=ops)
            sage: Y1.__ne__(Y2)
            True
            sage: Y2.__ne__(Y2)
            False
            sage: Y1.__ne__(Y1)
            False
            sage: Y3.__ne__(Y1)
            False
            sage: Y3.__ne__(Y2)
            True
        """
        return not self == other

    def __iter__(self):
        r"""
        Returns an iterator of the vertices in ``self``.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(4)]
            sage: Y = YangBaxterGraph(root=(1,0,2,1,0), operators=ops)
            sage: uniq(Y.__iter__())
            [(1, 0, 2, 1, 0), (1, 2, 0, 1, 0), (1, 2, 1, 0, 0), (2, 1, 0, 1, 0), (2, 1, 1, 0, 0)]
        """
        return self._digraph.vertex_iterator()

    def __len__(self):
        r"""
        Returns the number of vertices in ``self``.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(4)]
            sage: Y = YangBaxterGraph(root=(1,0,2,1,0), operators=ops)
            sage: Y.__len__()
            5
            sage: ops = [SwapIncreasingOperator(i) for i in range(5)]
            sage: Y = YangBaxterGraph(root=(0,1,0,2,1,0), operators=ops)
            sage: Y.__len__()
            16
        """
        return self._digraph.num_verts()

    def __copy__(self):
        r"""
        Returns a copy of ``self``.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(3)]
            sage: Y = YangBaxterGraph(root=(1,0,2,1,0), operators=ops); Y
            Yang-Baxter graph with root vertex (1, 0, 2, 1, 0)
            sage: B = copy(Y); B
            Yang-Baxter graph with root vertex (1, 0, 2, 1, 0)
            sage: Y is B
            False
            sage: Y == B
            True
        """
        from copy import copy
        Y = self.__class__(self._root, self._operators)
        Y._digraph = copy(self._digraph)
        return Y

    def _edges_in_bfs(self):
        r"""
        Returns an iterator of the edges of the digraph traversed in a
        breadth-first search of the vertices beginning at ``self.root()``.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(4)]
            sage: Y = YangBaxterGraph(root=(1,0,2,1,0), operators=ops)
            sage: list(Y._edges_in_bfs())
            [((1, 0, 2, 1, 0), (1, 2, 0, 1, 0), Swap-if-increasing at position 1), ((1, 2, 0, 1, 0), (1, 2, 1, 0, 0), Swap-if-increasing at position 2), ((1, 2, 0, 1, 0), (2, 1, 0, 1, 0), Swap-if-increasing at position 0), ((2, 1, 0, 1, 0), (2, 1, 1, 0, 0), Swap-if-increasing at position 2)]
        """
        digraph = self._digraph
        seen = {}
        queue = [self._root]
        seen[self._root] = True
        while queue:
            u = queue.pop()
            for w in digraph.neighbor_out_iterator(u):
                if w not in seen:
                    seen[w] = True
                    queue.append(w)
                    yield (u,w,digraph.edge_label(u,w))

    def root(self):
        r"""
        Returns the root vertex of ``self``.

        If ``self`` is the Yang-Baxter graph of the partition
        `[p_1,p_2,\dots,p_k]`, then this is the vertex
        `(p_k-1,p_k-2,\dots,0,\dots,p_1-1,p_1-2,\dots,0)`.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(4)]
            sage: Y = YangBaxterGraph(root=(1,0,2,1,0), operators=ops)
            sage: Y.root()
            (1, 0, 2, 1, 0)
            sage: Y = YangBaxterGraph(root=(0,1,0,2,1,0), operators=ops)
            sage: Y.root()
            (0, 1, 0, 2, 1, 0)
            sage: Y = YangBaxterGraph(root=(1,0,3,2,1,0), operators=ops)
            sage: Y.root()
            (1, 0, 3, 2, 1, 0)
            sage: Y = YangBaxterGraph(partition=[3,2])
            sage: Y.root()
            (1, 0, 2, 1, 0)
        """
        return self._root

    def successors(self, v):
        r"""
        Return the successors of the vertex ``v``.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(4)]
            sage: Y = YangBaxterGraph(root=(1,0,2,1,0), operators=ops)
            sage: Y.successors(Y.root())
            [(1, 2, 0, 1, 0)]
            sage: Y.successors((1, 2, 0, 1, 0))
            [(1, 2, 1, 0, 0), (2, 1, 0, 1, 0)]
        """
        return [a for (a,b) in self._successors(v)]

    def plot(self, *args, **kwds):
        r"""
        Plots ``self`` as a digraph.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(4)]
            sage: Y = YangBaxterGraph(root=(1,0,2,1,0), operators=ops)
            sage: Y.plot()
            Graphics object consisting of 16 graphics primitives
            sage: Y.plot(edge_labels=False)
            Graphics object consisting of 11 graphics primitives
        """
        if "edge_labels" not in kwds:
            kwds["edge_labels"] = True
        if "vertex_labels" not in kwds:
            kwds["vertex_labels"] = True
        return self._digraph.plot(*args, **kwds)

    def vertices(self):
        r"""
        Returns the vertices of ``self``.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(3)]
            sage: Y = YangBaxterGraph(root=(0,2,1,0), operators=ops)
            sage: Y.vertices()
            [(2, 0, 1, 0), (2, 1, 0, 0), (0, 2, 1, 0)]
        """
        return list(self)

    def edges(self):
        r"""
        Returns the (labelled) edges of ``self``.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(3)]
            sage: Y = YangBaxterGraph(root=(0,2,1,0), operators=ops)
            sage: Y.edges()
            [((0, 2, 1, 0), (2, 0, 1, 0), Swap-if-increasing at position 0), ((2, 0, 1, 0), (2, 1, 0, 0), Swap-if-increasing at position 1)]
        """
        return self._digraph.edges()

    def vertex_relabelling_dict(self, v, relabel_operator):
        r"""
        Return a dictionary pairing vertices ``u`` of ``self`` with
        the object obtained from ``v`` by applying the
        ``relabel_operator`` along a path from the root to ``u``. Note
        that the root is paired with ``v``.

        INPUT:

        - ``v`` -- an object

        - ``relabel_operator`` -- function mapping a vertex and a label to
          the image of the vertex

        OUTPUT:

        - dictionary pairing vertices with the corresponding image of ``v``

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(3)]
            sage: Y = YangBaxterGraph(root=(0,2,1,0), operators=ops)
            sage: def relabel_operator(op, u):
            ...       i = op.position()
            ...       return u[:i] + u[i:i+2][::-1] + u[i+2:]
            sage: Y.vertex_relabelling_dict((1,2,3,4), relabel_operator)
            {(0, 2, 1, 0): (1, 2, 3, 4),
             (2, 0, 1, 0): (2, 1, 3, 4),
             (2, 1, 0, 0): (2, 3, 1, 4)}
        """
        relabelling = {self._root:v}
        for (u,w,i) in self._edges_in_bfs():
            relabelling[w] = relabel_operator(i, relabelling[u])
        return relabelling

    def relabel_vertices(self, v, relabel_operator, inplace=True):
        r"""
        Relabel the vertices ``u`` of ``self`` by the object obtained
        from ``u`` by applying the ``relabel_operator`` to ``v`` along
        a path from ``self.root()`` to ``u``.

        Note that the ``self.root()`` is paired with ``v``.

        INPUT:

        - ``v`` -- tuple, Permutation, CombinatorialObject

        - ``inplace`` -- if ``True``, modifies ``self``; otherwise returns a
          modified copy of ``self``.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(3)]
            sage: Y = YangBaxterGraph(root=(0,2,1,0), operators=ops)
            sage: def relabel_op(op, u):
            ...       i = op.position()
            ...       return u[:i] + u[i:i+2][::-1] + u[i+2:]
            sage: d = Y.relabel_vertices((1,2,3,4), relabel_op, inplace=False); d
            Yang-Baxter graph with root vertex (1, 2, 3, 4)
            sage: Y.vertices()
            [(2, 0, 1, 0), (2, 1, 0, 0), (0, 2, 1, 0)]
            sage: e = Y.relabel_vertices((1,2,3,4), relabel_op); e
            sage: Y.vertices()
            [(2, 1, 3, 4), (1, 2, 3, 4), (2, 3, 1, 4)]
        """
        from copy import copy
        relabelling = self.vertex_relabelling_dict(v, relabel_operator)
        Y = self if inplace else copy(self)
        Y._root = relabelling[Y._root]
        Y._digraph.relabel(relabelling, inplace=True)
        if inplace is False:
            return Y

    def relabel_edges(self, edge_dict, inplace=True):
        r"""
        Relabel the edges of ``self``.

        INPUT:

        - ``edge_dict`` -- a dictionary keyed by the (unlabelled) edges.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: ops = [SwapIncreasingOperator(i) for i in range(3)]
            sage: Y = YangBaxterGraph(root=(0,2,1,0), operators=ops)
            sage: def relabel_op(op, u):
            ...       i = op.position()
            ...       return u[:i] + u[i:i+2][::-1] + u[i+2:]
            sage: Y.edges()
            [((0, 2, 1, 0), (2, 0, 1, 0), Swap-if-increasing at position 0), ((2, 0, 1, 0), (2, 1, 0, 0), Swap-if-increasing at position 1)]
            sage: d = {((0,2,1,0),(2,0,1,0)):17, ((2,0,1,0),(2,1,0,0)):27}
            sage: Y.relabel_edges(d, inplace=False).edges()
            [((0, 2, 1, 0), (2, 0, 1, 0), 17), ((2, 0, 1, 0), (2, 1, 0, 0), 27)]
            sage: Y.edges()
            [((0, 2, 1, 0), (2, 0, 1, 0), Swap-if-increasing at position 0), ((2, 0, 1, 0), (2, 1, 0, 0), Swap-if-increasing at position 1)]
            sage: Y.relabel_edges(d, inplace=True)
            sage: Y.edges()
            [((0, 2, 1, 0), (2, 0, 1, 0), 17), ((2, 0, 1, 0), (2, 1, 0, 0), 27)]
        """
        if inplace:
            Y = self
        else:
            from copy import copy
            Y = copy(self)
        digraph = Y._digraph
        for (u,v,i) in digraph.edges():
            digraph.set_edge_label(u,v,edge_dict[u,v])
        if not inplace:
            return Y

##### Yang-Baxter Graphs defined by a partition ###########################

class YangBaxterGraph_partition(YangBaxterGraph_generic):
    def __init__(self, partition):
        r"""
        A class to model the Yang-Baxter graph of a partition.

        The Yang-Baxter graph defined by a partition `[p_1,\dots,p_k]`
        is the labelled directed graph with vertex set obtained by
        bubble-sorting `(p_k-1,p_k-2,\dots,0,\dots,p_1-1,p_1-2,\dots,0)`;
        there is an arrow from `u` to `v` labelled by `i` if `v` is
        obtained by swapping the `i`-th and `(i+1)`-th elements of `u`.

        .. note::

           This is a lazy implementation: the digraph is only computed
           when it is needed.

        EXAMPLES::

            sage: Y = YangBaxterGraph(partition=[3,2,1]); Y
            Yang-Baxter graph of [3, 2, 1], with top vertex (0, 1, 0, 2, 1, 0)
            sage: loads(dumps(Y)) == Y
            True

        AUTHORS:

        - Franco Saliola (2009-04-23)
        """
        self._partition = partition
        beta = sorted(self._partition, reverse=True)
        root = tuple(sum(map(range,beta), []))[::-1]
        operators = [SwapIncreasingOperator(i) for i in range(sum(partition)-1)]
        super(YangBaxterGraph_partition, self).__init__(root, operators)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: Y = YangBaxterGraph(partition=[3,2])
            sage: Y.__repr__()
            'Yang-Baxter graph of [3, 2], with top vertex (1, 0, 2, 1, 0)'
        """
        return "Yang-Baxter graph of %s, with top vertex %s" % (self._partition, self._root)

    def __copy__(self):
        r"""
        Returns a copy of ``self``.

        EXAMPLES::

            sage: Y = YangBaxterGraph(partition=[3,2]); Y
            Yang-Baxter graph of [3, 2], with top vertex (1, 0, 2, 1, 0)
            sage: B = copy(Y); B
            Yang-Baxter graph of [3, 2], with top vertex (1, 0, 2, 1, 0)
            sage: Y is B
            False
            sage: Y == B
            True
        """
        from copy import copy
        Y = self.__class__(self._partition)
        Y._digraph = copy(self._digraph)
        return Y

    @lazy_attribute
    def _digraph(self):
        r"""
        Constructs the underlying digraph and stores the result as an
        attribute.

        EXAMPLES::

            sage: Y = YangBaxterGraph(partition=[2,1])
            sage: Y._digraph
            Digraph on 2 vertices
            sage: Y.edges()
            [((0, 1, 0), (1, 0, 0), Swap positions 0 and 1)]
        """
        digraph = super(YangBaxterGraph_partition, self)._digraph
        for (u,v,op) in digraph.edges():
            digraph.set_edge_label(u,v,SwapOperator(op.position()))
        return digraph

    @lazy_attribute
    def _vertex_ordering(self):
        r"""
        Returns a list of the vertices of ``self``, sorted using
        Pythons ``sorted`` method.

        EXAMPLES::

            sage: Y = YangBaxterGraph(partition=[3,2])
            sage: Y._vertex_ordering
            [(1, 0, 2, 1, 0), (1, 2, 0, 1, 0), (1, 2, 1, 0, 0), (2, 1, 0, 1, 0), (2, 1, 1, 0, 0)]
        """
        return sorted(self._digraph.vertices())

    def __iter__(self):
        r"""
        Iterate over the vertices ``self``.

        .. NOTE::

            The vertices are first sorted using Python's sorted command.

        EXAMPLES::

            sage: Y = YangBaxterGraph(partition=[3,2])
            sage: list(Y.__iter__())
            [(1, 0, 2, 1, 0), (1, 2, 0, 1, 0), (1, 2, 1, 0, 0), (2, 1, 0, 1, 0), (2, 1, 1, 0, 0)]
        """
        for v in self._vertex_ordering:
            yield v

    def _swap_operator(self, operator, u):
        r"""
        Return the image of ``u`` under ``operator``.

        INPUT:

        - ``i`` -- positive integer between 1 and len(u)-1, inclusive

        - ``u`` -- tuple, list, permutation, CombinatorialObject, ....

        EXAMPLES::

            sage: Y = YangBaxterGraph(partition=[3,1])
            sage: from sage.combinat.yang_baxter_graph import SwapOperator
            sage: ops = [SwapOperator(i) for i in range(3)]
            sage: [Y._swap_operator(op, (1,2,3,4)) for op in ops]
            [(2, 1, 3, 4), (1, 3, 2, 4), (1, 2, 4, 3)]
            sage: [Y._swap_operator(op, [4,3,2,1]) for op in ops]
            [[3, 4, 2, 1], [4, 2, 3, 1], [4, 3, 1, 2]]
            sage: [Y._swap_operator(op, Permutation([1,2,3,4])) for op in ops]
            [[2, 1, 3, 4], [1, 3, 2, 4], [1, 2, 4, 3]]
        """
        return operator(u)

    def vertex_relabelling_dict(self, v):
        r"""
        Return a dictionary pairing vertices ``u`` of ``self`` with the object
        obtained from ``v`` by applying transpositions corresponding to the
        edges labels along a path from the root to ``u``.

        Note that the root is paired with ``v``.

        INPUT:

        - ``v`` -- an object

        OUTPUT:

        - dictionary pairing vertices with the corresponding image of ``v``

        EXAMPLES::

            sage: Y = YangBaxterGraph(partition=[3,1])
            sage: Y.vertex_relabelling_dict((1,2,3,4))
            {(0, 2, 1, 0): (1, 2, 3, 4),
             (2, 0, 1, 0): (2, 1, 3, 4),
             (2, 1, 0, 0): (2, 3, 1, 4)}
            sage: Y.vertex_relabelling_dict((4,3,2,1))
            {(0, 2, 1, 0): (4, 3, 2, 1),
             (2, 0, 1, 0): (3, 4, 2, 1),
             (2, 1, 0, 0): (3, 2, 4, 1)}
        """
        return super(YangBaxterGraph_partition, self).vertex_relabelling_dict(v, self._swap_operator)

    def relabel_vertices(self, v, inplace=True):
        r"""
        Relabel the vertices of ``self`` with the object obtained from
        ``v`` by applying the transpositions corresponding to the edge
        labels along some path from the root to the vertex.

        INPUT:

        - ``v`` -- tuple, Permutation, CombinatorialObject

        - ``inplace`` -- if ``True``, modifies ``self``; otherwise
          returns a modified copy of ``self``.

        EXAMPLES::

            sage: Y = YangBaxterGraph(partition=[3,1]); Y
            Yang-Baxter graph of [3, 1], with top vertex (0, 2, 1, 0)
            sage: d = Y.relabel_vertices((1,2,3,4), inplace=False); d
            Digraph on 3 vertices
            sage: Y.vertices()
            [(0, 2, 1, 0), (2, 0, 1, 0), (2, 1, 0, 0)]
            sage: e = Y.relabel_vertices((1,2,3,4)); e
            sage: Y.vertices()
            [(1, 2, 3, 4), (2, 1, 3, 4), (2, 3, 1, 4)]
        """
        relabelling = self.vertex_relabelling_dict(v)
        if inplace:
            Y = self
            Y._root = relabelling[Y._root]
            Y._digraph.relabel(relabelling, inplace=inplace)
            Y._vertex_ordering = sorted(Y._digraph.vertices())
            return
        else:
            from copy import copy
            Y = copy(self)
            Y._root = relabelling[Y._root]
            return Y._digraph.relabel(relabelling, inplace=inplace)

##### Some Yang-Baxter operators ##########################################

class SwapOperator(SageObject):
    def __init__(self, i):
        r"""
        The operator that swaps the items in positions ``i`` and ``i+1``.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapOperator
            sage: s3 = SwapOperator(3)
            sage: s3 == loads(dumps(s3))
            True
        """
        self._position = i

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.combinat.yang_baxter_graph import SwapOperator
            sage: s = [SwapOperator(i) for i in range(3)]
            sage: map(hash, s)
            [0, 1, 2]
        """
        return hash(self._position)

    def __cmp__(self, other):
        r"""
        Compare two swap operators. The comparison is done by comparing the
        positions.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapOperator
            sage: s = [SwapOperator(i) for i in range(3)]
            sage: s[0] == s[0]
            True
            sage: s[1] < s[0]
            False
            sage: s[1] < s[2]
            True
        """
        if type(self) is type(other):
            return cmp(self._position, other._position)
        else:
            return cmp(type(self), type(other))

    def __repr__(self):
        r"""
        Representation string.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapOperator
            sage: s3 = SwapOperator(3)
            sage: s3.__repr__()
            'Swap positions 3 and 4'
        """
        return "Swap positions %s and %s" % (self._position, self._position+1)

    def __str__(self):
        r"""
        A short str representation (used, for example, in labelling edges of
        graphs).

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapOperator
            sage: s3 = SwapOperator(3)
            sage: s3.__str__()
            '3'
        """
        return "%s" % self._position

    def __call__(self, u):
        r"""
        Return the object obtained from swapping the items in positions
        ``i`` and ``i+1`` of ``u``.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapOperator
            sage: s3 = SwapOperator(3)
            sage: s3((1,2,3,4,5))
            (1, 2, 3, 5, 4)
            sage: s3([1,2,3,4,5])
            [1, 2, 3, 5, 4]
        """
        i = self._position
        if isinstance(u, Permutation):
            return Permutation(u[:i] + u[i:i+2][::-1] + u[i+2:])
        return type(u)(u[:i] + u[i:i+2][::-1] + u[i+2:])

    def position(self):
        r"""
        ``self`` is the operator that swaps positions ``i`` and ``i+1``. This
        method returns ``i``.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapOperator
            sage: s3 = SwapOperator(3)
            sage: s3.position()
            3
        """
        return self._position

class SwapIncreasingOperator(SwapOperator):
    def __repr__(self):
        r"""
        Representation string.

        EXAMPLES::

            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: s3 = SwapIncreasingOperator(3)
            sage: s3.__repr__()
            'Swap-if-increasing at position 3'
        """
        return "Swap-if-increasing at position %s" % self._position

    def __call__(self, u):
        r"""
        Returns a copy of ``u`` with ``u[i-1]`` and ``u[i]`` swapped if
        ``u[i-1] > u[i]``; otherwise returns ``u``.

        INPUT:

        - ``i`` -- positive integer between ``1`` and ``len(u)-1``, inclusive

        - ``u`` -- tuple, list, permutation, CombinatorialObject, ....

        EXAMPLES::

            sage: Y = YangBaxterGraph(partition=[2,2])
            sage: from sage.combinat.yang_baxter_graph import SwapIncreasingOperator
            sage: operators = [SwapIncreasingOperator(i) for i in range(3)]
            sage: [op((1,2,3,4)) for op in operators]
            [(2, 1, 3, 4), (1, 3, 2, 4), (1, 2, 4, 3)]
            sage: [op([4,3,2,1]) for op in operators]
            [[4, 3, 2, 1], [4, 3, 2, 1], [4, 3, 2, 1]]
            sage: [op(Permutation([1,3,2,4])) for op in operators]
            [[3, 1, 2, 4], [1, 3, 2, 4], [1, 3, 4, 2]]
        """
        i = self._position
        j = i+1
        if u[i] < u[j]:
            v = list(u)
            (v[j], v[i]) = (v[i], v[j])
            if isinstance(u, Permutation):
                return Permutation(v)
            return type(u)(v)
        else:
            return u

