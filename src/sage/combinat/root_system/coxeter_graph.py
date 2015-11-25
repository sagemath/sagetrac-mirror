"""
Coxeter graphs
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#       Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#       Copyright (C) 2015 Jean-Philippe Labbe <labbe@math.fu-berlin.de>
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
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.matrix.matrix import is_Matrix
from sage.graphs.graph import Graph
from sage.rings.all import ZZ
from sage.combinat.root_system.coxeter_type import CoxeterType 
from sage.combinat.root_system.coxeter_matrix import CoxeterMatrix

"""
def CoxeterGraph(*args, **kwds):
    
    Return the Coxeter graph corresponding to the input.

    INPUT:

    The input can be one of the following:

    - empty to obtain an empty Coxeter graph
    - a Coxeter type
    - a Coxeter matrix
    - a generalized Coxeter matrix
    - a (generalized) Coxeter matrix and an indexing set


    EXAMPLES::

        sage: CoxeterGraph(['A', 4])
        Graph on 4 vertices

        sage: CoxeterGraph(['A',1],['A',1])
        Graph on 2 vertices

        sage: R = RootSystem("A2xB2xF4")
        sage: CoxeterGraph(R)
        Graph on 8 vertices

        sage: R = RootSystem("A2xB2xF4")
        sage: CM = R.coxeter_matrix(); CM
        [1 3 2 2 2 2 2 2]
        [3 1 2 2 2 2 2 2]
        [2 2 1 4 2 2 2 2]
        [2 2 4 1 2 2 2 2]
        [2 2 2 2 1 3 2 2]
        [2 2 2 2 3 1 4 2]
        [2 2 2 2 2 4 1 3]
        [2 2 2 2 2 2 3 1]
        sage: DD = CoxeterGraph(CM); DD
        Graph on 8 vertices
        sage: DD.coxeter_matrix()
        [ 2 -1  0  0  0  0  0  0]
        [-1  2  0  0  0  0  0  0]
        [ 0  0  2 -1  0  0  0  0]
        [ 0  0 -2  2  0  0  0  0]
        [ 0  0  0  0  2 -1  0  0]
        [ 0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0 -2  2 -1]
        [ 0  0  0  0  0  0 -1  2]

    We can also create Coxeter graphs from arbitrary Coxeter matrices::

        sage: C = CoxeterMatrix([[1, -1], [-1, 1]])
        sage: CoxeterGraph(C)
        Coxeter graph of rank 2
        sage: C.index_set()
        (0, 1)
        sage: CI = CoxeterMatrix([[1, -2], [-2, 1]])
        sage: DI = CoxeterGraph(CI)
        sage: DI.index_set()
        (3, 5)
        sage: CII = CoxeterMatrix([[1, 5], [5, 1]])
        sage: DII = CoxeterGraph(CII, ('y', 'x'))
        sage: DII.index_set()
        ('x', 'y')

    .. SEEALSO::

        :func:`CoxeterType` for a general discussion on Coxeter
        types and in particular node labeling conventions.

    TESTS:

    Check that :trac:`15277` is fixed by not having edges from 0's::

        sage: CM =
        CoxeterMatrix([[1,-1,0,0],[-3,1,-2,-2],[0,-1,1,-1],[0,-1,-1,1]])
        sage: CM
        [ 2 -1  0  0]
        [-3  2 -2 -2]
        [ 0 -1  2 -1]
        [ 0 -1 -1  2]
        sage: CM.coxeter_graph().edges()
        [(0, 1, 3),
         (1, 0, 1),
         (1, 2, 1),
         (1, 3, 1),
         (2, 1, 2),
         (2, 3, 1),
         (3, 1, 2),
         (3, 2, 1)]
    
    if len(args) == 0:
        return CoxeterGraph_class()
    mat = args[0]
    if is_Matrix(mat):
        mat = CoxeterMatrix(*args)
    if isinstance(mat, CoxeterMatrix):
        if mat.coxeter_type() is not mat:
            try:
                return mat.coxeter_type().coxeter_graph()
            except AttributeError:
                ct = CoxeterType(*args)
                raise ValueError("Coxeter graph data not yet hardcoded for type %s"%ct)
        if len(args) > 1:
            index_set = tuple(args[1])
        elif "index_set" in kwds:
            index_set = tuple(kwds["index_set"])
        else:
            index_set = mat.index_set()
        D = CoxeterGraph_class(index_set=index_set)
        for (i, j) in mat.nonzero_positions():
            if i != j:
                D.add_edge(index_set[i], index_set[j], -mat[j, i])
        return D
    ct = CoxeterType(*args)
    try:
        return ct.coxeter_graph()
    except AttributeError:
        raise ValueError("Coxeter graph data not yet hardcoded for type %s"%ct)
"""

class CoxeterGraph(CoxeterType):
    r"""
    A generalized Coxeter graph.

    INPUT:

    The input can be one of the following:

    - empty to obtain an empty Coxeter graph
    - a Coxeter type
    - a Coxeter matrix
    - a generalized Coxeter matrix
    - a (generalized) Coxeter matrix and an indexing set
    - ``t`` -- a Coxeter type, a generalized Coxeter matrix, or ``None``

    A generalized Coxeter graph is a graph encoding
    a Coxeter system `(W, S)` and a bilinear form, where the relations are given by
    `(s_i s_j)^{m_{ij}}`. Thus `M` is symmetric and has entries
    in `\{1, 2, 3, \ldots, \infty\}` with `m_{ij} = 1` if and only
    if `i = j`.

    We represent `m_{ij} = \infty` by any number `m_{ij} \leq -1`. In
    particular, we can construct a bilinear form `B = (b_{ij})_{i,j \in I}`
    from `M` by

    .. MATH::

        b_{ij} = \begin{cases}
        m_{ij} & m_{ij} < 0\ (\text{i.e., } m_{ij} = \infty), \\
        -\cos\left( \frac{\pi}{m_{ij}} \right) & \text{otherwise}.
        \end{cases}
        
    The vertices of the generalized Coxeter graph are indexed by the elements of `S`.
    If `m_{ij} = 2` there is no edge between `i` and `j`.
    Otherwise, we label the edge between `i` and `j` by the value `m_{ij}`.

    EXAMPLES::

        sage: CoxeterGraph(['A', 3])
        O---O---O
        1   2   3
        A3
        sage: C = CoxeterMatrix([[1, -3], [-3, 1]])
        sage: CoxeterGraph(C)
        Coxeter graph of rank 2
        sage: C.coxeter_graph().coxeter_matrix() == C
        True
        
        

        sage: CoxeterGraph(['A', 4])
        Graph on 4 vertices

        sage: CoxeterGraph(['A',1],['A',1])
        Graph on 2 vertices

        sage: R = RootSystem("A2xB2xF4")
        sage: CoxeterGraph(R)
        Graph on 8 vertices

        sage: R = RootSystem("A2xB2xF4")
        sage: CM = R.coxeter_matrix(); CM
        [1 3 2 2 2 2 2 2]
        [3 1 2 2 2 2 2 2]
        [2 2 1 4 2 2 2 2]
        [2 2 4 1 2 2 2 2]
        [2 2 2 2 1 3 2 2]
        [2 2 2 2 3 1 4 2]
        [2 2 2 2 2 4 1 3]
        [2 2 2 2 2 2 3 1]
        sage: DD = CoxeterGraph(CM); DD
        Graph on 8 vertices
        sage: DD.coxeter_matrix()
        [ 2 -1  0  0  0  0  0  0]
        [-1  2  0  0  0  0  0  0]
        [ 0  0  2 -1  0  0  0  0]
        [ 0  0 -2  2  0  0  0  0]
        [ 0  0  0  0  2 -1  0  0]
        [ 0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0 -2  2 -1]
        [ 0  0  0  0  0  0 -1  2]

    We can also create Coxeter graphs from arbitrary Coxeter matrices::

        sage: C = CoxeterMatrix([[1, -1], [-1, 1]])
        sage: CoxeterGraph(C)
        Coxeter graph of rank 2
        sage: C.index_set()
        (0, 1)
        sage: CI = CoxeterMatrix([[1, -2], [-2, 1]])
        sage: DI = CoxeterGraph(CI)
        sage: DI.index_set()
        (3, 5)
        sage: CII = CoxeterMatrix([[1, 5], [5, 1]])
        sage: DII = CoxeterGraph(CII, ('y', 'x'))
        sage: DII.index_set()
        ('x', 'y')


    TESTS:

    Check that the correct type is returned when copied::

        sage: d = CoxeterGraph(['A', 3])
        sage: type(copy(d))
        <class 'sage.combinat.root_system.coxeter_graph.CoxeterGraph_class'>

    We check that :trac:`14655` is fixed::

        sage: cd = copy(d)
        sage: cd.add_vertex(4)
        sage: d.vertices() != cd.vertices()
        True
        
    Check that :trac:`15277` is fixed by not having edges from 0's::

        sage: CM =
        CoxeterMatrix([[1,-1,0,0],[-3,1,-2,-2],[0,-1,1,-1],[0,-1,-1,1]])
        sage: CM
        [ 2 -1  0  0]
        [-3  2 -2 -2]
        [ 0 -1  2 -1]
        [ 0 -1 -1  2]
        sage: CM.coxeter_graph().edges()
        [(0, 1, 3),
         (1, 0, 1),
         (1, 2, 1),
         (1, 3, 1),
         (2, 1, 2),
         (2, 3, 1),
         (3, 1, 2),
         (3, 2, 1)]

    Implementation note: if a Coxeter type is given, then the nodes
    are initialized from the index set of this Coxeter type.
    
    .. SEEALSO::

        :func:`CoxeterType` for a general discussion on Coxeter
        types and in particular node labeling conventions.

    """
    __metaclass__ = ClasscallMetaclass
    
    @staticmethod
    def __classcall_private__(cls, data=None, index_set=None, coxeter_type=None,
                              cartan_type=None, coxeter_type_check=True):
        r"""
        A Coxeter graph can we created via a graph, a Coxeter type, or
        a matrix.

        .. NOTE::

            To disable the Coxeter type check, use the optional argument
            ``coxeter_type_check = False``.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',1,1],['a','b'])
            sage: C2 = CoxeterMatrix([[1, -1], [-1, 1]])
            sage: C3 = CoxeterMatrix(matrix([[1, -1], [-1, 1]]), [0, 1])
            sage: C == C2 and C == C3
            True

        Check with `\infty` because of the hack of using `-1` to represent
        `\infty` in the Coxeter matrix::

            sage: G = Graph([(0, 1, 3), (1, 2, oo)])
            sage: W1 = CoxeterMatrix([[1, 3, 2], [3, 1, -1], [2, -1, 1]])
            sage: W2 = CoxeterMatrix(G)
            sage: W1 == W2
            True
            sage: CoxeterMatrix(W1.coxeter_graph()) == W1
            True

        The base ring of the matrix depends on the entries given::

            sage: CoxeterMatrix([[1,-1],[-1,1]])._matrix.base_ring()
            Integer Ring
            sage: CoxeterMatrix([[1,-3/2],[-3/2,1]])._matrix.base_ring()
            Rational Field
            sage: CoxeterMatrix([[1,-1.5],[-1.5,1]])._matrix.base_ring()
            Real Field with 53 bits of precision
        """
        if not data:
            if coxeter_type:
                data = CoxeterType(coxeter_type)
            elif cartan_type:
                data = CoxeterType(CartanType(cartan_type))

        # Special cases with no arguments passed
        if not data:
            graph = Graph()
            index_set = tuple()
            coxeter_type = None
            graph = typecall(cls, graph, coxeter_type, index_set)

            return graph

        if isinstance(data, CoxeterGraph):  # Initiate from itself
            return data

        # Initiate from a graph:
        if isinstance(data, Graph):
            return cls._from_graph(data, coxeter_type, coxeter_type_check)

        # Get the Coxeter type
        coxeter_type = None
        from sage.combinat.root_system.cartan_type import CartanType_abstract
        if isinstance(data, CartanType_abstract):
            coxeter_type = data.coxeter_type()
        else:
            try:
                coxeter_type = CoxeterType(data)
            except (TypeError, ValueError, NotImplementedError):
                pass

        # Initiate from a Coxeter type
        if coxeter_type:
            return cls._from_coxetertype(coxeter_type)

        # TODO:: remove when oo is possible in matrices.
        n = len(data[0])
        data = [x if x != infinity else -1 for r in data for x in r]
        data = matrix(n, n, data)
        # until here

        # Get the index set
        if index_set:
            index_set = tuple(index_set)
        else:
            index_set = tuple(range(1,n+1))
        if len(set(index_set)) != n:
                raise ValueError("the given index set is not valid")

        return cls._from_matrix(data, index_set, coxeter_type_check)

    def __init__(self, data, coxeter_type, index_set):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: d = CoxeterGraph(["A", 3])
            sage: TestSuite(d).run()
        """
        
        self._graph = Graph(data).copy(immutable=True)      
        self._coxeter_type = coxeter_type
        self._index_set = index_set
        self._rank = self._graph.num_verts()

        if self._coxeter_type is not None:
            if self._coxeter_type.is_finite():
                self._is_finite = True
                self._is_affine = False
            elif self._coxeter_type.is_affine():
                self._is_finite = False
                self._is_affine = True
            else:
                self._is_finite = False
                self._is_affine = False
        else:
            self._is_finite = False
            self._is_affine = False

    @classmethod
    def _from_graph(cls, graph, coxeter_type, coxeter_type_check):
        """
        Initiate the Coxeter graph from a graph.

        TESTS::

            sage: CoxeterMatrix(CoxeterMatrix(['A',4,1]).coxeter_graph())
            [1 3 2 2 3]
            [3 1 3 2 2]
            [2 3 1 3 2]
            [2 2 3 1 3]
            [3 2 2 3 1]
            sage: CoxeterMatrix(CoxeterMatrix(['B',4,1]).coxeter_graph())
            [1 2 3 2 2]
            [2 1 3 2 2]
            [3 3 1 3 2]
            [2 2 3 1 4]
            [2 2 2 4 1]
            sage: CoxeterMatrix(CoxeterMatrix(['F',4]).coxeter_graph())
            [1 3 2 2]
            [3 1 4 2]
            [2 4 1 3]
            [2 2 3 1]

            sage: G=Graph()
            sage: G.add_edge([0,1,oo])
            sage: CoxeterMatrix(G)
            [ 1 -1]
            [-1  1]
            sage: H = Graph()
            sage: H.add_edge([0,1,-1.5])
            sage: CoxeterMatrix(H)
            [ 1.00000000000000 -1.50000000000000]
            [-1.50000000000000  1.00000000000000]
        """
        check_coxeter_graph(graph)
        
        verts = sorted(graph.vertices())
        index_set = tuple(verts)
        n = len(index_set)
        
        if not coxeter_type:
            if n == 1:
                coxeter_type = CoxeterType(['A', 1])
            elif coxeter_type_check:
                coxeter_type = recognize_coxeter_type_from_graph(graph, index_set)
            else:
                coxeter_type = None
        graph = typecall(cls, graph, coxeter_type, index_set)

        return graph

    @classmethod
    def _from_matrix(cls, data, index_set, coxeter_type_check):
        """
        Initiate the Coxeter graph from a matrix.

        TESTS::

            sage: CM = CoxeterMatrix([[1,2],[2,1]]); CM
            [1 2]
            [2 1]
            sage: CM = CoxeterMatrix([[1,-1],[-1,1]]); CM
            [ 1 -1]
            [-1  1]
            sage: CM = CoxeterMatrix([[1,-1.5],[-1.5,1]]); CM
            [ 1.00000000000000 -1.50000000000000]
            [-1.50000000000000  1.00000000000000]
            sage: CM = CoxeterMatrix([[1,-3/2],[-3/2,1]]); CM
            [   1 -3/2]
            [-3/2    1]
            sage: CM = CoxeterMatrix([[1,-3/2,5],[-3/2,1,-1],[5,-1,1]]); CM
            [   1 -3/2    5]
            [-3/2    1   -1]
            [   5   -1    1]
            sage: CM = CoxeterMatrix([[1,-3/2,5],[-3/2,1,oo],[5,oo,1]]); CM
            [   1 -3/2    5]
            [-3/2    1   -1]
            [   5   -1    1]
        """
        # Check that the data is valid
        from sage.combinat.root_system.coxeter_matrix import check_coxeter_matrix
        check_coxeter_matrix(data)

        M = matrix(data)
        n = M.ncols()
        
        graph = Graph(multiedges=False)
        graph.add_vertices(index_set)
        for i in range(n):
            for j in range(i+1,n):
                label = M[i,j]
                if label != 2:
                    graph.add_edge(index_set[i],index_set[j],M[i,j])

        return cls._from_graph(graph, None, index_set, coxeter_type_check)

    @classmethod
    def _from_coxetertype(cls, coxeter_type):
        """
        Initiate the Coxeter graph from a Coxeter type.

        TESTS::

            sage: CoxeterMatrix(['A',4]).coxeter_type()
            Coxeter type of ['A', 4]
            sage: CoxeterMatrix(['A',4,1]).coxeter_type()
            Coxeter type of ['A', 4, 1]
            sage: CoxeterMatrix(['D',4,1]).coxeter_type()
            Coxeter type of ['D', 4, 1]
        """
        index_set = coxeter_type.index_set()
        from sage.rings.infinity import infinity # necessary???
        scalarproducts_to_order = { 0: 2,  1: 3,  2: 4,  3: 6, 4: infinity }
        graph = Graph(multiedges=False)
        a = coxeter_type._cartan_type.dynkin_diagram()
        graph.add_vertices(index_set)
        for i in index_set:
            for j in a.neighbors_out(i):
                # avoid adding the edge twice
                if not graph.has_edge(i,j):
                    graph.add_edge(i,j, scalarproducts_to_order[a[i,j]*a[j,i]])            

        return cls._from_graph(graph, coxeter_type, False)

    def _repr_(self):
        """
        EXAMPLES::

            sage: CoxeterGraph(['G',2])     # indirect doctest
              3
            O=<=O
            1   2
            G2
        """
        ct = self.coxeter_type()
        result = ct.ascii_art() +"\n" if hasattr(ct, "ascii_art") else ""

        if ct is None or isinstance(ct, CoxeterMatrix) or isinstance(ct, CoxeterGraph):
            return result+"Coxeter graph of rank %s"%self.rank()
        else:
            return result+"%s"%ct._repr_(compact=True)
            #return result+"Coxeter graph of type %s"%self.coxeter_type()._repr_(compact = True)

    def _latex_(self, scale=0.5):
        r"""
        Return a latex representation of this Coxeter graph

        EXAMPLES::

            sage: latex(CoxeterGraph(['A',3,1]))
            \begin{tikzpicture}[scale=0.5]
            \draw (-1,0) node[anchor=east] {$A_{3}^{(1)}$};
            \draw (0 cm,0) -- (4 cm,0);
            \draw (0 cm,0) -- (2.0 cm, 1.2 cm);
            \draw (2.0 cm, 1.2 cm) -- (4 cm, 0);
            \draw[fill=white] (0 cm, 0) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (2.0 cm, 1.2 cm) circle (.25cm) node[anchor=south east]{$0$};
            \end{tikzpicture}
        """
        if self.coxeter_type() is None:
            return "Coxeter graph of rank %s"%self.rank()
        ret = "\\begin{tikzpicture}[scale=%s]\n"%scale
        ret += "\\draw (-1,0) node[anchor=east] {$%s$};\n"%self.coxeter_type()._latex_()
        ret += self.coxeter_type()._latex_coxeter_graph()
        ret += "\n\\end{tikzpicture}"
        return ret

    def _matrix_(self):
        """
        Return a regular matrix from ``self``.

        EXAMPLES::

            sage: M = CoxeterGraph(['C',3])._matrix_(); M
            [ 2 -1  0]
            [-1  2 -2]
            [ 0 -1  2]
            sage: type(M)
            <class 'sage.combinat.root_system.cartan_matrix.CoxeterMatrix'>
        """
        return self.coxeter_matrix()._matrix_()

    def add_edge(self, i, j, label=1):
        """
        EXAMPLES::

            sage: from sage.combinat.root_system.coxeter_graph import CoxeterGraph_class
            sage: d = CoxeterGraph_class(CoxeterType(['A',3]))
            sage: list(sorted(d.edges()))
            []
            sage: d.add_edge(2, 3)
            sage: list(sorted(d.edges()))
            [(2, 3, 1), (3, 2, 1)]
        """
        Graph.add_edge(self, i, j, label)
        if not self.has_edge(j,i):
            self.add_edge(j,i,1)
    
    def edge_iterator(self):
        """
        EXAMPLES::
        
        
        """
        return self._graph.edge_iterator()

    def __hash__(self):
        """
        EXAMPLES::

            sage: d = CoxeterType(['A',3]).coxeter_graph()
            sage: hash(d) == hash((d.coxeter_type(), tuple(d.vertices()), tuple(d.edge_iterator(d.vertices()))))
            True
        """
        # Should assert for immutability!

        #return hash(self.coxeter_type(), self.vertices(), tuple(self.edges()))
        # FIXME: self.edges() currently tests at some point whether
        # self is a vertex of itself which causes an infinite
        # recursion loop. Current workaround: call self.edge_iterator directly
        return hash((self.coxeter_type(), tuple(self.vertices()), tuple(self.edge_iterator(self.vertices()))))

    @staticmethod
    def an_instance():
        """
        Returns an example of Coxeter graph

        EXAMPLES::

            sage: from sage.combinat.root_system.coxeter_graph import CoxeterGraph_class
            sage: g = CoxeterGraph_class.an_instance()
            sage: g
            Coxeter graph of rank 3
            sage: g.coxeter_matrix()
            [ 2 -1 -1]
            [-2  2 -1]
            [-1 -1  2]

        """
        # hyperbolic Dynkin diagram of Exercise 4.9 p. 57 of Kac Infinite Dimensional Lie Algebras.
        g = CoxeterGraph()
        g.add_vertices([1,2,3])
        g.add_edge(1,2,2)
        g.add_edge(1,3)
        g.add_edge(2,3)
        return g

    ##########################################################################
    # Coxeter type methods

    @cached_method
    def index_set(self):
        """
        EXAMPLES::

            sage: CoxeterGraph(['C',3]).index_set()
            (1, 2, 3)
            sage: CoxeterGraph("A2","B2","F4").index_set()
            (1, 2, 3, 4, 5, 6, 7, 8)
        """
        return tuple(self.vertices())

    def coxeter_type(self):
        """
        EXAMPLES::

            sage: CoxeterGraph("A2","B2","F4").coxeter_type()
            A2xB2xF4
        """
        return self._coxeter_type

    def rank(self):
        r"""
        Returns the index set for this Coxeter graph

        EXAMPLES::

            sage: CoxeterGraph(['C',3]).rank()
            3
            sage: CoxeterGraph("A2","B2","F4").rank()
            8
        """
        return self.num_verts()

    def coxeter_graph(self):
        """
        EXAMPLES::

            sage: CoxeterGraph(['C',3]).coxeter_graph()
            O---O=<=O
            1   2   3
            C3
        """
        return self

    @cached_method
    def coxeter_matrix(self):
        r"""
        Returns the Coxeter matrix for this Coxeter graph

        EXAMPLES::

            sage: CoxeterGraph(['C',3]).coxeter_matrix()
            [ 2 -1  0]
            [-1  2 -2]
            [ 0 -1  2]
        """
        return CoxeterMatrix(self)

    def is_finite(self):
        """
        Check if ``self`` corresponds to a finite Coxeter group.

        EXAMPLES::

            sage: CoxeterGroup(['F',4]).coxeter_graph().is_finite()
            True
            sage: D = CoxeterGroup(CoxeterMatrix([[1, -2], [-2, 1]]))
            sage: D.is_finite()
            False
        """
        if self._coxeter_type is not None:
            return self._coxeter_type.is_finite()
        return self.coxeter_matrix().is_finite()

    def is_affine(self):
        """
        Check if ``self`` corresponds to an affine Coxeter group.

        EXAMPLES::

            sage: CoxeterType(['F',4]).coxeter_graph().is_affine()
            False
            sage: D = CoxeterGraph(CoxeterMatrix([[1, -2], [-2, 1]]))
            sage: D.is_affine()
            False
        """
        if self._coxeter_type is not None:
            return self._coxeter_type.is_affine()
        return self.coxeter_matrix().is_affine()

    def is_irreducible(self):
        """
        Check if ``self`` corresponds to an irreducible Coxeter group.

        EXAMPLES::

            sage: CoxeterGroup(['F',4]).coxeter_graph().is_irreducible()
            True
        """
        return self._coxeter_type.is_irreducible()

    def is_crystallographic(self):
        """
        Implements :meth:`CoxeterType.is_crystallographic`???

        Checks if ``self`` corresponds to a crystallographic Coxeter group.

        EXAMPLES::

            sage: CoxeterGroup(['F',4]).coxeter_graph().is_crystallographic()
            True
            sage: CoxeterGroup(['H',4]).coxeter_graph().is_crystallographic()
            False

        """
        return True

    def __getitem__(self, i):
        r"""
        With a tuple (i,j) as argument, returns the scalar product
        `\langle
                \alpha^\vee_i, \alpha_j\rangle`.

        Otherwise, behaves as the usual DiGraph.__getitem__

        EXAMPLES: We use the `C_4` Coxeter graph as a coxeter
        matrix::

            sage: g = CoxeterGraph(['C',4])
            sage: matrix([[g[i,j] for j in range(1,5)] for i in range(1,5)])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -2]
            [ 0  0 -1  2]

        The neighbors of a node can still be obtained in the usual way::

            sage: [g[i] for i in range(1,5)]
            [[2], [1, 3], [2, 4], [3]]
        """
        if not isinstance(i, tuple):
            return Graph.__getitem__(self,i)
        [i,j] = i
        if i == j:
            return 2
        elif self.has_edge(j, i):
            return -self.edge_label(j, i)
        else:
            return 0

    def column(self, j):
        """
        Returns the `j^{th}` column `(a_{i,j})_i` of the
        Coxeter matrix corresponding to this Coxeter graph, as a container
        (or iterator) of tuples `(i, a_{i,j})`

        EXAMPLES::

            sage: g = Coxetergraph(["B",4])
            sage: [ (i,a) for (i,a) in g.column(3) ]
            [(3, 2), (2, -1), (4, -2)]
        """
        return [(j,2)] + [(i,-m) for (j1, i, m) in self.outgoing_edges(j)]

    def row(self, i):
        """
        Returns the `i^{th}` row `(a_{i,j})_j` of the
        Coxeter matrix corresponding to this Coxeter graph, as a container
        (or iterator) of tuples `(j, a_{i,j})`

        EXAMPLES::

            sage: g = CoxeterGraph(["C",4])
            sage: [ (i,a) for (i,a) in g.row(3) ]
            [(3, 2), (2, -1), (4, -2)]
        """
        return [(i,2)] + [(j,-m) for (j, i1, m) in self.incoming_edges(i)]

#####################################################################
## Other functions

def recognize_coxeter_type_from_graph(graph):
    """
    Return the Coxeter type of ``coxeter_graph`` if known,
    otherwise return ``None``.

    EXAMPLES:

    Some infinite ones::

        sage: C = CoxeterMatrix([[1,3,2],[3,1,-1],[2,-1,1]])
        sage: C.is_finite()  # indirect doctest
        False
        sage: C = CoxeterMatrix([[1,-1,-1],[-1,1,-1],[-1,-1,1]])
        sage: C.is_finite()  # indirect doctest
        False

    Some finite ones::

        sage: m = matrix(CoxeterMatrix(['D', 4]))
        sage: CoxeterMatrix(m).is_finite()  # indirect doctest
        True
        sage: m = matrix(CoxeterMatrix(['H', 4]))
        sage: CoxeterMatrix(m).is_finite()  # indirect doctest
        True

        sage: CoxeterMatrix(CoxeterType(['A',10]).coxeter_graph()).coxeter_type()
        Coxeter type of ['A', 10]
        sage: CoxeterMatrix(CoxeterType(['B',10]).coxeter_graph()).coxeter_type()
        Coxeter type of ['B', 10]
        sage: CoxeterMatrix(CoxeterType(['C',10]).coxeter_graph()).coxeter_type()
        Coxeter type of ['B', 10]
        sage: CoxeterMatrix(CoxeterType(['D',10]).coxeter_graph()).coxeter_type()
        Coxeter type of ['D', 10]
        sage: CoxeterMatrix(CoxeterType(['E',6]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 6]
        sage: CoxeterMatrix(CoxeterType(['E',7]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 7]
        sage: CoxeterMatrix(CoxeterType(['E',8]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 8]
        sage: CoxeterMatrix(CoxeterType(['F',4]).coxeter_graph()).coxeter_type()
        Coxeter type of ['F', 4]
        sage: CoxeterMatrix(CoxeterType(['G',2]).coxeter_graph()).coxeter_type()
        Coxeter type of ['G', 2]
        sage: CoxeterMatrix(CoxeterType(['H',3]).coxeter_graph()).coxeter_type()
        Coxeter type of ['H', 3]
        sage: CoxeterMatrix(CoxeterType(['H',4]).coxeter_graph()).coxeter_type()
        Coxeter type of ['H', 4]
        sage: CoxeterMatrix(CoxeterType(['I',100]).coxeter_graph()).coxeter_type()
        Coxeter type of ['I', 100]

    Some affine graphs::

        sage: CoxeterMatrix(CoxeterType(['A',1,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['A', 1, 1]
        sage: CoxeterMatrix(CoxeterType(['A',10,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['A', 10, 1]
        sage: CoxeterMatrix(CoxeterType(['B',10,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['B', 10, 1]
        sage: CoxeterMatrix(CoxeterType(['C',10,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['C', 10, 1]
        sage: CoxeterMatrix(CoxeterType(['D',10,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['D', 10, 1]
        sage: CoxeterMatrix(CoxeterType(['E',6,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 6, 1]
        sage: CoxeterMatrix(CoxeterType(['E',7,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 7, 1]
        sage: CoxeterMatrix(CoxeterType(['E',8,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 8, 1]
        sage: CoxeterMatrix(CoxeterType(['F',4,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['F', 4, 1]
        sage: CoxeterMatrix(CoxeterType(['G',2,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['G', 2, 1]

    TESTS:

    Check that we detect relabellings::

        sage: M = CoxeterMatrix([[1,2,3],[2,1,6],[3,6,1]], index_set=['a', 'b', 'c'])
        sage: M.coxeter_type()
        Coxeter type of ['G', 2, 1] relabelled by {0: 'a', 1: 'b', 2: 'c'}

        sage: from sage.combinat.root_system.coxeter_matrix import recognize_coxeter_type_from_matrix
        sage: for C in CoxeterMatrix.samples():
        ....:     relabelling_perm = Permutations(C.index_set()).random_element()
        ....:     relabelling_dict = {C.index_set()[i]: relabelling_perm[i] for i in range(C.rank())}
        ....:     relabeled_matrix = C.relabel(relabelling_dict)._matrix
        ....:     recognized_type = recognize_coxeter_type_from_matrix(relabeled_matrix, relabelling_perm)
        ....:     if C.is_finite() or C.is_affine():
        ....:         assert recognized_type == C.coxeter_type()
    """

    types = []
    for S in graph.connected_components_subgraphs():
        r = S.num_verts()
        # Handle the special cases first
        if r == 1:
            types.append(CoxeterType(['A',1]).relabel({1: S.vertices()[0]}))
            continue
        if r == 2: # Type B2, G2, or I_2(p)
            e = S.edge_labels()[0]
            if e == 3: # Can't be 2 because it is connected
                ct = CoxeterType(['B',2])
            elif e == 4:
                ct = CoxeterType(['G',2])
            elif e > 0 and e < float('inf'): # Remaining non-affine types
                ct = CoxeterType(['I',e])
            else: # Otherwise it is infinite dihedral group Z_2 \ast Z_2
                ct = CoxeterType(['A',1,1])
            if not ct.is_affine():
                types.append(ct.relabel({1: S.vertices()[0], 2: S.vertices()[1]}))
            else:
                types.append(ct.relabel({0: S.vertices()[0], 1: S.vertices()[1]}))
            continue

        test = [['A',r], ['B',r], ['A',r-1,1]]
        if r >= 3:
            if r == 3:
                test += [['G',2,1], ['H',3]]
            test.append(['C',r-1,1])
        if r >= 4:
            if r == 4:
                test += [['F',4], ['H',4]]
            test += [['D',r], ['B',r-1,1]]
        if r >= 5:
            if r == 5:
                test.append(['F',4,1])
            test.append(['D',r-1,1])
        if r == 6:
            test.append(['E',6])
        elif r == 7:
            test += [['E',7], ['E',6,1]]
        elif r == 8:
            test += [['E',8], ['E',7,1]]
        elif r == 9:
            test.append(['E',8,1])

        found = False
        for ct in test:
            ct = CoxeterType(ct)
            print ct
            T = ct.coxeter_graph()._graph
            iso, match = T.is_isomorphic(S, certify=True, edge_labels=True)
            if iso:
                types.append(ct.relabel(match))
                found = True
                break
        if not found:
            return None

    return CoxeterType(types)

def check_coxeter_graph(graph):
    """
    Check if ``graph`` represents a generalized Coxeter graph and raise
    and error if not.

    EXAMPLES::

    """
    for e in graph.edges():
        label = e[2]
        if label is not None:
            if label not in ZZ and label > -1:
                raise ValueError("invalid Coxeter graph label")
            elif label == 0 or label == 1:
                raise ValueError("invalid Coxeter graph label")

def precheck(t, letter=None, length=None, affine=None, n_ge=None, n=None):
    """
    EXAMPLES::

        sage: from sage.combinat.root_system.coxeter_graph import precheck
        sage: ct = CoxeterType(['A',4])
        sage: precheck(ct, letter='C')
        Traceback (most recent call last):
        ...
        ValueError: t[0] must be = 'C'
        sage: precheck(ct, affine=1)
        Traceback (most recent call last):
        ...
        ValueError: t[2] must be = 1
        sage: precheck(ct, length=3)
        Traceback (most recent call last):
        ...
        ValueError: len(t) must be = 3
        sage: precheck(ct, n=3)
        Traceback (most recent call last):
        ...
        ValueError: t[1] must be = 3
        sage: precheck(ct, n_ge=5)
        Traceback (most recent call last):
        ...
        ValueError: t[1] must be >= 5
    """
    if letter is not None:
        if t[0] != letter:
            raise ValueError("t[0] must be = '%s'"%letter)

    if length is not None:
        if len(t) != length:
            raise ValueError("len(t) must be = %s"%length)

    if affine is not None:
        try:
            if t[2] != affine:
                raise ValueError("t[2] must be = %s"%affine)
        except IndexError:
            raise ValueError("t[2] must be = %s"%affine)

    if n_ge is not None:
        if t[1] < n_ge:
            raise ValueError("t[1] must be >= %s"%n_ge)

    if n is not None:
        if t[1] != n:
            raise ValueError("t[1] must be = %s"%n)
