"""
Coxeter graphs

AUTHORS:

- Jean-Philippe Labbe (2014-10-22): Created the class based on Dynkin diagrams.

"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#       Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#       Copyright (C) 2014 Jean-Philippe Labbe <labbe@math.fu-berlin.de>
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
from sage.matrix.matrix import is_Matrix
from sage.graphs.graph import Graph
from sage.combinat.root_system.coxeter_type import CoxeterType 
from sage.combinat.root_system.coxeter_matrix import CoxeterMatrix

def CoxeterGraph(*args, **kwds):
    r"""
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
        sage: CM = R.cartan_matrix(); CM
        [ 2 -1| 0  0| 0  0  0  0]
        [-1  2| 0  0| 0  0  0  0]
        [-----+-----+-----------]
        [ 0  0| 2 -1| 0  0  0  0]
        [ 0  0|-2  2| 0  0  0  0]
        [-----+-----+-----------]
        [ 0  0| 0  0| 2 -1  0  0]
        [ 0  0| 0  0|-1  2 -1  0]
        [ 0  0| 0  0| 0 -2  2 -1]
        [ 0  0| 0  0| 0  0 -1  2]
        sage: DD = CoxeterGraph(CM); DD
        O---O
        1   2
        O=>=O
        3   4
        O---O=>=O---O
        5   6   7   8
        A2xB2xF4
        sage: DD.cartan_matrix()
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
    """
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


class CoxeterGraph_class(Graph, CoxeterType):
    """
    A generalized Coxeter graph.

    .. SEEALSO::

        :func:`CoxeterGraph()`

    INPUT:

    - ``t`` -- a Coxeter type, a generalized Coxeter matrix, or ``None``

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

    Implementation note: if a Coxeter type is given, then the nodes
    are initialized from the index set of this Coxeter type.
    """
    def __init__(self, t=None, index_set=None, **options):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: d = CoxeterGraph(["A", 3])
            sage: TestSuite(d).run()
        """
        if isinstance(t, Graph):
            if isinstance(t, CoxeterGraph_class):
                self._coxeter_type = t._coxeter_type
            else:
                self._coxeter_type = None
            Graph.__init__(self, data=t, **options)
            return

        Graph.__init__(self, **options)
        self._coxeter_type = t
        if index_set is not None:
            self.add_vertices(index_set)
        elif t is not None:
            self.add_vertices(t.index_set())

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

        if ct is None or isinstance(ct, CoxeterMatrix):
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
