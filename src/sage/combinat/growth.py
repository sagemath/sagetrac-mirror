r"""
Growth diagrams and dual graded graphs

AUTHORS:

- Martin Rubey (2016-09): Initial version

.. TODO::

    - when shape is given, check that it is compatible with filling or labels
    - implement backward rules for :class:`GrowthDiagramDomino`
    - optimise rules, mainly for :class:`GrowthDiagramRSK` and :class:`GrowthDiagramBurge`
    - make semistandard extension generic
    - implement rules from [LamShi2007]_

Growth diagrams, invented by Sergey Fomin [Fom1995]_, provide a vast
generalisation of the Robinson-Schensted-Knuth correspondence between
matrices with non-negative integer entries and pairs of semistandard
Young tableaux of the same shape.

The main fact is that many correspondences similar to RSK can be
defined by providing a so-called 'forward' rule: a function whose
input are three vertices x, y and t of a certain directed graph (in
the case of Robinson-Schensted: the directed graph corresponding to
Young's lattice) and an integer (in the case of Robinson-Schensted:
zero or one), and whose output is a fourth vertex z.  This rule
should be invertible in the following sense: there is a so-called
'backward' rule that recovers the integer and t given z, x and y.

The classical Robinson-Schensted-Knuth correspondence is provided by
:class:`GrowthDiagramRSK`.  Note that a growth diagram is printed
with matrix coordinates, the origin being in the top-left corner::

    sage: w = [2,3,6,1,4,5]; G = GrowthDiagramRSK(w); G
    0  0  0  1  0  0
    1  0  0  0  0  0
    0  1  0  0  0  0
    0  0  0  0  1  0
    0  0  0  0  0  1
    0  0  1  0  0  0

The 'forward' rule just mentioned assigns partitions to the corners
of each of the 36 cells of this matrix - with the exception of the
corners on the left and top boundary, which are initialized with the
empty partition.  The partitions along the boundary opposite of the
origin are obtained by using the method
:meth:`GrowthDiagramRSK.out_labels()`::

    sage: G.out_labels()
    [[],
     [1],
     [2],
     [3],
     [3, 1],
     [3, 2],
     [4, 2],
     [4, 1],
     [3, 1],
     [2, 1],
     [1, 1],
     [1],
     []]

However, in the case of a rectangular filling, it is more practical
to split this sequence of labels in two.  We then obtain the `P` and
`Q` symbols::

    sage: [G.P_symbol(), G.Q_symbol()]
    [[[1, 3, 4, 5], [2, 6]], [[1, 2, 3, 6], [4, 5]]]
    sage: RSK(w)
    [[[1, 3, 4, 5], [2, 6]], [[1, 2, 3, 6], [4, 5]]]

A great advantage of growth diagrams is that we immediately have
access also to the skew version of the RSK-correspondence, by
providing different initialisation for the labels near the origin.
We reproduce the original example of Bruce Sagan and Richard Stanley,
see also Tom Roby's thesis [Rob1991]_.  We can represent the
generalised permutation::

    1 2 4
    4 2 3

as a dictionary of coordinates, subtracting `1` from all entries
because lists in SageMath are zero-based::

    sage: w = {(1-1,4-1):1, (2-1,2-1):1, (4-1,3-1):1}
    sage: T = SkewTableau([[None, None],[None,5],[1]])
    sage: U = SkewTableau([[None, None],[None,3],[5]])
    sage: G = GrowthDiagramRSK(filling = w, shape = [5]*5, labels = T.to_chain()[::-1]+U.to_chain()[1:]); G
    0  0  0  0  0
    0  1  0  0  0
    0  0  0  1  0
    1  0  0  0  0
    0  0  0  0  0
    sage: G.P_symbol(), G.Q_symbol()
    ([[None, None, 2, 3], [None, None], [None, 4], [1], [5]],
     [[None, None, 1, 4], [None, None], [None, 2], [3], [5]])

Moreover, we are not forced to use rectangular fillings.  For
example, consider the Stanley-Sundaram correspondence between (skew)
oscillating tableaux and (partial) perfect matchings.  Again, from
Tom Roby's thesis::

    sage: o = [[2,1],[2,2],[3,2],[4,2],[4,1],[4,1,1],[3,1,1],[3,1],[3,2],[3,1],[2,1]]
    sage: l = [None]*(2*len(o)-1)
    sage: l[::2] = [Partition(la) for la in o]
    sage: l[1::2] = [l[2*i] if l[2*i].size() < l[2*i+2].size() else l[2*i+2] for i in range(len(o)-1)]
    sage: G = GrowthDiagramRSK(labels=l[1:-1], shape=[i for i in range(len(o)-2,0,-1)]); G
    0  0  0  0  0  0  0  1  0
    0  1  0  0  0  0  0  0
    0  0  0  0  0  0  0
    0  0  0  0  0  0
    0  0  1  0  0
    0  0  0  0
    0  0  0
    0  0
    0
    sage: ascii_art(SkewTableau(chain=G.in_labels()[len(o)-2:]), SkewTableau(chain=G.in_labels()[:len(o)-1][::-1]))
    .  1  .  7
    5     4

As mentioned at the beginning, the Robinson-Schensted-Knuth
correspondence is just a special case of growth diagrams.  In
particular, we have implemented local rules for the variation of RSK
originally due to Burge (:class:`GrowthDiagramBurge`), a
correspondence producing binary words originally due to Viennot
(:class:`GrowthDiagramBinWord`), and a correspondence producing
domino tableaux (:class:`GrowthDiagramDomino`) originally due to
Barbasch and Vogan.

REFERENCES:

.. [Fom1995] Sergey V. Fomin.
   *Schensted algorithms for dual graded graphs*.
   Journal of Algebraic Combinatorics Volume 4, Number 1 (1995), pp. 5-45

.. [Rob1991] Tom Roby.
   *Applications and extensions of Fomin's generalization of the Robinson-Schensted correspondence to differential posets*.
   M.I.T., Cambridge, Massachusetts

.. [Kra2006] Christian Krattenthaler.
   *Growth diagrams, and increasing and decreasing chains in fillings of Ferrers shapes*.
   Advances in Applied Mathematics Volume 37, Number 3 (2006), pp. 404-431

.. [Lam2004] Thomas Lam.
   *Growth diagrams, domino insertion and sign-imbalance*.
   Journal of Combinatorial Theory, Series A Volume 107, Number 1 (2004), pp. 87-115

.. [LamShi2007] Thomas Lam and Mark Shimozono.
   *Dual graded graphs for Kac-Moody algebras*.
   Algebra & Number Theory 1.4 (2007): pp. 451-488

.. [LLMSSZ2013] Thomas Lam, Luc Lapointe, Jennifer Morse, Anne Schilling, Mark Shimozono and Mike Zabrocki.
   *k-Schur functions and affine Schubert calculus*.
   https://arxiv.org/pdf/1301.3569.pdf

.. [HivNze] Florent Hivert and Janvier Nzeutchap.
   *Dual Graded Graphs in Combinatorial Hopf Algebras*.
   https://www.lri.fr/~hivert/PAPER/commCombHopfAlg.pdf

.. [vanLee1996] Marc van Leeuwen.
   *The Robinson-Schensted and Schützenberger algorithms, an elementary approach*.
   Electronic Journal of Combinatorics 3, no. 2 (1996): Research Paper 15, approx. 32 pp. (electronic)

.. [Vie1983] Xavier G. Viennot.
   *Maximal chains of subwords and up-down sequences of permutations*.
   Journal of Combinatorial Theory, Series A Volume 34, (1983), pp. 1-14

.. [Nze2007] Janvier Nzeutchap.
   *Binary Search Tree insertion, the Hypoplactic insertion, and Dual Graded Graphs*.
   :arXiv:`0705.2689` (2007)

.. [Sag1987] Bruce E. Sagan.
   *Shifted tableaux, Schur Q-functions, and a conjecture of R. Stanley*.
   Journal of Combinatorial Theory, Series A Volume 45 (1987), pp. 62-103

"""
from sage.structure.sage_object import SageObject
from sage.combinat.posets.posets import Poset
from sage.combinat.words.word import Word
from sage.combinat.words.words import Words
from sage.combinat.binary_tree import BinaryTree
from sage.combinat.binary_tree import BinaryTrees
from sage.combinat.composition import Compositions
from sage.combinat.partition import Partition, Partitions
from sage.combinat.skew_partition import SkewPartition
from sage.combinat.skew_tableau import SkewTableau
from sage.combinat.core import Core, Cores
from sage.combinat.k_tableau import WeakTableau, StrongTableau
from copy import copy
from sage.misc.functional import is_odd, is_even
from sage.rings.integer_ring import ZZ
from sage.graphs.digraph import DiGraph

class GrowthDiagram(SageObject):
    r"""
    The base class all variants of growth diagrams inherit from.

    Inheriting classes should provide an ``__init__`` method that
    checks that ``labels``, when provided, are of the correct
    type, and sets the following attributes:

    - ``self._zero``, the zero element of the vertices of the
      graphs.

    - ``self._rank_function``, the rank function of the dual
      graded graphs.

    - ``self._has_multiple_edges``, if the dual graded graph has
      multiple edges and therefore edges are triples consisting of
      two vertices and a label.  If not present, it is assumed to be
      ``False``.

    - ``self._zero_edge``, the zero label of the edges of the graphs.
      If not present, it is assumed to be the integer 0.

    - ``self._is_Q_edge``, ``self._is_P_edge``, functions that take
      two vertices as arguments and return ``True`` or ``False``, or,
      if multiple edges are allowed, the list of edge labels of the
      edges from the second to the first in the respective graded
      graph.  These are only used for checking user input and
      providing the dual graded graph, and are therefore not
      mandatory.

    It should then call the ``__init__`` method of this class.

    EXAMPLES::

        sage: w = [3,3,2,4,1]; G = GrowthDiagramRSK(w)
        sage: [G.P_symbol(), G.Q_symbol()]
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
        sage: RSK(w)
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]

    TESTS::

        sage: G = GrowthDiagramRSK()
        Traceback (most recent call last):
        ...
        ValueError: Please provide a filling or a sequence of labels.
    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        r"""
        Initialise a generalized Schensted growth diagram in the sense of
        Fomin.

        An instance of the class is a growth diagram consisting of a
        filling and labels on the boundary.  To initialise it, it is
        necessary to provide either a filling and labels for the
        border on the side of the origin (forward growth), or labels
        for the boundary opposite of the origin (backward growth).

        Coordinates are of the form (col, row) where the origin is in
        the upper left, to be consistent with permutation matrices
        and skew tableaux.  This is different from Fomin's
        convention, who uses a Cartesian coordinate system.

        INPUT:

        - ``filling`` is ``None``, if the growth diagram should be
          determined by applying the backward rules to the ``labels``
          decorating the boundary opposite of the origin of the
          ``shape``.  Otherwise it can be a dictionary with keys
          being coordinates and integer values, a sequence of
          sequences of integers (including matrices), or a word with
          integer letters (including permutations).  In the latter
          case, words with negative letters but without repetitions
          are allowed and interpreted as coloured permutations.

        - ``shape`` is a (possibly skew) partition or ``None``.  In
          the latter case it is determined as the Ferrers shape given
          by ``filling``, if that is a sequence of sequences, the
          bounding rectangle (including the origin) of ``filling``,
          or, if possible, inferred from ``labels``, if ``filling``
          is ``None``.

        - ``labels`` is ``None`` or a list.  If it is a list, it
          specifies a path whose length is the half-perimeter of
          ``shape``.  Thus, when ``self._has_multiple_edges`` is
          ``True``, it is of the form `(v_1, e_1,..., e_{n-1}, v_n)`,
          where `n` is the half-perimeter of ``shape``, and
          `(v_{i-1}, e_i, v_i)` is an edge in the dual graded graph
          for all `i`.  Otherwise, it is a list of `n` vertices.

          If ``filling`` is ``None``, its elements are the labels of
          the boundary opposite of the origin.  Otherwise its
          elements are the labels on the boundary on the side of the
          origin.  If ``labels`` is ``None`` (in which case
          ``filling`` must not be ``None``) the value of
          ``self._zero`` is used to initialise ``labels``.  """
        try:
            self._has_multiple_edges
        except AttributeError:
            self._has_multiple_edges = False
        try:
            self._zero_edge
        except AttributeError:
            self._zero_edge = 0
        try:
            self._is_Q_edge
        except AttributeError:
            if self._has_multiple_edges:
                self._is_Q_edge = lambda a,b: [self._zero_edge]
            else:
                self._is_Q_edge = lambda a,b: True
        try:
            self._is_P_edge
        except AttributeError:
            if self._has_multiple_edges:
                self._is_P_edge = lambda a,b: [self._zero_edge]
            else:
                self._is_P_edge = lambda a,b: True

        if filling is None:
            if labels is None:
                raise ValueError("Please provide a filling or a sequence of labels.")

            if shape is None:
                shape = self._shape_from_labels(labels)

            self._lambda, self._mu = self._init_shape_from_various_input(shape)
            self._out_labels = labels
            self._check_labels(self._out_labels)
            self._shrink()
        else:
            self._filling, (self._lambda, self._mu) = self._init_filling_and_shape_from_various_input(filling, shape)
            self._in_labels = self._init_labels_forward_from_various_input(labels)
            self._check_labels(self._in_labels)
            self._grow()

    @classmethod
    def _check_duality(cls, n, r=1):
        """
        Raise an error if the graphs are not r-dual at level n.
        """
        try:
            cls._has_multiple_edges
        except AttributeError:
            cls._has_multiple_edges = False
        if cls._has_multiple_edges:
            def check_vertex(w, P, Q):
                DUw = [v[0] for uw in P.outgoing_edges(w) for v in Q.incoming_edges(uw[1])]
                UDw = [v[1] for lw in Q.incoming_edges(w) for v in P.outgoing_edges(lw[0])]
                UDw.extend([w]*r)
                assert sorted(DUw) == sorted(UDw), "D U - U D differs from %s I for vertex %s!"%(r, w)
        else:
            def check_vertex(w, P, Q):
                DUw = [v for uw in P.upper_covers(w) for v in Q.lower_covers(uw)]
                UDw = [v for lw in Q.lower_covers(w) for v in P.upper_covers(lw)]
                UDw.extend([w]*r)
                assert sorted(DUw) == sorted(UDw), "D U - U D differs from %s I for vertex %s!"%(r, w)

        P = cls.P_graph(n+2)
        Q = cls.Q_graph(n+2)
        for w in cls.vertices(n):
            check_vertex(w, Q, P)

    @classmethod
    def P_graph(cls, n):
        r"""
        Return the first n levels of the first dual graded graph.
        """
        try:
            cls._has_multiple_edges
        except AttributeError:
            cls._has_multiple_edges = False
        if cls._has_multiple_edges:
            D = DiGraph([(x,y,e) for k in range(n-1)
                         for x in cls.vertices(k)
                         for y in cls.vertices(k+1)
                         for e in cls._is_P_edge(y, x)], multiedges=True)
            D.layout_default = D.layout_acyclic
            return D
        else:
            return Poset(([w for k in range(n) for w in cls.vertices(k)],
                          lambda x,y: cls._is_P_edge(y, x)),
                         cover_relations=True)

    @classmethod
    def Q_graph(cls, n):
        r"""
        Return the first n levels of the second dual graded graph.
        """
        try:
            cls._has_multiple_edges
        except AttributeError:
            cls._has_multiple_edges = False
        if cls._has_multiple_edges:
            D = DiGraph([(x,y,e) for k in range(n-1)
                            for x in cls.vertices(k)
                            for y in cls.vertices(k+1)
                            for e in cls._is_Q_edge(y, x)], multiedges=True)
            D.layout_default = D.layout_acyclic
            return D
        else:
            return Poset(([w for k in range(n) for w in cls.vertices(k)],
                          lambda x,y: cls._is_Q_edge(y, x)),
                         cover_relations=True)

    def filling(self):
        r"""
        Return the filling of the diagram as a dictionary.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: G.filling()
            {(0, 1): 1, (1, 0): 1, (2, 1): 2}
        """
        return self._filling

    def conjugate(self):
        r"""
        Return the conjugate growth diagram of ``self``. This is
        the growth diagram with the filling reflected over the
        main diagonal.

        When the filling is a permutation, the conjugate filling
        corresponds to its inverse.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: Gc = G.conjugate()
            sage: (Gc.P_symbol(), Gc.Q_symbol()) == (G.Q_symbol(), G.P_symbol())
            True
        """
        F = {(j,i): v for (i,j),v in self._filling.items()}
        return self.parent()(filling = F)

    def rotate(self):
        r"""
        Return the growth diagram with the filling rotated by 180 degrees.

        For RSK-growth diagrams and rectangular fillings, this
        corresponds to evacuation of the P- and the Q-symbol.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: Gc = G.rotate()
            sage: ascii_art(Gc.P_symbol(), Gc.Q_symbol())
            1  1  1  1  1  2
            2        3

            sage: ascii_art(SemistandardTableau(Tableau(G.P_symbol())).evacuation(), SemistandardTableau(Tableau(G.Q_symbol())).evacuation())
            1  1  1  1  1  2
            2        3

        """
        max_row = max(i for i, _ in self._filling)
        max_col = max(j for _, j in self._filling)
        F = {(max_row-i,max_col-j): v for (i,j),v in self._filling.items()}
        return self.parent()(filling = F)

    def shape(self):
        r"""
        Return the shape of the growth diagram as a skew partition.

        .. WARNING::

            In the literature the label on the corner opposite of the
            origin of a rectangular filling is often called the shape
            of the filling.  This method returns the shape of the
            region instead.

        EXAMPLES::

            sage: GrowthDiagramRSK([1]).shape()
            [1] / []

        """
        return SkewPartition([self._lambda, self._mu])

    def out_labels(self):
        r"""
        Return the labels along the boundary opposite of the origin.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: G.out_labels()
            [[], [1], [1, 1], [3, 1], [1], []]
        """
        return self._out_labels

    def in_labels(self):
        r"""
        Return the labels along the boundary on the side of the origin.

        EXAMPLES::

            sage: G = GrowthDiagramRSK(labels=[[2,2],[3,2],[3,3],[3,2]]); G
            1 0
            sage: G.in_labels()
            [[2, 2], [2, 2], [2, 2], [3, 2]]
        """
        return self._in_labels

    def P_chain(self):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram.

        EXAMPLES::

            sage: G = GrowthDiagramBinWord([4, 1, 2, 3])
            sage: G.P_chain()
            [word: , word: 1, word: 11, word: 111, word: 1011]
        """
        if self.is_rectangular():
            if self._has_multiple_edges:
                return self._out_labels[(2*self._lambda[0]):][::-1]
            else:
                return self._out_labels[self._lambda[0]:][::-1]
        else:
            raise ValueError("The P symbol is only defined for rectangular shapes.")

    def Q_chain(self):
        r"""
        Return the labels along the horizontal boundary of a rectangular
        growth diagram.

        EXAMPLES::

            sage: G = GrowthDiagramBinWord([[0,1,0,0], [0,0,1,0], [0,0,0,1], [1,0,0,0]])
            sage: G.Q_chain()
            [word: , word: 1, word: 10, word: 101, word: 1011]
        """
        if self.is_rectangular():
            if self._has_multiple_edges:
                return self._out_labels[:(2*self._lambda[0]+1)]
            else:
                return self._out_labels[:self._lambda[0]+1]
        else:
            raise ValueError("The Q symbol is only defined for rectangular shapes.")

    def is_rectangular(self):
        r"""
        Return ``True`` if the shape of the growth diagram is rectangular.

        EXAMPLES::

            sage: GrowthDiagramRSK([2,3,1]).is_rectangular()
            True
            sage: GrowthDiagramRSK([[1,0,1],[0,1]]).is_rectangular()
            False
        """
        return all(x == 0 for x in self._mu) and all(x == self._lambda[0] for x in self._lambda)

    def to_word(self):
        r"""
        Return the filling as a word, if the shape is rectangular and
        there is at most one nonzero entry in each column, which must
        be 1.

        EXAMPLES::

            sage: w = [3,3,2,4,1]; G = GrowthDiagramRSK(w)
            sage: G
            0  0  0  0  1
            0  0  1  0  0
            1  1  0  0  0
            0  0  0  1  0
            sage: G.to_word()
            [3, 3, 2, 4, 1]

        """
        if self.is_rectangular():
            w = [0]*self._lambda[0]
            for ((i,j), v) in self._filling.items():
                if v != 0:
                    if v == 1:
                        if w[i] == 0:
                            w[i] = j+1
                        else:
                            raise ValueError("Can only convert fillings with at most one entry per column to words.")
                    elif v == -1:
                        if w[i] == 0:
                            w[i] = -(j+1)
                        else:
                            raise ValueError("Can only convert fillings with at most one entry per column to words.")
                    else:
                        raise ValueError("Can only convert 0-1 fillings to words.  Try ``to_biword``.")
            return w
        else:
            raise ValueError("Can only convert fillings of rectangular shapes to words.")

    def to_biword(self):
        r"""
        Return the filling as a biword, if the shape is rectangular.

        EXAMPLES::

            sage: P = Tableau([[1,2,2],[2]]); Q = Tableau([[1,3,3],[2]])
            sage: bw = RSK_inverse(P, Q); bw
            [[1, 2, 3, 3], [2, 1, 2, 2]]
            sage: G = GrowthDiagramRSK(labels = Q.to_chain()[:-1] + P.to_chain()[::-1]); G
            0  1  0
            1  0  2

            sage: P = SemistandardTableau([[1, 1, 2], [2]])
            sage: Q = SemistandardTableau([[1, 2, 2], [2]])
            sage: G = GrowthDiagramRSK(labels = Q.to_chain()[:-1] + P.to_chain()[::-1]); G
            0  2
            1  1
            sage: G.to_biword()
            ([1, 2, 2, 2], [2, 1, 1, 2])
            sage: RSK([1, 2, 2, 2], [2, 1, 1, 2])
            [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]

        """
        if self.is_rectangular():
            w1 = []
            w2 = []
            for ((i,j), v) in sorted(self._filling.items()):
                if v >= 0:
                    w1.extend([i+1]*v)
                    w2.extend([j+1]*v)
                else:
                    raise ValueError("Can only convert fillings with nonnegative entries to words.")
            return (w1, w2)
        else:
            raise ValueError("Can only convert fillings of rectangular shapes to words.")

    def __iter__(self):
        r"""
        Return the rows of the filling.

        TESTS::

            sage: G = GrowthDiagramRSK({(0,1):1, (1,0):1}, SkewPartition([[2,1],[1]]))
            sage: list(G)
            [[None, 1], [1]]

            sage: pi = Permutation([2,3,1,6,4,5])
            sage: G = GrowthDiagramRSK(pi)
            sage: list(G)
            [[0, 0, 1, 0, 0, 0],
             [1, 0, 0, 0, 0, 0],
             [0, 1, 0, 0, 0, 0],
             [0, 0, 0, 0, 1, 0],
             [0, 0, 0, 0, 0, 1],
             [0, 0, 0, 1, 0, 0]]
        """
        return ([None]*self._mu[r] + [self._filling.get((self._mu[r]+j,r), 0)
                                      for j in range(self._lambda[r]-self._mu[r])]
                for r in range(len(self._lambda)))

    def __repr__(self):
        r"""
        Print the filling of the growth diagram as a skew tableau.

        TESTS::

            sage: GrowthDiagramRSK({(0,1):1, (1,0):1}, SkewPartition([[2,1],[1]]))
            .  1
            1

            sage: GrowthDiagramRSK({(0,1):1, (2,0):1}, SkewPartition([[3,1],[1]]))
            .  0  1
            1
        """
        return SkewTableau(expr=[self._mu,
                                 [[self._filling.get((self._mu[r]+j,r), 0)
                                   for j in range(self._lambda[r]-self._mu[r])]
                                  for r in range(len(self._lambda))][::-1]])._repr_diagram()

    def __eq__(self, other):
        r"""
        Return ``True`` if the growth diagram ``other`` has the same
        shape and the same filling as ``self``.

        EXAMPLES:

        Equality ignores zeros in fillings::

            sage: G1 = GrowthDiagramRSK({(0, 1): 1, (1, 0): 1})
            sage: G2 = GrowthDiagramRSK({(0, 0): 0, (0, 1): 1, (1, 0): 1})
            sage: G1 == G2
            True

        Growth diagrams with different shapes are different::

            sage: G1 = GrowthDiagramRSK([[0,1,0],[1,0]])
            sage: G2 = GrowthDiagramRSK([[0,1,0],[1]])
            sage: G1 == G2
            False

        Growth diagrams with different rules are different::

            sage: G1 = GrowthDiagramRSK({(0, 1): 1, (1, 0): 1})
            sage: G2 = GrowthDiagramBinWord({(0, 1): 1, (1, 0): 1})
            sage: G1 == G2
            False
        """
        return (self.parent() == other.parent() and
                self._lambda == other._lambda and
                self._mu == other._mu and
                self._filling == other._filling)

    def __ne__(self, other):
        r"""
        Return ``True`` if the growth diagram ``other`` does not have the
        same shape and the same filling as ``self``.
        """
        return not self == other

    def _half_perimeter(self):
        r"""
        Return half the perimeter of the shape of the growth diagram.

        Assumes that ``self._lambda`` is already set.

        TESTS::

            sage: G = GrowthDiagramRSK({(0,1):1, (2,0):1}, SkewPartition([[3,1],[1]])); G
            .  0  1
            1
            sage: G._half_perimeter()
            6

        """
        if len(self._lambda) == 0:
            return 1
        else:
            return self._lambda[0]+len(self._lambda)+1

    def _shape_from_labels(self, labels):
        r"""
        Determine the shape of the growth diagram given a list of labels
        during initialisation.

        The shape can be determined from the labels if the size of
        each label differs from the size of its successor.

        Assumes that ``self._rank_function`` is set.

        Otherwise raise an error.

        TESTS::

            sage: labels = [[],[2],[1],[],[1],[]]
            sage: G = GrowthDiagramRSK(labels=labels); G
            0 1
            1
            1
            sage: G._shape_from_labels(G.out_labels())
            [2, 1, 1]

        """
        if self._has_multiple_edges:
            def right_left(la, mu, e):
                if self._rank_function(la) < self._rank_function(mu):
                    assert e in self._is_Q_edge(mu, la), "%s has smaller rank than %s but there is no edge of color %s in Q!" %(la, mu, e)
                    return 1
                elif self._rank_function(la) > self._rank_function(mu):
                    assert e in self._is_P_edge(la, mu), "%s has smaller rank than %s but there is no edge of color %s in in P!" %(mu, la, e)
                    return 0
                else:
                    raise ValueError("Can only determine the shape of the growth diagram if ranks of successive labels differ.")
            return Partitions().from_zero_one([right_left(labels[i], labels[i+2], labels[i+1]) for i in range(0, len(labels)-2, 2)])
        else:
            def right_left(la, mu):
                if self._rank_function(la) < self._rank_function(mu):
                    assert self._is_Q_edge(mu, la), "%s has smaller rank than %s but isn't covered by it in Q!" %(la, mu)
                    return 1
                elif self._rank_function(la) > self._rank_function(mu):
                    assert self._is_P_edge(la, mu), "%s has smaller rank than %s but isn't covered by it in P!" %(mu, la)
                    return 0
                else:
                    raise ValueError("Can only determine the shape of the growth diagram if ranks of successive labels differ.")
            return Partitions().from_zero_one([right_left(labels[i], labels[i+1]) for i in range(len(labels)-1)])

    def _check_labels(self, labels):
        r"""
        Check sanity of the parameter ``labels``.

        Assumes that ``self._lambda`` is already set.

        TESTS::

            sage: GrowthDiagramRSK(shape=[1], labels=[[], [1]])                 # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: The number of labels is 2, but for this shape we need 3.

            sage: GrowthDiagramRSK(labels=[[], [1], [2], [2,1]])                # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: The number of labels is 4, but for this shape we need 1.

        .. TODO::

            Can we do something more sensible when the chain of labels is strictly increasing?
        """
        half_perimeter = self._half_perimeter()
        if self._has_multiple_edges:
            assert is_odd(len(labels)), "Only a list of odd length can specify a path, but %s has even length."%s
            path_length = (len(labels)+1)/2
        else:
            path_length = len(labels)

        if path_length != half_perimeter:
            raise ValueError("The number of labels is %s, but for this shape we need %s."
                             %(path_length, half_perimeter))

    def _init_labels_forward_from_various_input(self, labels):
        r"""
        Return a list of labels decorating the boundary near the origin.

        Assumes that ``self._lambda`` is already set.

        TESTS:

        ``labels`` is None::

            sage: filling = []
            sage: shape = SkewPartition([[4,2,1,1],[2,1,1]])
            sage: G = GrowthDiagramRSK(filling, shape)
            sage: G.in_labels()                                                 # indirect doctest
            [[], [], [], [], [], [], [], [], []]

        ``labels`` is a list of partitions::

            sage: filling = []
            sage: labels = [[],[1],[],[1],[]]
            sage: shape = SkewPartition([[2,1],[1]])
            sage: G = GrowthDiagramRSK(filling=filling, shape=shape, labels=labels)
            sage: G.in_labels()                                                 # indirect doctest
            [[], [1], [], [1], []]

        """
        if labels is None:
            if self._has_multiple_edges:
                return [self._zero, self._zero_edge]*(self._half_perimeter()-1) + [self._zero]
            else:
                return [self._zero]*(self._half_perimeter())
        else:
            return labels

    def _init_shape_from_various_input(self, shape):
        r"""
        Return a pair of partitions describing the region of the growth
        diagram.

        Assumes that ``self._filling`` is already set.

        TESTS:

        ``shape`` is a skew partition::

            sage: filling = []
            sage: shape = SkewPartition([[4,2,1,1],[2,1,1]])
            sage: G = GrowthDiagramRSK(filling, shape)
            sage: G._lambda, G._mu                                              # indirect doctest
            ([4, 2, 1, 1], [2, 1, 1, 0])

        ``shape`` is a partition::

            sage: filling = []
            sage: shape = Partition([3,2,1,1])
            sage: G = GrowthDiagramRSK(filling, shape)
            sage: G._lambda, G._mu                                              # indirect doctest
            ([3, 2, 1, 1], [0, 0, 0, 0])

        """
        try:
            shape = Partition(shape)
            return (list(shape), [0]*len(shape))
        except ValueError:
            try:
                shape = SkewPartition(shape)
                return (list(shape[0]),
                        list(shape[1]) + [0]*(len(shape[0])-len(shape[1])))
            except ValueError:
                raise ValueError("Cannot make sense of shape %s" %shape)


    def _init_filling_and_shape_from_various_input(self, filling, shape):
        r"""
        Return a dict ``F``, such that ``F[(i,j)]`` is the element in row ``i``
        and column ``j``, and a pair of partitions describing the
        region of the growth diagram.

        TESTS:

        ``filling`` is a dict of coordinates::

            sage: pi = Permutation([2,3,1,6,4,5])
            sage: G = GrowthDiagramRSK({(i,pi[i]-1):1 for i in range(len(pi))})
            sage: G._filling                                                    # indirect doctest
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        ``filling`` is a dict of dicts::

            sage: G = GrowthDiagramRSK({i:{pi[i]-1:1} for i in range(len(pi))})
            sage: G._filling                                                    # indirect doctest
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        ``filling`` is a matrix::

            sage: G = GrowthDiagramRSK(pi.to_matrix())
            sage: G._filling                                                    # indirect doctest
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        ``filling`` is a permutation::

            sage: G = GrowthDiagramRSK(pi)
            sage: G._filling                                                    # indirect doctest
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        ``filling`` is a list::

            sage: G = GrowthDiagramRSK([3,1,4,1,5])
            sage: G._filling                                                    # indirect doctest
            {(0, 2): 1, (1, 0): 1, (2, 3): 1, (3, 0): 1, (4, 4): 1}
            sage: G.shape()
            [5, 5, 5, 5, 5] / []

        ``filling`` is a list of lists::

            sage: G = GrowthDiagramRSK([[1,0,1],[0,1]])
            sage: G._filling                                                    # indirect doctest
            {(0, 0): 1, (1, 1): 1, (2, 0): 1}
            sage: G.shape()
            [3, 2] / []

        ``filling`` is a list of lists and shape is given::

            sage: G = GrowthDiagramRSK([[1,0,1],[0,1]], shape=SkewPartition([[3,2],[1]]))
            sage: G._filling                                                    # indirect doctest
            {(0, 0): 1, (1, 1): 1, (2, 0): 1}
            sage: G.shape()
            [3, 2] / [1]

        ``filling`` is empty and shape is ``None``::

            sage: G = GrowthDiagramRSK({})
            sage: (G.filling(), G.shape())
            ({}, [] / [])

        """
        if isinstance(filling, dict):
            try:
                from six import itervalues
                v = next(itervalues(filling))
                if isinstance(v, dict):
                    # it is a dict of dicts
                    F = dict()
                    for (i, row) in filling.items():
                        for (j, v) in row.items():
                            if v != 0:
                                F[(i,j)] = int(v)
                else:
                    # it is dict of coordinates
                    F = {(i,j): v for ((i,j), v) in filling.items()
                         if v != 0}
            except StopIteration:
                # it is an empty dict of coordinates
                F = filling

        else:
            # it is a sequence
            F = dict()
            try:
                # it is a sequence of sequences
                for i, row in enumerate(filling):
                    for j, v in enumerate(row):
                        if v != 0:
                            F[(j,i)] = int(v)
                if shape is None:
                    shape = [len(row) for row in filling]

            except TypeError:
                # it is a word - for convenience we allow signed words
                for i, l in enumerate(filling):
                    if l > 0:
                        F[(i, l-1)] = 1
                    else:
                        F[(i, -l-1)] = -1

        if shape is None:
            if F == {}:
                shape = []
            else:
                # find bounding rectangle of ``filling``
                max_row = max(i for i, _ in F)+1
                max_col = max(j for _, j in F)+1
                shape = [max_row]*max_col

        return (F, self._init_shape_from_various_input(shape))

    def _grow(self):
        r"""
        Compute the labels on the boundary opposite of the origin, given
        the filling.

        TESTS::

            sage: pi = Permutation([1])
            sage: G = GrowthDiagramRSK(pi)
            sage: G._out_labels                                                 # indirect doctest
            [[], [1], []]

            sage: pi = Permutation([1,2])
            sage: G = GrowthDiagramRSK(pi)
            sage: G._out_labels                                                 # indirect doctest
            [[], [1], [2], [1], []]

            sage: pi = Permutation([2,1])
            sage: G = GrowthDiagramRSK(pi)
            sage: G._out_labels                                                 # indirect doctest
            [[], [1], [1, 1], [1], []]

            sage: G = GrowthDiagramRSK({(0,1):1, (1,0):1}, SkewPartition([[2,1],[1]]))
            sage: G._out_labels                                                 # indirect doctest
            [[], [1], [], [1], []]

            sage: G = GrowthDiagramRSK({(1,1):1}, SkewPartition([[2,2],[1]]), labels=[[],[],[1],[],[]])
            sage: G._out_labels                                                 # indirect doctest
            [[], [1], [2], [1], []]

            sage: G = GrowthDiagramBinWord({(1,1):1}, SkewPartition([[2,2],[1]]), labels=[[],[],[1],[],[]])
            sage: G._out_labels                                                 # indirect doctest
            [word: , word: 1, word: 11, word: 1, word: ]
        """
        labels = copy(self._in_labels)
        l = len(self._lambda)
        if self._has_multiple_edges:
            for r in range(l):
                for c in range(self._mu[r]+l-r, self._lambda[r]+l-r):
                    j = r
                    i = c-l+r
                    (labels[2*c-1],
                     labels[2*c],
                     labels[2*c+1]) = self._forward_rule(labels[2*c-2],
                                                         labels[2*c-1],
                                                         labels[2*c],
                                                         labels[2*c+1],
                                                         labels[2*c+2],
                                                         self._filling.get((i,j), 0))
        else:
            for r in range(l):
                for c in range(self._mu[r]+l-r, self._lambda[r]+l-r):
                    j = r
                    i = c-l+r
                    labels[c] = self._forward_rule(labels[c-1],
                                                   labels[c],
                                                   labels[c+1],
                                                   self._filling.get((i,j), 0))

        self._out_labels = labels

    def _shrink(self):
        r"""
        Compute the labels on the boundary near the origin, and the filling.

        TESTS::

            sage: filling = [[0,0,1,0,0,0,0], [0,1,0,0,0,0,0], [1,0,0,0,0,0,0], [0,0,0,1,0,0,0], [0,0,0,0,0,0,1], [0,0,0,0,0,1,0], [0,0,0,0,1,0,0]]
            sage: G = GrowthDiagramRSK(filling)
            sage: list(GrowthDiagramRSK(labels=G._out_labels)) == filling
            True

            sage: labels = [[], [1], []]
            sage: G = GrowthDiagramRSK(labels=labels)
            sage: G._filling
            {(0, 0): 1}
            sage: G._in_labels                                                  # indirect doctest
            [[], [], []]

            sage: labels = [[], [1], [2], [2,1], [1,1], [1], []]
            sage: G = GrowthDiagramRSK(labels=labels)
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 0): 1}
            sage: G._in_labels                                                  # indirect doctest
            [[], [], [], [], [], [], []]

            sage: labels = [[], [1], [2], [3], [3, 1], [3, 2], [4, 2], [4, 1], [3, 1], [2, 1], [1, 1], [1], []]
            sage: G = GrowthDiagramRSK(labels=labels)
            sage: G._filling                                                    # indirect doctest
            {(0, 1): 1, (1, 2): 1, (2, 5): 1, (3, 0): 1, (4, 3): 1, (5, 4): 1}

            sage: labels = [[],[1],[1],[2],[2],[2,1],[2]]
            sage: G = GrowthDiagramRSK(labels=labels)
            Traceback (most recent call last):
            ...
            ValueError: Can only determine the shape of the growth diagram if ranks of successive labels differ.
            sage: G = GrowthDiagramRSK(shape=[3,2,1], labels=labels)
            sage: G._filling                                                    # indirect doctest
            {(1, 0): 1}
            sage: G._in_labels
            [[], [], [], [], [1], [1], [2]]

            sage: labels = [[], [1],[1],[2],[2],[2,1],[2],[2,1],[1,1],[2,1],[1,1]]
            sage: G = GrowthDiagramRSK(shape=[5,4,3,2,1], labels=labels)
            sage: G._filling
            {(1, 2): 1, (2, 1): 1, (4, 0): 1}
            sage: G._in_labels                                                  # indirect doctest
            [[], [], [], [], [], [], [1], [1], [1], [1, 1], [1, 1]]

            sage: labels = [[], [1],[1],[2],[2],[2,1],[2],[2,1],[1,1],[2,1],[1,1]]
            sage: G = GrowthDiagramRSK(shape=SkewPartition([[5,4,3,2,1],[3,2,1]]), labels=labels)
            sage: G._filling
            {(1, 2): 1, (2, 1): 1, (4, 0): 1}
            sage: G._in_labels                                                  # indirect doctest
            [[], [], [], [1], [1], [1], [1], [1], [1], [1, 1], [1, 1]]

        """
        F = dict()
        labels = copy(self._out_labels)
        l = len(self._lambda)
        if self._has_multiple_edges:
            for r in range(l):
                for c in range(self._lambda[l-r-1]+r, self._mu[l-r-1]+r, -1):
                    j = l-r-1
                    i = c-r-1
                    (labels[2*c-1],
                     labels[2*c],
                     labels[2*c+1], v) = self._backward_rule(labels[2*c-2],
                                                             labels[2*c-1],
                                                             labels[2*c],
                                                             labels[2*c+1],
                                                             labels[2*c+2])
                    if v != 0:
                        F[(i,j)] = v

        else:
            for r in range(l):
                for c in range(self._lambda[l-r-1]+r, self._mu[l-r-1]+r, -1):
                    j = l-r-1
                    i = c-r-1
                    labels[c], v = self._backward_rule(labels[c-1],
                                                       labels[c],
                                                       labels[c+1])
                    if v != 0:
                        F[(i,j)] = v

        self._in_labels = labels
        self._filling = F

######################################################################
# Specific classes of growth diagrams
######################################################################

class GrowthDiagramShiftedShapes(GrowthDiagram):
    r"""
    A class modelling the Schensted correspondence for shifted
    shapes, which agrees with Sagan and Worley's and Haiman's
    insertion algorithms.

    EXAMPLES::

        sage: GrowthDiagramShiftedShapes([3,4,1,2]).out_labels()
        [[], 1, [1], 2, [2], 3, [3], 1, [3, 1], 0, [2, 1], 0, [2], 0, [1], 0, []]

    Check example just before Corollary 3.2 in [Sag1987]_::

        sage: G = GrowthDiagramShiftedShapes([2,6,5,1,7,4,3])
        sage: G.P_chain()
        [[], 0, [1], 0, [2], 0, [3], 0, [3, 1], 0, [3, 2], 0, [4, 2], 0, [5, 2]]
        sage: G.Q_chain()
        [[], 1, [1], 2, [2], 1, [2, 1], 3, [3, 1], 2, [4, 1], 3, [4, 2], 3, [5, 2]]

    .. automethod:: _forward_rule
    .. automethod:: _backward_rule

    TESTS::

        sage: SY = GrowthDiagramShiftedShapes
        sage: SY._zero
        []

        sage: SY._check_duality(4)

    Check that the rules are bijective::

        sage: all(SY(labels=SY(pi).out_labels()).to_word() == pi for pi in Permutations(5))
        True
        sage: pi = Permutations(10).random_element()
        sage: G = SY(pi)
        sage: list(SY(labels=G.out_labels())) == list(G)
        True
    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        # TODO: should check that the filling is standard
        if labels is not None:
            labels = [Partition(labels[i]) if is_even(i) else ZZ(labels[i]) for i in range(len(labels))]
        super(GrowthDiagramShiftedShapes, self).__init__(filling = filling,
                                                   shape = shape,
                                                   labels = labels)
    __init__.__doc__ = GrowthDiagram.__init__.__doc__

    _zero = Partition([])
    _has_multiple_edges = True

    @staticmethod
    def vertices(n):
        if n == 0:
            return [GrowthDiagramShiftedShapes._zero]
        else:
            return Partitions(n, max_slope=-1)

    @staticmethod
    def _rank_function(w):
        return w.size()

    @staticmethod
    def _is_Q_edge(w, v):
        try:
            l = SkewPartition([w, v]).cells()
        except ValueError:
            return []
        else:
            if l[0][1] == 0:
                return [1]   # black
            else:
                return [2,3] # blue, red

    @staticmethod
    def _is_P_edge(w, v):
        return [0] if w.contains(v) else []

    @staticmethod
    def _forward_rule(y, e, t, f, x, content):
        r"""
        Return the output path given two incident edges and the content.

        See [Fom1995]_ Lemma 4.5.1, page 38.

        INPUT:

        - ``y, e, t, f, x`` -- a path of three partitions and two
          colors from a cell in a growth diagram, labelled as::

              t f x
              e
              y

        - ``content`` -- 0 or 1, the content of the cell.

        OUTPUT:

        The two colors and the fourth partition g, z, h according to
        Sagan - Worley insertion.

        TESTS::

            sage: G = GrowthDiagramShiftedShapes
            sage: G._forward_rule([], 0, [], 0, [], 1)
            (1, [1], 0)

            sage: G._forward_rule([1], 0, [1], 0, [1], 1)
            (2, [2], 0)

        if ``x != y``::

            sage: G._forward_rule([3], 0, [2], 1, [2,1], 0)
            (1, [3, 1], 0)

            sage: G._forward_rule([2,1], 0, [2], 2, [3], 0)
            (2, [3, 1], 0)

        if ``x == y != t``::

            sage: G._forward_rule([3], 0, [2], 2, [3], 0)
            (1, [3, 1], 0)

            sage: G._forward_rule([3,1], 0, [2,1], 2, [3,1], 0)
            (2, [3, 2], 0)

            sage: G._forward_rule([2,1], 0, [2], 1, [2,1], 0)
            (3, [3, 1], 0)

            sage: G._forward_rule([3], 0, [2], 3, [3], 0)
            (3, [4], 0)

        """
        assert e == 0, "The P-graph should not be colored"
        h = 0
        if x == t == y:
            assert f == 0, "Degenerate edge f should have color 0"
            if content == 0:
                g, z = 0, x
            elif content == 1:
                if len(x) == 0:
                    g, z = 1, Partition(x).add_cell(0) # black
                else:
                    g, z = 2, Partition(x).add_cell(0) # blue
            else:
                raise NotImplementedError
        elif content != 0:
            raise ValueError("For y=%s, t=%s, x=%s, the content should be 0 but is %s" %(y, t, x, content))
        elif x != t == y:
            g, z = f, x
        elif x == t != y:
            assert f == 0, "Degenerate edge f should have color 0"
            g, z = f, y
        else:
            if x != y:
                row = SkewPartition([x, t]).cells()[0][0]
                g, z = f, Partition(y).add_cell(row)
            elif x == y != t and f == 2: # blue
                row = 1+SkewPartition([x, t]).cells()[0][0]
                if row == len(y):
                    g, z = 1, Partition(y).add_cell(row) # black
                else:
                    g, z = 2, Partition(y).add_cell(row) # blue
            elif x == y != t and f in [1, 3]: # black or red
                c = SkewPartition([x, t]).cells()[0]
                col = c[0] + c[1] + 1
                # print y, t, x, c, col
                for i in range(len(y)):
                    if i+y[i] == col:
                        z = y[:i] + [y[i]+1] + y[i+1:]
                        break
                g = 3
            else:
                raise NotImplementedError
        return g, Partition(z), h

    @staticmethod
    def _backward_rule(y, g, z, h, x):
        r"""
        Return the input path and the content given two incident edges.

        See [Fom1995]_ Lemma 4.5.1, page 38.

        INPUT:

        - ``y, g, z, h, x`` -- a path of three partitions and two
          colors from a cell in a growth diagram, labelled as::

                  x
                  h
              y g z

        OUTPUT:

        A tuple ``(e, t, f, content)`` consisting of the shape ``t``
        of the fourth word, the colours of the incident edges and the
        content of the cell acording to Sagan - Worley insertion.
        """
        assert h == 0, "The P-graph should not be colored"

        if x == y == z:
            assert g == 0, "Degenerate edge g should have color 0"
            return (0, x, 0, 0)
        elif x == z != y:
            return (0, y, g, 0)
        elif x != z == y:
            assert g == 0, "Degenerate edge g should have color 0"
            return (0, x, 0, 0)
        else:
            if x != y:
                row = SkewPartition([z, x]).cells()[0][0]
                return (0, Partition(y).remove_cell(row), g, 0)
            else:
                row, col = SkewPartition([z, x]).cells()[0]
                if row > 0 and g in [1, 2]: # black or blue
                    return (0, Partition(y).remove_cell(row-1), 2, 0)
                elif row == 0 and g in [1, 2]: # black or blue
                    return (0, y, 0, 1)
                else:
                    # find last cell in column col-1
                    for i in range(len(y)-1,-1,-1):
                        if i+y[i] == col + row:
                            if y[i] == 1:
                                t = y[:i]
                                return (0, t, 1, 0)
                            else:
                                t = y[:i] + [y[i]-1] + y[i+1:]
                                return (0, t, 3, 0)
                    raise ValueError("This should not happen.")

class GrowthDiagramLLMSClass(GrowthDiagram):
    pass

def GrowthDiagramLLMS(k):
    """
    """
    class GrowthDiagramLLMS(GrowthDiagramLLMSClass):
        _k = k
        _zero = Core([], k)
        _zero_edge = None # to prevent confusion with the edge labelled with content 0

    return GrowthDiagramLLMS


class GrowthDiagramLLMSClass(GrowthDiagram):
    r"""
    A class modelling the Schensted correspondence for affine
    permutations.

    EXAMPLES::

    Check example of Figure 1 in [LamShi2007]_::

        sage: G = GrowthDiagramLLMS(3)([4,1,2,6,3,5])
        sage: G.P_symbol().pp()
        -1 -2 -3 -5
         3  5
        -4 -6
         5
         6

        sage: G.Q_symbol().pp()
        1  3  4  5
        2  5
        3  6
        5
        6

    Check Example 6.2 in [LLMSSZ2013]_::

        sage: G = GrowthDiagramLLMS(3)([4,1,3,2])
        sage: G.P_symbol().pp()
        -1 -2  3
        -3
        -4

        sage: G.Q_symbol().pp()
        1  3  4
        2
        3

    .. automethod:: _forward_rule
    .. automethod:: _backward_rule

    TESTS::

        sage: G = GrowthDiagramLLMS(3)
        sage: G._zero
        []

        sage: G._check_duality(4)

    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        # TODO: should check that the filling is standard
        if labels is not None:
            labels = [Core(labels[i], self._k) if is_even(i) else labels[i] for i in range(len(labels))]
        super(GrowthDiagramLLMSClass, self).__init__(filling = filling,
                                                     shape = shape,
                                                     labels = labels)
    __init__.__doc__ = GrowthDiagram.__init__.__doc__

    _has_multiple_edges = True

    @staticmethod
    def _rank_function(w):
        return w.length()

    @classmethod
    def vertices(cls, n):
        return Cores(cls._k, length=n)

    @staticmethod
    def _is_Q_edge(w, v):
        return [None] if w in v.weak_covers() else []

    @staticmethod
    def _is_P_edge(w, v):
        """For two k-cores v and w containing v, there are as many edges as
        there are components in the skew partition w/v.  These
        components are ribbons, and therefore contain a unique cell
        with maximal content.  We index the edge with this content.
        """
        if w in v.strong_covers():
            T = SkewPartition([w.to_partition(), v.to_partition()])
            return [max([j-i for i,j in c]) for c in T.cell_poset().connected_components()]
        else:
            return []

    def P_symbol(self):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram as a skew tableau.

        EXAMPLES::

            sage: G = GrowthDiagramLLMS(4)([3,4,1,2])
            sage: G.P_symbol().pp()
            -1 -2
            -3 -4

        """
        C = self.P_chain()
        T = SkewTableau(chain = C[::2])
        S = T.to_list()
        for entry, content in enumerate(C[1::2], 1):
            for i,j in T.cells_containing(entry):
                if j-i == content:
                    S[i][j] = -S[i][j]
                    break
        return StrongTableau(S, self._k-1)

    def Q_symbol(self):
        r"""
        Return the labels along the horizontal boundary of a rectangular
        growth diagram as a skew tableau.

        EXAMPLES::

            sage: G = GrowthDiagramLLMS(4)([3,4,1,2])
            sage: G.Q_symbol().pp()
            1 2
            3 4

        """
        return WeakTableau(SkewTableau(chain = self.Q_chain()[::2]), self._k-1)

    @classmethod
    def _forward_rule(cls, y, e, t, f, x, content):
        r"""
        Return the output path given two incident edges and the content.

        See [LamShi2007]_ Section 3.4 and [LLMSSZ2013]_ Section 6.3.

        INPUT:

        - ``y, e, t, f, x`` -- a path of three partitions and two
          colors from a cell in a growth diagram, labelled as::

              t f x
              e
              y

        - ``content`` -- 0 or 1, the content of the cell.

        OUTPUT:

        The two colors and the fourth partition g, z, h according to
        LLMS insertion.

        TESTS::

            sage: G = GrowthDiagramLLMS(3)
            sage: H = GrowthDiagramLLMS(4)
            sage: Z = G._zero
            sage: G._forward_rule(Z, None, Z, None, Z, 0)
            (None, [], None)

            sage: G._forward_rule(Z, None, Z, None, Z, 1)
            (None, [1], 0)

            sage: Y = Core([3,1,1], 3);
            sage: G._forward_rule(Y, None, Y, None, Y, 1)
            (None, [4, 2, 1, 1], 3)

        if ``x != y``::

            sage: Y = Core([1,1], 3); T = Core([1], 3); X = Core([2], 3)
            sage: G._forward_rule(Y, -1, T, None, X, 0)
            (None, [2, 1, 1], -1)

            sage: Y = Core([2], 4); T = Core([1], 4); X = Core([1,1], 4)
            sage: H._forward_rule(Y, 1, T, None, X, 0)
            (None, [2, 1], 1)

            sage: Y = Core([2,1,1], 3); T = Core([2], 3); X = Core([3,1], 3)
            sage: G._forward_rule(Y, -1, T, None, X, 0)
            (None, [3, 1, 1], -2)


        if ``x == y != t``::

            sage: Y = Core([1], 3); T = Core([], 3); X = Core([1], 3)
            sage: G._forward_rule(Y, 0, T, None, X, 0)
            (None, [1, 1], -1)

            sage: Y = Core([1], 4); T = Core([], 4); X = Core([1], 4)
            sage: H._forward_rule(Y, 0, T, None, X, 0)
            (None, [1, 1], -1)

            sage: Y = Core([2,1], 4); T = Core([1,1], 4); X = Core([2,1], 4)
            sage: H._forward_rule(Y, 1, T, None, X, 0)
            (None, [2, 2], 0)

        """
        assert f == None, "The Q-graph should not be colored"
        g = None
        if x == t == y:
            assert e == None, "Degenerate edge e should have color None"
            if content == 0:
                z, h = x, None
            elif content == 1:
                if t.size() == 0:
                    z = t.affine_symmetric_group_simple_action(0)
                else:
                    z = t.affine_symmetric_group_simple_action((t[0])%(cls._k))
                h = z[0]-1
            else:
                raise ValueError("Should not happen.")
        elif content != 0:
            raise ValueError("For y=%s, t=%s, x=%s, the content should be 0 but is %s" %(y, t, x, content))
        elif x != t == y:
            assert e == None, "Degenerate edge e should have color None"
            z, h = x, e
        elif x == t != y:
            z, h = y, e
        else: #  x != t and y != t
            qx = SkewPartition([x.to_partition(), t.to_partition()])
            qy = SkewPartition([y.to_partition(), t.to_partition()])
            if not all(c in qx.cells() for c in qy.cells()):
                res = [(j-i)%(cls._k) for i,j in qx.cells()]
                assert len(set(res)) == 1
                r = res[0]
                z = y.affine_symmetric_group_simple_action(r)
                if e%(cls._k) == r:
                    h = e-1
                else:
                    h = e
            elif x == y != t:
                # the addable cell with largest content at most e
                cprime = sorted([c for c in y.to_partition().addable_cells()
                                 if c[1]-c[0] <= e],
                                key = lambda c: -(c[1]-c[0]))[0]
                h = cprime[1]-cprime[0]
                z = y.affine_symmetric_group_simple_action(h%(cls._k))

        return g, z, h

class GrowthDiagramBinWord(GrowthDiagram):
    r"""
    A class modelling a Schensted-like correspondence for binary
    words.

    EXAMPLES::

        sage: pi = Permutation([4,1,8,3,6,5,2,7,9]); G = GrowthDiagramBinWord(pi); G
        0  1  0  0  0  0  0  0  0
        0  0  0  0  0  0  1  0  0
        0  0  0  1  0  0  0  0  0
        1  0  0  0  0  0  0  0  0
        0  0  0  0  0  1  0  0  0
        0  0  0  0  1  0  0  0  0
        0  0  0  0  0  0  0  1  0
        0  0  1  0  0  0  0  0  0
        0  0  0  0  0  0  0  0  1
        sage: G.out_labels()[9]
        word: 101010011

    The Kleitman Greene invariant is the descent word::

        sage: pi.descents(from_zero=False)
        [1, 3, 5, 6]

    .. automethod:: _forward_rule
    .. automethod:: _backward_rule

    TESTS::

        sage: BW = GrowthDiagramBinWord
        sage: BW._zero
        word:
        sage: G = BW(labels = [[1,1],[1,1,0],[0,1]])
        Traceback (most recent call last):
        ...
        AssertionError: 01 has smaller rank than 110 but isn't covered by it in P!

        sage: G = BW(labels = [[1,1],[1,0,1],[0,1]])
        Traceback (most recent call last):
        ...
        AssertionError: 11 has smaller rank than 101 but isn't covered by it in Q!

    Check duality::

        sage: BW._check_duality(4)

    Check that the rules are bijective::

        sage: all(BW(labels=BW(pi).out_labels()).to_word() == pi for pi in Permutations(4))
        True
        sage: pi = Permutations(10).random_element()
        sage: G = BW(pi)
        sage: list(BW(labels=G.out_labels())) == list(G)
        True

    Test that the Kleitman Greene invariant is indeed the descent word::

        sage: r=4; all(Word([0 if i in w.descents(from_zero=False) else 1 for i in range(r)]) == BW(w).out_labels()[r] for w in Permutations(r))
        True

    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        # TODO: should check that the filling is standard
        if labels is not None:
            labels = [Word(la, alphabet=[0,1]) for la in labels]
        super(GrowthDiagramBinWord, self).__init__(filling = filling,
                                                   shape = shape,
                                                   labels = labels)
    __init__.__doc__ = GrowthDiagram.__init__.__doc__

    _zero = Word([], alphabet=[0,1])

    @staticmethod
    def vertices(n):
        if n == 0:
            return [GrowthDiagramBinWord._zero]
        else:
            w1 = Word([1], [0,1])
            return [w1+w for w in Words([0,1], n-1)]

    @staticmethod
    def _rank_function(w):
        return len(w)

    @staticmethod
    def _is_Q_edge(w, v):
        return w[:-1] == v

    @staticmethod
    def _is_P_edge(w, v):
        return len(w) == len(v) + 1 and v.is_subword_of(w)

    @staticmethod
    def _forward_rule(y, t, x, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Fom1995]_ Lemma 4.6.1, page 40.

        INPUT:

        - ``y, t, x`` -- three binary words from a cell in a growth
          diagram, labelled as::

              t x
              y

        - ``content`` -- 0 or 1, the content of the cell.

        OUTPUT:

        The fourth binary word z according to Viennot's bijection [Vie1983]_.

        TESTS::

            sage: G = GrowthDiagramBinWord

            sage: G._forward_rule([], [], [], 1)
            word: 1

            sage: G._forward_rule([1], [1], [1], 1)
            word: 11

        if ``x != y`` append last letter of ``x`` to ``y``::

            sage: G._forward_rule([1,0], [1], [1,1], 0)
            word: 101

        if ``x == y != t`` append ``0`` to ``y``::

            sage: G._forward_rule([1,1], [1], [1,1], 0)
            word: 110

        """
        if x == t == y:
            if content == 0:
                z = x
            elif content == 1:
                z = Word(list(y) + [1], alphabet=[0,1])
            else:
                raise NotImplementedError
        elif content != 0:
            raise ValueError("For y=%s, t=%s, x=%s, the content should be 0 but is %s" %(y, t, x, content))
        elif x != t == y:
            z = x
        elif x == t != y:
            z = y
        else:
            if x != y:
                z = Word(list(y) + [x[-1]], alphabet=[0,1])
            elif x == y != t:
                z = Word(list(y) + [0], alphabet=[0,1])
            else:
                raise NotImplementedError
        return z

    @staticmethod
    def _backward_rule(y, z, x):
        r"""
        Return the content and the input shape.

        See [Fom1995]_ Lemma 4.6.1, page 40.

        - ``y, z, x`` -- three binary words from a cell in a growth diagram,
          labelled as::

                x
              y z

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        word and the content of the cell acording to Viennot's
        bijection [Vie1983]_.

        TESTS::

            sage: w = [4,1,8,3,6,5,2,7,9]; G = GrowthDiagramBinWord(w);
            sage: GrowthDiagramBinWord(labels=G.out_labels()).to_word() == w    # indirect doctest
            True

        """
        if x == y == z:
            return (x, 0)
        elif x == z != y:
            return (y, 0)
        elif x != z == y:
            return (x, 0)
        else:
            if x == y and len(z) > 0 and z[-1] == 1:
                return (x, 1)
            else:
                return (x[:-1], 0)

class GrowthDiagramSylvester(GrowthDiagram):
    r"""
    A class modelling a Schensted-like correspondence for binary
    trees.

    EXAMPLES:

    From [Nze2007]_::

        sage: pi = Permutation([3,5,1,4,2,6]); G = GrowthDiagramSylvester(pi); G
        0  0  1  0  0  0
        0  0  0  0  1  0
        1  0  0  0  0  0
        0  0  0  1  0  0
        0  1  0  0  0  0
        0  0  0  0  0  1
        sage: ascii_art(G.out_labels()[len(pi)])
          __o__
         /     \
        o       o
         \     / \
          o   o   o

    .. automethod:: _forward_rule

    TESTS::

        sage: SY = GrowthDiagramSylvester
        sage: SY._zero
        .

        sage: for pi in Permutations(5):
        ....:     G = SY(pi)
        ....:     R = LabelledBinaryTree(None)
        ....:     for i in pi:
        ....:         R = R.binary_search_insert(i)
        ....:     assert BinaryTree(R) == G.P_chain()[-1]

        sage: B = BinaryTree; R = B([None,[]]); L = B([[],None]); T = B([[],[]]); S = B([L,None])
        sage: G = SY(labels = [R, T, R])
        Traceback (most recent call last):
        ...
        AssertionError: [., [., .]] has smaller rank than [[., .], [., .]] but isn't covered by it in P!

        sage: G = SY(labels = [R, S, R])
        Traceback (most recent call last):
        ...
        AssertionError: [., [., .]] has smaller rank than [[[., .], .], .] but isn't covered by it in Q!

    Check duality::

        sage: SY._check_duality(4)

    Check that the rules are bijective::

        sage: all(SY(labels=SY(pi).out_labels()).to_word() == pi for pi in Permutations(4))
        True
        sage: pi = Permutations(10).random_element()
        sage: G = SY(pi)
        sage: list(SY(labels=G.out_labels())) == list(G)
        True

    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        # TODO: should check that the filling is standard
        if labels is not None:
            labels = [BinaryTree(la) for la in labels]
        super(GrowthDiagramSylvester, self).__init__(filling = filling,
                                                     shape = shape,
                                                     labels = labels)
    __init__.__doc__ = GrowthDiagram.__init__.__doc__

    _zero = BinaryTree()

    @staticmethod
    def vertices(n):
        return BinaryTrees(n)

    @staticmethod
    def _rank_function(w):
        return w.node_number()

    @staticmethod
    def _is_Q_edge(w, v):
        def is_subtree(T1, T2):
            if T2.is_empty():
                return False
            elif T2[0].is_empty() and T2[1].is_empty():
                return T1.is_empty()
            elif T1.is_empty():
                return False
            else:
                return ((T1[0] == T2[0] and is_subtree(T1[1], T2[1])) or
                        (T1[1] == T2[1] and is_subtree(T1[0], T2[0])))
        return is_subtree(v, w)

    @staticmethod
    def _is_P_edge(w, v):
        if w.is_empty():
            return False
        else:
            return v == GrowthDiagramSylvester._delete_right_most_node(w)

    @staticmethod
    def _delete_right_most_node(b):
        """
        Return the tree obtained by deleting the right most node from ``b``.
        """
        if b.is_empty():
            raise ValueError("Cannot delete right most node from empty tree")
        elif b[1].is_empty():
            return b[0]
        else:
            return BinaryTree([b[0], GrowthDiagramSylvester._delete_right_most_node(b[1])])

    @staticmethod
    def _forward_rule(x, t, y, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Nze2007]_, page 9.

        INPUT:

        - ``y, t, x`` -- three binary trees from a cell in a growth
          diagram, labelled as::

              t y
              x

        - ``content`` -- 0 or 1, the content of the cell.

        OUTPUT:

        The fourth binary tree z.

        TESTS::

            sage: G = GrowthDiagramSylvester
            sage: B = BinaryTree; E = B(); N = B([]); L = B([[],None]); R = B([None,[]]); T = B([[],[]])

            sage: ascii_art(G._forward_rule(E, E, E, 1))
            o
            sage: ascii_art(G._forward_rule(N, N, N, 1))
            o
             \
              o
            sage: ascii_art(G._forward_rule(L, L, L, 1))
              o
             / \
            o   o
            sage: ascii_art(G._forward_rule(R, R, R, 1))
            o
             \
              o
               \
                o

        if ``x != y``, obtain ``z`` from ``x`` adding a node such
        that deleting the right most gives ``y``::

            sage: ascii_art(G._forward_rule(R, N, L, 0))
              o
             / \
            o   o

            sage: ascii_art(G._forward_rule(L, N, R, 0))
              o
             /
            o
             \
              o

        if ``x == y != t``, obtain ``z`` from ``y`` by adding a node
        as left child to the right most node::

            sage: ascii_art(G._forward_rule(N, E, N, 0))
              o
             /
            o
            sage: ascii_art(G._forward_rule(T, L, T, 0))
              _o_
             /   \
            o     o
                 /
                o
            sage: ascii_art(G._forward_rule(L, N, L, 0))
                o
               /
              o
             /
            o
            sage: ascii_art(G._forward_rule(R, N, R, 0))
            o
             \
              o
             /
            o

        """
        def successors(b):
            """
            Return all trees obtained from ``b`` by adding a node.
            """
            if b.is_empty():
                yield BinaryTree([])
            else:
                for t in successors(b[0]):
                    yield BinaryTree([t, b[1]])
                for t in successors(b[1]):
                    yield BinaryTree([b[0], t])

        def union(x, y):
            """
            Return the unique tree obtained by adding a node to ``x`` such
            that deleting the right most node gives ``y``.
            """
            for t in successors(x):
                if GrowthDiagramSylvester._delete_right_most_node(t) == y:
                    return t
            raise ValueError("Couldn't find union of %s and %s" %(x,y))

        if x == t == y:
            if content == 0:
                z = x
            elif content == 1:
                z = t.over(BinaryTree([]))
            else:
                raise NotImplementedError
        elif content != 0:
            raise ValueError("For y=%s, t=%s, x=%s, the content should be 0 but is %s" %(y, t, x, content))
        elif x != t == y:
            z = x
        elif x == t != y:
            z = y
        else:
            z = union(x, y)
        return z

    @staticmethod
    def _backward_rule(y, z, x):
        r"""
        Return the output shape given three shapes and the content.

        See [Nze2007]_, page 9.

        INPUT:

        - ``y, z, x`` -- three binary trees from a cell in a growth
          diagram, labelled as::

                y
              x z

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        binary tree t and the content of the cell.

        TESTS::

            sage: G = GrowthDiagramSylvester
            sage: B = BinaryTree; E = B(); N = B([]); L = B([[],None]); R = B([None,[]]); T = B([[],[]])

            sage: ascii_art(G._backward_rule(E, E, E))
            ( , 0 )
            sage: ascii_art(G._backward_rule(N, N, N))
            ( o, 0 )

        """
        if x == y == z:
            return (x, 0)
        elif x == z != y:
            return (y, 0)
        elif x != z == y:
            return (x, 0)
        else:
            if x == y and z == x.over(BinaryTree([])):
                return (x, 1)
            else:
                t = GrowthDiagramSylvester._delete_right_most_node(y)
                return (t, 0)


class GrowthDiagramYoungFibonacci(GrowthDiagram):
    r"""

    A class modelling a Schensted-like correspondence for
    Young-Fibonacci-tableaux.

    EXAMPLES::

        sage: G = GrowthDiagramYoungFibonacci([4,1,8,3,6,5,2,7,9]); G
        0  1  0  0  0  0  0  0  0
        0  0  0  0  0  0  1  0  0
        0  0  0  1  0  0  0  0  0
        1  0  0  0  0  0  0  0  0
        0  0  0  0  0  1  0  0  0
        0  0  0  0  1  0  0  0  0
        0  0  0  0  0  0  0  1  0
        0  0  1  0  0  0  0  0  0
        0  0  0  0  0  0  0  0  1
        sage: G.out_labels()[9]
        word: 122121

    The Kleitman Greene invariant is: take the last letter and the
    largest letter of the permutation and remove them.  If they
    coincide write 1, otherwise write 2.

    .. automethod:: _forward_rule
    .. automethod:: _backward_rule

    TESTS::

        sage: YF = GrowthDiagramYoungFibonacci
        sage: YF._zero
        word:

    Check duality::

        sage: YF._check_duality(4)

        sage: G = YF(labels = [[1],[1,0],[1]])
        Traceback (most recent call last):
        ...
        ValueError: 0 not in alphabet!

        sage: G = YF(labels = [[1,1],[1,2]])
        Traceback (most recent call last):
        ...
        AssertionError: 11 has smaller rank than 12 but isn't covered by it in Q!

        sage: G = YF(labels = [[1,2],[1,1]])
        Traceback (most recent call last):
        ...
        AssertionError: 11 has smaller rank than 12 but isn't covered by it in P!


        sage: all(YF(labels=YF(pi).out_labels()).to_word() == pi for pi in Permutations(4))
        True
        sage: pi = Permutations(10).random_element()
        sage: G = YF(pi)
        sage: list(YF(labels=G.out_labels())) == list(G)
        True

    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        # TODO: should check that the filling is standard
        if labels is not None:
            labels = [Word(la, alphabet=[1,2]) for la in labels]
        super(GrowthDiagramYoungFibonacci, self).__init__(filling = filling,
                                                          shape = shape,
                                                          labels = labels)
    __init__.__doc__ = GrowthDiagram.__init__.__doc__

    _zero = Word([], alphabet=[1,2])

    @staticmethod
    def vertices(n):
        if n == 0:
            return [GrowthDiagramYoungFibonacci._zero]
        else:
            return [Word(list(w), [1,2]) for w in Compositions(n, max_part=2)]

    @staticmethod
    def _rank_function(w):
        return sum(w)

    @staticmethod
    def _is_P_edge(w, v):
        def covers(c):
            for i in range(len(c)+1):
                d = list(c)
                d.insert(i, 1)
                yield Word(d, alphabet=[1,2])
                if i < len(c) and c[i] == 1:
                    d = list(c)
                    d[i] = 2
                    yield Word(d, alphabet=[1,2])
                    break
        return sum(w) == sum(v) + 1 and w in covers(v)

    _is_Q_edge = _is_P_edge

    @staticmethod
    def _forward_rule(shape3, shape2, shape1, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Fom1995]_ Lemma 4.4.1, page 35.

        INPUT:

        - ``shape3, shape2, shape1`` -- three Fibonacci words from a
          cell in a growth diagram, labelled as::

              shape2 shape1
              shape3

        - ``content`` -- 0 or 1, the content of the cell.

        OUTPUT:

        The fourth Fibonacci word.

        TESTS::

            sage: G = GrowthDiagramYoungFibonacci

            sage: G._forward_rule([], [], [], 1)
            word: 1

            sage: G._forward_rule([1], [1], [1], 1)
            word: 11

            sage: G._forward_rule([1,2], [1], [1,1], 0)
            word: 21

            sage: G._forward_rule([1,1], [1], [1,1], 0)
            word: 21

        """
        if shape1 == shape2 == shape3:
            if content == 0:
                r = shape1
            elif content == 1:
                r = Word([1]+list(shape3), alphabet=[1,2])
            else:
                raise NotImplementedError
        elif content != 0:
            raise ValueError("For shape3=%s, shape2=%s, shape1=%s, the content should be 0 but is %s" %(shape3, shape2, shape1, content))
        elif shape1 == shape2:
            r = shape3
        elif shape3 == shape2:
            r = shape1
        else:
            if shape1 != shape2 != shape3:
                r = Word([2] + list(shape2), alphabet=[1,2])
            else:
                raise NotImplementedError("For shape3=%s, shape2=%s, shape1=%s, content %s we have no rule." %(shape3, shape2, shape1, content))
        return r

    @staticmethod
    def _backward_rule(y, z, x):
        r"""
        Return the content and the input shape.

        See [Fom1995]_ Lemma 4.4.1, page 35.

        - ``y, z, x`` -- three Fibonacci words from a cell in a
          growth diagram, labelled as::

                x
              y z

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        word and the content of the cell.

        TESTS::

            sage: w = [4,1,8,3,6,5,2,7,9]; G = GrowthDiagramYoungFibonacci(w);
            sage: GrowthDiagramYoungFibonacci(labels=G.out_labels()).to_word() == w    # indirect doctest
            True

        """
        if x == y == z:
            return (x, 0)
        elif x == z != y:
            return (y, 0)
        elif x != z == y:
            return (x, 0)
        else:
            if z[0] == 1:
                return (z[1:], 1)
            elif z[0] == 2:
                return (z[1:], 0)

class GrowthDiagramOnPartitions(GrowthDiagram):
    r"""
    A class for growth diagrams on Young's lattice on integer
    partitions graded by size.

    TESTS::

        sage: G = GrowthDiagramBurge([])                                    # indirect doctest
        sage: G._zero
        []
        sage: G = GrowthDiagramBurge(labels = [[1],[1]])                    # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: Can only determine the shape of the growth diagram if ranks of successive labels differ.

    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        if labels is not None:
            labels = [Partition(la) for la in labels]
        self._zero = Partition([])
        self._rank_function = lambda p: p.size()
        super(GrowthDiagramOnPartitions, self).__init__(filling = filling,
                                                        shape = shape,
                                                        labels = labels)
    __init__.__doc__ = GrowthDiagram.__init__.__doc__

    def P_symbol(self):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram as a skew tableau.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: G.P_symbol().pp()
            1  2  2
            2
        """
        return SkewTableau(chain = self.P_chain())

    def Q_symbol(self):
        r"""
        Return the labels along the horizontal boundary of a rectangular
        growth diagram as a skew tableau.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: G.Q_symbol().pp()
            1  3  3
            2
        """
        return SkewTableau(chain = self.Q_chain())

class GrowthDiagramRSK(GrowthDiagramOnPartitions):
    r"""
    A class modelling Robinson-Schensted-Knuth insertion.

    EXAMPLES::

        sage: pi = Permutation([2,3,6,1,4,5])
        sage: G = GrowthDiagramRSK(pi)
        sage: G.P_symbol(), G.Q_symbol()
        ([[1, 3, 4, 5], [2, 6]], [[1, 2, 3, 6], [4, 5]])
        sage: RSK(pi)
        [[[1, 3, 4, 5], [2, 6]], [[1, 2, 3, 6], [4, 5]]]

    .. automethod:: _forward_rule
    .. automethod:: _backward_rule

    """
    @staticmethod
    def _forward_rule(shape3, shape2, shape1, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Kra2006]_ `(F^1 0)-(F^1 2)`.

        INPUT:

        - ``shape3, shape2, shape1`` -- three partitions from a cell in a
          growth diagram, labelled as::

              shape2 shape1
              shape3

        - ``content`` -- a non-negative integer, the content of the cell.

        OUTPUT:

        The fourth partition according to the Robinson-Schensted-Knuth
        correspondence.


        TESTS::

            sage: G = GrowthDiagramRSK
            sage: G._forward_rule([2,1],[2,1],[2,1],1)
            [3, 1]

            sage: G._forward_rule([1],[],[2],2)
            [4, 1]
        """
        carry = content
        shape4 = []
        while True:
            if shape1 == []:
                row1 = 0
            else:
                row1 = shape1[0]
            if shape3 == []:
                row3 = 0
            else:
                row3 = shape3[0]
            newPart = max(row1, row3) + carry
            if newPart == 0:
                # returning this as a Partition costs a lot of time
                return shape4[::-1]
            else:
                shape4 = [newPart] + shape4
                if shape2 == []:
                    carry = min(row1, row3)
                else:
                    carry = min(row1, row3) - shape2[0]
                shape1 = shape1[1:]
                shape2 = shape2[1:]
                shape3 = shape3[1:]

    @staticmethod
    def _backward_rule(shape3, shape4, shape1):
        r"""
        Return the content and the input shape.

        See [Kra2006]_ `(B^1 0)-(B^1 2)`.

        INPUT:

        - ``shape3, shape4, shape1`` -- three partitions from a cell in a
          growth diagram, labelled as::

                   shape1
            shape3 shape4

        OUTPUT:

        A pair ``(shape2, content)`` consisting of the shape of the
        fourth word acording to the Robinson-Schensted-Knuth
        correspondence and the content of the cell.

        TESTS::

            sage: w = [4,1,8,3,6,5,2,7,9]; G = GrowthDiagramRSK(w)
            sage: GrowthDiagramRSK(labels=G._out_labels).to_word() == w         # indirect doctest
            True

        """
        carry = 0
        i = len(shape4)
        shape2 = []
        while i > 0:
            if len(shape1) < i:
                row1 = 0
            else:
                row1 = shape1[i-1]
            if len(shape3) < i:
                row3 = 0
            else:
                row3 = shape3[i-1]
            shape2 = [min(row1, row3) - carry] + shape2
            carry = shape4[i-1] - max(row1, row3)
            i = i-1
        return (Partition(shape2), carry)

class GrowthDiagramBurge(GrowthDiagramOnPartitions):
    r"""
    A class modelling Burge insertion.

    EXAMPLES::

        sage: GrowthDiagramBurge(labels=[[],[1,1,1],[2,1,1,1],[2,1,1],[2,1],[1,1],[]])
        1  1
        0  1
        1  0
        1  0

    .. automethod:: _forward_rule
    .. automethod:: _backward_rule
    """
    @staticmethod
    def _forward_rule(shape3, shape2, shape1, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Kra2006]_ `(F^4 0)-(F^4 2)`.

        INPUT:

        - ``shape3, shape2, shape1`` -- three  from a cell in a growth diagram,
          labelled as::

              shape2 shape1
              shape3

        - ``content`` -- a non-negative integer, the content of the cell.

        OUTPUT:

        The fourth partition according to the Burge correspondence.

        TESTS::

            sage: G = GrowthDiagramBurge
            sage: G._forward_rule([2,1],[2,1],[2,1],1)
            [3, 1]

            sage: G._forward_rule([1],[],[2],2)
            [2, 1, 1, 1]

        """
        carry = content
        shape4 = []
        while True:
            if shape1 == []:
                row1 = 0
            else:
                row1 = shape1[0]
            if shape2 == []:
                row2 = 0
            else:
                row2 = shape2[0]
            if shape3 == []:
                row3 = 0
            else:
                row3 = shape3[0]
            newPart = max(row1, row3) + min(int(row1 == row2 == row3), carry)
            if newPart == 0:
                return Partition(shape4[::-1])
            else:
                shape4 = [newPart] + shape4
                carry = carry - min(int(row1 == row2 == row3), carry) + min(row1, row3) - row2
                shape1 = shape1[1:]
                shape2 = shape2[1:]
                shape3 = shape3[1:]

    @staticmethod
    def _backward_rule(shape3, shape4, shape1):
        r"""
        Return the content and the input shape.

        See [Kra2006]_ `(B^4 0)-(B^4 2)`.  There is a typo in the
        computation of carry in `(B^4 2)` in the arXiv version of the
        article, `\rho` must be replaced by `\lambda`.

        INPUT:

        - ``shape3, shape4, shape1`` -- three partitions from a cell in a
          growth diagram, labelled as::

                   shape1         mu
            shape3 shape4      nu lambda

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        partition acording to the Burge correspondence and the content of
        the cell.

        TESTS::

            sage: w = [4,1,8,3,6,5,2,7,9]; G = GrowthDiagramBurge(w)
            sage: GrowthDiagramBurge(labels=G._out_labels).to_word() == w       # indirect doctest
            True

            sage: G = GrowthDiagramBurge
            sage: G._backward_rule([1,1,1],[2,1,1,1],[2,1,1])
            ([1, 1], 0)

        """
        carry = 0
        shape2 = []
        i = len(shape4)
        while i > 0:
            mu_i = 0 if len(shape1) < i else shape1[i-1]
            la_i = 0 if len(shape4) < i else shape4[i-1]
            nu_i = 0 if len(shape3) < i else shape3[i-1]

            shape2 = [min(mu_i, nu_i) - min(int(mu_i == nu_i == la_i), carry)] + shape2
            carry = carry - min(int(mu_i == nu_i == la_i), carry) + la_i - max(mu_i, nu_i)
            i = i-1
        return (Partition(shape2), carry)

class GrowthDiagramDomino(GrowthDiagram):
    r"""
    A class modelling domino insertion.

    EXAMPLES:

    Figure 3 in [Lam2004]_::

        sage: G = GrowthDiagramDomino([[0,0,0,-1],[0,0,1,0],[-1,0,0,0],[0,1,0,0]]); G
         0  0  0 -1
         0  0  1  0
        -1  0  0  0
         0  1  0  0

        sage: ascii_art(G.P_symbol(), G.Q_symbol())
        1  2  4  1  2  2
        1  2  4  1  3  3
        3  3     4  4

    TESTS::

        sage: G = GrowthDiagramDomino([[0,1,0],[0,0,-1],[1,0,0]]); G
        0  1  0
        0  0 -1
        1  0  0

        sage: ascii_art(G.P_symbol(), G.Q_symbol())
        1  1  1  1
        2  3  2  2
        2  3  3  3

        sage: l = {pi: GrowthDiagramDomino(pi) for pi in SignedPermutations(4)}
        sage: len(Set([(G.P_symbol(), G.Q_symbol()) for G in l.values()]))
        384

        The spin of a domino tableau is half the number of vertical dominos:

        sage: def spin(T):
        ....:     return sum(2*len(set(row)) - len(row) for row in T)/4

        According to [Lam2004]_, the number of negative entries in
        the signed permutation equals the sum of the spins of the two
        associated tableaux:

        sage: all(G.filling().values().count(-1) == spin(G.P_symbol()) + spin(G.Q_symbol()) for G in l.values())
        True

        Negating all signs transposes all the partitions:

        sage: all(l[pi].P_symbol() == l[SignedPermutations(4)([-e for e in pi])].P_symbol().conjugate() for pi in l)
        True

        Check part of Theorem 4.2.3 in [vanLee1996]_:

        sage: def to_permutation(pi):
        ....:     pi1 = list(pi)
        ....:     n = len(pi1)
        ....:     pi2 = [-e for e in pi][::-1] + pi1
        ....:     return Permutation([e+n+1 if e<0 else e+n for e in pi2])
        sage: def good(pi):
        ....:     return GrowthDiagramDomino(pi).P_chain()[-1] == GrowthDiagramRSK(to_permutation(pi)).P_chain()[-1]
        sage: all(good(pi) for pi in SignedPermutations(4))
        True

        sage: G = GrowthDiagramDomino(labels = [[1],[2,1]])
        Traceback (most recent call last):
        ...
        AssertionError: [1] has smaller rank than [2, 1] but isn't covered by it in Q!

        sage: G = GrowthDiagramDomino(labels = [[2,1],[1]])
        Traceback (most recent call last):
        ...
        AssertionError: [1] has smaller rank than [2, 1] but isn't covered by it in P!

    .. automethod:: _forward_rule
    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        if labels is not None:
            labels = [Partition(la) for la in labels]
        self._zero = Partition([])
        self._rank_function = lambda p: p.size()
        def covering(w, v):
            try:
                (row_1, col_1), (row_2, col_2) = SkewPartition([w, v]).cells()
                return row_1 == row_2 or col_1 == col_2
            except ValueError:
                return False

        self._is_Q_edge = self._is_P_edge = lambda w, v: covering(w, v)
        super(GrowthDiagramDomino, self).__init__(filling = filling,
                                                  shape = shape,
                                                  labels = labels)
    __init__.__doc__ = GrowthDiagram.__init__.__doc__

    def P_symbol(self):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram as a skew tableau.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: G.P_symbol().pp()
            1  2  2
            2
        """
        return SkewTableau(chain = self.P_chain())

    def Q_symbol(self):
        r"""
        Return the labels along the horizontal boundary of a rectangular
        growth diagram as a skew tableau.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: G.Q_symbol().pp()
            1  3  3
            2
        """
        return SkewTableau(chain = self.Q_chain())

    @staticmethod
    def _forward_rule(shape3, shape2, shape1, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Lam2004]_ Section 3.1.

        INPUT:

        - ``shape3, shape2, shape1`` -- three partitions from a cell in a
          growth diagram, labelled as::

              shape2 shape1
              shape3

        - ``content`` -- -1, 0 or 1, the content of the cell.

        OUTPUT:

        The fourth partition according to domino insertion.

        TESTS::

            sage: G = GrowthDiagramDomino

        Rule 1::

            sage: G._forward_rule([], [], [], 1)
            [2]

            sage: G._forward_rule([1,1], [1,1], [1,1], 1)
            [3, 1]

        Rule 2::

            sage: G._forward_rule([1,1], [1,1], [1,1], -1)
            [1, 1, 1, 1]

        Rule 3::

            sage: G._forward_rule([1,1], [1,1], [2,2], 0)
            [2, 2]

        Rule 4::

            sage: G._forward_rule([2,2,2], [2,2], [3,3], 0)
            [3, 3, 2]

            sage: G._forward_rule([2], [], [1,1], 0)
            [2, 2]

            sage: G._forward_rule([1,1], [], [2], 0)
            [2, 2]

            sage: G._forward_rule([2], [], [2], 0)
            [2, 2]

            sage: G._forward_rule([4], [2], [4], 0)
            [4, 2]

            sage: G._forward_rule([1,1,1,1], [1,1], [1,1,1,1], 0)
            [2, 2, 1, 1]

            sage: G._forward_rule([2,1,1], [2], [4], 0)
            [4, 1, 1]
        """
        def union(la, mu):
            """
            Return the union of the two partitions.
            """
            from six.moves import zip_longest
            return [max(p,q) for (p,q) in zip_longest(la, mu, fillvalue=0)]

        if content not in [0,1,-1]:
            raise ValueError("Domino: The content of the filling must be in {-1,0,1}")

        if content == 1:
            assert shape1 == shape2 == shape3
            if shape2 == []:
                shape4 = [2]
            else:
                shape4 = [shape2[0] + 2] + shape2[1:]

        elif content == -1:
            assert shape1 == shape2 == shape3
            shape4 = shape2 + [1,1]

        elif content == 0 and (shape2 == shape1 or shape2 == shape3):
            shape4 = union(shape1, shape3)

        else:
            # content == 0 and shape2 differs from shape1 and shape3 by
            # dominos gamma1 and gamma3

            # the following is certainly very slow
            gamma3 = set(SkewPartition([shape3, shape2]).cells())
            gamma1 = set(SkewPartition([shape1, shape2]).cells())
            diff = gamma1.intersection(gamma3)
            cell1, cell2 = gamma3
            if len(diff) == 0:
                shape4 = union(shape1, shape3)

            elif len(diff) == 1:
                shape4 = copy(shape1)
                # diff is a single cell
                (k,l) = diff.pop()
                # add (k+1, l+1) to shape1
                # either (k, l+1) or (k+1, l) must also be added
                if shape4[k] <= l+1:
                    shape4[k] += 1
                    shape4[k+1] += 1
                else:
                    if len(shape4) <= k+1:
                        shape4 += [2]
                    else:
                        shape4[k+1] += 2

            # diff has size 2, that is shape1 == shape3
            elif cell1[0] == cell2[0]:
                shape4 = copy(shape1)
                # a horizontal domino - add 2 to row below of gamma
                if len(shape4) <= cell1[0]+1:
                    shape4 += [2]
                else:
                    shape4[cell1[0]+1] += 2

            else:
                shape4 = copy(shape1)
                # a vertical domino - add 2 to column right of gamma
                # find first row shorter than cell1[1]+1
                for r, p in enumerate(shape4):
                    if p <= cell1[1]+1:
                        shape4[r] += 1
                        shape4[r+1] += 1
                        break
                else:
                    raise NotImplementedError("Domino: cannot call forward rule with shapes %s and content %s"
                                              %((shape3, shape2, shape1), content))

        return shape4
