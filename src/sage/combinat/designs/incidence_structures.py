"""
Incidence structures.

An incidence structure is specified by a list of points, blocks, and
an incidence matrix ([1]_, [2]_).

REFERENCES:

.. [1] Block designs and incidence structures from wikipedia,
  :wikipedia:`Block_design`
  :wikipedia:`Incidence_structure`

.. [2] E. Assmus, J. Key, Designs and their codes, CUP, 1992.

AUTHORS:

- Peter Dobcsanyi and David Joyner (2007-2008)

  This is a significantly modified form of part of the module block_design.py
  (version 0.6) written by Peter Dobcsanyi peter@designtheory.org.

- Vincent Delecroix (2014): major rewrite
"""
#***************************************************************************
#                              Copyright (C) 2007                          #
#                                                                          #
#                Peter Dobcsanyi       and         David Joyner            #
#           <peter@designtheory.org>          <wdjoyner@gmail.com>         #
#                                                                          #
#                                                                          #
#    Distributed under the terms of the GNU General Public License (GPL)   #
#    as published by the Free Software Foundation; either version 2 of     #
#    the License, or (at your option) any later version.                   #
#                    http://www.gnu.org/licenses/                          #
#***************************************************************************


from sage.misc.superseded import deprecated_function_alias
from sage.misc.cachefunc import cached_method

from sage.rings.all import ZZ
from sage.rings.integer import Integer

def IncidenceStructureFromMatrix(M, name=None):
    """
    Deprecated function that builds an incidence structure from a matrix.

    You should now use ``designs.IncidenceStructure(incidence_matrix=M)``.

    INPUT:

    - ``M`` -- a binary matrix. Creates a set of "points" from the rows and a
      set of "blocks" from the columns.

    EXAMPLES::

        sage: BD1 = designs.IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
        sage: M = BD1.incidence_matrix()
        sage: BD2 = IncidenceStructureFromMatrix(M)
        doctest:1: DeprecationWarning: IncidenceStructureFromMatrix is deprecated.
        Please use designs.IncidenceStructure(incidence_matrix=M) instead.
        See http://trac.sagemath.org/16553 for details.
        sage: BD1 == BD2
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(16553, 'IncidenceStructureFromMatrix is deprecated. Please use designs.IncidenceStructure(incidence_matrix=M) instead.')
    return IncidenceStructure(incidence_matrix=M, name=name)

class IncidenceStructure(object):
    """
    A base class for incidence structure (or block design) with explicit ground
    set and blocks.

    INPUT:

    - ``points`` -- the underlying set. If it is an integer `v`, then the set is
      considered to be `{0, ..., v-1}`.

    - ``blocks`` -- the blocks (might be any iterable)

    - ``incidence_matrix`` -- the incidence matrix

    - ``name`` (a string, such as "Fano plane").

    - ``check`` -- whether to check the input

    EXAMPLES::

    An incidence structure can be constructed by giving the number of points and the list of
    blocks::

        sage: designs.IncidenceStructure(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
        Incidence structure with 7 points and 7 blocks

    Or by its adjacency matrix (a `{0,1}`-matrix in which rows are indexed by
    points and columns by blocks)::

        sage: m = matrix([[0,1,0],[0,0,1],[1,0,1],[1,1,1]])
        sage: designs.BlockDesign(m)
        Incidence structure with 4 points and 3 blocks

    The order of the points and blocks do not matter (see :trac:`11333`)::

        sage: A = designs.BlockDesign([0,1,2], [[0],[0,2]])
        sage: B = designs.BlockDesign([1,0,2], [[0],[2,0]])
        sage: B == A
        True

        sage: C = designs.BlockDesign(2, [[0], [1,0]])
        sage: D = designs.BlockDesign(2, [[0,1], [0]])
        sage: C == D
        True
    """
    def __init__(self, points=None, blocks=None, incidence_matrix=None,
            name=None, check=True, test=None):
        """
        TESTS::

            sage: designs.IncidenceStructure(3, [[4]])
            Traceback (most recent call last):
            ...
            ValueError: Point 4 is not in the base set.

            sage: designs.IncidenceStructure(3, [[0,1],[0,2]], test=True)
            doctest:...: DeprecationWarning: the keyword test is deprecated,
            use check instead
            See http://trac.sagemath.org/16553 for details.
            Incidence structure with 3 points and 2 blocks

            sage: designs.IncidenceStructure(2, [[0,1,2,3,4,5]], test=False)
            Incidence structure with 2 points and 1 blocks
        """
        if test is not None:
            from sage.misc.superseded import deprecation
            deprecation(16553, "the keyword test is deprecated, use check instead")
            check = test

        from sage.matrix.constructor import matrix
        from sage.structure.element import Matrix

        if isinstance(points, Matrix):
            incidence_matrix = points
            points = None

        if points is None:
            assert incidence_matrix is not None
            M = matrix(incidence_matrix)
            v = M.nrows()
            self._points = range(v)
            self._blocks = sorted(M.nonzero_positions_in_column(i) for i in range(M.ncols()))
            self._blocks.sort()

        else:
            assert incidence_matrix is None
            if isinstance(points, (int,Integer)):
                self._points = range(points)
            else:
                self._points = sorted(points)
            v = len(self._points)
            self._blocks = [sorted(block) for block in blocks]
            self._blocks.sort()
            if check:
                for block in self._blocks:
                    for x in block:
                        if x not in self._points:
                            raise ValueError('Point {} is not in the base set.'.format(x))
                    if len(block) != len(set(block)):
                        raise ValueError("Repeated element in block {}".format(block))

        self._block_sizes = frozenset(len(b) for b in self._blocks)
        self._incidence_matrix = None
        self._name = str(name) if name is not None else 'IncidenceStructure'
        self._degrees = []

    def __iter__(self):
        """
        Iterator over the blocks.

        EXAMPLE::

            sage: sts = designs.steiner_triple_system(9)
            sage: list(sts)
            [[0, 1, 5], [0, 2, 4], [0, 3, 6], [0, 7, 8], [1, 2, 3], [1, 4, 7], [1, 6, 8], [2, 5, 8], [2, 6, 7], [3, 4, 8], [3, 5, 7], [4, 5, 6]]
        """
        return iter(self._blocks)

    def __repr__(self):
        """
        A print method.

        EXAMPLES::

            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD
            Incidence structure with 7 points and 7 blocks
        """
        return 'Incidence structure with {} points and {} blocks'.format(
                self.num_points(), self.num_blocks())

    def __str__(self):
        """
        A print method.

        EXAMPLES::

            sage: BD = designs.IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: print BD
            IncidenceStructure<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>
            sage: BD = designs.IncidenceStructure(range(7),[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: print BD
            IncidenceStructure<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>
        """
        return '{}<points={}, blocks={}>'.format(
                self._name, self._points, self._blocks)

    def __eq__(self, other):
        """
        Returns true if their points and blocks are equal (resp.).

        EXAMPLES::

            sage: BD1 = designs.BlockDesign(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: M = BD1.incidence_matrix()
            sage: BD2 = designs.BlockDesign(incidence_matrix=M)
            sage: BD1 == BD2
            True
        """
        if not isinstance(other, IncidenceStructure):
            return False
        return self._points == other._points and self._blocks == other._blocks

    def __ne__(self, other):
        r"""
        Difference test.

        EXAMPLES::

            sage: BD1 = designs.BlockDesign(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: M = BD1.incidence_matrix()
            sage: BD2 = designs.BlockDesign(incidence_matrix=M)
            sage: BD1 != BD2
            False
        """
        return not self.__eq__(other)

    def points(self):
        r"""
        Return the list of points.

        EXAMPLES::

            sage: designs.IncidenceStructure(3, [[0,1],[0,2]]).points()
            [0, 1, 2]
        """
        from copy import deepcopy
        return deepcopy(self._points)

    def num_points(self):
        r"""
        The number of points in that design.

        EXAMPLES::

            sage: designs.DesarguesianProjectivePlaneDesign(2).num_points()
            7
            sage: B = designs.IncidenceStructure(4, [[0,1],[0,2],[0,3],[1,2], [1,2,3]])
            sage: B.num_points()
            4
        """
        return len(self._points)

    def num_blocks(self):
        r"""
        The number of blocks.

        EXAMPLES::

            sage: designs.DesarguesianProjectivePlaneDesign(2).num_blocks()
            7
            sage: B = designs.BlockDesign(4, [[0,1],[0,2],[0,3],[1,2], [1,2,3]])
            sage: B.num_blocks()
            5
        """
        return len(self._blocks)

    def blocks(self):
        """
        Return the list of blocks.

        EXAMPLES::

            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.blocks()
            [[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]

        It is safe to modify the output::

            sage: b = BD.blocks()
            sage: b[0][0] = 18
            sage: BD.blocks()[0]
            [0, 1, 2]
        """
        from copy import deepcopy
        return deepcopy(self._blocks)

    def block_sizes(self):
        r"""
        Return the set of block sizes.

        EXAMPLES::

            sage: BD = designs.BlockDesign(8, [[0,1,3],[1,4,5,6],[1,2],[5,6,7]])
            sage: BD.block_sizes()
            frozenset([2, 3, 4])
            sage: BD = designs.BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.block_sizes()
            frozenset([3])
        """
        return self._block_sizes

    def is_connected(self):
        r"""
        Test whether the design is connected.

        EXAMPLES::

            sage: designs.IncidenceStructure(3, [[0,1],[0,2]]).is_connected()
            True
            sage: designs.IncidenceStructure(4, [[0,1],[2,3]]).is_connected()
            False
        """
        return self.incidence_graph().is_connected()

    def is_simple(self):
        r"""
        Test whether this design is simple (i.e. no repeated block).

        EXAMPLES::

            sage: designs.IncidenceStructure(3, [[0,1],[1,2],[0,2]]).is_simple()
            True
            sage: designs.IncidenceStructure(3, [[0],[0]]).is_simple()
            False
        """
        B = self._blocks
        for b in B:
            if B.count(b) > 1:
                return False
        return True

    def _gap_(self):
        """
        Return the GAP string describing the design.

        EXAMPLES::

            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD._gap_()
            'BlockDesign(7,[[1, 2, 3], [1, 4, 5], [1, 6, 7], [2, 4, 6], [2, 5, 7], [3, 4, 7], [3, 5, 6]])'
        """
        B = self.blocks()
        v = self.num_points()
        gB = []
        for b in B:
            gB.append([x+1 for x in b])
        return "BlockDesign("+str(v)+","+str(gB)+")"

    def incidence_matrix(self):
        """
        Return the incidence matrix `A` of the design. A is a `(v \times b)`
        matrix defined by: ``A[i,j] = 1`` if ``i`` is in block ``B_j`` and 0
        otherwise.

        EXAMPLES::

            sage: BD = designs.BlockDesign(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.block_sizes()
            frozenset([3])
            sage: BD.incidence_matrix()
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [1 0 0 0 0 1 1]
            [0 1 0 1 0 1 0]
            [0 1 0 0 1 0 1]
            [0 0 1 1 0 0 1]
            [0 0 1 0 1 1 0]
        """
        if self._incidence_matrix is None:
            from sage.matrix.constructor import Matrix
            from sage.rings.all import ZZ
            A = Matrix(ZZ, self.num_points(), self.num_blocks(), sparse=True)
            for j, b in enumerate(self.blocks()):
                for i in b:
                    A[i, j] = 1
            self._incidence_matrix = A
        return self._incidence_matrix

    def incidence_graph(self):
        """
        Returns the incidence graph of the design, where the incidence
        matrix of the design is the adjacency matrix of the graph.

        EXAMPLE::

            sage: BD = designs.BlockDesign(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.incidence_graph()
            Bipartite graph on 14 vertices
            sage: A = BD.incidence_matrix()
            sage: Graph(block_matrix([[A*0,A],[A.transpose(),A*0]])) == BD.incidence_graph()
            True

        REFERENCE:

        - Sage Reference Manual on Graphs
        """
        from sage.graphs.bipartite_graph import BipartiteGraph
        A = self.incidence_matrix()
        return BipartiteGraph(A)

    @cached_method
    def _degree_iterator(self):
        r"""
        Manage cached degree computation.

        EXAMPLES::

            sage: I = designs.IncidenceStructure(4, [[0,1,2],[0,1,3],[0,2,3],[1,2,3]])
            sage: I.degrees(None)      # indirect doctest
            [frozenset([4]), frozenset([3]), frozenset([2]), frozenset([1])]

            sage: I = designs.IncidenceStructure(4, [[0,1],[1,2],[0,1,2],[0,1,3]])
            sage: I.degrees(None)      # indirect doctest
            [frozenset([4]), frozenset([1, 2, 3, 4]), frozenset([0, 1, 2, 3]), frozenset([0, 1])]
        """
        points = frozenset(self._points)
        p_to_blocks = {p:[] for p in points}
        for i,b in enumerate(self._blocks):
            for j in b:
                p_to_blocks[j].append(i)
        for p in self._points:
            p_to_blocks[p] = frozenset(p_to_blocks[p])
        s_to_blocks = {}
        for i in p_to_blocks:
            s_to_blocks[frozenset([i])] = frozenset(p_to_blocks[i])

        no_empty = True
        yield frozenset([self.num_blocks()])
        yield frozenset(len(s) for s in p_to_blocks.itervalues())
        k = max(self._block_sizes)
        for i in xrange(2,k+1):
            s_to_blocks2 = {}
            for s1,B1 in s_to_blocks.iteritems():
                for p in points - s1:
                    B2 = p_to_blocks[p] & B1
                    if B2:
                        s_to_blocks2[s1.union([p])] = B2
                    else:
                        no_empty = False
            s_to_blocks = s_to_blocks2
            #print i
            #for s,b in s_to_blocks.iteritems():
            #    print "  {}: {}".format(sorted(s), sorted(b))
            d = frozenset(len(b) for b in s_to_blocks.itervalues())
            if not no_empty:
                d = d.union([0])
            yield d

    def degrees(self, t):
        r"""
        Return the set of ``t``-degrees of ``self``.

        Given a subset of `t` points, its *degree* is the number of blocks that
        contain it. Note that for `t=0` it is the number of blocks.

        EXAMPLES::

            sage: BD = designs.BlockDesign(8, [[0,1,3],[1,4,5,6],[1,2],[5,6,7]])
            sage: BD.points()
            [0, 1, 2, 3, 4, 5, 6, 7]
            sage: BD.degrees(0)
            frozenset([4])
            sage: BD.degrees(1)
            frozenset([1, 2, 3])
            sage: BD.degrees(2)
            frozenset([0, 1, 2])
            sage: BD.degrees(3)
            frozenset([0, 1])
            sage: BD.degrees(4)
            frozenset([0, 1])
            sage: BD.degrees(5)
            frozenset([])

            sage: P = designs.DesarguesianProjectivePlaneDesign(2)
            sage: P.degrees(2)
            frozenset([1])
            sage: P.degrees(3)
            frozenset([0, 1])
        """
        if t is None:
            for d in self._degree_iterator():
                self._degrees.append(d)
            return self._degrees[:]

        if t > max(self._block_sizes):
            return frozenset()

        while len(self._degrees) <= t:
            self._degrees.append(self._degree_iterator().next())
        return self._degrees[t]

    @cached_method
    def t_design_parameters(self):
        r"""
        Return a 4-tuple `(t,v,k,l)` such that the design is a `t-(v,k,l)`
        design with `t` maximum.

        If the design is not a `t`-design then an exception occurs.

        EXAMPLES::

        Many affine geometry design are examples of `t`-designs::

            sage: A = designs.AffineGeometryDesign(3, 1, GF(2))
            sage: A.t_design_parameters()
            (2, 8, 2, 1)
            sage: A = designs.AffineGeometryDesign(4, 2, GF(2))
            sage: A.t_design_parameters()
            (3, 16, 4, 1)

        Steiner triple and quadruple systems are other names for `2-(v,3,1)` and
        `3-(v,4,1)` designs::

            sage: S3_9 = designs.steiner_triple_system(9)
            sage: S3_9.t_design_parameters()
            (2, 9, 3, 1)

            sage: blocks = designs.steiner_quadruple_system(8)
            sage: S8 = designs.IncidenceStructure(8, blocks)
            sage: S8.t_design_parameters()
            (3, 8, 4, 1)

            sage: blocks = designs.steiner_quadruple_system(14)
            sage: S14 = designs.IncidenceStructure(14, blocks)
            sage: S14.t_design_parameters()
            (3, 14, 4, 1)

        Some examples of Witt designs that need the gap database::

            sage: BD = designs.WittDesign(9)     # optional - gap_packages
            sage: BD.t_design_parameters()       # optional - gap_packages
            (2, 9, 3, 1)
            sage: W12 = designs.WittDesign(12)    # optional - gap_packages
            sage: W12.t_design_parameters()       # optional - gap_packages
            (5, 12, 6, 1)

        Further examples::

            sage: D = designs.IncidenceStructure(4,[[],[]])
            sage: D.t_design_parameters()
            (0, 4, 0, 2)

            sage: D = designs.IncidenceStructure(4, [[0,1],[0,2],[0,3]])
            sage: D.t_design_parameters()
            (0, 4, 2, 3)

            sage: D = designs.IncidenceStructure(4, [[0],[1],[2],[3]])
            sage: D.t_design_parameters()
            (1, 4, 1, 1)

            sage: D = designs.IncidenceStructure(4,[[0,1],[2,3]])
            sage: D.t_design_parameters()
            (1, 4, 2, 1)

            sage: D = designs.IncidenceStructure(4, [range(4)])
            sage: D.t_design_parameters()
            (4, 4, 4, 1)

        Bad cases::

            sage: designs.IncidenceStructure(2, []).t_design_parameters()
            Traceback (most recent call last):
            ...
            ValueError: there should be at least one point and one block in a t-design

            sage: designs.IncidenceStructure(2, [[0],[0,1]]).t_design_parameters()
            Traceback (most recent call last):
            ...
            ValueError: in a t-design all blocks have the same size
        """
        from sage.rings.arith import binomial

        v = self.num_points()
        b = self.num_blocks()

        if v == 0 or b == 0:
            raise ValueError("there should be at least one point and one block in a t-design")

        k = len(self._blocks[0])
        if any(len(block) != k for block in self._blocks):
            raise ValueError("in a t-design all blocks have the same size")

        # here we already know that we are a 0-design
        deg = self.degrees(1)
        if len(deg) != 1:
            return (0,v,k,b)
        l, = deg

        # here we know that it is a 1-design
        if k == v:
            # this is the design (X, {X,X,X,...})
            return (v,v,v,b)

        t = 1

        for tt in range(2,k+1):
            # is lambda an integer?
            if (b*binomial(k,tt)) % binomial(v,tt) != 0:
                break

            ll = b*binomial(k,tt) // binomial(v,tt)
            if len(self.degrees(tt)) == 1:
                t = tt
                l = ll
            else:
                break

        return (t,v,k,l)

    def is_t_design(self, t=None, v=None, k=None, l=None):
        """
        Test whether ``self`` is a ``t-(v,k,l)` design.

        A `t-(v,k,\lambda)` (sometimes called `t`-design for short) is a block
        design in which:
        - the underlying set has cardinality `v`
        - the blocks have size `k`
        - each `t`-subset of points is covered by `\lambda` blocks

        INPUT:

        - ``t``, ``v``, ``k``, ``l`` -- optional parameters

        - ``verbose`` (boolean) -- prints useful information when the answer is
          negative.

        EXAMPLES::

            sage: fano_blocks = [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]]
            sage: BD = designs.BlockDesign(7, fano_blocks)
            sage: BD.is_t_design()
            True
            sage: BD.t_design_parameters()
            (2, 7, 3, 1)
            sage: BD.is_t_design(2, 7, 3, 1)
            True
            sage: BD.is_t_design(1, 7, 3, 3)
            True
            sage: BD.is_t_design(0, 7, 3, 7)
            True

            sage: BD.is_t_design(0,6,3,7) or BD.is_t_design(0,7,4,7) or BD.is_t_design(0,7,3,8)
            False

            sage: BD = designs.AffineGeometryDesign(3, 1, GF(2))
            sage: BD.is_t_design(2)
            True

        Steiner triple and quadruple systems are other names for `2-(v,3,1)` and
        `3-(v,4,1)` designs::

            sage: S3_9 = designs.steiner_triple_system(9)
            sage: S3_9.is_t_design(2,9,3,1)
            True

            sage: blocks = designs.steiner_quadruple_system(8)
            sage: S4_8 = designs.IncidenceStructure(8, blocks)
            sage: S4_8.is_t_design(3,8,4,1)
            True

            sage: blocks = designs.steiner_quadruple_system(14)
            sage: S4_14 = designs.IncidenceStructure(14, blocks)
            sage: S4_14.is_t_design(3,14,4,1)
            True

        Some examples of Witt designs that need the gap database::

            sage: BD = designs.WittDesign(9)         # optional - gap_packages
            sage: BD.is_t_design(2,9,3,1)            # optional - gap_packages
            True
            sage: W12 = designs.WittDesign(12)       # optional - gap_packages
            sage: W12.is_t_design(5,12,6,1)          # optional - gap_packages
            True
            sage: W12.is_t_design(4)                 # optional - gap_packages
            True

        TESTS::

            sage: blocks = designs.steiner_quadruple_system(8)
            sage: S4_8 = designs.IncidenceStructure(8, blocks)
            sage: R = range(15)
            sage: [(v,k,l) for v in R for k in R for l in R if S4_8.is_t_design(3,v,k,l)]
            [(8, 4, 1)]
            sage: [(v,k,l) for v in R for k in R for l in R if S4_8.is_t_design(2,v,k,l)]
            [(8, 4, 3)]
            sage: [(v,k,l) for v in R for k in R for l in R if S4_8.is_t_design(1,v,k,l)]
            [(8, 4, 7)]
            sage: [(v,k,l) for v in R for k in R for l in R if S4_8.is_t_design(0,v,k,l)]
            [(8, 4, 14)]
        """
        try:
            tt,vv,kk,ll = self.t_design_parameters()
        except ValueError:
            return False

        if t is None:
            t = tt

        if t > tt:
            return False

        # Handbook of combinatorial design theorem II.4.8: a t-(v,k,lambda) is
        # also a s-(v,k,lambda_s) design with:
        # lambda_s = lambda binomial(v-s,t-s) / binomial(k-s,t-s)

        from sage.rings.arith import binomial
        return ((v is None or v == vv) and
                (k is None or k == kk) and
                (l is None or l == ll * binomial(v-t,tt-t) // binomial(k-t,tt-t)))

    #####################
    # real computations #
    #####################

    def dual(self, algorithm=None):
        """
        Returns the dual of the incidence structure.

        INPUT:

        - ``algorithm`` -- whether to use Sage's implementation
          (``algorithm=None``, default) or use GAP's (``algorithm="gap"``).

          .. NOTE::

              The ``algorithm="gap"`` option requires GAP's Design package
              (included in the gap_packages Sage spkg).

        EXAMPLES:

        The dual of a projective plane is a projective plane::

            sage: PP = designs.DesarguesianProjectivePlaneDesign(4)
            sage: PP.dual().t_design_parameters()
            (2, 21, 5, 1)

        TESTS::

            sage: D = designs.IncidenceStructure(4, [[0,2],[1,2,3],[2,3]])
            sage: D
            Incidence structure with 4 points and 3 blocks
            sage: D.dual()
            Incidence structure with 3 points and 4 blocks
            sage: print D.dual(algorithm="gap")       # optional - gap_packages
            IncidenceStructure<points=[0, 1, 2], blocks=[[0], [0, 1, 2], [1], [1, 2]]>
            sage: blocks = [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]]
            sage: BD = designs.IncidenceStructure(7, blocks, name="FanoPlane");
            sage: BD
            Incidence structure with 7 points and 7 blocks
            sage: print BD.dual(algorithm="gap")         # optional - gap_packages
            IncidenceStructure<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>
            sage: BD.dual()
            Incidence structure with 7 points and 7 blocks

        REFERENCE:

        - Soicher, Leonard, Design package manual, available at
          http://www.gap-system.org/Manuals/pkg/design/htm/CHAP003.htm
        """
        if algorithm == "gap":
            from sage.interfaces.gap import gap
            gap.load_package("design")
            gD = self._gap_()
            gap.eval("DD:=DualBlockDesign("+gD+")")
            v = eval(gap.eval("DD.v"))
            gblcks = eval(gap.eval("DD.blocks"))
            gB = []
            for b in gblcks:
                gB.append([x-1 for x in b])
            return IncidenceStructure(range(v), gB, name=None, check=False)
        else:
            from block_design import BlockDesign
            return IncidenceStructure(
                          incidence_matrix=self.incidence_matrix().transpose(),
                          check=False)

    def automorphism_group(self):
        r"""
        Returns the subgroup of the automorphism group of the incidence graph
        which respects the P B partition. It is (isomorphic to) the automorphism
        group of the block design, although the degrees differ.

        EXAMPLES::

            sage: P = designs.DesarguesianProjectivePlaneDesign(2); P
            Incidence structure with 7 points and 7 blocks
            sage: G = P.automorphism_group()
            sage: G.is_isomorphic(PGL(3,2))
            True
            sage: G
            Permutation Group with generators [(2,3)(4,5), (2,4)(3,5), (1,2)(4,6), (0,1)(4,5)]

        A non self-dual example::

            sage: IS = designs.IncidenceStructure(range(4), [[0,1,2,3],[1,2,3]])
            sage: IS.automorphism_group().cardinality()
            6
            sage: IS.dual().automorphism_group().cardinality()
            1
        """
        from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        M1 = self.incidence_matrix().transpose()
        M2 = MatrixStruct(M1)
        M2.run()
        gens = M2.automorphism_group()[0]
        return PermutationGroup(gens, domain=self.points())

    ###############
    # Deprecation #
    ###############

    parameters = deprecated_function_alias(16553, t_design_parameters)
    dual_design = deprecated_function_alias(16553, dual)
    dual_incidence_structure = deprecated_function_alias(16553, dual)

    def block_design_checker(self, t, v, k, lmbda, type=None):
        """
        This method is deprecated and will soon be removed (see :trac:`16553`).
        You could use :meth:`is_t_design` or :meth:`t_design_parameters` instead.

        This is *not* a wrapper for GAP Design's IsBlockDesign. The GAP
        Design function IsBlockDesign
        http://www.gap-system.org/Manuals/pkg/design/htm/CHAP004.htm
        apparently simply checks the record structure and no mathematical
        properties. Instead, the function below checks some necessary (but
        not sufficient) "easy" identities arising from the identity.

        INPUT:

        - ``t`` - the t as in "t-design"

        - ``v`` - the number of points

        - ``k`` - the number of blocks incident to a point

        - ``lmbda`` - each t-tuple of points should be incident with
          lmbda blocks

        - ``type`` - can be 'simple' or 'binary' or 'connected'
          Depending on the option, this wraps IsBinaryBlockDesign,
          IsSimpleBlockDesign, or IsConnectedBlockDesign.

          - Binary: no block has a repeated element.

          - Simple: no block is repeated.

          - Connected: its incidence graph is a connected graph.

        WARNING: This is very fast but can return false positives.

        EXAMPLES::

            sage: BD = designs.IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.t_design_parameters()
            (2, 7, 3, 1)
            sage: BD.block_design_checker(2, 7, 3, 1)
            doctest:...: DeprecationWarning: .block_design_checker(v,t,k,lmbda) is deprecated; please use
            .is_t_design(v,t,k,lmbda) instead
            See http://trac.sagemath.org/16553 for details.
            True

            sage: BD.block_design_checker(2, 7, 3, 1,"binary")
            doctest:1: DeprecationWarning: .block_design_checker(type='binary') is
            deprecated; use .is_binary() instead
            See http://trac.sagemath.org/16553 for details.
            True

            sage: BD.block_design_checker(2, 7, 3, 1,"connected")
            doctest:1: DeprecationWarning: block_design_checker(type='connected') is
            deprecated, please use .is_connected() instead
            See http://trac.sagemath.org/16553 for details.
            True

            sage: BD.block_design_checker(2, 7, 3, 1,"simple")
            doctest:1: DeprecationWarning: .block_design_checker(type='simple')
            is deprecated; all designs here are simple!
            See http://trac.sagemath.org/16553 for details.
            True
        """
        from sage.misc.superseded import deprecation

        ans = self.is_t_design(t,v,k,lmbda)

        if type is None:
            deprecation(16553, ".block_design_checker(v,t,k,lmbda) is deprecated; please use .is_t_design(v,t,k,lmbda) instead")
            return ans

        if type == "binary":
            deprecation(16553, ".block_design_checker(type='binary') is deprecated; use .is_binary() instead")
            return True
        if type == "simple":
            deprecation(16553, ".block_design_checker(type='simple') is deprecated; all designs here are simple!")
            return True
        if type == "connected":
            deprecation(16553, "block_design_checker(type='connected') is deprecated, please use .is_connected() instead")
            return self.incidence_graph().is_connected()
