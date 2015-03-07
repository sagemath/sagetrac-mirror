r"""
Cluster Triangulations

A *cluster triangulation* (see [FominShapiroThurston]_) is a subset
of an ideal triangulation, i.e. a maximal collection of distinct
non-crossing arcs. An ideal triangulation of a surface with marked
points `(S,M)` is associated to a seed of a cluster algebra arising
from `(S,M)`.

.. SEEALSO::

    Cluster triangulations closely interact with
    :class:`~sage.combinat.cluster_algebra_quiver.cluster_seed.ClusterSeed`,
    :class:`~sage.combinat.cluster_algebra_quiver.quiver.ClusterQuiver`

REFERENCES:

[MSW_Positivity]_ [MSW_Bases]_ [MW_MatrixFormulae]_ [FominShapiroThurston]_
"""
from sage.structure.sage_object import SageObject


class ClusterTriangulation(SageObject):
    """
    An initial *ideal triangulation* associated to a surface

    INPUT:

    - ``data`` -- can be any of the following::

        * List of triangles - must be the list of 3-tuples (i.e. edge
          labels of a triangle) from an ideal triangulation (see
          Examples)

        * ClusterTriangulation

    .. TODO::

    - ``data`` -- can also be any of the following::

        * Matrix - a skew-symmetrizable matrix arising from a tagged
          triangulation (todo: not yet implemented)

        * DiGraph - must be the input data for a quiver from a tagged
          triangulation (todo: not yet implemented)

        * List of edges - must be the edge list of a digraph for a
          quiver from a tagged triangulation (todo: not yet
          implemented)

        * Objects that Theodosios Douvropoulos is currently working
          with (not related to cluster algebras from surfaces) (todo:
          not yet implemented)

    EXAMPLES:

    From a list of ideal triangles (forming an ideal triangulation of
    a surface)::

        sage: Triangles = [('nu','gamma','delta'), ('nu','alpha','beta')]
        sage: T = ClusterTriangulation(Triangles)
        sage: ClusterSeed(T).mutation_type()
        ['A', 5]

        sage: once_punctured_square = [('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')]
        sage: T = ClusterTriangulation(once_punctured_square, boundary_edges=['c','f','e','d'])
        sage: ClusterSeed(T).mutation_type()
        ['D', 4]

        sage: once_punctured_square = [('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')]
        sage: T = ClusterTriangulation(once_punctured_square, boundary_edges=[])
        sage: ClusterSeed(T).mutation_type()
        ['D', 8]


    Figure 6 of [MW_MatrixFormulae]_ ::

        sage: annulus22 = [('bd1','tau1','tau2'),('tau2','tau3','bd4'),('tau1','tau4','bd2'),('tau3','bd3','tau4')]
        sage: T = ClusterTriangulation(annulus22, boundary_edges=['bd3','bd2','bd1','bd4'])
        sage: T
        An ideal triangulation associated with cluster algebra of rank 4 with 4 boundary edges
        sage: ClusterSeed(T).mutation_type()
        ['A', [2, 2], 1]
        sage: T.triangulation_dictionary()
        [('tau1', x0),
        ('tau2', x1),
        ('tau3', x2),
        ('tau4', x3),
        ('bd1', b4),
        ('bd2', b5),
        ('bd3', b6),
        ('bd4', b7)]

    Figure 15 of [MSW_Positivity]_ ::

        sage: thrice_punctured_square = [(2,2,1), (1,3,11), (3,12,4), (4,5,14), (5,6,10), (6,7,9), (8,10,9), (7,13,8)]
        sage: T = ClusterTriangulation(thrice_punctured_square, boundary_edges=[14,12,13,11])
        sage: T
        An ideal triangulation associated with cluster algebra of rank 10 with 4 boundary edges
        sage: ClusterSeed(T).mutation_type()
        'undetermined finite mutation type'
        sage: T.triangulation_dictionary()
        [(1, x0*x1),
        (2, x1),
        (3, x2),
        (4, x3),
        (5, x4),
        (6, x5),
        (7, x6),
        (8, x7),
        (9, x8),
        (10, x9),
        (11, b10),
        (12, b11),
        (13, b12),
        (14, b13)]

        sage: twice_punctured_bigon = [('e','d','a'), ('a','r','b'), ('r','d','g'), ('g','n','b')]
        sage: T = ClusterTriangulation(twice_punctured_bigon, boundary_edges=['e','n'])
        sage: T
        An ideal triangulation associated with cluster algebra of rank 5 with 2 boundary edges
        sage: ClusterSeed(T).mutation_type()
        'undetermined finite mutation type'

        sage: T = ClusterTriangulation ( [[4, 5, 1], [4, 3, 2], [3, 7, 2], [2, 1, 6], [1, 4, 5]], boundary_edges=[1])
        sage: T
        An ideal triangulation associated with cluster algebra of rank 6 with 1 boundary edges

        sage: T = ClusterTriangulation ( [(4, 5, 1), (4, 3, 2), (3, 7, 2), (2, 1, 6), (1, 4, 5)])
        sage: T
        An ideal triangulation associated with cluster algebra of rank 7

        sage: Triangles = [(2,3,11),(2,1,1),(4,3,12),(0,4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,13)]
        sage: T = ClusterTriangulation(Triangles, boundary_edges=[11,12,13,0])
        sage: T
        An ideal triangulation associated with cluster algebra of rank 10 with 4 boundary edges

        sage: once_punctured_torus = ClusterTriangulation([(1,2,3),(3,1,2)])
        sage: S = ClusterSeed(once_punctured_torus).mutation_type()
        sage: S
        'undetermined finite mutation type'

    Triangles do not need to be connected::

        sage: T = ClusterTriangulation([(1,2,3),(3,4,5),(1,5,6),(10,20,30),(30,40,50)])
        sage: ClusterSeed(T).mutation_type()
        [ ['D', 6], ['A', 5] ]
    """
    def __init__(self, data, boundary_edges=None):
        """
        .. TODO::

            See my data for the TestSuite. Is this good enough?

        TESTS::

            sage: CT = ClusterTriangulation([('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')], boundary_edges=['c','d','e','f'])
            sage: TestSuite(CT).run()
        """
        from sage.combinat.cluster_algebra_quiver.surface import remove_duplicate_triangles, _triangulation_to_arrows, _surface_edge_list_to_matrix, _get_user_arc_labels, _get_triangulation_dictionary, _get_user_label_triangulation, _get_weighted_triangulation
        from sage.rings.all import QQ
        from sage.rings.all import FractionField, PolynomialRing
        from copy import copy

        # if data is a list of triangles (forming a triangulation)
        if isinstance(data,list) and \
        all(type(triangle) in [list,tuple] and len(triangle)==3 for triangle in data):
            self._boundary_edges = boundary_edges if boundary_edges else []
            self._arcs = _get_user_arc_labels(data)

            if not set(self._boundary_edges).issubset(self._arcs):
                raise ValueError(boundary_edges, ' are not a subset of ', self._arcs, ' .Optional parameter boundary_edges must be a list of edges that are part of the triangulation')
            data = remove_duplicate_triangles (data, boundary_edges)

            self._n = len(self._arcs) - len(self._boundary_edges) #if boundary_edges else len(self._arcs)
            self._triangles = data

            all_arrows = _triangulation_to_arrows(self._triangles)
            M = _surface_edge_list_to_matrix (all_arrows, self._arcs, self._boundary_edges, self._n)
            self._M = M[:self._n,:self._n] # In this implementation, we ignore the boundary edges (TODO)
            self._is_cluster_algebra = True
            self._description = 'An ideal triangulation associated with cluster algebra of rank %d'  %(self._n)
            if boundary_edges:
                self._description += ' with %d boundary edges' %(len(self._boundary_edges))
            self._R = FractionField(PolynomialRing(QQ,['x%s'%i for i in range(0,self._n)]+['b%s'%i for i in range(self._n,len(self._boundary_edges)+self._n)]))
            self._cluster = list(self._R.gens()[0:self._n])
            self._boundary_edges_vars = list(self._R.gens()[self._n:]) if boundary_edges else []
            self._triangulation_dictionary = _get_triangulation_dictionary (self._triangles, self._cluster, self._boundary_edges, self._boundary_edges_vars)
            self._triangulation = _get_user_label_triangulation(self._triangles)
            self._weighted_triangulation = _get_weighted_triangulation (self._triangles, self._triangulation_dictionary)

        # Construct data from a cluster triangulation
        elif isinstance(data, ClusterTriangulation):
            self._boundary_edges = copy(data._boundary_edges)
            self._arcs = copy(data._arcs)
            self._n = data._n
            self._triangles = copy( data._triangles )
            self._M = data._M[:data._n,:data._n]
            self._is_cluster_algebra = data._is_cluster_algebra
            self._description = copy( data._description )
            self._R = copy(data._R)
            self._cluster = copy(data._cluster)
            self._boundary_edges_vars = copy(data._boundary_edges_vars)
            self._triangulation_dictionary = copy(data._triangulation_dictionary)
            self._triangulation = copy(data._triangulation)
            self._weighted_triangulation = copy(data._weighted_triangulation)

        else:
            raise ValueError('Input must be a list of three-tuples')

    def __eq__(self, other):
        r"""
        Returns True iff ``self`` represent the same cluster triangulation as ``other``.

        EXAMPLES::

            sage: CT = ClusterTriangulation([('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')], boundary_edges=['c','d','e','f'])
            sage: CT1 = ClusterTriangulation([('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')])
            sage: CT.__eq__(CT1)
            False

            sage: CT2 = ClusterTriangulation([('aa','dd','cc'), ('aa','ll','bb'), ('rr','rr','ll'),('bb','ff','ee')], boundary_edges=['ee','dd','cc','ff'])
            sage: CT.__eq__(CT2)
            True
        """
        from sage.combinat.cluster_algebra_quiver.surface import _rearrange_triangle_small_first

        Tself, Tother = [], []
        for triangle in self._weighted_triangulation:
            Tself.append(_rearrange_triangle_small_first(triangle))
        for triangle in other._weighted_triangulation:
            Tother.append(_rearrange_triangle_small_first(triangle))
        Tself.sort()
        Tother.sort()
        return isinstance(other, ClusterTriangulation) and Tself == Tother

    def _repr_(self):
        """
        Return the description of ``self``.

        EXAMPLES::

            sage: T = ClusterTriangulation([(4, 5, 1), (4, 3, 2), (3, 7, 2), (2, 1, 6), (1, 4, 5)])
            sage: T._repr_()
            'An ideal triangulation associated with cluster algebra of rank 7'
        """
        name = self._description
        return name

    def triangles(self):
        """
        Return the underlying triangulation of ``self``.

        EXAMPLES::

            sage: Triangles = [(1,4,7),(1,2,5),(2,0,3),(0,6,3),(6,3,0)]
            sage: T = ClusterTriangulation(Triangles)
            sage: T.triangles()
            [(1, 4, 7), (1, 2, 5), (2, 0, 3), (0, 6, 3)]
        """
        return self._triangles

    def b_matrix(self):
        """
        Return the B-matrix of the coefficient-free cluster.

        The conventions for B-matrices are the opposite of
        [FominShapiroThurston]_.

        EXAMPLES:

        Twice-punctured monogon with 3 (non-ordinary) ideal triangles (affine D)::

            sage: twice_punctured_monogon = [[4,5,5],[1,2,3],[2,4,3]]
            sage: T = ClusterTriangulation(twice_punctured_monogon)
            sage: B = T.b_matrix()
            sage: B
            [ 0  1 -1  0  0]
            [-1  0  0  1  1]
            [ 1  0  0 -1 -1]
            [ 0 -1  1  0  0]
            [ 0 -1  1  0  0]
            sage: twice_punctured_monogon_mu2 = [(4,5,5),(2,3,3),(1,4,2)]

        2 self-folded triangles and 1 triangle with one vertex
        (affine D) Figure 9 (right) of [FominShapiroThurston]_ ::

            sage: Tmu2 = ClusterTriangulation(twice_punctured_monogon_mu2)
            sage: Bmu2 = Tmu2.b_matrix()
            sage: Bmu2
            [ 0 -1 -1  1  1]
            [ 1  0  0 -1 -1]
            [ 1  0  0 -1 -1]
            [-1  1  1  0  0]
            [-1  1  1  0  0]
            sage: B.mutate(1)
            sage: Bmu2 == B
            True

            sage: Qmu2 = ClusterQuiver(Bmu2)
            sage: Qmu2.mutation_type()
            'undetermined finite mutation type'

        Figure 10 (bottom) of [FominShapiroThurston]_ ::

            sage: twice_punctured_monogon = [(1,1,2), (4,4,3), ('boundary',2,3)]
            sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=['boundary'])
            sage: B = T.b_matrix()
            sage: B
            [ 0  0  1  1]
            [ 0  0  1  1]
            [-1 -1  0  0]
            [-1 -1  0  0]

        Figure 10 (top) of [FominShapiroThurston]_ ::

            sage: twice_punctured_monogon_mu3 = [(1,1,2), (2,4,3), (3,4,'boundary')]
            sage: Tmu3 = ClusterTriangulation(twice_punctured_monogon_mu3, boundary_edges=['boundary'])
            sage: Bmu3 = Tmu3.b_matrix()
            sage: Bmu3
            [ 0  0 -1  1]
            [ 0  0 -1  1]
            [ 1  1  0  0]
            [-1 -1  0  0]
            sage: B.mutate(2)
            sage: Bmu3 == B
            True

        Four-punctured sphere with 3 self-folded triangle::

            sage: four_punc_sphere = [(1,0,5),(3,5,4),(2,1,3),(0,2,4)]
            sage: T = ClusterTriangulation(four_punc_sphere)
            sage: B = T.b_matrix()
            sage: Tmu0 = ClusterTriangulation([(0,1,2),(3,5,4),(2,1,3),(0,4,5)])
            sage: B.mutate(0)
            sage: B == Tmu0.b_matrix()
            True
            sage: Tmu04 = ClusterTriangulation([(0,1,2),(4,5,5),(2,1,3),(3,4,0)])
            sage: B.mutate(4)
            sage: B == Tmu04.b_matrix()
            True
            sage: Tmu040 = ClusterTriangulation([(1,0,4),(4,5,5),(2,1,3),(3,0,2)])
            sage: B.mutate(0)
            sage: B == Tmu040.b_matrix()
            True
            sage: Tmu0402 = Tmu040.mutate(2,inplace=False)
            sage: Tmu0402 == ClusterTriangulation([(1,0,4),(4,5,5),(1,2,0),(3,2,3)])
            True
        """
        return self._M

    def mutate(self, sequence, inplace=True, user_labels=True):
        """
        Mutates ``self`` at a diagonal or a sequence of diagonals

        See :meth:`ClusterSeed.mutate`

        INPUT:

        - ``sequence`` -- a diagonal label (if user_labels is True)
            or a diagonal position (if user_labels is not True) of self
            or an iterator of diagonal labels/positions of self.
        - ``inplace`` -- (default: True) if False, the result is returned,
            otherwise ``self`` is modified.
        - ``user_labels`` -- (default: True) if False, user should input position/s
            of diagonal/s in input ``sequence``

        EXAMPLES:

        Twice-punctured monogon with 3 (non-ordinary) ideal triangles (affine D)::

            sage: T = ClusterTriangulation([('i1','i4','i2'),('i3','i4','i3'),('i2','i0','i1')])

        Two self-folded triangles and 1 triangle with one vertex (affine D)::

            sage: Tmu2 = ClusterTriangulation([('j1','j1','j2'),('j3','j4','j3'),('j2','j4','j0')])
            sage: T.mutate('i2',inplace=False).b_matrix() == Tmu2.b_matrix()
            True
            sage: Tmu2.mutate(2, user_labels=False)
            sage: T.b_matrix() == Tmu2.b_matrix()
            True

        Four-punctured sphere::

            sage: four_punc_sphere = [('i1','i0','i5'),('i3','i5','i4'),('i2','i1','i3'),('i0','i2','i4')]
            sage: T = ClusterTriangulation(four_punc_sphere)
            sage: T.mutate([0,4,0,2],user_labels=False)
            sage: Tmu0402 = ClusterTriangulation([(1,0,4),(4,5,5),(1,2,0),(3,2,3)])
            sage: T.b_matrix() == Tmu0402.b_matrix()
            True
            sage: Tmu04020 = ClusterTriangulation([(1,0,1),(4,5,5),(4,0,2),(3,2,3)])
            sage: T.mutate('i0',inplace=False).b_matrix() == Tmu04020.b_matrix()
            True
            sage: T.mutate('i0',inplace=False) == Tmu04020
            True
        """
        from sage.combinat.cluster_algebra_quiver.surface import _triangles_mutate, \
        _get_triangulation_dictionary, _get_user_label_triangulation, _get_weighted_triangulation

        if inplace:
            ct = self
        else:
            ct = ClusterTriangulation( self )

        n = ct._n
        V = range(n)
        arcs = self._arcs

        if user_labels:
            if sequence in arcs:
                seq = [sequence]
            else:
                seq = sequence
            if isinstance(seq, tuple):
                seq = list( seq )
        else:
            if sequence in V:
                seq = [sequence]
            else:
                seq = sequence
            if isinstance(seq, tuple):
                seq = list( seq )

        if not isinstance(seq, list):
            raise ValueError('The ideal triangulation can only be mutated at a diagonal or at a sequence of diagonals')
        if not isinstance(inplace, bool):
            raise ValueError('The second parameter must be boolean. To mutate at a sequence of length 2, input it as a list.')
        if not isinstance(user_labels, bool):
            raise ValueError('The third parameter must be boolean. To mutate at a sequence of length 3, input it as a list.')

        if user_labels:
            if any( a not in arcs for a in seq ):
                a = filter( lambda a: a not in arcs, seq )[0]
                raise ValueError(str( a ) + ' is not an arc in the triangulation')
            if any( b in ct._boundary_edges for b in seq ):
                b = filter( lambda b: b in ct._boundary_edges, seq )[0]
                raise ValueError('The triangulation cannot be mutated at a boundary edge ' + str( b ))
        else:
            if any( v not in V for v in seq ):
                v = filter( lambda v: v not in V, seq )[0]
                raise ValueError('The triangulation cannot be mutated at the diagonal position ' + str( v ))

        for diagonal in seq:
            if user_labels:
                pos = ct.get_edge_position(diagonal)
                ct._triangles = _triangles_mutate(ct._triangles, diagonal)
                ct._M.mutate(pos)
            else:
                pos = diagonal
                diagonal_label = ct._arcs[pos]
                ct._triangles = _triangles_mutate(ct._triangles, diagonal_label)
                ct._M.mutate(pos)

        ct._triangulation_dictionary = _get_triangulation_dictionary (ct._triangles, ct._cluster, ct._boundary_edges, ct._boundary_edges_vars)
        ct._triangulation = _get_user_label_triangulation(ct._triangles)
        ct._weighted_triangulation = _get_weighted_triangulation (ct._triangles, ct._triangulation_dictionary)

        if not inplace:
            return ct

    def arcs(self):
        """
        Return the sorted list of labels of diagonals and boundary edges of
        ``self`` given by user.

        EXAMPLES::

            sage: annulus22 = [('bd1','tau1','tau2'),('tau2','tau3','bd4'),('tau1','tau4','bd2'),('tau3','bd3','tau4')]
            sage: T = ClusterTriangulation(annulus22, boundary_edges=['bd3','bd2','bd1','bd4'])
            sage: T.arcs()
            ['bd1', 'bd2', 'bd3', 'bd4', 'tau1', 'tau2', 'tau3', 'tau4']
        """
        return self._arcs

    def boundary_edges(self):
        """
        Return the labels of boundary edges (not diagonals) of
        ``self`` given by user.

        EXAMPLES::

            sage: annulus22 = [('bd1','tau1','tau2'),('tau2','tau3','bd4'),('tau1','tau4','bd2'),('tau3','bd3','tau4')]
            sage: T = ClusterTriangulation(annulus22, boundary_edges=['bd3','bd2','bd1','bd4'])
            sage: T.boundary_edges()
            ['bd1', 'bd2', 'bd3', 'bd4']
        """
        return self._boundary_edges

    def cluster(self):
        """
        Return the cluster seed list (i.e. [x_0, x_1, ..., x_{n-1}])
        corresponding to ``self``.

        EXAMPLES::

            sage: once_punctured_square = [('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')]
            sage: T = ClusterTriangulation(once_punctured_square, boundary_edges=['c','f','e','d'])
            sage: T.cluster()
            [x0, x1, x2, x3]
        """
        return self._cluster

    def boundary_edges_vars(self):
        """
        Return the cluster seed list (i.e. [x_0, x_1, ..., x_{n-1}])
        corresponding to ``self``.

        EXAMPLES::

            sage: once_punctured_square = [('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')]
            sage: T = ClusterTriangulation(once_punctured_square, boundary_edges=['c','f','e','d'])
            sage: T.boundary_edges_vars()
            [b4, b5, b6, b7]
        """
        return self._boundary_edges_vars

    def triangulation_dictionary(self):
        """
        Return the correspondence between user-given labels (integers)
        and variables of ``self``.

        EXAMPLES:

         2 self-folded triangles and 1 triangle with one vertex (affine D),
         Figure 10 (bottom) of [FominShapiroThurston]_::

            sage: twice_punctured_monogon = [(4,5,5),(2,3,3),(1,4,2)]
            sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=[1])
            sage: T.triangulation_dictionary()
            [(2, x0*x1), (3, x1), (4, x2*x3), (5, x3), (1, b4)]
        """
        return self._triangulation_dictionary

    def triangulation(self):
        """
        Return the list of triangles of ``self`` where a self-folded
        triangle (r,r,ell) is replaced by (r, 'counterclockwise'), (r,
        'clockwise'), ell).

        EXAMPLES::

            sage: Triangles = [(1,4,7),(1,2,5),(2,0,3),(0,6,3)]
            sage: CT = ClusterTriangulation(Triangles)
            sage: CT.triangulation()
            [(1, 4, 7), (1, 2, 5), (2, 0, 3), (0, 6, 3)]

        Figure 10 (bottom) of [FominShapiroThurston]_::

            sage: twice_punctured_monogon = [(4,5,5),(2,3,3),(1,4,2)]
            sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=[1]) # 2 self-folded triangles and 1 triangle with one vertex (affine D)
            sage: T.triangulation()
            [((5, 'counterclockwise'), (5, 'clockwise'), 4),
            ((3, 'counterclockwise'), (3, 'clockwise'), 2),
            (1, 4, 2)]
        """
        return self._triangulation

    def weighted_triangulation(self):
        """
        Return the list of triangles of ``self`` where user-given
        labels are replaced by the variables corresponding to them.

        EXAMPLES::

            sage: Triangles = [(1,4,7),(1,2,5),(2,0,3),(0,6,3)]
            sage: CT = ClusterTriangulation(Triangles)
            sage: CT.weighted_triangulation()
            [(x1, x4, x7), (x1, x2, x5), (x2, x0, x3), (x0, x6, x3)]

        Two self-folded triangles and 1 triangle with one vertex (affine D). See Figure 10 (bottom) of [FominShapiroThurston]_::

            sage: twice_punctured_monogon = [(4,5,5),(2,3,3),(1,4,2)]
            sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=[1])
            sage: T.weighted_triangulation()
            [((x3, 'counterclockwise'), (x3, 'clockwise'), x2*x3),
            ((x1, 'counterclockwise'), (x1, 'clockwise'), x0*x1),
            (b4, x2*x3, x0*x1)]
        """
        return self._weighted_triangulation

    def get_edge_var(self,a):
        """
        Return the variable corresponding to the label (given by user)
        of an arc of boundary edge.

        .. SEEALSO::

            :meth:`get_edge_position`

        EXAMPLES::

            sage: T = ClusterTriangulation([(1,7,4),(1,5,2),(6,0,3),(2,3,0),(0,3,6),[7,4,1]], boundary_edges=[4,5,6,7])
            sage: T.get_edge_var(0)
            x0
            sage: T.get_edge_var(6)
            b6
        """
        from sage.combinat.cluster_algebra_quiver.surface import _get_weighted_edge
        return _get_weighted_edge(a,self._triangulation_dictionary)

    def get_edge_position(self,a):
        """
        Return the position of the variable corresponding to the label (given by user)
        of an arc or a boundary edge.

        .. SEEALSO::

            :meth:`get_edge_var`

        EXAMPLES::

            sage: T = ClusterTriangulation([('i1','i7','i4'),('i1','i5','i2'),\
            ('i6','i0','i3'),('i2','i3','i0'),('i0','i3','i6'),['i7','i4','i1']], \
            boundary_edges=['i4','i5','i6','i7'])
            sage: T.get_edge_position('i0')
            0
            sage: T.get_edge_position('i6')
            6
        """
        from sage.combinat.cluster_algebra_quiver.surface import _get_weighted_edge
        arcs = list(self._arcs)
        boundary_edges = list(self._boundary_edges)
        if a in arcs:
            return arcs.index(a)
        elif a in boundary_edges:
            return boundary_edges.index(a)
        else:
            raise ValueError(a, " is not a user-given label of an arc/boundary edge.")

    def get_edge_user_label(self,var):
        """
        Return the label (given by user) of an arc of boundary edge
        corresponding to a variable x_i or b_i (or a product x_i*x_j).

        INPUT:
        - ``var`` -- a variable or a product of two variables from
            self.triangulation_dictionary()

        EXAMPLES::

            sage: T = ClusterTriangulation([(1,7,4),(1,5,2),(6,0,3),(2,3,0),(0,3,6),[7,4,1]], boundary_edges=[4,5,6,7])
            sage: T.get_edge_user_label(T._cluster[0])
            0
            sage: T.get_edge_user_label(T._boundary_edges_vars[2])
            6
            sage: TT = ClusterTriangulation([('j1','j1','j2'),('j3','j4','j3'),('j2','j4','j0')])
            sage: TT.get_edge_user_label(TT._cluster[1]*TT._cluster[2])
            'j2'
        """
        from sage.combinat.cluster_algebra_quiver.surface import _get_edge_user_label
        edge_user_label = _get_edge_user_label(var,self._triangulation_dictionary)
        if edge_user_label is not None:
            return edge_user_label
        else:
            raise ValueError(var, ' is not a cluster variable (or a product of cluster variables) from self.triangulation_dictionary():', self._triangulation_dictionary)

    def list_snake_graph(self, crossed_arcs, first_triangle=None,
                         final_triangle=None, first_tile_orientation=1,
                         user_labels=True):
        """
        Return the snake graph description for a list of crossed arcs.

        .. SEEALSO::

            :meth:`draw_snake_graph`, :meth:`list_band_graph`

        INPUT:

        - ``crossed_arcs`` --  labels/variables corredsponding to arcs that cross curve

            * labels from self.cluster_triangulation().triangulation() (if ``user_labels`` is ``True``)
            * variables from self.cluster_triangulation().cluster() (if ``user_labels`` is ``False``)
            * If curve crosses a self-folded triangle (ell,r,ell), then specify:
                * ``ell, (r, 'counterclockwise), ell`` (if, as gamma is about to cross r, the puncture is to the right of gamma)
                * or ``ell, (r, clockwise), ell`` (if, as gamma is about to cross r, the puncture is to the left of gamma)

        - ``first_triangle`` -- (default:``None``) the first triangle
          crossed by curve

        - ``final_triangle`` -- (default:``None``) the last triangle
          crossed by curve

        - ``first_tile_orientation`` -- (default:1) the orientation
          (either +1 or -1) for the first tile of the snake graph

        - ``fig_size`` -- (default:``None``) image size

        - ``user_labels`` -- (default:``True``) whether or not
          ``crossed_arcs`` is a list of labels

        .. NOTE::

            The direction ('RIGHT' and 'ABOVE') for the last tile of
            the snake graph description does not mean anything.

        EXAMPLES:

        Figure 8 of Musiker - Schiffler - Williams' "Bases for Cluster Algebras from Surfaces" [MSW_Bases]_ such that:
        - the arc gamma crosses arcs  arelabeled 1,2,3,4,1;
        - outer boundary edges, clockwise starting from starting point of gamma are labeled: 7,8,9, 10;
        - inner boundary edges, clockwise starting from ending point of gamma are labeled:11, 0::

            sage: T = ClusterTriangulation([(8,7,5),(5,4,1),(1,0,2),(2,11,3),(3,4,6),(6,10,9)], boundary_edges=[7,8,9,10,11,0]) # Counterclockwise
            sage: c = [item for item in T.cluster()]
            sage: ClusterSeed(T).arc_laurent_expansion([c[0],c[1],c[2],c[3],c[0]], user_labels=False) == ClusterSeed(T).arc_laurent_expansion([c[0],c[3],c[2],c[1],c[0]], user_labels=False)
            True
            sage: T.list_snake_graph([c[0],c[1],c[2],c[3],c[0]],user_labels=False)
            [[(1, (x3, x0, x4)), (2, (x1, x0, b6), 'RIGHT')],
            [(-1, (x0, x1, b6)), (-2, (b11, x1, x2), 'ABOVE')],
            [(1, (b11, x2, x1)), (2, (x5, x2, x3), 'ABOVE')],
            [(-1, (x5, x3, x2)), (-2, (x0, x3, x4), 'RIGHT')],
            [(1, (x3, x0, x4)), (2, (x1, x0, b6), 'ABOVE')]]

        Figure 10 and 11 of Musiker - Schiffler - Williams "Positivity for Cluster Algebras from Surfaces" [MSW_Positivity]_ ::

            sage: T = ClusterTriangulation([(2,3,11),(2,1,1),(4,3,12),(0,4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,13)], boundary_edges=[11,12,13,0]) # Counterclockwise
            sage: c = [item for item in T.cluster()]
            sage: r = c[1-1]
            sage: ell = c[2-1]* r
            sage: ClusterSeed(T).arc_laurent_expansion([ell,(r,'counterclockwise'), ell, c[3-1],c[4-1],c[5-1],c[6-1]], user_labels=False)== ClusterSeed(T).arc_laurent_expansion([c[6-1],c[5-1],c[4-1],c[3-1], ell,(r,'clockwise'), ell], user_labels=False) # Gamma_1
            True
            sage: T.list_snake_graph([ell,(r,'counterclockwise'), ell, c[3-1],c[4-1],c[5-1],c[6-1]], first_tile_orientation=-1, user_labels=False) # Gamma_1
            [[(-1, (x2, x0*x1, b11)),
            (-2, ((x0, 'counterclockwise'), x0*x1, (x0, 'clockwise')), 'RIGHT')],
            [(1, (x0*x1, (x0, 'counterclockwise'), (x0, 'clockwise'))),
            (2, ((x0, 'counterclockwise'), (x0, 'clockwise'), x0*x1), 'ABOVE')],
            [(-1, ((x0, 'counterclockwise'), x0*x1, (x0, 'clockwise'))),
            (-2, (x2, x0*x1, b11), 'RIGHT')],
            [(1, (x0*x1, x2, b11)), (2, (x3, x2, b12), 'RIGHT')],
            [(-1, (x2, x3, b12)), (-2, (x4, x3, b10), 'RIGHT')],
            [(1, (x3, x4, b10)), (2, (x9, x4, x5), 'ABOVE')],
            [(-1, (x9, x5, x4)), (-2, (x6, x5, x8), 'ABOVE')]]

            sage: T.list_snake_graph([c[4],c[5],c[6],c[7],c[8],c[5],c[4]],first_tile_orientation=1, user_labels=False) # Ell_p
            [[(1, (x3, x4, b10)), (2, (x9, x4, x5), 'ABOVE')],
            [(-1, (x9, x5, x4)), (-2, (x6, x5, x8), 'RIGHT')],
            [(1, (x5, x6, x8)), (2, (x7, x6, b13), 'RIGHT')],
            [(-1, (x6, x7, b13)), (-2, (x9, x7, x8), 'ABOVE')],
            [(1, (x9, x8, x7)), (2, (x6, x8, x5), 'ABOVE')],
            [(-1, (x6, x5, x8)), (-2, (x9, x5, x4), 'ABOVE')],
            [(1, (x9, x4, x5)), (2, (x3, x4, b10), 'ABOVE')]]

            sage: ClusterSeed(T).arc_laurent_expansion([c[4],c[5],c[6],c[7],c[8],c[5],c[4]], user_labels=False) == ClusterSeed(T).arc_laurent_expansion([c[4],c[5],c[8],c[7],c[6],c[5],c[4]], user_labels=False) # Ell_p
            True

        Thrice-punctured square of Figure 10 of 'Positivity for Cluster Algebras from Surfaces', [MSW_Positivity]_::

            sage: thrice_punctured_square = [('r','r','ell'),(11,'ell',3),(3,12,4),(4,5,14),(5,6,10),(6,7,9),(8,10,9),(7,13,8)]
            sage: T = ClusterTriangulation(thrice_punctured_square, boundary_edges=[11,12,13,14])
            sage: T.list_snake_graph(['ell', ('r','counterclockwise'), 'ell', 3, 4, 5, 6],first_tile_orientation=-1,user_labels=True)
            [[(-1, (3, 'ell', 11)),
            (-2, (('r', 'counterclockwise'), 'ell', ('r', 'clockwise')), 'RIGHT')],
            [(1, ('ell', ('r', 'counterclockwise'), ('r', 'clockwise'))),
            (2, (('r', 'counterclockwise'), ('r', 'clockwise'), 'ell'), 'ABOVE')],
            [(-1, (('r', 'counterclockwise'), 'ell', ('r', 'clockwise'))),
            (-2, (3, 'ell', 11), 'RIGHT')],
            [(1, ('ell', 3, 11)), (2, (4, 3, 12), 'RIGHT')],
            [(-1, (3, 4, 12)), (-2, (5, 4, 14), 'RIGHT')],
            [(1, (4, 5, 14)), (2, (10, 5, 6), 'ABOVE')],
            [(-1, (10, 6, 5)), (-2, (7, 6, 9), 'ABOVE')]]
            sage: T.list_snake_graph([5,6,7,8,9,6,5], first_tile_orientation=1, user_labels=True)
            [[(1, (4, 5, 14)), (2, (10, 5, 6), 'ABOVE')],
            [(-1, (10, 6, 5)), (-2, (7, 6, 9), 'RIGHT')],
            [(1, (6, 7, 9)), (2, (8, 7, 13), 'RIGHT')],
            [(-1, (7, 8, 13)), (-2, (10, 8, 9), 'ABOVE')],
            [(1, (10, 9, 8)), (2, (7, 9, 6), 'ABOVE')],
            [(-1, (7, 6, 9)), (-2, (10, 6, 5), 'ABOVE')],
            [(1, (10, 5, 6)), (2, (4, 5, 14), 'ABOVE')]]
        """
        from sage.combinat.cluster_algebra_quiver.surface import _snake_graph

        if user_labels:
            return _snake_graph(self._triangulation, crossed_arcs,
                                first_triangle, final_triangle, is_arc=True,
                                first_tile_orientation=first_tile_orientation,
                                boundary_edges=self._boundary_edges)

        return _snake_graph(self._weighted_triangulation, crossed_arcs,
                            first_triangle, final_triangle, is_arc=True,
                            first_tile_orientation=first_tile_orientation,
                            boundary_edges=self._boundary_edges_vars)

    def list_band_graph(self, crossed_arcs, first_triangle=None,
                        final_triangle=None, first_tile_orientation=1,
                        user_labels=True):
        """
        Return the band graph description for a list of crossed arcs.

        .. SEEALSO::

            :meth:`draw_band_graph`, :meth:`list_snake_graph`

        INPUT:

        - ``crossed_arcs`` -- labels from
          self.cluster_triangulation().triangulation() (if
          ``user_labels`` is ``True``), and variables from
          self.cluster_triangulation().cluster() if (``user_labels``
          is ``False``) corresponding to arcs that are crossed by
          curve

        If curve crosses a self-folded triangle (ell,r,ell), then specify:

        * ``ell, (r, 'counterclockwise), ell`` (if, as gamma is about
          to cross r, the puncture is to the right of gamma)

        * or ``ell, (r, clockwise), ell`` (if, as gamma is about to
          cross r, the puncture is to the left of gamma)

        - ``first_triangle`` -- (default:``None``) the first triangle
          crossed by curve

        - ``final_triangle`` -- (default:``None``) the last triangle
          crossed by curve

        - ``first_tile_orientation`` -- (default:1) the orientation
          (either +1 or -1) for the first tile of the band graph

        - ``fig_size`` -- (default:``None``) image size

        - ``user_labels`` -- (default:``True``) whether or not
          ``crossed_arcs`` is a list of labels

        .. NOTE::

            In contrast to the method :meth:`list_snake_graph` (where
            the direction for the final tile of the snake graph does
            not mean anything), the direction ('RIGHT' and 'ABOVE')
            for the final tile of the band graph means that the 'cut'
            edge label is at the right (if 'RIGHT') or the ceiling (if
            'ABOVE') of the final tile

        EXAMPLES:

        Figure 6 of Musiker and Williams' "Matrix Formulae and Skein
        Relations for Cluster Algebras from Surfaces"
        [MW_MatrixFormulae]_ where tau_4, tau_1, tau_2, tau_3 =
        `0`, `1`, `2`, `3` and b1,b2,b3,b4=`b4`, `b5`, `b6`, `b7`::

            sage: T = ClusterTriangulation([(1,2,'b4'),(1,0,'b5'),(0,3,'b6'),(2,3,'b7')], boundary_edges=['b4','b5','b6','b7'])
            sage: c = [item for item in T.cluster()]
            sage: T.list_band_graph( [c[1],c[2],c[3],c[0],c[1]], user_labels=False) # Pick cut edge to be tau_1 = 1, go clockwise
            [[(1, (b5, x1, x0)), (2, (b4, x1, x2), 'ABOVE')],
            [(-1, (b4, x2, x1)), (-2, (x3, x2, b7), 'RIGHT')],
            [(1, (x2, x3, b7)), (2, (x0, x3, b6), 'RIGHT')],
            [(-1, (x3, x0, b6)), (-2, (b5, x0, x1), 'ABOVE')]]
        """
        from sage.combinat.cluster_algebra_quiver.surface import _snake_graph

        if user_labels:
            return _snake_graph(self._triangulation, crossed_arcs,
                                first_triangle, final_triangle, is_arc=False,
                                is_loop=True,
                                first_tile_orientation=first_tile_orientation,
                                boundary_edges=self._boundary_edges)

        return _snake_graph(self._weighted_triangulation, crossed_arcs,
                            first_triangle, final_triangle, is_arc=False,
                            is_loop=True,
                            first_tile_orientation=first_tile_orientation,
                            boundary_edges=self._boundary_edges_vars)

    def draw_snake_graph(self, crossed_arcs, first_triangle=None,
                         final_triangle=None, first_tile_orientation=1,
                         fig_size=None, user_labels=True):
        """
        Display the snake graph for a list of crossed arcs.

        .. SEEALSO::

            :meth:`ClusterSeed.arc_laurent_expansion`,
            :meth:`draw_band_graph`, :meth:`list_snake_graph`

        INPUT:

        - ``crossed_arcs`` -- labels from
          self.cluster_triangulation().triangulation() (if
          ``user_labels`` is ``True``) and variables from
          self.cluster_triangulation().cluster() if (``user_labels`` is
          ``False``) corresponding to arcs that are crossed by curve If
          curve crosses a self-folded triangle (ell,r,ell), then specify
          ``ell, (r, 'counterclockwise), ell`` (if, as gamma is about to
          cross r, the puncture is to the right of gamma) or ``ell, (r,
          clockwise), ell`` (if, as gamma is about to cross r, the
          puncture is to the left of gamma)

        - ``first_triangle`` -- (default:``None``) the first triangle
          crossed by curve

        - ``final_triangle`` -- (default:``None``) the last triangle
          crossed by curve

        - ``first_tile_orientation`` -- (default:1) the orientation
          (either +1 or -1) for the first tile of the snake graph

        - ``fig_size`` -- (default:``None``) image size

        - ``user_labels`` -- (default:``True``) whether or not
          ``crossed_arcs`` is a list of labels

        EXAMPLES:

        Figure 10 of Positivity for Cluster Algebras from Surfaces, [MSW_Positivity]_::

            sage: thrice_punctured_square = [('r','r','ell'),(11,'ell',3),(3,12,4),(4,5,14),(5,6,10),(6,7,9),(8,10,9),(7,13,8)]
            sage: T = ClusterTriangulation(thrice_punctured_square, boundary_edges=[11,12,13,14])
            sage: S = ClusterSeed(T)
            sage: r = T.get_edge_var('r')
            sage: ell = T.get_edge_var('ell')
            sage: three = T.get_edge_var(3)
            sage: four = T.get_edge_var(4)
            sage: five = T.get_edge_var(5)
            sage: six = T.get_edge_var(6)
            sage: crossed_arcs = [ell, (r,'counterclockwise'), ell, three, four, five, six]
            sage: T.draw_snake_graph(crossed_arcs,first_tile_orientation=-1,user_labels=False)
            sage: T.draw_snake_graph(['ell', ('r','counterclockwise'), 'ell', 3, 4, 5, 6],first_tile_orientation=-1,user_labels=True)
            sage: T.draw_snake_graph([5,6,7,8,9,6,5], first_tile_orientation=1, user_labels=True)
        """
        from sage.combinat.cluster_algebra_quiver.surface import _snake_graph, _draw_snake_graph
        drawing = _draw_snake_graph (self.list_snake_graph(crossed_arcs, first_triangle=first_triangle, final_triangle=final_triangle, first_tile_orientation=first_tile_orientation, user_labels=user_labels), print_user_labels=user_labels)
        drawing.show ( axes=False, figsize=fig_size )

    def draw_band_graph (self, crossed_arcs, first_triangle=None, final_triangle=None, first_tile_orientation=1, fig_size=None, user_labels=True):
        """
        Display the band graph for a list of crossed arcs.

        .. SEEALSO::

            :meth:`ClusterSeed.loop_laurent_expansion`,
            :meth:`draw_snake_graph`,  :meth:`list_band_graph`

        INPUT:

        - ``crossed_arcs`` -- labels from
          self.cluster_triangulation().triangulation() (if
          ``user_labels`` is ``True``) and variables from
          self.cluster_triangulation().cluster() if (``user_labels``
          is ``False``) corresponding to arcs that are crossed by
          curve If curve crosses a self-folded triangle (ell,r,ell),
          then specify ``ell, (r, 'counterclockwise), ell`` (if, as
          gamma is about to cross r, the puncture is to the right of
          gamma) or ``ell, (r, clockwise), ell`` (if, as gamma is
          about to cross r, the puncture is to the left of gamma)

        - ``first_triangle`` -- (default:``None``) the first triangle
          crossed by curve

        - ``final_triangle`` -- (default:``None``) the last triangle
          crossed by curve

        - ``first_tile_orientation`` -- (default:1) the orientation
          (either +1 or -1) for the first tile of the band graph

        - ``fig_size`` -- (default:``None``) image size

        - ``user_labels`` -- (default:``True``) whether or not
          ``crossed_arcs`` is a list of labels

        EXAMPLES:

        Figure 6 of Musiker and Williams' "Matrix Formulae and Skein
        Relations for Cluster Algebras from Surfaces" [MW_MatrixFormulae]_
        where tau_4, tau_1, tau_2, tau_3 = ``0``, ``1``, ``2``, ``3`` and
        b1,b2,b3,b4 = ``b4``, ``b5``, ``b6``, ``b7``, Pick `tau_1` (or
        ``c`` edge) to be the edge labeled 1, and go clockwise::

            sage: T = ClusterTriangulation([(1,2,'b4'),(1,0,'b5'),(0,3,'b6'),(2,3,'b7')], boundary_edges=['b4','b5','b6','b7'])
            sage: c = [item for item in T.cluster()]
            sage: T.draw_band_graph([c[1],c[2],c[3],c[0],c[1]], user_labels=False)
            sage: T.draw_band_graph([1,2,3,0,1], user_labels=True)
        """
        from sage.combinat.cluster_algebra_quiver.surface import _snake_graph, _draw_snake_graph
        drawing = _draw_snake_graph(self.list_band_graph(crossed_arcs, first_triangle=first_triangle, final_triangle=final_triangle, first_tile_orientation=first_tile_orientation, user_labels=user_labels), print_user_labels=user_labels)
        drawing.show ( axes=False, figsize=fig_size )

    def draw_lifted_arc(self, crossed_arcs, first_triangle=None,
                        final_triangle=None, fig_size=None, verbose=False,
                        user_labels=True):
        """
        Display illustration of the lifted triangulated polygon and lifted arc.

        See [MSW_Positivity]_ section 7.

        .. SEEALSO::

            :meth:`ClusterSeed.arc_laurent_expansion`, :meth:`draw_lifted_loop`

        .. TODO::

            Make more sophisticated graphics

        INPUT:

        - ``crossed_arcs`` -- the list of variables or user-given
          labels labeling the arcs crossed by curve

        - ``first_triangle`` -- (default:``None``) a 3-tuple with
          labels of the first triangle crossed by curve

        - ``final_triangle`` -- (default:``None``) a 3-tuple with
          labels of the final triangle crossed by curve

        - ``fig_size`` -- (default:``None``) size of graphics

        - ``verbose`` -- (default ``False``) if ``True``, description
          of lifted triangles is returned

        - ``user_labels`` -- (default ``True``) if ``True``, user
          should enter ``crossed_arcs`` from the list of
          user-specified labels from
          :meth:`ClusterTriangulation.arcs`, if ``False``, user should
          enter ``crossed_arcs`` from :meth:`ClusterTriangulation.cluster`

        EXAMPLES:

        Figure 10 and 11 of Musiker - Schiffler - Williams "Positivity
        for Cluster Algebras from Surfaces" [MSW_Positivity]_
        where boundary edges are labeled 11,12,13,0; and 1=r,
        2=ell, and 3,4,5,...,10 are other arcs::

            sage: T = ClusterTriangulation([(2,3,11),(2,1,1),(4,3,12),(0,4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,13)],
            ....:   boundary_edges=[11,12,13,0])
            sage: c = [item for item in T.cluster()]
            sage: r = c[1-1]
            sage: ell = c[2-1]*c[1-1]
            sage: T.draw_lifted_arc ([ell,(r,'counterclockwise'),ell, c[3-1],c[4-1],c[5-1],c[6-1]], user_labels=False)
            sage: T.draw_lifted_arc ([2,(1,'counterclockwise'),2,3,4,5,6], user_labels=True)

        Figure 8 of Musiker - Schiffler - Williams "Bases for Cluster
        Algebras from Surfaces" [MSW_Bases]_ where: the loop gamma
        crosses the arcs labeled 1,2,3,4,1 (1 is counted twice) outer
        boundary edges, clockwise starting from starting point of
        gamma: 7,8,9, 10 inner boundary edges, clockwise starting from
        ending point of gamma:11, 0::

            sage: T = ClusterTriangulation([(8,7,5),(5,4,1),(1,0,2),(2,11,3),(3,4,6),(6,10,9)], boundary_edges=[7,8,9,10,11,0])
            sage: c = [item for item in T.cluster()]
            sage: T.draw_lifted_arc([c[1-1],c[2-1],c[3-1],c[4-1],c[1-1]], user_labels=False)
            sage: T.draw_lifted_arc([1,2,3,4,1], user_labels=True)
        """
        from sage.combinat.cluster_algebra_quiver.surface import _lifted_polygon, _draw_lifted_curve

        if user_labels:
            lifted_polygon = _lifted_polygon(self._triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=True, is_loop=False)
        else:
            lifted_polygon = _lifted_polygon(self._weighted_triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=True, is_loop=False)
        drawing = _draw_lifted_curve(lifted_polygon, is_arc=True, is_loop=False)
        drawing.show(axes=False, figsize=fig_size)
        if verbose:
            return lifted_polygon

    def draw_lifted_loop(self, crossed_arcs, first_triangle=None,
                         final_triangle=None, fig_size=None, verbose=False,
                         user_labels=True):
        """
        Display the lifted triangulated polygon and lifted loop.

        See Musiker - Schiffler - Williams "Positivity for Cluster
        Algebras from Surfaces" [MSW_Positivity]_ section 7.

        .. SEEALSO::

            :meth:`ClusterSeed.loop_laurent_expansion`, :meth:`draw_lifted_arc`

        .. TODO::

            Make more sophisticated graphics

        INPUT:

        - ``crossed_arcs`` -- the list of variables or user-given labels labeling the arcs crossed by curve

        - ``first_triangle`` -- (default:``None``) a 3-tuple with labels of the first triangle crossed by curve

        - ``final_triangle`` -- (default:``None``) a 3-tuple with labels of the final triangle crossed by curve

        - ``fig_size`` -- (default:``None``) size of graphics

        - ``verbose`` -- (default:``False``) if ``True``, description of lifted triangles is returned

        - ``user_labels`` -- (default:``True``) if ``True``, user
            should enter ``crossed_arcs`` from the list of
            user-specified labels from
            :meth:`ClusterTriangulation.arcs`, if ``False``, user
            should enter ``crossed_arcs`` from
            :meth:`ClusterTriangulation.cluster`

        EXAMPLES:

        Figure 10 and 11 of Musiker - Schiffler - Williams "Positivity for Cluster Algebras from Surfaces" [MSW_Positivity]_, where
        the boundary edges are labeled 11,12,13,0, and the arcs are labeled 1 (a radius), 2 (a noose), and 3,4,5,...,10::

            sage: T = ClusterTriangulation([(2,3,'b11'),(2,1,1),(4,3,'b12'),('b0',4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,'b13')], boundary_edges=['b11','b12','b13','b0'])
            sage: c = [item for item in T.cluster()]
            sage: r = c[1-1]
            sage: ell = c[2-1]*c[1-1]
            sage: crossed_vars = [c[3-1],ell,(r,'counterclockwise'),ell, c[3-1],c[4-1],c[5-1],c[6-1],c[7-1],c[8-1],c[9-1],c[6-1],c[5-1],c[4-1],c[3-1]]
            sage: T.draw_lifted_loop(crossed_vars, user_labels=False)
            sage: crossed_userlabels = [3,2,(1,'counterclockwise'),2, 3,4,5,6,7,8,9,6,5,4,3]
            sage: T.draw_lifted_loop(crossed_userlabels)

        Figure 8 of Musiker - Schiffler - Williams "Bases for Cluster Algebras from Surfaces" [MSW_Bases]_
        with triangulation having arcs 1, ..., 6, and gamma crosses 1,2,3,4,1 (1 is listed twice); and
        the outer boundary edges, clockwise starting from starting point of gamma are 7,8,9, 10; and
        the inner boundary edges, clockwise starting from ending point of gamma:11, 0::

            sage: T = ClusterTriangulation([(8,7,5),(5,4,1),(1,0,2),(2,11,3),(3,4,6),(6,10,9)], boundary_edges=[7,8,9,10,11,0])
            sage: c = [item for item in T.cluster()]
            sage: T.draw_lifted_loop([c[1-1],c[2-1],c[3-1],c[4-1],c[1-1]], user_labels=False)
            sage: T.draw_lifted_loop([1,2,3,4,1], user_labels=True)
        """
        from sage.combinat.cluster_algebra_quiver.surface import _lifted_polygon, _draw_lifted_curve
        if user_labels:
            lifted_polygon = _lifted_polygon(self._triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=False, is_loop=True)
        else:
            lifted_polygon = _lifted_polygon(self._weighted_triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=False, is_loop=True)
        drawing = _draw_lifted_curve(lifted_polygon, is_arc=False, is_loop=True)
        drawing.show ( axes=False, figsize=fig_size)
        if verbose:
            return lifted_polygon
