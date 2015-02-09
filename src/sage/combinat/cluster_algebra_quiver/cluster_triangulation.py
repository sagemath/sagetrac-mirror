r"""
ClusterTriangulation

An *cluster triangulation* (see :arXiv:math/`0608367`) is an ideal triangulation, i.e. a maximal collection of distinct non-crossing arcs.
An ideal triangulation of a surface with marked points (S,M) is associated to a seed of a cluster algebra arising from (S,M).

.. SEEALSO:: Cluster triangulations are closely related to :meth:`~sage.combinat.cluster_algebra_quiver.cluster_seed.ClusterSeed`
and :meth:`~sage.combinat.cluster_algebra_quiver.quiver.ClusterQuiver`

"""

from sage.structure.sage_object import SageObject

######################################################################################################
############# begins: CREATING CLUSTER ALGEBRA FROM INITIAL TRIANGULATION INPUT ###########
######################################################################################################

class ClusterTriangulation(SageObject):
    r"""
    An initial *ideal triangulation* associated to a surface

    INPUT:

    - ``data`` -- can be any of the following::
        * List of triangles - must be the list of 3-tuples (i.e. edge labels of a triangle) from an ideal triangulation (see Examples)

    .. TODO::

    - ``data`` -- can also be any of the following::
        * Matrix - a skew-symmetrizable matrix arising from a tagged triangulation (todo: not yet implemented)
        * DiGraph - must be the input data for a quiver from a tagged triangulation (todo: not yet implemented)
        * List of edges - must be the edge list of a digraph for a quiver from a tagged triangulation (todo: not yet implemented)
        * Objects that Theodosios Douvropoulos is currently working with (not related to cluster algebras from surfaces) (todo: not yet implemented)

    EXAMPLES::

        from a List of ideal triangles (forming an ideal triangulation of a surface)::

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

            sage: annulus22 = [('bd1','tau1','tau2'),('tau2','tau3','bd4'),('tau1','tau4','bd2'),('tau3','bd3','tau4')]
            sage: T = ClusterTriangulation(annulus22, boundary_edges=['bd3','bd2','bd1','bd4'])
            sage: T
            An ideal triangulation associated with cluster algebra of rank 4 with 4 boundary edges
            sage: ClusterSeed(T).mutation_type() # Figure 6 of arXiv:1108.3382
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

            sage: thrice_punctured_square = [(2,2,1), (1,3,11), (3,12,4), (4,5,14), (5,6,10), (6,7,9), (8,10,9), (7,13,8)]
            sage: T = ClusterTriangulation(thrice_punctured_square, boundary_edges=[14,12,13,11]) # Figure 15 of arXiv:0906.0748
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

            sage: T = ClusterTriangulation ( [[4, 5, 1], [4, 3, 2], [3, 7, 2], [2, 1, 6], [1, 4, 5]], boundary_edges=[1] )
            sage: T
            An ideal triangulation associated with cluster algebra of rank 6 with 1 boundary edges

            sage: T = ClusterTriangulation ( [(4, 5, 1), (4, 3, 2), (3, 7, 2), (2, 1, 6), (1, 4, 5)] )
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

    """
    def __init__(self, data, boundary_edges=None):
        r"""
        TESTS::

            sage: CT = ClusterTriangulation([('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')], boundary_edges=['c','d','e','f'])
            sage: TestSuite(CT).run() # what CT should I put here? (todo)

        """
        from sage.combinat.cluster_algebra_quiver.surface import remove_duplicate_triangles, _triangulation_to_arrows, _surface_edge_list_to_matrix, _get_user_arc_labels, _get_triangulation_dictionary, _get_weighted_triangulation
        from sage.rings.all import QQ
        from sage.rings.all import FractionField, PolynomialRing

        self._boundary_edges = boundary_edges if boundary_edges else []

        # if data is a list of triangles (forming a triangulation)
        if isinstance(data,list) and \
        all(type(triangle) in [list,tuple] and len(triangle)==3 for triangle in data):
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
            self._weighted_triangulation = _get_weighted_triangulation (self._triangles, self._triangulation_dictionary)

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
            sage: CT.__eq__( CT2 )
            True
        """
        return isinstance(other, ClusterTriangulation) and self._weighted_triangulation == other._weighted_triangulation

    def _repr_(self):
        """
        Returns the description of ``self``.

        EXAMPLES::

            sage: T = ClusterTriangulation ( [(4, 5, 1), (4, 3, 2), (3, 7, 2), (2, 1, 6), (1, 4, 5)] )
            sage: T._repr_()
            'An ideal triangulation associated with cluster algebra of rank 7'
        """
        name = self._description
        return name

    def triangles(self):
        """
        Returns the underlying triangulation of ``self``.

        EXAMPLES::
            sage: Triangles = [(1,4,7),(1,2,5),(2,0,3),(0,6,3),(6,3,0)]
            sage: T = ClusterTriangulation(Triangles)
            sage: T.triangles()
            [(1, 4, 7), (1, 2, 5), (2, 0, 3), (0, 6, 3)]
        """
        return self._triangles

    def b_matrix(self):
        """
        Returns the B-matrix of the coefficient-free cluster.  The conventions
        for B-matrices are the opposite of arXiv:math/0608367.

        EXAMPLES::

            #sage: twice_punctured_monogon = [[1,4,2],[3,4,3],[2,0,1]]
            sage: twice_punctured_monogon = [[4,5,5],[1,2,3],[2,4,3]]
            sage: T = ClusterTriangulation(twice_punctured_monogon) # twice-punctured monogon with 3 (non-ordinary) ideal triangles (affine D)
            sage: B = T.b_matrix()
            sage: B
            [ 0  1 -1  0  0]
            [-1  0  0  1  1]
            [ 1  0  0 -1 -1]
            [ 0 -1  1  0  0]
            [ 0 -1  1  0  0]
            sage: twice_punctured_monogon_mu2 = [(4,5,5),(2,3,3),(1,4,2)]
            sage: Tmu2 = ClusterTriangulation(twice_punctured_monogon_mu2) # 2 self-folded triangles and 1 triangle with one vertex (affine D)
            sage: Bmu2 = Tmu2.b_matrix() #Figure 9 (right) of arXiv:math/0608367
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

            sage: twice_punctured_monogon = [(1,1,2), (4,4,3), ('boundary',2,3)] # Figure 10 (bottom) of arXiv:math/0608367
            sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=['boundary'])
            sage: B = T.b_matrix(); B
            [ 0  0  1  1]
            [ 0  0  1  1]
            [-1 -1  0  0]
            [-1 -1  0  0]
            sage: twice_punctured_monogon_mu3 = [(1,1,2), (2,4,3), (3,4,'boundary')] # Figure 10 (top) of arXiv:math/0608367
            sage: Tmu3 = ClusterTriangulation(twice_punctured_monogon_mu3, boundary_edges=['boundary'])
            sage: Bmu3 = Tmu3.b_matrix(); Bmu3
            [ 0  0 -1  1]
            [ 0  0 -1  1]
            [ 1  1  0  0]
            [-1 -1  0  0]
            sage: B.mutate(2)
            sage: Bmu3 == B
            True

        """
        return self._M

    def arcs(self):
        """
        Return the labels of diagonals (not boundary edges) of ``self`` given by user.

        EXAMPLES::

            sage: annulus22 = [('bd1','tau1','tau2'),('tau2','tau3','bd4'),('tau1','tau4','bd2'),('tau3','bd3','tau4')]
            sage: T = ClusterTriangulation(annulus22, boundary_edges=['bd3','bd2','bd1','bd4'])
            sage: T.arcs()
            ['bd1', 'bd2', 'bd3', 'bd4', 'tau1', 'tau2', 'tau3', 'tau4']
        """
        return self._arcs

    def boundary_edges(self):
        """
        Return the labels of boundary edges (not diagonals) of ``self`` given by user.

        EXAMPLES::
            sage: annulus22 = [('bd1','tau1','tau2'),('tau2','tau3','bd4'),('tau1','tau4','bd2'),('tau3','bd3','tau4')]
            sage: T = ClusterTriangulation(annulus22, boundary_edges=['bd3','bd2','bd1','bd4'])
            sage: T.boundary_edges()
            ['bd1', 'bd2', 'bd3', 'bd4']

        """
        return self._boundary_edges

    def cluster(self):
        """
        Return the cluster seed list (i.e. [x_0, x_1, ..., x_{n-1}]) corresponding to ``self``.

        EXAMPLES::

            sage: once_punctured_square = [('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')]
            sage: T = ClusterTriangulation(once_punctured_square, boundary_edges=['c','f','e','d'])
            sage: T.cluster()
            [x0, x1, x2, x3]
        """
        return self._cluster

    def boundary_edges_vars(self):
        """
        Return the cluster seed list (i.e. [x_0, x_1, ..., x_{n-1}]) corresponding to ``self``.

        EXAMPLES::

            sage: once_punctured_square = [('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')]
            sage: T = ClusterTriangulation(once_punctured_square, boundary_edges=['c','f','e','d'])
            sage: T.boundary_edges_vars()
            [b4, b5, b6, b7]

        """
        return self._boundary_edges_vars

    def triangulation_dictionary(self):
        """
        Return the correspondence between user-given labels (integers) and variables of ``self``.

        EXAMPLES::

            sage: twice_punctured_monogon = [(4,5,5),(2,3,3),(1,4,2)] #Figure 10 (bottom) of arXiv:math/0608367
            sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=[1]) # 2 self-folded triangles and 1 triangle with one vertex (affine D)
            sage: T.triangulation_dictionary()
            [(2, x0*x1), (3, x1), (4, x2*x3), (5, x3), (1, b4)]
        """
        return self._triangulation_dictionary

    def weighted_triangulation(self):
        """
        Return the list of triangles of ``self`` where user-given labels are replaced by the variables corresponding to them.

        EXAMPLES::

            sage: Triangles = [(1,4,7),(1,2,5),(2,0,3),(0,6,3)]
            sage: CT = ClusterTriangulation(Triangles)
            sage: CT.weighted_triangulation()
            [(x1, x4, x7), (x1, x2, x5), (x2, x0, x3), (x0, x6, x3)]

            sage: twice_punctured_monogon = [(4,5,5),(2,3,3),(1,4,2)] #Figure 10 (bottom) of arXiv:math/0608367
            sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=[1]) # 2 self-folded triangles and 1 triangle with one vertex (affine D)
            sage: T.weighted_triangulation()
            [((x3, 'counterclockwise'), (x3, 'clockwise'), x2*x3),
            ((x1, 'counterclockwise'), (x1, 'clockwise'), x0*x1),
            (b4, x2*x3, x0*x1)]
        """
        return self._weighted_triangulation

    def get_edge_var(self,a):
        """
        Return the variable corresponding to the label (given by user) of an arc of boundary edge

        EXAMPLES::

            sage: T = ClusterTriangulation([(1,7,4),(1,5,2),(6,0,3),(2,3,0),(0,3,6),[7,4,1]], boundary_edges=[4,5,6,7])
            sage: T.get_edge_var(0)
            x0
            sage: T.get_edge_var(6)
            b6
        """
        from sage.combinat.cluster_algebra_quiver.surface import _get_weighted_edge
        return _get_weighted_edge(a,self._triangulation_dictionary)