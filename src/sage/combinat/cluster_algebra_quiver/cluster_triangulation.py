r"""
ClusterTriangulation

An *cluster triangulation* (see :arxiv:`math/0608367`) is a subset of an ideal triangulation, i.e. a maximal collection of distinct non-crossing arcs.
An ideal triangulation of a surface with marked points (S,M) is associated to a seed of a cluster algebra arising from (S,M).

.. SEEALSO::

    Cluster triangulations closely interact with
    :class:`~sage.combinat.cluster_algebra_quiver.cluster_seed.ClusterSeed`,
    :class:`~sage.combinat.cluster_algebra_quiver.cluster_seed.ClusterQuiver`,
"""

from sage.structure.sage_object import SageObject

class ClusterTriangulation(SageObject):
    """
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
            sage: ClusterSeed(T).mutation_type() # Figure 6 of :arxiv:`1108.3382`
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

             Figure 15 of :arxiv:`0906.0748`::

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
        TESTS::

            sage: CT = ClusterTriangulation([('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')], boundary_edges=['c','d','e','f'])
            sage: TestSuite(CT).run() # what CT should I put here? (todo)
        """
        from sage.combinat.cluster_algebra_quiver.surface import remove_duplicate_triangles, _triangulation_to_arrows, _surface_edge_list_to_matrix, _get_user_arc_labels, _get_triangulation_dictionary, _get_user_label_triangulation, _get_weighted_triangulation
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
            self._triangulation = _get_user_label_triangulation(self._triangles)
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
            sage: CT.__eq__(CT2)
            True
        """
        return isinstance(other, ClusterTriangulation) and self._weighted_triangulation == other._weighted_triangulation

    def _repr_(self):
        """
        Returns the description of ``self``.

        EXAMPLES::

            sage: T = ClusterTriangulation([(4, 5, 1), (4, 3, 2), (3, 7, 2), (2, 1, 6), (1, 4, 5)])
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
        for B-matrices are the opposite of :arxiv:`math/0608367`.

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
            sage: Bmu2 = Tmu2.b_matrix() #Figure 9 (right) of :arxiv:`math/0608367`
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

            sage: twice_punctured_monogon = [(1,1,2), (4,4,3), ('boundary',2,3)] # Figure 10 (bottom) of :arxiv:`math/0608367`
            sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=['boundary'])
            sage: B = T.b_matrix()
            sage: B
            [ 0  0  1  1]
            [ 0  0  1  1]
            [-1 -1  0  0]
            [-1 -1  0  0]
            sage: twice_punctured_monogon_mu3 = [(1,1,2), (2,4,3), (3,4,'boundary')] # Figure 10 (top) of :arxiv:`math/0608367`
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

             2 self-folded triangles and 1 triangle with one vertex (affine D), Figure 10 (bottom) of :arxiv:`math/0608367`::

                sage: twice_punctured_monogon = [(4,5,5),(2,3,3),(1,4,2)]
                sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=[1])
                sage: T.triangulation_dictionary()
                [(2, x0*x1), (3, x1), (4, x2*x3), (5, x3), (1, b4)]
        """
        return self._triangulation_dictionary

    def triangulation(self):
        """
        Return the list of triangles of ``self`` where a self-folded triangle (r,r,ell) is replaced by (r, 'counterclockwise'), (r, 'clockwise'), ell)

        EXAMPLES::

            sage: Triangles = [(1,4,7),(1,2,5),(2,0,3),(0,6,3)]
            sage: CT = ClusterTriangulation(Triangles)
            sage: CT.triangulation()
            [(1, 4, 7), (1, 2, 5), (2, 0, 3), (0, 6, 3)]

            Figure 10 (bottom) of :arxiv:`math/0608367`::

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
        Return the list of triangles of ``self`` where user-given labels are replaced by the variables corresponding to them.

        EXAMPLES::

            sage: Triangles = [(1,4,7),(1,2,5),(2,0,3),(0,6,3)]
            sage: CT = ClusterTriangulation(Triangles)
            sage: CT.weighted_triangulation()
            [(x1, x4, x7), (x1, x2, x5), (x2, x0, x3), (x0, x6, x3)]

            Two self-folded triangles and 1 triangle with one vertex (affine D). See Figure 10 (bottom) of :arxiv:`math/0608367`::

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

    def get_edge_user_label(self,var):
        """
        Return the label (given by user) of an arc of boundary edge corresponding to a variable x_i or b_i

        EXAMPLES::

            sage: T = ClusterTriangulation([(1,7,4),(1,5,2),(6,0,3),(2,3,0),(0,3,6),[7,4,1]], boundary_edges=[4,5,6,7])
            sage: T.get_edge_user_label(T._cluster[0])
            0
            sage: T.get_edge_user_label(T._boundary_edges_vars[2])
            6
        """
        from sage.combinat.cluster_algebra_quiver.surface import _get_edge_user_label
        return _get_edge_user_label(var,self._triangulation_dictionary)


    def snake_graph(self, crossed_arcs, first_triangle=None, final_triangle=None, first_tile_orientation=1, user_labels=False):
        """
        INPUT:
        crossed_arcs = x0, x1, ... etc
        (optional) first_triangle = (a,b,c) -> the first triangle crossed by arc
        (optional) final_triangle = (d,e,f) -> the final triangle crossed by arc

        Algorithm will try to find the first triangle and final triangle crossed by arc, and
        will warn user to enter these triangles if the algorithm fails

        1 labels the bottom triangle of a positively-oriented tile,
        2 labels the top triangle of a positively-oriented tile,
        -1 labels the bottom triangle of a positively-oriented tile,
        -2 labels the top triangle of a negatively-oriented tile.

        The direction (RIGHT or ABOVE) that is attached to the top triangle (labeled -2 or 2)
        indicates the location of the tile after the current tile.

        If this is a snake graph (not a band graph),
        then the direction for the last tile does not mean anything

        EXAMPLES::

            Figure 8 of Musiker - Schiffler - Williams' "Bases for Cluster Algebras from Surfaces" :arxiv:`1110.4364` such that:
            the arc gamma crosses arcs  arelabeled 1,2,3,4,1;
            outer boundary edges, clockwise starting from starting point of gamma are labeled: 7,8,9, 10;
            inner boundary edges, clockwise starting from ending point of gamma are labeled:11, 0::

                sage: T = ClusterTriangulation([(8,7,5),(5,4,1),(1,0,2),(2,11,3),(3,4,6),(6,10,9)], boundary_edges=[7,8,9,10,11,0]) # Counterclockwise
                sage: c = [item for item in T.cluster()]
                sage: ClusterSeed(T).arc ([c[0],c[1],c[2],c[3],c[0]], user_labels=False) == ClusterSeed(T).arc([c[0],c[3],c[2],c[1],c[0]], user_labels=False)
                True
                sage: T.snake_graph([c[0],c[1],c[2],c[3],c[0]],user_labels=False)
                [[(1, (x3, x0, x4)), (2, (x1, x0, b6), 'RIGHT')],
                [(-1, (x0, x1, b6)), (-2, (b11, x1, x2), 'ABOVE')],
                [(1, (b11, x2, x1)), (2, (x5, x2, x3), 'ABOVE')],
                [(-1, (x5, x3, x2)), (-2, (x0, x3, x4), 'RIGHT')],
                [(1, (x3, x0, x4)), (2, (x1, x0, b6), 'ABOVE')]]

            Figure 10 and 11 of Musiker - Schiffler - Williams "Positivity for Cluster Algebras from Surfaces" :arxiv:`0906.0748`::

                sage: T = ClusterTriangulation([(2,3,11),(2,1,1),(4,3,12),(0,4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,13)], boundary_edges=[11,12,13,0]) # Counterclockwise
                sage: c = [item for item in T.cluster()]
                sage: r = c[1-1]
                sage: ell = c[2-1]* r
                sage: ClusterSeed(T).arc([ell,(r,'counterclockwise'), ell, c[3-1],c[4-1],c[5-1],c[6-1]], user_labels=False)== ClusterSeed(T).arc([c[6-1],c[5-1],c[4-1],c[3-1], ell,(r,'clockwise'), ell], user_labels=False) # Gamma_1
                True
                sage: T.snake_graph([ell,(r,'counterclockwise'), ell, c[3-1],c[4-1],c[5-1],c[6-1]], first_tile_orientation=-1, user_labels=False) # Gamma_1
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

                sage: T.snake_graph([c[4],c[5],c[6],c[7],c[8],c[5],c[4]],first_tile_orientation=1, user_labels=False) # Ell_p
                [[(1, (x3, x4, b10)), (2, (x9, x4, x5), 'ABOVE')],
                [(-1, (x9, x5, x4)), (-2, (x6, x5, x8), 'RIGHT')],
                [(1, (x5, x6, x8)), (2, (x7, x6, b13), 'RIGHT')],
                [(-1, (x6, x7, b13)), (-2, (x9, x7, x8), 'ABOVE')],
                [(1, (x9, x8, x7)), (2, (x6, x8, x5), 'ABOVE')],
                [(-1, (x6, x5, x8)), (-2, (x9, x5, x4), 'ABOVE')],
                [(1, (x9, x4, x5)), (2, (x3, x4, b10), 'ABOVE')]]

                sage: ClusterSeed(T).arc([c[4],c[5],c[6],c[7],c[8],c[5],c[4]], user_labels=False) == ClusterSeed(T).arc([c[4],c[5],c[8],c[7],c[6],c[5],c[4]], user_labels=False) # Ell_p
                True

            Thrice-punctured square of Figure 10 of 'Positivity for Cluster Algebras from Surfaces', :arxiv:`0906.0748`::

                sage: thrice_punctured_square = [('r','r','ell'),(11,'ell',3),(3,12,4),(4,5,14),(5,6,10),(6,7,9),(8,10,9),(7,13,8)]
                sage: T = ClusterTriangulation(thrice_punctured_square, boundary_edges=[11,12,13,14])
                sage: T.snake_graph(['ell', ('r','counterclockwise'), 'ell', 3, 4, 5, 6],first_tile_orientation=-1,user_labels=True)
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
                sage: T.snake_graph([5,6,7,8,9,6,5], first_tile_orientation=1, user_labels=True)
                [[(1, (4, 5, 14)), (2, (10, 5, 6), 'ABOVE')],
                [(-1, (10, 6, 5)), (-2, (7, 6, 9), 'RIGHT')],
                [(1, (6, 7, 9)), (2, (8, 7, 13), 'RIGHT')],
                [(-1, (7, 8, 13)), (-2, (10, 8, 9), 'ABOVE')],
                [(1, (10, 9, 8)), (2, (7, 9, 6), 'ABOVE')],
                [(-1, (7, 6, 9)), (-2, (10, 6, 5), 'ABOVE')],
                [(1, (10, 5, 6)), (2, (4, 5, 14), 'ABOVE')]]

        REFERENCES:

        .. Musiker-Schiffler-Williams, Positivity for Cluster Algebras from Surfaces, :arxiv:`0906.0748`
        and Bases for Cluster Algebras from Surfaces :arxiv:`1110.4364`
        """
        from sage.combinat.cluster_algebra_quiver.surface import _snake_graph

        if user_labels:
            return _snake_graph (self._triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=True, first_tile_orientation=first_tile_orientation, boundary_edges=self._boundary_edges)
        else:
            return _snake_graph (self._weighted_triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=True, first_tile_orientation=first_tile_orientation, boundary_edges=self._boundary_edges_vars)

    def band_graph(self, crossed_arcs, first_triangle=None, final_triangle=None, first_tile_orientation=1, user_labels=False):
       """
       Return the band graph description for a list of crossed arcs

       EXAMPLES::

            Figure 6 of Musiker and Williams' "Matrix Formulae and Skein Relations for Cluster Algebras from Surfaces" :arXiv:`1108.3382`
            where tau_4, tau_1, tau_2, tau_3 = ``0``,``1``,``2``,``3`` and b1,b2,b3,b4=``b4``,``b5``,``b6``,``b7``::

                sage: T = ClusterTriangulation([(1,2,'b4'),(1,0,'b5'),(0,3,'b6'),(2,3,'b7')], boundary_edges=['b4','b5','b6','b7'])
                sage: c = [item for item in T.cluster()]
                sage: T.band_graph( [c[1],c[2],c[3],c[0],c[1]], user_labels=False) # Pick cut edge to be tau_1 = 1, go clockwise
                [[(1, (b5, x1, x0)), (2, (b4, x1, x2), 'ABOVE')],
                [(-1, (b4, x2, x1)), (-2, (x3, x2, b7), 'RIGHT')],
                [(1, (x2, x3, b7)), (2, (x0, x3, b6), 'RIGHT')],
                [(-1, (x3, x0, b6)), (-2, (b5, x0, x1), 'ABOVE')]]
       """
       from sage.combinat.cluster_algebra_quiver.surface import _snake_graph

       if user_labels:
           return _snake_graph (self._triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=False, is_loop=True, first_tile_orientation=first_tile_orientation, boundary_edges=self._boundary_edges)
       else:
           return _snake_graph (self._weighted_triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=False, is_loop=True, first_tile_orientation=first_tile_orientation, boundary_edges=self._boundary_edges_vars)

    def draw_snake_graph(self, crossed_arcs, first_triangle=None, final_triangle=None, first_tile_orientation=1, fig_size=None, user_labels=True):
        """
        Display the snake graph for a list of crossed arcs

        EXAMPLES::

            Figure 10 of Positivity for Cluster Algebras from Surfaces, :arXiv:`0906.0748`::

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
        drawing = _draw_snake_graph (self.snake_graph(crossed_arcs, first_triangle=first_triangle, final_triangle=final_triangle, first_tile_orientation=first_tile_orientation, user_labels=user_labels), user_labels=user_labels)
        drawing.show ( axes=False, figsize=fig_size )

    def draw_band_graph (self, crossed_arcs, first_triangle=None, final_triangle=None, first_tile_orientation=1, fig_size=None, user_labels=True):
        """
        Display the band graph for a list of crossed arcs

        EXAMPLES::

            Figure 6 of Musiker and Williams' "Matrix Formulae and Skein Relations for Cluster Algebras from Surfaces" :arXiv:`1108.3382`
            where tau_4, tau_1, tau_2, tau_3 = ``0``,``1``,``2``,``3`` and b1,b2,b3,b4=``b4``,``b5``,``b6``,``b7``,
            Pick `tau_1` (or ``c`` edge) to be the edge labeled 1, and go clockwise::

                sage: T = ClusterTriangulation([(1,2,'b4'),(1,0,'b5'),(0,3,'b6'),(2,3,'b7')], boundary_edges=['b4','b5','b6','b7'])
                sage: c = [item for item in T.cluster()]
                sage: T.draw_band_graph([c[1],c[2],c[3],c[0],c[1]], user_labels=False)
                sage: T.draw_band_graph([1,2,3,0,1], user_labels=True)
        """
        from sage.combinat.cluster_algebra_quiver.surface import _snake_graph, _draw_snake_graph
        drawing = _draw_snake_graph ( self.band_graph (crossed_arcs, first_triangle=first_triangle, final_triangle=final_triangle, first_tile_orientation=first_tile_orientation, user_labels=user_labels), user_labels=user_labels)
        drawing.show ( axes=False, figsize=fig_size )

    def draw_lifted_polygon_arc (self, crossed_arcs, first_triangle=None, final_triangle=None, fig_size=None, return_value=False, user_labels=True):
        """
        Display illustration of the lifted triangulatied polygon and lifted curve (see :arxiv:`0906.0748` section 7)

        .. SEEALSO:: :meth:`ClusterSeed.arc`, :meth:``draw_lifted_polygon_loop``

        EXAMPLES::

            Figure 10 and 11 of Musiker - Schiffler - Williams "Positivity for Cluster Algebras from Surfaces" :arxiv:`0906.0748`
            where boundary edges are labeled 11,12,13,0; and  1=r, 2=ell, and 3,4,5,...,10 are other arcs::

                sage: T = ClusterTriangulation([(2,3,11),(2,1,1),(4,3,12),(0,4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,13)], boundary_edges=[11,12,13,0])
                sage: c = [item for item in T.cluster()]
                sage: r = c[1-1]
                sage: ell = c[2-1]*c[1-1]
                sage: T.draw_lifted_polygon_arc ([ell,(r,'counterclockwise'),ell, c[3-1],c[4-1],c[5-1],c[6-1]], user_labels=False)
                sage: T.draw_lifted_polygon_arc ([2,(1,'counterclockwise'),2,3,4,5,6], user_labels=True)

            Figure 8 of Musiker - Schiffler - Williams "Bases for Cluster Algebras from Surfaces" :arxiv:`1110.4364` where:
            the loop gamma crosses the arcs labeled 1,2,3,4,1 (1 is counted twice)
            outer boundary edges, clockwise starting from starting point of gamma: 7,8,9, 10
            inner boundary edges, clockwise starting from ending point of gamma:11, 0::

                sage: T = ClusterTriangulation([(8,7,5),(5,4,1),(1,0,2),(2,11,3),(3,4,6),(6,10,9)], boundary_edges=[7,8,9,10,11,0])
                sage: c = [item for item in T.cluster()]
                sage: T.draw_lifted_polygon_arc([c[1-1],c[2-1],c[3-1],c[4-1],c[1-1]], user_labels=False)
                sage: T.draw_lifted_polygon_arc([1,2,3,4,1], user_labels=True)
        """
        from sage.combinat.cluster_algebra_quiver.surface import _lifted_polygon, _draw_lifted_polygon

        #if self._weighted_boundary_edges is not None:
        #    boundary_edges_crossed = set(self._weighted_boundary_edges).intersection(crossed_arcs)
        #    if boundary_edges_crossed != set():
        #        raise ValueError (list[boundary_edges_crossed], ' are both boundary edges and crossed arcs.')

        if user_labels:
            lifted_polygon = _lifted_polygon(self._triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=True, is_loop=False)
        else:
            lifted_polygon = _lifted_polygon(self._weighted_triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=True, is_loop=False)
        drawing = _draw_lifted_polygon(lifted_polygon, is_arc=True, is_loop=False)
        drawing.show ( axes=False, figsize=fig_size)
        if return_value:
            return lifted_polygon

    def draw_lifted_polygon_loop(self, crossed_arcs, first_triangle=None, final_triangle=None, fig_size=None, return_value=False, user_labels=True):
        """
        Display the lifted triangulated polygon and lifted curve (see :arxiv:`0906.0748` section 7)

        .. SEEALSO:: :meth:`ClusterSeed.loop`, :meth:``draw_lifted_polygon_arc``

        EXAMPLES::

            Figure 10 and 11 of Musiker - Schiffler - Williams "Positivity for Cluster Algebras from Surfaces" :arxiv:`0906.0748`, where
            the boundary edges are labeled 11,12,13,0, and the arcs are labeled 1 (a radius), 2 (a noose), and 3,4,5,...,10::

                sage: T = ClusterTriangulation([(2,3,11),(2,1,1),(4,3,12),(0,4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,13)], boundary_edges=[11,12,13,0])
                sage: c = [item for item in T.cluster()]
                sage: r = c[1-1]
                sage: ell = c[2-1]*c[1-1]
                sage: T.draw_lifted_polygon_loop([c[3-1],ell,(r,'counterclockwise'),ell, c[3-1],c[4-1],c[5-1],c[6-1],c[7-1],c[8-1],c[9-1],c[6-1],c[5-1],c[4-1],c[3-1]], user_labels=False)

            Figure 8 of Musiker - Schiffler - Williams "Bases for Cluster Algebras from Surfaces" :arxiv:`1110.4364`
            with triangulation having arcs 1, ..., 6, and gamma crosses 1,2,3,4,1 (1 is listed twice); and
            the outer boundary edges, clockwise starting from starting point of gamma are 7,8,9, 10; and
            the inner boundary edges, clockwise starting from ending point of gamma:11, 0::

                sage: T = ClusterTriangulation([(8,7,5),(5,4,1),(1,0,2),(2,11,3),(3,4,6),(6,10,9)], boundary_edges=[7,8,9,10,11,0])
                sage: c = [item for item in T.cluster()]
                sage: T.draw_lifted_polygon_loop([c[1-1],c[2-1],c[3-1],c[4-1],c[1-1]], user_labels=False)
                sage: T.draw_lifted_polygon_loop([1,2,3,4,1], user_labels=True)
        """
        from sage.combinat.cluster_algebra_quiver.surface import _lifted_polygon, _draw_lifted_polygon
        if user_labels:
            lifted_polygon = _lifted_polygon(self._triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=False, is_loop=True)
        else:
            lifted_polygon = _lifted_polygon(self._weighted_triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=False, is_loop=True)
        drawing = _draw_lifted_polygon(lifted_polygon, is_arc=False, is_loop=True)
        drawing.show ( axes=False, figsize=fig_size)
        if return_value:
            return lifted_polygon
