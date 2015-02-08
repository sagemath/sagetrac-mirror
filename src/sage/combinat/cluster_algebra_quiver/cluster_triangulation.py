r"""
ClusterTriangulation

An Ideal triangulation associated to a cluster algebra (arXiv:math/0608367) or other triangulations (not implemented yet)

Cluster triangulations are closely related to :meth:`~sage.combinat.cluster_algebra_quiver.cluster_seed.ClusterSeed`
and :meth:`~sage.combinat.cluster_algebra_quiver.quiver.ClusterQuiver`

"""

from sage.structure.sage_object import SageObject
#from copy import copy
#from sage.structure.unique_representation import UniqueRepresentation

######################################################################################################
############# begins: CREATING CLUSTER ALGEBRA FROM INITIAL TRIANGULATION INPUT ###########
######################################################################################################

class ClusterTriangulation(SageObject):
    r"""
    An initial *ideal triangulation* associated to a surface

    INPUT:

    - ``data`` -- can be any of the following::

        * Matrix - a skew-symmetrizable matrix arising from a tagged triangulation (todo: not yet implemented)
        * DiGraph - must be the input data for a quiver from a tagged triangulation (todo: not yet implemented)
        * List of triangles - must be the list of triangles (counterclockwise) from an ideal triangulation (see Examples)
        * List of edges - must be the edge list of a digraph for a quiver from a tagged triangulation (todo: not yet implemented)
        * Objects that Theodosios Douvropoulos is currently working with (not related to cluster algebras from surfaces) (todo: not yet implemented)

    EXAMPLES::

        from a List of ideal triangles (forming an ideal triangulation of a surface)::

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

    """
    def __init__(self, data, boundary_edges=None):
        r"""
        TESTS::

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
            data = remove_duplicate_triangles (data)

            self._n = len(self._arcs) - len(self._boundary_edges) #if boundary_edges else len(self._arcs)
            self._triangles = data

            all_arrows = _triangulation_to_arrows(self._triangles)
            M = _surface_edge_list_to_matrix (all_arrows, self._arcs, self._boundary_edges, self._n)
            self._M = M[:self._n,:self._n] # In this implementation, we ignore the boundary edges (TODO)
            self._is_cluster_algebra = True
            self._description = 'An ideal triangulation associated with cluster algebra of rank %d'  %self._n
            if boundary_edges:
                self._description += ' with %d boundary edges' %len(self._boundary_edges)
            self._R = FractionField(PolynomialRing(QQ,['x%s'%i for i in range(0,self._n)]+['b%s'%i for i in range(self._n,len(self._boundary_edges)+self._n)]))
            self._cluster = list(self._R.gens()[0:self._n])
            self._boundary_edges_vars = self._R.gens()[self._n:] if boundary_edges else []
            self._triangulation_dictionary = _get_triangulation_dictionary (self._triangles, self._cluster, self._boundary_edges, self._boundary_edges_vars)
            self._weighted_triangulation = _get_weighted_triangulation (self._triangles, self._triangulation_dictionary)

        else:
            raise ValueError('Input must be a list of three-tuples')

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
        Returns the B-matrix of the standard quiver of self.  The conventions
        for B-matrices agree with Fomin-Zelevinsky (up to a reordering of
        the simple roots).

        EXAMPLES::

            sage: T = ClusterTriangulation([[1,4,2],[3,4,3],[2,0,1]]) # twice-punctured monogon with 3 (non-ordinary) ideal triangles (affine D)
            sage: B = T.b_matrix()
            sage: B
            [ 0  1 -1  0  0]
            [-1  0  0  1  1]
            [ 1  0  0 -1 -1]
            [ 0 -1  1  0  0]
            [ 0 -1  1  0  0]

            sage: Tmu2 = ClusterTriangulation([(1,1,2),(3,4,3),(2,4,0)]) # 2 self-folded triangles and 1 triangle with one vertex (affine D)
            sage: Bmu2 = Tmu2.b_matrix()
            sage: Bmu2
            [ 0  1  1 -1 -1]
            [-1  0  0  1  1]
            [-1  0  0  1  1]
            [ 1 -1 -1  0  0]
            [ 1 -1 -1  0  0]
            sage: B.mutate(2)
            sage: Bmu2 == B
            True

            sage: Qmu2 = ClusterQuiver(Bmu2)
            sage: Qmu2.mutation_type()
            'undetermined finite mutation type'

        """
        return self._M

    def arcs(self):
        """
        Return the labels of diagonals (not boundary edges) given by user
        """
        return self._arcs

    def boundary_edges(self):
        """
        Return the labels of boundary edges (not diagonals) given by user
        """
        return self._boundary_edges

    def cluster(self):
        """
        Return the cluster corresponding to the initial triangulation
        """
        return self._cluster

    def boundary_edges_vars(self):
        """
        Return the labels b_i of the bounday edges specified by user
        """
        return self.boundary_edges_vars

    def triangulation_dictionary(self):
        """
        Return the correspondence between user-given labels (integers) and variables
        """
        return self._triangulation_dictionary

    def weighted_triangulation(self):
        """
        Return the triangulation given by user with weights (e.g. [(x1, x2, x0),(x1,x3,x5), ...)
        EXAMPLES::
            sage: T = [(1,4,7),(1,2,5),(2,0,3),(0,6,3)]
            sage: S = ClusterTriangulation(T)
            sage: S.weighted_triangulation()
            [(x1, x4, x7), (x1, x2, x5), (x2, x0, x3), (x0, x6, x3)]
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