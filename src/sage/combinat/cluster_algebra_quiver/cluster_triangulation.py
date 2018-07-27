r"""
Cluster Triangulations

A *cluster triangulation* (see [FominShapiroThurston]_) is a subset of
an ideal triangulation, i.e. a maximal collection of distinct
non-crossing arcs. An ideal triangulation of a surface with marked
points `(S, M)` is associated to a seed of a cluster algebra arising
from `(S, M)`.

REFERENCES:

.. [MSW_Positivity] Musiker - Schiffler - Williams,
   *Positivity for Cluster Algebras from Surfaces*,
   :arxiv:`0906.0748`

.. [MSW_Bases] Musiker - Schiffler - Williams,
   *Bases for Cluster Algebras from Surfaces*,
   :arxiv:`1110.4364`

.. [MW_MatrixFormulae] Musiker and Williams,
   *Matrix Formulae and Skein Relations for Cluster Algebras
   from Surfaces*,
   :arxiv:`1108.3382`

.. [FominShapiroThurston] Fomin - Shapiro - Thurston,
   *Cluster algebras and triangulated surfaces. part I: Cluster
   complexes*,
   :arxiv:`math/0608367`

.. [SchifflerThomas] Shiffler - Thomas,
   *On cluster algebras arising from unpunctured surfaces*,
   :arxiv:`abs/0712.4131`

.. [DupontThomas] Dupont - Thomas,
   *Atomic Basis in Types A and Affine A*,
   :arxiv:`1106.3758`

.. SEEALSO::

    Cluster triangulations closely interact with
    :class:`~sage.combinat.cluster_algebra_quiver.cluster_seed.ClusterSeed`,
    :class:`~sage.combinat.cluster_algebra_quiver.quiver.ClusterQuiver`

"""
from sage.combinat.cluster_algebra_quiver.cluster_seed import ClusterSeed

class ClusterTriangulation(ClusterSeed):
    r"""
    An initial *ideal triangulation* associated to a surface

    INPUT:

    - ``data`` -- can be any of the following

        * List of triangles - must be the list of 3-tuples (i.e. edge
          labels of a triangle) from an ideal triangulation (see
          Examples)

        * ClusterTriangulation

    .. TODO::

        - ``data`` -- can also be any of the following

            * Matrix - a skew-symmetrizable matrix arising from a tagged
              triangulation (not yet implemented)

            * DiGraph - must be the input data for a quiver from a tagged
              triangulation (not yet implemented)

            * List of edges - must be the edge list of a digraph for a
              quiver from a tagged triangulation (not yet
              implemented)

            * Objects that Theodosios Douvropoulos is currently
              working with (not related to cluster algebras from
              surfaces) (not yet implemented)

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
        A seed for a cluster algebra associated with an ideal triangulation of rank 4 with 4 boundary edges
        sage: ClusterSeed(T).mutation_type()
        ['A', [2, 2], 1]
        sage: T.map_label_to_variable()
        {'bd1': b4,
        'bd2': b5,
        'bd3': b6,
        'bd4': b7,
        'tau1': x0,
        'tau2': x1,
        'tau3': x2,
        'tau4': x3}

    Figure 15 of [MSW_Positivity]_ ::
    
        sage: thrice_punctured_square = [(2,2,1), (1,3,11), (3,12,4), (4,5,14), (5,6,10), (6,7,9), (8,10,9), (7,13,8)]
        sage: T = ClusterTriangulation(thrice_punctured_square, boundary_edges=[14,12,13,11])
        sage: T
        A seed for a cluster algebra associated with an ideal triangulation of rank 10 with 4 boundary edges
        sage: ClusterSeed(T).mutation_type()
        'undetermined finite mutation type from a surface'
        sage: T.map_label_to_variable()
        {1: x0*x1,
        2: x1,
        3: x2,
        4: x3,
        5: x4,
        6: x5,
        7: x6,
        8: x7,
        9: x8,
        10: x9,
        11: b10,
        12: b11,
        13: b12,
        14: b13}

        sage: twice_punctured_bigon = [('e','d','a'), ('a','r','b'), ('r','d','g'), ('g','n','b')]
        sage: T = ClusterTriangulation(twice_punctured_bigon, boundary_edges=['e','n'])
        sage: T
        A seed for a cluster algebra associated with an ideal triangulation of rank 5 with 2 boundary edges
        sage: ClusterSeed(T).mutation_type()
        'undetermined finite mutation type from a surface'

        sage: T = ClusterTriangulation ( [[4, 5, 1], [4, 3, 2], [3, 7, 2], [2, 1, 6], [1, 4, 5]], boundary_edges=[1])
        sage: T
        A seed for a cluster algebra associated with an ideal triangulation of rank 6 with 1 boundary edges

        sage: T = ClusterTriangulation ( [(4, 5, 1), (4, 3, 2), (3, 7, 2), (2, 1, 6), (1, 4, 5)])
        sage: T
        A seed for a cluster algebra associated with an ideal triangulation of rank 7

        sage: Triangles = [(2,3,11),(2,1,1),(4,3,12),(0,4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,13)]
        sage: T = ClusterTriangulation(Triangles, boundary_edges=[11,12,13,0])
        sage: T
        A seed for a cluster algebra associated with an ideal triangulation of rank 10 with 4 boundary edges

        sage: once_punctured_torus = ClusterTriangulation([(1,2,3),(3,1,2)])
        sage: S = ClusterSeed(once_punctured_torus).mutation_type()
        sage: S
        'undetermined finite mutation type from a surface'

    Triangles do not need to be connected::

        sage: T = ClusterTriangulation([(1,2,3),(3,4,5),(1,5,6),(10,20,30),(30,40,50)])
        sage: ClusterSeed(T).mutation_type()
        [ ['D', 6], ['A', 5] ]

    Surfaces that are not allowed::

        sage: one_selffolded_triangle = [(1,1,2)]
        sage: ClusterTriangulation(one_selffolded_triangle)
        Traceback (most recent call last):
        ...
        ValueError: The following surfaces are not allowed: a sphere with 1, 2, or 3 punctures; a monogon with zero or 1 puncture; a bigon or triangle without punctures.

        sage: one_triangle = [(0,1,2)]
        sage: ClusterTriangulation(one_triangle)
        Traceback (most recent call last):
        ...
        ValueError: The following surfaces are not allowed: a sphere with 1, 2, or 3 punctures; a monogon with zero or 1 puncture; a bigon or triangle without punctures.

        sage: thrice_punctured_sphere = [(0,1,2),(2,1,0)]
        sage: ClusterTriangulation(thrice_punctured_sphere)
        Traceback (most recent call last):
        ...
        ValueError: The following surfaces are not allowed: a sphere with 1, 2, or 3 punctures; a monogon with zero or 1 puncture; a bigon or triangle without punctures.

        sage: selffolded_triangle_and_ordinary_triangle = [(0,0,1),(2,3,4)]
        sage: ClusterTriangulation(selffolded_triangle_and_ordinary_triangle)
        Traceback (most recent call last):
        ...
        ValueError: A noose of a self-folded triangle must be a side of another triangle.
    """
    def __init__(self, data, frozen=None, is_principal=None, from_surface=False, boundary_edges=None):
        r"""
        .. TODO::

            See my data for the TestSuite. Is this good enough?

        TESTS::

            sage: CT = ClusterTriangulation([('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')], boundary_edges=['c','d','e','f'])
            sage: TestSuite(CT).run()
        """
        from sage.combinat.cluster_algebra_quiver.surface import remove_duplicate_triangles, _triangulation_to_arrows, _surface_edge_list_to_matrix, _get_user_arc_labels, produce_dict_label_to_variable, produce_dict_variable_to_label, _get_user_label_triangulation, _get_weighted_triangulation
        from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        from sage.rings.all import QQ
        from sage.rings.all import FractionField, PolynomialRing
        from copy import copy
        from sage.matrix.all import identity_matrix

        self._from_surface = True # Maybe try to remove this? This should not be necessary for class ClusterTriangulation

        # if data is a list of triangles (forming a triangulation)
        if isinstance(data,list) and \
        all(type(triangle) in [list,tuple] and len(triangle)==3 for triangle in data):
            self._boundary_edges = boundary_edges if boundary_edges else []
            self._edges = _get_user_arc_labels(data)
            self._arcs = list(self._edges)
            for b in self._boundary_edges:
                self._arcs.remove(b)

            if not set(self._boundary_edges).issubset(self._edges):
                raise ValueError(boundary_edges, ' are not a subset of ', self._edges, ' .Optional parameter boundary_edges must be a list of edges that are part of the triangulation')
            data = remove_duplicate_triangles (data, boundary_edges)

            self._n = len(self._edges) - len(self._boundary_edges) #if boundary_edges else len(self._edges)
            self._nlist = range(self._n)
            self._triangles = data

            all_arrows = _triangulation_to_arrows(self._triangles)
            M = _surface_edge_list_to_matrix (all_arrows, self._edges, self._boundary_edges, self._n)
            self._M = M[:self._n,:self._n]

            #self._is_cluster_algebra = True
            self._description = 'A seed for a cluster algebra associated with an ideal triangulation of rank %d'  %(self._n)
            if boundary_edges:
                self._description += ' with %d boundary edges' %(len(self._boundary_edges))

            if is_principal:
                #if self._M.nrows()==self._n:
                #    self._M = self._M.stack(identity_matrix(self._n))
                self._M = self._M.stack(identity_matrix(self._n))
                #if 'with principal coefficients' not in data._description:
                self._description += ' with principal coefficients'

            self._m = self._M.nrows() - self._n
            self._mlist = range(self._m)
            #print 'm :', self._m
            self._quiver = ClusterQuiver(self._M, from_surface=True)

            if is_principal:
                self._R = FractionField(PolynomialRing(QQ,['x%s'%i for i in range(0,self._n)]+['y%s'%i for i in range(0,self._n)]+['b%s'%i for i in range(self._n,len(self._boundary_edges)+self._n)]))
            else:
                self._R = FractionField(PolynomialRing(QQ,['x%s'%i for i in range(0,self._n)]+['b%s'%i for i in range(self._n,len(self._boundary_edges)+self._n)]))
            #self._boundary_edges_vars = FractionField(PolynomialRing(QQ,['b%s'%i for i in range(self._n,len(self._boundary_edges)+self._n)]))
            #self._cluster = list(self._R.gens()[0:self._n])
            self._cluster = list(self._R.gens()[0:self._n+self._m])
            self._boundary_edges_vars = list(self._R.gens()[self._n+self._m:]) if boundary_edges else []
            self._map_label_to_variable = produce_dict_label_to_variable(self._triangles, self._cluster[0:self._n], self._boundary_edges, self._boundary_edges_vars)
            self._map_variable_to_label = produce_dict_variable_to_label(self._map_label_to_variable)
            self._triangulation = _get_user_label_triangulation(self._triangles)
            self._weighted_triangulation = _get_weighted_triangulation (self._triangles, self._map_label_to_variable)
            self._mutation_type = self._quiver.mutation_type()
            if self._mutation_type == 'undetermined finite mutation type':
                self._mutation_type += ' from a surface'
            self._is_principal = is_principal

			#sets appropriate booleans to 'False'            
            self._use_g_vec = False
            self._use_c_vec = False
            self._use_d_vec = False
            self._use_fpolys = False
            self._bot_is_c = False
            
            #sets up ability to track mutations
            self._mut_path = []
            self._track_mut = True
            
            #sets currently unused dictionaries, names, data to 'None'
            self._user_labels = None
            self._user_labels_prefix = None
            self._init_vars = None
            self._init_exch = None
            self._F = None
            self._y = None
            self._yhat = None
            #self._mut_path = None
            
            #The initial B-matrix is set to be the B-matrix corresponding to the input triangulation
            self._b_initial = copy(self._M)
            
            #Constructs the appropriate coefficient ring
            self._U = PolynomialRing(QQ,['y%s' % i for i in range(self._n)])

        # Construct data from a cluster triangulation
        elif isinstance(data, ClusterTriangulation):
            self._boundary_edges = copy(data._boundary_edges)
            self._edges = copy(data._edges)
            self._arcs = copy(data._arcs)
            self._n = data._n
            self._nlist = range(self._n)
            self._triangles = copy( data._triangles )
            self._M = copy(data._M)
            self._m = self._M.nrows() - self._n
            self._mlist = range(self._m)
            #self._quiver = ClusterQuiver(data._quiver)
            self._quiver = ClusterQuiver( data._quiver ) if data._quiver else None
            #self._is_cluster_algebra = data._is_cluster_algebra
            self._description = copy( data._description )

            self._R = copy(data._R)
            self._cluster = copy(data._cluster)
            self._boundary_edges_vars = copy(data._boundary_edges_vars)
            self._map_label_to_variable = copy(data._map_label_to_variable)
            self._map_variable_to_label = copy(data.map_variable_to_label)
            self._triangulation = copy(data._triangulation)
            self._weighted_triangulation = copy(data._weighted_triangulation)
            self._mutation_type = data._mutation_type
            self._is_principal = copy(data._is_principal)
            
			#sets appropriate booleans to 'False'            
            self._use_g_vec = False
            self._use_c_vec = False
            self._use_d_vec = False
            self._use_fpolys = True
            self._bot_is_c = False
            
            #sets up ability to track mutations
            self._mut_path = []
            self._track_mut = True
            
            #sets currently unused dictionaries, names, data to 'None'
            self._user_labels = None
            self._user_labels_prefix = None
            self._init_vars = None
            self._init_exch = None
            self._F = None
            self._y = None
            self._yhat = None
            
            #The initial B-matrix is set to be the B-matrix corresponding to the input ClusterTriangulation
            self._b_initial = copy(self._M)
            
            #Constructs the appropriate coefficient ring
            self._U = PolynomialRing(QQ,['y%s' % i for i in range(self._n)])
        elif isinstance(data,ClusterSeed):
        	#copy data that we want to retain in the ClusterTriangulation
        	#(some of this is not necessary, so we might want to just ignore it)
        	self._n = data._n
        	self._nlist = data._nlist
        	self._m = data._m
        	self._mlist = data._mlist
        	self._M = data._M
        	self._R = data._R
        	self._cluster = copy(data._cluster)
        	self._mutation_type = data._mutation_type
        	self._is_principal = data._is_principal
        	self._use_g_vec = data._use_g_vec
        	self._use_c_vec = data._use_c_vec
        	self._use_d_vec = data._use_d_vec
        	self._use_fpolys = data._use_fpolys
        	self._bot_is_c = data._bot_is_c
        	self._mut_path = data._mut_path
        	self._track_mut = data._track_mut
        	self._user_labels = data._user_labels
        	self._user_labels_prefix = data._user_labels_prefix
        	self._init_vars = data._init_vars
        	self._init_exch = data._init_exch
        	self._F = data._F
        	self._y = data._y
        	self._yhat = data._yhat
        	
        	if data._mutation_type.letter() == 'A':
        		#By default, I'll have all the boundary edge variables be b1, b2, etc
        		#And the default user labels for boundaries will just be 1?
        		# Maybe there's another choice that makes more sense
        		self._boundary_edges = ['b%i' % i for i in range(self._n+3)]
        		self._boundary_edges_vars = [1]*(self._n+3)
        		#print("It's type A!")
        		
        		#Constructing a list of triangles that we want in our triangulation
        		#Note that the triangulation comes from the initial cluster, so we use the initial B matrix
        		B = data._b_initial
        		triangles = []
        		avail_boundary_edges = self._boundary_edges
        		
        		#print(avail_boundary_edges)
        		
        		#The first triangle will always have sides (0,b_0,b_1)
        		triangles += [(avail_boundary_edges[0], avail_boundary_edges[1], 0)]
        		avail_boundary_edges = avail_boundary_edges[2:]
        		
        		#print(avail_boundary_edges)
        		#print(triangles)
        		#print(avail_boundary_edges)
        		
        		#creates intermediate triangles based on quiver orientation
        		for i in range(self._n)[:-1]:
        			if B[i,i+1] > 0:
        				#print("+")
        				triangles += [(i+1,i,avail_boundary_edges[0])]
        				avail_boundary_edges = avail_boundary_edges[1:]
        				#print(triangles)
        				#print(avail_boundary_edges)
        			if B[i,i+1] < 0:
        				#print("-")
        				triangles += [(i,i+1,avail_boundary_edges[-1])]
        				avail_boundary_edges = avail_boundary_edges[:-1]
        				#print(triangles)
        				#print(avail_boundary_edges)
        			
        		#Adds final triangle with side n and remaining boundary sides	        		
        		triangles += [(self._n-1,avail_boundary_edges[0],avail_boundary_edges[1])]
        		
        		print triangles
        		#print(avail_boundary_edges)
        		
        	if data._mutation_type.letter() == 'D':
        		#Same default choices, same comment about there probably being a better choice
        		self._boundary_edges = ['b%i' % i for i in range(self._n)]
        		self._boundary_edges_vars = [1]*(self._n)
        		#print("It's type D!")
        		
        		B = data._b_initial
        		triangles = []
        		avail_boundary_edges = self._boundary_edges
        		
        		#The first triangle always has sides (0,b_0,b_1)
        		triangles += [(avail_boundary_edges[0],avail_boundary_edges[1],0)]
        		avail_boundary_edges = avail_boundary_edges[2:]
        		
        		#creates intermediate vertices (prior to last two)
        		
        		for i in range(self._n)[:-3]:
        			if B[i,i+1] > 0:
        				triangles += [(i+1,i,avail_boundary_edges[0])]
        				avail_boundary_edges = avail_boundary_edges[1:]
        			if B[i,i+1] < 0:
        					triangles += [(i,i+1,avail_boundary_edges[-1])]
        					avail_boundary_edges = avail_boundary_edges[:-1]
        					
        		
        		#creates the final two triangles (possibly self-folded)
        		#there are four cases, based on the four possible orientations for this part of D_n
        		if B[self._n -2, self._n - 3] > 0 and B[self._n -1, self._n -3] > 0:
        			triangles += [(self._n - 2, self._n - 3, avail_boundary_edges[0]),(self._n - 2, self._n -1, self._n -1)]
        		elif B[self._n - 2, self._n -3] > 0 and B[self._n -1, self._n -3] < 0:
        			triangles += [(self._n -1, self._n -2, self._n -3), (self._n - 1, self._n -2, avail_boundary_edges[0])]
        		elif B[self._n -2, self._n -3] < 0 and B[self._n - 1, self._n - 3] > 0:
        			triangles += [(self._n - 3, self._n - 1, self._n - 2),(self._n - 1, avail_boundary_edges[0], self._n - 2)]
        		else:
        			triangles += [(self._n - 3, self._n - 2, avail_boundary_edges[0]),(self._n - 2, self._n - 2, self._n -1)]
        					
        		print triangles
        		#print avail_boundary_edges
        	
        	
        	#print("Input was a ClusterSeed that's NOT a ClusterTriangulation")
        	
        else:
            raise ValueError('Input must be a list of three-tuples or a ClusterTriangulation class. You entered data: ', data)

        #self._m = 0
        #self._quiver = ClusterQuiver(self._M)
        #if not self._quiver:
        #    self._quiver = ClusterQuiver(self._M, from_surface=True)
        #ClusterSeed.__init__(self, self._quiver)
        ##ClusterSeed.__init__(self, self._quiver, from_surface=True)
        #ClusterSeed.__init__(self, self, from_surface=True, is_principal=is_principal)
        #super(ClusterSeed, self).__init__(self._M)

    def __eq__(self, other):
        r"""
        Return ``True`` iff ``self`` represent the same cluster triangulation
        as ``other``.

        EXAMPLES::

            sage: CT = ClusterTriangulation([('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')], boundary_edges=['c','d','e','f'])
            sage: CT1 = ClusterTriangulation([('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')])
            sage: CT.__eq__(CT1)
            False
            sage: ClusterTriangulation(CT).__eq__(CT)
            True

            sage: CT2 = ClusterTriangulation([('aa','dd','cc'), ('aa','ll','bb'), ('rr','rr','ll'),('bb','ff','ee')], boundary_edges=['ee','dd','cc','ff'])
            sage: CT.__eq__(CT2)
            True
            sage: ClusterSeed(CT) == ClusterSeed(CT2)
            True
        """
        from sage.combinat.cluster_algebra_quiver.surface import _rearrange_triangle_small_first

        Tself, Tother = [], []
        if not isinstance(other, ClusterTriangulation):
            return False
        else:
            for triangle in self._weighted_triangulation:
                Tself.append(_rearrange_triangle_small_first(triangle))
            for triangle in other._weighted_triangulation:
                Tother.append(_rearrange_triangle_small_first(triangle))
            Tself.sort()
            Tother.sort()
            return self._M == other._M and Tself == Tother

    def _repr_(self):
        """
        Return the description of ``self``.

        EXAMPLES::

            sage: T = ClusterTriangulation([(4, 5, 1), (4, 3, 2), (3, 7, 2), (2, 1, 6), (1, 4, 5)])
            sage: T._repr_()
            'A seed for a cluster algebra associated with an ideal triangulation of rank 7'
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
            [ 0 -1  1  0  0]
            [ 1  0  0 -1 -1]
            [-1  0  0  1  1]
            [ 0  1 -1  0  0]
            [ 0  1 -1  0  0]

        2 self-folded triangles and 1 triangle with one vertex
        (affine D) Figure 9 (right) of [FominShapiroThurston]_ ::

            sage: twice_punctured_monogon_mu2 = [(4,5,5),(2,3,3),(1,4,2)]
            sage: Tmu2 = ClusterTriangulation(twice_punctured_monogon_mu2)
            sage: Bmu2 = Tmu2.b_matrix()
            sage: Bmu2
            [ 0  1  1 -1 -1]
            [-1  0  0  1  1]
            [-1  0  0  1  1]
            [ 1 -1 -1  0  0]
            [ 1 -1 -1  0  0]
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
            [ 0  0 -1 -1]
            [ 0  0 -1 -1]
            [ 1  1  0  0]
            [ 1  1  0  0]

        Figure 10 (top) of [FominShapiroThurston]_ ::

            sage: twice_punctured_monogon_mu3 = [(1,1,2), (2,4,3), (3,4,'boundary')]
            sage: Tmu3 = ClusterTriangulation(twice_punctured_monogon_mu3, boundary_edges=['boundary'])
            sage: Bmu3 = Tmu3.b_matrix()
            sage: Bmu3
            [ 0  0  1 -1]
            [ 0  0  1 -1]
            [-1 -1  0  0]
            [ 1  1  0  0]
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
            sage: Tmu0402.triangulation() == \
            ....: ClusterTriangulation([(1,0,4),(4,5,5),(1,2,0),(3,2,3)]).triangulation()
            True
        """
        return self._M

    def mutate(self, sequence, inplace=True, user_labels=True):
        """
        Mutate ``self`` at a diagonal or a sequence of diagonals.

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

        Twice-punctured monogon (labeled 'i0') with 3 (non-ordinary) ideal triangles (affine D)::

            sage: T = ClusterTriangulation([('i1','i4','i2'),('i3','i4','i3'),('i2','i0','i1')])
            sage: T.mutate('i1', inplace=False).cluster()
            [x0, (x3*x4 + x0)/x1, x2, x3, x4]
            sage: T.mutate('i2', inplace=False).cluster()
            [x0, x1, (x3*x4 + x0)/x2, x3, x4]

        Mutating at a self-folded triangle (r,r,ell)'s radius r is an involution except that
        mutating at r twice returns the original triangulation with the label r and ell switched::

            sage: T.mutate('i3', inplace=False).cluster()
            [x0, x1, x2, (x1 + x2)/x3, x4]
            sage: ClusterSeed(T._M).mutate(T.get_edge_position('i3'), inplace=False).cluster() == T.mutate('i3', inplace=False).cluster()
            True
            sage: T.mutate(['i3','i3'], inplace=False).cluster()
            [x0, x1, x2, x3, x4]
            sage: T.mutate(['i3','i3'], inplace=False).weighted_triangulation()
            [(x2, x0, x1),
            (x1, x3*x4, x2),
            ((x4, 'counterclockwise'), (x4, 'clockwise'), x3*x4)]
            sage: T.weighted_triangulation()
            [(x1, x3*x4, x2),
            ((x3, 'counterclockwise'), (x3, 'clockwise'), x3*x4),
            (x2, x0, x1)]
            sage: T.mutate(['i4','i4'], inplace=False) == T
            True

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
            sage: T.mutate('i0',inplace=False).quiver() == Tmu04020.quiver()
            True

        A once-punctured square's triangulation with self-folded
        triangle, border edges are labeled 4,5,6,7, 2nd triangulation
        in oral paper ell-loop is labeled 3, radius is labeled 0::

            sage: T = ClusterTriangulation([(1,7,4),(1,5,2),(2,3,6),(3,0,0)], boundary_edges=[4,5,6,7])
            sage: S = ClusterSeed(T); S.cluster()
            [x0, x1, x2, x3]
            sage: S.mutate(T.get_edge_position(0))
            sage: S.cluster()
            [(x2 + 1)/x0, x1, x2, x3]
            sage: T.mutate(0)
            sage: T.triangles()
            [(1, 7, 4), (1, 5, 2), (2, 3, 0), (0, 3, 6)]
            sage: T.cluster()
            [(x2 + 1)/x0, x1, x2, x3]
            sage: ClusterSeed(T).cluster() # Cluster is reset to initial cluster
            [x0, x1, x2, x3]
            sage: T.quiver() == S.quiver()
            True

        A once-punctured torus::

            sage: T = ClusterTriangulation([('a','b','c'),('c','a','b')])
            sage: S = ClusterSeed(T)
            sage: S.mutate(T.get_edge_position('a'))
            sage: T.mutate('a')
            sage: T.triangulation()
            [('b', 'a', 'c'), ('b', 'a', 'c')]
            sage: T.cluster()
            [(x1^2 + x2^2)/x0, x1, x2]
            sage: S.cluster() == T.cluster()
            True
            sage: S.b_matrix() == T.b_matrix()
            True
            sage: S.quiver() == T.quiver()
            True
            sage: T.mutate('a')
            sage: T == ClusterTriangulation([('a','b','c'),('c','a','b')])
            True

        An edge that is only contained in one triangle and is not a
        self-folded triangle's radius is not mutable::

            sage: twice_punc_monogon = ClusterTriangulation([('i1','i4','i2'),('i3','i4','i3'),('i2','i0','i1')])
            sage: twice_punc_monogon.mutate('i0')
            Traceback (most recent call last):
            ...
            ValueError: ('The ideal triangulation cannot be mutated at ', 'i0', '.There is only one triangle ', ('i2', 'i0', 'i1'), ', not a self-folded triangle, with side ', 'i0')
        """
        from sage.combinat.cluster_algebra_quiver.surface import _triangles_mutate, \
        produce_dict_label_to_variable, produce_dict_variable_to_label, _get_user_label_triangulation, _get_weighted_triangulation

        if inplace:
            ct = self
        else:
            ct = ClusterTriangulation( self )

        S = ClusterSeed(ct)

        #print "I am in ClusterTriangulation.mutate with ct: ", ct # TODO ERASE

        n = ct._n
        V = range(n)
        edges = ct._edges

        if user_labels:
            if sequence in edges:
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
            raise ValueError('The second parameter must be boolean. To mutate at a sequence of length 2, input list of diagonals.')
        if not isinstance(user_labels, bool):
            raise ValueError('The third parameter must be boolean. To mutate at a sequence of length 3, input a list of diagonals.')

        if user_labels:
            if any( a not in edges for a in seq ):
                a = filter( lambda a: a not in edges, seq )[0]
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
                #ct._M.mutate(pos)
                #ct._quiver.mutate(pos)
                #S = S.mutate(pos, inplace=False)
            else:
                pos = diagonal
                diagonal_label = ct._arcs[pos]
                ct._triangles = _triangles_mutate(ct._triangles, diagonal_label)
                #ct._M.mutate(pos)
                #ct._quiver.mutate(pos)
                #S = S.mutate(pos, inplace=False)

            ct._M.mutate(pos)
            S = S.mutate(pos, inplace=False)

        ct._cluster = S.cluster()
        ct._map_label_to_variable = produce_dict_label_to_variable (ct._triangles, ct._cluster[0:ct._n], ct._boundary_edges, ct._boundary_edges_vars)
        ct._map_variable_to_label = produce_dict_variable_to_label (ct._map_label_to_variable)
        ct._triangulation = _get_user_label_triangulation(ct._triangles)
        ct._weighted_triangulation = _get_weighted_triangulation (ct._triangles, ct._map_label_to_variable)

        #ct._M = S._M
        ct._quiver = None

        if not inplace:
            return ct

    def arcs(self):
        """
        Return the sorted list of labels of diagonals of
        ``self`` given by user.

        EXAMPLES::

            sage: annulus22 = [('bd1','tau1','tau2'),('tau2','tau3','bd4'),('tau1','tau4','bd2'),('tau3','bd3','tau4')]
            sage: T = ClusterTriangulation(annulus22, boundary_edges=['bd3','bd2','bd1','bd4'])
            sage: T.arcs()
            ['tau1', 'tau2', 'tau3', 'tau4']
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

    #def triangulation_dictionary(self):
    def map_label_to_variable(self):
        """
        Return a dictionary with keys user-given labels (numbers or strings)
        and items variables x_i/b_i of ``self``.

        EXAMPLES:

         2 self-folded triangles and 1 triangle with one vertex (affine D),
         Figure 10 (bottom) of [FominShapiroThurston]_::

            sage: twice_punctured_monogon = [(4,5,5),(2,3,3),(1,4,2)]
            sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=[1])
            sage: T.map_label_to_variable()
            {1: b4, 2: x0*x1, 3: x1, 4: x2*x3, 5: x3}
            sage: T.mutate(3)
            sage: T.cluster()
            [x0, (x2*x3 + 1)/x1, x2, x3]
            sage: T.map_label_to_variable()
            {1: b4, 2: x0, 3: (x2*x3 + 1)/x1, 4: x2*x3, 5: x3}
        """
        return self._map_label_to_variable

    #def triangulation_dictionary_variable_to_label(self):
    def map_variable_to_label(self):
        """
        Return a dictionary with keys variables x_i/b_i
        and items user-given labels (numbers or strings) of ``self``.

        EXAMPLES:

         2 self-folded triangles and 1 triangle with one vertex (affine D),
         Figure 10 (bottom) of [FominShapiroThurston]_::

            sage: twice_punctured_monogon = [(4,5,5),(2,3,3),(1,4,2)]
            sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=[1])
            sage: T.map_variable_to_label()
            {b4: 1, x3: 5, x1: 3, x2*x3: 4, x0*x1: 2}
        """
        return self._map_variable_to_label

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

    def _get_map_label_to_variable(self,a):
        """
        Return the variable corresponding to the label (given by user)
        of an arc of boundary edge.

        .. SEEALSO::

            :meth:`get_edge_position`

        EXAMPLES::

            sage: T = ClusterTriangulation([(1,7,4),(1,5,2),(6,0,3),(2,3,0),(0,3,6),[7,4,1]], boundary_edges=[4,5,6,7])
            sage: T._get_map_label_to_variable(0)
            x0
            sage: T._get_map_label_to_variable(6)
            b6
            sage: once_punctured_square = [('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')]
            sage: T = ClusterTriangulation(once_punctured_square, boundary_edges=['c','f','e','d'])
            sage: T._get_map_label_to_variable('r')
            x3
            sage: T._get_map_label_to_variable('ll')
            x2*x3
            sage: T.mutate('a')
            sage: T._get_map_label_to_variable('a')
            (x2*x3 + x1)/x0
        """
        from sage.combinat.cluster_algebra_quiver.surface import _get_weighted_edge
        return _get_weighted_edge(a, self._map_label_to_variable)

    def get_edge_position(self, a):
        """
        Return the position of the variable corresponding to the label
        (given by user) of an arc or a boundary edge.

        .. SEEALSO::

            :meth:`_get_map_label_to_variable`

        EXAMPLES::

            sage: T = ClusterTriangulation([('i1','i7','i4'),('i1','i5','i2'),\
            ('i6','i0','i3'),('i2','i3','i0'),('i0','i3','i6'),['i7','i4','i1']], \
            boundary_edges=['i4','i5','i6','i7'])
            sage: T.get_edge_position('i0')
            0
            sage: T.get_edge_position('i6')
            2
            sage: once_punctured_square = [('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')]
            sage: T = ClusterTriangulation(once_punctured_square, boundary_edges=['c','f','e','d'])
            sage: T.get_edge_position('r')
            3
            sage: T.get_edge_position('c')
            0
        """
        from sage.combinat.cluster_algebra_quiver.surface import _get_weighted_edge
        arcs = list(self._arcs)
        boundary_edges = list(self._boundary_edges)
        if a in arcs:
            return arcs.index(a)
        elif a in boundary_edges:
            return boundary_edges.index(a)
        raise ValueError(a, " is not a user-given label of an arc/boundary edge.")

    def _get_map_variable_to_label(self,var):
        """
        Return the label (given by user) of an arc of boundary edge
        corresponding to a variable x_i or b_i (or a product x_i*x_j).

        INPUT:

        - ``var`` -- a variable or a product of two variables from
          self.map_label_to_variable()

        EXAMPLES::

            sage: T = ClusterTriangulation([(1,7,4),(1,5,2),(6,0,3),(2,3,0),(0,3,6),[7,4,1]], boundary_edges=[4,5,6,7])
            sage: T._get_map_variable_to_label(T._cluster[0])
            0
            sage: T._get_map_variable_to_label(T._boundary_edges_vars[2])
            6
            sage: T.mutate(0)
            sage: T.cluster()[0]*T.cluster_variable(3)
            (x2*x3 + x3)/x0
            sage: T._get_map_variable_to_label(T.cluster()[0]*T.cluster_variable(3))
            0
            sage: T._get_map_label_to_variable(0)
            (x2*x3 + x3)/x0
            sage: T.cluster()
            [(x2 + 1)/x0, x1, x2, x3]

            sage: TT = ClusterTriangulation([('j1','j1','j2'),('j3','j4','j3'),('j2','j4','j0')])
            sage: TT._get_map_variable_to_label(TT._cluster[1]*TT._cluster[2])
            'j2'
        """
        from sage.combinat.cluster_algebra_quiver.surface import _get_edge_user_label
        edge_user_label = _get_edge_user_label(var,self._map_variable_to_label)
        if edge_user_label is not None:
            return edge_user_label
        else:
            raise ValueError(var, ' is not a cluster variable (or a product of cluster variables) from self.map_variable_to_label():', self._map_variable_to_label)

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
            sage: T.arc_laurent_expansion([c[0],c[1],c[2],c[3],c[0]], user_labels=False) == T.arc_laurent_expansion([c[0],c[3],c[2],c[1],c[0]], user_labels=False)
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
            sage: T.arc_laurent_expansion([ell,(r,'counterclockwise'), ell, c[3-1],c[4-1],c[5-1],c[6-1]], user_labels=False)== T.arc_laurent_expansion([c[6-1],c[5-1],c[4-1],c[3-1], ell,(r,'clockwise'), ell], user_labels=False) # Gamma_1
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

            sage: T.arc_laurent_expansion([c[4],c[5],c[6],c[7],c[8],c[5],c[4]], user_labels=False) == T.arc_laurent_expansion([c[4],c[5],c[8],c[7],c[6],c[5],c[4]], user_labels=False) # Ell_p
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
        else:
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
        else:
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
            sage: r = T._get_map_label_to_variable('r')
            sage: ell = T._get_map_label_to_variable('ell')
            sage: three = T._get_map_label_to_variable(3)
            sage: four = T._get_map_label_to_variable(4)
            sage: five = T._get_map_label_to_variable(5)
            sage: six = T._get_map_label_to_variable(6)
            sage: crossed_arcs = [ell, (r,'counterclockwise'), ell, three, four, five, six]
            sage: T.draw_snake_graph(crossed_arcs,first_tile_orientation=-1,user_labels=False)
            Graphics object consisting of 43 graphics primitives
            sage: T.draw_snake_graph(['ell', ('r','counterclockwise'), 'ell', 3, 4, 5, 6],first_tile_orientation=-1,user_labels=True)
            Graphics object consisting of 43 graphics primitives
            sage: T.draw_snake_graph([5,6,7,8,9,6,5], first_tile_orientation=1, user_labels=True)
            Graphics object consisting of 43 graphics primitives
        """
        from sage.combinat.cluster_algebra_quiver.surface import _draw_snake_graph
        drawing = _draw_snake_graph (self.list_snake_graph(crossed_arcs, first_triangle=first_triangle, final_triangle=final_triangle, first_tile_orientation=first_tile_orientation, user_labels=user_labels), print_user_labels=user_labels)
        #return drawing.plot( axes=False, figsize=fig_size )
        return drawing.plot()

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
            Graphics object consisting of 25 graphics primitives
            sage: T.draw_band_graph([1,2,3,0,1], user_labels=True)
            Graphics object consisting of 25 graphics primitives
        """
        from sage.combinat.cluster_algebra_quiver.surface import _draw_snake_graph
        drawing = _draw_snake_graph(self.list_band_graph(crossed_arcs, first_triangle=first_triangle, final_triangle=final_triangle, first_tile_orientation=first_tile_orientation, user_labels=user_labels), print_user_labels=user_labels)
        return drawing.plot( axes=False, figsize=fig_size )

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
            Graphics object consisting of 64 graphics primitives
            sage: T.draw_lifted_arc ([2,(1,'counterclockwise'),2,3,4,5,6], user_labels=True)
            Graphics object consisting of 64 graphics primitives

        Figure 8 of Musiker - Schiffler - Williams "Bases for Cluster
        Algebras from Surfaces" [MSW_Bases]_ where: the loop gamma
        crosses the arcs labeled 1,2,3,4,1 (1 is counted twice) outer
        boundary edges, clockwise starting from starting point of
        gamma: 7,8,9, 10 inner boundary edges, clockwise starting from
        ending point of gamma:11, 0::

            sage: T = ClusterTriangulation([(8,7,5),(5,4,1),(1,0,2),(2,11,3),(3,4,6),(6,10,9)], boundary_edges=[7,8,9,10,11,0])
            sage: c = [item for item in T.cluster()]
            sage: T.draw_lifted_arc([c[1-1],c[2-1],c[3-1],c[4-1],c[1-1]], user_labels=False)
            Graphics object consisting of 48 graphics primitives
            sage: T.draw_lifted_arc([1,2,3,4,1], user_labels=True)
            Graphics object consisting of 48 graphics primitives
        """
        from sage.combinat.cluster_algebra_quiver.surface import _lifted_polygon, _draw_lifted_curve

        if user_labels:
            lifted_polygon = _lifted_polygon(self._triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=True, is_loop=False)
        else:
            lifted_polygon = _lifted_polygon(self._weighted_triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=True, is_loop=False)
        drawing = _draw_lifted_curve(lifted_polygon, is_arc=True, is_loop=False)
        #if verbose:
        #    print lifted_polygon
        return drawing.plot(axes=False, figsize=fig_size)

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
            Graphics object consisting of 120 graphics primitives
            sage: crossed_userlabels = [3,2,(1,'counterclockwise'),2, 3,4,5,6,7,8,9,6,5,4,3]
            sage: T.draw_lifted_loop(crossed_userlabels)
            Graphics object consisting of 120 graphics primitives

        Figure 8 of Musiker - Schiffler - Williams "Bases for Cluster Algebras from Surfaces" [MSW_Bases]_
        with triangulation having arcs 1, ..., 6, and gamma crosses 1,2,3,4,1 (1 is listed twice); and
        the outer boundary edges, clockwise starting from starting point of gamma are 7,8,9, 10; and
        the inner boundary edges, clockwise starting from ending point of gamma:11, 0::

            sage: T = ClusterTriangulation([(8,7,5),(5,4,1),(1,0,2),(2,11,3),(3,4,6),(6,10,9)], boundary_edges=[7,8,9,10,11,0])
            sage: c = [item for item in T.cluster()]
            sage: T.draw_lifted_loop([c[1-1],c[2-1],c[3-1],c[4-1],c[1-1]], user_labels=False)
            Graphics object consisting of 40 graphics primitives
            sage: T.draw_lifted_loop([1,2,3,4,1], user_labels=True)
            Graphics object consisting of 40 graphics primitives
        """
        from sage.combinat.cluster_algebra_quiver.surface import _lifted_polygon, _draw_lifted_curve
        if user_labels:
            lifted_polygon = _lifted_polygon(self._triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=False, is_loop=True)
        else:
            lifted_polygon = _lifted_polygon(self._weighted_triangulation, crossed_arcs, first_triangle, final_triangle, is_arc=False, is_loop=True)
        drawing = _draw_lifted_curve(lifted_polygon, is_arc=False, is_loop=True)
        # if verbose:
        #    print lifted_polygon
        return drawing.plot( axes=False, figsize=fig_size)

    def arc_laurent_expansion(self, crossed_arcs, first_triangle=None,
                              final_triangle=None, verbose=False,
                              fig_size=1, user_labels=True, return_labels=False):
        """Return the Laurent expansion of a generalized arc gamma
        (i.e. a curve between marked point/s) which crosses the arc in
        crossed_arcs

        INPUT:

        - ``crossed_arcs`` -- labels from
          self.cluster_triangulation().triangulation() (if
          ``user_labels`` is ``True``) and variables from
          ``self.cluster_triangulation().cluster()`` if
          (``user_labels`` is ``False``) corresponding to arcs that
          are crossed by gamma If curve crosses a self-folded triangle
          (ell,r,ell), then specify ``ell, (r, 'counterclockwise),
          ell`` (if, as gamma is about to cross r, the puncture is to
          the right of gamma) or ``ell, (r, clockwise), ell`` (if, as
          gamma is about to cross r, the puncture is to the left of
          gamma)

        - ``first_triangle`` -- (default:``None``) the first triangle
          crossed by curve

        - ``final_triangle`` -- (default:``None``) the last triangle
          crossed by curve

        - ``verbose`` -- (default:``False``) display the image of the
          perfect matchings of the snake graph if ``True``

        - ``fig_size`` -- (default:1) image size

        - ``user_labels`` -- (default:``True``) whether or not
          ``crossed_arcs`` is a list of labels
          
        - ``return_labels`` -- (default: ``False``) allows the Laurent expansion
          to be returned in terms of labels rather than variables

        .. SEEALSO::

            :meth:`loop_laurent_expansion`,
            :meth:`ClusterTriangulation.draw_lifted_arc`

        ALGORITHM:

        See the perfect matching formula from in [MSW_Positivity]_ (section 4).

        EXAMPLES:

        An ideal triangulation of a once-punctured square having 2 radii, and the boundary edges are labeled 4,5,6,7::

            sage: once_punctured_square = [(1,7,4),(1,5,2),(6,0,3),(2,3,0),(0,3,6),[7,4,1]]
            sage: T = ClusterTriangulation(once_punctured_square, boundary_edges=[4,5,6,7])
            sage: Q = ClusterQuiver(T)
            sage: S = ClusterSeed(Q)
            sage: c = [item for item in T.cluster()]
            sage: T.arc_laurent_expansion([S.x(1),S.x(2),S.x(3)], user_labels=False)
            (x0*x2^2 + 2*x0*x2 + x1*x3 + x0)/(x1*x2*x3)
            sage: T.arc_laurent_expansion([c[1],c[2],c[3]], first_triangle=[c[1],T._get_map_label_to_variable(7),T._get_map_label_to_variable(4)], \
            ....: final_triangle=( c[0],c[3], T._get_map_label_to_variable(6) ), \
            ....: user_labels=False) == T.arc_laurent_expansion([S.x(1),S.x(2),S.x(3)], user_labels=False)
            True
            sage: T.arc_laurent_expansion([1,2,3],user_labels=True) == T.arc_laurent_expansion([S.x(1),S.x(2),S.x(3)],user_labels=False)
            True
            sage: T.arc_laurent_expansion([1,2,3],user_labels=True) == S.mutate([3,2,1],inplace=False).cluster_variable(1)
            True

            sage: TP = T.principal_extension()
            sage: TP.arc_laurent_expansion([1,2,3],user_labels=True)
            (x1*x3*y1*y2*y3 + x0*x2^2 + x0*x2*y1 + x0*x2*y3 + x0*y1*y3)/(x1*x2*x3)
            sage: SP = S.principal_extension()
            sage: SP.mutate([3,2,1],inplace=False).cluster_variable(1) == TP.arc_laurent_expansion([1,2,3],user_labels=True)
            True

        A once-punctured square's triangulation with self-folded
        triangle, border edges are labeled 4,5,6,7, 2nd triangulation
        in oral paper ell-loop is labeled 3, radius is labeled 0::

            sage: T = ClusterTriangulation([(1,7,4),(1,5,2),(2,3,6),(3,0,0)], boundary_edges=[4,5,6,7])
            sage: S = ClusterSeed(T)
            sage: S.mutation_type()
            ['D', 4]
            sage: c = [item for item in T.cluster()]
            sage: r=c[0]
            sage: ell=c[3]*r
            sage: T.arc_laurent_expansion([c[1],c[2],ell,(r,'counterclockwise'),ell], user_labels=False)
            (x2^3 + x0*x1*x3 + 3*x2^2 + 3*x2 + 1)/(x0*x1*x2*x3)
            sage: gamma = T.arc_laurent_expansion([1,2,3,(0,'counterclockwise'),3], user_labels=True)
            sage: gamma == T.arc_laurent_expansion([c[1],c[2],ell,(r,'counterclockwise'),ell], user_labels=False)
            True
            sage: gamma == T.mutate([0,3,2,1],inplace=False).cluster_variable(1)
            True
            sage: gamma == T.arc_laurent_expansion([3,(0,'clockwise'),3,2,1], user_labels=True)
            True
            sage: TP = T.principal_extension()
            sage: SP = S.principal_extension()
            sage: TP.arc_laurent_expansion([1,2,3,(0,'counterclockwise'),3], user_labels=True) == \
            ....: SP.mutate([0,3,2,1],inplace=False).cluster_variable(1)
            True

        An 8-gon triangulation from Figure 2 of [SchifflerThomas]_
        where tau_i = i and tau_13 is labeled 0::

            sage: T = ClusterTriangulation([('1','7','8'), ('2','9','10'), ('3','1','2'), ('5','4','3'), ('11','12','5'), ('4','0','6')],\
            ....: boundary_edges=['6','7','8','9','10','11','12','0'])
            sage: S = ClusterSeed(T)
            sage: c = [item for item in S.cluster()]
            sage: gamma = T.arc_laurent_expansion([c[1-1],c[3-1],c[5-1]], user_labels=False)
            sage: gamma == T.arc_laurent_expansion(['5','3','1'])
            True
            sage: S.mutate([1-1,3-1,5-1], inplace=True)
            sage: S.cluster_variable(5-1) == gamma
            True
            sage: TP = T.principal_extension()
            sage: TP.mutate(['1','3','5'],user_labels=True,inplace=False).cluster_variable(4) == TP.arc_laurent_expansion(['5','3','1'],user_labels=True)
            True
            sage: TP.arc_laurent_expansion(['1','3','5']) == TP.arc_laurent_expansion(['5','3','1'])
            True

        Affine A(2,2) triangulation from Figure 3 of
        [SchifflerThomas]_ where tau_i = i and tau_8 is labeled 0::

            sage: T = ClusterTriangulation([(7,4,3),(4,1,5),(3,6,2),(2,1,0)], boundary_edges=[5,6,7,0])
            sage: S = ClusterSeed(T)
            sage: c = [item for item in T.cluster()]
            sage: gamma = T.arc_laurent_expansion([c[1-1],c[2-1],c[3-1],c[4-1],c[1-1]], user_labels=False)
            sage: gamma == T.arc_laurent_expansion([1,4,3,2,1])
            True
            sage: S.mutate([1-1,3-1,4-1,2-1,3-1], inplace=True)
            sage: S.cluster_variable(2) == gamma
            True
            sage: S.mutation_type()
            ['A', [2, 2], 1]
            sage: TP = T.principal_extension()
            sage: TP.mutate([1-1,3-1,4-1,2-1,3-1], user_labels=False, inplace=False).cluster_variable(2) == TP.arc_laurent_expansion([1,4,3,2,1])
            True

            sage: T = ClusterTriangulation([(0,2,1),(0,4,3),(1,6,5)])
            sage: S = ClusterSeed(T)
            sage: S1 = S.mutate(0, inplace=False)
            sage: S1.cluster_variable(0) == T.arc_laurent_expansion ([S.cluster_variable(0)], user_labels=False)
            True
            sage: Tp = T.principal_extension()
            sage: Sp = S.principal_extension()
            sage: Sp.mutate(0, inplace=False).cluster_variable(0) == Tp.arc_laurent_expansion ([0], user_labels=True)
            True

            sage: once_punctured_torus = ClusterTriangulation([(0,1,2),(2,0,1)])
            sage: S = ClusterSeed(once_punctured_torus)
            sage: c = S.cluster()
            sage: once_punctured_torus.arc_laurent_expansion([c[0],c[1],c[0],c[2],c[0],c[1],c[0]], user_labels=False) == \
            ....: S.mutate([0,1,2], inplace=False).cluster_variable(2)
            True
            sage: once_punctured_torus.principal_extension().arc_laurent_expansion([c[0],c[1],c[0],c[2],c[0],c[1],c[0]], user_labels=False) == \
            ....: S.principal_extension().mutate([0,1,2], inplace=False).cluster_variable(2)
            True

        An example that Salomon Dominguez wrote::

            sage: affineD = ClusterTriangulation([('5','6','4'),('6','7','10'),('7','1','11'),('1','5','2'),('2','3','8'),('3','4','9')], boundary_edges=['8','9','10','11'])
            sage: gamma = affineD.arc_laurent_expansion(['6','7','1','2','3'],verbose=False)
            sage: gamma == affineD.mutate(['3','2','1','7','6'],user_labels=True, inplace=False).cluster()[5]
            True

        Test bug for when a generalized arc's first cross and last cross are the same::

            sage: Annulus41 = ClusterTriangulation([(1,6,7),(1,3,2),(3,5,4),(5,0,2),(0,8,9)], boundary_edges=[6,7,8,9])
            sage: GeneralizedArc = Annulus41.arc_laurent_expansion([3,5,2])
            sage: Z = Annulus41.loop_laurent_expansion([2,3,5,2])
            sage: GeneralizedArc == Z * Annulus41.cluster_variable(1) + Annulus41.cluster_variable(0)
            True

        A markov quiver example which was a bug in a previous version::

            sage: once_punc_torus = ClusterTriangulation([(0,1,2),(0,1,2)])
            sage: x2_times_x2WithOneNotching = once_punc_torus.arc_laurent_expansion([0,1,2,0,1])
            sage: x2_times_x2WithOneNotching
            (2*x0^2*x2 + 2*x1^2*x2 + 2*x2^3)/(x0*x1)
        """
        from sage.combinat.cluster_algebra_quiver.surface import _get_weighted_edges, LaurentExpansionFromSurface
        CT = self#._cluster_triangulation

        if user_labels:
            crossed_arcs = _get_weighted_edges(crossed_arcs,
                                               CT._map_label_to_variable)
            first_triangle = _get_weighted_edges(first_triangle,
                                                 CT._map_label_to_variable)
            final_triangle = _get_weighted_edges(final_triangle,
                                                 CT._map_label_to_variable)
        if not return_labels: 
        	return LaurentExpansionFromSurface(CT, crossed_arcs, first_triangle, final_triangle, True, False, verbose, CT._boundary_edges_vars, fig_size=fig_size)
        else:
        	expansion = LaurentExpansionFromSurface(CT, crossed_arcs, first_triangle, final_triangle, True, False, verbose, CT._boundary_edges_vars, fig_size = fig_size)
        	# reverses the label and variable dictionary so we can use it for substitution
        	labelDict = {v:k for k,v in self._map_label_to_variable.items()}
        	expansion = str(expansion)
        	for key in labelDict.keys():
        		expansion = expansion.replace(str(key),labelDict.get(key))
        	return expansion
        	
        	#expansion = LaurentExpansionFromSurface(CT, crossed_arcs, first_triangle, final_triangle, True, False, verbose, CT._boundary_edges_vars, fig_size=fig_size)
        	#print expansion
		        # reverses the label -> variable dictionary so we can use it for substitution
		        #labelDict = {v:k for k,v in T._map_label_to_variable.items()}
		        #expansion = str(expansion)
		        #for key in labelDict.keys():
                        #    expansion = expansion.replace(str(key),labelDict.get(key))
                        #return expansion

    def loop_laurent_expansion(self, crossed_arcs, first_triangle=None,
                               final_triangle=None, verbose=False,
                               fig_size=1, user_labels=True):
        """Return the Laurent expansion of a loop (living in the
        surface's interior) in the variables of
        ``self.cluster_triangulation()._cluster``.

        See algorithm in [MSW_Bases]_ sections 3.1-3.2::

            #. Pick an orientation of the closed loop gamma and an ideal triangle Delta0
            crossed by gamma. Delta0 should be chosen so that it has 3 distinct edges.
            #. Let ``tau1`` be the second edge of Delta0 that is crossed by gamma
            (in Figure 9 of [MSW_Bases]_, this edge is labeled ``c``).
            #. Let input ``crossed_arcs`` be the list of arcs that are crossed by gamma in order,
            where ``tau1`` is counted twice, so that crossed_arcs[0]=crossed_arcs[-1]=``tau1``

        .. SEEALSO::

            :meth:`arc_laurent_expansion`,
            :meth:`ClusterTriangulation.draw_lifted_loop`

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
          crossed by loop

        - ``final_triangle`` -- (default:``None``) the last triangle
          crossed by loop

        - ``verbose`` -- (default:``False``) display the image of the
          perfect matchings of the band graph if ``True``

        - ``fig_size`` -- (default:1) image size

        - ``user_labels`` -- (default:``True``) whether or not
          ``crossed_arcs`` is a list of labels

        EXAMPLES:

        Figure 6 of [MW_MatrixFormulae]_ where tau_4, tau_1, tau_2,
        tau_3 = ``0``,``1``,``2``,``3`` and
        b1,b2,b3,b4=``b4``,``b5``,``b6``,``b7``.  We pick tau_1 to
        be 1, and go clockwise, so that crossed_arcs = [1,2,3,0,1]::

            sage: T = ClusterTriangulation([(1,2,'b4'),(1,0,'b5'),(0,3,'b6'),(2,3,'b7')], boundary_edges=['b4','b5','b6','b7'])
            sage: c = [item for item in T.cluster()]
            sage: T.loop_laurent_expansion([c[1],c[2],c[3],c[0],c[1]], user_labels=False)
            (x0*x1^2*x2 + x0*x2*x3^2 + x1^2 + 2*x1*x3 + x3^2)/(x0*x1*x2*x3)
            sage: T.loop_laurent_expansion([1,2,3,0,1], user_labels=True) == T.loop_laurent_expansion([c[1],c[2],c[3],c[0],c[1]], user_labels=False)
            True

        Example 3.6 from [DupontThomas]_::

            sage: T = ClusterTriangulation([(0,1,2),(0,1,3)], boundary_edges=[2,3])
            sage: c = [item for item in T.cluster()]
            sage: T.loop_laurent_expansion([1,0,1,0,1],first_triangle=(0,1,2)) == T.loop_laurent_expansion([0,1,0,1,0],first_triangle=(0,1,2))
            True
            sage: T.loop_laurent_expansion([1,0,1,0,1],first_triangle=(0,1,2)) == T.loop_laurent_expansion([0,1,0,1,0],first_triangle=(0,1,3))
            True

            sage: crossed_arcs = [c[0], c[1], c[0]] # loop z_1 with no self-intersection
            sage: T.loop_laurent_expansion( crossed_arcs, first_triangle=(c[0],c[1], T._get_map_label_to_variable(2)), user_labels=False)
            (x0^2 + x1^2 + 1)/(x0*x1)
            sage: crossed_arcs = [c[0], c[1], c[0], c[1], c[0]] # loop z_2 with 1 self-intersection
            sage: T.loop_laurent_expansion(crossed_arcs, first_triangle = (c[0],c[1], T._get_map_label_to_variable(2)), user_labels=False)
            (x0^4 + x1^4 + 2*x0^2 + 2*x1^2 + 1)/(x0^2*x1^2)

        Once-punctured square with 2 radii and boundary edges labeled
        4,5,6,7 where the loop is contractible to the puncture::

            sage: once_punctured_square = ClusterTriangulation([(1,7,4),(1,5,2),(6,0,3),(2,3,0),(0,3,6),[7,4,1]], boundary_edges=[4,5,6,7])
            sage: once_punctured_square.loop_laurent_expansion(crossed_arcs=[0,3,0], first_triangle = [0,3,6], user_labels=True)
            2
        """
        from sage.combinat.cluster_algebra_quiver.surface import _get_weighted_edges, LaurentExpansionFromSurface
        CT = self#._cluster_triangulation

        if user_labels:
            crossed_arcs = _get_weighted_edges(crossed_arcs,
                                               CT._map_label_to_variable)
            first_triangle = _get_weighted_edges(first_triangle,
                                                 CT._map_label_to_variable)
            final_triangle = _get_weighted_edges(final_triangle,
                                                 CT._map_label_to_variable)
        return LaurentExpansionFromSurface(CT, crossed_arcs, first_triangle, final_triangle, False, True, verbose, CT._boundary_edges_vars, fig_size=fig_size)

    def principal_extension(self,ignore_coefficients=False):
        r"""
        Returns the principal extension of self, yielding a 2n-by-n matrix.  Raises an error if the input seed has a non-square exchange matrix,
        unless 'ignore_coefficients=True' is set.  In this case, the method instead adds n frozen variables to any previously frozen variables.
        I.e., the seed obtained by adding a frozen variable to every exchangeable variable of ``self``.

        Note that this function will not do anything if called more than once, in contrast to :meth:`ClusterSeed.principal_extension`

        EXAMPLES::

            sage: T = ClusterTriangulation([(1,2,'b4'),(1,0,'b5'),(0,3,'b6'),(2,3,'b7')], boundary_edges=['b4','b5','b6','b7']); T
            A seed for a cluster algebra associated with an ideal triangulation of rank 4 with 4 boundary edges

            sage: TP = T.principal_extension(); TP
            A seed for a cluster algebra associated with an ideal triangulation of rank 4 with 4 boundary edges with principal coefficients

            sage: TP.b_matrix()
            [ 0  1  0 -1]
            [-1  0 -1  0]
            [ 0  1  0 -1]
            [ 1  0  1  0]
            [ 1  0  0  0]
            [ 0  1  0  0]
            [ 0  0  1  0]
            [ 0  0  0  1]

            sage: TP.mutate(0)
            sage: TP.cluster()
            [(x3*y0 + x1)/x0, x1, x2, x3]
            sage: TP._cluster
            [(x3*y0 + x1)/x0, x1, x2, x3]
            sage: SP = ClusterSeed(TP); SP
            A seed for a cluster algebra associated with an ideal triangulation of rank 4 with 4 boundary edges with principal coefficients of type ['A', [2, 2], 1] with principal coefficients

            sage: TP._cluster == SP._cluster
            True

            sage: TP.principal_extension()
            Traceback (most recent call last):
            ...
            ValueError: The b-matrix is not square. Use ignore_coefficients to ignore this.

            sage: T2 = TP.principal_extension(ignore_coefficients=True)
            sage: TP.b_matrix()
            [ 0 -1  0  1]
            [ 1  0 -1 -1]
            [ 0  1  0 -1]
            [-1  1  1  0]
            [-1  1  0  0]
            [ 0  1  0  0]
            [ 0  0  1  0]
            [ 0  0  0  1]
            """
        from sage.matrix.all import identity_matrix
        from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        if not ignore_coefficients and self._m != 0:
            raise ValueError("The b-matrix is not square. Use ignore_coefficients to ignore this.")
        #M = self._M.stack(identity_matrix(self._n))
        #Q = ClusterQuiver(M)
        #is_principal = (self._m == 0)
        is_principal = True
        #seed = ClusterTriangulation( self, is_principal=is_principal )
        seed = ClusterTriangulation( self._triangles, is_principal=is_principal, boundary_edges=self._boundary_edges ) # this will give initial cluster (x0,x1,..,xn) even if self has been mutated
        seed._mutation_type = self._mutation_type
        return seed
