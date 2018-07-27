r"""
Surface

This file contains helper functions for producing an initial surface
ideal triangulation for the :class:`ClusterTriangulation` class and
for computing the Laurent expansion for cluster algebra elements not
belonging to the initial ideal triangulation.

See :meth:`ClusterSeed.arc_laurent_expansion`, :meth:`ClusterSeed.loop_laurent_expansion`

See [MSW_Positivity]_, [MSW_Bases]_, [MW_MatrixFormulae]_
"""

#########################################################################################
############ begins: CREATING CLUSTER ALGEBRA FROM INITIAL TRIANGULATION INPUT ##########
#########################################################################################


RIGHT = 'RIGHT'
ABOVE = 'ABOVE'


def are_triangles_equal(triangleA, triangleB):
    """
    If triangles are equal (including orientation), return
    ``True``. Otherwise return ``False``.

    INPUT:

    - ``triangleA`` -- a 3-tuple (or list of length 3) of labels of a
      triangle (with at least 2 unique labels)

    - ``triangleB`` -- a 3-tuple (or list of length 3) of labels of a
      triangle (with at least 2 unique labels)

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import are_triangles_equal
        sage: are_triangles_equal((1,2,3),[2,3,1])
        True
        sage: are_triangles_equal((1,2,3),(1,3,2))
        False
        sage: are_triangles_equal([1,1,2],[1,2,1])
        True
    """
    if (triangleA[0], triangleA[1], triangleA[2]) == (triangleB[0], triangleB[1], triangleB[2]):
        return True
    elif (triangleA[0], triangleA[1], triangleA[2]) == (triangleB[1], triangleB[2], triangleB[0]):
        return True
    elif (triangleA[0], triangleA[1], triangleA[2]) == (triangleB[2], triangleB[0], triangleB[1]):
        return True
    else:
        return False


def is_selffolded(t):
    """
    If triangle t only two distinct entries, return the triangle t in
    the form (a,a,b).

    Otherwise, return ``False``.

    INPUT:

    - ``t`` -- a 3-tuple (or list of length 3) of labels of a triangle
      (with at least 2 unique labels)

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import is_selffolded
        sage: is_selffolded(('ell','r','r'))
        ('r', 'r', 'ell')
        sage: is_selffolded([1,1,2])
        (1, 1, 2)
        sage: is_selffolded([1,2,1])
        (1, 1, 2)
        sage: is_selffolded((1,2,3))
        False
    """
    if t[0]==t[1] and t[0]!=t[2]:
        ell = t[2]
        radius = t[0]
        return (radius, radius, ell)
    elif t[0]==t[2] and t[0]!=t[1]:
        ell = t[1]
        radius = t[0]
        return (radius, radius, ell)
    elif t[1]==t[2] and t[1]!=t[0]:
        ell = t[0]
        radius = t[1]
        return (radius, radius, ell)
    else:
        return False

def remove_duplicate_triangles(data,boundary_edges=None):
    """
    In case user accidentally inputs a duplicate triangle, we remove
    the duplicate.

    The only exception: The once-punctured torus' triangulation has
    two triangles with identical labels.

    This function also checks for triangulations that are not allowed.

    See :class:`ClusterQuiver`

    INPUT:

    - ``data`` -- a list of 3-tuples (or lists of length 3) of triangles

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import remove_duplicate_triangles
        sage: remove_duplicate_triangles([(1,2,3),(1,2,4),[3,1,2],[2,3,1]])
        [(1, 2, 3), (1, 2, 4)]
        sage: once_punctured_torus = [(1,2,3),(3,1,2)]
        sage: remove_duplicate_triangles(once_punctured_torus)
        [(1, 2, 3), (3, 1, 2)]
    """
    list_triangles = []

    # Triangle/once-punctured monogon and thrice-punctured sphere are not allowed
    if len(data) == 1 or (len(data)==2 and are_triangles_equal(data[0],(data[1][::-1]))):
        raise ValueError('The following surfaces are not allowed: a sphere with 1, 2, or 3 punctures; a monogon with zero or 1 puncture; a bigon or triangle without punctures.')

    if len(data) == 2 and are_triangles_equal(data[0],data[1]):
        if boundary_edges:
            raise ValueError('Two triangles sharing the same edges form a triangulation of a once-punctured torus, which has no boundary')
        else:
            return data

    for triangle in data:
        is_duplicate_triangle = False
        for t in list_triangles:
            if are_triangles_equal(t, triangle):
                is_duplicate_triangle = True
                break
        if not is_duplicate_triangle:
            list_triangles.append(triangle)
    return list_triangles


def _triangulation_to_arrows(list_triangles):
    """
    Return the list of directed edges [a,b] corresponding to the list
    ``list_triangles`` of ideal triangles.

    .. NOTE::

        The order of ``list_triangles`` does not matter. The output
        may contain 2-cycles.

    For each triangle t in list_triangles:

    #. if t = [a, b, c] has distinct edges, add cycles a -> b -> c -> a
    #. if t = [r, r, ell] is a self-folded triangle and there is an
       edge between ell and b, add the same directed edge between r
       and b.
    #. if t = [a, b, c] has distinct edges and a and b are both nooses (of
       two self-folded triangles [ra,ra,a] and [rb,rb,b]), add an extra
       arrow between ra and rb.

    See :func:`~sage.combinat.cluster_algebra_quiver.surface._surface_edge_list_to_matrix` and :class:`ClusterQuiver`

    INPUT:

    - ``list_triangles`` -- a list of triangles, i.e. a 3-tuple (or a
      list of length 3) with at least 2 unique labels

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _triangulation_to_arrows, _surface_edge_list_to_matrix, _get_user_arc_labels
        sage: annulus2_2 = [('bd1','tau1','tau2'),('tau2','tau3','bd4'),('tau1','tau4','bd2'),('tau3','bd3','tau4')]
        sage: _triangulation_to_arrows(annulus2_2)
        [['bd1', 'tau1', None],
         ['tau1', 'tau2', None],
         ['tau2', 'bd1', None],
         ['tau2', 'tau3', None],
         ['tau3', 'bd4', None],
         ['bd4', 'tau2', None],
         ['tau1', 'tau4', None],
         ['tau4', 'bd2', None],
         ['bd2', 'tau1', None],
         ['tau3', 'bd3', None],
         ['bd3', 'tau4', None],
         ['tau4', 'tau3', None]]

        sage: annulus1_1 = [[2, 1, 0], [3, 2, 1]]
        sage: _surface_edge_list_to_matrix(_triangulation_to_arrows(annulus1_1),_get_user_arc_labels(annulus1_1),[],4)
        [ 0  1 -1  0]
        [-1  0  2 -1]
        [ 1 -2  0  1]
        [ 0  1 -1  0]

        sage: T = [[1,4,2],[3,4,3],[2,0,1]] # twice-punctured monogon with 3 (non-ordinary) ideal triangles (affine D)
        sage: B = _surface_edge_list_to_matrix(_triangulation_to_arrows(T),_get_user_arc_labels(T),[],5)
        sage: Q = ClusterQuiver(B)
        sage: Q.b_matrix()
        [ 0 -1  1  0  0]
        [ 1  0  0 -1 -1]
        [-1  0  0  1  1]
        [ 0  1 -1  0  0]
        [ 0  1 -1  0  0]

        sage: Tmu2 = [(1,1,2),(3,4,3),(2,4,0)] # 2 self-folded triangles and 1 triangle with one vertex (affine D)
        sage: Bmu2 = _surface_edge_list_to_matrix(_triangulation_to_arrows(Tmu2),_get_user_arc_labels(Tmu2),[],5)
        sage: Bmu2
        [ 0 -1 -1  1  1]
        [ 1  0  0 -1 -1]
        [ 1  0  0 -1 -1]
        [-1  1  1  0  0]
        [-1  1  1  0  0]

        sage: Bmu2.mutate(2)
        sage: Bmu2 == B
        True

        sage: Qmu2 = ClusterQuiver(Bmu2)
        sage: Qmu2.mutation_type()
        'undetermined finite mutation type'
        sage: Qmu2.is_mutation_finite()
        True

        sage: arrows = _triangulation_to_arrows([[4, 5, 1], [4, 3, 2], [3, 7, 2], [2, 1, 6], [1, 4, 5]])
        sage: arrows
        [[4, 5, None],
        [5, 1, None],
        [1, 4, None],
        [4, 3, None],
        [3, 2, None],
        [2, 4, None],
        [3, 7, None],
        [7, 2, None],
        [2, 3, None],
        [2, 1, None],
        [1, 6, None],
        [6, 2, None],
        [1, 4, None],
        [4, 5, None],
        [5, 1, None]]
    """
    digraph_edges = []
    selffolded_triangles = []
    nooses = []

    for t in list_triangles:
        selffolded = is_selffolded(t)
        if t[0] != t[1] and t[1] != t[2] and t[0] != t[2]:
            # if t has three distinct edges, add 3 edges for those
            digraph_edges.extend([[t[0],t[1],None], [t[1],t[2],None], [t[2],t[0],None]])
        elif selffolded != False:
            radius = selffolded[0]
            ell = selffolded[2]
            selffolded_triangles.append([radius,radius, ell])
            if _is_arc_mutable(list_triangles, ell)[0]:
                nooses.append(ell)
            else:
                raise ValueError('A noose of a self-folded triangle must be a side of another triangle.')
        else:
            raise ValueError('An ideal triangle has to have 3 distinct edges or 2 distinct edges')
            break

    selffolded_edges = []
    radius_to_radius_edges = []
    for t in selffolded_triangles:
        radius = t[0]
        ell = t[2]
        flagFoundSharedNooseA = False
        flagFoundSharedNooseB = False
        for e in digraph_edges:
            if e[0]==ell:
                selffolded_edges.append([radius, e[1], None])
                if e[1] in nooses:
                    radius_e = _get_radius(e[1],selffolded_triangles)
                    if [radius, radius_e, None] not in radius_to_radius_edges:
                        radius_to_radius_edges.append([radius, radius_e, None])
                flagFoundSharedNooseA = True
            elif e[1]==ell:
                selffolded_edges.append([e[0],radius, None])
                if e[0] in nooses:
                    radius_e = _get_radius(e[0],selffolded_triangles)
                    if [radius_e, radius, None] not in radius_to_radius_edges:
                        radius_to_radius_edges.append([radius_e, radius, None])
                flagFoundSharedNooseB = True
            if flagFoundSharedNooseA and flagFoundSharedNooseB:
                break

    return digraph_edges + selffolded_edges + radius_to_radius_edges

def _get_radius(noose,selffolded_triangles):
    """
    Return the radius of the self-folded triangle with input noose from the list of all self-folded triangles.

    See :func:`~sage.combinat.cluster_algebra_quiver.surface._triangulation_to_arrows`

    INPUT:

    - ``noose`` -- a string or number which labels a noose arc.
    - ``self-folded triangles`` -- a list of 3-tuples in the form of [[radiusA, radiusA, nooseA], [radiusB, radiusB, nooseB], ...]

    EXAMPLES:

    A triangulation with 2 self-folded triangles and 1 triangle with one vertex (affine D)::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _get_radius
        sage: _get_radius(4,[(1,1,2),(3,3,4)])
        3
    """
    for t in selffolded_triangles:
        if noose == t[2]:
            return t[0]

    raise ValueError('There is no self-folded triangles with noose ', noose)

def _surface_edge_list_to_matrix(arrows, arcs_and_boundary_edges, boundary_edges, arcs_count):
    """
    Returns the matrix obtained from the edge list of a quiver (possibly with loops, two-cycles, multiply-listed edges).

    .. NOTE:

    This function does not ignore duplicates. Compare this to
    :func:`sage.combinat.cluster_algebra_quiver.quiver_mutation_type._edge_list_to_matrix`
    which ignores duplicates::

        sage: from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import _edge_list_to_matrix
        sage: _edge_list_to_matrix([(0,1,(1,-1)),(1,2,(1,-1)),(2,0,(1,-1)),(2,0,(1,-1))],3,0)
        [ 0  1 -1]
        [-1  0  1]
        [ 1 -1  0]

        sage: from sage.combinat.cluster_algebra_quiver.surface import _surface_edge_list_to_matrix
        sage: _surface_edge_list_to_matrix([(0,1,None),(1,2,None),(2,0,None),(2,0,None)],[0,1,2],[],3)
        [ 0 -1  2]
        [ 1  0 -1]
        [-2  1  0]

        sage: _edge_list_to_matrix([(0,1,(1,-1)),(1,0,(1,-1))],2,0)
        [ 0 -1]
        [ 1  0]

        sage: _surface_edge_list_to_matrix([(0,1,None),(1,0,None),(0,2,None),\
        (2,1,None),(1,3,None),(3,0,None)],[0,1,2,3],[2,3],2)
        [0 0]
        [0 0]

    INPUT:

    - ``arrows``: the list of edges from the quiver associated to a triangulation.
    - ``arcs_and_boundary_edges``: the list of user-given labels of arcs and boundary edges
    - ``boundary_edges``: the list of user-given boundary edges
    - ``arcs_count``: the number of surface edges (that are not boundary edges)

    OUTPUT:

    - An `arcs_count \times arcs_count` skew-symmetric matrix corresponding to the arrows-list (ignoring the boundary edges).

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _surface_edge_list_to_matrix, _triangulation_to_arrows
        sage: annulus2_2 = [('bd1','tau1','tau2'),('tau2','tau3','bd4'),('tau1','tau4','bd2'),('tau3','bd3','tau4')]
        sage: _surface_edge_list_to_matrix(_triangulation_to_arrows(annulus2_2),['bd1', 'bd2', 'bd3', 'bd4', 'tau1', 'tau2', 'tau3', 'tau4'],['bd1', 'bd2', 'bd3', 'bd4'],4)
        [ 0 -1  0 -1]
        [ 1  0 -1  0]
        [ 0  1  0  1]
        [ 1  0 -1  0]

        sage: _surface_edge_list_to_matrix([[0,1,None],[1,2,None],[2,0,None],[2,3,None],[0,3,None],[3,0,None]],[0,1,2,3],[],4)
        [ 0 -1  1  0]
        [ 1  0 -1  0]
        [-1  1  0 -1]
        [ 0  0  1  0]
    """
    from sage.rings.all import ZZ
    from sage.matrix.all import matrix

    arcs = []
    if len(boundary_edges) == 0:
        arcs = arcs_and_boundary_edges
    else:
        for pos in range(0,len(arcs_and_boundary_edges)):
            tau = arcs_and_boundary_edges[pos]
            if tau not in boundary_edges:
                arcs.append(tau)

    #dic = []
    if len(arcs) == arcs_count:
        dic = dict((arcs[pos],pos) for pos in range(0,arcs_count))
        #for pos in range(0,arcs_count):
        #    dic.append((arcs[pos], pos))

    M = matrix(ZZ,arcs_count,arcs_count,sparse=True)

    for user_edge in arrows:
        if user_edge[0] in boundary_edges or user_edge[1] in boundary_edges:
            continue
        if user_edge[2] is None: # todo: arrows coming from surfaces will not have weight, so the elif should be removed
            #edge = (_get_weighted_edge(user_edge[0],dic), _get_weighted_edge(user_edge[1],dic), (1,-1))
            edge = (dic[user_edge[0]],dic[user_edge[1]],(1,-1))
        elif user_edge[2] in ZZ:
            #edge = (_get_weighted_edge(user_edge[0],dic), _get_weighted_edge(user_edge[2],dic), (user_edge[2],-user_edge[2]))
            edge = (dic[user_edge[0]],dic[user_edge[1]],(2,-2))
        #print 'edge: ', edge
        v1,v2,(a,b) = edge

        if v1 < arcs_count:
            M[v2,v1] += b
        if v2 < arcs_count:
            M[v1,v2] += a
    return -1*M

def _get_user_arc_labels(T):
    """
    Return the sorted list of all entries from a list of labels in a list of 3-tuples

    INPUT:

    - ``T``: a list of 3-tuples

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _get_user_arc_labels
        sage: annulus2_2 = [('bd1','tau1','tau2'),('tau2','tau3','bd4'),('tau1','tau4','bd2'),('tau3','bd3','tau4')]
        sage: _get_user_arc_labels(annulus2_2)
        ['bd1', 'bd2', 'bd3', 'bd4', 'tau1', 'tau2', 'tau3', 'tau4']
    """
    T_user_labels = list(set([edge for triangle in T for edge in triangle]))
    T_user_labels.sort()
    return T_user_labels

def produce_dict_label_to_variable(T, cluster_xs, boundary_edges, boundary_edges_vars):
    """
    Return a dictionary of tuples ``{a:b, ...}`` where ``a`` is a user-given label from input ``T``
    and ``b`` is the variable x_i or b_i corresponding to ``a1``.

    See :class:`ClusterTriangulation`

    INPUT:

    - ``T``: triangulation, a list of 3-tuples, see ``ClusterTriangulation.triangulation``
    - ``cluster``: see :meth:`~ClusterTriangulation.cluster`
    - ``boundary_edges``: a list of user-given boundary edges, possibly an empty list, see :meth:`~ClusterTriangulation.boundary_edges`
    - ``boundary_edges_vars``: see :meth:`~ClusterTriangulation.boundary_edges_vars`

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import produce_dict_label_to_variable
        sage: Triangles = [(1, 4, 7), (1, 2, 5), (6, 3, 0), (2, 0, 3), (0, 6, 3), [7, 1, 4]]
        sage: T = ClusterTriangulation(Triangles)
        sage: produce_dict_label_to_variable(T._triangles, T.cluster(), T._boundary_edges, T._boundary_edges_vars)
        {0: x0, 1: x1, 2: x2, 3: x3, 4: x4, 5: x5, 6: x6, 7: x7}

        sage: twice_punctured_bigon = [(1,1,2),(3,4,3),(2,4,0),(0,6,7)]
        sage: T = ClusterTriangulation(twice_punctured_bigon)
        sage: produce_dict_label_to_variable(T._triangles, T.cluster(), T._boundary_edges, T._boundary_edges_vars)
        {0: x0, 1: x1, 2: x1*x2, 3: x3, 4: x3*x4, 6: x5, 7: x6}

        sage: twice_punctured_bigon = [('e','d','a'), ('a','r','b'), ('r','d','g'), ('g','n','b')]
        sage: T = ClusterTriangulation(twice_punctured_bigon, boundary_edges=['e','n'])
        sage: produce_dict_label_to_variable(T._triangles, T._cluster, T._boundary_edges, T._boundary_edges_vars)
        {'a': x0, 'b': x1, 'd': x2, 'e': b5, 'g': x3, 'n': b6, 'r': x4}
    """
    dic = []
    list_radius = []
    list_ell = []

    T_user_labels = _get_user_arc_labels(T)
    arcs = []
    boundary_edges.sort()

    #print 'T_user_labels: ',T_user_labels

    for l in T_user_labels:
        if l not in boundary_edges:
            arcs.append(l)

    #print 'arcs: ', arcs
    #print 'boundary_edges: ', boundary_edges
    #print 'cluster_xs: ', cluster_xs
    #print 'boundary_edges_vars: ', boundary_edges_vars
    #if len(arcs) + len(boundary_edges)==len(cluster) + len(boundary_edges_vars):
    dic_x = dict((arcs[pos],cluster_xs[pos]) for pos in range(0,len(cluster_xs)))
    dic_b = dict((boundary_edges[pos],boundary_edges_vars[pos]) for pos in range(0,len(boundary_edges_vars)))
    dic = dict(list(dic_x.items()) + list(dic_b.items()))

        #for pos in range(0,len(cluster)):
        #    dic.append((arcs[pos], cluster[pos]))
        #for pos in range(0,len(boundary_edges_vars)):
        #    dic.append((boundary_edges[pos], boundary_edges_vars[pos]))

    #print 'dic: ', dic

    for triangle in T: # Keep track of any radius and ell-loop of a self-folded triangle
        selffolded = is_selffolded(triangle)
        if type(selffolded) in [tuple,list]:
            radius_label = selffolded[0]
            ell_label = selffolded[2]
            radius = _get_weighted_edge(radius_label, dic)
            ell = _get_weighted_edge(ell_label, dic)
            list_radius.append(radius)
            list_ell.append(ell)

    #for pos in range(0,len(dic.keys())): # Assign a noose to the product of r * r\notch
    for user_label in dic:
        cluster_var = dic.get(user_label) #dic[pos][1]
        if list_ell.count(cluster_var):
            ind = list_ell.index(cluster_var)
            dic[user_label] =  list_ell[ind] * list_radius[ind]

    #print 'I am in _get_map_label_to_variable, dic: ', dic
    return dic

def produce_dict_variable_to_label(td):
    """
    Return a dictionary of tuples ``{a:b, ...}`` where is a variable x_i or b_i and
    ``a`` is the user-given label corresponding to ``b``.

    See :class:`ClusterTriangulation` and :func:`sage.combinat.cluster_algebra_quiver.quiver_mutation_type.produce_dict_label_to_label`

    INPUT:

    - ``td``: a dictionary from edge labels to variables x_i or b_i

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import produce_dict_label_to_variable, produce_dict_variable_to_label
        sage: Triangles = [(1, 4, 7), (1, 2, 5), (6, 3, 0), (2, 0, 3), (0, 6, 3), [7, 1, 4]]
        sage: T = ClusterTriangulation(Triangles)
        sage: produce_dict_variable_to_label(produce_dict_label_to_variable(T._triangles, T._cluster, T._boundary_edges, T._boundary_edges_vars))
        {x7: 7, x6: 6, x5: 5, x4: 4, x3: 3, x2: 2, x1: 1, x0: 0}

        sage: twice_punctured_bigon = [(1,1,2),(3,4,3),(2,4,0),(0,6,7)]
        sage: T = ClusterTriangulation(twice_punctured_bigon)
        sage: produce_dict_variable_to_label(produce_dict_label_to_variable(T._triangles, T._cluster, T._boundary_edges, T._boundary_edges_vars))
        {x6: 7, x5: 6, x3: 3, x1: 1, x0: 0, x3*x4: 4, x1*x2: 2}

        sage: twice_punctured_bigon = [('e','d','a'), ('a','r','b'), ('r','d','g'), ('g','n','b')]
        sage: T = ClusterTriangulation(twice_punctured_bigon, boundary_edges=['e','n'])
        sage: produce_dict_variable_to_label(produce_dict_label_to_variable(T._triangles, T._cluster, T._boundary_edges, T._boundary_edges_vars))
        {b6: 'n', b5: 'e', x4: 'r', x3: 'g', x2: 'd', x1: 'b', x0: 'a'}
    """
    return {v: k for k, v in td.items()}

def _get_edge_user_label(edge_var, triangulation_dictionary_variable_to_label):
    """
    Return tuple[0] if the input ``edge`` if tuple[1] is in the the input
    ``triangulation dictionary``

    INPUT:

    - ``edge_var`` -- a variable x_i or b_i returned by :
        meth:`ClusterTriangulation.cluster` or :meth:`ClusterTriangulation.boundary_edges_vars`
    - ``triangulation_dictionary_variable_to_label`` -- see :meth:`ClusterTriangulation.map_variable_to_label`

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _get_edge_user_label
        sage: S = ClusterSeed(['A',4])
        sage: _get_edge_user_label(S.x(3), {S.x(0):1,S.x(1):6,S.x(2):7,S.x(3):-9})
        -9

        sage: twice_punctured_monogon = [(4,5,5),(2,3,3),(1,4,2)]
        sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=[1])
        sage: T.map_variable_to_label()
        {b4: 1, x3: 5, x1: 3, x2*x3: 4, x0*x1: 2}
        sage: _get_edge_user_label(T._cluster[2]*T._cluster[3], T._map_variable_to_label)
        4

        sage: TT = ClusterTriangulation([('j1','j1','j2'),('j3','j4','j3'),('j2','j4','j0')])
        sage: _get_edge_user_label(TT._cluster[1]*TT._cluster[2], TT._map_variable_to_label)
        'j2'
    """
    if triangulation_dictionary_variable_to_label.has_key(edge_var):
        return triangulation_dictionary_variable_to_label[edge_var]
    #for td in triangulation_dictionary:
    #    if td[1] == edge_var:
    #        return td[0]

def _get_weighted_edge(arc, triangulation_dictionary):
    """
    Return the variable corresponding to the input arc label.

    See :func:`~sage.combinat.cluster_algebra_quiver.surface._get_weighted_edges`

    INPUT:

    - ``arc`` -- user-given label for an arc or boundary edge, possibly a tuple (label,'counterclockwise') or (label,'clockwise')
    - ``triangulation_dictionary`` -- see :meth:`ClusterTriangulation.map_label_to_variable`

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _get_weighted_edge
        sage: twice_punctured_monogon = [(4,5,5),(2,3,3),(1,4,2)]
        sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=[1])
        sage: T.map_label_to_variable()
        {1: b4, 2: x0*x1, 3: x1, 4: x2*x3, 5: x3}
        sage: _get_weighted_edge((5, 'clockwise'), T._map_label_to_variable)
        (x3, 'clockwise')
        sage: _get_weighted_edge(2, T._map_label_to_variable)
        x0*x1
        sage: _get_weighted_edge(1, T._map_label_to_variable)
        b4
    """
    if arc is None:
        return arc

    dir = None

    if isinstance(arc,tuple):
        arc_label, dir = arc
    else:
        arc_label = arc

    if triangulation_dictionary.has_key(arc_label):
        return (triangulation_dictionary[arc_label], dir) if dir else triangulation_dictionary[arc_label]

    #for td in triangulation_dictionary:
    #    if td[0] == arc_label:
    #        return (td[1], dir) if dir else td[1]

def _get_weighted_edges(edges, triangulation_dictionary):
    """
    Return the variables corresponding to the list of input arc label/s

    See :func:`~sage.combinat.cluster_algebra_quiver.surface._get_weighted_edge`

    INPUT:

    - ``edges`` -- list of user-given labels for arcs/boundary edges
    - ``triangulation_dictionary`` -- see :meth:`ClusterTriangulation.map_label_to_variable`

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _get_weighted_edges
        sage: twice_punctured_monogon = [(4,5,5),(2,3,3),(1,4,2)]
        sage: T = ClusterTriangulation(twice_punctured_monogon, boundary_edges=[1])
        sage: T.map_label_to_variable()
        {1: b4, 2: x0*x1, 3: x1, 4: x2*x3, 5: x3}
        sage: _get_weighted_edges([4, (5, 'clockwise'), 4, 2], T._map_label_to_variable)
        [x2*x3, (x3, 'clockwise'), x2*x3, x0*x1]
        sage: _get_weighted_edges([2,1], T._map_label_to_variable)
        [x0*x1, b4]
    """
    if edges is None:
        return edges

    weighted_edges = []
    for e in edges:
        weighted_e = _get_weighted_edge(e, triangulation_dictionary)
        weighted_edges.append(weighted_e)
    return weighted_edges

def _get_weighted_triangulation(T, triangulation_dictionary):
    """
    Return the triangulation ``T`` given by user such that user-given labels are replaced by variables x_i and b_i.
    See :class:`ClusterTriangulation` and :meth:`ClusterTriangulation.map_label_to_variable`

    INPUT:

    - ``T`` -- list of triangles with labels given by user
    - ``triangulation_dictionary`` -- see :meth:`ClusterTriangulation.map_label_to_variable`

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _get_weighted_triangulation
        sage: T = ClusterTriangulation([(1, 4, 7), (1, 2, 5), (6, 3, 0), (2, 0, 3), (0, 6, 3), [7, 1, 4]])
        sage: T._map_label_to_variable
        {0: x0, 1: x1, 2: x2, 3: x3, 4: x4, 5: x5, 6: x6, 7: x7}
        sage: _get_weighted_triangulation(T._triangles, T._map_label_to_variable)
        [(x1, x4, x7), (x1, x2, x5), (x6, x3, x0), (x2, x0, x3)]
    """
    weighted_T = []
    for pos in range(0,len(T)):
        a = _get_weighted_edge(T[pos][0], triangulation_dictionary)
        b = _get_weighted_edge(T[pos][1], triangulation_dictionary)
        c = _get_weighted_edge(T[pos][2], triangulation_dictionary)

        false_or_r_r_ell = is_selffolded ((a,b,c))
        if isinstance(false_or_r_r_ell, tuple):
            a = (false_or_r_r_ell[0],'counterclockwise')
            b = (false_or_r_r_ell[1],'clockwise')
            c = false_or_r_r_ell[2]
        weighted_T.append((a,b,c))
    return weighted_T

def _get_user_label_triangulation(T):
    """
    Return the same input list ``T`` (of triangles with user-given edge labels)
    where each triangle (a,a,b) or (a,b,a) or (b,a,a) is replaced by ((a, 'counterclockwise'),(a, 'clockwise'),b)

    INPUT:

    - ``T`` -- list of triangles with labels given by user

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _get_user_label_triangulation
        sage: twice_punctured_monogon = [(1,1,2), (4,4,3), ('boundary',2,3)]
        sage: _get_user_label_triangulation(twice_punctured_monogon)
        [((1, 'counterclockwise'), (1, 'clockwise'), 2),
        ((4, 'counterclockwise'), (4, 'clockwise'), 3),
        ('boundary', 2, 3)]

        sage: once_punctured_square = [('a','d','c'), ('a','ll','b'), ('r','r','ll'),('b','f','e')]
        sage: _get_user_label_triangulation(once_punctured_square)
        [('a', 'd', 'c'),
        ('a', 'll', 'b'),
        (('r', 'counterclockwise'), ('r', 'clockwise'), 'll'),
        ('b', 'f', 'e')]
    """
    informative_T = []
    for pos in range(0,len(T)):
        a, b, c = T[pos]

        false_or_r_r_ell = is_selffolded ((a,b,c))
        if isinstance(false_or_r_r_ell, tuple):
            a = (false_or_r_r_ell[0],'counterclockwise')
            b = (false_or_r_r_ell[1],'clockwise')
            c = false_or_r_r_ell[2]
        informative_T.append((a,b,c))
    return informative_T

############# ENDING: CREATING CLUSTER ALGEBRA FROM INITIAL TRIANGULATION INPUT ##########
##########################################################################################

##########################################################################################
########################################### BEGINS: MUTATE A TRIANGULATION ###############

def _is_arc_mutable(T, diagonal):
    """
    Return True if ``diagonal`` is the side of two triangles
    or if ``diagonal`` is the radius of a self-folded triangle which noose
    is the side of another triangle

    INPUT:

    - ``diagonal`` -- a side of one of the triangles
    - ``T`` -- list of triangles (3-tuples of user-given labels)

    EXAMPLES::
        sage: from sage.combinat.cluster_algebra_quiver.surface import _is_arc_mutable

        sage: once_punctured_disk = [('r','r','ell'),('a','ell','b')]
        sage: _is_arc_mutable(once_punctured_disk, 'r')
        (True, [('r', 'r', 'ell'), ('a', 'ell', 'b')], True)
        sage: _is_arc_mutable(once_punctured_disk, 'ell')
        (True, [('r', 'r', 'ell'), ('a', 'ell', 'b')], False)
        sage: _is_arc_mutable(once_punctured_disk, 'b')
        (False,
        ('The ideal triangulation cannot be mutated at ',
        'b',
        '.There is only one triangle ',
        ('a', 'ell', 'b'),
        ', not a self-folded triangle, with side ',
        'b'),
        False)

        sage: once_punctured_torus = [(0,1,2),(0,1,2)]
        sage: _is_arc_mutable(once_punctured_torus, 2)
        (True, [(0, 1, 2), (0, 1, 2)], False)

        sage: bad_surface = [(0,0,1)] + once_punctured_disk
        sage: _is_arc_mutable(bad_surface, 1)
        (False,
        ('The ideal triangulation cannot be mutated at ',
        1,
        '. There is only one triangle, ',
        (0, 0, 1),
        ', with side ',
        1,
        ', but ',
        1,
        ' is the noose, not the radius'),
        False)

        sage: _is_arc_mutable(bad_surface, 0)
        (False,
        ('The ideal triangulation cannot be mutated at radius ',
        0,
        '. The noose surrounding it is not the side of any triangle with three distinct sides.'),
        True)
    """
    two_triangles = _get_triangle(T, diagonal, None)
    if len(two_triangles) == 2:
        return (True, two_triangles, False)
    elif len(two_triangles) == 1:
        selffolded_tri = two_triangles[0]
        rrl = is_selffolded(selffolded_tri)
        if not rrl:
            #raise ValueError(
            ErrorMessage = 'The ideal triangulation cannot be mutated at ', diagonal, \
            '.There is only one triangle ', selffolded_tri, \
            ', not a self-folded triangle, with side ', diagonal
            return (False, ErrorMessage, False)

        r = rrl[0]
        ell = rrl[2]

        #selffolded_rrl = (r,r,ell)
        if r != diagonal:
            ErrorMessage = 'The ideal triangulation cannot be mutated at ', diagonal, \
            '. There is only one triangle, ', selffolded_tri, \
            ', with side ', diagonal, \
            ', but ', diagonal, ' is the noose, not the radius'
            return (False, ErrorMessage, False)
                #raise ValueError('The ideal triangulation cannot be mutated at ', diagonal, \
                #'. There is only one triangle ', selffolded_tri, \
                #', and ', diagonal, ' is the noose, not the radius')

        two_triangles_with_side_ell = _get_triangle(T, ell, None)
        if len(two_triangles_with_side_ell) < 2:
            #raise ValueError(
            ErrorMessage = 'The ideal triangulation cannot be mutated at radius ', diagonal, \
            '. The noose surrounding it is not the side of any triangle with three distinct sides.'
            return (False, ErrorMessage, True)

        if len(two_triangles_with_side_ell) > 2:
            ErrorMessage = 'The ideal triangulation cannot be mutated at radius ', diagonal, \
            '. The noose surrounding it is the side of too many triangles: ', two_triangles_with_side_ell
            return (False, ErrorMessage, True)

        if len(two_triangles_with_side_ell) == 2:
            two_triangles_with_side_ell.remove(selffolded_tri)
            three_distinct_vertex_triangle = two_triangles_with_side_ell[0]
            return (True, [selffolded_tri, three_distinct_vertex_triangle], True)
    elif len(two_triangles) > 2:
        raise ValueErro('TODO: is this user error or bug? There should not be more than two triangles with side ', diagonal)

def _triangles_mutate(T,diagonal):
    """
    Return a list of triangles after we mutate at ``diagonal``, i.e.
    we remove the two triangles (before_tau_A, diagonal, after_tau_A) and
    (before_tau_B, diagonal, after_tau_B) which have ``diagonal`` as an edge from the
    list ``T``, and replace them with two new triangles
    (after_tau_A, diagonal, before_tau_B) and (after_tau_B, diagonal, before_tau_A).
    Note that the input list ``T`` will be modified.
    Note that the new arc k' is labeled k.
    Note that if diagonal is a self-folded triangle's radius, then applying the
    function twice will give a symmetric but different list of triangles,
    which may not be a desirable outcome TODO

    See :meth:`ClusterTriangulation.mutate`

    INPUT:

    - ``triangles`` -- list of triangles (3-tuples of user-given labels)
    - ``diagonal`` -- the edge that is to be mutated

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _triangles_mutate

        sage: once_punctured_disk = [('r','r','ell'),('a','ell','b')]
        sage: _triangles_mutate(once_punctured_disk, 'r')
        [('a', 'ell', 'r'), ('r', 'ell', 'b')]
        sage: _triangles_mutate(once_punctured_disk, 'r')
        [('a', 'r', 'b'), ('ell', 'r', 'ell')]
        sage: once_punctured_disk
        [('a', 'r', 'b'), ('ell', 'r', 'ell')]
        sage: _triangles_mutate(once_punctured_disk, 'a')
        Traceback (most recent call last):
        ...
        ValueError: ('The ideal triangulation cannot be mutated at ', 'a', '.There is only one triangle ', ('a', 'r', 'b'), ', not a self-folded triangle, with side ', 'a')

        sage: twice_punctured_disk = [(0,0,1),(1,3,2),(2,5,4),(4,6,3)]
        sage: _triangles_mutate(twice_punctured_disk,0)
        [(2, 5, 4), (4, 6, 3), (2, 1, 0), (0, 1, 3)]

        sage: four_punc_sphere = [(1,0,5),(3,5,4),(2,1,3),(0,2,4)]
        sage: Tmu0 = _triangles_mutate(four_punc_sphere,0)
        sage: Tmu0
        [(3, 5, 4), (2, 1, 3), (5, 0, 4), (2, 0, 1)]
        sage: four_punc_sphere
        [(3, 5, 4), (2, 1, 3), (5, 0, 4), (2, 0, 1)]
        sage: Tmu04 = _triangles_mutate(four_punc_sphere,4)
        sage: Tmu04
        [(2, 1, 3), (2, 0, 1), (3, 4, 0), (5, 4, 5)]
        sage: Tmu040 = _triangles_mutate(four_punc_sphere,0)
        sage: Tmu040
        [(2, 1, 3), (5, 4, 5), (1, 0, 4), (3, 0, 2)]
        sage: Tmu0402 = _triangles_mutate(four_punc_sphere,2)
        sage: Tmu0402
        [(5, 4, 5), (1, 0, 4), (1, 2, 0), (3, 2, 3)]
        sage: Tmu04020 = _triangles_mutate(four_punc_sphere,0)
        sage: Tmu04020
        [(5, 4, 5), (3, 2, 3), (4, 0, 2), (1, 0, 1)]
        sage: Tmu040203 = _triangles_mutate(four_punc_sphere,3)
        sage: Tmu040203
        [(5, 4, 5), (1, 0, 1), (0, 2, 3), (3, 2, 4)]
    """
    #two_triangles = _get_triangle(T, diagonal, None)

    #if len(two_triangles) == 2:
    #    triangleA, triangleB = two_triangles[0], two_triangles[1]
    is_mutable, two_triangles, is_selffolded_radius = _is_arc_mutable(T, diagonal)
    if not is_mutable:
        ErrorMessage = two_triangles
        raise ValueError(ErrorMessage)

    if is_mutable and not is_selffolded_radius:
        triangleA, triangleB = two_triangles
        T.remove(triangleA)
        T.remove(triangleB)
        triangleA = _rearrange_triangle_for_snakegraph(triangleA, diagonal, 1)
        triangleB = _rearrange_triangle_for_snakegraph(triangleB, diagonal, 1)
        before_tau_A, after_tau_A = triangleA[0], triangleA[2]
        before_tau_B, after_tau_B = triangleB[0], triangleB[2]
        T.extend([(after_tau_A, diagonal, before_tau_B),(after_tau_B, diagonal, before_tau_A)])
    elif is_mutable and is_selffolded_radius:
    #elif len(two_triangles) == 1:
    #    selffolded_tri = two_triangles[0]
    #    if is_selffolded(selffolded_tri):
    #        r,r,ell = is_selffolded(selffolded_tri)
    #        if r != diagonal:
    #            raise ValueError('The ideal triangulation cannot be mutated at ', diagonal, \
    #            '. There is only one triangle ', selffolded_tri, \
    #            ', and ', diagonal, ' is the noose, not the radius')
    #    else:
    #        raise ValueError('The ideal triangulation cannot be mutated at ', diagonal, \
    #        '.There is only one triangle ', selffolded_tri, \
    #        ', not a self-folded triangle, with side ', diagonal)
    #    two_triangles_with_side_ell = _get_triangle(T, ell, None)
    #    if len(two_triangles_with_side_ell) == 2:
    #        two_triangles_with_side_ell.remove(selffolded_tri)
    #        three_distinct_vertex_triangle = two_triangles_with_side_ell[0]
    #    else:
    #        raise ValueError('The ideal triangulation cannot be mutated at radius ', diagonal, \
    #        '. The noose surrounding it is not the side of any triangle with three distinct sides.')
        selffolded_tri, three_distinct_vertex_triangle = two_triangles
        T.remove(selffolded_tri)
        T.remove(three_distinct_vertex_triangle)

        rrl = is_selffolded(selffolded_tri)
        ell = rrl[2]#; print 'ell :', ell #TODO erase
        r = rrl[0]#; print 'r :', r #TODO erase
        three_distinct_vertex_triangle = _rearrange_triangle_for_snakegraph(three_distinct_vertex_triangle, ell, 1)
        #print 'three_distinct_vertex_triangle after rearrange: ', three_distinct_vertex_triangle #TODO Erase
        before_ell, after_ell = three_distinct_vertex_triangle[0], three_distinct_vertex_triangle[2]

        T.extend([(before_ell, ell, r), (r, ell, after_ell)])
    #else:
    #    raise ValueError(ErrorMessage) #('The search for the two triangles with edge ', diagonal, ' returns ', two_triangles)

    return T

############# ENDING: MUTATE A TRIANGULATION #############################################
##########################################################################################



##########################################################################################
########################################### BEGINS: LAURENT EXPANSION ####################

def LaurentExpansionFromSurface(CT, crossed_arcs, first_triangle=None,
                                final_triangle=None, is_arc=None,
                                is_loop=None, verbose=False,
                                boundary_edges=None, fig_size=1):
    """
    Return the Laurent expansion of the cluster algebra element
    (corresponding to the curve that crosses the arcs of ``crossed_arcs``)
    with respect to the initial seed corresponding to the triangulation ``T``.

    See :meth:`ClusterSeed.arc_laurent_expansion` and :meth:`ClusterSeed.loop`

    INPUT:

    - ``CT`` -- ClusterTriangulation class

    - ``crossed_arcs`` -- cluster variables corresponding to arcs that
      are crossed by gamma

    - ``first_triangle`` -- (default:``None``) the first triangle
      (a,b,c) crossed by curve

    - ``final_triangle`` -- (default:``None``) the last triangle
      (d,e,f) crossed by curve

    - ``is_arc`` -- (default:``None``) True if curve is between marked point/s

    - ``is_loop`` -- (default:``None``) True if curve is a loop in the
      interior of the surface

    - ``verbose`` -- (default:``False``) display the image of all
      perfect matchings of the snake graph if ``True``

    - ``fig_size`` -- (default:1) image size

    - ``is_principal`` -- (default:``False``) ``True`` if the seed corresponding to ``T``
      has principal coefficients (otherwise the seed is coefficient-free)

    ALGORITHM:

    See the perfect matching formula from Musiker-Schiffler-Williams'
    "Positivity for cluster algebras from surfaces" [MSW_Positivity]_
    (section 4) and "Bases for cluster algebras from surfaces"
    [MSW_Bases]_ (sections 3.1-3.2)

    To compute the Laurent polynomial expansion:

    #. Compute the sum ``all_sum`` of all weights of all perfect matchings
       in ``all_perfect_matchings``, i.e. the weight of a perfect matching
       is the product of all weights of edges in the matchings.
    #. The Laurent polynomial expansion corresponding to the curve is equal
       to ``all_sum`` over prod[``crossed_arcs``]

    .. TODO::

        In case parameter ``verbose`` is set to ``True``, user may
        want to specify that the perfect matchings should be labeled
        by user-given labels (and not the weights from the variables)

        Currently the number of perfect matchings that can be displayed
        without error messages are very limited.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import LaurentExpansionFromSurface
        sage: CT = ClusterTriangulation([(0,2,1),(0,4,3),(1,6,5)])
        sage: S = ClusterSeed(CT)
        sage: S1=S.mutate(0,inplace=False)
        sage: S1.cluster()
        [(x1*x3 + x2*x4)/x0, x1, x2, x3, x4, x5, x6]
        sage: S1.cluster_variable(0) ==\
        ....: LaurentExpansionFromSurface(CT,[CT.cluster_variable(0)],None,None,\
        ....: is_arc=True,is_loop=None,verbose=None,boundary_edges=None)
        True

        sage: SP = S.principal_extension()
        sage: SP.mutate(0)
        sage: SP.cluster()
        [(x2*x4*y0 + x1*x3)/x0, x1, x2, x3, x4, x5, x6, y0, y1, y2, y3, y4, y5, y6]

        sage: CTP = CT.principal_extension()
        sage: LaurentExpansionFromSurface(CTP,[CTP.cluster_variable(0)],None,None,True,None,None,None,None)
        (x2*x4*y0 + x1*x3)/x0

    Figure 6 of Musiker and Williams' "Matrix Formulae and Skein
    Relations for Cluster Algebras from Surfaces" [MW_MatrixFormulae]_
    which is a loop where tau_4, tau_1, tau_2, tau_3 = ``0``, ``1``, ``2``, ``3`` and
    b1, b2, b3, b4 = ``b4``, ``b5``, ``b6``, ``b7``.

    We pick tau_1 to be 1, and go clockwise, so that crossed_arcs =
    [1,2,3,0,1] ::

        sage: T = ClusterTriangulation([(1,2,'b4'),(1,0,'b5'),(0,3,'b6'),(2,3,'b7')], boundary_edges=['b4','b5','b6','b7'])
        sage: S = ClusterSeed(T)
        sage: c = [item for item in S.cluster()]
        sage: LaurentExpansionFromSurface(T,[c[1],c[2],c[3],c[0],c[1]],None,None,None,True,None,T._boundary_edges_vars,None)
        (x0*x1^2*x2 + x0*x2*x3^2 + x1^2 + 2*x1*x3 + x3^2)/(x0*x1*x2*x3)

        sage: TP = T.principal_extension()
        sage: LaurentExpansionFromSurface(TP,[c[1],c[2],c[3],c[0],c[1]],\
        ....: None,None,is_arc=None,is_loop=True,verbose=True,boundary_edges=T._boundary_edges_vars)
        **************** Perfect Matchings and Their Weights: *****************
        (x0*x2*x3^2*y0*y1*y2*y3 + x3^2*y0*y2*y3 + x0*x1^2*x2 + x1*x3*y0*y3 + x1*x3*y2*y3 + x1^2*y3)/(x0*x1*x2*x3)
    """
    #from copy import deepcopy
    if not isinstance(crossed_arcs,list) or len(crossed_arcs) < 1:
        raise ValueError('crossed_arcs should be a non-empty list object '
                         'of cluster variable/s')
    is_principal = CT._is_principal
    T = list(CT.weighted_triangulation())
    G_x = _snake_graph(T,crossed_arcs,first_triangle, final_triangle,
                     is_arc, is_loop, 1, boundary_edges)

    if is_principal:
        # Replace all x_i with the correct height monomial (see [MSW_Positivity]_ Thm 4.10)
        G_y = replace_x_with_y(CT, G_x)
    else:
        G_y = None

    all_matchings = GetAllMatchings(G_x)

    all_matchings_symmetricdifference_minpm = None
    if is_principal: # To compute principal coefficients
        all_matchings_symmetricdifference_minpm = []
        min_pm = GetMinimalMatching(G_x)
        for pm in all_matchings:
            minpm_symmetric_difference_pm = matching_symmetric_difference(min_pm, pm)

            #print ' '
            #print 'pm: ', pm
            #print 'minpm_symmetric_difference_pm: ', minpm_symmetric_difference_pm
            #print ' '

            all_matchings_symmetricdifference_minpm.append(minpm_symmetric_difference_pm)

    if verbose:
        print "**************** Perfect Matchings and Their Weights: *****************"
        from sage.plot.graphics import Graphics
        drawing = Graphics()
        xy=(0,0)
        draw_G = _draw_snake_graph(G_x, print_user_labels=False, xy=xy)

        for pos in range(0,len(all_matchings)):
            PM = all_matchings[pos]
            matching_weight = GetMonomialTerm(G_x, PM, None, is_loop)
            if is_principal:
                SD = all_matchings_symmetricdifference_minpm[pos]
                symmetric_difference_weight = GetCoefficientTerm(G_y, SD)
                matching_weight *= symmetric_difference_weight
            else:
                SD = None
            draw_G = _draw_snake_graph(G_x, print_user_labels=False,xy=xy)
            matching_drawing, xy = _draw_matching(perfect_matching=PM, \
            matching_weight=matching_weight,pos=pos,xy=xy, white_space=1, \
            print_user_labels=False, symmetric_difference=SD)
            drawing += matching_drawing + draw_G

        drawing.set_aspect_ratio(1)
        #fig_size = fig_size*2*n/3
        fig_size = 0.8*fig_size
        drawing.show(axes=False,figsize=[fig_size*len(crossed_arcs)*(len(all_matchings)+1), fig_size*len(crossed_arcs)])

    return SumOfMonomialTerms(G_x, all_matchings, boundary_edges, is_principal, G_y, all_matchings_symmetricdifference_minpm, is_loop)/ GetDenominator(G_x)

def replace_x_with_y(CT, G_x):
    """
    Replace the diagonal x_i in snake graph ``G_x`` with y_i replace_x_with_y

    EXAMPLES:

    Figure 10 and 11 of Musiker - Schiffler - Williams "Positivity for Cluster Algebras from Surfaces" [MSW_Positivity]_ ::

        sage: from sage.combinat.cluster_algebra_quiver.surface import replace_x_with_y

        sage: T = ClusterTriangulation([(2,3,11),(2,1,1),(4,3,12),(0,4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,13)], boundary_edges=[11,12,13,0])
        sage: T = T.principal_extension()
        sage: c = [item for item in T.cluster()]
        sage: r = c[1-1]
        sage: ell = c[2-1]* r
        sage: snakegraph_x = T.list_snake_graph([ell,(r,'counterclockwise'), ell, c[3-1],c[4-1],c[5-1],c[6-1]], first_tile_orientation=-1, user_labels=False)
        sage: snakegraph_x
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

        sage: replace_x_with_y(T, snakegraph_x)
        [[(-1, (x2, y1, b11)),
        (-2, ((x0, 'counterclockwise'), y1, (x0, 'clockwise')), 'RIGHT')],
        [(1, (x0*x1, (y0/y1, 'counterclockwise'), (x0, 'clockwise'))),
        (2, ((x0, 'counterclockwise'), (y0/y1, 'clockwise'), x0*x1), 'ABOVE')],
        [(-1, ((x0, 'counterclockwise'), y1, (x0, 'clockwise'))),
        (-2, (x2, y1, b11), 'RIGHT')],
        [(1, (x0*x1, y2, b11)), (2, (x3, y2, b12), 'RIGHT')],
        [(-1, (x2, y3, b12)), (-2, (x4, y3, b10), 'RIGHT')],
        [(1, (x3, y4, b10)), (2, (x9, y4, x5), 'ABOVE')],
        [(-1, (x9, y5, x4)), (-2, (x6, y5, x8), 'ABOVE')]]

    """
    G_y = []
    selffoldedtriangles = []
    ell_loops = []
    for tri in CT.weighted_triangulation():
        if isinstance(tri[0],tuple):
            #selffoldedtriangles.append(tri)
            ell_loops.append(tri[2])

    for tile_pos in range(0,len(G_x)):
        triangle_bottom = G_x[tile_pos][0][1]
        triangle_top = G_x[tile_pos][1][1]
        tile_dir = G_x[tile_pos][1][2]
        diagonal_x = triangle_bottom[1]
        if isinstance(diagonal_x, tuple):
            #print 'diagonal_x tuple: ', diagonal_x
            label_of_radius_x = CT._get_map_variable_to_label(diagonal_x[0])
            position_of_radius_x = CT.get_edge_position(label_of_radius_x)
            radius_y = CT.cluster()[CT._n+position_of_radius_x]

            triangle_rrl = _get_triangle(CT.triangles(), label_of_radius_x, None)[0]
            r,r,label_of_ell_x = is_selffolded(triangle_rrl)
            position_of_ell_x = CT.get_edge_position(label_of_ell_x)

            r_notched_y = CT.cluster()[CT._n+position_of_ell_x]
            diagonal_y = radius_y/r_notched_y
            diagonal_y_bottom = (diagonal_y, diagonal_x[1])
            if diagonal_x[1] == 'clockwise':
                clockwise_or_counterclockwise = 'counterclockwise'
            else:
                clockwise_or_counterclockwise = 'clockwise'
            diagonal_y_top = (diagonal_y, clockwise_or_counterclockwise)
            #print 'diagonal_y = (diagonal_y, diagonal_x[1]): ', diagonal_y
        elif diagonal_x in ell_loops:
            #pos_of_ell_in_list = ell_loops.index(diagonal_x)
            #triangle_rrl = selffoldedtriangles[pos_of_ell_in_list]
            x_ell = diagonal_x
            #x_r = triangle_rrl[0][0]
            label_of_x_ell = CT._get_map_variable_to_label(x_ell)
            #label_of_x_r = CT._get_map_variable_to_label(x_r)
            position_of_x_rnotched = CT.get_edge_position(label_of_x_ell)
            #position_of_x_r = CT.get_edge_position(label_of_x_r)
            y_rnotched = CT.cluster()[CT._n+position_of_x_rnotched]
            #y_r = CT.cluster()[CT._n+position_of_x_r]
            diagonal_y_bottom = diagonal_y_top = y_rnotched
        else:
            label_of_diagonal_x = CT._get_map_variable_to_label(diagonal_x)
            position_of_diagonal_x = CT.get_edge_position(label_of_diagonal_x)
            diagonal_y_bottom = diagonal_y_top = CT.cluster()[CT._n+position_of_diagonal_x]

        G_y_triangle_bottom = (G_x[tile_pos][0][0],(triangle_bottom[0], diagonal_y_bottom, triangle_bottom[2]))
        G_y_triangle_top = (G_x[tile_pos][1][0],(triangle_top[0], diagonal_y_top, triangle_top[2]), tile_dir)

        G_y.append([G_y_triangle_bottom, G_y_triangle_top])
        #print 'G_y: ', G_y
    if len(G_y) == 0:
        raise ValueError('Bug in replace_x_with_y. Snake graph G_y is empty')
    return G_y

def matching_symmetric_difference(min_pm, pm):
    """
    Return the symmetric difference of the (snake graph) perfect matchings
    ``min_pm`` and ``pm``

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import matching_symmetric_difference, GetMinimalMatching, _snake_graph
        sage: T = ClusterTriangulation([(0,2,1),(0,4,3),(1,6,5)])
        sage: G = _snake_graph(T._weighted_triangulation,[T.cluster_variable(0)],
        ....: first_triangle=None, final_triangle=None, is_arc=True, is_loop=False,
        ....: first_tile_orientation=1, boundary_edges=None)
        sage: Pmin = GetMinimalMatching(G)
        sage: Pmin
        [['minimal PM'], [[(1, 0, 1, 0), 'ABOVE']]]
        sage: matching_symmetric_difference(Pmin,Pmin)
        [['minimal PM'], [[(0, 0, 0, 0), 'ABOVE']]]
        sage: matching_symmetric_difference(Pmin,[[0], [[(0, 1, 0, 1), 'ABOVE']]])
        [[0], [[(1, 1, 1, 1), 'ABOVE']]]
    """
    minpm_union_pm = matching_union(min_pm, pm)
    minpm_intersect_pm = matching_intersection(min_pm, pm)
    minpm_symmetric_difference_pm = matching_subtract(minpm_union_pm, minpm_intersect_pm)
    return minpm_symmetric_difference_pm

def matching_union(pmA, pmB):
    """
    Return the union of the (snake graph) perfect matchings
    ``min_pm`` and ``pm``

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import matching_union, GetMinimalMatching, _snake_graph
        sage: T = ClusterTriangulation([(0,2,1),(0,4,3),(1,6,5)])
        sage: G = _snake_graph(T._weighted_triangulation,[T.cluster_variable(0)],
        ....: first_triangle=None, final_triangle=None, is_arc=True, is_loop=False,
        ....: first_tile_orientation=1, boundary_edges=None)
        sage: Pmin = GetMinimalMatching(G)
        sage: Pmin
        [['minimal PM'], [[(1, 0, 1, 0), 'ABOVE']]]
        sage: matching_union(Pmin,Pmin) == Pmin
        True
        sage: matching_union(Pmin,[[0], [[(0, 1, 0, 1), 'ABOVE']]])
        [[0], [[(1, 1, 1, 1), 'ABOVE']]]
    """
    tile_count = len(pmA[1]) # Always put pmA as the minimal matching
    pmC = []
    tileC_matching = [-1,-1,-1,-1]
    for tile_pos in range(0,tile_count):
        tileA = pmA[1][tile_pos]
        tileB = pmB[1][tile_pos]
        #print 'tileA:', tileA
        for pos in range(0,4):
            #print 'tileA[0][pos]: ', tileA[0][pos]
            tileC_matching[pos] = tileA[0][pos] or tileB[0][pos]
        tileC = [(tileC_matching[0],tileC_matching[1],tileC_matching[2],tileC_matching[3]), tileA[1]]
        pmC.append(tileC)

    #print [pmB[0], pmC]
    return [pmB[0], pmC]

def matching_intersection(pmA, pmB):
    """
    Return the intersection of the (snake graph) perfect matchings
    ``min_pm`` and ``pm``

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import matching_intersection, GetMinimalMatching, _snake_graph
        sage: T = ClusterTriangulation([(0,2,1),(0,4,3),(1,6,5)])
        sage: G = _snake_graph(T._weighted_triangulation,[T.cluster_variable(0)],
        ....: first_triangle=None, final_triangle=None, is_arc=True, is_loop=False,
        ....: first_tile_orientation=1, boundary_edges=None)
        sage: Pmin = GetMinimalMatching(G)
        sage: Pmin
        [['minimal PM'], [[(1, 0, 1, 0), 'ABOVE']]]
        sage: matching_intersection(Pmin,Pmin) == Pmin
        True
        sage: matching_intersection(Pmin,[[0], [[(0, 1, 0, 1), 'ABOVE']]])
        [[0], [[(0, 0, 0, 0), 'ABOVE']]]
    """
    tile_count = len(pmA[1]) # Always put pmA as the minimal matching
    pmC = []
    tileC_matching = [-1,-1,-1,-1]
    for tile_pos in range(0,tile_count):
        tileA = pmA[1][tile_pos]
        tileB = pmB[1][tile_pos]
        #print 'tileA:', tileA
        for pos in range(0,4):
            #print 'tileA[0][pos]: ', tileA[0][pos]
            tileC_matching[pos] = tileA[0][pos] and tileB[0][pos]
        tileC = [(tileC_matching[0],tileC_matching[1],tileC_matching[2],tileC_matching[3]), tileA[1]]
        pmC.append(tileC)
    return [pmB[0], pmC]

def matching_subtract(bigger_pm, smaller_pm):
    """
    Return subtraction of (snake graph) perfect matchings, ``bigger_pm`` minus ``smaller_pm``.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import matching_subtract, GetMinimalMatching, matching_intersection, matching_union, GetAllMatchings
        sage: affineD = ClusterTriangulation([('5','6','4'),('6','7','10'),('7','1','11'),('1','5','2'),('2','3','8'),('3','4','9')], boundary_edges=['8','9','10','11'])
        sage: bandgraph = affineD.list_band_graph(['1','2','3','4','6','7','1'])

        sage: Pmin = GetMinimalMatching(bandgraph)
        sage: P24 = [[2, 4],\
        ....: [[(1, 0, 1, 0), 'RIGHT'],\
        ....: [(0, 1, 0, 0), 'RIGHT'],\
        ....: [(0, 1, 0, 1), 'ABOVE'],\
        ....: [(0, 0, 0, 0), 'ABOVE'],\
        ....: [(0, 1, 0, 1), 'ABOVE'],\
        ....: [(0, 0, 1, 0), 'RIGHT']]]
        sage: Pmin_intersect_P24 = matching_intersection(Pmin,P24)
        sage: Pmin_intersect_P24
        [[2, 4],
        [[(1, 0, 1, 0), 'RIGHT'],
        [(0, 0, 0, 0), 'RIGHT'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 0), 'RIGHT']]]
        sage: Pmin_union_P24 = matching_union(Pmin,P24)
        sage: Pmin_union_P24
        [[2, 4],
        [[(1, 0, 1, 0), 'RIGHT'],
        [(0, 1, 0, 0), 'RIGHT'],
        [(1, 1, 0, 1), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 1, 1, 1), 'RIGHT']]]
        sage: matching_subtract(Pmin_union_P24,Pmin_intersect_P24)
        [[2, 4],
        [[(0, 0, 0, 0), 'RIGHT'],
        [(0, 1, 0, 0), 'RIGHT'],
        [(1, 1, 0, 1), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 1, 1, 1), 'RIGHT']]]

        sage: P025 = [[0,2,5], \
        ....: [[(0,1,0,1), 'RIGHT'],\
        ....: [(0,1,0,1), 'RIGHT'],\
        ....: [(0,1,0,1), 'ABOVE'],\
        ....: [(0,0,1,0), 'ABOVE'],\
        ....: [(1,0,1,0), 'ABOVE'],\
        ....: [(1,0,1,0), 'RIGHT']]]
        sage: Pmin_intersect_P025 = matching_intersection(Pmin,P025)
        sage: Pmin_intersect_P025
        [[0, 2, 5],
        [[(0, 0, 0, 0), 'RIGHT'],
        [(0, 0, 0, 0), 'RIGHT'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 0), 'RIGHT']]]
        sage: Pmin_union_P025 = matching_union(Pmin,P025)
        sage: Pmin_union_P025
        [[0, 2, 5],
        [[(1, 1, 1, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(1, 1, 0, 1), 'ABOVE'],
        [(0, 1, 1, 1), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE'],
        [(1, 1, 1, 1), 'RIGHT']]]
        sage: matching_subtract(Pmin_union_P025,Pmin_intersect_P025) == Pmin_union_P025
        True
    """
    pmA = bigger_pm
    pmB = smaller_pm
    tile_count = len(pmA[1]) # Always put pmA as the minimal matching
    pmC = []
    tileC_matching = [-1,-1,-1,-1]
    for tile_pos in range(0,tile_count):
        tileA = pmA[1][tile_pos]
        tileB = pmB[1][tile_pos]
        #print 'tileA:', tileA
        for pos in range(0,4):
            #print 'tileA[0][pos]: ', tileA[0][pos]
            tileC_matching[pos] = tileA[0][pos] - tileB[0][pos]
        tileC = [(tileC_matching[0],tileC_matching[1],tileC_matching[2],tileC_matching[3]), tileA[1]]
        pmC.append(tileC)
    return [pmB[0], pmC]

def GetAllMatchings(G):
    """
    Return all perfect matchings of the snake/band graph ``G``.

    INPUT:

    - ``G`` -- snake/band graph, see
      :meth:`ClusterTriangulation.list_snake_graph` and
      :meth:`ClusterTriangulation.list_band_graph`

    ALGORITHM:

    To produce ``all_perfect_matchings`` of a snake graph/band graph:

    #. Produce the minimal matching min_pm, i.e. the perfect matching
       containing:

        #. the floor edge of the first tile (with +1 orientation),
        #. and only boundary edges of the snake graph.

    #. Produce more perfect matching/s by flipping (one at a time) every
       flippable tile (i.e. (1,0,1,0) to (0,1,0,1)) of min_pm.
    #. Continue flipping (one at a time) every flippable tile of current
       matchings.
    #. When this process stops producing new perfect matchings, quit the loop.

    EXAMPLES:

    An ideal triangulation of a once-punctured square having 2 radii, and the boundary edges are labeled 4,5,6,7::

        sage: from sage.combinat.cluster_algebra_quiver.surface import GetAllMatchings
        sage: once_punctured_square = [(1,'b7','b4'),(1,'b5',2),('b6',0,3),(2,3,0),(0,3,'b6'),['b7','b4',1]]
        sage: T = ClusterTriangulation(once_punctured_square, boundary_edges=['b4','b5','b6','b7'])
        sage: S = ClusterSeed(T)
        sage: crossed = [S.x(1),S.x(2),S.x(3)]
        sage: gamma = S.mutate([3,2,1],inplace=False).cluster_variable(1)
        sage: snakegraph = T.list_snake_graph(crossed,user_labels=False)
        sage: sum(gamma.numerator().coefficients()) == len(GetAllMatchings(snakegraph))
        True

        sage: T=ClusterTriangulation([[1,2,3],[1,2,3]])
        sage: G=T.list_snake_graph([T.cluster()[0]],user_labels=False)
        sage: G
        [[(1, (x2, x0, x1)), (2, (x2, x0, x1), 'ABOVE')]]
        sage: GetAllMatchings(G)
        [[['minimal PM'], [[(1, 0, 1, 0), 'ABOVE']]],
        [['maximal PM'], [[(0, 1, 0, 1), 'ABOVE']]]]
    """
    from sage.combinat.combinat import fibonacci # todo: eventually remove this after we are sure we don't need the upper bound

    MinMatching = GetMinimalMatching(G)  # Return [['minimal PM'], [minimal matching with directions]]
    #print 'MinMatching = GetMinimalMatching(G) : ', MinMatching

    tile_flip_max = len(G) # The maximal number of flips is equal to the number of tiles of G.
    # The maximal number of flips produces the maximal matching, i.e. the opposite of the minimal matching.
    # Other matchings are produced by flips fewer than the maximal number.

    old_matchings = []
    current_matchings = [MinMatching]

    if len(G) == 1:
        horizontal_PM = [FlipTile(MinMatching[1][0])]
        all_matchings =[MinMatching, [['maximal PM'], horizontal_PM]]
    else:
        for loop_count in range(0,tile_flip_max):
            new_matchings = FlipAllFlippableTilesInList(current_matchings) # All tiles (except the ones that were flipped last) are flipped
            if new_matchings != []:
                old_current_new_matchings = UniqueList(old_matchings + current_matchings + new_matchings)

                # The same matching can be produced via different flips, so we merge all of them into one item
                new_matchings_corrected_indices = \
                GetMoreLastFlippedTilesInList(new_matchings, old_current_new_matchings)

                old_matchings = UniqueList(old_matchings + current_matchings)
                current_matchings = new_matchings_corrected_indices
            else:
                print 'WARNING: This else break should be removed.'
                print 'All the perfect matchings of G are found when the maximal matching is produced (after tile_flip_max) flips.'
                print 'loop count: ', loop_count
                break

        all_matchings = UniqueList(old_matchings + current_matchings)

    return all_matchings


def UniqueList(in_list):
    """
    Return the list ``in_list`` with any multiple entries removed

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import UniqueList
        sage: UniqueList([5,4,3,2,1,'a','b','c','d','e',4,3,2,1,'d','c','c','c'])
        [5, 4, 3, 2, 1, 'a', 'b', 'c', 'd', 'e']
    """
    new_list = []
    for elt in in_list:
        if elt not in new_list:
            new_list.append(elt)
    return new_list

########################################### END: LAURENT EXPANSION ####################
##########################################################################################



##########################################################################################
####### BEGINS: FUNCTIONS THAT EXTRACT WEIGHTS OF PERFECT MATCHINGS ######################
##########################################################################################

def PartitionIntoTuples(L):
    """
    Partition a list into a list of tuples. This function acts like Partition[L,2] in Mathematica.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import PartitionIntoTuples
        sage: L = [1,2,3,4,5,6,7,8,9,10]
        sage: PartitionIntoTuples(L)
        [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]]
        sage: L = [2,3,4,5,6,7,8,9,10]
        sage: PartitionIntoTuples(L)
        [[2, 3], [4, 5], [6, 7], [8, 9], [10]]
    """
    ListTuples = []
    if len(L) % 2 == 1:
        for pos in range(0,len(L)-1,2):
            ListTuples.append([L[pos],L[pos+1]])
        ListTuples.append([L[-1]])
    else:
        for pos in range(0,len(L),2):
            ListTuples.append([L[pos],L[pos+1]])
    return ListTuples

def GetDenominator(G):
    """
    Return the product of the variables corresponding to arcs that are crossed by curve gamma
    (this is the denominator of the Laurent polynomial expansion returned by
    :func:`~sage.combinat.cluster_algebra_quiver.surface.LaurentExpansionFromSurface`)

    To debug, original Mathematica function: todaslasdiagonales

    INPUT:

    - ``G`` -- snake/band graph information, see :meth:`ClusterTriangulation.list_band_graph` or :meth:`ClusterTriangulation.list_snake_graph`

    EXAMPLES:

    Figure 6 of Musiker and Williams' "Matrix Formulae and Skein
    Relations for Cluster Algebras from Surfaces" [MW_MatrixFormulae]_
    tau_4, tau_1, tau_2, tau_3 = 0,1,2,3, and we pick tau_1=1 to be
    the first and last arc crossed::

        sage: from sage.combinat.cluster_algebra_quiver.surface import GetDenominator
        sage: T = ClusterTriangulation([(1,2,'b1'),(1,0,'b2'),(0,3,'b3'),(2,3,'b4')],boundary_edges=['b1','b2','b3','b4']) # Counterclockwise triangulation
        sage: c = [item for item in T.cluster()]
        sage: BG = T.list_band_graph([c[1], c[2], c[3], c[0], c[1]],user_labels=False)
        sage: GetDenominator(BG)
        x0*x1*x2*x3
    """
    denom = 1
    for pos in range(0,len(G)):
        diagonal = G[pos][0][1][1]
        if type(diagonal) in [list,tuple]:
            diagonal = diagonal[0]
        denom *= diagonal
    return denom

def SumOfMonomialTerms(snakegraph, all_matchings, boundary_edges=None, is_principal=False, snakegraph_y=None, all_symmetric_differences=None, is_loop=False):
    """
    Return sum of all monomial terms (corresponding to all perfect
    matchings of ``snakegraph``) which is the numerator of the Laurent
    polynomial expansion returned by
    :func:`~sage.combinat.cluster_algebra_quiver.surface.LaurentExpansionFromSurface`

    To debug, original Mathematica function: terminoPolinomioSobreListaDeConfig

    INPUT:

    - ``snakegraph`` -- see
      :meth:`ClusterTriangulation.list_band_graph` or
      :meth:`ClusterTriangulation.list_snake_graph`

    - ``all_matchings`` -- all perfect matchings for ``snakegraph``

    - ``boundary_edges`` -- (default:``None``) variables b_i (which is
      set to 1) corresponding to boundary edges

    EXAMPLES:

    Affine A(2,2) triangulation from Figure 3 of Shiffler-Thomas' paper :arxiv:`abs/0712.4131` where tau_i = i and tau_8 is labeled 0::

        sage: from sage.combinat.cluster_algebra_quiver.surface import SumOfMonomialTerms, GetDenominator, GetAllMatchings, _get_weighted_edge
        sage: T = ClusterTriangulation([('b7',4,3),(4,1,'b5'),(3,'b6',2),(2,1,'b8')], boundary_edges=['b5','b6','b7','b8'])
        sage: S = ClusterSeed(T)
        sage: crossed = [S.x(0),S.x(1),S.x(2),S.x(3),S.x(0)]
        sage: boundary_variables = [_get_weighted_edge(b,T._map_label_to_variable) for b in ['b5','b6','b7','b8']]
        sage: G = T.list_snake_graph(crossed,user_labels=False)
        sage: all_matchings = GetAllMatchings(G)
        sage: gamma_numerator = SumOfMonomialTerms(G, all_matchings, boundary_edges=boundary_variables)
        sage: gamma_denominator = GetDenominator(G)
        sage: S.mutate([0,2,3,1,2], inplace=False).cluster_variable(2) == gamma_numerator/gamma_denominator
        True
    """
    sumTerms = 0
    pos = 0
    for matching in all_matchings:
        if is_principal:
            symmetric_difference = all_symmetric_differences[pos]
            sumTerms += GetMonomialTerm(snakegraph, matching, boundary_edges, is_loop) * \
            GetCoefficientTerm(snakegraph_y, symmetric_difference)
        else:
            sumTerms += GetMonomialTerm(snakegraph, matching, boundary_edges, is_loop)
        pos += 1
    return sumTerms

def ExtractWeight(tile, abcd, is_final_tile):
    """
    Return the list of cluster variables which label the marked edges
    of ``tile`` (such that an interior edge is ignored if it borders
    the next tile).

    To debug, see Mathematica function: rescatePositivo4

    INPUT:

    - ``tile`` -- a tile from a band/snake graph in the format
      [[1,(x,y,z)],[2,(b,y,a),DIR]]

    - ``abcd`` -- whether or not a matching contains an edge of the
      tile, in the form (bottom,right,top,left)

    - ``is_final_tile`` -- True if the tile is the final tile of the
      snake/band graph

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import ExtractWeight
        sage: T = ClusterTriangulation([(1,0,4),(3,0,2)])
        sage: G = T.list_snake_graph([T.cluster()[0]], user_labels=False)
        sage: G[0]
        [(1, (x1, x0, x4)), (2, (x3, x0, x2), 'ABOVE')]
        sage: ExtractWeight(G[0],(1,0,1,0),True)
        [x1, x3]
        sage: ExtractWeight(G[0],(0,1,0,1),True)
        [x2, x4]
    """
    x = tile[0][1][0]
    if type(x) in [tuple, list]: x=x[0]
    z = tile[0][1][2]
    if type(z) in [tuple, list]: z=z[0]
    b = tile[1][1][0]
    if type(b) in [tuple, list]: b=b[0]
    a = tile[1][1][2]
    if type(a) in [tuple, list]: a=a[0]

    DIR = tile[1][2]

    (i1,i2,i3,i4) = (abcd[0], abcd[1], abcd[2], abcd[3])

    weights = [x*i1, a*i2, b*i3, z*i4]

    if i2 == 1 and DIR == RIGHT and not is_final_tile:
        weights = [x*i1, b*i3, z*i4]

    if i3 == 1 and DIR == ABOVE and not is_final_tile:
        weights = [x*i1, a*i2, z*i4]

    while 0 in weights:
        weights.remove(0)

    return weights

def GetMonomialTerm(snakegraph, PM, boundary_edges, is_loop):
    """
    Return the monomial term corresponding to the input perfect matching ``PM``.

    To debug, see Mathematica function: terminoPolinomio

    Called by :func:`~sage.combinat.cluster_algebra_quiver.surface.SumOfMonomialTerms and :func:`~sage.combinat.cluster_algebra_quiver.surface.LaurentExpansionFromSurface`

    INPUT:

    - ``snakegraph`` -- snake/band graph

    - ``PM`` -- a perfect matching of a band/snake graph

    EXAMPLES:

    Figure 10 of [MSW_Positivity]_::

        sage: from sage.combinat.cluster_algebra_quiver.surface import GetMonomialTerm, FlipAllFlippableTilesInList, FlipAllFlippableTiles, GetMinimalMatching
        sage: thrice_punctured_square = [('r','r','ell'),('11','ell','3'),('3','12','4'),\
        ....: ('4','5','14'),('5','6','10'),('6','7','9'),('8','10','9'),('7','13','8')]
        sage: T = ClusterTriangulation(thrice_punctured_square, boundary_edges=['11','12','13','14'])
        sage: crossed = [T._get_map_label_to_variable(e) for e in ['5','6','7','8','9','6','5']]
        sage: G = T.list_snake_graph(crossed, first_tile_orientation=1, user_labels=False)
        sage: pm_a, pm_b = FlipAllFlippableTilesInList([GetMinimalMatching(G)])
        sage: pm_a
        [[2],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE']]]
        sage: GetMonomialTerm(G,pm_a,boundary_edges=T._boundary_edges_vars,is_loop=False)
        x2^2*x3^2*x7^3
        sage: pm_b
        [[5],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 1), 'RIGHT'],
        [(1, 0, 1, 0), 'RIGHT'],
        [(0, 1, 0, 0), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE']]]
        sage: GetMonomialTerm(G,pm_b,boundary_edges=T._boundary_edges_vars, is_loop=None)
        x0*x2^2*x3*x4*x5*x6*x7
    """
    tile_weights = []
    if boundary_edges is None:
        boundary_edges = []

    if len(snakegraph) == len(PM[1]):  # this should always be equal
        total_weight = []
        is_final_tile = False
        for pos in range(0,len(snakegraph)):
            if pos == len(snakegraph)-1 :
                is_final_tile = True
            abcd = PM[1][pos][0] # abcd = (bottom,right,top,left)
            tile_weight = ExtractWeight(snakegraph[pos], abcd, is_final_tile)
            tile_weights.append(tile_weight)
    else:
        print 'warning: see GetMonomialTerm'

    #matching_weight = Multiply elements in tile_weights
    matching_weight = 1
    #print 'tile_weights:', tile_weights
    for var in [item for sublist in tile_weights for item in sublist]:
        if var not in boundary_edges:
            matching_weight = matching_weight * var

    first_triangle = snakegraph[0][0][1] # the very first triangle [x,y,z]
    final_triangle = snakegraph[-1][1][1] # the very last triangle [x,y,z]

    # a,b,c and x,y,z are counterclockwise and b,y are diagonals so that
    # a = bottom, b= diagonal, c = left
    # x = top, y=diagonal, z=right
    a1 = first_triangle[0]
    b1 = first_triangle[1]
    c1 = first_triangle[2]

    xn = final_triangle[0]
    yn = final_triangle[1]
    zn = final_triangle[2]

    first_tile_matching = PM[1][0][0] #e.g. {1,0,1,0}
    final_tile_matching = PM[1][-1][0]
    myarray = 1

    if is_loop:
        #4 cases:
        if a1==xn and b1 == zn and c1 == yn: # second case, bottom of first tile == top of final tile, and diagonal of first tile == right of final tile
            if (first_tile_matching[0] == 0 # bottom of first tile is not marked
            and final_tile_matching[2] == 0): # top of final tile is not marked
                myarray = 0
            else:
                if a1 in boundary_edges:
                    myarray = 1
                else:
                    myarray = 1/a1

        if a1 == zn and b1 == xn and c1 == yn: # first case
            if (first_tile_matching[0] == 0 # bottom of the first tile is not marked
            and final_tile_matching[1] == 0): # right of final tile is not marked
                myarray = 0
            else:
                if a1 in boundary_edges:
                    myarray = 1
                else:
                    myarray = 1/a1

        if a1 == yn and b1 == zn and c1 == xn: # third case
            if (first_tile_matching[3] == 0 # left of first tile is not marked
            and final_tile_matching[2] == 0): # top of final tile is not marked
                myarray = 0
            else:
                if c1 in boundary_edges:
                    myarray = 1
                else:
                    myarray = 1/c1

        if a1 == yn and b1 == xn and c1 == zn: # fourth case
            if (first_tile_matching[3] == 0 # left of first tile is not marked
            and final_tile_matching[1] == 0): # right of final tile is not marked
                myarray = 0
            else:
                if c1 in boundary_edges:
                    myarray = 1
                else:
                    myarray = 1/c1

    return matching_weight  * myarray

def GetCoefficientTerm(snakegraph, SD):
    """
    Return the coefficient product corresponding to the symmetric difference ``SD``
    is a union of disjoint cycle/s.

    Called by :func:`~sage.combinat.cluster_algebra_quiver.surface.SumOfMonomialTerms and :func:`~sage.combinat.cluster_algebra_quiver.surface.LaurentExpansionFromSurface`

    INPUT:

    - ``snakegraph`` -- snake/band graph

    - ``SD`` -- a union of disjoint cycles along edges of a band/snake graph

    EXAMPLES:

    Figure 10 of [MSW_Positivity]_::

        sage: from sage.combinat.cluster_algebra_quiver.surface import GetCoefficientTerm, \
        FlipAllFlippableTilesInList, FlipAllFlippableTiles, GetMinimalMatching, \
        matching_symmetric_difference, replace_x_with_y
        sage: thrice_punctured_square = [('r','r','ell'),('11','ell','3'),('3','12','4'),\
        ('4','5','14'),('5','6','10'),('6','7','9'),('8','10','9'),('7','13','8')]
        sage: TP = ClusterTriangulation(thrice_punctured_square, boundary_edges=['11','12','13','14']).principal_extension()
        sage: crossed = [TP._get_map_label_to_variable(e) for e in ['5','6','7','8','9','6','5']]
        sage: G = TP.list_snake_graph(crossed, first_tile_orientation=1, user_labels=False)
        sage: pm_min = GetMinimalMatching(G)
        sage: pm_a, pm_b = FlipAllFlippableTilesInList([pm_min])
        sage: G_y = replace_x_with_y(TP, G)
        sage: sd_a = matching_symmetric_difference(pm_min,pm_a)
        sage: sd_a
        [[2],
        [[(0, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 0), 'RIGHT'],
        [(1, 1, 1, 1), 'RIGHT'],
        [(0, 0, 0, 1), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE']]]
        sage: GetCoefficientTerm(G_y,sd_a)
        y5
        sage: sd_b = matching_symmetric_difference(pm_min,pm_b)
        sage: sd_b
        [[5],
        [[(0, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 0), 'RIGHT'],
        [(0, 0, 0, 0), 'RIGHT'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE'],
        [(1, 1, 1, 1), 'ABOVE'],
        [(1, 0, 0, 0), 'ABOVE']]]
        sage: GetCoefficientTerm(G_y,sd_b)
        y4
   """
    diagonal_weights = []
    #if boundary_edges is None:
    #    boundary_edges = []

    if len(snakegraph) == len(SD[1]):  # this should always be equal
        total_weight = []
        is_final_tile = False
        for pos in range(0,len(snakegraph)):
            DIR = snakegraph[pos][1][2]
            #print snakegraph[pos]
            if pos == len(snakegraph)-1 :
                is_final_tile = True
            abcd = SD[1][pos][0] # abcd = (bottom,right,top,left)
            bottom, right, top, left = abcd[0:4]
            #print 'abcd:', abcd
            #tile_weight = ExtractWeight(snakegraph[pos], abcd, is_final_tile)
            if (bottom and right) or (bottom and right) or (top and right) or (top and left):
                diagonal_weight = snakegraph[pos][1][1][1]
            elif (bottom and top) and (DIR == RIGHT or right or left):
                diagonal_weight = snakegraph[pos][1][1][1]
            elif (left and right) and (DIR == ABOVE or bottom or top):
                diagonal_weight = snakegraph[pos][1][1][1]
            else:
                diagonal_weight = 1
            if isinstance(diagonal_weight,tuple):
                diagonal_weight = diagonal_weight[0]
            diagonal_weights.append(diagonal_weight)
    else:
        print 'warning: see GetCoefficientTerm'

    #diagonal_weight = Multiply elements in diagonal_weights
    diagonal_weight = 1
    #print 'diagonal_weights: ', diagonal_weights
    for var in diagonal_weights:
        diagonal_weight = diagonal_weight * var

    first_triangle = snakegraph[0][0][1] # the very first triangle [x,y,z]
    final_triangle = snakegraph[-1][1][1] # the very last triangle [x,y,z]

    # a,b,c and x,y,z are counterclockwise and b,y are diagonals so that
    # a = bottom, b= diagonal, c = left
    # x = top, y=diagonal, z=right
    #a1 = first_triangle[0]
    #b1 = first_triangle[1]
    #c1 = first_triangle[2]

    #xn = final_triangle[0]
    #yn = final_triangle[1]
    #zn = final_triangle[2]

    first_tile_matching = SD[1][0][0] #e.g. {1,0,1,0}
    final_tile_matching = SD[1][-1][0]
    myarray = 1

    return diagonal_weight  * myarray

##########################################################################################
####### ENDS: FUNCTIONS THAT EXTRACT WEIGHTS OF PERFECT MATCHINGS ########################
##########################################################################################


#########################################################################################
################### BEGIN: FUNCTIONS FOR CONSTRUCTING BAND/SNAKE GRAPH ##################
#########################################################################################

def _get_first_final_triangles(T,crossed_arcs, first_triangle, final_triangle, is_arc, is_loop):
    """
    Search for and return [first_triangle, last_triangle], the first and final triangles crossed by curve.

    Called by :func:`~sage.combinat.cluster_algebra_quiver.surface_snake_graph`
    and :func:`~sage.combinat.cluster_algebra_quiver.surface_lifted_polygon`

    EXAMPLES:

    Example 3.6 from Dupont - Thomas' Atomic Basis in Types A and Affine A :arXiv:`1106.3758`::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _get_first_final_triangles
        sage: T = ClusterTriangulation([(0,1,'b2'),(0,1,'b3')], boundary_edges=['b2','b3'])
        sage: _get_first_final_triangles(T._triangulation, [0,1,0,1,0], first_triangle=(0,1,'b2'), final_triangle=None, is_arc=False, is_loop=True)
        [(0, 1, 'b2'), (0, 1, 'b2')]

    A once-punctured square's triangulation with self-folded triangle,
    border edges are labeled 4,5,6,7, 2nd triangulation in oral
    paper ell-loop is labeled 3, radius is labeled 0::

        sage: T = ClusterTriangulation([(1,7,4),(1,5,2),(2,3,6),(3,0,0)], boundary_edges=[4,5,6,7])
        sage: _get_first_final_triangles(T._triangulation, [1,2,3,(0,'counterclockwise'),3], None, None, is_arc=True, is_loop=False)
        [(1, 7, 4), (2, 3, 6)]
    """
    if (is_arc == is_loop):
        raise ValueError('is_arc and is_loop cannot have the same value')
    if is_loop:
        if len(crossed_arcs) < 2:
            raise ValueError('Since gamma is a loop, user needs to specify at least two crossing arcs.',\
            'If gamma only crosses one arc tau, then enter [tau,tau]')
        if crossed_arcs[0]!=crossed_arcs[-1]:
            raise ValueError('Since gamma is a loop, user needs to specify a sequence of tau_1, tau_2, ..., tau_d where tau_1=tau_d')

    # If only one arc tau of T is crossed, then we assign the two triangles
    # that have an edge tau to be first_triangle and final_triangle
    if len(crossed_arcs)==1 or (is_loop==True and len(crossed_arcs)==2):
        [first_triangle, final_triangle] = \
        try_to_find_end_triangles_for_one_crossed_arc(T, crossed_arcs[0], first_triangle, final_triangle, is_arc, is_loop)

    # We try to determine first_triangle and final_triangle by tau_1, tau_2 and tau_{d-1} and tau_{d}
    else:
        if is_loop:
            if not(first_triangle is None) and final_triangle is None:
                final_triangle = first_triangle
            elif first_triangle is None and not(final_triangle is None):
                first_triangle = final_triangle

        #if first_triangle == None:
        first_triangle = try_to_find_end_triangle(T,crossed_arcs, 'First', is_arc, is_loop, first_triangle)

        #if is_loop == True and final_triangle == None:
        if is_loop and final_triangle is None:
            final_triangle = first_triangle

        #if final_triangle == None:
        final_triangle = try_to_find_end_triangle(T,crossed_arcs, 'Final', is_arc, is_loop, final_triangle)

        if is_loop:
            #if are_triangles_equal(first_triangle, final_triangle) == False:
            if not are_triangles_equal(first_triangle, final_triangle):
                raise ValueError('Input error. Gamma is a loop, but first_triangle = ', first_triangle, ' is not equal',\
                ' final_triangle = ', final_triangle)
    if first_triangle is None or final_triangle is None:
        raise ValueError('Error. [first_triangle,final_triangle] = ', [first_triangle,final_triangle])
    return [first_triangle,final_triangle]

def _list_of_tau_k_and_tau_kplus1(T, crossed_arcs):
    """
    Return a list [(None,tau_1), (tau_1,tau_2), (tau_2,tau_3), ... ,(final_tau, None)]
    where tau_k is the k-th arc that is crossed by the curve

    If tau_k is a radius that is crossed in a counterclockwise direction (respectively, clockwise direction),
    then return [..., (tau_{k-1}, (tau_k,'counterclockwise)), ((tau_k, 'clockwise'), tau_{k+1}), ...]
    (respectively, return (tau_{k-1}, (tau_k,'clockwise)), ((tau_k, 'counterclockwise'),)

    Called by :func:`~sage.combinat.cluster_algebra_quiver.surface_snake_graph`
    and :func:`~sage.combinat.cluster_algebra_quiver.surface_lifted_polygon`

    EXAMPLES:

    Figure 10 and 11 of [MSW_Positivity]_::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _list_of_tau_k_and_tau_kplus1
        sage: T = ClusterTriangulation([(2,3,'b11'),(2,1,1),(4,3,'b12'),('b0',4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,'b13')], boundary_edges=['b11','b12','b13','b0']) # Counterclockwise
        sage: _list_of_tau_k_and_tau_kplus1(T._triangulation,[2,(1,'counterclockwise'),2,3,4,5,6])
        [(None, 2),
        (2, (1, 'counterclockwise')),
        ((1, 'clockwise'), 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, None)]
        sage: _list_of_tau_k_and_tau_kplus1(T._triangulation,[6,5,4,3,2,(1,'clockwise'),2])
        [(None, 6),
        (6, 5),
        (5, 4),
        (4, 3),
        (3, 2),
        (2, (1, 'clockwise')),
        ((1, 'counterclockwise'), 2),
        (2, None)]
    """
    edges = [(None, crossed_arcs[0])]
    for k in range(0,len(crossed_arcs)-1): # k from 1 to d-1
        tau_k = crossed_arcs[k]
        tau_kplus1 = crossed_arcs[k+1]
        if type(tau_k) in [list,tuple] and len(tau_k)==2: # If tau_k is a radius of a self-folded triangle
            tau_k_a_dir = edges[k][1][1]
            if tau_k_a_dir == 'clockwise':
                tau_k = (tau_k[0],'counterclockwise')
                #edges.append((tau_k, tau_kplus1))
            elif tau_k_a_dir == 'counterclockwise':
                tau_k = (tau_k[0],'clockwise')
        edges.append((tau_k, tau_kplus1)) # Gamma crosses triangle_k from tau_k to tau_kplus1
    tau_final = crossed_arcs[-1]
    edges.append((tau_final, None))
    return edges

def _list_triangles_crossed_by_curve(T, crossed_arcs, first_triangle, final_triangle, edges):
    """
    Return the list of ideal triangles Delta_k (in order) which are crossed by curve

    Called by :func:`~sage.combinat.cluster_algebra_quiver.surface_snake_graph`
    and :func:`~sage.combinat.cluster_algebra_quiver.surface_lifted_polygon`

    INPUT:

    - ``T`` -- triangulation (with edges labeled by user-labels or by
      variables)

    - ``crossed_arcs`` -- list of arcs crossed by curve (with arcs
      labeled by user-labels or by variables)

    - ``first_triangle`` -- the first triangle crossed by curve

    - ``final_triangle`` -- the final triangle crossed by curve

    - ``edges`` -- list of (triangle) edges crossed by curve,
      retrieved by
      :func:`~sage.combinat.cluster_algebra_quiver._list_of_tau_k_and_tau_kplus1`

    EXAMPLES:

    Figure 10 and 11 of [MSW_Positivity]_::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _list_triangles_crossed_by_curve, _list_of_tau_k_and_tau_kplus1
        sage: T = ClusterTriangulation([(2,3,'b11'),(2,1,1),(4,3,'b12'),('b0',4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,'b13')], boundary_edges=['b11','b12','b13','b0']) # Counterclockwise
        sage: tau_k_and_kplus1_list = _list_of_tau_k_and_tau_kplus1(T._triangulation,[2,(1,'counterclockwise'),2,3,4,5,6])
        sage: delta_k_list = _list_triangles_crossed_by_curve(T._triangulation, [2,(1,'counterclockwise'),2,3,4,5,6], (2, 3, 'b11'), (6, 7, 9), tau_k_and_kplus1_list)
        sage: delta_k_list
        [(2, 3, 'b11'),
        ((1, 'counterclockwise'), (1, 'clockwise'), 2),
        ((1, 'counterclockwise'), (1, 'clockwise'), 2),
        (2, 3, 'b11'),
        (4, 3, 'b12'),
        ('b0', 4, 5),
        (5, 6, 10),
        (6, 7, 9)]

        sage: delta_k_list.reverse()
        sage: tau_k_and_kplus1_list = _list_of_tau_k_and_tau_kplus1(T._triangulation,[6,5,4,3,2,(1,'clockwise'),2])
        sage: delta_k_list == _list_triangles_crossed_by_curve(T._triangulation, [6,5,4,3,2,(1,'clockwise'),2], (6, 7, 9), (2, 3, 'b11'), tau_k_and_kplus1_list)
        True
    """
    triangles = [first_triangle]
    for k in range(1,len(edges)-1): # Get each triangle (triangle_k) with edges tau_k and tau_{k+1}
        tau_k = edges[k][0]
        tau_kplus1 = edges[k][1]
        triangle_k = _get_triangle(T,tau_k,tau_kplus1)
        if len(triangle_k)>1: # If there are two triangles with the same two edges, compare with previous triangle
            if are_triangles_equal(triangle_k[0],triangles[k-1]):
                triangle_k = [triangle_k[1]]
        elif len(triangle_k)!=1:
            raise ValueError('Error. _get_triangle for ', tau_k, ' and ', tau_kplus1, ' returns ', triangle_k)
        triangle_k = triangle_k[0]
        triangles.append(triangle_k)
    triangles.append(final_triangle)
    return triangles

def _snake_graph(T,crossed_arcs, first_triangle=None, final_triangle=None, is_arc=True, is_loop=False, first_tile_orientation=1, boundary_edges=None):
    """
    Return descriptions of the snake/band graph for the curve crossing the arcs in input ``crossed_arcs``

    See :meth:`ClusterTriangulation.list_snake_graph`, :meth:`ClusterTriangulation.list_band_graph`,
    :meth:`ClusterTriangulation.draw_snake_graph`, and :meth:`ClusterTriangulation.draw_band_graph`

    To debug, see Mathematica function: banda

    INPUT:

    - ``T`` -- triangulation (with edges labeled by user or variables), see :meth:`ClusterTriangulation.triangulation` or :meth:`ClusterTriangulation.weighted_triangulation`
    - ``crossed_arcs`` -- A list of labels/variables corresponding to arcs that cross curve gamma, in order.
                        If curve crosses a self-folded triangle (ell,r,ell), then specify
                        ``ell, (r, 'counterclockwise), ell`` (if, as gamma is about to cross r, the puncture is to the right of gamma)
                        or ``ell, (r, clockwise), ell`` (if, as gamma is about to cross r, the puncture is to the left of gamma)
    - ``first_triangle`` -- (default:``None``) The first triangle crossed by curve
    - ``final_triangle`` -- (default:``None``) The final triangle crossed by curve (if is_loop is True, ``first_triangle`` = ``final_triangle``)
    - ``is_arc`` -- (default:True) (TODO: change this to None) Return snake graph if True
    - ``is_loop`` -- (default:False) (TODO: change this to None) Return band graph if True
    - ``first_tile_orientation`` -- (default:1) The orientation of the first tile (either +1 or -1)
    - ``boundary_edges`` -- (default:``None``) Labels or variables corresponding to boundary edges, i.e.

        *:meth:`ClusterTriangulation.boundary_edges` and
        *:meth:`ClusterTriangulation.boundary_edges_vars`

    .. NOTE::

        * 1 labels the bottom triangle of a positively-oriented tile,
        * 2 labels the top triangle of a positively-oriented tile,
        * -1 labels the bottom triangle of a positively-oriented tile,
        * -2 labels the top triangle of a negatively-oriented tile.

        * The direction (RIGHT or ABOVE) that is attached to the top triangle (labeled -2 or 2)
            indicates the location of the tile after the current tile.
            * If RIGHT, the next tile is glued to the right of the current tile
            * If ABOVE, the next tile is glued on top of the current tile

        * If this is a snake graph, then the direction for the final tile does not mean anything.
        * If this is a band graph, the direction for the final tile indicates that
            the 'cut' edge is the ceiling (if ABOVE) or to the right (if RIGHT) of the final tile.

    EXAMPLES:

    Thrice-punctured square of Figure 10 of [MSW_Positivity]_::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _snake_graph
        sage: thrice_punctured_square = [('r','r','ell'),(11,'ell',3),(3,12,4),(4,5,14),(5,6,10),(6,7,9),(8,10,9),(7,13,8)]
        sage: T = ClusterTriangulation(thrice_punctured_square, boundary_edges=[11,12,13,14])
        sage: _snake_graph(T._triangulation, ['ell', ('r','counterclockwise'), 'ell', 3, 4, 5, 6],first_tile_orientation=-1)
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
    if not(boundary_edges is None) and boundary_edges != []:
        for edge in boundary_edges:
            if edge in crossed_arcs:
                raise ValueError(edge, ' is both a boundary edge and a crossed arc.')

    end_triangles = _get_first_final_triangles(T,crossed_arcs, first_triangle, final_triangle, is_arc, is_loop)
    first_triangle = end_triangles[0]
    final_triangle = end_triangles[1]

    if is_loop:
        crossed_arcs = crossed_arcs[0:len(crossed_arcs)-1] # If loop, remove tau_final, which is equal to tau_1
    edges = _list_of_tau_k_and_tau_kplus1(T, crossed_arcs)
    triangles = _list_triangles_crossed_by_curve(T, crossed_arcs, first_triangle, final_triangle, edges)

    tau_first_a = edges[0][1] # The edge tau_k of triangle_0 (as opposed to triangle_1)

    tile_orientation = first_tile_orientation
    tile_1_bottom = (tile_orientation, _rearrange_triangle_for_snakegraph(first_triangle,tau_first_a, tile_orientation))

    out_snakegraph = [tile_1_bottom]
    triangle_kmin1 = first_triangle # The very first triangle (bottom half of the first tile)

    for k in range(1,len(edges)-1):
        tau_k_b = edges[k][0] # The edge tau_k of triangle_k (as opposed to triangle_{k-1} which also has edge tau_k)
        tau_kplus1_a = edges[k][1] #The edge tau_{k+1} of triangle_k (as opposed to triangle_{k+1} which also has edge tau_{k+1})
        #tau_kplus1_b = edges[k+1][0] # The edge tau_{k+1} of triangle_{k+1} (as opposed to triangle_k which also has edge tau_{k+1})

        rearranged_top_triangle_with_orientation = \
        _rearrange_triangle_for_snakegraph(triangles[k],tau_k_b, tile_orientation)

        # We glue the next triangle to the RIGHT or ABOVE this current tile
        if rearranged_top_triangle_with_orientation[0] == tau_kplus1_a: #actual_diagonal:
            tile_direction = RIGHT
        else:
            tile_direction = ABOVE

        top_triangle_of_tile_k = (2*tile_orientation, rearranged_top_triangle_with_orientation, tile_direction)

        # Tile k+1 will have opposite orientation as tile k
        tile_orientation = tile_orientation * (-1)

        bottom_triangle_of_tile_kplus1 = \
        (tile_orientation, _rearrange_triangle_for_snakegraph(triangles[k], tau_kplus1_a ,tile_orientation))

        # tau_kplus1_b
        out_snakegraph.extend([top_triangle_of_tile_k, bottom_triangle_of_tile_kplus1])

    #tau_d_a = edges[-2][1] # For arc, edge of the penultimate triangle (as opposed to the final triangle which also has edge tau)
    tau_d_b = edges[-1][0] # edge of the final triangle (as opposed to the penultimate triangle which also has the same edge)
    #tau_1_b = edges[-1][0] # For band

    rearranged_top_triangle_with_orientation = \
    _rearrange_triangle_for_snakegraph(final_triangle,tau_d_b, tile_orientation)

    if is_loop and k==len(edges)-2 and rearranged_top_triangle_with_orientation[0] == tau_first_a:
        tile_direction = RIGHT
        #if isinstance(tau_1_b,tuple):
        #    tile_direction = ABOVE
    else:
        tile_direction = ABOVE
        #if isinstance(tau_kplus1_a,tuple):
        #    tile_direction = RIGHT

    top_triangle_of_tile_d = (2*tile_orientation, rearranged_top_triangle_with_orientation, tile_direction)
    out_snakegraph.append(top_triangle_of_tile_d)

    return PartitionIntoTuples(out_snakegraph)

def _rearrange_triangle_for_snakegraph(triangle, diagonal, s):
    """
    Return a rearrangement of the 3-tuple input ``triangle``
    so that the ``diagonal`` is in the middle and
    the original orientation is kept (if ``s`` is +1) and reversed (if ``s`` is -1)

    To debug, see Mathematica function: rot

    INPUT:

    - ``triangle`` -- a 3-tuple triangle
    - ``diagonal`` -- one of the edges of ``triangle`` which is to go in the middle
    - ``s`` -- the orientation of the current tile (either +1 or -1)

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _rearrange_triangle_for_snakegraph
        sage: _rearrange_triangle_for_snakegraph(('a','b','c'),'b',1)
        ('a', 'b', 'c')
        sage: _rearrange_triangle_for_snakegraph((1,2,3),1,1)
        (3, 1, 2)
        sage: _rearrange_triangle_for_snakegraph((1,2,3),3,1)
        (2, 3, 1)
        sage: _rearrange_triangle_for_snakegraph(('a','b','c'),'b',-1)
        ('c', 'b', 'a')
        sage: _rearrange_triangle_for_snakegraph((1,2,3),1,-1)
        (2, 1, 3)
        sage: _rearrange_triangle_for_snakegraph((1,2,3),3,-1)
        (1, 3, 2)
    """
    x = triangle[0]
    y = triangle[1]
    z = triangle[2]

    if s == 1 and x == diagonal:
        return (z,x,y)
    if s == 1 and z == diagonal:
        return (y,z,x)
    if s == -1 and x == diagonal:
        return (y,x,z)
    if s == -1 and z == diagonal:
        return (x,z,y)
    if s == -1 and y == diagonal:
        return (z,y,x)
    return (x,y,z)

def _rearrange_triangle_small_first(triangle):
    """
    Return a rearrangement of the 3-tuple input ``triangle``
    so that the smallest edge label is in the first entry

    INPUT:

    - ``triangle`` -- a 3-tuple triangle

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _rearrange_triangle_small_first
        sage: _rearrange_triangle_small_first(('a','b','c'))
        ('a', 'b', 'c')
        sage: _rearrange_triangle_small_first((3,1,2))
        (1, 2, 3)
        sage: _rearrange_triangle_small_first((2,3,1))
        (1, 2, 3)
        sage: _rearrange_triangle_small_first((1,3,2))
        (1, 3, 2)
        sage: _rearrange_triangle_small_first((3,2,1))
        (1, 3, 2)
        sage: _rearrange_triangle_small_first((2,1,3))
        (1, 3, 2)
    """
    x, y, z = triangle[0], triangle[1], triangle[2]
    smallest = min(x,y,z)
    if smallest == x:
        return (x,y,z)
    if smallest == y:
        return (y,z,x)
    if smallest == z:
        return (z,x,y)

def _get_triangle(T, tau_k, tau_k1=None):
    """
    Return the triangle/s that share an edge with tau_k (and tau_k1 , if given).

    INPUT:

    - ``T`` -- triangulations with edges labeled by either user-labels or variables
    - ``tau_k`` -- the label or variable for an arc of T
    - ``tau_k1`` -- (default=``None``) the label of variable for the arc ``tau_{k+1}`` which is crossed after ``tau_k``

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _get_triangle
        sage: T = [(1, 4, 7), (1, 2, 5), (2, 0, 3), (0, 6, 3)]
        sage: _get_triangle(T,0)
        [(2, 0, 3), (0, 6, 3)]
        sage: _get_triangle(T,0,3)
        [(2, 0, 3), (0, 6, 3)]
        sage: _get_triangle(T,0,2)
        [(2, 0, 3)]
        sage: _get_triangle(T,0,1)
        []
    """
    triangle_k = []
    if isinstance(tau_k,list): tau_k = (tau_k[0],tau_k[1])

    if tau_k1 is not None:
        if isinstance(tau_k1,list): tau_k1 = (tau_k1[0],tau_k1[1])
        for t in T:
            if tau_k in t and tau_k1 in t:
                triangle_k.append(t)
    else:
        for t in T:
            if tau_k in t:
                triangle_k.append(t)

    return triangle_k

def try_to_find_end_triangles_for_one_crossed_arc(T, tau, first_triangle, final_triangle, is_arc, is_loop):
    """
    We assume there is only one crossed arc, or the self-folded triangle (ell, r, ell) is the only triangle crossed.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import try_to_find_end_triangles_for_one_crossed_arc
        sage: Tri = ClusterTriangulation([(2, 3, 11),(2, 1, 1),(4, 3, 12),(0, 4, 5),(5, 6, 10),(6, 7, 9),(9, 8, 10),(8, 7, 13)], boundary_edges=[11,12,13,0])
        sage: Tri.map_label_to_variable()
        {0: b10,
        1: x0,
        2: x0*x1,
        3: x2,
        4: x3,
        5: x4,
        6: x5,
        7: x6,
        8: x7,
        9: x8,
        10: x9,
        11: b11,
        12: b12,
        13: b13}

        sage: try_to_find_end_triangles_for_one_crossed_arc (Tri._weighted_triangulation,Tri.cluster()[3], None,None,True,False)
        [(x3, x2, b12), (b10, x3, x4)]

        sage: try_to_find_end_triangles_for_one_crossed_arc (Tri._triangulation,4, None,None,True,False)
        [(4, 3, 12), (0, 4, 5)]
    """
    triangle0 = _get_triangle(T, tau, None)

    if is_loop:
        raise ValueError('Input error. Gamma crossing only one arc of T means Gamma is contractible to a puncture.')

    # The case when gamma is an arc
    if len(triangle0) == 2:
        if not(first_triangle is None) and final_triangle is None:
            if are_triangles_equal(triangle0[0], first_triangle):
                final_triangle = triangle0[1]
            elif are_triangles_equal(triangle0[1], first_triangle):
                final_triangle = triangle0[0]
            else:
                raise ValueError('Input error. User gives first_triangle = ', first_triangle, ' that does not exist')
        elif first_triangle is None and not(final_triangle is None):
            if are_triangles_equal(triangle0[0], final_triangle):
                fist_triangle = triangle0[1]
            elif are_triangles_equal(triangle0[1], final_triangle):
                first_triangle = triangle0[0]
            else:
                raise ValueError('Input error. User gives final_triangle = ', final_triangle, ' that does not exist')
        elif first_triangle is None and final_triangle is None:
            first_triangle = triangle0[0]
            final_triangle = triangle0[1]

    elif len(triangle0) == 1:
        raise ValueError('Input error. Only one triangle ', triangle0[0], ' has edge tau_1=', tau)
    elif len(triangle0)==0:
        raise ValueError('Input error. No triangle has edge tau_1=', tau)
    elif len(triangle0)>2:
        raise ValueError('Input error. More than 2 triangles ', triangle0[0], ' has edge tau_1=', tau)

    return [first_triangle, final_triangle]

def try_to_find_end_triangle(T,crossed_arcs, first_or_final, is_arc, is_loop, input_triangle=None):
    """
    We assume there are at least two crossed arcs. Return the first/final triangle crossed by curve.

    See :func:`~sage.combinat.cluster_algebra_quiver.surface._get_first_final_triangles`, :func:`~sage.combinat.cluster_algebra_quiver.surface._snake_graph`

    INPUT:

    - ``T`` -- a list of triangles
    - ``crossed_arcs`` -- a list (of length at least 2) of arcs crossed by a curve
    - ``first_or_final`` -- ``First`` or ``Final``
    - ``is_arc`` -- True if curve is between marked point/s
    - ``is_loop`` -- True if curve is not between marked point/s
    - ``input_triangle`` -- (default:None) the first/final triangle specified by user


    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import try_to_find_end_triangle

    """
    if first_or_final == 'First':
        tau_1 = crossed_arcs[0]
        tau_2 = crossed_arcs[1]
    else:
        tau_1 = crossed_arcs[-1]
        tau_2 = crossed_arcs[-2]

    triangle1 = _get_triangle(T, tau_1, tau_2)
    N = len(triangle1)

    if N == 1:
        triangle1_A = triangle1[0]
        triangle0 = _get_triangle(T, tau_1, None)
        if len(triangle0) == 1:
            raise ValueError('Incorrect input. Only one triangle ', triangle0,' has edge ', tau_1)
        elif len(triangle0)==2:
            triangle0.remove(triangle1_A)
            if is_loop and first_or_final == 'Final': # Maybe erase
                triangle0 = [triangle1_A] #Maybe erase
            if not(input_triangle is None):
                if are_triangles_equal(input_triangle,triangle0[0]):
                    return input_triangle
                else:
                    raise ValueError('Incorrect input. User inputs an end triangle ', first_or_final, ' triangle: ', input_triangle, ' but it should be ', triangle0[0])
            return triangle0[0]
        elif len(triangle0) > 2 :
            raise ValueError('Incorrect input. There are more than 2 triangles with edge =', tau_1)

    elif N == 2:
        if not(input_triangle is None):
            triangle0 = _get_triangle(T, tau_1, None)
            if are_triangles_equal(input_triangle,triangle0[0]) or are_triangles_equal(input_triangle,triangle0[1]):
                return input_triangle
            else:
                raise ValueError('Incorrect input. User inputs an end triangle ', first_or_final, ' triangle: ', input_triangle, ' which does not exist.')
        if are_triangles_equal(triangle1[0],triangle1[1]) and len(T)==2: # If T is a triangulation of a once-punctured torus
            return triangle1[0]
        if first_or_final == 'First':
            raise ValueError('In this case user must specify the first triangle = triangle_0 crossed by gamma. '\
            ,'since there are two triangles ', triangle1, ' with edges ', tau_1,' and ', tau_2)
        else:
            raise ValueError('In this case user must specify the final triangle = triangle_d crossed by gamma. '\
            ,'since there are two triangles ', triangle1, ' with edges ', tau_1,' and ', tau_2)

    elif N == 0:
        raise ValueError('Incorrect input. There are no triangle with edges ', tau_1, ' and ', tau_2)

    elif N > 3:
        raise ValueError('Incorrect input. There are more than 2 triangles ', triangle1 \
        ,  ' with edges ', tau_1, ' and ', tau_2)

#########################################################################################
################### ENDS: FUNCTIONS FOR CONSTRUCTING BAND/SNAKE GRAPH ###################
#########################################################################################

##################################
### BEGIN: DRAWING SNAKE GRAPH ###
##################################

def _draw_matching(perfect_matching, matching_weight=None, pos=None, xy=(0,0), white_space=1, print_user_labels=False, symmetric_difference=None):
    """
    Returns drawing of a single snake graph with perfect matchings and symmetric difference.

    INPUT:

    - ``perfect_matching`` -- perfect_matching
    - ``matching_weight`` -- matching_weight
    - ``pos`` -- pos
    - ``xy`` -- current position's Cartesian coordinate (x,y)
    - ``white_space`` -- white_space
    - ``print_user_labels`` -- whether to print user labels or to print variables x_i for labeling the edges of the perfect matching
    - ``symmetric_difference`` -- symmetric difference of the minimal matching and perfect_matching

    Matching of tile (a,b,c,d) = (floor, right, ceiling, left) where an entry
    is equal to 1 if it is part of the matching and 0 otherwise

    .. TODO::

        Currently the parameter ``print_user_labels`` is always ``False``.
        We may want to enable user to print out perfect matchings of snake/band graphs with edges labeled by user-labels from
        :meth:`ClusterTriangulation.arcs` and :meth:`ClusterTriangulation.boundary_edges`
        (instead of variables from :meth:`ClusterTriangulation.cluster` and :meth:`ClusterTriangulation.boundary_edges_vars`)

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _draw_matching, GetMonomialTerm, GetMinimalMatching, FlipAllFlippableTiles
        sage: T = ClusterTriangulation([(1,2,'b4'),(1,0,'b5'),(0,3,'b6'),(2,3,'b7')], boundary_edges=['b4','b5','b6','b7'])
        sage: c = T.cluster()
        sage: G = T.list_band_graph([c[1],c[2],c[3],c[0],c[1]], user_labels=False)
        sage: min_pm = GetMinimalMatching(G)
        sage: min_pm
        [['minimal PM'],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 1), 'RIGHT'],
        [(1, 0, 1, 0), 'RIGHT'],
        [(0, 1, 0, 0), 'ABOVE']]]
        sage: PM = FlipAllFlippableTiles(min_pm)[0]
        sage: monomial_term = GetMonomialTerm(G, PM, boundary_edges=T._boundary_edges_vars, is_loop=False)
        sage: PM
        [[2],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'ABOVE']]]
        sage: _draw_matching(PM, monomial_term, print_user_labels=False)
        (Graphics object consisting of 8 graphics primitives, (4, 0))
    """
    from sage.plot.graphics import Graphics
    from sage.plot.line import line
    from sage.plot.text import text

    drawing = Graphics()
    PM_color = 'black'
    PM_thickness = 6
    SD_color = 'red'
    SD_thickness = 4

    (x,y)=xy

    if not(matching_weight is None):
        if print_user_labels:
            drawing = drawing + text('$'+str(matching_weight)+'$', (x+0.5*white_space,y-0.6), rgbcolor = 'gray')
        else:
            drawing = drawing + text('$'+str(matching_weight).replace('x','x_{').replace('y','y_{').replace('b','b_{').replace('*','}').replace('/','}/')+'}$', (x+0.5*white_space,y-0.6), rgbcolor = 'black')
    if not(pos is None):
        pos_str = '$(' + str(pos) + ')$. '
        drawing = drawing + text(pos_str, (x ,y-0.2), rgbcolor=(0,1,0.2))

    for PM in perfect_matching[1]:
        if PM[0][0] == 1: # floor
            drawing = drawing + line([(x+0,y+0),(x+1,y+0)], rgbcolor=PM_color, thickness=PM_thickness)
        elif PM[0][1] == 1: # right
            drawing = drawing + line([(x+1,y+0),(x+1,y+1)], rgbcolor=PM_color, thickness=PM_thickness)
        if PM[0][2] == 1: # ceiling
            drawing = drawing + line([(x+0,y+1),(x+1,y+1)], rgbcolor=PM_color, thickness=PM_thickness)
        elif PM[0][3] == 1: # left
            drawing = drawing + line([(x+0,y+0),(x+0,y+1)], rgbcolor=PM_color, thickness=PM_thickness)
        DIR = PM[1]
        if DIR == ABOVE:
            y=y+1
        else:
            x=x+1

    if symmetric_difference:
        (x,y)=xy
        for SD in symmetric_difference[1]:
            if SD[0][0] == 1: # floor
                drawing = drawing + line([(x+0,y+0),(x+1,y+0)], rgbcolor=SD_color, thickness=SD_thickness)
            if SD[0][1] == 1: # right
                drawing = drawing + line([(x+1,y+0),(x+1,y+1)], rgbcolor=SD_color, thickness=SD_thickness)
            if SD[0][2] == 1: # ceiling
                drawing = drawing + line([(x+0,y+1),(x+1,y+1)], rgbcolor=SD_color, thickness=SD_thickness)
            if SD[0][3] == 1: # left
                drawing = drawing + line([(x+0,y+0),(x+0,y+1)], rgbcolor=SD_color, thickness=SD_thickness)

            DIR = SD[1]
            if DIR == ABOVE:
                y=y+1
            else:
                x=x+1

    return drawing, (x+ 2*white_space,0)


def _draw_snake_graph(G, print_user_labels, xy=(0, 0)):
    """
    Return the plot of the snake graph `G`.

    INPUT:

    - ``G`` -- snake graph description, see :meth:`~sage.combinat.cluster_algebra_quiver.ClusterTriangulation.list_snake_graph`
    - ``user_labels`` -- whether the labels of the snake graph are user-given labels or variables (x_i and b_i) corresponding to arcs/boundary edges
    - ``xy`` -- (default:(0,0)) snake graph should be plotted at xy=(a,b)

    EXAMPLES:

    Figure 10 of [MSW_Positivity]_::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _draw_snake_graph
        sage: thrice_punctured_square = [('r','r','ell'),(11,'ell',3),(3,12,4),(4,5,14),(5,6,10),(6,7,9),(8,10,9),(7,13,8)]
        sage: T = ClusterTriangulation(thrice_punctured_square, boundary_edges=[11,12,13,14])
        sage: G_user_labels = T.list_snake_graph(['ell', ('r','counterclockwise'), 'ell', 3, 4, 5, 6],first_tile_orientation=-1,user_labels=True)
        sage: _draw_snake_graph(G_user_labels,print_user_labels=True)
        Graphics object consisting of 43 graphics primitives
    """
    from sage.plot.graphics import Graphics
    from sage.plot.line import line
    from sage.plot.text import text

    drawing = Graphics()
    x, y = 0,0
    (x,y)=xy

    for pos in range(0,len(G)):
        tile = G[pos]

        #print tile

        tile_drawing = line([(x+1,y+0),(x+0,y+0),(x+0,y+1),(x+1,y+1),(x+1,y+0),(x+0,y+1)])
        floor = tile[0][1][0]
        if type(floor) in [tuple, list]: floor=floor[0]
        floor = str(floor)
        if not print_user_labels: # a noose represent a product of two cluster variables (parallel arcs with different notchings at a puncture)
            floor = '$' + floor.replace('*','}').replace('x','x_{').replace('y','y_{').replace('b','b_{') + '}$'

        diagonal = tile[0][1][1]
        if type(diagonal) in [tuple, list]: diagonal=diagonal[0]
        diagonal = str(diagonal)
        if not print_user_labels:
            diagonal = '$' + diagonal.replace('*','}').replace('x','x_{').replace('y','y_{').replace('b','b_{') + '}$'

        left_side = tile[0][1][2]
        if type(left_side) in [tuple, list]: left_side=left_side[0]
        left_side = str(left_side)
        if not print_user_labels:
            left_side = '$' + left_side.replace('*','}').replace('x','x_{').replace('y','y_{').replace('b','b_{') + '}$'

        right_side = tile[1][1][2]
        if type(right_side) in [tuple, list]: right_side=right_side[0]
        right_side = str(right_side)
        if not print_user_labels:
            right_side = '$' + right_side.replace('*','}').replace('x','x_{').replace('y','y_{').replace('b','b_{') + '}$'

        ceiling = tile[1][1][0]
        if type(ceiling) in [tuple, list]: ceiling=ceiling[0]
        ceiling = str(ceiling)
        if not print_user_labels:
            ceiling = '$' + ceiling.replace('*','}').replace('x','x_{').replace('y','y_{').replace('b','b_{') + '}$'

        orientation = tile[0][0]
        if orientation == 1: orientation='$+$'
        else: orientation='$-$'

        text_color = (1,0,0) # red
        labels = text(diagonal,(x+0.5,y+0.5),vertical_alignment='bottom', rgbcolor=text_color)\
        + text(right_side,(x+1,y+0.5),horizontal_alignment='left', rgbcolor=text_color)\
        + text(ceiling,(x+0.5,y+1),vertical_alignment='bottom', rgbcolor=text_color) \
        + text(orientation,(x+0.8, y+0.8))

        if pos>0:
            PREVIOUS_DIR = G[pos-1][1][-1]
            if PREVIOUS_DIR == RIGHT: # Then draw the label of the bottom edge
                labels = labels + text(floor,(x+0.5,y+0),vertical_alignment='bottom', rgbcolor=text_color)
            elif PREVIOUS_DIR == ABOVE: # Then draw the label of the left edge
                labels = labels + text(left_side,(x+0,y+0.5),horizontal_alignment='left', rgbcolor=text_color)
        else: # For the first tile, draw labels for both bottom and left edges
            labels = labels + \
            text(floor,(x+0.5,y+0),vertical_alignment='bottom', rgbcolor=text_color)\
            + text(left_side,(x+0,y+0.5),horizontal_alignment='left', rgbcolor=text_color)

        DIR = tile[1][-1]
        if DIR == ABOVE:
            y=y+1
        else:
            x=x+1

        drawing = drawing + tile_drawing + labels
        drawing.axes(False)
        drawing.set_aspect_ratio(1)

    return drawing


#########################################################################################
################### BEGIN: FUNCTIONS FOR COMPUTING MINIMAL MATCHING #####################
#########################################################################################

def GetMinimalMatching(G):
    """
    Return the minimal matching of the snake/band graph ``G``.

    To debug, see Mathematica function: MachingInicial[listaDirecciones]

    INPUT:

    - ``G`` -- a band/snake graph

    EXAMPLES:

    The minimal matching of a graph with one tile is such that the
    floor and ceiling edges of the only tile are marked (1,0,1,0).

    In addition, when a tile is the final tile, we write DIR='ABOVE'
    direction by default::

        sage: from sage.combinat.cluster_algebra_quiver.surface import GetMinimalMatching,_snake_graph
        sage: T = ClusterTriangulation([(0,2,1),(0,4,3),(1,6,5)])
        sage: G_varlabel = _snake_graph(T._weighted_triangulation,[T.cluster()[0]],is_arc=True)
        sage: GetMinimalMatching(G_varlabel)
        [['minimal PM'], [[(1, 0, 1, 0), 'ABOVE']]]
        sage: G_userlabel = _snake_graph(T._triangulation,[0],is_arc=True)
        sage: GetMinimalMatching(G_userlabel)
        [['minimal PM'], [[(1, 0, 1, 0), 'ABOVE']]]

    Figure 10 of [MSW_Positivity]_::

        sage: thrice_punctured_square = [('r','r','ell'),(11,'ell',3),(3,12,4),(4,5,14),(5,6,10),(6,7,9),(8,10,9),(7,13,8)]
        sage: T = ClusterTriangulation(thrice_punctured_square, boundary_edges=[11,12,13,14])
        sage: G = T.list_snake_graph([5,6,7,8,9,6,5], first_tile_orientation=1, user_labels=True)
        sage: GetMinimalMatching(G)
        [['minimal PM'],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 1), 'RIGHT'],
        [(1, 0, 1, 0), 'RIGHT'],
        [(0, 1, 0, 0), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE']]]

        sage: G_opp = T.list_snake_graph([5,6,7,8,9,6,5], first_tile_orientation=-1, user_labels=True)
        sage: GetMinimalMatching(G_opp)
        [['minimal PM'],
        [[(0, 0, 0, 1), 'RIGHT'],
        [(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'RIGHT'],
        [(0, 0, 0, 0), 'RIGHT'],
        [(1, 0, 1, 0), 'RIGHT'],
        [(0, 1, 0, 0), 'ABOVE']]]
    """
    # The list of (for each tile in the band/snake graph, the
    # direction of the next tile)
    graph_directions = snake_graph_tile_directions(G)
    #print 'graph_directions = snake_graph_tile_directions(G) and G is :', G

    if len(graph_directions) == 1:  # If curve crosses the triangulation once
        return [['minimal PM'], [[(1, 0, 1, 0), graph_directions[0]]]]

    #print 'graph_directions: ', graph_directions
    orientation_of_first_tile = G[0][0][0]
    initial_matching = [[_minimal_matching_first_tile(graph_directions[0],orientation_of_first_tile), graph_directions[0]]]

    # Continue assigning minimal matching for the rest of the tiles, except the final tile
    for pos in range(1,len(graph_directions)-1):
       # The edges that we have marked so far (from the previous tile)
       last_marking_in_list = initial_matching[-1][0]
       current_tile_mark = _minimal_matching_current_tile(graph_directions[pos-1],graph_directions[pos], last_marking_in_list)
       initial_matching.append([current_tile_mark, graph_directions[pos]])

    last_marking_in_list = initial_matching[-1][0]
    penultimate_direction = graph_directions[-2]

    final_tile = _minimal_matching_final_tile(penultimate_direction, last_marking_in_list)
    initial_matching.append([final_tile, graph_directions[-1]])
    initial_matching =  [['minimal PM'],initial_matching]
    return initial_matching

def snake_graph_tile_directions(G):
    """
    Return [the positions (RIGHT or ABOVE) of all tiles in the band/snake graph including the last tile]

    INPUT:

    - ``G`` -- a snake/band graph

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import snake_graph_tile_directions
        sage: G = [[(1, ('b2', 1, 0)), (2, ('b1', 1, 2), 'ABOVE')],\
        [(-1, ('b1', 2, 1)), (-2, (3, 2, 'b4'), 'RIGHT')],\
        [(1, (2, 3, 'b4')), (2, (0, 3, 'b3'), 'RIGHT')],\
        [(-1, (3, 0, 'b3')), (-2, ('b2', 0, 1), 'ABOVE')]]
        sage: snake_graph_tile_directions(G)
        ['ABOVE', 'RIGHT', 'RIGHT', 'ABOVE']
    """
    directions = []
    for pos in range(0,len(G)):
        direction = G[pos][1][2]
        directions.append(direction)
    return directions

def _minimal_matching_first_tile(DIR, orientation_of_first_tile):
    """
    Return the minimal matching of the first tile of the band/snake graph

    To debug, see Mathematica function: funAuxA

    INPUT:

    - ``DIR`` -- ABOVE (resp, RIGHT) if the second tile is above (resp, to the right) the first tile

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _minimal_matching_first_tile
        sage: _minimal_matching_first_tile('ABOVE',1)
        (1, 0, 0, 0)
        sage: _minimal_matching_first_tile('RIGHT',1)
        (1, 0, 1, 0)
        sage: _minimal_matching_first_tile('ABOVE',-1)
        (0, 1, 0, 1)
        sage: _minimal_matching_first_tile('RIGHT',-1)
        (0, 0, 0, 1)
    """
    if orientation_of_first_tile == 1:
        if DIR == ABOVE:
            mark = (1,0,0,0)
            return mark
        elif DIR == RIGHT:
            mark = (1,0,1,0)
            return mark
    elif orientation_of_first_tile == -1:
        if DIR == ABOVE:
            mark = (0,1,0,1)
            return mark
        elif DIR == RIGHT:
            mark = (0,0,0,1)
            return mark
    raise ValueError("Bug. DIR should be ABOVE or RIGHT. orientation_of_first_tile should be 1 or -1")

def _minimal_matching_current_tile(previous_DIR, current_DIR, last_marking_in_list):
    """
    INPUT:
    [previous_DIR, current_DIR]
    - ``last_marking_in_list`` -> the last edges that we have marked so far on the previous tile

    Returns the marking for the current tile (to achieve the minimal matching).
    We assume the current tile is not the final tile in the band/snake graph

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _minimal_matching_current_tile
        sage: _minimal_matching_current_tile('ABOVE','RIGHT',(0, 1, 0, 0))
        (0, 0, 1, 0)
    """
     # case: current tile is above previous tile, and next tile is above the current tile
    if [previous_DIR,current_DIR] == [ABOVE, ABOVE]:
        if last_marking_in_list == (1, 0, 0, 0):
            return (0, 1, 0, 1)
        if last_marking_in_list == (0, 1, 0, 1):
            return (0, 0, 0, 0)
        if last_marking_in_list == (0, 0, 0, 0):
            return (0, 1, 0, 1)
        if last_marking_in_list == (0, 1, 0, 0):
            return (0, 0, 0, 0)

    # case: current tile is above previous tile, and next tile is to the right of the current tile
    if [previous_DIR,current_DIR] == [ABOVE, RIGHT]:
        if last_marking_in_list == (1, 0, 0, 0):
            return (0, 0, 0, 1)
        if last_marking_in_list == (0, 1, 0, 1):
            return (0, 0, 1, 0)
        if last_marking_in_list == (0, 1, 0, 0):
            return (0, 0, 1, 0)
        if last_marking_in_list == (0, 0, 0, 0):
            return (0, 0, 0, 1)

    # case: current tile is to the right of previous tile, and next tile is to the right of the current tile
    if [previous_DIR,current_DIR] == [RIGHT, RIGHT]:
        if last_marking_in_list == (0, 0, 0, 1):
            return (1, 0, 1, 0)
        if last_marking_in_list == (0, 0, 1, 0):
            return (0, 0, 0, 0)
        if last_marking_in_list == (1, 0, 1, 0):
            return (0, 0, 0, 0)
        if last_marking_in_list == (0, 0, 0, 0):
            return (1, 0, 1, 0)

    # case: current tile is to the right of previous tile, and next tile is above the current tile
    if [previous_DIR,current_DIR] == [RIGHT, ABOVE]:
        if last_marking_in_list == (0, 0, 0, 1):
            return (1, 0, 0, 0)
        if last_marking_in_list == (1, 0, 1, 0):
            return (0, 1, 0, 0)
        if last_marking_in_list == (0, 0, 1, 0):
            return (0, 1, 0, 0)
        if last_marking_in_list == (0, 0, 0, 0):
            return (1, 0, 0, 0)

def _minimal_matching_final_tile(penultimate_direction, penultimate_tile_mark):
    """
    Return the matching for the final tile

    To debug, see Mathematica function: funAuxFinDEBanda

    INPUT:

    - ``penultimate_direction`` -- whether the final tile is to the RIGHT or ABOVE the penultimate tile,
    - ``penultimate_tile_mark`` -- the edges that we have marked on the penultimate tile

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _minimal_matching_final_tile
        sage: _minimal_matching_final_tile('ABOVE',(0,1,0,1)) == (0, 0, 1, 0)
        True
    """
    if penultimate_direction == ABOVE:
        if penultimate_tile_mark == (0,1,0,1):
            return (0,0,1,0)
        if penultimate_tile_mark == (0,0,0,0):
            return (0,1,0,1)
        if penultimate_tile_mark == (0,1,0,0):
            return (0,0,1,0)
        if penultimate_tile_mark == (1,0,0,0):
            return (0,1,0,1)

    if penultimate_direction == RIGHT:
        if penultimate_tile_mark == (0,0,0,1):
            return (1,0,1,0)
        if penultimate_tile_mark == (1,0,1,0):
            return (0,1,0,0)
        if penultimate_tile_mark == (0,0,0,0):
            return (1,0,1,0)
        if penultimate_tile_mark == (0,0,1,0):
            return (0,1,0,0)

#########################################################################################
################### ENDS: FUNCTIONS FOR COMPUTING MINIMAL MATCHING #####################
#########################################################################################

##########################################################################################
### BEGINS: FUNCTIONS THAT CREATE ALL PERFECT MATCHINGS FROM THE MINIMAL MATCHING ########
##########################################################################################

def FlipAllFlippableTilesInList(input_list_matchings):
    """
    For each perfect matching ``PM`` in ``input_list_matchings``,
    create more perfect matchings by flipping every flippable tile.

    Return the list ``list_new_matchings`` of perfect matchings
    created by this method.

    The list ``list_new_matchings`` does not purposely contain any
    perfect matching info from the input ``input_list_matchings``, but
    ``list_new_matchings`` is likely to most perfect matchings from
    ``input_list_matchings``.

    See :func:`~sage.combinat.cluster_algebra_quiver.surface.GetAllMatchings`

    To debug, see Mathematica function: confiSelectivoSobrelista

    INPUT:

    - ``input_list_matchings`` -- list of matchings of a band/snake graph

    EXAMPLES:

    Figure 10 of [MSW_Positivity]_::

        sage: from sage.combinat.cluster_algebra_quiver.surface import FlipAllFlippableTilesInList, FlipAllFlippableTiles, GetMinimalMatching
        sage: thrice_punctured_square = [('r','r','ell'),(11,'ell',3),(3,12,4),(4,5,14),(5,6,10),(6,7,9),(8,10,9),(7,13,8)]
        sage: T = ClusterTriangulation(thrice_punctured_square, boundary_edges=[11,12,13,14])
        sage: G = T.list_snake_graph([5,6,7,8,9,6,5], first_tile_orientation=1, user_labels=True)
        sage: current_matchings = FlipAllFlippableTilesInList([GetMinimalMatching(G)])
        sage: current_matchings
        [[[2],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE']]],
        [[5],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 1), 'RIGHT'],
        [(1, 0, 1, 0), 'RIGHT'],
        [(0, 1, 0, 0), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE']]]]
        sage: FlipAllFlippableTilesInList(current_matchings)
        [[[1],
        [[(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'RIGHT'],
        [(0, 1, 0, 0), 'RIGHT'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE']]],
        [[3],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 0, 0, 1), 'RIGHT'],
        [(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE']]],
        [[5],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE']]],
        [[2],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE']]],
        [[6],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 1), 'RIGHT'],
        [(1, 0, 1, 0), 'RIGHT'],
        [(0, 1, 0, 0), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE'],
        [(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE']]]]

    """
    list_new_matchings = []
    for matching in input_list_matchings:
        list_new_matchings = UniqueList(list_new_matchings + FlipAllFlippableTiles(matching))

    return list_new_matchings


def FlipAllFlippableTiles(input_tiles):
    """
    Return a list ``out_list_new_matchings`` of new perfect matchings
    produced by flipping a flippable tile of ``input_tiles``.

    The output list does not include the input perfect matching ``input_tiles``.

    See :func:`~sage.combinat.cluster_algebra_quiver.surface.FlipAllFlippableTilesInList`

    To debug, see Mathematica function: confiSelectivo

    .. NOTE::

        Record the indices of the tiles that we flip, i.e.
        if we flip ``[[indices],[perfect matching PM]]`` at tile ``k``, we produce [[k],[PM flipped at k]].
        If a tile indexed by ``k`` is flippable but ``k`` is in ``input_tiles[0]``, don't flip tile ``k``.

    INPUT:

    - ``input_tiles`` -- [[the indices of tiles that were last flipped],[a perfect matching of a snake/band graph]] if ``input_tiles`` is not the minimal matching. [['minimal PM'], [the minimal matching]] if ``input_tiles`` is the minimal matching.

    EXAMPLES:

    Figure 10 of [MSW_Positivity]_::

        sage: from sage.combinat.cluster_algebra_quiver.surface import FlipAllFlippableTiles, GetMinimalMatching
        sage: thrice_punctured_square = [('r','r','ell'),(11,'ell',3),(3,12,4),(4,5,14),(5,6,10),(6,7,9),(8,10,9),(7,13,8)]
        sage: T = ClusterTriangulation(thrice_punctured_square, boundary_edges=[11,12,13,14])
        sage: G = T.list_snake_graph([5,6,7,8,9,6,5], first_tile_orientation=1, user_labels=True)
        sage: MinPM = GetMinimalMatching(G)
        sage: MinPM
        [['minimal PM'],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 1), 'RIGHT'],
        [(1, 0, 1, 0), 'RIGHT'],
        [(0, 1, 0, 0), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE']]]
        sage: FlipAllFlippableTiles(MinPM)
        [[[2],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE']]],
        [[5],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 0, 0, 1), 'RIGHT'],
        [(1, 0, 1, 0), 'RIGHT'],
        [(0, 1, 0, 0), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE']]]]
    """
    last_flipped_tiles = input_tiles[0] # Tiles that we will not flip this time
    tiles = input_tiles[1]

    # The list of tile markings [(a, b, c, d), (e, f, g, h), ...], without RIGHT or ABOVE information
    list_tile_marks_without_directions = [tile[0] for tile in tiles]

    #horizonal_tiles_indices = [positions of list_tile_marks_without_directions that match [1,0,1,0]]
    if (1,0,1,0) in list_tile_marks_without_directions:
        horizonal_indices = [list_tile_marks_without_directions.index((1,0,1,0))]
        for pos in range(horizonal_indices[0]+1,len(list_tile_marks_without_directions)):
            if list_tile_marks_without_directions[pos]==(1,0,1,0):
                horizonal_indices.append(pos)
    else:
        horizonal_indices = []

    #vertical_tiles_indices = [positions of list_tile_marks_without_directions that match (0,1,0,1)]
    if (0,1,0,1) in list_tile_marks_without_directions:
        vertical_indices = [list_tile_marks_without_directions.index((0,1,0,1))]
        for pos in range(vertical_indices[0]+1,len(list_tile_marks_without_directions)):
            if list_tile_marks_without_directions[pos]==(0,1,0,1):
                vertical_indices.append(pos)
    else:
        vertical_indices = []

    flippable_tiles_indices = horizonal_indices + vertical_indices
    flippable_tiles_indices.sort()
    out_list_new_matchings = []

    for j in range(0,len(flippable_tiles_indices)):
        current_tile_pos = flippable_tiles_indices[j]

        if current_tile_pos not in last_flipped_tiles:
            # then we flip this tile
            current_tile = tiles[current_tile_pos] # [[a,b,c,d], ABOVE/RIGHT]
            flipped_current_tile = FlipTile(current_tile)

            # new_matching = Replace position current_tile_pos of "tiles" with flipped_current_tile
            new_matching = list(tiles)
            new_matching[current_tile_pos] = flipped_current_tile

            # We may need to change the marks of edges that the current tile shares with the tiles before or after
            if current_tile_pos < len(tiles)-1:
                next_tile = tiles[current_tile_pos + 1]
                next_tile_direction = next_tile[1]
                NextTileMark = next_tile[0]
                NextTileMark = GetNextTileMarking(current_tile[0],current_tile[1], NextTileMark)
                flipped_next_tile = [NextTileMark, next_tile_direction]
                new_matching[current_tile_pos+1] = flipped_next_tile

            if current_tile_pos > 0:
                previous_tile = tiles[current_tile_pos-1]
                previous_tile_direction = previous_tile [1]
                PreviousTileMark = previous_tile[0]
                PreviousTileMark = GetPreviousTileMarking(previous_tile_direction, current_tile[0], PreviousTileMark)
                flipped_previous_tile =[PreviousTileMark, previous_tile_direction]
                new_matching[current_tile_pos-1] = flipped_previous_tile

            new_matching_with_last_flipped_tile = [[current_tile_pos], new_matching]
            out_list_new_matchings.append(new_matching_with_last_flipped_tile)
    return out_list_new_matchings

def GetMoreLastFlippedTilesInList(list_matchings, list_other_matchings):
    """
    Return a (possibly non-unique) list of perfect matchings of a snake/band graph
    (each perfect matching P includes info on all possible indices of tiles that were flipped to get to P)

    To debug, see Mathematica function: comparoNuevosContraVarios

    INPUT:

    - ``list_matchings`` -- a (possibly with repeats) list of perfect matchings of a snake/band graph
    - ``list_other_matchings`` -- a list of perfect matchings of a snake/band graph which is to be compared with ``list_matchings``

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import GetMoreLastFlippedTilesInList
        sage: list_matchings = [[[3], [[(0, 1, 0, 1), 'ABOVE'], [(0, 0, 1, 0), 'RIGHT'], [(0, 0, 0, 0), 'RIGHT'], [(1, 0, 1, 0), 'ABOVE'], [(1, 0, 1, 0), 'ABOVE']]], \
        [[0], [[(0, 1, 0, 1), 'ABOVE'], [(0, 0, 1, 0), 'RIGHT'], [(0, 0, 0, 0), 'RIGHT'], [(1, 0, 1, 0), 'ABOVE'], [(1, 0, 1, 0), 'ABOVE']]], \
        [[4], [[(1, 0, 1, 0), 'ABOVE'], [(1, 0, 1, 0), 'RIGHT'], [(0, 0, 0, 0), 'RIGHT'], [(1, 0, 0, 0), 'ABOVE'], [(0, 1, 0, 1), 'ABOVE']]], \
        [[1], [[(1, 0, 1, 0), 'ABOVE'], [(1, 0, 1, 0), 'RIGHT'], [(0, 0, 0, 0), 'RIGHT'], [(1, 0, 0, 0), 'ABOVE'], [(0, 1, 0, 1), 'ABOVE']]]]

        sage: list_other_matchings = [[['minimal PM'], [[(1, 0, 0, 0), 'ABOVE'], [(0, 0, 0, 1), 'RIGHT'], [(1, 0, 1, 0), 'RIGHT'], [(0, 1, 0, 0), 'ABOVE'], [(0, 0, 1, 0), 'ABOVE']]], \
        [[2], [[(1, 0, 0, 0), 'ABOVE'], [(0, 1, 0, 1), 'RIGHT'], [(0, 1, 0, 1), 'RIGHT'], [(0, 1, 0, 1), 'ABOVE'], [(0, 0, 1, 0), 'ABOVE']]], \
        [[1], [[(1, 0, 1, 0), 'ABOVE'], [(1, 0, 1, 0), 'RIGHT'], [(0, 1, 0, 0), 'RIGHT'], [(0, 1, 0, 1), 'ABOVE'], [(0, 0, 1, 0), 'ABOVE']]], \
        [[3], [[(1, 0, 0, 0), 'ABOVE'], [(0, 1, 0, 1), 'RIGHT'], [(0, 0, 0, 1), 'RIGHT'], [(1, 0, 1, 0), 'ABOVE'], [(1, 0, 1, 0), 'ABOVE']]], \
        [[0], [[(0, 1, 0, 1), 'ABOVE'], [(0, 0, 1, 0), 'RIGHT'], [(0, 1, 0, 0), 'RIGHT'], [(0, 1, 0, 1), 'ABOVE'], [(0, 0, 1, 0), 'ABOVE']]], \
        [[1, 3], [[(1, 0, 1, 0), 'ABOVE'], [(1, 0, 1, 0), 'RIGHT'], [(0, 0, 0, 0), 'RIGHT'], [(1, 0, 1, 0), 'ABOVE'], [(1, 0, 1, 0), 'ABOVE']]], \
        [[4], [[(1, 0, 0, 0), 'ABOVE'], [(0, 1, 0, 1), 'RIGHT'], [(0, 0, 0, 1), 'RIGHT'], [(1, 0, 0, 0), 'ABOVE'], [(0, 1, 0, 1), 'ABOVE']]], \
        [[3], [[(0, 1, 0, 1), 'ABOVE'], [(0, 0, 1, 0), 'RIGHT'], [(0, 0, 0, 0), 'RIGHT'], [(1, 0, 1, 0), 'ABOVE'], [(1, 0, 1, 0), 'ABOVE']]], \
        [[0], [[(0, 1, 0, 1), 'ABOVE'], [(0, 0, 1, 0), 'RIGHT'], [(0, 0, 0, 0), 'RIGHT'], [(1, 0, 1, 0), 'ABOVE'], [(1, 0, 1, 0), 'ABOVE']]], \
        [[4], [[(1, 0, 1, 0), 'ABOVE'], [(1, 0, 1, 0), 'RIGHT'], [(0, 0, 0, 0), 'RIGHT'], [(1, 0, 0, 0), 'ABOVE'], [(0, 1, 0, 1), 'ABOVE']]], \
        [[1], [[(1, 0, 1, 0), 'ABOVE'], [(1, 0, 1, 0), 'RIGHT'], [(0, 0, 0, 0), 'RIGHT'], [(1, 0, 0, 0), 'ABOVE'], [(0, 1, 0, 1), 'ABOVE']]]]

        sage: GetMoreLastFlippedTilesInList(list_matchings, list_other_matchings)
        [[[0, 3],
        [[(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'RIGHT'],
        [(0, 0, 0, 0), 'RIGHT'],
        [(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE']]],
        [[0, 3],
        [[(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'RIGHT'],
        [(0, 0, 0, 0), 'RIGHT'],
        [(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'ABOVE']]],
        [[1, 4],
        [[(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'RIGHT'],
        [(0, 0, 0, 0), 'RIGHT'],
        [(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE']]],
        [[1, 4],
        [[(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'RIGHT'],
        [(0, 0, 0, 0), 'RIGHT'],
        [(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE']]]]
    """
    out_new_matchings_corrected_indices = []
    for outer_ct in range(0,len(list_matchings)):
        current_matching = list_matchings[outer_ct]
        for inner_ct in range(0,len(list_other_matchings)):
            temp_current_matching = \
            GetMoreLastFlippedTiles(current_matching, list_other_matchings[inner_ct])
            current_matching = temp_current_matching
        out_new_matchings_corrected_indices.append(current_matching)

    return out_new_matchings_corrected_indices

def GetMoreLastFlippedTiles(current_matching, to_compare_matching):
    """
    Return input ``current_matching`` if input ``to_compare_matching`` contains a different perfect matching, and
    return [[indices of last flipped tiles which may include new indices if ``to_compare_matching`` is the same matching],
    [current_matching[1]]]

    To debug, see Mathematica function: comparo1A1

    INPUT:

    - ``current_matching`` -- [[last flipped tiles of current matching], [current matching]]
    - ``to_compare_matching`` -- another matching that may or may not be equal to current_matching

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import GetMoreLastFlippedTiles
        sage: PM = [[1],\
        [[(1, 0, 1, 0), 'ABOVE'],\
        [(1, 0, 1, 0), 'RIGHT'],\
        [(0, 1, 0, 0), 'RIGHT'],\
        [(0, 1, 0, 1), 'ABOVE'],\
        [(0, 0, 0, 0), 'ABOVE'],\
        [(0, 1, 0, 1), 'ABOVE'],\
        [(0, 0, 1, 0), 'ABOVE']]]
        sage: PM_to_compare = [[0,3], PM[1]]
        sage: GetMoreLastFlippedTiles(PM, PM_to_compare)
        [[0, 1, 3],
        [[(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 1, 0), 'RIGHT'],
        [(0, 1, 0, 0), 'RIGHT'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE']]]

        sage: PM = [[3],\
        [[(1, 0, 0, 0), 'ABOVE'],\
        [(0, 1, 0, 1), 'RIGHT'],\
        [(0, 0, 0, 1), 'RIGHT'],\
        [(1, 0, 1, 0), 'ABOVE'],\
        [(1, 0, 0, 0), 'ABOVE'],\
        [(0, 1, 0, 1), 'ABOVE'],\
        [(0, 0, 1, 0), 'ABOVE']]]
        sage: PM_to_compare = [[5],\
        [[(1, 0, 0, 0), 'ABOVE'],\
        [(0, 1, 0, 1), 'RIGHT'],\
        [(0, 1, 0, 1), 'RIGHT'],\
        [(0, 1, 0, 1), 'ABOVE'],\
        [(0, 0, 1, 0), 'ABOVE'],\
        [(1, 0, 1, 0), 'ABOVE'],\
        [(1, 0, 1, 0), 'ABOVE']]]
        sage: PM[1] == PM_to_compare[1]
        False
        sage: GetMoreLastFlippedTiles(PM, PM_to_compare)
        [[3],
        [[(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'RIGHT'],
        [(0, 0, 0, 1), 'RIGHT'],
        [(1, 0, 1, 0), 'ABOVE'],
        [(1, 0, 0, 0), 'ABOVE'],
        [(0, 1, 0, 1), 'ABOVE'],
        [(0, 0, 1, 0), 'ABOVE']]]
    """
    last_flipped_tiles = current_matching[0]
    if current_matching[1] == to_compare_matching[1]: # If the two input matchings have the same marks {a,b,c,d}
        last_flipped_tiles = UniqueList(current_matching[0] + to_compare_matching[0])
        last_flipped_tiles.sort()
    return [last_flipped_tiles, current_matching[1]]

def GetNextTileMarking(tile,DIR, NextTileMark):
    """
    (1, 0, 1, 0), RIGHT, (a, b, c, d) -> return (a, b, c, 1)

    (1, 0, 1, 0), ABOVE, (a, b, c, d) -> return (0, b, c, d)

    (0, 1, 0, 1), RIGHT, (a, b, c, d) -> return (a, b, c, 0)

    (0, 1, 0, 1), ABOVE, (a, b, c, d) -> return (1, b, c, d)

    INPUT:

    - ``tile`` -- tile marking (either (1,0,1,0) or (0,1,0,1))

    - ``DIR`` -- whether the next tile is glued to the ``RIGHT`` or
      ``ABOVE`` the current tile

    - ``NextTileMark`` -- the marks of the next tile

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import GetNextTileMarking
        sage: GetNextTileMarking((1, 0, 1, 0), 'RIGHT', (0, 1, 0, 0))
        (0, 1, 0, 1)
        sage: GetNextTileMarking((0, 1, 0, 1), 'ABOVE', (0, 0, 0, 0))
        (1, 0, 0, 0)
        sage: GetNextTileMarking((0, 1, 0, 1), 'RIGHT', (0, 1, 0, 1))
        (0, 1, 0, 0)
        sage: GetNextTileMarking((0, 1, 0, 1), 'ABOVE', (0, 0, 1, 0))
        (1, 0, 1, 0)
        sage: GetNextTileMarking((1, 0, 1, 0), 'ABOVE', (1, 0, 1, 0))
        (0, 0, 1, 0)
    """
    a = NextTileMark[0]
    b = NextTileMark[1]
    c = NextTileMark[2]
    d = NextTileMark[3]
    if [tile,DIR] == [(1,0,1,0), RIGHT]:
        return (a, b, c, 1)
    if [tile,DIR] == [(1,0,1,0), ABOVE]:
        return (0, b, c, d)
    if [tile,DIR] == [(0, 1, 0, 1), RIGHT]:
        return (a, b, c, 0)
    if [tile,DIR] == [(0, 1, 0, 1), ABOVE]:
        return (1, b, c, d)

def GetPreviousTileMarking(DIR, tile, PreviousTileMark):
    """
    RIGHT, (1, 0, 1, 0), (a, b,c, d) -> return (a, 1, c, d)

    ABOVE, (1, 0, 1, 0), (a, b,c, d) -> return (0, b, c, d)

    RIGHT, (0, 1, 0, 1), (a, b,c, d) -> return (a, 0, c, d)

    ABOVE, (0, 1, 0, 1), (a, b,c, d) -> return (a, b, 1, d)

    INPUT:

    - ``DIR`` -- whether the next tile is glued to the ``RIGHT`` or
      ``ABOVE`` the current tile

    - ``tile`` -- tile marking (either (1,0,1,0) or (0,1,0,1))

    - ``PreviousTileMark`` -- the marks of the previous tile

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import GetPreviousTileMarking
        sage: GetPreviousTileMarking('RIGHT', (1, 0, 1, 0), (0, 0, 0, 1))
        (0, 1, 0, 1)
        sage: GetPreviousTileMarking('RIGHT', (0, 1, 0, 1), (0, 1, 0, 1))
        (0, 0, 0, 1)
        sage: GetPreviousTileMarking('ABOVE', (0, 1, 0, 1), (0, 0, 0, 0))
        (0, 0, 1, 0)
        sage: GetPreviousTileMarking('ABOVE', (0, 1, 0, 1), (1, 0, 0, 0))
        (1, 0, 1, 0)
        sage: GetPreviousTileMarking('ABOVE', (1, 0, 1, 0), (1, 0, 1, 0))
        (1, 0, 0, 0)
    """
    a = PreviousTileMark[0]
    b = PreviousTileMark[1]
    c = PreviousTileMark[2]
    d = PreviousTileMark[3]
    if [DIR, tile] == [RIGHT, (1,0,1,0)]:
        return (a, 1, c, d)
    if [DIR, tile] == [ABOVE, (1,0,1,0)]:
        return (a, b, 0, d)
    if [DIR, tile] == [RIGHT, (0,1,0,1)]:
        return (a, 0, c, d)
    if [DIR, tile] == [ABOVE, (0,1,0,1)]:
        return (a, b, 1, d)

def FlipTile(tile_info):
    """
    [(1, 0, 1, 0),DIR] -> return [(0, 1, 0, 1), DIR]

    [(0, 1, 0, 1), DIR] -> return [(1, 0, 1, 0), DIR]

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.surface import FlipTile
        sage: FlipTile([(1,0,1,0),'ABOVE'])
        [(0, 1, 0, 1), 'ABOVE']
        sage: FlipTile([(0, 1, 0, 1), 'ABOVE'])
        [(1, 0, 1, 0), 'ABOVE']
    """
    mark = tile_info[0]
    DIR = tile_info[1]
    if mark == (1, 0, 1, 0):
        return [(0, 1, 0, 1), DIR]
    if mark == (0, 1, 0, 1):
        return [(1, 0, 1, 0), DIR]

##########################################################################################
####### ENDS: FUNCTIONS THAT CREATE ALL PERFECT MATCHINGS FROM THE MINIMAL MATCHING ######
##########################################################################################

######################################################################
######### BEGIN: FUNCTIONS THAT DRAW LIFTED POLYGONS #################
######################################################################

def _triangle_to_draw(triangle, triangle_type, tau, tau_placement, test_k=None):
    """
    Returns the follow information (dictionary) on how the triangle should be displayed and glued to the next triangle::

        {'bottom'/'top': edge label or variable info,
        'left': edge label or variable info,
        'right': edge label or variable info,
        'glued_on': 'glued_to_the_left'/'glued_to_the_right'/'glued_to_the_top'/'glued_to_the_bottom'}

    If 'glued_on': 'glued_to_the_bottom', the current input (upside-down) triangle ``triangle`` is glued BELOW the previous (pyramid) triangle.
    If 'glued_on': 'glued_to_the_top', the current input (pyramid) triangle ``triangle`` is glued ABOVE the previous (upside-down) triangle.
    If 'glued_on': 'glued_to_the_left', this current input triangle ``triangle`` is glued to the LEFT of the previous triangle.
    If 'glued_on': 'glued_to_the_right', this current input triangle ``triangle`` is glued to the RIGHT of the previous triangle.

    If 'bottom' is returned, the input triangle ``triangle`` is a pyramid.
    If 'top' is returned, the input triangle ``triangle`` is an upside-down pyramid.

    INPUT:

    - ``triangle`` -- a 3-tuple of user-given labels or variables corresponding to triangle edges

    - ``tau`` -- the first edge of ``triangle`` crossed by curve gamma
        (i.e. if ``triangle`` is Delta_k, then tau is tau_k)

    - ``tau_placement`` -- (options:``bottom``,``top``,``left``,``right``)

        *  If ``tau_placement`` = ``bottom``, this current input ``triangle`` to be glued below the previous triangle.
        *  If ``tau_placement`` = ``top``, this current input ``triangle`` to be glued above the previous triangle.
        *  If ``tau_placement`` = ``left``, this current input ``triangle`` to be glued to the left the previous triangle.
        *  If ``tau_placement`` = ``right``, this current input ``triangle`` to be glued to the right the previous triangle.

    - ``test_k`` -- include the current index of the triangle (for debugging purposes)

    EXAMPLES:

    Figure 10 and 11 of [MSW_Positivity]_::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _triangle_to_draw
        sage: from sage.combinat.cluster_algebra_quiver.surface import _list_triangles_crossed_by_curve, _list_of_tau_k_and_tau_kplus1
        sage: T = ClusterTriangulation([(2,3,'b11'),(2,1,1),(4,3,'b12'),('b0',4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,'b13')], boundary_edges=['b11','b12','b13','b0']) # Counterclockwise
        sage: tau_k_and_kplus1_list = _list_of_tau_k_and_tau_kplus1(T._triangulation,[2,(1,'counterclockwise'),2,3,4,5,6])
        sage: delta_k_list = _list_triangles_crossed_by_curve(T._triangulation, [2,(1,'counterclockwise'),2,3,4,5,6], (2, 3, 'b11'), (6, 7, 9), tau_k_and_kplus1_list)
        sage: delta_k_list
        [(2, 3, 'b11'),
        ((1, 'counterclockwise'), (1, 'clockwise'), 2),
        ((1, 'counterclockwise'), (1, 'clockwise'), 2),
        (2, 3, 'b11'),
        (4, 3, 'b12'),
        ('b0', 4, 5),
        (5, 6, 10),
        (6, 7, 9)]
        sage: first_triangle_to_draw = _triangle_to_draw(delta_k_list[0], 'pyramid', tau_k_and_kplus1_list[0][1], 'right')
        sage: first_triangle_to_draw
        {'bottom': 'b11', 'glued_on': 'glued_to_the_left', 'left': 3, 'right': 2}
        sage: _triangle_to_draw(delta_k_list[2], 'upsidedown triangle', tau_k_and_kplus1_list[2][0], 'top', test_k=None)
        {'glued_on': 'glued_to_the_bottom',
        'left': 2,
        'right': (1, 'counterclockwise'),
        'top': (1, 'clockwise')}
        sage: _triangle_to_draw(delta_k_list[5], 'pyramid', tau_k_and_kplus1_list[5][0], 'left', test_k=None)
        {'bottom': 5, 'glued_on': 'glued_to_the_right', 'left': 4, 'right': 'b0'}
    """
    alpha, beta = None, None
    if triangle[0] == tau: # tau_k of the triangle_k (with edges tau_k and tau_{k+1})
        alpha = triangle[1]
        beta = triangle[2]
    elif triangle[1] == tau:
        alpha = triangle[2]
        beta = triangle[0]
    elif triangle[2] == tau:
        alpha = triangle[0]
        beta = triangle[1]
    else:
        raise ValueError ('Error!! tau=', tau, ' but triangle: ', triangle, 'test_k: ', test_k)

    if triangle_type == 'pyramid':
        if tau_placement == 'bottom':
            return dict(bottom=tau, right=alpha, left=beta, glued_on='glued_to_the_top') # this (pyramid) triangle is glued above the previous (upside down) triangle
        elif tau_placement == 'right':
            return dict(right=tau, left=alpha, bottom=beta, glued_on='glued_to_the_left') # this (pyramid) triangle is glued to the left of the previous (upside down) triangle
        elif tau_placement == 'left':
            return dict(left=tau, bottom=alpha, right=beta, glued_on='glued_to_the_right') # this (pyramid) triangle is glued to the right of the previous (upside down) triangle

    elif triangle_type == 'upsidedown triangle':
        if tau_placement == 'top':
            return dict(top=tau, left=alpha, right=beta, glued_on='glued_to_the_bottom') # this (upside down) triangle is glued to the bottom of the previous (pyramid) triangle
        elif tau_placement == 'left':
            return dict(left=tau, right=alpha, top=beta, glued_on='glued_to_the_right') # this (upside down) triangle is glued to the right of the previous (pyramid) triangle
        elif tau_placement  == 'right':
            return dict(right=tau, top=alpha, left=beta, glued_on='glued_to_the_left') # this (upside down) triangle is glued to the left of the previous (pyramid) triangle
    else:
        raise ValueError ('Input error: ', triangle, triangle_type, tau, tau_placement)

def _lifted_polygon(T, crossed_arcs, first_triangle, final_triangle ,is_arc, is_loop):
    """Return list of copies of the triangles which are crossed by gamma.
    See Section 7 of [MSW_Positivity]_

    INPUT:

    - ``T`` -- list of triangles, see :meth:`ClusterTriangulation.weighted_triangulation` and :meth:`ClusterTriangulation.triangulation`
    - ``crossed_arcs`` -- the list of variables or user-given labels labeling the arcs crossed by curve
    - ``first_triangle`` -- a 3-tuple with labels of the first triangle crossed by curve
    - ``final_triangle`` -- a 3-tuple with labels of the final triangle crossed by curve
    - ``is_arc`` -- True if curve is between marked point/s
    - ``is_loop`` -- True if curve is a loop in the interior of surface

    EXAMPLES:

    Figure 10 and 11 of [MSW_Positivity]_, where the boundary edges
    are labeled 11,12,13,0, and the arcs are labeled 1 (a radius), 2
    (a noose), and 3,4,5,...,10::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _lifted_polygon
        sage: T = ClusterTriangulation([(2,3,11),(2,1,1),(4,3,12),(0,4,5),(5,6,10),(6,7,9),(9,8,10),(8,7,13)], boundary_edges=[11,12,13,0])
        sage: c = T.cluster()
        sage: r = c[0]
        sage: ell = c[1] * r

        sage: _lifted_polygon (T._weighted_triangulation, [ell,(r,'counterclockwise'),ell, c[3-1],c[4-1],c[5-1],c[6-1]], None, None, True, False)
        [{'bottom': b11, 'glued_on': 'glued_to_the_left', 'left': x2, 'right': x0*x1},
        {'glued_on': 'glued_to_the_right',
        'left': x0*x1,
        'right': (x0, 'counterclockwise'),
        'top': (x0, 'clockwise')},
        {'bottom': (x0, 'clockwise'),
        'glued_on': 'glued_to_the_top',
        'left': (x0, 'counterclockwise'),
        'right': x0*x1},
        {'glued_on': 'glued_to_the_right', 'left': x0*x1, 'right': x2, 'top': b11},
        {'bottom': b12, 'glued_on': 'glued_to_the_right', 'left': x2, 'right': x3},
        {'glued_on': 'glued_to_the_right', 'left': x3, 'right': x4, 'top': b10},
        {'bottom': x5, 'glued_on': 'glued_to_the_right', 'left': x4, 'right': x9},
        {'glued_on': 'glued_to_the_bottom', 'left': x6, 'right': x8, 'top': x5}]

        sage: _lifted_polygon (T._triangulation, [2,(1,'counterclockwise'),2,3,4,5,6], None, None, True, False)
        [{'bottom': 11, 'glued_on': 'glued_to_the_left', 'left': 3, 'right': 2},
        {'glued_on': 'glued_to_the_right',
        'left': 2,
        'right': (1, 'counterclockwise'),
        'top': (1, 'clockwise')},
        {'bottom': (1, 'clockwise'),
        'glued_on': 'glued_to_the_top',
        'left': (1, 'counterclockwise'),
        'right': 2},
        {'glued_on': 'glued_to_the_right', 'left': 2, 'right': 3, 'top': 11},
        {'bottom': 12, 'glued_on': 'glued_to_the_right', 'left': 3, 'right': 4},
         {'glued_on': 'glued_to_the_right', 'left': 4, 'right': 5, 'top': 0},
        {'bottom': 6, 'glued_on': 'glued_to_the_right', 'left': 5, 'right': 10},
        {'glued_on': 'glued_to_the_bottom', 'left': 7, 'right': 9, 'top': 6}]

        sage: crossed_vars = [c[2],ell,(r,'counterclockwise'),ell, c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[5],c[4],c[3],c[2]]
        sage: _lifted_polygon(T._weighted_triangulation, crossed_vars, first_triangle=None, final_triangle=None, is_arc=False, is_loop=True)
        [{'bottom': x3, 'glued_on': 'glued_to_the_left', 'left': b12, 'right': x2},
        {'glued_on': 'glued_to_the_right', 'left': x2, 'right': b11, 'top': x0*x1},
        {'bottom': x0*x1,
        'glued_on': 'glued_to_the_top',
        'left': (x0, 'clockwise'),
        'right': (x0, 'counterclockwise')},
        {'glued_on': 'glued_to_the_right',
        'left': (x0, 'clockwise'),
        'right': x0*x1,
        'top': (x0, 'counterclockwise')},
        {'bottom': x2, 'glued_on': 'glued_to_the_right', 'left': x0*x1, 'right': b11},
        {'glued_on': 'glued_to_the_bottom', 'left': b12, 'right': x3, 'top': x2},
        {'bottom': x4, 'glued_on': 'glued_to_the_right', 'left': x3, 'right': b10},
        {'glued_on': 'glued_to_the_bottom', 'left': x5, 'right': x9, 'top': x4},
        {'bottom': x8, 'glued_on': 'glued_to_the_left', 'left': x6, 'right': x5},
        {'glued_on': 'glued_to_the_left', 'left': x7, 'right': x6, 'top': b13},
        {'bottom': x8, 'glued_on': 'glued_to_the_left', 'left': x9, 'right': x7},
        {'glued_on': 'glued_to_the_bottom', 'left': x5, 'right': x6, 'top': x8},
        {'bottom': x4, 'glued_on': 'glued_to_the_left', 'left': x9, 'right': x5},
        {'glued_on': 'glued_to_the_bottom', 'left': b10, 'right': x3, 'top': x4},
        {'bottom': x2, 'glued_on': 'glued_to_the_right', 'left': x3, 'right': b12},
        {'glued_on': 'glued_to_the_bottom', 'left': b12, 'right': x3, 'top': x2}]

        sage: crossed_userlabels = [3,2,(1,'counterclockwise'),2, 3,4,5,6,7,8,9,6,5,4,3]
        sage: _lifted_polygon(T._triangulation, crossed_userlabels, first_triangle=None, final_triangle=None, is_arc=False, is_loop=True)
        [{'bottom': 4, 'glued_on': 'glued_to_the_left', 'left': 12, 'right': 3},
        {'glued_on': 'glued_to_the_right', 'left': 3, 'right': 11, 'top': 2},
        {'bottom': 2,
        'glued_on': 'glued_to_the_top',
        'left': (1, 'clockwise'),
        'right': (1, 'counterclockwise')},
        {'glued_on': 'glued_to_the_right',
        'left': (1, 'clockwise'),
        'right': 2,
        'top': (1, 'counterclockwise')},
        {'bottom': 3, 'glued_on': 'glued_to_the_right', 'left': 2, 'right': 11},
        {'glued_on': 'glued_to_the_bottom', 'left': 12, 'right': 4, 'top': 3},
        {'bottom': 5, 'glued_on': 'glued_to_the_right', 'left': 4, 'right': 0},
        {'glued_on': 'glued_to_the_bottom', 'left': 6, 'right': 10, 'top': 5},
        {'bottom': 9, 'glued_on': 'glued_to_the_left', 'left': 7, 'right': 6},
        {'glued_on': 'glued_to_the_left', 'left': 8, 'right': 7, 'top': 13},
        {'bottom': 9, 'glued_on': 'glued_to_the_left', 'left': 10, 'right': 8},
        {'glued_on': 'glued_to_the_bottom', 'left': 6, 'right': 7, 'top': 9},
        {'bottom': 5, 'glued_on': 'glued_to_the_left', 'left': 10, 'right': 6},
        {'glued_on': 'glued_to_the_bottom', 'left': 0, 'right': 4, 'top': 5},
        {'bottom': 3, 'glued_on': 'glued_to_the_right', 'left': 4, 'right': 12},
        {'glued_on': 'glued_to_the_bottom', 'left': 12, 'right': 4, 'top': 3}]

    Figure 8 of [MSW_Bases]_ with triangulation having arcs 1, ..., 6,
    and gamma crosses 1,2,3,4,1 (1 is listed twice); and the outer
    boundary edges, clockwise starting from starting point of
    gamma are 7,8,9, 10; and the inner boundary edges, clockwise
    starting from ending point of gamma:11, 0::

        sage: T = ClusterTriangulation([(8,7,5),(5,4,1),(1,0,2),(2,11,3),(3,4,6),(6,10,9)], boundary_edges=[7,8,9,10,11,0])
        sage: c = [item for item in T.cluster()]

        sage: _lifted_polygon(T._weighted_triangulation, [c[1-1],c[2-1],c[3-1],c[4-1],c[1-1]],None, None, True, False)
        [{'bottom': x3, 'glued_on': 'glued_to_the_left', 'left': x4, 'right': x0},
        {'glued_on': 'glued_to_the_right', 'left': x0, 'right': b6, 'top': x1},
        {'bottom': x1, 'glued_on': 'glued_to_the_top', 'left': x2, 'right': b11},
        {'glued_on': 'glued_to_the_left', 'left': x5, 'right': x2, 'top': x3},
        {'bottom': x3, 'glued_on': 'glued_to_the_top', 'left': x4, 'right': x0},
        {'glued_on': 'glued_to_the_right', 'left': x0, 'right': b6, 'top': x1}]

        sage: _lifted_polygon(T._triangulation, [1,2,3,4,1], None, None, True, False)
        [{'bottom': 4, 'glued_on': 'glued_to_the_left', 'left': 5, 'right': 1},
        {'glued_on': 'glued_to_the_right', 'left': 1, 'right': 0, 'top': 2},
        {'bottom': 2, 'glued_on': 'glued_to_the_top', 'left': 3, 'right': 11},
        {'glued_on': 'glued_to_the_left', 'left': 6, 'right': 3, 'top': 4},
        {'bottom': 4, 'glued_on': 'glued_to_the_top', 'left': 5, 'right': 1},
        {'glued_on': 'glued_to_the_right', 'left': 1, 'right': 0, 'top': 2}]

        sage: _lifted_polygon(T._weighted_triangulation, [c[1-1],c[2-1],c[3-1],c[4-1],c[1-1]], None, None, False, True)
        [{'bottom': x3, 'glued_on': 'glued_to_the_left', 'left': x4, 'right': x0},
        {'glued_on': 'glued_to_the_right', 'left': x0, 'right': b6, 'top': x1},
        {'bottom': x1, 'glued_on': 'glued_to_the_top', 'left': x2, 'right': b11},
        {'glued_on': 'glued_to_the_left', 'left': x5, 'right': x2, 'top': x3},
        {'bottom': x3, 'glued_on': 'glued_to_the_top', 'left': x4, 'right': x0},
        {'glued_on': 'glued_to_the_right', 'left': x0, 'right': x4, 'top': x3}]

        sage: _lifted_polygon(T._triangulation, [1,2,3,4,1], None, None, False, True)
        [{'bottom': 4, 'glued_on': 'glued_to_the_left', 'left': 5, 'right': 1},
        {'glued_on': 'glued_to_the_right', 'left': 1, 'right': 0, 'top': 2},
        {'bottom': 2, 'glued_on': 'glued_to_the_top', 'left': 3, 'right': 11},
        {'glued_on': 'glued_to_the_left', 'left': 6, 'right': 3, 'top': 4},
        {'bottom': 4, 'glued_on': 'glued_to_the_top', 'left': 5, 'right': 1},
        {'glued_on': 'glued_to_the_right', 'left': 1, 'right': 5, 'top': 4}]

    .. TODO::

        Can the output be more readable?
    """
    end_triangles = _get_first_final_triangles(T,crossed_arcs, first_triangle,
                                               final_triangle, is_arc, is_loop)
    first_triangle = end_triangles[0]
    final_triangle = end_triangles[1]

    edges = _list_of_tau_k_and_tau_kplus1 (T, crossed_arcs)
    triangles = _list_triangles_crossed_by_curve (T, crossed_arcs, first_triangle, final_triangle, edges)
    triangles_to_draw = [_triangle_to_draw(triangles[0], 'pyramid', edges[0][1], 'right')]

    for k in range(1,len(triangles)):
        tau = edges[k][0]
        placement = None
        other_side_of_radius = None
        triangle_kmin1 = triangles_to_draw[k-1]

        if type(tau) in [tuple,list]: # If tau is a radius, then switch clockwise/ counterclockwise
            other_side_of_radius = edges[k-1][1]

        if k % 2 == 1:
            if triangle_kmin1 ['bottom'] in [tau, other_side_of_radius]:
                placement = 'top'
            elif triangle_kmin1 ['right'] in [tau, other_side_of_radius]:
                placement = 'left'
            elif triangle_kmin1 ['left'] in [tau, other_side_of_radius]:
                placement = 'right'
            else:
                raise ValueError ('Bug? odd k= ', k ,' tau_k = ', tau, 'triangle_kmin1: ', triangle_kmin1)
            triangles_to_draw.append(_triangle_to_draw(triangles[k], 'upsidedown triangle', tau, placement, test_k=k))

        elif k % 2 == 0:
            if triangle_kmin1 ['top'] in [tau, other_side_of_radius]:
                placement = 'bottom'
            elif triangle_kmin1 ['left'] in [tau, other_side_of_radius]:
                placement = 'right'
            elif triangle_kmin1 ['right'] in [tau, other_side_of_radius]:
                placement = 'left'
            else:
                raise ValueError ('Bug? even k= ', k, ' tau_k = ', tau, 'triangle_kmin1: ', triangle_kmin1)
            triangles_to_draw.append(_triangle_to_draw(triangles[k], 'pyramid', tau, placement, test_k=k))
    return triangles_to_draw

def _draw_lifted_curve(lifted_polygon, is_arc, is_loop):
    """
    Returns the graphics of a lifted triangulated polygon/annulus and the lifted curve.
    See :meth:`ClusterTriangulation.draw_lifted_arc` and :meth:`ClusterTriangulation.draw_lifted_loop`

    EXAMPLES:

    Figure 8 of [MSW_Bases]_ with triangulation having arcs 1, ..., 6,
    and gamma crosses 1,2,3,4,1 (1 is listed twice); and the outer
    boundary edges, clockwise starting from starting point of gamma
    are 7,8,9, 10; and the inner boundary edges, clockwise starting
    from ending point of gamma:11, 0::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _draw_lifted_curve, _lifted_polygon
        sage: T = ClusterTriangulation([(8,7,5),(5,4,1),(1,0,2),(2,11,3),(3,4,6),(6,10,9)], boundary_edges=[7,8,9,10,11,0])
        sage: lifted_polygon = _lifted_polygon(T._triangulation, [1,2,3,4,1], None, None, is_arc=False, is_loop=True)
        sage: _draw_lifted_curve(lifted_polygon, is_arc=False, is_loop=True)
        Graphics object consisting of 40 graphics primitives
    """
    from sage.plot.graphics import Graphics
    from sage.plot.line import line
    from sage.plot.arrow import arrow2d
    #from sage.plot.text import text
    drawing = Graphics()
    x, y = 0,0
    k=0
    curve_style = 'dashed'
    curve_thickness = 2
    drawing = drawing + _draw_triangle (lifted_polygon[0],x,y, 0)

    if is_arc:
        #drawing = drawing + line([(x-2,y-2),(x+1,y-1)],linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # from left corner to right edge of pyramid
        drawing = drawing + arrow2d((x-2,y-2),(x+1,y-1),linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # from left corner to right edge of pyramid
    elif is_loop:
        drawing = drawing + line([(x,y-1),(x+1,y-1)],linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # from center of pyramid to right edge of pyramid

    triangle_count = len(lifted_polygon)
    if is_loop:
        triangle_count = triangle_count - 1 # If a loop, don't draw the last triangle in lifted_polygon

    for k in range(1,triangle_count):
        triangle = lifted_polygon [k]
        if k<len(lifted_polygon)-1:
            next_triangle = lifted_polygon[k+1]

        if k % 2 == 1: # draw an upside-down triangle
            if triangle['glued_on'] == 'glued_to_the_bottom':
                # glue to the bottom of the previous (pyramid) triangle_{k-1}
                y = y-4 # new x,y is at the bottom of the upside-down triangle
                if k<len(lifted_polygon)-1:
                    if next_triangle['glued_on'] == 'glued_to_the_right':
                        drawing = drawing + arrow2d((x,y+2),(x+1,y+1), linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # from top to right edge
                    if next_triangle['glued_on'] == 'glued_to_the_left':
                        drawing = drawing + arrow2d((x,y+2),(x-1,y+1), linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # from top to left edge
                else:
                    drawing = drawing + arrow2d((x,y+2),(x,y), linestyle = curve_style, rgbcolor='red') # from top to bottom corner

            elif triangle['glued_on'] == 'glued_to_the_right':
                # glue to the right of the previous (pyramid) triangle_{k-1}
                x = x+2
                y = y-2 # new x,y is at the bottom of the upside-down triangle
                if k<len(lifted_polygon)-1:
                    if next_triangle['glued_on'] == 'glued_to_the_right':
                        drawing = drawing + arrow2d((x-1,y+1),(x+1,y+1),linestyle = curve_style,  rgbcolor='red')#, thickness=curve_thickness) # from left edge to right edge
                    elif next_triangle['glued_on'] == 'glued_to_the_top':
                        drawing = drawing + arrow2d((x-1,y+1),(x,y+2), linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # from left edge to top edge
                else:
                    drawing = drawing + arrow2d((x-1,y+1),(x+2,y+2), linestyle = curve_style, rgbcolor='red') # from left edge to corner

            elif triangle['glued_on'] == 'glued_to_the_left':
                # glue to the left of the previous (pyramid) triangle_{k-1}
                x = x-2
                y = y-2 # new x,y is at the bottom of the upside-down triangle
                if k<len(lifted_polygon)-1:
                    if next_triangle['glued_on'] == 'glued_to_the_top':
                        drawing = drawing + arrow2d((x+1,y+1),(x,y+2), linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # from right side to top edge
                    elif next_triangle['glued_on'] == 'glued_to_the_left':
                        drawing = drawing + arrow2d((x+1,y+1),(x-1,y+1), linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # from right side to left
                else:
                    drawing = drawing + arrow2d((x+1,y+1),(x-2,y+2), linestyle = curve_style, rgbcolor='red') # from right side to corner

        elif k % 2 == 0: # draw a pyramid
            if triangle['glued_on'] == 'glued_to_the_top':
                # glue to the top of the previous (upside-down) triangle_{k-1}
                y = y+4 # new x,y is at the bottom of the upside-down triangle
                if k<len(lifted_polygon)-1:
                    if next_triangle['glued_on'] == 'glued_to_the_left':
                        drawing = drawing + arrow2d((x,y-2),(x-1,y-1), linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # from bottom edge to left edge
                    elif next_triangle['glued_on'] == 'glued_to_the_right':
                        drawing = drawing + arrow2d((x,y-2),(x+1,y-1), linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # from bottom edge to right corner
                else:
                    drawing = drawing + arrow2d((x,y-2),(x,y),linestyle = curve_style, rgbcolor='red') # from bottom to top corner

            if triangle['glued_on'] == 'glued_to_the_left': #### TODO CHECK
                # glue to the left of the previous (upside) triangle_{k-1}
                x = x-2
                y = y+2 # new x,y is at the bottom of the upside-down triangle
                if k<len(lifted_polygon)-1:
                    if next_triangle['glued_on'] == 'glued_to_the_left':
                        drawing = drawing + arrow2d((x+1,y-1),(x-1,y-1), linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # From right edge to left edge
                    elif next_triangle['glued_on'] == 'glued_to_the_bottom':
                        drawing = drawing + arrow2d((x+1,y-1),(x,y-2), linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # From right edge to bottom edge
                else:
                    drawing = drawing + arrow2d((x+1,y-1),(x-2,y-2), linestyle = curve_style, rgbcolor='red') # From right edge to corner
            elif triangle['glued_on'] == 'glued_to_the_right':
                # glue to the right of the previous (upside) triangle_{k-1}
                x = x+2
                y = y+2 # new x,y is at the bottom of the upside-down triangle
                if k<len(lifted_polygon)-1:
                    if next_triangle['glued_on'] == 'glued_to_the_right':
                        drawing = drawing + arrow2d((x-1,y-1),(x+1,y-1), linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # From left edge to right edge
                    if next_triangle['glued_on'] == 'glued_to_the_bottom':
                        drawing = drawing + arrow2d((x-1,y-1),(x,y-2), linestyle = curve_style, rgbcolor='red')#, thickness=curve_thickness) # From left edge to bottom edge
                else:
                    drawing = drawing + arrow2d((x-1,y-1),(x+2,y-2), linestyle = curve_style, rgbcolor='red') # From left edge to corner

        drawing = drawing + _draw_triangle (triangle,x,y, k)

    drawing.axes(False)
    drawing.set_aspect_ratio(1)
    return drawing

def _draw_triangle(triangle,x,y, k):
    """
    Returns a triangle starting at (x,y) in the shape of a pyramid (resp, an upside-down triangle) if k is even (resp, k is odd).
    See :func:`~sage.combinat.cluster_algebra_quiver.surface._draw_lifted_curve`,
    :meth:`ClusterTriangulation.draw_lifted_arc` and :meth:`ClusterTriangulation.draw_lifted_loop`

    INPUT:

    - ``triangle`` -- a 3-tuple (or list of length 3) of labels of a triangle
    - ``x,y`` -- (x,y) indicates where triangle should be drawn
    - ``k`` -- a natural number which indexes the triangles

    EXAMPLES:

    Figure 8 of [MSW_Bases]_ with triangulation having arcs 1, ..., 6,
    and gamma crosses 1,2,3,4,1 (1 is listed twice); and the outer
    boundary edges, clockwise starting from starting point of gamma
    are 7,8,9,1 and 0; and the inner boundary edges, clockwise
    starting from ending point of gamma are 11 and 0::

        sage: from sage.combinat.cluster_algebra_quiver.surface import _draw_triangle
        sage: from sage.combinat.cluster_algebra_quiver.surface import _draw_lifted_curve, _lifted_polygon
        sage: T = ClusterTriangulation([(8,7,5),(5,4,1),(1,0,2),(2,11,3),(3,4,6),(6,10,9)], boundary_edges=[7,8,9,10,11,0])
        sage: lifted_polygon = _lifted_polygon(T._triangulation, [1,2,3,4,1], None, None, is_arc=False, is_loop=True)
        sage: _draw_lifted_curve(lifted_polygon, is_arc=False, is_loop=True)
        Graphics object consisting of 40 graphics primitives
    """
    from sage.plot.graphics import Graphics
    from sage.plot.line import line
    from sage.plot.text import text
    lines = Graphics()
    texts = Graphics()

    text_color = 'green'

    for side in triangle:
        if type(triangle[side]) in [tuple, list]: # If the side is (radius, 'clockwise/counterclockwise)
            triangle[side] = triangle[side][0] # assign as simple radius

    if k % 2 == 0:
        lines = lines + line ([(x,y),(x-2,y-2)]) # left side of pyramid
        lines = lines + line ([(x-2,y-2),(x+2,y-2)]) # bottom of pyramid
        lines = lines + line ([(x,y),(x+2,y-2)]) # right side of pyramid

        left = str(triangle['left'])
        if not isinstance(triangle['left'], str):
            left = '$' + left.replace('*','}').replace('x','x_{').replace('b','b_{') + '}$'
        texts = texts + text (left,(x-1,y-1), horizontal_alignment='left', rgbcolor=text_color) # left side of pyramid

        bottom = str(triangle['bottom'])
        if not isinstance(triangle['bottom'], str):
            bottom = '$' + bottom.replace('*','}').replace('x','x_{').replace('b','b_{') + '}$'
        texts = texts + text (bottom,(x,y-2), vertical_alignment='bottom', rgbcolor=text_color) # bottom of pyramid

        right = str(triangle['right'])
        if not isinstance(triangle['right'], str):
            right = '$' + right.replace('*','}').replace('x','x_{').replace('b','b_{') + '}$'
        texts = texts + text (right, (x+1,y-1), horizontal_alignment='left', rgbcolor=text_color) # right side of pyramid

        texts = texts + text ('$\\Delta_{'+str(k)+'}$', (x+0.2, y-0.7), rgbcolor=(0,1,0))

    elif k % 2 == 1:
        lines = lines + line ([(x,y),(x-2,y+2)]) # left side of the upside-down triangle
        lines = lines + line ([(x-2,y+2),(x+2,y+2)]) # top of the upside-down triangle
        lines = lines + line ([(x,y),(x+2,y+2)]) # right side of the upside-down triangle

        if not isinstance(triangle['left'], str):
            texts = texts + text('$'+str(triangle['left']).replace('x','x_{').replace('b','b_{').replace('*','}')+'}$',(x-1,y+1), horizontal_alignment='left', rgbcolor=text_color) # left side of the upside-down triangle
        else:
            texts = texts + text(triangle['left'], (x-1,y+1), horizontal_alignment='left', rgbcolor=text_color) # left side of the upside-down triangle

        if not isinstance(triangle['top'], str):
            texts = texts + text('$'+str(triangle['top']).replace('x','x_{').replace('b','b_{').replace('*','}')+'}$',(x,y+2), vertical_alignment='bottom', rgbcolor=text_color) # top of the upside-down triangle
        else:
            texts = texts + text(triangle['top'], (x,y+2), vertical_alignment='bottom', rgbcolor=text_color) # top of the upside-down triangle

        if not isinstance(triangle['right'], str):
            texts = texts + text('$'+str(triangle['right']).replace('x','x_{').replace('b','b_{').replace('*','}')+'}$', (x+1,y+1), horizontal_alignment='left', rgbcolor=text_color) # right side of the upside-down triangle
        else:
            texts = texts + text(triangle['right'], (x+1,y+1), horizontal_alignment='left', rgbcolor=text_color) # right side of the upside-down triangle

        texts = texts + text ('$\\Delta_{'+str(k)+'}$', (x+0.2, y+0.7), rgbcolor=(0,1,0))

    return lines + texts
