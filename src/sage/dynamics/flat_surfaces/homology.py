r"""
Simplicial complex, homology of surfaces and translation surfaces.

EXAMPLES::

To initialize a ribbon graph with angles, you have to input the standard data to
initialize a ribbon graph plus a list of positive rational numbers which
corresponds to the angles between darts (more precisely, the number at position
i is the angle between i and v(i))::

    sage: e = '(0,1)(2,3)'
    sage: f = '(0,2,1,3)'
    sage: a = [1/2,1/2,1/2,1/2]
    sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
    sage: r.spin_parity()
    1
"""
from sage.all import FreeModule
from sage.matrix.constructor import matrix, identity_matrix

from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.misc.cachefunc import cached_method

################################
# permutation temporary stuff  #
# (to be moved somewhere else) #
################################
# we implement two kinds of permutations based on datatypes
#  list : dense permutation (p[i] is the image of i)
#  dictionnary : sparse permutation (p[i] is the image of i)

def perm_cycle_tuples(p):
    r"""
    Return the cycle decomposition of a permutation with partial domain
    """
    seen = [True] * len(p)
    cycles = []
    for i in xrange(len(p)):
        if p[i] is not None and seen[i]:
            cycle = []
            while seen[i]:
                seen[i] = False
                cycle.append(i)
                i = p[i]
            cycles.append(tuple(cycle))
    return cycles

def perm_invert(p):
    r"""
    TESTS::

        sage: from sage.dynamics.flat_surfaces.homology import perm_invert
        sage: perm_invert([])
        []
        sage: perm_invert([None, None])
        [None, None]
        sage: perm_invert([None, 3, None, 4, 1])
        [None, 4, None, 1, 3]
    """
    if isinstance(p,list):
        q = [None] * len(p)
        for i in xrange(len(p)):
            if p[i] is not None:
                q[p[i]] = i
        return q
    elif isinstance(p,dict):
        return dict((j,i) for i,j in q.iteritems())

    raise ValueError, "q should be either a list or a dict"

def perm_check(p):
    r"""
    TESTS::

        sage: from sage.dynamics.flat_surfaces.homology import perm_check
        sage: perm_check([])
        True
        sage: perm_check([None,None,3,2])
        True
        sage: perm_check([None,None])
        True
        sage: perm_check([None,0])
        False
        sage: perm_check([None,3])
        False
        sage: perm_check([1,1])
        False
    """
    if isinstance(p,list):
        N = len(p)
        seen = [False] * N
        for i in xrange(len(p)):
            if p[i] is not None:
                j = int(p[i])
                if j < 0 or j > N:
                    return False
                if seen[j] or p[j] is None:
                    return False
                seen[j] = True
        return True

    elif isinstance(p,dict):
        return set(p) == set(p.values())

    else:
        raise ValueError, "p should be either a list or a dict"

def perm_compose(p,q,sparse=None):
    r"""
    Return the product pq q of partial permutations ``p`` and ``q``

    the product is from left to right. In other words this function return the
    composition `q \circ p`.

    TESTS::

        sage: from sage.dynamics.flat_surfaces.homology import perm_compose
        sage: perm_compose([],[])
        []
        sage: perm_compose([None,2,3,1],[None,2,1,3])
        [None, 1, 3, 2]
    """
    if isinstance(p,list):
        if sparse is None or sparse is False:
            r = [None] * len(p)
        else:
                r = {}
        if isinstance(q,list):
            for i in xrange(len(p)):
                if p[i] is not None and p[i] < len(q):
                    r[i] = q[p[i]]
        elif isinstance(q,dict):
            for i in xrange(len(p)):
                if p[i] is not None and p[i] in q:
                    r[i] = q[p[i]]
        else:
            raise ValueError, "q should be either list or dict"
        return r

    elif isinstance(p,dict):
        if sparse is None or sparse is True:
            r = {}
        if isinstance(q,list):
            for i in p:
                if p[i] is not None and p[i] < len(q):
                    r[i] = q[p[i]]
        elif isinstance(q,dict):
            for i in p:
                if p[i] is not None and p[i] in q:
                    r[i] = q[p[i]]
        else:
            raise ValueError, "q should be either list or dict"
        return r

    else:
        raise ValueError, "p should be either list or dict"

def perm_domain(p):
    r"""
    Returns the domain of a partial permutation.

    TESTS::

        sage: from sage.dynamics.flat_surfaces.homology import perm_domain
        sage: perm_domain([])
        []
        sage: perm_domain([None,None])
        []
        sage: perm_domain([None,2,1,None])
        [1, 2]
    """
    if isinstance(p,list):
        return [i for i in xrange(len(p)) if p[i] is not None]
    elif isinstance(p,dict):
        return p.keys()
    else:
        raise ValueError, "p should be either list or dict"

def perm_range(p):
    r"""
    TESTS::

        sage: from sage.dynamics.flat_surfaces.homology import perm_range
        sage: perm_range([])
        []
        sage: perm_range([None,None])
        []
        sage: perm_range([None,2,1,None])
        [1, 2]
    """
    if isinstance(p,list):
        return sorted([p[i] for i in xrange(len(p)) if p[i] is not None])
    elif isinstance(p,dict):
        return sorted(p.values())
    else:
        raise ValueError, "p should be either list or dict"

def perm_orbit(p,i):
    if isinstance(p,list):
        if 0 <= i and i < len(p):
            j = p[i]
            if j is None:
                return None
            c = [i]
            while i != j:
                c.append(j)
                j = p[j]
    elif isinstance(p,dict):
        if i not in p:
            return None
        c = [i]
        j = p[i]
        while i != j:
            c.append(j)
            j = p[j]
    else:
        raise ValueError, "p should be either a list or a dict"

    return c

def from_cycles(c, sparse=False):
    r"""
    Build a partial permutation from cycles

    EXAMPLES::

        sage: from sage.dynamics.flat_surfaces.homology import from_cycles
        sage: from_cycles([])
        []
        sage: from_cycles([],sparse=True)
        {}
        sage: from_cycles([(1,6,2),(5,8)])
        [None, 6, 1, None, None, 8, 2, None, 5]
        sage: from_cycles([(1,6,2),(5,8)],sparse=True)
        {8: 5, 1: 6, 2: 1, 5: 8, 6: 2}
    """
    if sparse:
        p = {}
    else:
        N = 0
        for cc in c:
            N = max(N,max(cc)+1)
        p = [None] * (N)
    for cc in c:
        for i in xrange(len(cc)-1):
            p[cc[i]] = cc[i+1]
        p[cc[-1]] = cc[0]
    return p

def from_string(s, sparse=False):
    r"""
    Build a partial permutation from a string

    EXAMPLES::

        sage: from sage.dynamics.flat_surfaces.homology import from_string
        sage: from_string('()')
        []
        sage: from_string('()',sparse=True)
        {}
        sage: from_string('(1,6,2)(5,8)')
        [None, 6, 1, None, None, 8, 2, None, 5]
        sage: from_string('(1,6,2)(5,8)',sparse=True)
        {8: 5, 1: 6, 2: 1, 5: 8, 6: 2}
    """
    if s:
        if s[0] != '(' or s[-1] != ')':
            raise ValueError, "the syntax should be like (1,5)(3,6,2)*"
        if len(s) == 2:
            if sparse:
                return {}
            return []
        return from_cycles([map(int,c.split(',')) for c in s[1:-1].split(")(")],sparse=sparse)
    return []

def equalize_perms(p):
    r"""
    Put all permutations with the same length and set to identity where others
    are defined.
    """
    N = max(len(pp) for pp in p)
    for j in xrange(len(p)):
        p[j].extend([None]*(N-len(p[j])))
    for i in xrange(N):
        if any(p[j][i] is not None for j in xrange(len(p))):
            for j in xrange(len(p)):
                if p[j][i] is None:
                    p[j][i] = i

#############################
# Ribbon graph and homology #
#############################

def RibbonGraph(vertices=None,edges=None,faces=None,sparse=False,check=True):
    r"""
    Construct a ribbon graph from the given data.

    If only vertices or faces is provided then edges is initialized as the
    involution without fixed point which exchanges `2i` with `2i+1`.

    EXAMPLES::

        sage: RibbonGraph([],[],[])
        Ribbon graph with 1 vertex, 0 edge and 1 face
        sage: RibbonGraph('()','(0,1)','(0,1)')
        Ribbon graph with 2 vertices, 1 edge and 1 face

        sage: G = RibbonGraph('(0,3)(1,2)','(0,1)(2,3)','(0,2)(1,3)')
        sage: G
        Ribbon graph with 2 vertices, 2 edges and 2 faces
        sage: G.darts()
        [0, 1, 2, 3]
        sage: G.genus()
        0

        sage: G = RibbonGraph(edges='(1,3)(2,4)(5,7)(6,8)',faces='(1,2,3,4,5,6,7,8)')
        sage: G
        Ribbon graph with 1 vertex, 4 edges and 1 face
        sage: G.darts()
        [1, 2, 3, 4, 5, 6, 7, 8]
        sage: G.genus()
        2

        sage: G = RibbonGraph(vertices='(0,2,3,6)(1,4,5,7)')
        sage: G
        Ribbon graph with 2 vertices, 4 edges and 4 faces
        sage: G.edges()
        [(0, 1), (2, 3), (4, 5), (6, 7)]
        sage: G.genus()
        0
    """
    if sparse:
        raise NotImplementedError, "not yet implemented"


    if vertices is not None:
        if isinstance(vertices, str):
            vertices = from_string(vertices)
        elif isinstance(vertices, list):
            if len(vertices) != 0:
                if isinstance(vertices[0], tuple):
                    vertices = from_cycles(vertices)
                    if check: perm_check(vertices)

    if edges is not None:
        if isinstance(edges,str):
            edges = from_string(edges)
        elif isinstance(edges, list):
            if len(edges) != 0:
                if isinstance(edges[0],tuple):
                    edges = from_cycles(edges)
                    if check: perm_check(edges)

    if faces is not None:
        if isinstance(faces,str):
            faces = from_string(faces)
        elif isinstance(faces,list):
            if len(faces) != 0:
                if isinstance(faces[0],tuple):
                    faces = from_cycles(faces)
                    if check: perm_check(faces)
    if edges is None:
        if vertices is None and faces is None:
            raise ValueError, "at least vertices or faces should be not None"
        if not (vertices is None or faces is None):
            equalize_perms([vertices,faces])
            edges = perm_compose(perm_invert(vertices),perm_invert(faces))
        else:
            if vertices is None:
                n = len(faces)
            elif faces is None:
                n = len(vertices)
            if n%2:
                raise ValueError, "there should be an even number of darts"
            edges = []
            for i in xrange(n//2):
                edges.append(2*i+1)
                edges.append(2*i)
            if vertices is None:
                vertices = perm_compose(perm_invert(faces),perm_invert(edges))
            elif faces is None:
                faces = perm_compose(perm_invert(edges),perm_invert(vertices))
    elif vertices is None:
        if edges is None or faces is None:
            raise ValueError, "at least two of the entries should be not None"
        equalize_perms([edges,faces])
        vertices = perm_compose(perm_invert(faces),perm_invert(edges))
    elif faces is None:
        if vertices is None or edges is None:
            raise ValueError, "at least two of the entries should be not None"
        equalize_perms([vertices,edges])
        faces = perm_compose(perm_invert(edges),perm_invert(vertices))
    else:
        equalize_perms([vertices,edges,faces])

    return RibbonGraphDense(vertices,edges,faces,check=check)

class RibbonGraphGeneric(SageObject):
    r"""
    Generic class for Ribbon graph.

    A Ribbon graph is a graph embedded in a surface. This class implements many
    generic functionnalities as computation of homology. It should not be used
    directly as it does not manage any data.

    This class uses representation as a triple ``(v,e,f)`` of permutations such that
    `vef = 1` and the action of the group generated by  ̀`v,e,f` acts transitvely
    in the the domain. The cycles of ``v`` are considered as vertices, the ones
    of ``e`` are considered as edges and the ones of ``f`` as the faces. Each
    element of the domain is a half-edge which is called a *dart*. A dart is
    also associated to an oriented edge.

    Should be implemented in derived classes

    for XXX either ``vertex``, or ``edge`` or ``face``

     - XXXs(self) - cycles of XXX
     - dart_to_XXX(self,i)
     - num_XXX(self) - number of cycles of XXX
     - XXX_perm(self) - XXX as a permutation
     - XXX_orbit(self,i) - orbit of i for XXX

     also

     - darts - return the list of darts
     - num_darts - return the number of darts (cardinality of the domain)
    """
    def darts(self):
        raise NotImplementedError
    def num_darts(self):
        raise NotImplementedError
    def num_vertices(self):
        r"""
        Number of vertices
        """
        raise NotImplementedError
    def vertices(self):
        raise NotImplementedError
    def vertex_perm(self):
        raise NotImplementedError
    def vertex_orbit(self,i):
        raise NotImplementedError
    def num_edges(self):
        r"""
        Number of edges
        """
        raise NotImplementedError
    def edges(self):
        raise NotImplementedError
    def edge_perm(self):
        raise NotImplementedError
    def edge_orbit(self,i):
        raise NotImplementedError
    def num_faces(self):
        r"""
        Number of faces
        """
        raise NotImplementedError
    def faces(self):
        raise NotImplementedError
    def face_perm(self):
        raise NotImplementedError
    def face_orbit(self,i):
        raise NotImplementedError
    def dart_to_vertex(self,d):
        raise NotImplementedError

    def _repr_(self):
        r"""
        String representation.

        TESTS::

            sage: RibbonGraph(edges='(0,1)(2,3)',faces='(0,1,2,3)')._repr_()
            'Ribbon graph with 3 vertices, 2 edges and 1 face'
            sage: RibbonGraph(edges='(0,1)',faces='(0)(1)')._repr_()
            'Ribbon graph with 1 vertex, 1 edge and 2 faces'
            sage: RibbonGraph(edges='(0,1)(2,3)',faces='(0,2)(1,3)')._repr_()
            'Ribbon graph with 2 vertices, 2 edges and 2 faces'
        """
        n = self.num_vertices()
        if n <= 1:
            vert_str = "%d vertex" %n
        else:
            vert_str = "%d vertices" %n
        n = self.num_edges()
        if n <= 1:
            edge_str = "%d edge" %n
        else:
            edge_str = "%d edges" %n
        n = self.num_faces()
        if n <= 1:
            face_str = "%d face" %n
        else:
            face_str = "%d faces" %n
        return "Ribbon graph with %s, %s and %s" %(vert_str,edge_str,face_str)

    def dual(self):
        r"""
        Returns the dual Ribbon graph.

        The *dual* ribbon graph of `(v,e,f)` is `(f^{-1}, e, v^{-1})`.

        EXAMPLES::

            sage: r = RibbonGraph(edges='(0,1)',faces='(0)(1)'); r
            Ribbon graph with 1 vertex, 1 edge and 2 faces
            sage: r.dual()
            Ribbon graph with 2 vertices, 1 edge and 1 face
        """
        return RibbonGraph(
            vertices=perm_invert(self._faces),
            edges=self._edges,
            faces=perm_invert(self._vertices))

    #
    # euler characteristic
    #

    def euler_characteristic(self):
        r"""
        Returns the Euler characteristic of the embedded surface.

        The *Euler characteristic* of a surface complex is `V - E + F`, where
        `V` is the number of vertices, `E` the number of edges and `F` the
        number of faces.

        EXAMPLES::

            sage: r = RibbonGraph(edges='(0,1)(2,3)(4,5)',faces='(0,2,4)(3,5)')
            sage: r.euler_characteristic()
            2

            sage: r = RibbonGraph(edges='(0,1)(2,3)',faces='(0,2,1,3)')
            sage: r.euler_characteristic()
            0
        """
        return self.num_vertices() - self.num_edges() + self.num_faces()

    def is_plane(self):
        r"""
        Returns true if and only if the ribbon graph belongs in a sphere. In
        other words if it has genus 0.

        EXAMPLES::

            sage: r = RibbonGraph(vertices='(0)(1)',edges='(0,1)')
            sage: r.is_plane()
            True

            sage: r = RibbonGraph(vertices='(0,1)',edges='(0,1)')
            sage: r.is_plane()
            True

            sage: r = RibbonGraph(edges='(0,1)(2,3)',faces='(0,2)(1,3)')
            sage: r.is_plane()
            True

            sage: r = RibbonGraph(edges='(0,1)(2,3)',faces='(0,2,1,3)')
            sage: r.is_plane()
            False
        """
        return self.euler_characteristic() == 2

    def is_plane_tree(self):
        r"""
        Returns True if and only if the ribbon graph is a planar tree. In other
        words, it has genus 0 and only one face.

        EXAMPLES::

            sage: r = RibbonGraph(vertices='(0)(1)',edges='(0,1)')
            sage: r.is_plane_tree()
            True

            sage: r = RibbonGraph(vertices='(0)(1,2,4)(3)(5)',edges='(0,1)(2,3)(4,5)')
            sage: r.is_plane_tree()
            True

            sage: r = RibbonGraph(vertices='(0,1)',edges='(0,1)')
            sage: r.is_plane_tree()
            False
            sage: r.is_plane()
            True
        """
        return (self.num_faces() == 1 and self.genus() == 0)

    def is_triangulated(self):
        r"""
        Returns True if the surface is triangulated. In other words, faces
        consist only of the product of 3-cycles.

        EXAMPLES::

            sage: r = RibbonGraph(edges='(0,1)(2,3)(4,5)',faces='(0,2,4)(1,5,3)')
            sage: r.is_triangulated()
            True

            sage: r = RibbonGraph(edges='(0,1)(2,3)',faces='(0,2,1,3)')
            sage: r.is_triangulated()
            False
        """
        return all(len(c) == 3 for c in self.faces())

    def genus(self):
        r"""
        Return the genus of the surface associated to this Ribbon graph.

        EXAMPLES::

            sage: R = RibbonGraph(vertices='(1)(2)',edges='(1,2)')
            sage: R.genus()
            0

            sage: e='(1,3)(2,4)'
            sage: f='(1,2,3,4)'
            sage: RibbonGraph(edges=e,faces=f).genus()
            1

            sage: e='(1,3)(2,4)(5,7)(6,8)'
            sage: f='(1,2,3,4,5,6,7,8)'
            sage: RibbonGraph(edges=e,faces=f).genus()
            2

            sage: e='(1,3)(2,4)(5,7)(6,8)(9,11)(10,12)'
            sage: f='(1,2,3,4,5,6,7,8,9,10,11,12)'
            sage: RibbonGraph(edges=e,faces=f).genus()
            3
        """
        return 1 - self.euler_characteristic()//2

    #
    # cycles and fundamental group
    #

    def spanning_tree(self):
        r"""
        Return a spanning tree

        OUTPUT:

        - spanning tree as a DiGraph

        - remaining edges as 2-tuples ``(i,e[i])``

        EXAMPLES::

            sage: R = RibbonGraph('(1,2,3)','(1,2)(3,4)')
            sage: R
            Ribbon graph with 2 vertices, 2 edges and 2 faces
            sage: T,o = R.spanning_tree()
            sage: T
            Digraph on 2 vertices
            sage: T.edges()
            [(0, 1, (3, 4))]
            sage: o
            [(1, 2)]

            sage: R = RibbonGraph('(1,2,3)(4,5,6)','(1,2)(3,4)(5,6)')
            sage: R
            Ribbon graph with 2 vertices, 3 edges and 3 faces
            sage: T,o = R.spanning_tree()
            sage: T
            Digraph on 2 vertices
            sage: T.edges()
            [(0, 1, (3, 4))]
            sage: o
            [(1, 2), (5, 6)]

            sage: e = '(1,3)(5,7)(2,4)(6,8)'
            sage: f = '(1,2,3,4,5,6,7,8)'
            sage: R = RibbonGraph(edges=e, faces=f)
            sage: T,o = R.spanning_tree()
            sage: T
            Digraph on 1 vertex
            sage: o
            [(1, 3), (2, 4), (5, 7), (6, 8)]
        """
        from sage.graphs.digraph import DiGraph

        d = self.darts()
        v = self.vertices()
        e = self.edge_perm()

        if self.num_darts() == 0:
            return DiGraph(),[]
        if self.num_vertices() == 1:
            return DiGraph({0:[]}),self.edges()

        T = DiGraph()

        v0 = 0
        for root in v[0]:
            v1 = self.dart_to_vertex(e[root])
            if v1 != 0:
                break

        T.add_edge(v0,v1,(root,e[root]))
        o = []
        seen = set([v0,v1])       # seen vertices
        waiting = [e[root],root]  # waiting darts
        cc = []

        while waiting:
            ii = waiting.pop() # this is a dart
            v0 = self.dart_to_vertex(ii)
            seen.add(v0)
            for j in self.vertex_orbit(ii)[1:]:
                v1 = self.dart_to_vertex(e[j])
                if v1 in seen:
                    if j < e[j]:
                        o.append((j,e[j]))
                else:
                    T.add_edge(v0,v1,(j,e[j]))
                    waiting.append(e[j])
                    seen.add(v1)

        return T, sorted(o)

    def collapse(self,spanning_tree=None,sparse=False):
        r"""
        Return a ribbon graph callapsed along a spanning tree.

        The resulting graph is on the same surface as the preceding but has only
        one vertex. It could be used twice to provide a polygonal representation
        with one vertex and one face.

        EXAMPLES::

            sage: R = RibbonGraph(vertices='(0,1,2,5)(3,7)(4,10,9)(6,11,12)(8,13)')
            sage: R.genus()
            1
            sage: R.num_vertices()
            5
            sage: R.num_edges()
            7
            sage: R.num_faces()
            2
            sage: R2 = R.collapse()
            sage: R2
            Ribbon graph with 1 vertex, 3 edges and 2 faces
            sage: R
            Ribbon graph with 5 vertices, 7 edges and 2 faces
            sage: R3 = R2.dual().collapse().dual()
            sage: R3
            Ribbon graph with 1 vertex, 2 edges and 1 face

        """
        from copy import deepcopy

        if spanning_tree is None:
            spanning_tree,_ = self.spanning_tree()

        darts_to_kill = set([])
        for v0,v1,e in spanning_tree.edges():
            darts_to_kill.add(e[0])
            darts_to_kill.add(e[1])

        new_edges = []
        for e in self.edges():
            if e[0] not in darts_to_kill:
                new_edges.append(e)

        new_faces = []
        for f in self.faces():
            ff = tuple(i for i in f if i not in darts_to_kill)
            if ff:
                new_faces.append(ff)

        return RibbonGraph(edges=new_edges,faces=new_faces,sparse=sparse)

    def boundaries(self):
        r"""
        Return the list of cycles which are boundaries.

        A cycle is a *boundary* if it bounds a face.

        EXAMPLES::

            sage: r = RibbonGraph('(1,2,3)(4,5,6)','(1,2)(3,4)(5,6)')
            sage: r.boundaries()
            [[(1, 2)],  [(2, 1), (3, 4), (6, 5), (4, 3)], [(5, 6)]]

            sage: r = RibbonGraph('(1,2,3)(4,5)(6,7,8)',edges='(1,2)(3,4)(5,6)(7,8)')
            sage: r.boundaries()
            [[(1, 2)],  [(2, 1), (3, 4), (5, 6), (8, 7), (6, 5), (4, 3)], [(7, 8)]]
        """
        e = self.edge_perm()
        return sorted([[(i,e[i]) for i in f] for f in self.faces()])

    def cycle_basis(self, intersection=False, verbose=False):
        r"""
        Returns a base of oriented cycles of the Ribbon graph modulo boundaries.

        If ``intersection`` is set to True then the method also returns the
        intersection matrix of the cycles.

        EXAMPLES::

            sage: r = RibbonGraph('(1,2,3)(4,5,6)','(1,2)(3,4)(5,6)')
            sage: r.cycle_basis()
            []

            sage: r = RibbonGraph('(1,2,3)(4,5)(6,7,8)',edges='(1,2)(3,4)(5,6)(7,8)')
            sage: r.cycle_basis()
            []

            sage: r = RibbonGraph('(1,4,5)(2,3)(6,7,8)',edges='(1,2)(3,4)(5,6)(7,8)')
            sage: r.cycle_basis()
            []

            sage: e = '(1,3)(2,4)(5,7)(6,8)'
            sage: f = '(1,2,3,4,5,6,7,8)'
            sage: r = RibbonGraph(edges=e,faces=f)
            sage: r.cycle_basis()
            [[(1, 3)], [(2, 4)], [(5, 7)], [(6, 8)]]

            sage: f = '(0,10,13)(6,17,11)(2,14,7)(15,12,3)(16,20,19)(18,1,9)(4,22,21)(23,8,5)'
            sage: e = [(i,i+1) for i in xrange(0,24,2)]
            sage: r = RibbonGraph(edges=e,faces=f); r
            Ribbon graph with 2 vertices, 12 edges and 8 faces
            sage: c,m = r.cycle_basis(intersection=True)
            sage: c
            [[(0, 1), (4, 5)], [(8, 9)], [(12, 13)], [(14, 15), (1, 0)]]
            sage: m
            [ 0  1  0  0]
            [-1  0  0  0]
            [ 0  0  0  1]
            [ 0  0 -1  0]
        """
        T,o = self.spanning_tree()

        # build a Ribbon graph with one vertex and one face
        r = self.collapse(T).dual().collapse().dual()
        if T is None:
            return r.edges()

        if intersection:
            c = r.vertices()[0]
            M = len(c)
            I = []

        cycles = [] # the cycles
        for e in r.edges():
            if verbose:
                print "build cycle from edge %s between vertex v0=%d and v1=%d" %(str(e),self.dart_to_vertex(e[0]),self.dart_to_vertex(e[1]))

            # build the branch to the root from v0
            v0 = self.dart_to_vertex(e[0])
            if verbose:
                print " build branch from v0=%d" %v0
            p0 = []
            while v0 != 0:
                v0,_,e0 = T.incoming_edges(v0)[0] # (v_in,v_out,label)
                p0.append(e0)
                if verbose:
                    print " add %d" %v0
            if verbose:
                print " branch is %s" %str(p0)
            # build the branch to the root from v1

            v1 = self.dart_to_vertex(e[1])
            if verbose:
                print " build branch from v1=%d" %v1
            p1 = []
            while v1 != 0:
                v1,_,e1 = T.incoming_edges(v1)[0]
                p1.append(e1)
                if verbose:
                    print " add %d" %v1
            if verbose:
                print " branch is %s" %str(p1)
            # clean the branches by removing common part
            while p0 and p1 and p0[-1] == p1[-1]:
                if verbose:
                    print "find common element",p0[-1]
                p0.pop(-1)
                p1.pop(-1)

            # add the cycle to the list
            cycles.append((p0,e,p1))

            # compute algebraic intersection with preceding cycles
            if intersection:
                i = []
                for _,ee,_ in cycles:
                    if verbose:
                        print "compute intersection"
                    p_in = c.index(e[1])
                    p_out = (c.index(e[0]) - p_in) % M
                    q_in  = (c.index(ee[1]) - p_in) % M
                    q_out = (c.index(ee[0]) - p_in) % M
                    if verbose:
                        print "  after reduction: p_out = %d, q_in = %d, q_out = %d" %(p_out,q_in,q_out)

                    # compute intersection
                    # p_in = 0 and the others 3 are positive
                    if q_in < p_out and p_out < q_out:
                        i.append(1)
                    elif q_out < p_out and p_out < q_in:
                        i.append(-1)
                    else:
                        i.append(0)

                I.append(i)

        # make cycle as list
        cycles = [p0[::-1]+[e]+[c[::-1] for c in p1] for p0,e,p1 in cycles]

        if intersection:
            m = matrix(len(cycles))
            for j in xrange(len(I)):
                for jj in xrange(len(I[j])):
                    m[j,jj] = I[j][jj]
                    m[jj,j] = -I[j][jj]

            return cycles, m
        return cycles

    def is_cycle(self,c):
        r"""
        Test whether ``c`` is a cycle.

        A *path* is a sequence of oriented edges such that each edge starts
        where the preceding one ends. A *cycle* is a path which starts where it
        ends.
        """
        for i in xrange(len(c)-1):
            if self.dart_to_vertex(c[i][1]) != self.dart_to_vertex(c[i+1][0]):
                return False
        if self.dart_to_vertex(c[-1][1]) != self.dart_to_vertex(c[0][0]):
            return False
        return True

    # chain complex of homology

    @cached_method
    def chain_complex(self,ring=None):
        r"""
        Return the chain complex associated to self
        """
        return SurfaceChainComplex(self,ring)

class RibbonGraphDense(RibbonGraphGeneric):
    r"""
    Dense ribbon graph.

    The set of active darts is a subset of range(0,N)

    A dense ribbon graph has the following attributes

      - total_darts - non negative integer - the total number darts
      - num_darts - non negative integer - the number of active darts
      - active_darts - bitset - list of lengths _total_darts with True or
        False. The position i is True if i is an active dart.

      - vertices, vertices_inv - list - partial permutations of [0,N] which are
        inverse of each other
      - vertex_cycles - the cycles of the partial permutation vertices
      - dart_to_vertex_index

      - edges, edges_inv - list - partial permutations of [0,N] which are
        inverse of each other
      - edge_cycles - the cycles of the partial permutation edge
      - dart_to_edge_index

      - faces, faces_inv - list - partial permutations of [0,N] which are
        inverse of each other
      - face_cycles - the cycles of the partial permutation faces
      - dart_to_face_index

    TODO:

    - orientation of the edges
    """
    def __init__(self,vertices,edges,faces,check=True):
        self._total_darts = total_darts = len(vertices)
        self._active_darts = [False] * total_darts
        num_darts = 0
        assert(len(edges) == total_darts and len(faces) == total_darts)
        for i in xrange(total_darts):
            if vertices[i] is None:
                assert(edges[i] is None and faces[i] is None)
            else:
                assert(edges[i] is not None and faces[i] is not None)
                num_darts += 1
                self._active_darts[i] = True
        self._num_darts = num_darts

        self._vertices = vertices
        self._vertex_cycles = perm_cycle_tuples(vertices)
        self._dart_to_vertex_index = [None] * total_darts
        for i,c in enumerate(self._vertex_cycles):
            for j in c:
                self._dart_to_vertex_index[j] = i

        self._edges = edges
        self._edge_cycles = perm_cycle_tuples(edges)
        self._dart_to_edge_index = [None] * total_darts
        for i,e in enumerate(self._edge_cycles):
            self._dart_to_edge_index[e[0]] = i
            self._dart_to_edge_index[e[1]] = i

        self._faces = faces
        self._face_cycles = perm_cycle_tuples(faces)
        self._dart_to_face_index = [None] * total_darts
        for i,c in enumerate(self._face_cycles):
            for j in c:
                self._dart_to_face_index[j] = i

        self._vertices_inv = [None] * total_darts
        self._edges_inv = [None] * total_darts
        self._faces_inv = [None] * total_darts
        for i in xrange(total_darts):
            if self._vertices[i] is not None:
                self._vertices_inv[self._vertices[i]] = i
                self._edges_inv[self._edges[i]] = i
                self._faces_inv[self._faces[i]] = i

        if check:
            self._check()

    def is_sparse(self):
        r"""
        Returns False.
        """
        return False

    def add_extra_darts(self,n):
        r"""
        Add extra darts to the current vertex in order to support a total of
        ``n`` darts.
        """
        m = self._total_darts
        if n > m:
            self._total_darts = int(n)
            for p in [self._vertices, self._vertices_inv,
                      self._edges,self._edges_inv,
                      self._faces, self._faces_inv]:
                p.extend([None] * (n - m))
            self._active_darts.extend([False] * (n-m))
            self._check()

    def _check(self):
        r"""
        Check that the data of the Ribbon graph is coherent
        """
        from sage.graphs.graph import Graph
        G = Graph()

        if len(self._active_darts) != self._total_darts:
            raise ValueError, "the length of active darts is not total_darts"
        if self._active_darts.count(True) != self._num_darts:
            raise ValueError, "the number of darts do not coincide with active darts"

        for i in xrange(self._total_darts):
            if self._active_darts[i]:
                G.add_edge(i,self._vertices[i])
                G.add_edge(i,self._edges[i])
                G.add_edge(i,self._faces[i])
                if self._vertices[i] is None or self._vertices_inv[i] is None:
                    raise ValueError, "dart %d is active but has no vertex" %i
                if self._edges[i] is None or self._edges_inv[i] is None:
                    raise ValueError, "dart %d is active but has no edge" %i
                if self._faces[i] is None or self._faces_inv[i] is None:
                    raise ValueError, "dart %d is active but has no face" %i

                if self._vertices[self._vertices_inv[i]] != i:
                    raise ValueError, "vertices is not the inverse of vertices_inv"
                if self._edges[self._edges_inv[i]] != i:
                    raise ValueError, "edges is not the inverse of edges_inv"
                if self._faces[self._faces_inv[i]] != i:
                    raise ValueError, "faces is not the inverse of faces_inv"

                if self._faces[self._edges[self._vertices[i]]] != i:
                    raise ValueError, "the Ribbon graph condition vef=() is not satisfied for %d" %i
                if self._edges[i] == i or self._edges[self._edges[i]] != i:
                    raise ValueError, "edges is not an involution without fixed point for %d" %i

            else:
                if self._vertices[i] is not None or self._vertices_inv[i] is not None:
                    raise ValueError, "dart %d is not active but has a vertex" %i
                if self._edges[i] is not None or self._edges_inv[i] is not None:
                    raise ValueError, "dart %d is not active but has an edge" %i
                if self._faces[i] is not None or self._faces_inv[i] is not None:
                    raise ValueError, "dart %d is not active but has a face" %i

        if not G.is_connected():
            raise ValueError, "the graph is not connected"

    def relabel(self, perm=None):
        r"""
        perm is a of range(0,N)

        If ``perm`` is None, relabel the darts on 0,2M keeping the relative
        order of the darts.
        """
        if perm is None:
            perm=[None]*self.num_darts()
            k = 0
            for i in xrange(self.num_darts()):
                if self._active_darts[i]:
                    perm[i] = k
                    k += 1

        vertices = [None] * self.num_darts()
        edges = [None] * self.num_darts()
        faces = [None] * self.num_darts()
        for i in xrange(self.num_darts()):
            if self._active_darts[i]:
                vertices[perm[i]] = perm[self._vertices[i]]
                edges[perm[i]] = perm[self._edges[i]]
                faces[perm[i]] = perm[self._faces[i]]

        return RibbonGraphDense(vertices,edges,faces)

    #
    # Darts
    #

    def num_darts(self):
        r"""
        Returns the number of darts.
        """
        return self._num_darts

    def darts(self):
        r"""
        Return the list of darts
        """
        return [i for i in xrange(self._total_darts) if self._active_darts[i]]

    def num_vertices(self):
        r"""
        Returns the number of vertices.
        """
        return max(1,len(self._vertex_cycles))

    def vertex_perm(self):
        r"""
        Returns the permutation that define the vertices.
        """
        return self._vertices

    def vertex_orbit(self, i):
        r"""
        Return the orbit of ``i`` under the permutation that define the
        vertices.
        """
        if self._active_darts[i]:
            return perm_orbit(self._vertices,i)
        return None

    def vertices(self):
        r"""
        Return the list of vertices as cycles decomposition of the vertex
        permutation.
        """
        return self._vertex_cycles

    def dart_to_vertex(self,i):
        r"""
        Return the vertex on which the dart ``i`` is attached.
        """
        if self._active_darts[i]:
            return self._dart_to_vertex_index[i]
        raise ValueError, "dart %d is not active" %i

    #
    # Edges
    #

    def num_edges(self):
        r"""
        Returns the number of edges.
        """
        return len(self._edge_cycles)

    def edge_perm(self):
        r"""
        Return the permutation that define the edges.
        """
        return self._edges

    def edge_orbit(self, i):
        r"""
        Return the orbit of the dart ``i`` under the permutation that defines
        the edges.
        """
        if self._active_darts[i]:
            return perm_orbit(self._edges,i)
        return None

    def edges(self):
        r"""
        Return the set of edges.
        """
        return self._edge_cycles

    def dart_to_edge(self, i, orientation=False):
        r"""
        Returns the edge the darts ``i`` belongs to.

        If orientation is set to ``True`` then the output is a `2`-tuple
        ``(e,o)`` where ``e`` is the index of the edge and ``o`` is its
        orientation as ``+1`` or ``-1``.
        """
        if self._active_darts[i]:
            if not orientation:
                return self._dart_to_edge_index[i]
            j = self._dart_to_edge_index[i]
            if i == self._edge_cycles[j][0]:
                return (j,1)
            elif i == self._edge_cycles[j][1]:
                return (j,-1)
            else:
                raise ValueError, "this should not happen!"
        raise ValueError, "dart %d is not active" %i

    #
    # Faces
    #

    def num_faces(self):
        r"""
        Return the number of faces.
        """
        return max(1,len(self._face_cycles))

    def face_perm(self):
        r"""
        Return the permutation that defines the face.
        """
        return self._faces

    def face_orbit(self, i):
        r"""
        Return the orbit of ``i`` under the permutation associated to faces.
        """
        if self._active_darts[i]:
            return perm_orbit(self._faces,i)
        return None

    def faces(self):
        r"""
        Return the list of faces.
        """
        return self._face_cycles

    def dart_to_face(self, i):
        if self._active_darts[i]:
            return self._dart_to_face_index[i]
        raise ValueError, "dart %d is not active" %i

#TODO: doc + test + _check()
class RibbonGraphSparse(RibbonGraphGeneric):
    r"""
    Implementation using dictionnaries
    """
    def __init__(self,vertices,edges,faces,check=True):
        total_darts = len(vertices)
        assert(len(edges) == total_darts and len(faces) == total_darts)

        self._vertices = vertices
        self._vertex_cycles = perm_cycle_tuples(vertices)
        self._dart_to_vertex_index = {}
        for i,c in enumerate(self._vertex_cycles):
            for j in c:
                self._dart_to_vertex_index[j] = i

        self._edges = edges
        self._edge_cycles = perm_cycle_tuples(edges)
        self._dart_to_edge_index = {}
        for i,c in enumerate(self._edge_cycles):
            for j in c:
                self._dart_to_edge_index[j] = i

        self._faces = faces
        self._face_cycles = perm_cycle_tuples(faces)
        self._dart_to_face_index = {}
        for i,c in enumerate(self._face_cycles):
            for j in c:
                self._dart_to_edge_index[j] = i

        self._vertices_inv = dict((j,i) for i,j in self._vertices.iteritems())
        self._edges_inv = dict((j,i) for i,j in self._edges.iteritems())
        self._faces_inv = dict((j,i) for i,j in self._faces.iteritems())

        self._check()

    #TODO
    def _check():
        pass

    def is_sparse(self):
        return True

    def darts(self):
        return sorted(self._vertices.keys())

    def num_darts(self):
        return len(self._vertices)

    def dart_to_vertex(self,i):
        return self._dart_to_vertex[i]

    def num_vertices(self):
        return len(self._vertex_cycles)

    def vertices(self):
        return self._vertex_cycles

    def vertex_perm(self):
        return self._vertices

    def vertex_orbit(self,i):
        return perm_orbit(self._vertices,i)

    def dart_to_edge(self,i):
        return self._dart_to_edge(i)

    def num_edges(self):
        return len(self._edge_cycles)

    def edges(self):
        return self._edge_cycles

    def edge_perm(self):
        return self._edges

    def edge_orbit(self,i):
        return perm_orbit(self._edges,i)

    def dart_to_face(self,i):
        return self._dart_to_face(i)

    def num_faces(self):
        return len(self._face_cycles)

    def faces(self):
        return self._face_cycles

    def face_perm(self):
        return self._faces

    def face_orbit(self,i):
        return perm_orbit(self._faces,i)

class RibbonGraphWithAngles(RibbonGraphDense):
    r"""
    A Ribbon graph with angles between edges

    Currently angles can only be *rational* multiples of pi.

    TODO:

    - allows any kind of angles by providing a sum for the total and considering
      each angle as a (projective) portion of the total angle.
    """
    def __init__(self, vertices=None, edges=None, faces=None, angles=None):
        r = RibbonGraph(vertices,edges,faces)
        RibbonGraphDense.__init__(self,r.vertex_perm(),r.edge_perm(),r.face_perm())

        if len(angles) != self.num_darts():
            raise ValueError, "there are %d angles and %d darts" %(len(angles),self.num_darts())
        self._angles = map(QQ,angles)
          # angle between a dart and its vertex-neighbour
          # (rational number as multiple of pi)

        self._total_angle = []
          # total angle around vertices
          # (integer which corresponds to a multiple of pi)

        for v in self.vertices():
            self._total_angle.append(sum(angles[i] for i in v))

        for f in self.faces():
            a = sum(self._angles[i] for i in f)
            if a != len(f)-2:
                raise ValueError, "the angle of a face should be (nb_edges - 2) x pi"

    def angle_between_darts(self, d1, d2):
        r"""
        Return the angle between the darts ``d1`` and ``d2``
        """
        v = self.vertex_orbit(d1)
        if d2 not in v:
            raise ValueError, "d1=%s and d2=%s are not at the same vertex" %(str(d1),str(d2))

        a = 0
        i = 0
        while v[i] != d2:
            a += self._angles[v[i]]
            i += 1
        return a

    def angle_at_vertex(self,v):
        r"""
        Angle at a vertex (coefficient of pi)
        """
        return self._total_angle[v]

    def angle_at_vertices(self):
        r"""
        Return the list of angles at a vertex.
        """
        return self._total_angle

    def winding(self, c):
        r"""
        Return winding number along the cycle ``c``.

        This is NOT well defined because it depends on the way we choose to pass
        on the left or on the right at singularity.
        """
        a = 0
        for i in xrange(len(c)-1):
            d1 = c[i][1]
            d2 = c[i+1][0]
            if self.dart_to_vertex(d1) != self.dart_to_vertex(d2):
                raise ValueError, "c is not a cycle"
            a += self.angle_between_darts(d1,d2)-1
        d1 = c[-1][1]
        d2 = c[0][0]
        if self.dart_to_vertex(d1) != self.dart_to_vertex(d2):
            raise ValueError, "c is not a cycle"
        a += self.angle_between_darts(d1,d2)-1

        return a

    def holonomy_representation(self):
        r"""
        Return the holonomy representation in `SO(2)` as two lists.

        The first list correspond to cycles around vertices, while the second
        correspond to a cycle basis that generate homology.

        EXAMPLES::

            sage: e = '(0,1)(2,3)'
            sage: f = '(0,2,1,3)'
            sage: a = [1/2,1/2,1/2,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.holonomy_representation()
            ([0], [0, 0])

        The standard cube::

            sage: e = [(i,i+1) for i in xrange(0,24,2)]
            sage: f = '(0,20,7,10)(16,22,19,21)(2,9,5,23)(14,3,17,1)(12,8,15,11)(18,4,13,6)'
            sage: a = [1/2]*24
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.holonomy_representation()
            ([3/2, 3/2, 3/2, 3/2, 3/2, 3/2, 3/2, 3/2], [])

        Two copies of a triangle::

            sage: e = '(0,1)(2,3)(4,5)'
            sage: f = '(0,2,4)(1,5,3)'
            sage: a = [1/2,1/6,1/3,1/3,1/6,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.holonomy_representation()
            ([1, 1/2, 1/2], [])

            sage: a = [1/3,7/15,1/5,1/5,7/15,1/3]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.holonomy_representation()
            ([2/3, 2/3, 2/3], [])
        """
        from sage.functions.other import floor

        l1 = []
        for c in xrange(self.num_vertices()):
            w = self.angle_at_vertex(c)
            l1.append(w - 2*floor(w/2))

        l2 = []
        for c in self.cycle_basis():
            w = self.winding(c)
            l2.append(w - 2*floor(w/2))

        return l1,l2

    def has_trivial_holonomy(self):
        r"""
        Test whether self has trivial holonomy representation

        EXAMPLES::

            sage: e = '(0,1)(2,3)'
            sage: f = '(0,2,1,3)'
            sage: a = [1/2,1/2,1/2,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.has_trivial_holonomy()
            True

            sage: e = '(0,1)(2,3)(4,5)'
            sage: f = '(0,2,4)(1,5,3)'
            sage: a = [1/3,7/15,1/5,1/5,7/15,1/3]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.has_trivial_holonomy()
            False
        """
        l1,l2 = self.holonomy_representation()
        return all(i==0 for i in l1) and all(i==0 for i in l2)

    def spin_parity(self,check=True,verbose=False):
        r"""
        Return the spin parity of the Ribbon graph with angles.

        The surface should be holonomy free and with odd multiple of 2 pi
        angles.

        EXAMPLES:

        We first consider the case of the torus::

            sage: e = '(0,1)(2,3)'
            sage: f = '(0,2,1,3)'
            sage: a = [1/2,1/2,1/2,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.spin_parity()
            1

        Then the case of genus 2 surface (with an angle of 6pi)::

            sage: e = '(0,1)(2,3)(4,5)(6,7)'
            sage: f = '(0,2,4,3,6,1,7,5)'
            sage: a = [1/2,1/2,1,1/2,1/2,1,3/2,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.spin_parity()
            1

            sage: e = '(0,1)(2,3)(4,5)(6,7)'
            sage: f = '(0,2,4,6,1,3,5,7)'
            sage: a = [1/2,1/2,1,1,1,1,1/2,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.spin_parity()
            1

            sage: e = '(0,1)(2,3)(4,5)(6,7)'
            sage: f = '(0,2,4,6,1,3,5,7)'
            sage: a = [3/4]*8
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.spin_parity()
            1

        In genus 3 two spin parities occur for one conical angle 10pi::

            sage: e = '(0,1)(2,3)(4,5)(6,7)(8,9)(10,11)'
            sage: f1 = '(0,4,6,8,10,2,1,9,11,5,7,3)'
            sage: f2 = '(0,4,6,8,10,2,1,5,7,9,11,3)'
            sage: a = [1/2,1/2,1/2,1/2] + [1]*8
            sage: r1 = RibbonGraphWithAngles(edges=e,faces=f1,angles=a)
            sage: r1.spin_parity()
            1
            sage: r2 = RibbonGraphWithAngles(edges=e,faces=f2,angles=a)
            sage: r2.spin_parity()
            0
        """
        from sage.rings.finite_rings.constructor import GF
        # mod F2 we have: q(x+y) = B(x,y) + q(x) + q(y)

        if not self.has_trivial_holonomy():
            raise ValueError, "the surface does not have trivial holonomy"
        if any((i+2)%4 for i in self.angle_at_vertices()):
            raise ValueError, "each angle should be odd multiple of 2pi"

        GF2 = GF(2)

        c,M = self.cycle_basis(intersection=True)

        winding = []
        for cc in c:
            w = self.winding(cc)
            if w % 2 != 0:
                raise ValueError, "fatal error ! each winding should be a multiple of 2"
            winding.append(GF2(w//2))

        if verbose:
            print "cycles with winding"
            for i in xrange(len(c)):
                print c[i], winding[i]
            print "intersection matrix on Z"
            print M

        # compute a base change to get a symplectic basis
        _,P = M.symplectic_form()
        M = M.change_ring(GF2)
        P = P.change_ring(GF2)
        if verbose:
            print "base change for symplectic basis on GF(2)"
            print P

        g = self.genus()

        s = GF2(0)
        for i in xrange(g):
            # 1. computation of q(P.row(i))
            a = P.row(i)
            a_indices = [j for j in xrange(2*g) if a[j] != 0]
            ## winding + nb_components
            t_a = sum(winding[i]+1 for i in a_indices)
            ## self intersection
            for j1 in xrange(len(a_indices)):
                for j2 in xrange(j1+1,len(a_indices)):
                    t_a += M[a_indices[j1],a_indices[j2]]

            # 2. computation of q(P.row(g+i))
            b = P.row(g+i)
            b_indices = [j for j in xrange(2*g) if b[j] != 0]
            ## winding + nb_components
            t_b = sum(winding[i]+1 for i in b_indices)
            ## self intersection
            for j1 in xrange(len(b_indices)):
                for j2 in xrange(j1+1,len(b_indices)):
                    t_b += M[b_indices[j1],b_indices[j2]]

            # 3. add to s the contribution of the couple
            if verbose:
                print "contribution from %d is %d * %d = %d"%(i,t_a,t_b,t_a*t_b)
            s += t_a*t_b

        return s

##############################
# Chain complex and homology #
##############################

class SurfaceChainComplex(SageObject):
    r"""
    Surface chain complex associated to a Ribbon graph

    -> intersection of cycles

    EXAMPLES::

        sage: g = RibbonGraph(edges='(0,1)(2,3)(4,5)(6,7)',faces='(0,2,4,6,1,3,5,7)')
        sage: C = g.chain_complex()
        sage: C.chain_space(0).rank() == g.num_vertices()
        True
        sage: C.chain_space(1).rank() == g.num_edges()
        True
        sage: C.chain_space(2).rank() == g.num_faces()
        True

        sage: C.cycle_space(0).rank() - C.boundary_space(0).rank() == 0
        True
        sage: C.cycle_space(1).rank() - C.boundary_space(1).rank() == 2 * g.genus()
        True
        sage: C.cycle_space(2).rank() - C.boundary_space(2).rank() == 1
        True
    """
    def __init__(self, g, ring=None):
        self._ribbon_graph = g
        if ring is None:
            ring = ZZ
        self._base_ring = ring
        self._chain_spaces = [
            FreeModule(ring,g.num_vertices()),
            FreeModule(ring,g.num_edges()),
            FreeModule(ring,g.num_faces())]

        self._derivatives = []
        m = matrix(ring, [1]*g.num_vertices())
        self._derivatives.append(m)

        m = matrix(ring, g.num_vertices(), g.num_edges())
        for i,e in enumerate(g.edges()):
            m[g.dart_to_vertex(e[0]),i] -= 1
            m[g.dart_to_vertex(e[1]),i] += 1
        self._derivatives.append(m)

        m = matrix(ring, g.num_edges(), g.num_faces())
        for i,f in enumerate(g.faces()):
            for e in f:
                j,o = g.dart_to_edge(e,orientation=True)
                m[j,i] += o
        self._derivatives.append(m)

    def parent(self):
        r"""
        Returns the ribbon graph from which this chain complex is defined.
        """
        return self._ribbon_graph

    def _repr_(self):
        return "Chain complex over %s" %str(self._base_ring)

    def differential(self, degree=None):
        r"""
        Returns a differential of given degree.
        """
        if degree is None:
            return self._derivatives
        return self._derivatives[degree]

    def chain_space(self, degree=None):
        if degree is None:
            return self._chain_spaces
        return self._chain_spaces[degree]

    @cached_method
    def cycle_space(self,degree=None):
        r"""
        The space of cycles.
        """
        if degree is None:
            return [self.cycle_space(i) for i in xrange(3)]
        return self._derivatives[degree].right_kernel()

    @cached_method
    def boundary_space(self,degree=None):
        r"""
        Return the boundary space of the given ``degree``.
        """
        if degree is None:
            return [self.boundary_space(i) for i in xrange(3)]
        if degree == 2:
            return self.chain_space(2).submodule([])
        else:
            return self._derivatives[degree+1].column_space()


