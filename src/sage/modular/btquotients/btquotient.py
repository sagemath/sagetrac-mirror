#########################################################################
#       Copyright (C) 2011 Cameron Franc and Marc Masdeu
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################
r"""
Arithmetic quotients of the Bruhat-Tits tree

"""
from sage.rings.integer import Integer
from sage.structure.element import Element
from sage.matrix.constructor import Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.structure.sage_object import SageObject
from sage.rings.all import ZZ,Zmod,QQ
from sage.misc.latex import latex
from sage.plot import plot
from sage.rings.padics.precision_error import PrecisionError
from itertools import islice
import collections
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.rings.arith import gcd,xgcd,kronecker_symbol
from sage.rings.padics.all import Qp,Zp
from sage.algebras.quatalg.all import QuaternionAlgebra
from sage.quadratic_forms.all import QuadraticForm
from sage.graphs.all import Graph
from sage.libs.all import pari
from sage.interfaces.all import magma
from copy import copy
from sage.plot.colors import rainbow
from sage.rings.number_field.all import NumberField
from sage.modular.arithgroup.all import Gamma0


class DoubleCosetReduction(SageObject):
    r"""
    Edges in the Bruhat-Tits tree are represented as certain double 
    cosets in `\GL_2`. Given a matrix `x` in `\GL_2`, this class computes and stores 
    the data corresponding to its double coset representation. 

    More precisely:
    Initialized with an element `x` of `\GL_2(\ZZ)`, finds elements 
    `\gamma`, `t` and an edge `e` such that `get=x`. It stores these 
    values as members ``gamma``, ``label`` and functions ``self.sign()``, 
    ``self.t()`` and ``self.igamma()``, satisfying:
        if ``self.sign()==+1``:
            ``igamma()*edge_list[label].rep*t()==x``
        if ``self.sign()==-1``:
            ``igamma()*edge_list[label].opposite.rep*t()==x``

    It also stores a member called power so that:
        ``p**(2*power)=gamma.reduced_norm()``
   
    The usual decomposition ``get=x``, with would be:
        g=gamma/(p**power)
        e=edge_list[label]
        t'=t*p**power
    Here usual denotes that we've rescaled gamma to have unit
    determinant, and so that the result is honestly an element 
    of the arithmetic quarternion group under consideration. In 
    practice we store integral multiples and keep track of the 
    powers of `p`.

    INPUT:

    - ``Y`` -  BTQuotient object in which to work

    - ``x`` -  Something coercible into a matrix in `\GL_2(\ZZ)`. In 
       principle we should allow elements in `\GL_2(\QQ_p)`, but for what 
       we do it is enough to work with integral entries

    - ``extrapow`` - gets added to the power attribute, and it is 
       used for the Hecke action.

    EXAMPLES::

        sage: Y = BTQuotient(5,13)
        sage: x = Matrix(ZZ,2,2,[123,153,1231,1231])
        sage: d = DoubleCosetReduction(Y,x)
        sage: d.sign()
        -1
        sage: d.igamma()*Y._edge_list[d.label - len(Y.get_edge_list())].opposite.rep*d.t()==x
        True
        sage: x = Matrix(ZZ,2,2,[1423,113553,11231,12313])
        sage: d = DoubleCosetReduction(Y,x)
        sage: d.sign()
        1
        sage: d.igamma()*Y._edge_list[d.label].rep*d.t()==x
        True

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu

    """
    def __init__(self,Y,x,extrapow=0):
        r"""
        Initializes and computes the reduction as a double coset.
        """
        e1=Y._BT.edge(x)
        try:
            g,label,parity=Y._cached_decomps[e1]
        except KeyError:
            valuation=e1.determinant().valuation(Y._p)
            parity=valuation%2
            v1=Y._BT.target(e1)
            v=Y.fundom_rep(v1)
            g,e=Y._find_equivalent_edge(e1,v.entering_edges,valuation=valuation)
            label=e.label
            Y._cached_decomps[e1]=(g,label,parity)

        self._parent=Y
        self.parity=parity
        self._num_edges = len(Y.get_edge_list())
        self.label=label + parity * self._num_edges # The label will encode whether it is an edge or its opposite !
        self.gamma=g[0]
        self.x=x
        self.power=g[1]+extrapow
        self._t_prec=-1
        self._igamma_prec=-1


    def sign(self):
        r"""
        The direction of the edge.

        The BT quotients are directed graphs but we only store 
        half the edges (we treat them more like unordered graphs). 
        The sign tells whether the matrix self.x is equivalent to the 
        representative in the quotient (sign = +1), or to the
       opposite of one of the representatives (sign = -1).
        """
        if self.parity == 0:
            return 1
        else:
            return -1

    def igamma(self,embedding = None):
        r"""
        Image under gamma.

        Elements of the arithmetic group can be regarded as elements
        of the global quarterion order, and hence may be represented 
        exactly. This function computes the image of such an element 
        under the local splitting and returns the corresponding p-adic
        approximation.

        INPUT:

          - ``embedding`` - an integer, or a function (Default: none). If
            ``embedding`` is None, then the image of ``self.gamma`` under the local splitting
            associated to ``self.Y`` is used. If ``embedding`` is an integer, then the
            precision of the local splitting of self.Y is raised (if necessary)
            to be larger than this integer, and this new local splitting is
            used. If a function is passed, then map ``self.gamma``
            under ``embedding``.
        """
        Y=self._parent
        if embedding is None:
            prec=Y._prec
        else:
            try:
                # The user wants higher precision
                prec=ZZ(embedding)
            except TypeError:
                # The user knows what she is doing, so let it go
                return embedding(self.gamma)
        if prec > self._igamma_prec:
            self._igamma_prec=prec
            self._cached_igamma=Y.embed_quaternion(self.gamma,exact = False, prec = prec)
        return self._cached_igamma

    def t(self, prec = None):
        r"""
        The 't part' of the decomposition using the rest of the data.

        EXAMPLES::

            sage: Y = BTQuotient(5,13)
            sage: x = Matrix(ZZ,2,2,[123,153,1231,1232])
            sage: d = DoubleCosetReduction(Y,x)
            sage: t = d.t(20)
            sage: t[1,0].valuation() > 0
            True
        """
        Y=self._parent
        if prec is None:
            prec = Y._prec
        if(self._t_prec>=prec):
            return self._cached_t #.change_ring(self._parent._R)
        e=Y._edge_list[self.label % self._num_edges]
        self._t_prec=prec
        if(self.parity==0):
            self._cached_t=(self.igamma(prec)*e.rep).inverse()*self.x
        else:
            self._cached_t=(self.igamma(prec)*e.opposite.rep).inverse()*self.x
        return self._cached_t

class BruhatTitsTree(SageObject, UniqueRepresentation):
    r"""
    An implementation of the Bruhat-Tits tree for `\GL_2(\QQ_p)`.

    INPUT:

    - ``p`` - a prime number. The corresponding tree is then p+1 regular 

    EXAMPLES::

    Here we create the tree for `\GL_2(\QQ_5)`:

        sage: p = 5
        sage: T = BruhatTitsTree(p)
        sage: m = Matrix(ZZ,2,2,[p**5,p**2,p**3,1+p+p*3])
        sage: e = T.edge(m); e
        [  0  25]
        [625  21]
        sage: v0 = T.origin(e); v0
        [ 25   0]
        [ 21 125]
        sage: v1 = T.target(e); v1
        [ 25   0]
        [ 21 625]
        sage: T.origin(T.opposite(e)) == v1
        True
        sage: T.target(T.opposite(e)) == v0
        True

    A value error is raised if a prime is not passed:

        sage: T = BruhatTitsTree(4)
        Traceback (most recent call last):
        ...
        ValueError: Input (4) must be prime

    AUTHORS:

    - Marc Masdeu (2012-02-20)
    """
    def __init__(self,p):
        if not(ZZ(p).is_prime()):
            raise ValueError, 'Input (%s) must be prime'%p
        self._p=ZZ(p)
        self._Mat_22=MatrixSpace(ZZ,2,2)
        self._mat_p001=self._Mat_22([self._p,0,0,1])

    def target(self,e,normalized = False):
        r"""
        Returns the target vertex of the edge represented by the
        input matrix M.

        INPUT:

          - ``e`` - a 2x2 matrix with integer entries

          - ``normalized`` - boolean (default: false). If true 
            then the input matrix is assumed to be normalized.

        OUPUT:

          A 2x2 integer matrix

        """
        if normalized:
            #then the normalized target vertex is also M and we save some
            #row reductions with a simple return
            return e
        else:
            #must normalize the target vertex representative
            return self.vertex(e)

    def origin(self, e ,normalized = False):
        r"""
        Returns the origin vertex of the edge represented by the
        input matrix M.

        INPUT:

          - ``e`` - a 2x2 matrix with integer entries

          - ``normalized`` - boolean (default: false). If true 
            then the input matrix M is assumed to be normalized

        OUTPUT:

          A 2x2 integer matrix

        """
        if not normalized:
            #then normalize
            x=copy(self.edge(e))
        else:
            x=copy(M)
        x.swap_columns(0,1)
        x.rescale_col(0,self._p)
        return self.vertex(x)

    def edge(self,M):
        r"""
        Normalizes a matrix to the correct normalized edge
        representative.

        INPUT:

        - ``M`` - a 2x2 integer matrix

        OUTPUT:

          A 2x2 integer matrix

        EXAMPLES::

            sage: T = BruhatTitsTree(3)
            sage: T.edge( Matrix(ZZ,2,2,[0,-1,3,0]) )
            [0 1]
            [3 0]

        """
        p=self._p
        M = M.change_ring(ZZ)
        try:
            v=min([M[i,j].valuation(p) for i in range(2) for j in range(2)])
            val=lambda x:x.valuation(p)
        except TypeError:
            v=min([M[i,j].valuation() for i in range(2) for j in range(2)])
            val=lambda x:x.valuation()

        if v != 0:
            M=p**(-v)*M

        m00=val(M[0,0])
        m01=val(M[0,1])

        if (m00 <= m01):
            tmp=val(M.determinant())-m00
            bigpower=p**(1+tmp)
            r=M[0,0]
            if r!=0:
                r/=p**m00
                assert val(r) == 0
            g,s,_=xgcd(r,bigpower)
            assert g==1
            r=(M[1,0]*s)%bigpower
            newM=self._Mat_22([p**m00,0,r,bigpower/p])
        else:
            tmp=val(M.determinant())-m01
            bigpower=p**tmp
            r=M[0,1]
            if r!=0:
                r/=p**m01
                assert val(r)==0
            g,s,_ = xgcd(r,bigpower)
            assert g == 1
            r=(M[1,1]*s)%bigpower
            newM=self._Mat_22([0,p**m01,bigpower,r])
        newM.set_immutable()
        # assert self.belongs_to_group(newM.inverse()*M,as_edge = True)
        return newM

    def belongs_to_group(self,x,as_edge = False):
        det = x.determinant()
        det_padic = Qp(self._p,20)(det)
        x = (1/self._p**det_padic.valuation(self._p))*x
        assert x.determinant().valuation(self._p)==0
        v = min([xx.valuation(self._p) for xx in x.list()])
        if v<0:
            return False
        if as_edge == True:
            if x[1,0].valuation(self._p)>0:
                return True
            else:
                return False
        else:
            return True

    def vertex(self,M):
        r"""
        Normalizes a matrix to the corresponding normalized
        vertex representative

        INPUT:

        - ``M`` - 2x2 integer matrix

        OUTPUT:

        -  a 2x2 integer matrix

        EXAMPLES::

            sage: p = 5
            sage: T = BruhatTitsTree(p)
            sage: m = Matrix(ZZ,2,2,[p**5,p**2,p**3,1+p+p*3])
            sage: e = T.edge(m)
            sage: t = m.inverse()*e
            sage: scaling = Qp(p,20)(t.determinant()).sqrt()
            sage: t = 1/scaling * t
            sage: min([t[ii,jj].valuation(p) for ii in range(2) for jj in range(2)]) >= 0
            True
            sage: t[1,0].valuation(p) > 0
            True
        """
        p=self._p
        try:
            v=min([M[i,j].valuation(p) for i in range(2) for j in range(2)])
            val=lambda x:x.valuation(p)
        except TypeError:
            v=min([M[i,j].valuation() for i in range(2) for j in range(2)])
            val=lambda x:x.valuation()
        if v != 0:
            M=p**(-v)*M
        m00=val(M[0,0])
        m01=val(M[0,1])
        if m01<m00:
            M=copy(M)
            M.swap_columns(0,1)
            m00=m01
        m10=val(M[1,0])
        tmp=val(M.determinant())-m00
        bigpower=p**tmp
        r=M[0,0]
        if r!=0:
            r/=p**m00
            assert val(r) == 0
        g,s,_=xgcd(r,bigpower)
        assert g == 1
        r=(M[1,0]*s)%bigpower
        newM=self._Mat_22([p**m00,0,r,bigpower])
        newM.set_immutable()
        # assert self.belongs_to_group(newM.inverse()*M,as_edge = False)
        return newM


    def edges_leaving_origin(self):
        r"""
        Find normalized representatives for the `p+1` edges leaving the origin vertex 
        corresponding to the homothety class of `\ZZ_p^2`. These are cached.

        OUTPUT:

        -  A list of size `p+1` of 2x2 integer matrices

        EXAMPLES::

            sage: T = BruhatTitsTree(3)
            sage: T.edges_leaving_origin()
            [
            [0 1]  [3 0]  [0 1]  [0 1]
            [3 0], [0 1], [3 1], [3 2]
            ]

        """
        try: return self._edges_leaving_origin
        except:
            p=self._p
            self._edges_leaving_origin=[self.edge(self._Mat_22([0,-1,p,0]))]
            self._edges_leaving_origin.extend([self.edge(self._Mat_22([p,i,0,1])) for i in range(p)])
            return self._edges_leaving_origin

    def edge_between_vertices(self,v1,v2, normalized = False):
        r"""
        This function computes the normalized matrix rep. for the edge
        passing between two vertices.

        INPUT:

        - ``v1`` - 2x2 integer matrix
        - ``v2`` - 2x2 integer matrix
        - ``normalized`` - boolean (Default: False) Whether the vertices are normalized.

        OUTPUT:

        -  2x2 integer matrix, representing the edge from ``v1`` to ``v2``.
          If ``v1`` and ``v2`` are not at distance `1`, raise a ``ValueError``.

        EXAMPLES::

            sage: p = 7
            sage: T = BruhatTitsTree(p)
            sage: v1 = T.vertex(Matrix(ZZ,2,2,[p,0,0,1])); v1
            [7 0]
            [0 1]
            sage: v2 = T.vertex(Matrix(ZZ,2,2,[p,1,0,1])); v2
            [1 0]
            [1 7]
            sage: T.edge_between_vertices(v1,v2)
            Traceback (most recent call last):
            ...
            ValueError: Vertices are not adjacent.

            sage: v3 = T.vertex(Matrix(ZZ,2,2,[1,0,0,1])); v3
            [1 0]
            [0 1]
            sage: T.edge_between_vertices(v1,v3)
            [0 1]
            [1 0]
        """
        if normalized:
            v22=v2
        else:
            v22=self.vertex(v2)
        for e in self.leaving_edges(v1):
            if self.target(e)==v22:
                return e
        raise ValueError, 'Vertices are not adjacent.'

    def leaving_edges(self,M):
        r"""
        Edges leaving a vertex

        INPUT:

        - ``M`` - 2x2 integer matrix

        OUTPUT:

          List of size p+1 of 2x2 integer matrices

        EXAMPLES::

            sage: p = 7
            sage: T = BruhatTitsTree(p)
            sage: T.leaving_edges(Matrix(ZZ,2,2,[1,0,0,1]))
            [
            [0 1]  [7 0]  [0 1]  [0 1]  [0 1]  [0 1]  [0 1]  [0 1]
            [7 0], [0 1], [7 1], [7 4], [7 5], [7 2], [7 3], [7 6]
            ]

        """
        return [self.edge(M*A) for A in self.edges_leaving_origin()]

    def opposite(self,e):
        r"""
        This function returns the edge oriented oppositely to a given edge.

        INPUT:

        - ``e`` - 2x2 integer matrix

        OUPUT:

          2x2 integer matrix

        EXAMPLES::

            sage: p = 7
            sage: T = BruhatTitsTree(p)
            sage: e = Matrix(ZZ,2,2,[1,0,0,1])
            sage: T.opposite(e)
            [0 1]
            [7 0]
            sage: T.opposite(T.opposite(e)) == e
            True

        """
        x=copy(e)
        x.swap_columns(0,1)
        x.rescale_col(0,self._p)
        return self.edge(x)

    def entering_edges(self,v):
        r"""
        This function returns the edges entering a given vertex.

        INPUT:

        - ``v`` - 2x2 integer matrix

        OUTPUT:

          A list of size p+1 of 2x2 integer matrices

        EXAMPLES::

            sage: p = 7
            sage: T = BruhatTitsTree(p)
            sage: T.entering_edges(Matrix(ZZ,2,2,[1,0,0,1]))
            [
            [1 0]  [0 1]  [1 0]  [1 0]  [1 0]  [1 0]  [1 0]  [1 0]
            [0 1], [1 0], [1 1], [4 1], [5 1], [2 1], [3 1], [6 1]
            ]

        """
        return [self.opposite(e) for e in self.leaving_edges(v)]

    def subdivide(self,edgelist,level):
        r"""
        (Ordered) edges of self may be regarded as open balls in P_1(Qp).
        Given a list of edges, this function return a list of edges
        corresponding to the level-th subdivision of the corresponding opens. 
        That is, each open ball of the input is broken up into `p^\mbox{level}`
        subballs of equal radius.

        INPUT:

          - ``edgelist`` - a list of edges

          - ``level`` - an integer

        OUTPUT:

          A list of 2x2 integer matrices

        EXAMPLES::

            sage: p = 3
            sage: T = BruhatTitsTree(p)
            sage: T.subdivide([Matrix(ZZ,2,2,[p,0,0,1])],2)
            [
            [27  0]  [0 9]  [0 9]  [0 3]  [0 3]  [0 3]  [0 3]  [0 3]  [0 3]
            [ 0  1], [3 1], [3 2], [9 1], [9 4], [9 7], [9 2], [9 5], [9 8]
            ]

        """
        all_edges=[]
        if(level<0):
            return []
        if(level==0):
            return [self._Mat_22(edge) for edge in edgelist]
        else:
            newEgood=[]
            for edge in edgelist:
                edge=self._Mat_22(edge)
                origin=self.origin(edge)
                newE=self.leaving_edges(self.target(edge))
                newEgood.extend([e for e in newE if self.target(e)!=origin])
            return self.subdivide(newEgood,level-1)

    def get_balls(self,center=1,level=1):
        r"""
        Returns a decomposition of `\PP^1(\QQ_p)` into compact open balls.

        Each vertex in the Bruhat-Tits tree gives a decomposition of
        `\PP^1(\QQ_p)` into `p+1` open balls. Each of these balls may
        be further subdivided, to get a finer decomposition.

        This function returns the decompostion of `\PP^1(\QQ_p)` corresponding
        to ``center`` into `(p+1)p^\mbox{level}` balls.

        EXAMPLES::

            sage: p = 2
            sage: T = BruhatTitsTree(p)
            sage: T.get_balls(Matrix(ZZ,2,2,[p,0,0,1]),1)
            [
            [0 1]  [0 1]  [8 0]  [0 4]  [0 2]  [0 2]
            [2 0], [2 1], [0 1], [2 1], [4 1], [4 3]
            ]

        """
        return self.subdivide(self.leaving_edges(center),level)

    def find_path(self,v,boundary=None):
        r"""
        Computes a path from a vertex to a given set of so-called
        boundary vertices, whose interior must contain the origin vertex.
        In the case that the boundary is not specified, it
        computes the geodesic between the given vertex and the origin.
        In the case that the boundary contains more than one vertex,
        it computes the geodesic to some point of the boundary.

        INPUT:

        - ``v`` - a 2x2 matrix representing a vertex

        - ``boundary`` - a list of matrices (default: None). If ommitted, finds
          the geodesic from ``v`` to the central vertex.

        OUTPUT:

        An ordered list of edges describing the geodesic from ``v`` to
        ``boundary``, followed by the vertex in the boundary that is
        closest to ``v``.

        EXAMPLES::

            sage: p = 3
            sage: T = BruhatTitsTree(p)
            sage: T.find_path( Matrix(ZZ,2,2,[p^4,0,0,1]) )
            ([[81  0]
            [ 0  1], [27  0]
            [ 0  1], [9 0]
            [0 1], [3 0]
            [0 1]], [1 0]
            [0 1])
            sage: T.find_path( Matrix(ZZ,2,2,[p^3,0,134,p^2]) )
            ([[27  0]
            [ 8  9], [27  0]
            [ 2  3], [27  0]
            [ 0  1], [9 0]
            [0 1], [3 0]
            [0 1]], [1 0]
            [0 1])

        """
        if boundary is None:
            m=self._Mat_22(1)
            m.set_immutable()
            boundary = {m:m}
        m=self._mat_p001
        new_v=self.vertex(v)
        chain=[]
        while new_v[1,0]!=0 or new_v[0,0].valuation(self._p)<new_v[1,1].valuation(self._p):
            if boundary.has_key(new_v):
                return chain,boundary[new_v]
            chain.append(new_v)
            new_v=self.vertex(new_v*m)

        if boundary.has_key(new_v):
            return chain,boundary[new_v]

        while True:
            if boundary.has_key(new_v):
                return chain,boundary[new_v]
            chain.append(new_v)
            new_v=self._Mat_22([new_v[0,0]/self._p,0,0,1])
            new_v.set_immutable()
        raise RuntimeError

    def find_containing_affinoid(self,z):
        r"""
        Returns the vertex corresponding to the affinoid in 
        the `p`-adic upper half plane that a given (unramified!) point reduces to.

        INPUT:

          - ``z`` - an element of an unramified extension of `\QQ_p` that is not contained
            in `\QQ_p`.

        OUTPUT:

          A 2x2 integer matrix representing a vertex of ``self``.

        EXAMPLES::

            sage: T = BruhatTitsTree(5)
            sage: K.<a> = Qq(5^2,20)
            sage: T.find_containing_affinoid(a)
            [1 0]
            [0 1]
            sage: z = 5*a+3
            sage: v = T.find_containing_affinoid(z).inverse(); v
            [   1    0]
            [-2/5  1/5]

        Note that the translate of ``z`` belongs to the standard affinoid. That is,
        it is a `p`-adic unit and its reduction modulo `p` is not in `\FF_p`::

            sage: gz = (v[0,0]*z+v[0,1])/(v[1,0]*z+v[1,1]); gz
            (a + 1) + O(5^19)
            sage: gz.valuation() == 0
            True
        """
        #Assume z belongs to some extension of QQp.
        p=self._p
        if(z.valuation()<0):
            return self.vertex(self._Mat_22([0,1,p,0])*self.find_containing_affinoid(1/(p*z)))
        a=0
        pn=1
        val=z.valuation()
        L=[]
        for ii in range(val):
            L.append(0)
        L.extend(z.list())
        for n in range(len(L)):
            if(L[n]!=0):
                if(len(L[n])>1):
                    break
                if(len(L[n])>0):
                    a+=pn*L[n][0]
            pn*=p
        return self.vertex(self._Mat_22([pn,a,0,1]))

    def find_geodesic(self,v1,v2,normalized = True):
        r"""
        This function computes the geodesic between two vertices

        INPUT:

          - ``v1`` - 2x2 integer matrix representing a vertex

          - ``v2`` - 2x2 integer matrix representing a vertex

          - ``normalized`` - boolean (Default: True)

        OUTPUT:

          ordered list of 2x2 integer matrices representing edges

        EXAMPLES::

            sage: p = 3
            sage: T = BruhatTitsTree(p)
            sage: v1 = T.vertex( Matrix(ZZ,2,2,[p^3, 0, 1, p^1]) ); v1
            [27  0]
            [ 1  3]
            sage: v2 = T.vertex( Matrix(ZZ,2,2,[p,2,0,p]) ); v2
            [1 0]
            [6 9]
            sage: T.find_geodesic(v1,v2)
            [
            [27  0]  [27  0]  [9 0]  [3 0]  [1 0]  [1 0]  [1 0]
            [ 1  3], [ 0  1], [0 1], [0 1], [0 1], [0 3], [6 9]
            ]

        """
        if not normalized:
            v1,v2=self.vertex(v1),self.vertex(v2)
        gamma=v2
        vv=self.vertex(gamma.adjoint()*v1)
        chain,v0=self.find_path(vv)
        chain.append(v0)
        return [self.vertex(gamma*x) for x in chain]

    def find_covering(self,z1,z2):
        r"""
        This function computes a covering of P1(Qp) adapted to a 
        certain geodesic in self.

        More precisely, the `p`-adic upper half plane points ``z1``
        and ``z2`` reduce to vertices `v_1`, `v_2`.
        The returned covering consists of all the edges leaving the
        geodesic from `v_1` to `v_2`.

        INPUT:

          - ``z1``, ``z2`` - unramified algebraic points of h_p

        OUTPUT:

          a list of 2x2 integer matrices representing edges of self

        EXAMPLES::

            sage: p = 3
            sage: K.<a> = Qq(p^2)
            sage: T = BruhatTitsTree(p)
            sage: z1 = a + a*p
            sage: z2 = 1 + a*p + a*p^2 - p^6
            sage: T.find_covering(z1,z2)
            [
            [0 1]  [3 0]  [0 1]  [0 1]  [0 1]  [0 1]
            [3 0], [0 1], [3 2], [9 1], [9 4], [9 7]
            ]

        NOTES::

          This function is used to compute certain Coleman integrals on `\PP^1`. That's
          why the input consists of two points of the `p`-adic upper half plane, but
          decomposes `\PP^1(\QQ_p)`. This decomposition is what allows us to represent
          the relevant integrand as a locally analytic function. The ``z1`` and ``z2``
          appear in the integrand.

        """
        v1=self.find_containing_affinoid(z1)
        v2=self.find_containing_affinoid(z2)
        vertex_set=[self._Mat_22(0)]+self.find_geodesic(v1,v2)+[self._Mat_22(0)]
        E=[]
        for ii in range(1,len(vertex_set)-1):
            vv=vertex_set[ii]
            newE=self.leaving_edges(vv)
            for e in newE:
                targ=self.target(e)
                if targ!=vertex_set[ii-1] and targ!=vertex_set[ii+1]:
                    E.append(e)
        return E


class Vertex(SageObject):
    r"""
    This is a structure to represent vertices of quotients of the Bruhat-Tits tree. 
    It is useful to enrich the representation of the vertex as a matrix with extra
    data.

    INPUT:

     - ``label`` - An integer which uniquely identifies this vertex.
     - ``rep`` - A 2x2 matrix in reduced form representing this vertex.
     - ``leaving_edges`` - (Default: empty list) A list of edges leaving this vertex.
     - ``entering_edges`` - (Default: empty list) A list of edges entering this vertex.
     - ``determinant`` - (Default: None) The determinant of ``rep``, if known.
     - ``valuation`` - (Default: None) The valuation of the determinant of ``rep``, if known.


    AUTHORS:

    - Marc Masdeu (2012-02-20)
    """
    def __init__(self,owner,label,rep,leaving_edges=None,entering_edges=None,determinant=None,valuation=None):
        if leaving_edges is None:
            leaving_edges = []
        if entering_edges is None:
            entering_edges = []
        if determinant is None:
            determinant = rep.determinant()
        if valuation is None:
            valuation = determinant.valuation(owner._p)
        self.owner=owner
        self.label=label
        self.rep=rep
        self.rep.set_immutable()
        self.determinant=determinant
        self.valuation=valuation
        self.parity=valuation%2
        self.leaving_edges=leaving_edges
        self.entering_edges=entering_edges


class Edge(SageObject):
    r"""
    This is a structure to represent edges of quotients of the Bruhat-Tits tree. 
    It is useful to enrich the representation of an edge as a matrix with extra
    data.

    INPUT:

     - ``label`` - An integer which uniquely identifies this edge.
     - ``rep`` - A 2x2 matrix in reduced form representing this edge.
     - ``origin`` - The origin vertex of ``self``.
     - ``target`` - The target vertex of ``self``.
     - ``links`` - (Default: empty list) A list of elements of `\Gamma` which identify different
     edges in the Bruhat-Tits tree which are equivalent to ``self``.
     - ``opposite`` - (Default: None) The edge opposite to ``self``
     - ``determinant`` - (Default: None) The determinant of ``rep``, if known.
     - ``valuation`` - (Default: None) The valuation of the determinant of ``rep``, if known.


    AUTHORS:

    - Marc Masdeu (2012-02-20)
    """
    def __init__(self,owner,label,rep,origin,target,links = None,opposite = None,determinant = None,valuation = None):
        if links is None:
            links = []
        if determinant is None:
            determinant=rep.determinant()
        if valuation is None:
            valuation = determinant.valuation(owner._p)
        self.owner=owner
        self.label=label
        self.rep=rep
        self.rep.set_immutable()
        self.origin=origin
        self.target=target
        self.links=links
        self.opposite=opposite
        self.determinant=determinant
        self.valuation=valuation
        self.parity=valuation%2

class BTQuotient(SageObject, UniqueRepresentation):
    @staticmethod
    def __classcall__(cls,p,Nminus,Nplus=1,use_magma = False):
        return super(BTQuotient,cls).__classcall__(cls,p,Nminus,Nplus,use_magma)

    r"""
    This function computes the quotient of the Bruhat-Tits tree by an arithmetic
    quaternionic group. The group in question is the group of norm 1 elements in
    an eichler Z[1/p]-order of some (tame) level inside of a definite quaternion 
    algebra that is unramified at the prime p.

    INPUT:

     - ``p`` - a prime number
     - ``Nminus`` - squarefree integer divisible by an odd number of distinct primes and 
       relatively prime to p. This is the discriminant of the definite quaternion algebra 
       that one is quotienting by.
     - ``Nplus`` - an integer corpime to pNminus (Default: 1). This is the tame level. It need
       not be squarefree! If Nplus is not 1 then the user currently needs magma installed due
       to sage's inability to compute well with nonmaximal Eichler orders in rational (definite) 
       quaternion algebras.
     - ``use_magma`` - boolean (default: False). If True, uses magma for quaternion arithmetic.

    EXAMPLES::

        sage: X = BTQuotient(13,19)
        sage: X.genus()
        19
        sage: G = X.get_graph(); G
        Multi-graph on 4 vertices

    NOTES::

      A sage implementation of Eichler orders in rational quaternions algebras would remove the
      dependency on magma.

    AUTHORS::

    - Marc Masdeu (2012-02-20)
    """
    def __init__(self,p,Nminus,Nplus=1,use_magma = False):
        Nminus=Integer(Nminus)
        Nplus=Integer(Nplus)
        p=Integer(p)
        lev=p*Nminus
        if not p.is_prime():
            raise ValueError, "p must be a prime"
        if not lev.is_squarefree():
            raise ValueError, "level must be squarefree"
        if(gcd(lev,Nplus)>1):
            raise ValueError, "level and conductor must be coprime"

        if len(Nminus.factor())%2 != 1:
            raise ValueError, "Nminus should be divisible by an odd number of primes"

        self._pN=p
        self._p=p
        self._Nminus=Nminus
        self._Nplus=Nplus
        if(use_magma==True or not self._Nminus.is_prime() or self._Nplus!=1 or self._p==2):
            try:
                self._magma=magma
                magmap=self._magma(p)
                # print "Warning: this input needs magma to work..."
            except RuntimeError:
                raise NotImplementedError,'Sage does not know yet how to work with the kind of orders that you are trying to use. Try installing Magma first and set it up so that Sage can use it.'
            self._use_magma = True
        else:
            self._use_magma = False

        self._BT=BruhatTitsTree(p)
        self._prec=-1
        self._cached_vertices=dict()
        self._cached_edges=dict()
        self._cached_paths=dict()
        self._cached_decomps=dict()
        self._cached_equivalent=dict()
        self._CM_points=dict()

        self._V=(QQ**4).ambient_module().change_ring(ZZ)
        self._Mat_44=MatrixSpace(ZZ,4,4)
        self._Mat_22=MatrixSpace(ZZ,2,2)
        self._Mat_41=MatrixSpace(ZZ,4,1)
        self._Xv=[self._Mat_22([1,0,0,0]),self._Mat_22([0,1,0,0]),self._Mat_22([0,0,1,0]),self._Mat_22([0,0,0,1])]
        self._Xe=[self._Mat_22([1,0,0,0]),self._Mat_22([0,1,0,0]),self._Mat_22([0,0,self._p,0]),self._Mat_22([0,0,0,1])]

    def _repr_(self):
        r"""
        This returns the representation of self as a string.

        EXAMPLES::

            sage: X = BTQuotient(5,13); X
            Quotient of the Bruhat Tits tree of GL_2(QQ_5) with discriminant 13 and level 1

        """
        return "Quotient of the Bruhat Tits tree of GL_2(QQ_%s) with discriminant %s and level %s"%(self.prime(),self.Nminus().factor(),self.Nplus().factor())

    def _latex_(self):
        r"""
        Returns the LaTeX representation of self.

        EXAMPLES::

            sage: X = BTQuotient(5,13); latex(X)
            X(5 \cdot 13,1)\otimes_{\mathbb{Z}} \mathbb{F}_{5}

        """
        return "X(%s,%s)\\otimes_{\\mathbb{Z}} \\mathbb{F}_{%s}"%(latex(self.level().factor()),latex(self.Nplus().factor()),latex(self.prime()))

    def get_vertex_dict(self):
        r"""
        This function returns the vertices of the quotient viewed as
        a dict.

        OUTPUT:

        A python dict with the vertices of the quotient.

        """
        try: return self._boundary
        except AttributeError:
            self._compute_quotient()
            return self._boundary

    def get_vertex_list(self):
        r"""
        This function returns a list of the vertices of the quotient.

        OUTPUT:

        A list with the vertices of the quotient.

        """
        try: return self._vertex_list
        except AttributeError:
            self._compute_quotient()
            return self._vertex_list

    def get_edge_list(self):
        r"""
        This function returns a list of ``Edge``s which represent a fundamental domain
        inside the Bruhat-Tits tree for the quotient.

        OUTPUT:

          A list of ``Edge``s.

        """
        try: return self._edge_list
        except AttributeError:
            self._compute_quotient()
            return self._edge_list

    def get_list(self):
        r"""
        This function returns a list of ``Edge``s which represent a fundamental domain
        inside the Bruhat-Tits tree for the quotient, together with a list of the
        opposite edges. This is used to work with automorphic forms.

        OUTPUT:

          A list of ``Edge``s.

        """
        E = self.get_edge_list()
        return E + [e.opposite for e in E]

    def get_generators(self):
        r"""
        This function uses a fundamental domain in the Bruhat-Tits tree, and certain gluing
        data for boundary vertices, in order to compute a collection of generators for the
        arithmetic quaternionic group that one is quotienting by. This is analogous to using
        a polygonal rep. of a compact real surface to present its fundamental domain.

        OUTPUT:

          A list of elements of an arithmetic quaternionic group.

        EXAMPLES::

            sage: X = BTQuotient(3,13)
            sage: X.get_generators()
            [
            [ 2]  [-5]  [ 4]
            [-5]  [ 3]  [ 1]
            [ 1]  [ 1]  [-3]
            [ 0], [ 2], [-2]
            ]

        """
        try: return list(self._generators)
        except AttributeError:
            self._compute_quotient()
            return list(self._generators)

    @cached_method
    def genus(self):
        r"""
        This function computes the genus of the Shimura curve corresponding to this
        quotient via Cerednik-Drinfeld. It is computed via a formula and not in terms
        of the quotient graph.

        INPUT:

        - level: Integer (default: None) a level. By default, use that of ``self``.
        - Nplus: Integer (default: None) a conductor. By default, use that of ``self``.

        OUTPUT:

          An integer equal to the genus

        EXAMPLES::

            sage: X = BTQuotient(7,19)
            sage: X.genus()
            9
        """
        return self.dimension_harmonic_cocycles(2)

    @cached_method
    def dimension_harmonic_cocycles(self,k,lev = None,Nplus = None):
        r"""
        This function computes the dimension of the space of harmonic cocycles
        of weight `k` on ``self``.

        OUTPUT:

          An integer equal to the dimension

        EXAMPLES::

            sage: X = BTQuotient(3,7)
            sage: print [X.dimension_harmonic_cocycles(k) for k in range(2,20,2)]
            [1, 4, 4, 8, 8, 12, 12, 16, 16]
            sage: print [len(HarmonicCocycles(X,k,100).basis()) for k in range(2,20,2)] # long time
            [1, 4, 4, 8, 8, 12, 12, 16, 16]

            sage: X = BTQuotient(2,5)
            sage: print [X.dimension_harmonic_cocycles(k) for k in range(2,40,2)]
            [0, 1, 3, 1, 3, 5, 3, 5, 7, 5, 7, 9, 7, 9, 11, 9, 11, 13, 11]
            sage: print [len(HarmonicCocycles(X,k,100).basis()) for k in range(2,40,2)] # long time
            [0, 1, 3, 1, 3, 5, 3, 5, 7, 5, 7, 9, 7, 9, 11, 9, 11, 13, 11]

        """

        k = ZZ(k)
        if lev is None:
            lev = self._p * self._Nminus
        else:
            lev = ZZ(lev)
        if Nplus is None:
            Nplus = self._Nplus
        else:
            Nplus = ZZ(Nplus)

        if k == 0:
            return 0

        if lev == 1:
            return Gamma0(Nplus).dimension_cusp_forms(k = k)

        f = lev.factor()
        if any([l[1] != 1 for l in f]):
            raise NotImplementedError, 'The level should be squarefree for this function to work... Sorry!'

        divs = lev.divisors()

        return Gamma0(lev*Nplus).dimension_cusp_forms(k = k) - sum([len(ZZ(lev/d).divisors())*self.dimension_harmonic_cocycles(k,d,Nplus) for d in divs[:len(divs)-1]])

    def Nplus(self):
        r"""
        This function returns the tame level `N^+`.

        OUTPUT:

          An integer equal to `N^+`.

        EXAMPLES::

            sage: X = BTQuotient(5,7,1)
            sage: X.Nplus()
            1

        """
        return self._Nplus


    def Nminus(self):
        r"""
        This function returns the discriminant of the relevant definite
        quaternion algebra.

        OUTPUT:

          An integer equal to `N^-`.

        EXAMPLES::

            sage: X = BTQuotient(5,7)
            sage: X.Nminus()
            7

        """
        return self._Nminus

    @cached_method
    def level(self):
        r"""
        This function returns `p N^-`, which is the discriminant of the indefinite quaternion
        algebra that is uniformed by Cerednik-Drinfeld.

        OUTPUT:

          An integer equal to `p N^-`.

        EXAMPLES::

            sage: X = BTQuotient(5,7)
            sage: X.level()
            35

        """
        return self._Nminus*self._p

    def prime(self):
        r"""
        This function returns the prime one is working with.

        OUTPUT:

          An integer equal to the fixed prime p

        EXAMPLES::

            sage: X = BTQuotient(5,7)
            sage: X.prime()
            5

        """
        return self._p


    def get_graph(self):
        r"""
        This returns the quotient graph (and computes it if needed).

        OUTPUT:

          A graph representing the quotient of the Bruhat-Tits tree.

        EXAMPLES::

            sage: X = BTQuotient(11,5)
            sage: X.get_graph()
            Multi-graph on 2 vertices

        """
        try: return self._S
        except AttributeError:
            self._compute_quotient()
            return self._S

    def get_fundom_graph(self):
        r"""
        This returns the fundamental domain (and computes it if needed).

        OUTPUT:

          A fundamental domain for the action of `\Gamma`.

        EXAMPLES::

            sage: X = BTQuotient(11,5)
            sage: X.get_fundom_graph()
            Graph on 24 vertices

        """
        try: return self._Sfun
        except AttributeError:
            self._compute_quotient()
            return self._Sfun

    def plot(self,*args,**kwargs):
        r"""
        This function plots the quotient graph.

        OUTPUT:

          A plot of the quotient graph

        """
        S=self.get_graph()
        vertex_colors = {}
        v0 = Matrix(ZZ,2,2,[1,0,0,1])
        v0.set_immutable()
        rainbow_color = rainbow(len(self._vertex_list))
        for v in S.vertex_iterator():
            key =rainbow_color[S.get_vertex(v).label]
            if vertex_colors.has_key(key):
                vertex_colors[key].append(v)
            else:
                vertex_colors[key]=[v]

        my_args = dict()
        my_args['vertex_colors'] = vertex_colors
        my_args['color_by_label'] = True
        my_args['vertex_labels'] = False
        my_args.update(kwargs)
        return S.plot(*args,**my_args)
        return S.plot(*args,**kwargs)

    def plot_fundom(self,*args,**kwargs):
        r"""
        This function plots a fundamental domain.

        OUTPUT:

          A plot of the fundamental domain.

        """
        S=self.get_fundom_graph()
        vertex_colors = {}
        rainbow_color = rainbow(len(self._vertex_list))
        for v in S.vertex_iterator():
            key =rainbow_color[S.get_vertex(v).label]
            if vertex_colors.has_key(key):
                vertex_colors[key].append(v)
            else:
                vertex_colors[key]=[v]

        my_args = dict()
        my_args['vertex_colors'] = vertex_colors
        my_args['color_by_label'] = True
        my_args['vertex_labels'] = True
        my_args.update(kwargs)
        return S.plot(*args,**my_args)

    def is_admissible(self,D):
        r"""
        This function tests whether the imaginary quadratic field of
        discriminant `D` embeds in the quaternion algebra. It furthermore
        tests the Heegner hypothesis in this setting (e.g., is `p` inert
        in the field, etc).

        INPUT:

        - ``D`` - an integer whose squarefree part will define the quadratic field

        OUTPUT:

          A boolean describing whether the quadratic field is admissible

        EXAMPLES::

            sage: X = BTQuotient(5,7)
            sage: print [X.is_admissible(D) for D in range(-1,-20,-1)]
            [False, True, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, True, False]

        """
        D=Integer(D).squarefree_part()
        disc=D
        if(D%4!=1):
            disc*=4
        ff=self.level().factor()
        for f in ff:
            if kronecker_symbol(disc,f[0])!=-1:
                return False
        ff=self._Nplus.factor()
        for f in ff:
            if kronecker_symbol(disc,f[0])!=1:
                return False
        return True

    def _local_splitting_map(self,prec):
        r"""
        Returns an embedding of the definite quaternion algebra
        into the algebra of 2x2 matrices with coefficients in `\QQ_p`.

        INPUT:

            prec -- Integer. The precision of the splitting.

        OUTPUT:

            A function giving the splitting.

        EXAMPLES::

            sage: X = BTQuotient(11,3)
            sage: phi = X._local_splitting_map(10)
            sage: B.<i,j,k> = QuaternionAlgebra(3)
            sage: phi(i)**2 == QQ(i**2)*phi(B(1))
            True

        """
        I,J,K=self._local_splitting(prec)
        def phi(q):
            R=I.parent()
            v=q.coefficient_tuple()
            return R(v[0] + I*v[1] + J*v[2] + K*v[3])
        return phi

    def _local_splitting(self,prec):
        r"""
        Finds an embedding of the definite quaternion algebra
        into the algebra of 2x2 matrices with coefficients in `\QQ_p`.

        INPUT:

        - prec -- Integer. The precision of the splitting.

        OUTPUT:

        - Matrices I, J, K giving the splitting.

        EXAMPLES::

            sage: X = BTQuotient(11,3)
            sage: phi = X._local_splitting_map(10)
            sage: B.<i,j,k> = QuaternionAlgebra(3)
            sage: phi(i)**2 == QQ(i**2)*phi(B(1))
            True

        """
        assert(self._use_magma == False)
        if(prec<=self._prec):
            return self._II,self._JJ,self._KK
        self._prec=prec
        A=self.get_quaternion_algebra()

        ZZp=Zp(self._p,prec)
        v=A.invariants()
        a =ZZp(v[0])
        b = ZZp(v[1])
        if (A.base_ring() != QQ):
            raise ValueError, "must be rational quaternion algebra"
        if (A.discriminant() % self._p == 0):
            raise ValueError, "p (=%s) must be an unramified prime"%self._p
        M = MatrixSpace(ZZp, 2)

        if a.is_square():
            alpha=a.sqrt()
            self._II=M([alpha,0,2*alpha,-alpha])
            self._JJ=M([b,-b,b-1,-b])
        else:
            self._II = M([0,a,1,0])
            z=0
            self._JJ=0
            while(self._JJ==0):
                c=a*z*z+b
                if c.is_square():
                    x=c.sqrt()
                    self._JJ=M([x,-a*z,z,-x])
                else:
                    z+=1
        self._KK = self._II*self._JJ
        return self._II, self._JJ, self._KK

    def _compute_embedding_matrix(self,prec):
        r"""
        Returns a matrix representing the embedding with the given precision.

        INPUT:

            - prec - Integer. The precision of the embedding matrix.
        EXAMPLES:

        Note that the entries of the matrix are elements of Zmod:

            sage: X = BTQuotient(3,7)
            sage: A = X._compute_embedding_matrix(10); A
            [26830 29524 53659 59048]
            [29525 26829     1 53659]
            [29525 26830     1 53659]
            [32220 29525  5390     1]
            sage: R = A.base_ring()
            sage: B = X.get_eichler_order_basis()
            sage: R(B[0,0].reduced_trace()) == A[0,0]+A[3,0]
            True

        """
        if self._use_magma == True:
            try: return Matrix(Zmod(self._pN),4,4,self._cached_Iota0_matrix)
            except AttributeError: pass
            Ord=self.get_eichler_order()
            M,f,rho=self._magma.function_call('pMatrixRing',args=[Ord,self._p],params={'Precision':2000},nvals=3)
            OBasis=Ord.Basis()
            v=[f.Image(OBasis[i]) for i in [1,2,3,4]]
            self._cached_Iota0_matrix=[v[kk][ii,jj].sage() for ii in range(1,3) for jj in range(1,3) for kk in range(4)]
            return Matrix(Zmod(self._pN),4,4,self._cached_Iota0_matrix)
        else:
            phi=self._local_splitting_map(prec)
            B=self.get_eichler_order_basis()
            return Matrix(Zmod(self._p**prec),4,4,[phi(B[kk,0])[ii,jj] for ii in range(2) for jj in range(2) for kk in range(4)])

    def _increase_precision(self,amount=1):
        r"""
        Increase the working precision.

        INPUT:

           - ``amount`` Integer (Default: 1). The amount by which to increase the precision.

        """
        if amount >= 1:
            self.get_embedding_matrix(prec = self._prec+amount)
            return
        else:
            return

    def get_embedding_matrix(self, exact = False, prec = None):
        r"""
        Returns the matrix of the embedding.

        INPUT:

        - ``exact`` boolean (Default: False). If True, return an embedding into a matrix algebra with coefficients
        in a number field. Otherwise, embed into matrices over `p`-adic numbers.

        - ``prec`` Integer (Default: None). If specified, return the marix with precision ``prec``. Otherwise,
          return the the cached matrix (with the current working precision).

        OUTPUT:

        - A 4x4 matrix representing the embedding.
        """
        if exact:
            try:
                return self._Iota_exact
            except:
                raise RuntimeError, 'Exact splitting not available.'
        else:
            if(prec is None or prec is self._prec):
                try:
                    return self._Iota
                except AttributeError: pass

            self._pN=self._p**prec
            self._R=Qp(self._p,prec=prec)

            if(prec>self._prec):
                Iotamod=self._compute_embedding_matrix(prec)
                self._Iotainv_lift=Iotamod.inverse().lift()
                self._Iota=Matrix(self._R,4,4,[Iotamod[ii,jj] for ii in range(4) for jj in range(4)])

            self._prec=prec
            self._Iotainv=self._Mat_44([self._Iotainv_lift[ii,jj]%self._pN for ii in range(4) for jj in range(4)])
            return self._Iota

    def embed_quaternion(self, g, exact = False, prec=None):
        r"""
        Embeds the quaternion element ``g`` into a matrix algebra.

        INPUT:

          - ``g`` a column vector of size `4` whose entries represent a quaternion in our basis.
          - ``exact`` boolean (Default: False) - If True, tries to embed ``g`` into a matrix
          algebra over a number field. If False, the target is the matrix algebra over `\QQ_p`.

        OUTPUT:

          A 2x2 matrix with coefficients in `\QQ_p` if ``exact`` is False, or a number field
          if ``exact`` is True.

        """
        if exact==True:
            return Matrix(self.get_splitting_field(),2,2,(self.get_embedding_matrix(exact = True)*g).list())
        else:
            A=self.get_embedding_matrix(prec = prec)*Matrix(self._R,4,1,g)
            return Matrix(self._R,2,2,A.list())

    def get_embedding(self,prec=None):
        r"""
        Returns a function which embeds quaternions into a matrix algebra.

        EXAMPLES:

        ::

            sage: X = BTQuotient(5,3)
            sage: f = X.get_embedding(prec = 4)
            sage: b = Matrix(ZZ,4,1,[1,2,3,4])
            sage: f(b)
            [2 + 3*5 + 2*5^2 + 4*5^3 + O(5^4)       3 + 2*5^2 + 4*5^3 + O(5^4)]
            [        5 + 5^2 + 3*5^3 + O(5^4)           4 + 5 + 2*5^2 + O(5^4)]

        """
        A = self.get_embedding_matrix(prec = prec)
        return lambda g: Matrix(self._R,2,2,(A*g).list())

    def get_edge_stabs(self):
        r"""
        This function computes the stabilizers in the arithmetic group of all
        edges in the Bruhat-Tits tree within a fundamental domain for the 
        quotient graph.

        OUTPUT:

          A list of edge stabilizers. Each edge stabilizer is a finite cyclic subgroup,
          so we return generators for these subgroups.

        EXAMPLES:

        ::

            sage: X = BTQuotient(13,2)
            sage: S = X.get_edge_stabs()
            sage: gamma = X.embed_quaternion(S[0][0][0][0],prec = 20)
            sage: e = X.get_edge_list()[0].rep
            sage: X._BT.edge(gamma*e) == e
            True

        """
        try: return self._edge_stabs
        except AttributeError:
            self._edge_stabs=[self._stabilizer(e.rep,as_edge=True) for e in self.get_edge_list()]
            return self._edge_stabs

    def get_stabilizers(self):
        r"""
        This function computes the stabilizers in the arithmetic group of all
        edges in the Bruhat-Tits tree within a fundamental domain for the
        quotient graph, as well as those of opposite edges. Used for automorphic forms.

        OUTPUT:

          A list of edge stabilizers. Each edge stabilizer is a finite cyclic subgroup,
          so we return generators for these subgroups.

        EXAMPLES:

        ::

        """
        S = self.get_edge_stabs()
        return S + S

    def get_vertex_stabs(self):
        r"""
        This function computes the stabilizers in the arithmetic group of all
        vertices in the Bruhat-Tits tree within a fundamental domain for the 
        quotient graph.

        OUTPUT:

          A list of vertex stabilizers. Each vertex stabilizer is a finite cyclic 
          subgroup, so we return generators for these subgroups.

        EXAMPLES:

        ::

            sage: X = BTQuotient(13,2)
            sage: S = X.get_vertex_stabs()
            sage: gamma = X.embed_quaternion(S[0][0][0][0],prec = 20)
            sage: v = X.get_vertex_list()[0].rep
            sage: X._BT.vertex(gamma*v) == v
            True

        """
        try: return self._vertex_stabs
        except AttributeError:
            self._vertex_stabs=[self._stabilizer(e.rep,as_edge=False) for e in self.get_vertex_list()]
            return self._vertex_stabs

    def get_quaternion_algebra(self):
        r"""
        This function returns the underlying quaternion algebra.

        OUTPUT:

          The underlying definite quaternion algebra

        EXAMPLES:

        ::

            sage: X = BTQuotient(5,7)
            sage: X.get_quaternion_algebra()
            Quaternion Algebra (-1, -7) with base ring Rational Field

        """
        try: return self._A
        except AttributeError: pass
        self._init_order()
        return self._A

    def get_eichler_order(self):
        r"""
        This function returns the underlying Eichler order of level `N^+`.

        OUTPUT:

          Underlying Eichler order.

        EXAMPLES:

        ::

            sage: X = BTQuotient(5,7)
            sage: X.get_eichler_order()
            Order of Quaternion Algebra (-1, -7) with base ring Rational Field with basis (1/2 + 1/2*j, 1/2*i + 1/2*k, j, k)

        """
        try: return self._O
        except AttributeError: pass
        self._init_order()
        return self._O

    def get_splitting_field(self):
        r"""
        Returns an quadratic field splitting the quaternion algebra attached to ``self``. Currently requires Magma.

        EXAMPLES:

        ::

            sage: X = BTQuotient(5,11)
            sage: X.get_splitting_field()
            Traceback (most recent call last):
            ...
            NotImplementedError: Sage does not know yet how to work with the kind of orders that you are trying to use. Try installing Magma first and set it up so that Sage can use it.

        If we do have Magma installed, then it works:

        ::

            sage: X = BTQuotient(5,11,use_magma = True) # optional - magma
            sage: X.get_splitting_field() # optional - magma
            Number Field in a with defining polynomial X1^2 + 11

        """
        if self._use_magma == False:
            raise NotImplementedError,'Sage does not know yet how to work with the kind of orders that you are trying to use. Try installing Magma first and set it up so that Sage can use it.'
        try: return self._FF
        except AttributeError: pass
        self._compute_exact_splitting()
        return self._FF

    def get_eichler_order_basis(self):
        r"""
        This function returns an basis for the global Eichler order.

        OUTPUT:

          Basis for the underlying Eichler order of level Nplus.

        EXAMPLES:

        ::

            sage: X = BTQuotient(7,11)
            sage: X.get_eichler_order_basis()
            [  1/2 + 1/2*j]
            [1/2*i + 1/2*k]
            [            j]
            [            k]

        """
        try: return self._B
        except AttributeError: pass
        self._init_order()
        return self._B

    def get_eichler_order_quadform(self):
        r"""
        This function returns the norm form for the underlying Eichler order
        of level Nplus. Required for finding elements in the arithmetic 
        subgroup Gamma.

        OUTPUT:

          The norm form of the underlying Eichler order

        EXAMPLES:

        ::

            sage: X = BTQuotient(7,11)
            sage: X.get_eichler_order_quadform()
            Quadratic form in 4 variables over Integer Ring with coefficients: 
            [ 3 0 11 0 ]
            [ * 3 0 11 ]
            [ * * 11 0 ]
            [ * * * 11 ]

        """
        try: return self._OQuadForm
        except AttributeError: pass
        self._init_order()
        return self._OQuadForm

    def get_eichler_order_quadmatrix(self):
        r"""
        This function returns the matrix of the quadratic form of the underlying
        Eichler order in the fixed basis.

        OUTPUT:

          A 4x4 integral matrix describing the norm form.

        EXAMPLES:

        ::

            sage: X = BTQuotient(7,11)
            sage: X.get_eichler_order_quadmatrix()
            [ 6  0 11  0]
            [ 0  6  0 11]
            [11  0 22  0]
            [ 0 11  0 22]

        """
        try: return self._OM
        except AttributeError: pass
        self._init_order()
        return self._OM

    @cached_method
    def get_units_of_order(self):
        r"""
        This function returns the units of the underlying Eichler `\ZZ`-order. This is a finite
        group since the order lives in a definite quaternion algebra over `\QQ`.

        OUTPUT:

          A list of elements of the global Eichler `\ZZ`-order of level `N^+`.

        EXAMPLES:

        ::

            sage: X = BTQuotient(7,11)
            sage: X.get_units_of_order()
            [[ 0 -2  0  1], [-2  0  1  0]]

        """
        OM=self.get_eichler_order_quadmatrix()
        v=pari('qfminim(%s,2,0)'%(OM._pari_()))
        n_units=Integer(v[0].python()/2)
        v=pari('qfminim(%s,0,%s)'%((OM._pari_()),n_units))
        O_units=[]
        for jj in range(n_units):
            vec=Matrix(ZZ,1,4,[v[2][ii,jj].python() for ii in range(4)])
            O_units.append(vec)
        return O_units

    def get_CM_points(self,disc,prec):
        r"""
        Find the CM points corresponding to a given discriminant.

        INPUT:

        - ``disc``: An integer. Not necessarily a fundamental discriminant.
        - ``prec``: Integer. The precision to which compute the points.

        EXAMPLES::

            sage: p = 7
            sage: lev=2
            sage: X = BTQuotient(7,2,use_magma = True) # optional - magma
            sage: X.get_splitting_field() # optional - magma
            Number Field in a with defining polynomial X1^2 + 1
            sage: V = X.get_CM_points(5^4 * -11, 20) # optional - magma
            sage: len(V) # optional - magma
            20

        """
        p = self._p

        try:
            CM_points = self._CM_points[disc]
        except KeyError:
            self.get_splitting_field()
            self._magma.eval('rf := ReducedForms(%s)'%disc)
            h = ZZ(magma.eval(' # rf'))
            reduced_forms = [tuple(ZZ(magma.eval('rf[%s][%s]'%(ii,jj))) for jj in range(1,4)) for ii in range(1,h+1)]
            D0 = fundamental_discriminant(disc)
            cond = ZZ((disc/D0).sqrt())
            K = QuadraticField(D0, names = "delta")

            CM_points = []
            for f in reduced_forms:
                self._magma.eval('K<delta> := QuadraticField(%s)'%disc)
                self._magma.eval('O<w> := sub < MaximalOrder(K) | %s>'%cond)
                self._magma.eval('tau,_ := Embed(O,R)')
                CM_points.append(Matrix(ZZ,4,1,[ZZ(self._magma.eval('tau[%s]'%ii)) for ii in range(1,5)]))
            self._CM_points[disc] = CM_points

        out_list = []
        for elt in CM_points:
            a,b,c,d = self.embed_quaternion(elt,prec = prec).list()
            trace = a+d
            norm = a*d-b*c
            Kp=Qq(p**2,prec=prec,names='g')
            g=Kp.gen()
            # Compute the fixed points of the matrix m acting on the Kp points of Hp.
            A=Kp(a-d)
            D2=Kp(trace**2-4*norm)
            if D2==0:
                D=D2
            else:
                # Compute the square root of D in a naive way
                for a,b in product(range(p),repeat=2):
                    y0=a+b*g
                    if((y0**2-D2).valuation()>0):
                        break
                y1=y0
                D=0
                while(D!=y1):
                    D=y1
                    y1=(D**2+D2)/(2*D)
            out_list.append((A+D)/(2*c))
        return out_list

    @cached_method
    def _get_Up_data(self):
        r"""
        Returns (computes if necessary) Up data.

        This is a vector of length p,
        and each entry consists of the corresponding data for the matrix [p,a,0,1]
        where a varies from 0 to p-1.
        The data is a tuple (acter,edge_images), with edge images being of type
        ``DoubleCosetReduction``.

        EXAMPLES:

        """
        E=self.get_edge_list()
        vec_a=self._BT.subdivide([1],1)
        return [[alpha.inverse(),[DoubleCosetReduction(self,e.rep*alpha) for e in E]+[DoubleCosetReduction(self,e.opposite.rep*alpha) for e in E]] for alpha in vec_a]

    @cached_method
    def _get_atkin_lehner_data(self,q):
        r"""

        """
        E=self.get_edge_list()
        self._increase_precision(20)

        B=self.get_eichler_order_basis()
        BB=self._BB
        v=[]
        nninc=-2
        while(len(v)==0):
            nninc+=2
            #print 'Searching for norm', q*self._p**nninc
            v=self._find_elements_in_order(q*self._p**nninc)
        beta1=Matrix(QQ,4,1,v[0])
        success=False
        while(not success):
            try:
                x=self.embed_quaternion(beta1)
                nn=ceil(x.determinant().valuation())
                T=[beta1,[DoubleCosetReduction(self,x.adjoint()*e.rep,extrapow=nn) for e in E]]
                success=True
            except PrecisionError:
                self._increase_precision(10)
        return T

    @cached_method
    def _get_hecke_data(self,l):
        r"""

        """
        #print 'Finding representatives...'
        def enumerate_words(v):
            L=len(v)
            n=[]
            while(True):
                add_new=True
                for jj in range(len(n)):
                    n[jj]=(n[jj]+1)%L
                    if(n[jj]!=0):
                        add_new=False
                        break
                if(add_new):
                    n.append(0)
                wd=prod([v[x] for x in n])
                yield wd

        E=self.get_edge_list()
        self._increase_precision(20)

        if((self.level()*self.Nplus())%l==0):
            Sset=[]
        else:
            Sset=[self._p]
        B=self.get_eichler_order_basis()
        BB=self._BB
        T=[]
        T0=[]
        v=[]
        nninc=-2
        while(len(v)==0):
            nninc+=2
            #print 'Searching for norm', l*self._p**nninc
            v=self._find_elements_in_order(l*self._p**nninc)
        alpha1=v[0]
        alpha0=self._conv(alpha1)
        alpha=Matrix(QQ,4,1,alpha1)
        alphamat=self.embed_quaternion(alpha)
        A=self.get_quaternion_algebra()
        letters=self.get_generators()+[y[0] for Se in self.get_vertex_stabs() for y in Se if y[2]]
        I=enumerate_words([self._conv(x) for x in letters])
        n_tests=0
        while(len(T)<l+1):
            n_tests+=1
            v=I.next()
            v0=v*alpha0
            vinv=A((v0)**(-1))
            new=True
            for tt in T0:
                r=vinv*tt
                x=BB*Matrix(QQ,4,1,r.coefficient_tuple())
                if(all([x[jj,0].is_S_integral(Sset) for jj in range(4)])):
                    new=False
                    break
            if new:
                v1=BB*Matrix(QQ,4,1,v.coefficient_tuple())
                success = False
                while not success:
                    try:
                        x=self.embed_quaternion(v1)*alphamat
                        nn=ceil(x.determinant().valuation())
                        T.append([v1,[DoubleCosetReduction(self,x.adjoint()*e.rep,extrapow=nn) for e in E]])
                        success = True
                    except PrecisionError:
                        self._increase_precision(10)
                        alphamat=self.embed_quaternion(alpha)
                T0.append(v0)
        #print 'Done (used %s reps in total)'%(n_tests)
        return T,alpha

    def _find_equivalent_vertex(self,v0,V=None,valuation=None):
        r"""
        Finds a vertex in ``V`` equivalent to ``v0``.

        INPUT:
        - ``v0`` -- a 2x2 matrix in `\ZZ_p` representing a vertex in the Bruhat-Tits tree.

        - ``V`` -- list (Default: None) If a list of Vertex is given, restrict the search
            to the vertices in ``V``. Otherwise use all the vertices in a fundamental domain.

        - ``valuation`` -- an integer (Default: None): The valuation of the determinant
            of ``v0``, if known (otherwise it is calculated).

        OUTPUT:

        A pair ``g``, ``v``, where ``v`` is a Vertex in ``V`` equivalent to ``v0``, and
         ``g`` is such that `g\cdot v_0= v`.

        EXAMPLES::

        """
        try:
            return self._cached_vertices[v0]
        except KeyError: pass
        if V is None:
            V=self._vertex_list
        if valuation is None:
            valuation=v0.determinant().valuation(self._p)
        parity=valuation%2
        for v in filter(lambda v:v.parity==parity,V):
            g=self._are_equivalent(v0,v.rep,False,valuation+v.valuation)
            if g is not None:
                self._cached_vertices[v0]=(g,v)
                return g,v
        return 0,None

    def _find_equivalent_edge(self,e0,E=None,valuation=None):
        r"""
        Finds an edge in ``E`` equivalent to ``e0``.

        INPUT:
        - ``e0`` -- a 2x2 matrix in `\ZZ_p` representing an edge in the Bruhat-Tits tree.

        - ``E`` -- list (Default: None) If a list of Edge is given, restrict the search
            to the vertices in ``E``. Otherwise use all the edges in a fundamental domain.

        - ``valuation`` -- an integer (Default: None): The valuation of the determinant
            of ``e0``, if known (otherwise it is calculated).

        OUTPUT:

        A pair ``g``, ``e``, where ``e`` is an Edge in ``E`` equivalent to ``e0``, and
         ``g`` is such that `g\cdot e_0= e`.

        EXAMPLES::

        """
        try:
            return self._cached_edges[e0]
        except KeyError: pass
        if valuation is None:
            valuation=e0.determinant().valuation(self._p)
        parity=valuation%2
        if E is None:
            if parity == 0:
                E=self._edge_list
            else:
                E=[e.opposite for e in self._edge_list]
        for e in filter(lambda x:x.parity==parity,E):
            g=self._are_equivalent(e.rep,e0,True,valuation+e.valuation)
            if g is not None:
                self._cached_edges[e0]=(g,e)
                return g,e
        return 0,None

    def fundom_rep(self,v1):
        r"""
        Finds an equivalent vertex in the fundamental domain.

        INPUT:

        - ``v1`` - a 2x2 matrix representing a normalized vertex.

        OUTPUT:

        - A ``Vertex`` equivalent to ``v1``, in the fundamental domain.

        EXAMPLES:

        """
        try:
            tmp=self._cached_paths[v1]
            return tmp
        except KeyError: pass
        # print 'v1=',v1
        chain,v=self._BT.find_path(v1,self.get_vertex_dict())
        # print 'chain =', chain
        while(len(chain)>0):
            v0=chain.pop()
            g,v=self._find_equivalent_vertex(v0,V=[e.target for e in v.leaving_edges])
            assert not v is None
            self._cached_paths[v0]=v
        return v


    def _find_lattice(self,v1,v2,as_edges,m):
        r"""
        Find the lattice attached to the pair ``v1``,``v2``.

        INPUT:

        - ``v1``, ``v2`` - 2x2 matrices. They represent either a pair of normalized vertices or a pair of normalized edges.
        - ``as_edges`` - boolean. If True, the inputs will be considered as edges instead of vertices.
        - ``m`` - integer - The valuation of the determinant of ``v1``*``v2``.

        EXAMPLES:

        """
        if(as_edges):
            X=self._Xe
        else:
            X=self._Xv
        p=self._p
        if m+1>self._prec:
            self.get_embedding_matrix(prec = m+1)
        v1adj=v1.adjoint()
        R=self._Mat_44
        vecM=[v2*X[ii]*v1adj for ii in range(4)]
        M=(self._Iotainv*R([[vecM[ii][jj,kk] for ii in range(4) ] for jj in range(2) for kk in range(2)])).augment(R(self._pN)).transpose()
        E = M.echelon_form().submatrix(0,0,4,4)
        Et = E.transpose()
        return Et,E*self.get_eichler_order_quadmatrix()*Et

    def _stabilizer(self,e,as_edge=True):
        r"""
        Finds the stabilizer of an edge or vertex.

        EXAMPLES:

        """
        p=self._p
        m=e.determinant().valuation(p)
        twom=2*m
        E,A = self._find_lattice(e,e,as_edge,twom)
        n_units=len(self.get_units_of_order())
        ## Using PARI to get the shortest vector in the lattice (via LLL)
        v=pari('qfminim(%s,0,%s,flag = 0)'%(A._pari_(),2*n_units))

        mat=v[2].python().transpose()
        n_vecs=mat.nrows()
        stabs=[]
        for jj in range(n_vecs):
            vect = mat.row(jj).row()
            vec = vect.transpose()
            try: nrd=Integer((vec*A*vect)[0,0]/2)
            except TypeError: continue
            if nrd == p**twom:
                w=E*vect
                x=self._conv(w.transpose())
                w.set_immutable()
                stabs.append([w,m,x!=p**m])
        if len(stabs) <= 1:
            return [[self.B_one(),0,False]]
        else:
            return stabs

    def _are_equivalent(self,v1,v2,as_edges=False,twom=None,check_parity = False):
        r"""
        This function determines whether two vertices (or edges) of the Bruhat-Tits
        tree are equivalent under the arithmetic group in question. The computation
        boils down to an application of the LLL short-vector algorithm to a 
        particular lattice; for details see [FM].
        
        INPUT:

          - ``v1``, ``v2`` - two 2x2 integral matrices representing either vertices or edges
          
          - ``as_edges`` - boolean (Default: False). Tells whether the matrices should be
            interpreted as edges (if true), or as vertices (if false)
          
            - ``twom`` - integer (Default: None) If specified, indicates the valuation of the determinant of ``v1`` `\times` ``v2``.

        OUTPUT:

          If the objects are equivalent, this returns an element of the arithemtic
          group Gamma that takes v1 to v2. Otherwise is returns false.

        EXAMPLES:

        This example illustrates ...

        ::

        REFERENCES:

          [FM] "Computing quotients of the Bruhat-Tits tree...", Cameron Franc, Marc Masdeu.

        """
        try:
            return self._cached_equivalent[(v1,v2,as_edges)]
        except KeyError: pass
        p=self._p
        if twom is None:
            twom=v1.determinant().valuation(p)+v2.determinant().valuation(p)
        if check_parity:
            if twom % 2 != 0:
                self._cached_equivalent[(v1,v2,as_edges)]=None
                return None

        E,A=self._find_lattice(v1,v2,as_edges,twom)
        ## Using PARI to get the shortest vector in the lattice (via LLL)
        vec=pari('qfminim(%s,0,1,flag=0)'%(A._pari_()))[2].python()

        vect=vec.transpose()
        nrd=Integer((vect*A*vec)[0,0]/2)
        # print nrd
        # print p
        # print twom
        # print '--'
        if nrd == p**twom:
            m=Integer(twom/2)
            g=E*vec
            g.set_immutable()
            self._cached_equivalent[(v1,v2,as_edges)]=(g,m)
            return (g,m)
        else:
            self._cached_equivalent[(v1,v2,as_edges)]=None
            return None

    def _compute_exact_splitting(self):
        self._init_order()
        self._magma.eval('f:=MatrixRepresentation(R)')
        f=self._magma.function_call('MatrixRepresentation',args=[self._O],nvals=1)
        self._FF=NumberField(f.Codomain().BaseRing().DefiningPolynomial().sage(),'a')
        allmats=[]
        for kk in range(4):
           self._magma.eval('x:=f(R.%s)'%(kk+1))
           all_str=[]
           for ii in range(2):
               for jj in range(2):
                   self._magma.eval('v%s%s:=[Sage(z) : z in Eltseq(x[%s,%s])]'%(ii,jj,ii+1,jj+1))
           v=[self._FF(self._magma.eval('Sprint(v%s)'%tt)) for tt in ['00','01','10','11']]
           allmats.append(Matrix(self._FF,2,2,v))
        #print [(x.determinant(),x.trace()) for x in allmats]
        self._Iota_exact=Matrix(self._FF,4,4,[self._FF(allmats[kk][ii,jj]) for ii in range(2) for jj in range(2) for kk in range(4) ])

    def _init_order(self):
        r"""
        Initialize the order of the quaternion algebra. Here we
        possibly use Magma to split it.

        EXAMPLES:

            sage: X = BTQuotient(2,3,5,use_magma = False)
        """
        if(self._use_magma == True):
            A=self._magma.QuaternionAlgebra(self._Nminus)
            self._magma.eval('A:=QuaternionAlgebra(%s)'%(self._Nminus))
            self._magma.eval('R:=QuaternionOrder(A,%s)'%(self._Nplus))
            g=A.gens()
            # We store the order because we need to split it
            self._O=A.QuaternionOrder(self._Nplus)
            OBasis=self._O.Basis()
            self._A=QuaternionAlgebra((g[0]**2).sage(),(g[1]**2).sage())
            v=[1]+self._A.gens()
            self._B=Matrix(self._A,4,1,[sum([OBasis[tt+1][rr+1].sage()*v[rr] for rr in range(4)]) for tt in range(4)])

        else:
            # Note that we can't work with non-maximal orders in sage
            assert(self._Nplus==1)
            self._A=QuaternionAlgebra(self._Nminus)
            v=[1]+self._A.gens()
            self._O=self._A.maximal_order()
            OBasis=self._O.basis()

            self._B=Matrix(self._A,4,1,[OBasis[tt] for tt in range(4)])
        self._OQuadForm=QuadraticForm(self._Mat_44([(self._B[ii,0]*self._B[jj,0].conjugate()).reduced_trace() for ii in range(4) for jj in range(4)]))
        self._OM=self._OQuadForm.matrix()
        self._BB=Matrix(QQ,4,4,[[self._B[ii,0][jj] for ii in range(4)] for jj in range(4)]).inverse()

    def B_one(self):
        r"""
        Returns the coordinates of `1` in the basis for the quaternion order.

        EXAMPLES:

            sage: X = BTQuotient(7,11)
            sage: v,pow = X.B_one()
            sage: X._conv(v) == 1
            True

        """
        try: return self._B_one
        except AttributeError:
            V = self.get_units_of_order()
            for v in V:
                vt = v.transpose()
                vt.set_immutable()
                b = self._conv(v)
                if b == 1:
                    self._B_one = (vt,0)
                    break
                if b == -1:
                    self._B_one = (-vt,0)
                    break
            return self._B_one


    def _conv(self,v):
        r"""
        Returns a quaternion having coordinates in the fixed basis for the
        order given by ``v``.

        EXAMPLES:

            sage: X = BTQuotient(5,7)
            sage: A = X.get_quaternion_algebra()
            sage: i,j,k = A.gens()
            sage: B = X.get_eichler_order_basis()
            sage: X._conv([1,2,3,4]) == B[0,0]+2*B[1,0]+3*B[2,0]+4*B[3,0]
            True

        """
        if hasattr(v,"list"):
            v=v.list()
        return (Matrix(QQ,1,4,v)*self.get_eichler_order_basis())[0,0]

    @cached_method
    def _find_elements_in_order(self, norm, trace = None, primitive=False):
        r"""
        Returns elements in the order of the quaternion algebra of specified
        reduced norm, and possibly reduced trace.

        INPUT:

        - ``norm`` - integer. The required reduced norm.
        - ``trace`` - integer (Default: None). If specified, returns elements only
        reduced trace ``trace``.
        - ``primitive`` boolean (Default: False). If True, return only elements that
        cannot be divided by `p`.

        EXAMPLES:

            sage: X = BTQuotient(5,7)
            sage: X._find_elements_in_order(23)
            [[2, 9, -1, -5], [0, 8, 0, -5], [-2, 9, 1, -5], [6, 7, -3, -4], [2, 5, -1, -4], [0, 6, -1, -4], [0, 8, -1, -4], [2, 9, -1, -4], [-2, 5, 1, -4], [0, 6, 1, -4], [0, 8, 1, -4], [-2, 9, 1, -4], [-6, 7, 3, -4], [7, 6, -4, -3], [7, 6, -3, -3], [6, 7, -3, -3], [0, 8, 0, -3], [-7, 6, 3, -3], [-6, 7, 3, -3], [-7, 6, 4, -3], [0, 1, -1, -2], [0, 6, -1, -2], [0, 1, 1, -2], [0, 6, 1, -2], [9, 2, -5, -1], [6, 0, -4, -1], [8, 0, -4, -1], [5, 2, -4, -1], [9, 2, -4, -1], [1, 0, -2, -1], [6, 0, -2, -1], [0, -1, -1, -1], [-1, 0, -1, -1], [5, 2, -1, -1], [2, 5, -1, -1], [0, -1, 1, -1], [1, 0, 1, -1], [-5, 2, 1, -1], [-2, 5, 1, -1], [-6, 0, 2, -1], [-1, 0, 2, -1], [-8, 0, 4, -1], [-6, 0, 4, -1], [-9, 2, 4, -1], [-5, 2, 4, -1], [-9, 2, 5, -1], [8, 0, -5, 0], [8, 0, -3, 0]]
            sage: X._find_elements_in_order(23,1)
            [[1, 0, -2, -1], [1, 0, 1, -1]]

        """
        OQuadForm=self.get_eichler_order_quadform()
        V=OQuadForm.vectors_by_length(norm)[norm]
        W=V if not primitive else filter(lambda v: any((vi%self._p!=0 for vi in v)),V)
        return W if trace is None else filter(lambda v:self._conv(v).reduced_trace() == trace,W)

    def _compute_quotient(self):
        r"""
        Computes the quotient graph.

        EXAMPLES:

            sage: X = BTQuotient(11,2)
            sage: X.get_graph()
            Multi-graph on 2 vertices

            sage: X = BTQuotient(17,19)
            sage: X.get_graph()
            Multi-graph on 4 vertices

            sage: X = BTQuotient(5,7,12)
            sage: X.get_graph()
            Multi-graph on 24 vertices
            sage: len(X._edge_list)
            72

            sage: X = BTQuotient(2,3,5)
            sage: X.get_graph()
            Multi-graph on 4 vertices

            sage: X = BTQuotient(2,3,35)
            sage: X.get_graph()
            Multi-graph on 16 vertices

            sage: X = BTQuotient(53,11,2)
            sage: X.get_graph()
            Multi-graph on 6 vertices

            sage: X = BTQuotient(2,13,9)
            sage: X.get_graph()
            Multi-graph on 24 vertices

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu

        """
        generators=set([])
        genus=self.genus()
        num_edges=0
        self.get_embedding_matrix(prec = 1)
        p=self._p
        num_verts=0
        v0=Vertex(self,num_verts,self._Mat_22([1,0,0,1]),determinant = 1,valuation = 0)
        V=collections.deque([v0])
        S=Graph(0,multiedges=True,weighted=True)
        Sfun = Graph(0)
        edge_list=[]
        vertex_list=[v0]
        self._num_edges=0
        num_verts+=1
        n_units=len(self.get_units_of_order())
        while len(V)>0:
            v=V.popleft()
            # found_edges=sum([len(e.links) for e in v.leaving_edges])
            # v_stabilizer=len(self._stabilizer(v.rep,as_edge=False))
            # total_edges=(p+1)/v_stabilizer

            E=self._BT.leaving_edges(v.rep)

            # print 'V = %s, E = %s, G = %s (target = %s), lenV = %s'%(num_verts,num_edges,1+num_edges-num_verts,genus,len(V))
            for e in E: #ii in range(p+1):
                # if found_edges == total_edges:
                #     break

                #e=E[ii]
                edge_det=e.determinant()
                edge_valuation=edge_det.valuation(p)
                # if v_stabilizer == 1 and ii == len(v.leaving_edges):
                #     e1 = None
                # else:
                #     g,e1=self._find_equivalent_edge(e,v.leaving_edges,valuation=edge_valuation)

                g,e1=self._find_equivalent_edge(e,v.leaving_edges,valuation=edge_valuation)

                if e1 is not None: # The edge is old. We just update the links
                    e1.links.append(g)
                    target = self._BT.target(e)
                    if e1.parity == 0:
                        Sfun.add_edge(v.rep,target,label = e1.label)
                    else:
                        Sfun.add_edge(v.rep,target,label = e1.opposite.label)

                    Sfun.set_vertex(target,e1.target)
 
                else: # The edge is new.
                    target=self._BT.target(e)
                    target.set_immutable()
                    new_det=target.determinant()
                    new_valuation=new_det.valuation(p)
                    new_parity=new_valuation%2
                    g1,v1=self._find_equivalent_vertex(target,V,valuation=new_valuation)
                    if v1 is None:
                        #The vertex is also new
                        v1=Vertex(self,num_verts,target,determinant = new_det,valuation = new_valuation)
                        vertex_list.append(v1)
                        num_verts+=1
                        #Add the vertex to the list of pending vertices
                        V.append(v1)
                    else:
                        generators.add(g1[0])


                    # Add the edge to the list
                    new_e=Edge(self,num_edges,e,v,v1,determinant = edge_det,valuation = edge_valuation)
                    new_e.links.append(self.B_one())
                    Sfun.add_edge(v.rep,target,label = num_edges)
                    Sfun.set_vertex(target,v1)

                    # Add the edge to the graph
                    S.add_edge(v.rep,v1.rep,num_edges)
                    S.set_vertex(v.rep,v)
                    S.set_vertex(v1.rep,v1)

                    # Find the opposite edge
                    opp=self._BT.opposite(e)
                    opp_det=opp.determinant()
                    new_e_opp=Edge(self,num_edges,opp,v1,v,opposite = new_e)
                    new_e.opposite=new_e_opp

                    if new_e.parity == 0:
                        edge_list.append(new_e)
                    else:
                        edge_list.append(new_e_opp)

                    v.leaving_edges.append(new_e)
                    v.entering_edges.append(new_e_opp)
                    v1.entering_edges.append(new_e)
                    v1.leaving_edges.append(new_e_opp)
                    num_edges += 1
                    # found_edges+=ZZ(v_stabilizer/len(self._stabilizer(e,as_edge=True)))
            # if genus == 1 - len(vertex_list) +num_edges:
            #     break
        computed_genus=Integer(1- len(vertex_list)+num_edges)
        if computed_genus != genus:
            print 'You found a bug! Please report!'
            print 'Computed genus =',computed_genus
            print 'Theoretical genus =', genus
            raise RuntimeError

        self._generators=generators
        self._boundary= {v.rep:v for v in vertex_list}
        self._edge_list=edge_list
        self._vertex_list=vertex_list
        self._S=S
        self._Sfun=Sfun


