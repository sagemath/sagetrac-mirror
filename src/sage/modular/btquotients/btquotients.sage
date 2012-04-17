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
        sage: d.igamma()*Y._edge_list[d.label].opposite.rep*d.t()==x
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
        self.label=label
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
        if(self.parity==0):
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
        e=Y._edge_list[self.label]
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
    def genus_old(self):
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

        REFERENCE:

        Source: Rotger, V. "Non-elliptic Shimura curves of genus one", Proposition 5.2

        """
        Nplus=self._Nplus
        lev=self._p*self._Nminus
        # Compute the genus using the genus formulas
        e4=1
        e3=1
        mu=Nplus
        for f in lev.factor():
            e4*=(1-kronecker_symbol(-4,Integer(f[0])))
            e3*=(1-kronecker_symbol(-3,Integer(f[0])))
            mu*=Integer(f[0])-1
        for f in Nplus.factor():
            if (f[1]==1):
                e4*=(1+kronecker_symbol(-4,Integer(f[0])))
                e3*=(1+kronecker_symbol(-3,Integer(f[0])))
            else:
                if(kronecker_symbol(-4,Integer(f[0]))==1):
                    e4*=2
                else:
                    e4=0
                if(kronecker_symbol(-3,Integer(f[0]))==1):
                    e3*=2
                else:
                    e3=0
            mu*=1+1/Integer(f[0])
        if(lev==1):
            einf=sum([euler_phi(gcd(d,Integer(Nplus/d))) for d in Integer(Nplus).divisors()])
        else:
            einf=0
        return 1+Integer(mu/12-e3/3-e4/4-einf/2)

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
            [1 1 1 2]
            [2 0 1 1]
            [2 1 1 1]
            [0 2 2 1]
            sage: R = A.base_ring()
            sage: B = X.get_eichler_order_basis()
            sage: R(B[0,0].reduced_trace()) == A[0,0]+A[3,0]
            True

        """
        if self._use_magma==True:
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
            return Matrix(Zmod(self._pN),4,4,[phi(B[kk,0])[ii,jj] for ii in range(2) for jj in range(2) for kk in range(4)])

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

        This example illustrates ...

        ::

        """
        return lambda g:Matrix(self._R,2,2,(self.get_embedding_matrix(prec = prec)*g).list())

    def get_edge_stabs(self):
        r"""
        This function computes the stabilizers in the arithmetic group of all
        edges in the Bruhat-Tits tree within a fundamental domain for the 
        quotient graph.

        OUTPUT:

          A list of edge stabilizers. Each edge stabilizer is a finite cyclic subgroup,
          so we return generators for these subgroups.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        try: return self._edge_stabs
        except AttributeError:
            self._edge_stabs=[self._stabilizer(e.rep,as_edge=True) for e in self.get_edge_list()]
            return self._edge_stabs

    def get_vertex_stabs(self):
        r"""
        This function computes the stabilizers in the arithmetic group of all
        vertices in the Bruhat-Tits tree within a fundamental domain for the 
        quotient graph.

        OUTPUT:

          A list of vertex stabilizers. Each vertex stabilizer is a finite cyclic 
          subgroup, so we return generators for these subgroups.

        EXAMPLES:

        This example illustrates ...

        ::

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
        """
        try: return self._O
        except AttributeError: pass
        self._init_order()
        return self._O

    def get_splitting_field(self):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
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

    @cached_method
    def CM_point(self,f):
        r"""
        Find the CM point corresponding to a given binary quadratic form.

        INPUT:

            - ``f``: A tuple specifying a quadratic form.

        EXAMPLES:

        """
        p=self._p
        A,B,C = f
        D = B**2 - 4*A*C
        K.<delta> = QuadraticField(D)
        D0 = K.discriminant()
        cond = ZZ((D/D0).sqrt())

        # w = (-B + delta)/2
        self._magma.eval('K<delta> := QuadraticField(%s)'%D)
        self._magma.eval('O<w> := sub < MaximalOrder(K) | %s>'%cond)
        self._magma.eval('tau,_ := Embed(O,R)')
        elt = Matrix(ZZ,4,1,[ZZ(self._magma.eval('tau[%s]'%ii)) for ii in range(1,5)])
        return elt

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
        self.get_splitting_field()
        self._magma.eval('rf := ReducedForms(%s)'%disc)
        h = ZZ(magma.eval(' # rf'))
        reduced_forms = [(ZZ(magma.eval('rf[%s][%s]'%(ii,jj))) for jj in range(1,4)) for ii in range(1,h+1)]

        out_list = []
        for f in reduced_forms:
            elt = self.CM_point(f)
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

    # @cached_method
    # def get_CM_points(self,val,prec,ignore_units=False):
    #     r"""
    #     Find the CM points corresponding to a given order.
    # 
    #     INPUT:
    #         - ``val``: An integer or an order. If it is an integer, consider the maximal order `O` in the field `\QQ(\sqrt{val})`.
    # 
    #         - ``ignore_units``: boolean (Default: False). If True, will return points which differ by a unit.
    #         - ``prec``: Integer. The precision to which compute the points.
    # 
    #     EXAMPLES:
    # 
    #     """
    #     p=self._p
    #     try:
    #         disc=Integer(val)
    #         R = QQ['x']
    #         f = R([-disc, 0, 1])
    #         K=NumberField(f, 'sq', check=False)
    #         val=K.maximal_order()
    #     except TypeError:
    #         disc=val.discriminant()
    #         # K=val.fraction_field()
    #     w=val.ring_generators()[0]
    #     trace = w.trace()
    #     norm = w.norm()
    #     try: all_elts=self._CM_points[val]
    #     except KeyError:
    #         if(not self.is_admissible(disc)):
    #             return []
    #         all_elts=[]
    #         nn=-1
    #         while(len(all_elts)==0):
    #             nn+=1
    #             #print 'nn=',nn
    #             all_elts=self._find_elements_in_order(p**(2*nn)*norm,(p**nn)*trace)
    # 
    #         all_elts=[[v[ii]/(p**nn) for ii in range(4)] for v in all_elts]
    #         # Now we take into account the action of units
    #         all_elts0=[self._conv(v) for v in all_elts]
    # 
    #         if not ignore_units:
    #             units=[[u[ii]/(p**nn) for ii in range(4)] for u in self._find_elements_in_order(p**(2*nn))]
    #             units0=[self._conv(u) for u in units]
    #         else:
    #             units=[]
    #             units0=[]
    # 
    #         all_elts_purged0=[]
    #         all_elts_purged=[]
    # 
    #         while(len(all_elts0)>0):
    #             v0=all_elts0.pop(0)
    #             v1=all_elts.pop(0)
    #             new=True
    #             for tt in all_elts_purged0:
    #                 #compare v1 with tt
    #                 for u in units0:
    #                     if(tt*u==u*v0):
    #                         new=False
    #                         break
    #                 if(not new):
    #                     break
    #             if(new):
    #                 all_elts_purged0.append(v0)
    #                 all_elts_purged.append(v1)
    # 
    #         self._CM_points[val]=copy(all_elts_purged)
    # 
    #     V=self._CM_points[val]
    #     all_elts_split=[self.embed_quaternion(matrix(4,1,y),prec = prec) for y in V]
    #     Kp=Qq(p**2,prec=prec,names='g')
    #     g=Kp.gen()
    #     W=[]
    #     for m in all_elts_split:
    #         m.set_immutable()
    #         # Compute the fixed points of the matrix m acting on the Kp points of Hp.
    #         a,b,c,d = m.list()
    #         A=Kp(a-d)
    #         D2=Kp(trace**2-4*norm)
    #         if D2==0:
    #             D=D2
    #         else:
    #             # Compute the square root of D in a naive way
    #             for a,b in product(range(p),repeat=2):
    #                 y0=a+b*g
    #                 if((y0**2-D2).valuation()>0):
    #                     break
    #             y1=y0
    #             D=0
    #             while(D!=y1):
    #                 D=y1
    #                 y1=(D**2+D2)/(2*D)
    # 
    #         W.append( (A+D)/(2*c) )
    #     return W

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


#########################################################################
#       Copyright (C) 2011 Cameron Franc and Marc Masdeu
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################
from collections import namedtuple
from sage.structure.element import Element
from sage.structure.element import ModuleElement
from sage.modules.module import Module
from sage.rings.all import Integer
from sage.structure.element import Element
from sage.matrix.constructor import Matrix, zero_matrix
from sage.matrix.matrix_space import MatrixSpace_generic
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.all import Qp
from sage.rings.all import RationalField
from sage.rings.number_field.all import NumberField
from copy import copy
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.modular.hecke.all import (AmbientHeckeModule, HeckeSubmodule, HeckeModuleElement)
import sage.rings.arith as arith
import sage.modular.hecke.hecke_operator


class HarmonicCocycleElement(HeckeModuleElement):
    r"""
    Objects of this type are Gamma-invariant harmonic cocycles on the 
    Bruhat-Tits tree. Gamma-invariance is necessary so that the cocycle
    can be stored in terms of a finite amount of data.

    More precisely, given a BTQuotient T, we store harmonic cocycles as 
    a list of values in some coefficient module (e.g. for weight 2 forms
    can take Cp) indexed by edges of a fundamental domain for T in the 
    Bruhat-Tits tree. Evaluate the cocycle at other edges using Gamma
    invariance (although the values may not be equal over an orbit of
    edges as the coefficient module action may be nontrivial).

    INPUT:

     - ``vec`` - (default: None) 
     - ``from_values`` -  (default: False) 

    EXAMPLES:

    ::

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu
    """
    def __init__(self,_parent,vec=None,from_values=False):
        HeckeModuleElement.__init__(self,_parent,None)
        self._parent=_parent
        if from_values:
            self._R=_parent._U._R
            self._wt=_parent._k
            self._nE=len(_parent._E)
            self._F=copy(vec)
            return
        if isinstance(vec, self.__class__):
            #The first argument is just a modular form that we copy as is
            self._R=vec._R
            self._nE=vec._nE
            self._wt=vec._wt
            self._F=[_parent._U.element_class(_parent._U,vec._F[ii],quick=True) for ii in range(self._nE)]
            return
        # When vec contains coordinates for the basis
        self._R=self._parent._R
        try:
            v=[self._R(x) for x in vec.list()]
        except AttributeError:
            v=[self._R(vec) for ii in range(self._parent.dimension())]
        self._wt=_parent._k
        self._nE=len(_parent._E)
        vmat=Matrix(self._R,1,self._parent.dimension(),v)
        tmp=(vmat*self._parent.ambient_module().basis_matrix()).row(0)
        self._F=[_parent._U.element_class(_parent._U,Matrix(self._R,self._wt-1,1,tmp[e*(self._wt-1):(e+1)*(self._wt-1)]),quick=True) for e in range(self._nE)]
        return

    def _add_(self,g):
        r"""
        This function adds two cocycles componentwise.
        """
        #Should ensure that self and g are modular forms of the same weight and on the same curve
        C=self.__class__
        return C(self._parent,self.element()+g.element())
    
    def _sub_(self,g):
        r"""
        This function computes the difference of two cocycles.
        """
        #Should ensure that self and g are modular forms of the same weight and on the same curve
        C=self.__class__
        return C(self._parent,self.element()-g.element())

    def _rmul_(self,a):
        r"""
        This function multiplies a cocycle by a scalar.
        """
        #Should ensure that 'a' is a scalar
        C=self.__class__
        return C(self._parent,a*self.element())


    def _repr_(self):
        r"""
        Retuns a string representing the values of the cocycle on the edges.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        tmp='Element of '+str(self._parent)
        return tmp

    def __eq__(self,other):
        r"""
        Test for equality with another cocycle. Two cocycles are equal if they
        take the same values on all the edges in a fundamental domain.

        EXAMPLES:

        This example illustrates ...

        ::
        """
        return all([self._F[e].__eq__(other._F[e]) for e in range(self._nE)])

    def __ne__(self,other):
        r"""
        Test for non-equality with another cocycle. Two cocycles are non-equal if they
        take the different values on at least one edge in a fundamental domain.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        return any([self._F[e].__ne__(other._F[e]) for e in range(self._nE)])

    def __nonzero__(self):
        r"""
        Test for being non-zero.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        return any([self._F[e].__nonzero__() for e in range(self._nE)])

    def valuation(self):
        r"""
        Returns the valuation of the cocycle, defined as the minimum of the values
        it takes on a set of representatives.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        if not self.__nonzero__():
            return oo
        else:
            return min([self._F[e].valuation() for e in range(self._nE)])

    def _compute_element(self):
        r"""
    
        """
        R=self._R
        A=self._parent.basis_matrix().transpose()
        B=Matrix(R,self._nE*(self._parent._k-1),1,[self._F[e]._val[ii,0]  for e in range(self._nE) for ii in range(self._parent._k-1) ])
        res=(A.solve_right(B)).transpose()
        return self._parent.free_module()(res.row(0))

    #In HarmonicCocycle
    def evaluate(self,e1):
        r"""
        This function evaluates the cocycle on an edge of the Bruhat-Tits tree.

        INPUT:

        - ``e1`` - a matrix corresponding to an edge of the Bruhat-Tits tree

        OUTPUT

          An element of the coefficient module of the cocycle which describes 
          the value of the cocycle on e1

        EXAMPLES:
        """
        X=self._parent._X
        p=X._p
        u=DoubleCosetReduction(X,e1)
        return (u.sign()*self._F[u.label]).l_act_by(u.igamma(self._parent.embed_quaternion)*(p**(-u.power)))

    #In HarmonicCocycle
    def riemann_sum(self,f,center=1,level=0,E=None):
        r"""
        This function evaluates the integral of the funtion ``f`` with respect to the measure determined by ``self``.

        EXAMPLES:

        This example illustrates ...

        ::
        """
        R1=LaurentSeriesRing(f.base_ring(),'r1')
        R1.set_default_prec(self._parent._k-1)
        R2=PolynomialRing(f.base_ring(),'r2')

        if E is None:
            E=self._parent._X._BT.get_balls(center,level)
        else:
            E=self._parent._X._BT.subdivide(E,level)
        value=0
        ii=0
        for e in E:
            ii+=1
            exp=R1([e[1,1],e[1,0]])**(self._parent._k-2)*e.determinant()**(-(self._parent._k-2)/2)*f(R1([e[0,1],e[0,0]])/R1([e[1,1],e[1,0]])).truncate(self._parent._k-1)
            new=self.evaluate(e).l_act_by(e.inverse()).evaluate(exp)
            value+=new
        return value

    def modular_form(self,z=None,level=0):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        return self.derivative(z,level,order=0)

    # In HarmonicCocycle
    def derivative(self,z=None,level=0,order=1):
        r"""
        The derivative of the modular form attached to ``self``.

        EXAMPLES:
        ::
        """
        def F(z):
            R=PolynomialRing(z.parent(),'x,y').fraction_field()
            Rx=PolynomialRing(z.parent(),'x1').fraction_field()
            x1=Rx.gen()
            subst=R.hom([x1,z],codomain=Rx)
            x,y=R.gens()
            center=self._parent._X._BT.find_containing_affinoid(z)
            zbar=z.trace()-z
            f=R(1)/(x-y)
            k=self._parent._k
            V=[f]
            for ii in range(order):
                V=[v.derivative(y) for v in V]+[k/(y-zbar)*v for v in V]
                k+=2
            return sum([self.riemann_sum(subst(v),center,level) for v in V])
        if(z is None):
            return F
        else:
            return F(z)


class HarmonicCocycles(AmbientHeckeModule):
    Element=HarmonicCocycleElement
    r"""
    This object represents a space of Gamma invariant harmonic cocycles valued in
    a cofficient module.

    INPUT:

     - ``X`` - A BTQuotient object

     - ``k`` - integer - The weight.

     - ``prec`` - integer (Default: None). If specified, the precision for the coefficient module

     - ``basis_matrix`` - integer (Default: None)

     - ``base_field`` - (Default: None)

    EXAMPLES:

    ::

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu
    """
    def __init__(self,X,k,prec=None,basis_matrix=None,base_field=None):
        self._k=k
        self._X=X
        self._E=self._X.get_edge_list()
        self._V=self._X.get_vertex_list()

        if prec is None:
            self._prec=None
            if base_field is None:
                try:
                    self._R= X.get_splitting_field()
                except AttributeError:
                    raise ValueError, "It looks like you are not using Magma as backend...and still we don't know how to compute splittings in that case!"
            else:
                pol=X.get_splitting_field().defining_polynomial().factor()[0][0]
                self._R=base_field.extension(pol,pol.variable_name()).absolute_field(name='r')
            self._U=OCVn(self._k-2,self._R)
        else:
            self._prec=prec
            if base_field is None:
                self._R=Qp(self._X._p,prec=prec)
            else:
                self._R=base_field
            self._U=OCVn(self._k-2,self._R,self._k-1)
        self.__rank = self._X.dimension_harmonic_cocycles(self._k)
        if basis_matrix is not None:
            self.__matrix=basis_matrix
            self.__matrix.set_immutable()
            assert self.__rank == self.__matrix.nrows()

        AmbientHeckeModule.__init__(self, self._R, self.__rank, self._X.prime()*self._X.Nplus()*self._X.Nminus(), weight=self._k)
        self._populate_coercion_lists_()

    def base_extend(self,base_ring):
        r"""
        This function extends the base ring of the coefficient module.
        
        INPUT:

        - ``base_ring`` - a ring that has a coerce map from the current base ring

        EXAMPLES:

        ::
        """
        if not base_ring.has_coerce_map_from(self.base_ring()):
            raise ValueError, "No coercion defined"
        else:
            return self.change_ring(base_ring)

    def change_ring(self, new_base_ring):
        r"""
        This function changes the base ring of the coefficient module.

        INPUT:

        - ``new_base_ring'' - a ring that has a coerce map from the current base ring

        EXAMPLES:
        ::

        """
        if not new_base_ring.has_coerce_map_from(self.base_ring()):
            raise ValueError, "No coercion defined"
        else:
            return self.__class__(self._X,self._k,prec=self._prec,basis_matrix=self.basis_matrix().change_ring(base_ring),base_field=new_base_ring)

    def rank(self):
        r"""
        The rank (dimension) of ``self``.
        
        EXAMPLES:
        ::

        """
        return self.__rank

    def submodule(self,v,check=False):
        r"""
        Return the submodule of ``self`` spanned by ``v``.

        EXAMPLES:
        """
        return HarmonicCocyclesSubmodule(self,v,dual=None,check=check)

    def is_simple(self):
        r"""
        Whether ``self`` is irreducible.

        EXAMPLES:
        ::

        """
        return self.rank()==1

    def _repr_(self):
        r"""
        This returns the representation of self as a string.
        """
        return 'Space of harmonic cocycles of weight %s on %s'%(self._k,self._X)

    def _latex_(self):
        r"""
        A LaTeX representation of ``self``.

        EXAMPLES:
        ::
        """
        s='\\text{Space of harmonic cocycles of weight }'+latex(self._k)+'\\text{ on }'+latex(self._X)
        return s

    def _an_element_(self):
        r"""

        """
        return self.basis()[0]


    def _coerce_map_from_(self, S):
        r"""
        Can coerce from other HarmonicCocycles or from pAutomorphicForms
        """
        if isinstance(S,(HarmonicCocycles,pAutomorphicForms)):
            if(S._k!=self._k):
                return False
            if(S._X!=self._X):
                return False
            return True
        return False

    def __cmp__(self,other):
        r"""

        """
        try:
            res=(self.base_ring()==other.base_ring() and self._X==other._X and self._k==other._k)
            return res
        except:
            return False

    def _element_constructor_(self,x):
        r"""

        """
        #Code how to coherce x into the space
        #Admissible values of x?
        if isinstance(x,HarmonicCocycleElement):
            return HarmonicCocycleElement(self,x)
        elif isinstance(x,pAutomorphicForm):
            tmp=[self._U.element_class(_parent._U,x._F[ii]).l_act_by(self._E[ii].rep) for ii in range(self._nE)]
            return HarmonicCocycleElement(self,tmp,from_values=True)
        else:
            return HarmonicCocycleElement(self,x)


    def free_module(self):
        r"""
        This function returns the underlying free module

        EXAMPLES:
        ::
        """
        try: return self.__free_module
        except AttributeError: pass
        V = self.base_ring()**self.dimension()
        self.__free_module = V
        return V

    def character(self):
        r"""
        Only implemented the trivial character so far.

        EXAMPLES:

        """
        return lambda x:x

    def embed_quaternion(self,g):
        r"""
        Embed the quaternion element ``g`` into the matrix algebra.

        EXAMPLES:
        ::
        """
        return self._X.embed_quaternion(g,exact = self._R.is_exact(), prec = self._prec)

    def basis_matrix(self):
        r"""
        Returns a basis of ``self`` in matrix form.

        If the coefficient module M is of finite rank then the space of Gamma invariant
        M valued harmonic cocycles can be represented as a subspace of the finite rank
        space of all functions from the finitely many edges in the corresponding 
        BTQuotient into M. This function computes this representation of the space of
        cocycles.

        OUTPUT:

          A basis matrix describing the cocycles in the spaced of all M valued Gamma
          invariant functions on the tree

        EXAMPLES:

        This example illustrates ...

        ::

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu (2012-02-20)
        """
        try: return self.__matrix
        except AttributeError: pass
        nV=len(self._V)
        nE=len(self._E)
        stab_conds=[]
        S=self._X.get_edge_stabs()
        p=self._X._p
        d=self._k-1
        for e in self._E:
            try:
                g=filter(lambda g:g[2],S[e.label])[0]
                C=self._U.l_matrix_representation(self.embed_quaternion(g[0]))
                C-=self._U.l_matrix_representation(Matrix(QQ,2,2,p**g[1]))
                stab_conds.append([e.label,C])
            except IndexError: pass

        n_stab_conds=len(stab_conds)
        self._M=Matrix(self._R,(nV+n_stab_conds)*d,nE*d,0,sparse=True)
        for v in self._V:
            for e in filter(lambda e:e.parity==0,v.leaving_edges):
                C=sum([self._U.l_matrix_representation(self.embed_quaternion(x[0])) for x in e.links],Matrix(self._R,d,d,0))
                self._M.set_block(v.label*d,e.label*d,C)
            for e in filter(lambda e:e.parity==0,v.entering_edges):
                C=sum([self._U.l_matrix_representation(self.embed_quaternion(x[0])) for x in e.opposite.links],Matrix(self._R,d,d,0))
                self._M.set_block(v.label*d,e.opposite.label*d,C)

        for kk in range(n_stab_conds):
            v=stab_conds[kk]
            self._M.set_block((nV+kk)*d,v[0]*d,v[1])

        x1=self._M.right_kernel().matrix()
        K=[c for c in x1.rows()]

        if not self._R.is_exact():
            for ii in range(len(K)):
                s=min([t.valuation() for t in K[ii]])
                for jj in range(len(K[ii])):
                    K[ii][jj]=(p**(-s))*K[ii][jj]

        self.__matrix=Matrix(self._R,len(K),nE*d,K)
        self.__matrix.set_immutable()
        return self.__matrix

    def __apply_atkin_lehner(self,q,f):
        r"""
        This function applies an Atkin-Lehner involution to a harmonic cocycle

        INPUT:

          - ``q`` - an integer dividing the full level p*Nminus*Nplus

          - ``f`` - a harmonic cocycle

        OUTPUT:

          The harmonic cocycle obtained by hitting f with the Atkin-Lehner at q

        EXAMPLES:
        ::
        """
        R=self._R
        Data=self._X._get_atkin_lehner_data(q)
        p=self._X._p
        tmp=[self._U.element_class(self._U,zero_matrix(self._R,self._k-1,1),quick=True) for jj in range(len(self._E))]
        d1=Data[1]
        mga=self.embed_quaternion(Data[0])
        for jj in range(len(self._E)):
            t=d1[jj]
            tmp[jj]+=(t.sign()*f._F[t.label]).l_act_by(p**(-t.power)*mga*t.igamma(self.embed_quaternion))
        return HarmonicCocycleElement(self,tmp,from_values=True)

    def test_app_hecke(self,l,f):
        return self.__apply_hecke_operator(l,f)

    def __apply_hecke_operator(self,l,f):
        r"""
        This function applies a Hecke operator to a harmonic cocycle.

        INPUT:

          - ``l`` - an integer

          - ``f`` - a harmonic cocycle

        OUTPUT:

          A harmonic cocycle which is the result of applying the lth Hecke operator
          to f

        EXAMPLES:
        ::

        """
        R=self._R
        HeckeData,alpha=self._X._get_hecke_data(l)
        if(self.level()%l==0):
            factor=QQ(l**(Integer((self._k-2)/2))/(l+1))
        else:
            factor=QQ(l**(Integer((self._k-2)/2)))
        p=self._X._p
        alphamat=self.embed_quaternion(alpha)
        tmp=[self._U.element_class(self._U,zero_matrix(self._R,self._k-1,1),quick=True) for jj in range(len(self._E))]
        for ii in range(len(HeckeData)):
            d1=HeckeData[ii][1]
            mga=self.embed_quaternion(HeckeData[ii][0])*alphamat
            for jj in range(len(self._E)):
                t=d1[jj]
                tmp[jj]+=(t.sign()*f._F[t.label]).l_act_by(p**(-t.power)*mga*t.igamma(self.embed_quaternion))
        return HarmonicCocycleElement(self,[factor*x for x in tmp],from_values=True)

    def _compute_atkin_lehner_matrix(self,d):
        r"""
        When the underlying coefficient module is finite, this function computes the 
        matrix of an Atkin-Lehner involution in the basis provided by the function
        basis_matrix

        INPUT:

          - ``d`` - an integer dividing p*Nminus*Nplus

        OUTPUT:

          The matrix of the AL-involution at d in the basis given by self.basis_matrix

        EXAMPLES:
        ::

        """
        res=self.__compute_operator_matrix(lambda f:self.__apply_atkin_lehner(d,f))
        return res

    def _compute_hecke_matrix_prime(self,l):
        r"""
        When the underlying coefficient module is finite, this function computes the 
        matrix of a (prime) Hecke operator in the basis provided by the function
        basis_matrix

        INPUT:

          - ``l`` - an integer prime

        OUTPUT:

          The matrix of T_l acting on the cocycles in the basis given by 
          self.basis_matrix

        EXAMPLES:
        ::

        """
        res=self.__compute_operator_matrix(lambda f:self.__apply_hecke_operator(l,f))
        return res

    def __compute_operator_matrix(self,T):
        r"""
        Compute the matrix of the operator ``T``.

        EXAMPLES:
        ::

        """
        R=self._R
        A=self.basis_matrix().transpose()
        basis=self.basis()
        B=zero_matrix(R,len(self._E)*(self._k-1),self.dimension())
        for rr in range(len(basis)):
            g=T(basis[rr])
            B.set_block(0,rr,Matrix(R,len(self._E)*(self._k-1),1,[g._F[e]._val[ii,0]  for e in range(len(self._E)) for ii in range(self._k-1) ]))

        try:
            res=(A.solve_right(B)).transpose()
            res.set_immutable()
            return res
        except ValueError:
            print A
            print B
            raise ValueError

class HarmonicCocyclesSubmodule(sage.modular.hecke.submodule.HeckeSubmodule,HarmonicCocycles):
    r"""

    INPUT:

     - ``x`` - integer (default: 1) the description of the
       argument x goes here.  If it contains multiple lines, all
       the lines after the first need to be indented.

     - ``y`` - integer (default: 2) the ...

    EXAMPLES:

    AUTHOR:

    - Marc Masdeu (2012-02-20)
    """
    def __init__(self, ambient_module, submodule, dual=None, check=False):
        """
            ambient_module -- HarmonicCocycles
            submodule -- a submodule of the ambient space.
            dual_module -- (default: None) ignored
            check -- (default: False) whether to check that the
                     submodule is Hecke equivariant
        """
        A = ambient_module
        sage.modular.hecke.submodule.HeckeSubmodule.__init__(self, A, submodule, check=check)
        self.__rank=submodule.dimension()
        HarmonicCocycles.__init__(self,X=A._X,k=A._k,prec=A._prec,basis_matrix=submodule.basis_matrix()*A.basis_matrix())

    def rank(self):
        r"""
        Returns the rank (dimension) of the submodule.

        OUTPUT:

        integer - The rank of ``self``.

        """
        return self.__rank

    def _repr_(self):
        r"""
        Returns the representation of self as a string.
        """
        return "Subspace of %s of dimension %s"%(self.ambient(),self.dimension())


class pAutomorphicForm(ModuleElement):
    r"""
    This class is a rudimentary implementation of a class for a p-adic automorphic
    form on a definite quaternion algebra over Q. These are required in order to
    compute moments of measures associated to harmonic cocycles on the BT-tree
    using the overconvergent modules of Darmon-Pollack and Matt Greenberg. See
    Greenberg's thesis for more details.

    INPUT:
     - ``vec`` - 

     - ``quick`` - boolean (default: False) 

    EXAMPLES:

    This example illustrates ...

    ::

    ...

    REFERENCES:

    Matthew Greenberg's thesis (available on his webpage as of 02/12).

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu


    """
    def __init__(self,_parent,vec,quick=False):
        ModuleElement.__init__(self,_parent)
        self._parent=_parent
        self._nE=2*len(_parent._E) # We record the values at the opposite edges
        self._cached_values=dict()
        self._R=Qp(_parent._X._p,prec=_parent._prec)
        if(quick):
                self._F=[v for v in vec]
        else:
            if(isinstance(vec,pAutomorphicForm)):
                self._F=[self._parent._U(vec._F[ii]) for ii in range(self._nE)]
                self._make_invariant()

            elif(isinstance(vec,HarmonicCocycleElement)):
                assert(_parent._U.weight()==vec._wt-2)
                self._F=[]
                assert(2*len(vec._F)==self._nE)
                assert(isinstance(_parent._U,OCVn))
                E=self._parent._E


                MMM=vec._parent._U.element_class
                tmp=[]
                for ii in range(len(vec._F)):
                    newtmp=MMM(vec._parent._U,vec._F[ii]).l_act_by(E[ii].rep.inverse())
                    tmp.append(newtmp)
                    self._F.append(_parent._U(newtmp))
                A=Matrix(QQ,2,2,[0,-1/_parent._X._p,-1,0])
                for ii in range(len(vec._F)):
                    self._F.append(_parent._U(-1*tmp[ii].r_act_by(A)))
                self._make_invariant()

            elif(isinstance(vec,list) and len(vec)==self._nE):
                try:
                    self._F=[self._parent._U(v) for v in vec]
                except:
                    try:
                        veczp=_parent._U._R(vec)
                        self._parent=_parent
                        self._F=[self._parent._U(veczp) for ii in range(self._nE)]
                    except:
                        print vec
                        assert(0)
            else:
                try:
                    veczp=_parent._U._R(vec)
                    self._parent=_parent
                    self._F=[self._parent._U(veczp) for ii in range(self._nE)]
                except:
                    raise ValueError,"Cannot initialize a p-adic automorphic form with the given input="+str(vec)

    def _make_invariant(self):
        r"""
        EXAMPLES:

        ::
        """
        S=self._parent._X.get_edge_stabs()
        E=self._parent._E
        M=[e.rep for e in E]+[e.opposite.rep for e in E]
        lS=len(S)
        assert(2*len(S)==self._nE)
        MMM=self._parent._U.element_class
        for ii in range(self._nE):
            Si=S[ii%lS]
            if(any([v[2] for v in Si])):
                x=MMM(self._F[ii].parent(),self._F[ii]._val,quick=True)
                self._F[ii]=MMM(self._F[ii].parent(),0)
                s=QQ(0)
                m=M[ii]
                for v in Si:
                    s+=1
                    self._F[ii]+=x.r_act_by(m.adjoint()*self._parent.embed_quaternion(v[0])*m)
                self._F[ii]=(1/s)*self._F[ii]

    def precision(self):
        r"""
        The precision of ``self``, which is the minimum among the precision of the values on a fundamental domain.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        return min(x.precision() for x in self._F)

    def _add_(self,g):
        r"""
        This function adds two p-adic automorphic forms.
        
        INPUT:

          - ``g`` - a p-adic automorphic form

        OUTPUT:

          the result of adding g to self

        EXAMPLES:

        This example illustrates ...

        ::
        """
        #Should ensure that self and g are of the same weight and on the same curve
        vec=[self._F[e]+g._F[e] for e in range(self._nE)]
        return pAutomorphicForm(self._parent,vec,quick=True)
    
    def _sub_(self,g):
        r"""
        This function subtracts a p-adic automorphic form from another.

        INPUT:

          - ``g`` - a p-adic automorphic form

        OUTPUT:

          the result of subtracting g from self

        EXAMPLES:

        This example illustrates ...

        ::

        """
        #Should ensure that self and g are of the same weight and on the same curve
        vec=[self._F[e]-g._F[e] for e in range(self._nE)]
        return pAutomorphicForm(self._parent,vec,quick=True)

    def _getitem_(self,e1):
        r"""
        Evaluates a p-adic automorphic form on a matrix in GL2(Qp).

        INPUT:

          - ``e1`` - a matrix in GL2(Qp)

        OUTPUT:

          the value of self evaluated on e1

        EXAMPLES:

        This example illustrates ...

        ::

        """
        return self.evaluate(e1)

    def evaluate(self,e1):
        r"""
        Evaluates a p-adic automorphic form on a matrix in GL2(Qp).

        INPUT:

          - ``e1`` - a matrix in GL2(Qp)

        OUTPUT:

          the value of self evaluated on e1

        EXAMPLES:

        This example illustrates ...

        ::

        """
        X=self._parent._X
        u=DoubleCosetReduction(X,e1)
        return (self._F[u.label+u.parity*self._nE/2]).r_act_by((u.t())*X._p**(u.power))

    def _rmul_(self,a):
        r"""
        Multiplies the automorphic form by a scalar.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        #Should ensure that 'a' is a scalar
        return pAutomorphicForm(self._parent,[a*self._F[e] for e in range(self._nE)],quick=True)

    def left_act_by(self,alpha):
        r"""
        Act by ``alpha`` on the left.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        Y=self._parent._X
        E=self._parent._E
        p=Y._p
        Tf=[]
        for e in E:
            u=DoubleCosetReduction(Y,alpha*e.rep)
            r=u.t()*p**(-(u.power))
            if(u.parity==0):
                tmp=self._F[u.label].r_act_by(r)
            else:
                tmp=self._F[u.label+len(E)].r_act_by(r)
            Tf.append(tmp)
        for e in E:
            u=DoubleCosetReduction(Y,alpha*e.opposite.rep)
            r=u.t()*gg*p**(-(u.power))
            if(u.parity==0):
                tmp=self._F[u.label].r_act_by(r)
            else:
                tmp=self._F[u.label+len(E)].r_act_by(r)
            Tf.append(tmp)
        return pAutomorphicForm(self._parent,Tf,quick=True)


    def _repr_(self):
        r"""
        This returns the representation of self as a string.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        tmp='p-adic automorphic form on '+str(self._parent)+':\n'
        tmp+='   e   |   c(e)'
        tmp+='\n'
        for e in range(Integer(self._nE/2)):
            tmp+='  '+str(e)+' | '+str(self._F[e])+'\n'
        return tmp

    def valuation(self):
        r"""
        The valuation of ``self``, defined as the minimum of the valuations of the values that it takes on a set of edge representatives.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        if not self.__nonzero__():
            return oo
        else:
            return(min([self._F[e].valuation() for e in range(self._nE)]))

    def __nonzero__(self):
        r"""
        """
        return(any([self._F[e].__nonzero__() for e in range(self._nE)]))

    def improve(self, verbose = True):
        r"""
        Repeatedly applies the `U_p` operator to a p-adic automorphic form. This
        is used to compute moments of a measure associated to a rigid modular form in the
        following way: lift a rigid modular form to an ``overconvergent'' `p`-adic automorphic
        form in any way, and then repeatedly apply `U_p` to project to the ordinary part.
        The resulting form encodes the moments of the measure of the original rigid modular 
        form (assuming it is ordinary). 


        EXAMPLES:


        REFERENCES:

        For details see Matthew Greenberg's thesis (available on his webpage as of 02/12).
        Alternatively check out Darmon-Pollack for the analogous algorithm in the case of
        modular symbols.

        AUTHORS:


        - Cameron Franc (2012-02-20)
        - Marc Masdeu

        """
        MMM=self._parent
        if verbose == True:
            sys.stdout.flush()
        h2=MMM._apply_Up_operator(self,True)
        if verbose == True:
            sys.stdout.write("#")
            sys.stdout.flush()
        ii=0
        current_val=0
        old_val=-Infinity
        init_val=self.valuation()
        while(current_val>old_val):
            old_val=current_val
            ii+=1
            self._F=[self._parent._U(c) for c in h2._F]
            h2=MMM._apply_Up_operator(self,scale = True)
            current_val=(h2-self).valuation()-init_val
            # print 'val =',current_val
            if current_val is Infinity:
                break
            if verbose == True:
                sys.stdout.write("#")
                sys.stdout.flush()
        if verbose == True:
            print ''
        self._F=[self._parent._U(c) for c in h2._F]

    def integrate(self,f,center=1,level=0,method='moments'):
        r"""
        Calculate
        .. MATH::

            \int_{\PP^1(\QQ_p)} f(x)d\mu(x)

        were `\mu` is the measure associated to ``self``.

        INPUT:

         - ``f`` - An analytic function.

         - ``center`` - 2x2 matrix over Qp (default: 1)

         - ``level`` - integer (default: 0)

         - ``method`` - string (default: 'moments'). Which method of integration to use. Either 'moments' or 'riemann_sum'.


        EXAMPLES:

        ::


        AUTHORS:

        - Marc Masdeu (2012-02-20)
        - Cameron Franc (2012-02-20)

        """
        E=self._parent._X._BT.get_balls(center,level)
        R1=LaurentSeriesRing(f.base_ring(),'r1')
        R2=PolynomialRing(f.base_ring(),'r2')
        value=0
        ii=0
        if(method=='riemann_sum'):
            R1.set_default_prec(self._parent._U.weight()+1)
            for e in E:
                ii+=1
                #print ii,"/",len(E)
                exp=((R1([e[1,1],e[1,0]]))**(self._parent._U.weight())*e.determinant()**(-(self._parent._U.weight())/2))*f(R1([e[0,1],e[0,0]])/R1([e[1,1],e[1,0]]))
                #exp=R2([tmp[jj] for jj in range(self._parent._k-1)])
                new=self.evaluate(e).evaluate(exp.truncate(self._parent._U.weight()+1))
                value+=new
        elif(method=='moments'):
            R1.set_default_prec(self._parent._U.dimension())
            for e in E:
                ii+=1
                #print ii,"/",len(E)
                tmp=((R2([e[1,1],e[1,0]]))**(self._parent._U.weight())*e.determinant()**(-(self._parent._U.weight())/2))*f(R2([e[0,1],e[0,0]])/R2([e[1,1],e[1,0]]))
                exp=R1(tmp.numerator())/R1(tmp.denominator())
                new=self.evaluate(e).evaluate(exp)
                value+=new
        else:
            print 'The available methods are either "moments" or "riemann_sum". The latter is only provided for consistency check, and should never be used.'
            return False
        return value

    def modular_form(self,z=None,level=0,method='moments'):
        r"""
        Returns the modular form corresponding to ``self``.
        
        INPUT:

          - ``z`` - (default: None). If specified, returns the value of the form at the point ``zz`` in the `p`-adic upper half plane.

          - ``level`` - integer (default: 0). If ``method`` is 'riemann_sum', will use a covering of `\PP^1(\QQ_p)` with balls of size `p^-\mbox{level]`.
        
          - ``method`` - string (default: ``moments``). It must be either ``moments`` or ``riemann_sum``.

        OUTPUT:

        A function from the `p`-adic upper half plane to `\CC_p`. If an argument ``z`` was passed, returns instead the value at that point.

        """
        return self.derivative(z,level,method,order=0)

    def derivative(self,z=None,level=0,method='moments',order=1):
        r"""
        Returns the derivative of the modular form corresponding to
        ``self``.

        INPUT:

          - ``z`` - (Default: None). If specified, evaluates the derivative
              at the point ``z`` in the `p`-adic upper half plane.

          - ``level`` - integer (default: 0). If ``method`` is 'riemann_sum', will use a covering of `\PP^1(\QQ_p)` with balls of size `p^-\mbox{level]`.
        
          - ``method`` - string (default: ``moments``). It must be either ``moments`` or ``riemann_sum``.

          - ``order`` - integer (Default: 1). The order of the derivative to be computed.

        OUTPUT:

        A function from the `p`-adic upper half plane to `\CC_p`. If an argument ``z`` was passed, returns instead the value of the derivative at that point.

        EXAMPLES:

        """
        def F(z,level=0,method='moments'):
            R=PolynomialRing(z.parent(),'x,y').fraction_field()
            Rx=PolynomialRing(z.parent(),'x1').fraction_field()
            x1=Rx.gen()
            subst=R.hom([x1,z],codomain=Rx)
            x,y=R.gens()
            center=self._parent._X._BT.find_containing_affinoid(z)
            zbar=z.trace()-z
            f=R(1)/(x-y)
            k=self._parent._n+2
            V=[f]
            for ii in range(order):
                V=[v.derivative(y) for v in V]+[k/(y-zbar)*v for v in V]
                k+=2
            return sum([self.integrate(subst(v),center,level,method) for v in V])
        if z is None:
            return F

        try:
            try: tmp=self._cached_values[z,level,method,order,self.precision()]
            except KeyError:
                tmp=F(z,level,method)
                self._cached_values[z,level,method,order,self.precision()]=tmp
            return tmp
        except TypeError:
            return F(z,level,method)


    # So far we can't break it into two integrals because of the pole at infinity.
    def coleman(self,t1,t2,E=None,method='moments',mult=False,delta=-1,level=0):
        r"""
        If ``self`` is a `p`-adic automorphic form that corresponds to a rigid modular form,
        then this computes the coleman integral of this form between two points on
        the boundary `\PP^1(\QQ_p)` of the `p`-adic upper half plane.

        INPUT:

         - ``t1``, ``t2`` - elements of `\PP^1(\QQ_p)` (the endpoints of integration)

         - ``E`` - (Default: None). If specified, will not compute the covering adapted
            to ``t1`` and ``t2`` and instead use the given one. In that case, ``E``
            should be a list of matrices corresponding to edges describing the open
            balls to be considered.

         - ``method`` - string (Default: 'moments'). Tells which algorithm to use 
           (alternative is 'riemann_sum', which is unsuitable for computations 
           requiring high precision)

         - ``mult`` - boolean (Default: False). Whether to use the multiplicative version.

         - ``delta`` - integer (Default: -1)

         - ``level`` - integer (Default: 0)

        OUTPUT:

          The result of the coleman integral

        EXAMPLES:

        This example illustrates ...

        ::

        TESTS:

        ::
            sage: p = 7
            sage: lev = 2
            sage: prec = 20
            sage: X = BTQuotient(p,lev, use_magma = True) # optional - magma
            sage: k = 2 # optional - magma
            sage: M = HarmonicCocycles(X,k,prec) # optional - magma
            sage: B = M.basis() # optional - magma
            sage: f = 3*B[0] # optional - magma
            sage: MM = pAutomorphicForms(X,k,prec,overconvergent=True) # optional - magma
            sage: D = -11 # optional - magma
            sage: X.is_admissible(D) # optional - magma
            True
            sage: K.<a> = QuadraticField(D) # optional - magma
            sage: CM = X.get_CM_points(D,prec=prec) # optional - magma
            sage: Kp = CM[0].parent() # optional - magma
            sage: P=CM[0] # optional - magma
            sage: Q=P.trace()-P # optional - magma
            sage: F=MM.lift(f, verbose = False) # long time optional - magma
            sage: J0=F.coleman(P,Q,mult=True) # long time optional - magma
            sage: E=EllipticCurve([1,0,1,4,-6]) # optional - magma
            sage: T=E.tate_curve(p) # optional - magma
            sage: xx,yy=getcoords(T,J0,prec) # long time optional -magma
            sage: P = E.base_extend(K).lift_x(algdep(xx,1).roots(QQ)[0][0]); P # long time optional - magma
            (7/11 : 58/121*a - 9/11 : 1)

            sage: p = 13 # optional - magma
            sage: lev = 2 # optional - magma
            sage: prec = 20 # optional - magma
            sage: Y = BTQuotient(p,lev, use_magma = True) # optional - magma
            sage: k = 2 # optional - magma
            sage: M=HarmonicCocycles(Y,k,prec) # optional - magma
            sage: B=M.basis() # optional - magma

            sage: f = B[1] # optional - magma
            sage: g = -4*B[0]+3*B[1] # optional - magma
            sage: MM = pAutomorphicForms(Y,k,prec,overconvergent=True) # optional - magma
            sage: D = -11 # optional - magma
            sage: Y.is_admissible(D) # optional - magma
            True
            sage: K.<a> = QuadraticField(D) # optional - magma
            sage: CM = Y.get_CM_points(D,prec=prec) # optional - magma
            sage: Kp = parent(CM[0]) # optional - magma
            sage: P = CM[0] # optional - magma
            sage: Q = P.trace()-P # optional - magma
            sage: F = MM.lift(f, verbose = False) # long time optional - magma
            sage: J11 = F.coleman(P,Q,mult = True) # long time optional - magma
            sage: E = EllipticCurve('26a2') # optional - magma
            sage: T = E.tate_curve(p) # optional - magma
            sage: xx,yy = getcoords(T,J11,prec) # long time optional - magma
            sage: HP = E.base_extend(K).lift_x(algdep(xx,1).roots(QQ)[0][0]); HP # long time optional - magma
            (-137/11 : 2/121*a + 63/11 : 1)

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu (2012-02-20)
        """
        if(mult and delta>=0):
            raise NotImplementedError, "Need to figure out how to implement the multiplicative part."
        p=self._parent._X._p
        K=t1.parent()
        R=PolynomialRing(K,'x')
        x=R.gen()
        R1=LaurentSeriesRing(K,'r1')
        r1=R1.gen()
        if(E is None):
            E=self._parent._X._BT.find_covering(t1,t2)
            # print 'Got %s open balls.'%len(E)
        value=0
        ii=0
        value_exp=K(1)
        if(method=='riemann_sum'):
            R1.set_default_prec(self._parent._U.weight()+1)
            for e in E:
                ii+=1
                b=e[0,1]
                d=e[1,1]
                y=(b-d*t1)/(b-d*t2)
                poly=R1(y.log()) #R1(our_log(y))
                c_e=self.evaluate(e)
                new=c_e.evaluate(poly)
                value+=new
                if mult:
                    value_exp *= K.teichmuller(y)**Integer(c_e[0].rational_reconstruction())

        elif(method=='moments'):
            R1.set_default_prec(self._parent._U.dimension())
            for e in E:
                ii+=1
                f=(x-t1)/(x-t2)
                a,b,c,d=e.list()
                y0=f(R1([b,a])/R1([d,c])) #f( (ax+b)/(cx+d) )
                y0=p**(-y0(0).valuation())*y0
                mu=K.teichmuller(y0(0))
                y=y0/mu-1
                poly=R1(0)
                ypow=y
                for jj in range(1,R1.default_prec()+10):
                    poly+=(-1)**(jj+1)*ypow/jj
                    ypow*=y
                if(delta>=0):
                    poly*=((r1-t1)**delta*(r1-t2)**(self._parent._n-delta))
                c_e=self.evaluate(e)
                new=c_e.evaluate(poly)
                value+=new
                if mult:
                    value_exp *= K.teichmuller((b-d*t1)/(b-d*t2))**Integer(c_e[0].rational_reconstruction())

        else:
            print 'The available methods are either "moments" or "riemann_sum". The latter is only provided for consistency check, and should not be used in practice.'
            return False
        if mult:
            return value_exp * value.exp()
        return value


class pAutomorphicForms(Module):
    Element=pAutomorphicForm
    r"""
    The module of (quaternionic) `p`-adic automorphic forms.

    EXAMPLES:

    ::


    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu (2012-02-20)
    """
    def __init__(self,X,U,prec=None,t=None,R=None,overconvergent=False):
        if(R is None):
            if(not isinstance(U,Integer)):
                self._R=U.base_ring()
            else:
                if(prec is None):
                    prec=100
                self._R=Qp(X._p,prec)
        else:
            self._R=R
        #U is a CoefficientModuleSpace
        if(isinstance(U,Integer)):
            if(t is None):
                if(overconvergent):
                    t=prec-U+1
                else:
                    t=0
            self._U=OCVn(U-2,self._R,U-1+t)
        else:
            self._U=U
        self._X=X
        self._V=self._X.get_vertex_list()
        self._E=self._X.get_edge_list()
        self._prec=self._R.precision_cap()
        self._n=self._U.weight()
        Module.__init__(self,base=self._R)
        self._populate_coercion_lists_()

    def _repr_(self):
        r"""
        This returns the representation of self as a string.

        EXAMPLES:

        This example illustrates ...

        ::
        """
        s='Space of automorphic forms on '+str(self._X)+' with values in '+str(self._U)
        return s

    def _coerce_map_from_(self, S):
        r"""
        Can coerce from other HarmonicCocycles or from pAutomorphicForms
        """
        if(isinstance(S,HarmonicCocycles)):
            if(S._k-2!=self._n):
                return False
            if(S._X!=self._X):
                return False
            return True
        if(isinstance(S,pAutomorphicForms)):
            if(S._n!=self._n):
                return False
            if(S._X!=self._X):
                return False
            return True
        return False

    def _element_constructor_(self,x):
        r"""
        """
        #Code how to coherce x into the space
        #Admissible values of x?
        if isinstance(x,(HarmonicCocycleElement,pAutomorphicForm)):
            return pAutomorphicForm(self,x)

    def _an_element_(self):
        r"""
        Returns an element of the module.
        """
        return pAutomorphicForm(self,1)

    def embed_quaternion(self,g):
        r"""
        Returns the image of ``g`` under the embedding
        of the quaternion algebra into 2x2 matrices.
        """
        return self._X.embed_quaternion(g,prec = self._prec)

    def lift(self,f, verbose = True):
        r"""
        Lifts the harmonic cocycle ``f`` to an
        overconvegent automorphic form, thus computing
        all the moments.
        """
        F=self.element_class(self,f)
        F.improve(verbose = verbose)
        return F

    def _apply_Up_operator(self,f,scale=False, fix_lowdeg_terms = True):
        r"""
        Apply the Up operator to ``f``.

        EXAMPLES:

        ::
            sage: X = BTQuotient(3,11)
            sage: M = HarmonicCocycles(X,4,30)
            sage: A = pAutomorphicForms(X,4,10, overconvergent = True)
            sage: F = A.lift(M.basis()[0], verbose = False); F
            p-adic automorphic form on Space of automorphic forms on Quotient of the Bruhat Tits tree of GL_2(QQ_3) with discriminant 11 and level 1 with values in Overconvergent coefficient module of weight n=2 over the ring 3-adic Field with capped relative precision 10 and depth 10:
            e   |   c(e)
            0 | 3^2 + O(3^12) + O(3^32)*z + O(3^26)*z^2 + (2*3^2 + 3^3 + 2*3^5 + 3^7 + 3^8 + 2*3^9 + O(3^10))*z^3 + (2*3^5 + 2*3^6 + 2*3^7 + 3^9 + O(3^10))*z^4 + (3^2 + 3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 3^7 + 2*3^8 + 3^9 + O(3^10))*z^5 + (3^2 + 2*3^3 + 3^4 + 2*3^6 + O(3^10))*z^6 + (2*3^3 + 3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 3^8 + 3^9 + O(3^10))*z^7 + (3^2 + 3^3 + 2*3^6 + 3^7 + 3^8 + 3^9 + O(3^10))*z^8 + (2*3^2 + 2*3^3 + 2*3^5 + 2*3^7 + 3^8 + 2*3^9 + O(3^10))*z^9
            1 | 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + O(3^12) + (3^2 + O(3^12))*z + (2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + O(3^12))*z^2 + (2*3^2 + 2*3^3 + 3^4 + 2*3^5 + 3^6 + 2*3^7 + 2*3^8 + O(3^10))*z^3 + (2*3^3 + 3^5 + 3^7 + 3^8 + O(3^10))*z^4 + (2*3^3 + 3^6 + 3^7 + 3^9 + O(3^10))*z^5 + (3^3 + 2*3^4 + 2*3^5 + 2*3^7 + 3^8 + 3^9 + O(3^10))*z^6 + (3^7 + 2*3^8 + 2*3^9 + O(3^10))*z^7 + (3^3 + 2*3^4 + 3^7 + O(3^10))*z^8 + (2*3^2 + 3^4 + 3^6 + 2*3^7 + 3^8 + 2*3^9 + O(3^10))*z^9
            2 | 3^2 + 2*3^3 + 2*3^6 + 3^7 + 2*3^8 + O(3^12) + (3 + 2*3^2 + 2*3^3 + 3^5 + 2*3^6 + 3^7 + 3^10 + O(3^11))*z + (2*3 + 2*3^2 + 3^4 + 2*3^5 + 2*3^6 + 2*3^8 + 3^10 + O(3^11))*z^2 + (2*3 + 3^2 + 2*3^7 + 3^9 + O(3^10))*z^3 + (3 + 2*3^2 + 2*3^4 + 3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10))*z^4 + (3 + 3^2 + 3^4 + 2*3^9 + O(3^10))*z^5 + (3^3 + 2*3^5 + 3^6 + 3^8 + 2*3^9 + O(3^10))*z^6 + (3^5 + 2*3^7 + 3^9 + O(3^10))*z^7 + (2*3 + 3^3 + 3^4 + 2*3^6 + O(3^10))*z^8 + (2*3 + 2*3^3 + 2*3^4 + 2*3^6 + O(3^10))*z^9
            3 | 3^2 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + O(3^12) + (3^3 + 2*3^4 + 2*3^8 + O(3^13))*z + (3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^6 + 3^7 + 2*3^8 + 2*3^9 + 2*3^10 + O(3^11))*z^2 + (3^2 + 2*3^3 + 3^4 + 3^7 + 3^8 + 2*3^9 + O(3^10))*z^3 + (3 + 2*3^2 + 2*3^3 + 3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 3^8 + O(3^10))*z^4 + (3 + 3^3 + 3^4 + 2*3^5 + 2*3^6 + 3^7 + 2*3^8 + 2*3^9 + O(3^10))*z^5 + (3 + 3^4 + 3^5 + 3^6 + 2*3^7 + O(3^10))*z^6 + (2*3 + 3^2 + 2*3^3 + 3^4 + 2*3^6 + 3^8 + 3^9 + O(3^10))*z^7 + (3 + 3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 3^7 + 2*3^9 + O(3^10))*z^8 + (2*3^2 + 3^4 + 3^5 + 3^8 + 3^9 + O(3^10))*z^9

        """
        HeckeData=self._X._get_Up_data()
        if scale == False:
            factor=(self._X._p)**(self._U.weight()/2)
        else:
            factor=1
        Tf=[]
        for jj in range(2*len(self._E)):
            tmp=self._U(0)
            for d in HeckeData:
                gg=d[0] # acter
                u=d[1][jj] # edge_list[jj]
                r=self._X._p**(-(u.power))*(u.t()*gg)
                tmp+=f._F[u.label+u.parity*len(self._E)].r_act_by(r)
            tmp *= factor
            for ii in range(self._n+1):
                tmp[ii] = f._F[jj][ii]
            Tf.append(tmp)
        return pAutomorphicForm(self,Tf,quick=True)

#########################################################################
#       Copyright (C) 2011 Cameron Franc and Marc Masdeu
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

from sage.structure.element import ModuleElement
from sage.modules.module import Module
from sage.matrix.constructor import Matrix
from sage.matrix.matrix_space import MatrixSpace_generic
from sage.matrix.matrix_space import MatrixSpace
from copy import copy
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.power_series_ring import PowerSeriesRing
from sage.structure.unique_representation import UniqueRepresentation

class OCVnElement(ModuleElement):
    r"""
    This class represents elements in an overconvergent coefficient module.

    INPUT:

     - ``parent`` - An overconvergent coefficient module.

     - ``val`` - The value that it needs to store (default: 0). It can be another OCVnElement,
       in which case the values are copied. It can also be a column vector (or something
       coercible to a column vector) which represents the values of the element applied to
       the polynomials `1`, `x`, `x^2`, ... ,`x^n`.

     - ``quick`` - boolean (default: False). If set to true, no checks are done and ``val`` is
       assumed to be the a column vector.

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu (2012-02-20)
    """
    def __init__(self,parent,val=0,quick=False):
        ModuleElement.__init__(self,parent)
        self._parent=parent
        self._n=self._parent._n
        self._nhalf=Integer(self._n/2)
        self._depth=self._parent._depth
        if quick:
            self._val=copy(val)
        else:
            if isinstance(val,self.__class__):
                d=min([val._parent._depth,parent._depth])
                assert(val._parent.weight()==parent.weight())
                self._val=Matrix(self._parent._R,self._depth,1,0)
                for ii in range(d):
                    self._val[ii,0]=val._val[ii,0]
            else:
                try:
                    self._val=MatrixSpace(self._parent._R,self._depth,1)(val)
                except:
                    self._val=val*ones_matrix(self._parent._R,self._depth,1)

    def __getitem__(self,r):
        r"""
        Returns the value of ``self`` on the polynomial `x^r`.

        INPUT:
          - ``r`` - an integer. The power of `x`. 

        EXAMPLES:

        """
        return self._val[r,0]

    def __setitem__(self,r, val):
        r"""
        Sets the value of ``self`` on the polynomial `x^r` to ``val``.

        INPUT:
        - ``r`` - an integer. The power of `x`.
        - ``val`` - a value.

        EXAMPLES:

        """
        self._val[r,0] = val

    def element(self):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        tmp=self.matrix_rep()
        return [tmp[ii,0] for ii in range(tmp.nrows())]

    def list(self):
        r"""
        EXAMPLES:

        This example illustrates ...

        ::
        """
        return self.element()

    def matrix_rep(self,B=None):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        #Express the element in terms of the basis B
        if(B is None):
            B=self._parent.basis()
        A=Matrix(self._parent._R,self._parent.dimension(),self._parent.dimension(),[[b._val[ii,0] for b in B] for ii in range(self._depth)])
        tmp=A.solve_right(self._val)
        return tmp

    def _add_(self,y):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        val=self._val+y._val
        return self.__class__(self._parent,val,quick=True)

    def _sub_(self,y):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        val=self._val-y._val
        return self.__class__(self._parent,val,quick=True)

    def l_act_by(self,x):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        #assert(x.nrows()==2 and x.ncols()==2) #An element of GL2
        return self._l_act_by(x[0,0],x[0,1],x[1,0],x[1,1],extrafactor=x.determinant()**(-self._nhalf))

    def r_act_by(self,x):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        #assert(x.nrows()==2 and x.ncols()==2) #An element of GL2
        return self._l_act_by(x[1,1],-x[0,1],-x[1,0],x[0,0],extrafactor=x.determinant()**(-self._nhalf))

    def _l_act_by(self,a,b,c,d,extrafactor=1):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        R=self._parent._R
        if(self._parent.base_ring().is_exact()):
            factor=1
        else:
            t=min([R(x).valuation() for x in [a,b,c,d] if x!=0])
            factor=R.prime()**(-t)
        try:
            x=self._parent._powers[(factor*a,factor*b,factor*c,factor*d)]
            return self.__class__(self._parent,(extrafactor*factor**(-self._n))*(x*self._val),quick=True)
        except KeyError:
            tmp=self._parent._get_powers_and_mult(factor*a,factor*b,factor*c,factor*d,extrafactor*factor**(-self._n),self._val)

            return self.__class__(self._parent,tmp,quick=True)

    def _rmul_(self,a):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        #assume that a is a scalar
        return self.__class__(self._parent,a*self._val,quick=True)

    def precision_absolute(self):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        #This needs to be thought more carefully...
        if not self._parent.base_ring().is_exact():
            return [self._val[ii,0].precision_absolute() for ii in range(self._depth)]
        else:
            return oo

    def precision(self):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        #This needs to be thought more carefully...
        if not self._parent.base_ring().is_exact():
            return min([self._val[ii,0].precision_absolute() for ii in range(self._depth)])
        else:
            return oo

    def precision_relative(self):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        #This needs to be thought more carefully...
        if not self._parent.base_ring().is_exact():
            return [self._val[ii,0].precision_relative() for ii in range(self._depth)]
        else:
            return oo

    def _repr_(self):
        r"""
        This returns the representation of self as a string.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        R=PowerSeriesRing(self._parent._R,default_prec=self._depth,name='z')
        z=R.gen()
        s=str(sum([R(self._val[ii,0]*z**ii) for ii in range(self._depth)]))
        return s

    def __cmp__(self,other):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        return cmp(self._val,other._val)

    def __nonzero__(self):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        return self._val!=0

    def evaluate(self,P):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        p=self._parent._R.prime()
        try:
            r=min([P.degree()+1,self._depth])
            return sum([self._val[ii,0]*P[ii] for ii in range(r)])
        except:
            return self._val[0,0]*P

    def valuation(self,l=None):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        if not self._parent.base_ring().is_exact():
            if(not l is None and l!=self._parent._R.prime()):
                raise ValueError, "This function can only be called with the base prime"
            return min([self._val[ii,0].valuation() for ii in range(self._depth)])
        else:
            return min([self._val[ii,0].valuation(l) for ii in range(self._depth)])


class OCVn(Module,UniqueRepresentation):
    Element=OCVnElement
    r"""
    This class represents objects in the overconvergent approximation modules used to
    describe overconvergent p-adic automorphic forms. 

    INPUT:

     - ``n`` - integer 

     - ``R`` - ring

     - ``depth`` - integer (Default: None)

     - ``basis`` - (Default: None)


    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu (2012-02-20)
    """
    def __init__(self,n,R,depth=None,basis=None):
        Module.__init__(self,base=R)
        if basis is not None:
            self._basis=copy(basis)
        self._n=n
        self._R=R
        if R.is_exact():
            self._Rmod=self._R
        else:
            self._Rmod=Zmod(self._R.prime()**(self._R.precision_cap()))

        if depth is None:
            depth=n+1
        if depth != n+1:
            if R.is_exact(): raise ValueError, "Trying to construct an over-convergent module with exact coefficients, how do you store p-adics ??"
        self._depth=depth
        self._PowerSeries=PowerSeriesRing(self._Rmod,default_prec=self._depth,name='z')
        self._powers=dict()
        self._populate_coercion_lists_()

    def _an_element_(self):
        r"""
        """
        return OCVnElement(self,Matrix(self._R,self._depth,1,range(1,self._depth+1)),quick=True)

    def _coerce_map_from_(self, S):
        r"""

        EXAMPLES:

        ::

        """
        # Nothing coherces here, except OCVnElement
        return False

    def _element_constructor_(self,x):
        r"""

        EXAMPLES:

        """
        #Code how to coherce x into the space
        #Admissible values of x?
        return OCVnElement(self,x)

    def _get_powers_and_mult(self,a,b,c,d,lambd,vect):
        r"""
        Compute the action of a matrix on the basis elements.

        EXAMPLES:

        ::

        """
        R=self._PowerSeries
        r=R([b,a])
        s=R([d,c])
        n=self._n
        if(self._depth==n+1):
            rpows=[R(1)]
            spows=[R(1)]
            for ii in range(n):
                rpows.append(r*rpows[ii])
                spows.append(s*spows[ii])
            x=Matrix(self._Rmod,n+1,n+1,0)
            for ii in range(n+1):
                y=rpows[ii]*spows[n-ii]
                for jj in range(self._depth):
                    x[ii,jj]=y[jj]
        else:
            ratio=r*(s**(-1))
            y=s**n
            x=Matrix(self._Rmod,self._depth,self._depth,0)
            for jj in range(self._depth):
                x[0,jj]=y[jj]
            for ii in range(1,self._depth):
                y*=ratio
                for jj in range(self._depth):
                    x[ii,jj]=y[jj]
        if self._Rmod is self._R:
            xnew=x
        else:
            xnew=x.change_ring(self._R.base_ring())
            xnew=xnew.change_ring(self._R)
        self._powers[(a,b,c,d)]=xnew
        return self._R(lambd)*xnew*vect

    def _repr_(self):
        r"""
        This returns the representation of self as a string.

        EXAMPLES:

        """
        s='Overconvergent coefficient module of weight n='+str(self._n)+' over the ring '+ str(self._R)+' and depth '+str(self._depth)
        return s

    def basis(self):
        r"""
        A basis of the module.

        INPUT:

         - ``x`` - integer (default: 1) the description of the
           argument x goes here.  If it contains multiple lines, all
           the lines after the first need to be indented.

         - ``y`` - integer (default: 2) the ...

        OUTPUT:

        integer -- the ...

        EXAMPLES:


        """
        try: return self._basis
        except: pass
        self._basis=[OCVnElement(self,Matrix(self._R,self._depth,1,{(jj,0):1},sparse=False),quick=True) for jj in range(self._depth)]
        return self._basis

    def base_ring(self):
        r"""
        This function returns the base ring of the overconvergent element.

        EXAMPLES::

        This example illustrates ...

        ::

        """
        return self._R

    def depth(self):
        r"""
        Returns the depth of the module.
        """
        return self._depth

    def dimension(self):
        r"""
        Returns the dimension (rank) of the module.
        """
        return self._depth

    def weight(self):
        r"""
        Returns the cohomological weight of the automorphic form.
        """
        return self._n

    def l_matrix_representation(self,g,B=None):
        r"""
        Matrix representation of ``g`` in a given basis.

        """
        if B is None:
            B=self.basis()
        A=[(b.l_act_by(g)).matrix_rep(B) for b in B]
        d=self.dimension()
        return Matrix(self._R,d,d,[A[jj][ii,0] for ii in range(d) for jj in range(d)])


#########################################################################
#       Copyright (C) 2011 Cameron Franc and Marc Masdeu
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################


from itertools import product,chain
from sage.rings.all import Qp

def getcoords(E,u,prec=20):
    q = E.parameter(prec=prec)
    un = u * q**(-(u.valuation()/q.valuation()).floor())
    precn = (prec/q.valuation()).floor() + 4

    # formulas in Silverman II (Advanced Topics in the Arithmetic of Elliptic curves, p. 425)

    xx = un/(1-un)**2 + sum( [q**n*un/(1-q**n*un)**2 + q**n/un/(1-q**n/un)**2-2*q**n/(1-q**n)**2 for n in range(1,precn) ])

    yy = un**2/(1-un)**3 + sum( [q**(2*n)*un**2/(1-q**n*un)**3 - q**n/un/(1-q**n/un)**3+q**n/(1-q**n)**2 for n in range(1,precn) ])

    C,r,s,t = E._inverse_isomorphism(prec=prec)
    C2 = C**2
    return ( r + C2 * xx, t + s * C2 * xx + C * C2 * yy )


def our_sqrt(x,K):
    if(x==0):
        return x
    x=K(x)
    p=K.base_ring().prime()
    z=K.gen()
    found=False
    for a,b in product(range(p),repeat=2):
        y0=a+b*z
        if((y0**2-x).valuation()>0):
            found=True
            break
    y1=y0
    y=0
    while(y!=y1):
        y=y1
        y1=(y**2+x)/(2*y)
    return y

def our_log(x,prec=None):
    K=x.parent()
    if prec is None:
        prec=K.precision_cap()+10
    x0=x.unit_part()
    y=x0/K.teichmuller(x0)-1
    tmp=K(0)
    ypow=y
    for ii in range(1,prec+1):
        tmp+=(-1)**(ii+1)*ypow/ii
        ypow*=y
    return tmp

def our_exp(x,prec=None):
    K=x.parent()
    if prec is None:
        prec=K.precision_cap()+10
    tmp=K(1+x)
    xpow=x**2
    iifact=2
    for ii in range(3,prec):
        tmp+=xpow/iifact
        xpow*=x
        iifact*=ii
    return tmp


def fix_deg_monomials(v,n):
    return [reduce(lambda x,y:x*y,[v[ii]**(part[ii]-1) for ii in range(len(v))]) for part in OrderedPartitions(len(v)+n,len(v))]


#The list of elements elts must be in the form [a1,a1^-1,a2,a2^{-1}, etc]
def free_group_words(elts,op=None,init=[1]):
    if op is None:
        op=lambda x,y:x*y
    allwords=[]

    ii=0
    n=1
    # Generate words of length 1
    for i in range(len(elts)):
        wd=[i,op(elts[i],init),[i]]
        ii+=1
        if ii%10000==0:
            print ii
        yield wd[1]
        #yield wd[1],n,wd[2]
        allwords.append(wd)

    # Generate longer words
    while True:
        n+=1
        newwords = []
        for pairs in allwords:
            leftind = pairs[0]
            if leftind % 2 == 0:
                omit = leftind+1
            else:
                omit = leftind-1
            for i in range(omit)+range(omit+1,len(elts)):
                wd=[i,op(elts[i],pairs[1]),[i]+pairs[2]]
                ii+=1
                if ii%10000==0:
                    print ii
                yield wd[1]
                #yield wd[1],n,wd[2]
                newwords.append(wd)
        allwords=newwords


#Act by a fractional linear transformation on an element of the p-adic upper half plane
# The parameter twist corresponds to applying a change of variables given by the
# matrix [1,0,twist,1]
def act_by_flt(g,Z,twist = 0):
    bb=g[0,1]
    btwist=bb*twist
    aa, dd=g[0,0]+btwist,g[1,1]-btwist
    cc=g[1,0]+(g[1,1]-aa)*twist
    try:
        return [(aa*z + bb)/(cc*z + dd) for z in Z]
    except TypeError:
        return (aa*Z + bb)/(cc*Z + dd)


def get_action_flt(twist):
    return lambda g,Z:act_by_flt(g,Z,twist)

def find_good_monomial(f):
    d=max(f.degrees())
    for x in f.parent().gens():
        x2d=x**d
        print 'Trying monomial ',x
        print 'Appears in degree',f.degree(x)
        print 'and the other deg is',(f-f.coefficient(x2d)*x2d).degree(x)

        if f.degree(x)>0 and (f-f.coefficient(x2d)*x2d).degree(x)==0:
            return x2d
    return None

# Finds relations among the modular forms in X
# Up to a given degree
def find_relations(X,dmax,prec,generators,h=0):
    genus=len(X)
    p=X[0].parent()._X.prime()
    K.<g>=Qq(p^2,prec)
    g=K.gen()
    max_num_monomials=binomial(genus+dmax-1,dmax)

    sys.stdout.flush()
    CEP=[]
    for ii in range(max_num_monomials+h):
        Pt=g+p*ii
        sys.stdout.write("#")
        sys.stdout.flush()
        CEP.append([f.modular_form(Pt) for f in X])

    V=[]
    for d in range(2,dmax+1):
        num_monomials=binomial(genus+d-1,d)
        A=Matrix(K,num_monomials+h,num_monomials,[fix_deg_monomials(CEP[ii][:num_monomials],d) for ii in range(num_monomials+h)])
        for v in V:
            # Find a suitable monomial to cancel higher degrees
            d0=v[0]
            f0=sum([x[0] for x in v[1]])
            xi2d=find_good_monomial(f0)
            assert not xi2d is None
            tmons=fix_deg_monomials(generators,d-d0)
            degdmons=fix_deg_monomials(generators,d)
            pos=[(xi2d*t,degdmons.index(xi2d*t)) for t in tmons]
            A=A.stack(Matrix(K,len(pos),num_monomials,dict([((ii,pos[ii][1]),1) for ii in range(len(pos))])))
        B=A.right_kernel().matrix()
        assert(B.nrows()==1)
        mons=fix_deg_monomials(generators,d)
        tmp=B.row(0)
        newV=filter(lambda x:x[1]!=0,zip(mons,tmp))
        print 'newV=',newV
        V.append((d,newV))
    return V


def find_invariants(genus,V,P):
    generators=P.gens()
    goodMons=list(chain.from_iterable([v[1] for v in V]))
    assert all([x[1]!=0 for x in goodMons])

    A=copy(Matrix(ZZ,len(goodMons),genus,[tuple(x[0].degrees()) for x in goodMons]).kernel().matrix())

    n_invariants=A.nrows()
    goodcols=[]

    # Try to select columns to become dependent variables
    for ii in range(A.nrows()):
        found=False
        for jj in range(A.ncols()):
            if ZZ(A[ii,jj]).abs()==1 and all([all([A[i1,jj]*A[i1,j1]==0 for j1 in goodcols]) for i1 in range(ii+1,A.nrows())]):
                goodcols.append(jj)
                found=True
                break
        if not found: raise RuntimeError
        A.rescale_row(ii,A[ii,jj])
        assert(A[ii,jj]==1)
        for i0 in range(ii)+range(ii+1,A.nrows()):
            A.add_multiple_of_row(i0,ii,-A[i0,jj])

    badcols=range(A.ncols())
    for x in goodcols:
        badcols.remove(x)

    ################
    # Just to gather more information
    print 'goodcols=',goodcols
    print 'badcols=',badcols
    for ii in range(A.nrows()):
        r=A.row(ii)
        tmp=1
        for jj in range(A.ncols()):
            if(A[ii,jj]!=0):
                tmp*=goodMons[jj][1]**ZZ(A[ii,jj])
                if jj<5:
                    print 'a%s^(%s)'%(jj,ZZ(A[ii,jj])),
                else:
                    print 'b%s^(%s)'%(jj-5,ZZ(A[ii,jj])),
        print ''
        rat=algdep(tmp,1).roots(RationalField())[0][0]
        print 'rat=',rat
    ################

    S0=PolynomialRing(QQ,genus,names='a')
    S=S0.fraction_field()
    lst=[]
    for j0 in range(A.ncols()):
        try: lst.append(S.gen(badcols.index(j0)))
        except ValueError:
            ii=goodcols.index(j0)
            r=A.row(ii)
            tmp=1
            mon=1
            for jj in range(A.ncols()):
                if(A[ii,jj]!=0):
                    tmp*=goodMons[jj][1]**ZZ(A[ii,jj])
                    if jj!=j0:
                        mon*=S.gen(badcols.index(jj))**(-ZZ(A[ii,jj]))
            rat=algdep(tmp,1).roots(RationalField())[0][0]
            lst.append(S(rat*mon))
    PolyS=P.change_ring(S)
    F=[]
    ii=0
    for d,v in V:
        f=PolyS(0)
        for x in filter(lambda x:x[1]!=0,v):
            f+=PolyS(lst[ii])*PolyS(x[0])
            ii+=1
        F.append(f*f.denominator())
    PolyS0=P.change_ring(S0)
    return [PolyS0(f) for f in F]

def substitute(F,**args):
    R=F[0].parent()
    tmp=[R(f.subs(**args)) for f in F]
    return [lcm([x.denominator() for x in f.coefficients()])*f for f in tmp]

def find_divisor(F,x):
    R=F[0].parent()
    gens=R.gens()
    y=gens[(gens.index(x)+1)%len(gens)]
    F1=[f.subs(dict([(x,0),(y,1)])) for f in F]
    S.<y>=PolynomialRing(RationalField())
    others=[]
    for f in F1:
        if list(f.degrees()).count(0)==len(gens)-1:
            # It means that it is really a single variable polynomial
            ii=list(f.degrees()).index(f.degree())
            xi=gens[ii]
            lst=[]
            for jj in range(len(gens)):
                if jj==ii:
                    lst.append(S.gen(0))
                else:
                    lst.append(0)
            phi=R.hom(lst,codomain=S,check=False)
            fone=phi(f)
            S0=S.base_extend((fone/fone.leading_coefficient()).root_field('a'))
            a=S0(fone).roots()[0][0]
        else:
            others.append(f)
    others=[f.subs(dict([(f.parent().gen(ii),a)])) for f in others]
    return others
