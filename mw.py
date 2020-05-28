r"""
Graded Ring of Minkowski Weights

A Minkowski weight is a function on cones of the fan subject
to balancing conditions. A weight on cones of codimension $c$
is in graded degree $c$. With the order of cones fixed by 
the fan, we represent a weight as a list of values.

In [FS1994]_ it is shown that this ring is isomorphic to
the (operational) Chow ring of the associated toric variety
so this class provides a way to compute in this ring. For this
reason, the ranks and bases of graded pieces are of interest,
as well as the implementation of a product of weights. This
computation is done via fan displacement.

EXAMPLES:

The simplest example is P2 where this matches cohomology::

    sage: A = MW(toric_varieties.P2().fan()); A
    Graded ring of Minkowski weights with top degree 2
    
The ranks are as follows, ordered by codimension::

    sage: A.ranks()
    [1, 1, 1]
    
For a given codimension, we get a basis::

    sage: A.basis(1)
    [
    (1, 1, 1)
    ]
    
Finally, we can multiply weights to obtain new weights::

    sage: w1 = [2,2,2]
    sage: w2 = [3,3,3]
    sage: A.product(w1,1,w2,1)
    [6]


AUTHORS:
- Elise Villella (2019-07-25) Initial version
"""

#*****************************************************************************
#       Copyright (C) 2019      Elise Villella <elisevillella@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from __future__ import print_function
from sage.matrix.constructor import Matrix
from sage.misc.functional import rank
from sage.misc.prandom import random
from sage.geometry.cone import Cone
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors


class MW:
    def __init__(self,fan):
        r"""
        The graded ring of Minkowski weights on fan.
        
        A Minkowski weight is a function on cones of the fan subject
        to balancing conditions. A weight on cones of codimension `c`
        is in graded degree `c`. With the order of cones fixed by 
        the fan, we represent a weight as a list of values.
        
        In [FS1994]_ it is shown that this ring is isomorphic to
        the (operational) Chow ring of the associated toric variety
        so this class provides a way to compute in this ring. For this
        reason, the ranks and bases of graded pieces are of interest,
        as well as the implementation of a product of weights. This
        computation is done via fan displacement.
        
        INPUT:
        
        - ``fan`` -- A Fan (e.g. normal fan of a polytope)
        
        OUTPUT:
        
        - MW object
        
        EXAMPLES::
            
            sage: A = MW(toric_varieties.P2().fan()); A
            Graded ring of Minkowski weights with top degree 2
       
            sage: A.ranks()
            [1, 1, 1]
    
            sage: A.basis(1)
            [
            (1, 1, 1)
            ]
    
            sage: w1 = [2,2,2]
            sage: w2 = [3,3,3]
            sage: A.product(w1,1,w2,1)
            [6]
            
        
            
        REFERENCES:
        
        - [FS1994]_
        """
        self._fan = fan
        self._cones = fan.cones
        self._dim = fan.dim()
        self._rk, self._A = self._set_ranks_matrices() ## NOTE: ranks and A's are for COdimension
        self._A_RR = [M.echelon_form() for M in self._A]
        self._bases = [M.T.kernel().basis() for M in self._A]
        self._N = fan.lattice()
        
    def basis(self,k):
        """Returns a basis for MW^k.
        
        Computes a basis for the space of weights of codim `k`, so
        each basis element will be a balanced weight on cones of 
        dimension `k`. This is the easiest way to generate balanced
        weights on fans.
        
        INPUT:
        
        - ``k`` -- integer; This integer must be between 0 and the 
          dimension of the fan.
          
        OUTPUT: A list of basis elements for degree (codimension) `k`.
        
        EXAMPLES::
            
            sage: A = MW(toric_varieties.Cube_face_fan().fan())
            sage: A.basis(2)
            [
            (1, 0, 0, 0, 0, 1, 1, -1),
            (0, 1, 0, 0, 0, 1, 0, 0),
            (0, 0, 1, 0, 0, 0, 1, 0),
            (0, 0, 0, 1, 0, 0, 0, 1),
            (0, 0, 0, 0, 1, -1, -1, 1)
            ]
        """
        return self._bases[k]
    
    def ranks(self):
        """Returns list of ranks of MW graded by codimension.
        
        OUTPUT: List of integer ranks graded by codimension.
        
        EXAMPLES::
        
            sage: A = MW(toric_varieties.Cube_face_fan().fan())
            sage: A.ranks()
            [1, 1, 5, 1]
        """
        return self._rk
    
    def product(self,w1,cd1,w2,cd2):
        r"""
        Computes product of two weights.
        
        Takes weights w1 and w2 of codimensions cd1 and cd2 respectively.
        Returns their product; a balanced weight of codimension cd1+cd2
        Product formula is: $
        `w1*w2(\gamma)= \sum_{\sigma,\tau} m_{\sigma,\tau} w1(\sigma) w2(\tau)`
        where `m_{\sigma,\tau} = [N,N_\sigma + N_\tau]` as long as `\sigma` 
        meets `\tau + v` for generic fixed `v` and `\gamma \subset \sigma \cap \tau`. 
        See [FS1994]_.
        
        INPUT:
        
        - ``w1`` -- list or vector. This must be a balanced weight on the fan.
        - ``cd1`` -- integer. This must be the degree of `w1`.
        
        - ``w2`` -- list or vector. This must be a balanced weight on the fan.
        - ``cd2`` -- integer. This must be the degree of `w2`.
        
        OUTPUT: Balanced weight of degree (codimension) cd1+cd2. This is an empty list if cd1+cd2 > dim fan.
        
        EXAMPLES::
        
            sage: A = MW(toric_varieties.P2().fan())
            sage: w1 = [2,2,2]
            sage: w2 = [3,3,3]
            sage: A.product(w1,1,w2,1)
            [6]
            
        TESTS::
        
            sage: A = MW(toric_varieties.P2().fan())
            sage: w3 = [2]
            sage: w4 = [5]
            sage: A.product(w3,2,w4,2)
            []
            sage: w5 = [2,3]

        """
        cd = cd1+cd2
        d = self._dim - cd
        d1 = self._dim - cd1
        d2 = self._dim - cd2
        assert self._check_weight(w1,cd1),"Weight w1 is not a valid weight of codimension cd1"
        assert self._check_weight(w2,cd2),"Weight w2 is not a valid weight of codimension cd2"
        assert (len(w1)==len(self._fan(d1))),"Length of weight 1 does not match length of basis"##not needed anymore
        assert (len(w2)==len(self._fan(d2))),"Length of weight 2 does not match length of basis"##replaced by new check weight function
        if d < 0:
            return []
        else: #compute product and evaluate on all cones of dim d
            soln = []
            cones = self._fan(d)
            Sigmas = self._fan(d1)
            Taus = self._fan(d2)

            #get a generic v
            v = self._generic(self._fan)

            new_weight = []
            for c in cones:
                # only need to compute m_{sig,tau} for pairs which both contain c as a face
                sig_c = self._cones_containing(c,Sigmas)
                tau_c = self._cones_containing(c,Taus)

                relevant = self._check_generic(sig_c,tau_c,v)

                c_sum = 0

                #for each relevant ordered pair, compute coefficient m_{sig,tau} = [N:N_sig + N_tau]
                for (sig,tau) in relevant:
                    N_sig = sig.sublattice()
                    N_tau = tau.sublattice()
                    N_sum = self._N.submodule(N_sig.basis()+N_tau.basis()) # computes N_sig + N_tau

                    #N_sum.index_in(N) is m_{sig,tau}
                    # for c(sig), first find sig in list of dim d1 cones, then take that entry of w1
                    c_sum += N_sum.index_in(self._N)*w1[Sigmas.index(sig)]*w2[Taus.index(tau)]
                new_weight.append(c_sum)

            assert self._check_weight(new_weight,cd),"New weight fails to be balanced or the right size for codimension cd"
            return new_weight                    
        
    def __repr__(self):
        r"""
        Return a string representation of ``self``.
        
        OUTPUT: String describing graded ring MW.
        
        EXAMPLES::
        
            sage: MW(toric_varieties.P2().fan()).__repr__()
            'Graded ring of Minkowski weights with top degree 2'
            
        """
        return ('Graded ring of Minkowski weights with top degree %s' % (self._dim))

    def _balancing(self, tau):
        r"""
        Returns relations associated with cone `\tau`.
        
        Returns list of relations for cone $\tau$ for cones of dimension (dim tau + 1) containing $\tau$
        Relation is: $\sum_{\sigma > \tau, \dim \sigma - \dim \tau = 1} \langle u, n_{\sigma,\tau} \rangle c(\sigma) = 0$ 
        $n_{\sigma,\tau}$ is a lattice point in $\sigma$ which generates $N_\sigma/N_\tau$ (1D v.space) modulo $\tau$
        Returned list contains one relation for each basis element u of $\tau^\perp$
        
        TESTS::
        
            sage: vertices = [(0,1,0),(0,1,1),(0,2,0),(0,2,2),(1,1,1),(1,2,1),(1,2,2)]
            sage: fan = NormalFan(Polyhedron(vertices = vertices))
            sage: A = MW(fan)
            sage: r1 = fan(1)[2]
            sage: A._balancing(r1)
            [[0, 1, -1, 0, 0, 0, 0, 1, 0, 0, 0], [0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0]]
            sage: r2 = fan(0)[0]
            sage: A._balancing(r2)
            [[0, 0, 1, -1, 0, -1], [1, -1, 0, 0, 1, 0], [-1, 0, 0, 1, 0, 0]]
            sage: r3 = fan(3)[0]
            sage: A._balancing(r3)
            []
            
        """

        relns = []

        d = tau.dim()
        l = len(self._cones(d+1)) # dimension of vector relations to be returned

        # tuple of cones of dimension d+1 with tau as facet, i.e., relevant cones
        relevant = tau.facet_of()

        # choose basis for space perpendicular to tau, in list
        basis = tau.orthogonal_sublattice().basis()

        for b in basis: #each b gives relation, i.e., vector in R^l
            v = [] # vector holding coefficients of relation for given b
            for c in self._cones(d+1):# each c gives one component of this vector
                if c in tau.facet_of():
                    Q = c.relative_quotient(tau)
                    n = Q.gens()[0] #only one generator?
                    v.append(b*n)
                else:
                    v.append(0)

            relns.append(v)
        return relns

    def _set_rank(self,cd):
        r"""
        Returns the rank of degree `cd` as well as the matrix of relations.
        
        For given codimension cd, returns dimension of kernel of relations which is rank of MW^cd and matrix A of relations
        Any vector in the kernel of A is a balanced weight of codimension cd
        
        TESTS::
        
            sage: vertices = [(0,1,0),(0,1,1),(0,2,0),(0,2,2),(1,1,1),(1,2,1),(1,2,2)]
            sage: fan = NormalFan(Polyhedron(vertices = vertices))
            sage: A = MW(fan)
            sage: A._set_rank(2)
            (
               [ 0  0  1 -1  0 -1]
               [ 1 -1  0  0  1  0]
            3, [-1  0  0  1  0  0]
            )
            sage: A._set_rank(5)
            (0, [])

        
        """
        d = self._dim - cd
        ConeList = self._cones(d) #it's a tuple!
        n = len(ConeList)

        #generate balancing conditions for every cone of dimension d-1
        ConeRelns = self._cones(d-1)
        relations = []
        for c in ConeRelns:
            balance = self._balancing(c)
            for b in balance:
                relations.append(b)

        # remove redundant relations, rank = dimension of space perpendicular to the list of vectors
        A = Matrix(relations)
        r = rank(A)
        # kernel(A.T) gives a basis for MW^cd where coordinates designate coefficients of cones in that place
        return n-r, A

    def _set_ranks_matrices(self):
        r"""
        Iterates through all codimensions and initializes lists of ranks and relation matrices
        
        TESTS::
        
            sage: A = MW(toric_varieties.P2().fan())
            sage: ascii_art(sorted(A._set_ranks_matrices()))
            [ [ [ 1  0 -1]                  ]              ]
            [ [ [ 1 -1  0]  [ 1  0 -1]      ]              ]
            [ [ [ 0 -1  1], [ 0  1 -1], [0] ], [ 1, 1, 1 ] ]
            
        """
        rnks = []
        mtrx = []
        for cd in range(self._dim):
            r,M = self._set_rank(cd)
            rnks.append(r)
            mtrx.append(M)
        rnks.append(1)
        mtrx.append(Matrix(1))
        return rnks, mtrx

    def _check_weight(self,w,cd):
        r"""
        Determines whether `w` is a balanced weight of codimension `cd`.
        
        Takes weight `w` and codimension `cd`, returns True if `w` can be 
        evaluated on cones of codim `cd` AND balancing conditions are 
        satisfied on the fan of self.
        
        TESTS::
        
            sage: A = MW(toric_varieties.P2().fan())
            sage: w1 = [5,5,5]
            sage: w2 = [2,3,2]
            sage: A._check_weight(w1,1)
            True
            sage: A._check_weight(w2,1)
            False
            sage: A._check_weight(w1,2)
            False
            sage: A._check_weight(w2,13)
            False
            
        """
        d = self._dim-cd
        if len(w)!=len(self._fan(d)): #list lengths incompatible; can't multiply
            return False        
        res = self._A[cd]*Matrix(w).T ##res is a vector of length = num relations, should be zero vector
        return res.is_zero()
    
    def _generic(self,fan):
        r""" Returns a vector `v` which is generic with respect to the given fan .
        
        TESTS::
        
            sage: A = MW(toric_varieties.P2().fan())
            sage: A._generic(A._fan()) # random output
            [0.2573820339704006, 0.24526498452629364]
        
        """
        d = fan.dim()

        v = [random() for r in range(d)] #candidate for generic

        needToCheck = True

        while needToCheck:
            coneFlag = False #change to true if v in a cone of codim 1, that makes v NOT generic
            for c in fan.cones(d-1):
                if v in c:
                    coneFlag = True
            if coneFlag: #try new generic vector
                v = [random() for r in range(d)]
            else: #v is generic
                needToCheck = False
        return v

    def _check_generic(self,cones1,cones2,v):
        """ 
        Returns pairs of cones which intersect with respect to `v`.
        
        For every pair of cones `c1`, `c2` from respective lists, check whether 
        they intersect generically with respect to `v`.
        Return list of ordered pairs of cones which intersect with respect to `v`.
        This means the intersection of `c1` and `c2+v` is non-empty.
        
        EXAMPLES::
            
            sage: A = MW(toric_varieties.P2().fan())
            sage: cones = A._fan(1)
            sage: v = (.1,.3)
            sage: A._check_generic(cones,cones,v)
            [(1-d cone of Rational polyhedral fan in 2-d lattice N, 1-d cone of Rational polyhedral fan in 2-d lattice N)]
        
        """
        good_pairs = []
        for c1 in cones1:
            for c2 in cones2:
                # check whether v is in c1 - c2 (i.e. mink sum of c1 and -c2)
                C = Cone(rays = [r for r in c1.rays()]+[-1*r for r in c2.rays()])
                if C.contains(v):
                    good_pairs.append((c1,c2))
        return good_pairs

    def _cones_containing(self,cone,conelist):
        r""" Returns sublist of conelist whose cones have cone as a face. 
        
        EXAMPLES::
        
            sage: A = MW(toric_varieties.P(3).fan())
            sage: conelist = A._fan(3)
            sage: cone = A._fan(1)[0]
            sage: A._cones_containing(cone,conelist)
            [3-d cone of Rational polyhedral fan in 3-d lattice N, 3-d cone of Rational polyhedral fan in 3-d lattice N, 3-d cone of Rational polyhedral fan in 3-d lattice N]
            
        """
        toreturn = []
        for c in conelist:
            if cone.is_face_of(c):
                toreturn.append(c)
        return toreturn