r"""
Graded Ring of Minkowski Weights

EXAMPLE::
    sage: A = MW(toric_varieties.P2().fan()); A
    Graded ring of Minkowski weights with top degree 2
    


AUTHORS:
- Elise Villella (2019) Initial version
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


class MW:
    def __init__(self,fan):
        """Defines graded ring of Minkowski weights associated to given fan. Computes ranks of each graded piece, 
        computes basis for each space of weights, can multiply weights.
        
        ::
            
            sage: A = MW(toric_varieties.P2().fan()); A
            Graded ring of Minkowski weights with top degree 2
            
        """
        self.fan = fan
        self.cones = fan.cones
        self.dim = fan.dim()
        self.rk, self.A = self._set_ranks_matrices() ## NOTE: ranks and A's are for COdimension
        self.A_RR = [M.echelon_form() for M in self.A]
        self.bases = [M.T.kernel().basis() for M in self.A]
        self.N = fan.lattice()
        
    def basis(self,k):
        """Returns a basis (list of basis elements) for MW^k, space of weights of codim k"""
        return self.bases[k]
    
    def ranks(self):
        """Returns list of ranks of each graded pieces"""
        return self.rk
    
    def product(self,w1,cd1,w2,cd2):
        #takes in two weights, w1 and w2 of codimensions cd1 and cd2 resp
        #returns their product w of weight cd1+cd2
        #product evaluated on cone of correct dimension is:
        #sum_{sigma,tau of correct dimensions} m_{sigma,tau} c(sigma) \tilde c(tau)
        #where m_{sigma,tau} is [N,N_sigma+N_tau] as long as sigma meets tau+v (generic fixed v)
        #and gamma is contained in both sigma and tau
        cd = cd1+cd2
        d = self.dim - cd
        d1 = self.dim - cd1
        d2 = self.dim - cd2
        assert self._check_weight(w1,cd1),"Weight w1 is not a valid weight of codimension cd1"
        assert self._check_weight(w2,cd2),"Weight w2 is not a valid weight of codimension cd2"
        assert (len(w1)==len(self.fan(d1))),"Length of weight 1 does not match length of basis"##not needed anymore
        assert (len(w2)==len(self.fan(d2))),"Length of weight 2 does not match length of basis"##replaced by new check weight function
        if d < 0:
            return []
        else: #will have to compute product and evaluate on all cones of dim d
            soln = []
            cones = self.fan(d)
            Sigmas = self.fan(d1)
            Taus = self.fan(d2)

            #determine generic v
            v = self._generic(self.fan)

            new_weight = []
            for c in cones:
                # only need to compute m_{sig,tau} for pairs which both contain c as a face
                sig_c = self._cones_containing(c,Sigmas) #cones in Sigmas which contain c as a face
                tau_c = self._cones_containing(c,Taus)

                relevant = self._check_generic(sig_c,tau_c,v)

                c_sum = 0

                #for each ordered pair in relevant, compute coefficient m_{sig,tau} = [N:N_sig + N_tau]
                for (sig,tau) in relevant:
                    N_sig = sig.sublattice()
                    N_tau = tau.sublattice()
                    N_sum = self.N.submodule(N_sig.basis()+N_tau.basis()) # computes N_sig + N_tau

                    #N_sum.index_in(N) is m_{sig,tau}
                    # for c(sig), first find location of sig in list of cones of dimension d1, then take corresp. entry of w1
                    c_sum += N_sum.index_in(self.N)*w1[Sigmas.index(sig)]*w2[Taus.index(tau)]
                new_weight.append(c_sum)

            assert self._check_weight(new_weight,cd),"New weight fails to be balanced or the right size for codimension cd"
            return new_weight
        
    def __repr__(self):
        r"""
        Return a string representation of ``self``.
        EXAMPLES::
            sage: MW(toric_varieties.P2().fan())._repr_()
            "Graded ring of Minkowski weights with top degree 2"
        """
        return ('Graded ring of Minkowski weights with top degree %s' % (self.dim))

    def _balancing(self, tau):
        """returns list of relations for cone tau (coefficient vectors) for cones of dimension (dim tau + 1)
        relation is: sum_{sigma > tau, codim 1} < u, n_{sigma,tau} > c(sigma) = 0, where 
        n_{sigma,tau} is lattice pt in sigma which generates N_sigma/N_tau (1D v.space)
        get one relation for each basis element u for tau^perp"""

        relns = []

        d = tau.dim()
        l = len(self.cones(d+1)) # dimension of vector relations to be returned
        #print('number of cones of dimension '+str(d+1)+ ' is: ' + str(l))

        # tuple of cones of dimension d+1 with tau as facet, i.e., relevant cones
        relevant = tau.facet_of()

        # choose basis for space perpendicular to tau, in list
        basis = tau.orthogonal_sublattice().basis()

        for b in basis: #each b gives relation, i.e., vector in R^l
            v = [] # vector holding coefficients of relation for given b, at end of loop append to relns and reset
            for c in self.cones(d+1):# each c gives one component of this vector
                if c in tau.facet_of():
                    Q = c.relative_quotient(tau)
                    n = Q.gens()[0] #only one generator?
                    v.append(b*n)
                else:
                    v.append(0)

            relns.append(v)
        return relns

    def _set_rank(self,cd):
        # returns the rank of MW^cd, the space of functions on cones of codim cd, matrix A whose null space is space of weights
        d = self.dim - cd
        #print('dimension is: '+str(d))
        ConeList = self.cones(d) #it's a tuple!
        n = len(ConeList)
        #print(str(n)+' cones of that dimension')

        #generate balancing conditions for every cone of one dimension smaller than d
        ConeRelns = self.cones(d-1)
        relations = []
        for c in ConeRelns:
            balance = self._balancing(c)
            for b in balance:
                relations.append(b)
        #print(relations)

        #interpret relations, i.e., may be redundant, and determine rank
        #that is, dimension of space perpendicular to the list of vectors
        A = Matrix(relations)
        r = rank(A)
        # kernel(A.T) gives a basis for the space of MW^cd where coordinates designate coefficients of cones in that place
        return n-r, A

    def _set_ranks_matrices(self):
        # computes ranks, matrices whose null spaces give minkowski weights
        rnks = []
        mtrx = []
        for cd in range(self.dim):
            r,M = self._set_rank(cd)
            rnks.append(r)
            mtrx.append(M)
        rnks.append(1)
        mtrx.append(Matrix(1))
        return rnks, mtrx

    def _check_weight(self,w,cd):
        """Takes weight w and codimension cd, returns True if w can be evaluated on cones of codim cd AND balancing
        conditions are satisfied on the fan of self"""
        d = self.dim-cd
        res = self.A[cd]*Matrix(w).T ##res is a vector of length = num relations, should be zero vector
        if len(w)!=len(self.fan(d)):
            print("incompatible list lengths")
        if res.is_zero == False:
            print("relations not satisfied")
        return (len(w)==len(self.fan(d))) and res.is_zero()
    
    def _generic(self,fan):
        """ returns a vector v which is generic with respect to the given fan """
        d = fan.dim()

        #random candidate for genric vector
        v = [random() for r in range(d)]

        #check whether generic, i.e., v should NOT be contained in any cones of codimension 1
        needToCheck = True

        while needToCheck:
            coneFlag = False #change to true if v in a cone of codim 1
            for c in fan.cones(d-1):
                if v in c:
                    coneFlag = True
            if coneFlag:
                #generate new random vector and try again
                v = [random() for r in range(d)]
            else:
                #vector v is generic
                needToCheck = False
        return v

    def _check_generic(self,cones1,cones2,v):
        """ for every pair of cones c1, c2 from respective lists, check whether they intersect generically wrt v
        return list of ordered pairs of cones which intersect wrt v """
        good_pairs = []
        for c1 in cones1:
            for c2 in cones2:
                # check whether v is in c1 - c2 (i.e. mink sum of c1 and -c2)
                C = Cone(rays = [r for r in c1.rays()]+[-1*r for r in c2.rays()])
                if C.contains(v):
                    good_pairs.append((c1,c2))
        return good_pairs

    def _cones_containing(self,cone,conelist):
        """ returns sublist of conelist whose cones have cone as a face """
        toreturn = []
        for c in conelist:
            if cone.is_face_of(c):
                toreturn.append(c)
        return toreturn
