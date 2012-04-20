######################################################################
##
##  Code to create the Manin Relations class, which solves the "Manin
##  relations".  That is, a description of Div^0(P^1(Q)) as a
##  Z[Gamma_0(N)]-module in terms of generators and relations is
##  found.  The method used is geometric, constructing a nice
##  fundamental domain for Gamma_0(N) and reading the relevant Manin
##  relations off of that picture.  The algorithm follows the paper of
##  Pollack and Stevens "Overconvergent modular symbols and p-adic
##  L-functions"
##
##  Copyright (c) 2012, Rob Pollack and Jonathan Hanke
##      <rpollack@math.bu.edu>
##      <jonhanke@gmail.com>
##
##  Released under the GNU Public License, 2012.
##
######################################################################

from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
from sage.modular.modsym.all import P1List
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from sage.modules.free_module_element import zero_vector
from copy import deepcopy

M2Z = MatrixSpace_ZZ_2x2()
t00 = (0,0)
t10 = (1,0)
t01 = (0,1)
t11 = (1,1)

class PSModularSymbolsDomain(SageObject):
    def __init__(self, N, reps, indices, equiv_ind):
        self._N = N
        self._reps = reps
        self._indices = sorted(indices)
        self._gens = [reps[i] for i in self._indices]
        self._equiv_ind = equiv_ind
        self._equiv_rep = {}
        for ky in equiv_ind:
            self._equiv_rep[ky] = reps[equiv_ind[ky]]

    def __len__(self):
        return len(self._reps)

    def __getitem__(self, i):
        return self._reps[i]

    def __iter__(self):
        return iter(self._reps)

    def gens(self):
        return self._gens

    def gen(self, n=0):
        return self._gens[n]

    def ngens(self):
        return len(self._gens)

    def level(self):
        r"""
        Returns the level `N` of `\Gamma_0(N)` that we work with.

        OUTPUT:

        - The integer `N` of the group `\Gamma_0(N)` for which the
          Manin Relations are being computed.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.level()
            11
        """
        return self._N

    def indices(self):
        return self._indices

    def reps(self, n=None):
        r"""
        Returns the n-th coset rep associated with our fundamental
        domain or all coset reps if n is not specified.

        INPUT:

        - ``n`` -- integer (default: None)

        OUTPUT:

        - If n is given then the n-th coset representative is returned
          and otherwise all coset reps are returned.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.reps(0)
            [1 0]
            [0 1]
            sage: A.reps(1)
            [ 1  1]
            [-1  0]
            sage: A.reps(2)
            [ 0 -1]
            [ 1  3]
            sage: A.reps()
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]  [ 0 -1]  [ 1  0]
            [0 1], [-1  0], [ 1  3], [ 3  2], [ 2  3], [ 3  1], [ 1  2], [-2  1],
            <BLANKLINE>
            [ 0 -1]  [ 1  0]  [-1 -1]  [ 1 -1]
            [ 1  1], [-1  1], [ 2  1], [-1  2]
            ]
        """
        if n is None:
            return self._reps
        else:
            return self._reps[n]

######################################
##  Define the Manin Relation Class ##
######################################

class ManinRelations(SageObject):
    """
    This class gives a description of Div^0(P^1(QQ)) as a
    `\ZZ[\Gamma_0(N)]`-module.

    INPUT:

    - ``N`` -- a positive integer
    """
    def __init__(self, N):
        r"""
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
        """
        ## Store the level
        self._N = N

        ## Creates and stores the Sage representation of P^1(Z/NZ)
        P = P1List(N)
        self._P = P

        ## Creates a fundamental domain for Gamma_0(N) whose boundary is a union
        ## of unimodular paths (except in the case of 3-torsion).
        ## We will call the intersection of this domain with the real axis the
        ## collection of cusps (even if some are Gamma_0(N) equivalent to one another).
        cusps = self.form_list_of_cusps()

        ## Takes the boundary of this fundamental domain and finds SL_2(Z) matrices whose
        ## associated unimodular path gives this boundary.  These matrices form the
        ## beginning of our collection of coset reps for Gamma_0(N) / SL_2(Z).
        coset_reps = self.fd_boundary(cusps)

        ## Make the identity matrix
        Id = M2Z([1,0,0,1])

        ## Gives names to matrices of order 2 and 3 (in PSL_2)
        sig = M2Z([0,1,-1,0])
        tau = M2Z([0,-1,1,-1])

        ## Takes the bottom row of each of our current coset reps,
        ## thinking of them as distinct elements of P^1(Z/NZ)
        p1s = [(coset_reps[j])[1] for j in range(len(coset_reps))]

        ## Initializes relevant Manin data
        gens_index = []
        twotor_index = []
        twotorrels = []
        threetor_index = []
        threetorrels = []
        rels = [0 for i in range(0,len(coset_reps))]

        ## the list rels (above) will give Z[Gamma_0(N)] relations between
        ## the associated divisor of each coset representatives in terms
        ## of our chosen set of generators.
        ## entries of rel will be lists of elements of the form (c,A,r)
        ## with c a constant, A a Gamma_0(N) matrix, and r the index of a
        ## generator.  The meaning is that the divisor associated to the
        ## j-th coset rep will equal the sum of:
        ##
        ##   c * A^(-1) * (divisor associated to r-th coset rep)
        ##
        ## as one varies over all (c,A,r) in rel[j].
        ## (Here r must be in self.generator_indices().)
        ##
        ## This will be used for modular symbols as then the value of a
        ## modular symbol phi on the (associated divisor) of the j-th
        ## element of coset_reps will be the sum of c * phi (r-th genetator) | A
        ## as one varies over the tuples in rel[j]

        boundary_checked = [False for i in range(0,len(coset_reps))]

        ## The list boundary_checked keeps track of which boundary pieces of the
        ## fundamental domain have been already used as we are picking
        ## our generators

#        glue_data = [0 for i in range(0,len(coset_reps))]   ## ????


        ## The following loop will choose our generators by picking one edge
        ## out of each pair of edges that are glued to each other and picking
        ## each edge glued to itself (arising from two-torsion)
        ## ------------------------------------------------------------------
        for r in range(0,len(coset_reps)):
            if boundary_checked[r] == False:

                ## We now check if this boundary edge is glued to itself by
                ## Gamma_0(N)

                if P.index(p1s[r][0],p1s[r][1]) == P.index(-p1s[r][1],p1s[r][0]):
                    ## This edge is glued to itself and so coset_reps[r]
                    ## needs to be added to our generator list.
                    rels[r] = [(1,Id,r)] ## this relation expresses the fact
                                              ## that coset_reps[r] is one of our
                                                   ## basic generators

                    gens_index = gens_index + [r]    ## the index r is adding to our list
                                                  ## of indexes of generators
                    twotor_index = twotor_index + [r]  ## the index r is adding to our
                                                         ## list of indexes of generators
                                          ## which satisfy a 2-torsion relation

                    gam = coset_reps[r] * sig * coset_reps[r]**(-1)
                    ## gam is 2-torsion matrix and in Gamma_0(N).
                    ## if D is the divisor associated to coset_reps[r]
                    ## then gam * D = - D and so (1+gam)D=0.

                    ## This gives a restriction to the possible values of
                    ## modular symbols on D

                    twotorrels = twotorrels + [gam]

                    ## The 2-torsion matrix gam is recorded in our list of
                    ## 2-torsion relations.

#                    glue_data[r]=(r,gam)

                    ## this records that the edge associated to coset_reps[r]
                    ##is glued to itself by gam

                    boundary_checked[r] = True

                    ## We have now finished with this edge.

                else:
                    c = coset_reps[r][1,0]
                    d = coset_reps[r][1,1]
                    if (c**2+d**2+c*d)%N == 0:
                    ## If true, then the ideal triangle below the unimodular
                    ## path described by coset_reps[r] contains a point
                    ## fixed by a 3-torsion element.

                        gens_index = gens_index + [r]
                        ## the index r is adding to our list of indexes
                        ## of generators

                        rels[r] = [(1,Id,r)]
                        ## this relation expresses the fact that coset_reps[r]
                        ## is one of our basic generators

                        threetor_index = threetor_index + [r]
                        ## the index r is adding to our list of indexes of
                        ##generators which satisfy a 3-torsion relation

                        gam = coset_reps[r] * tau * coset_reps[r]**(-1)
                        ## gam is 3-torsion matrix and in Gamma_0(N).
                        ## if D is the divisor associated to coset_reps[r]
                        ## then (1+gam+gam^2)D=0.
                        ## This gives a restriction to the possible values of
                        ## modular symbols on D

                        threetorrels = threetorrels + [gam]

                        ## The 3-torsion matrix gam is recorded in our list of
                        ## 2-torsion relations.
                        ###  DID I EVER ADD THE GLUE DATA HERE?  DO I NEED IT?

                        ## The reverse of the unimodular path associated to
                        ## coset_reps[r] is not Gamma_0(N) equivalent to it, so
                        ## we need to include it in our list of coset
                        ## representatives and record the relevant relations.

                        a = coset_reps[r][0][0]
                        b = coset_reps[r][0][1]
                        A = M2Z([-b,a,-d,c])
                        A.set_immutable()
                        coset_reps.append(A)
                        ## A (representing the reversed edge) is included in
                        ## our list of coset reps

                        rels = rels + [[(-1,Id,r)]]
                        ## This relation means that phi on the reversed edge
                        ## equals -phi on original edge

                        boundary_checked[r] = True
                        ## We have now finished with this edge.

                    else:

                        ## This is the generic case where neither 2 or
                        ## 3-torsion intervenes.
                        ## The below loop searches through the remaining edges
                        ## and finds which one is equivalent to the reverse of
                        ##  coset_reps[r]
                        ## ---------------------------------------------------
                        found = False

                        s = r + 1
                        ## s is the index we use to run through the
                        ## remaining edges

                        while (not found):
                            if P.index(p1s[s][0],p1s[s][1]) == P.index(-p1s[r][1],p1s[r][0]):
                                ## the reverse of coset_reps[r] is
                                ## Gamma_0(N)-equivalent to coset_reps[s]
                                ## coset_reps[r] will now be made a generator
                                ## and we need to express phi(coset_reps[s])
                                ## in terms of phi(coset_reps[r])

                                gens_index = gens_index + [r]
                                ## the index r is adding to our list of
                                ## indexes of generators

                                rels[r] = [(1,Id,r)]
                                ## this relation expresses the fact that
                                ## coset_reps[r] is one of our basic generators

                                A = coset_reps[s] * sig
                                ## A corresponds to reversing the orientation
                                ## of the edge corr. to coset_reps[r]

                                gam = coset_reps[r] * A**(-1)
                                ## gam is in Gamma_0(N) (by assumption of
                                ## ending up here in this if statement)

                                rels[s] = [(-1,gam,r)]
                                ## this relation means that phi evaluated on
                                ## coset_reps[s] equals -phi(coset_reps[r])|gam
                                ## To see this, let D_r be the divisor
                                ## associated to coset_reps[r] and D_s to
                                ## coset_reps[s]. Then gam D_s = -D_r and so
                                ## phi(gam D_s) = - phi(D_r) and thus
                                ## phi(D_s) = -phi(D_r)|gam
                                ## since gam is in Gamma_0(N)

#                                glue_data[r] = (s,gam)  ## ????
                                found = True
                                boundary_checked[r] = 1
                                boundary_checked[s] = 1

                            else:
                                s = s + 1  ## moves on to the next edge

        ## We now need to complete our list of coset representatives by
        ## finding all unimodular paths in the interior of the fundamental
        ## domain, as well as express these paths in terms of our chosen set
        ## of generators.
        ## -------------------------------------------------------------------

        for r in range(len(cusps)-2):

        ## r is the index of the cusp on the left of the path.  We only run
        ## thru to the number of cusps - 2 since you can't start an interior
        ## path on either of the last two cusps

            for s in range(r+2,len(cusps)):
            ## s is in the index of the cusp on the the right of the path
                cusp1 = cusps[r]
                cusp2 = cusps[s]
                if self.is_unimodular_path(cusp1,cusp2):
                    A,B = self.unimod_to_matrices(cusp1,cusp2)
                    ## A and B are the matrices whose associated paths
                    ## connect cusp1 to cusp2 and cusp2 to cusp1 (respectively)
                    coset_reps = coset_reps + [A,B]
                    ## A and B are added to our coset reps
                    vA = []
                    vB = []

                    ## This loop now encodes the relation between the
                    ## unimodular path A and our generators.  This is done
                    ## simply by accounting for all of the edges that lie
                    ## below the path attached to A (as they form a triangle)
                    ## Similarly, this is also done for B.

                    for j in range(r+2,s+2):
                    ## Running between the cusps between cusp1 and cusp2
                        vA = vA + rels[j]  ## Edge relation added
                        t = (-rels[j][0][0],rels[j][0][1],rels[j][0][2]) ## This is simply the negative of the above edge relation.
                        vB = vB + [t]  ## Negative of edge relation added
                    rels = rels + [vA,vB]  ## Relations for A and B adding to relations list
####        return [coset_reps,gens_index,twotor_index,twotorrels,threetor_index,threetorrels,rels,glue_data]


        ## Store the data coming from solving the Manin Relations
        ## ======================================================

        self._mats = coset_reps

        ## Coset representatives of Gamma_0(N) coming from the geometric
        ## fundamental domain algorithm

        ## Make the translation table between the Sage and Geometric
        ## descriptions of P^1

        equiv_ind = {}
        for i, rep in enumerate(coset_reps):
            ky = P.normalize(rep[t10],rep[t11])
            equiv_ind[ky] = i

        self._gens_index = gens_index
        ## This is a list of indices of the (geometric) coset representatives
        ## whose values (on the associated degree zero divisors) determine the
        ## modular symbol.

        self._ngens = len(self._gens_index)
        self._gens = [self._mats[i] for i in self._gens_index]

        self._twotor_index = twotor_index
        ## A list of indices of the (geometric) coset representatives whose
        ## paths are identified by some 2-torsion element (which switches the
        ## path orientation)
        self._twotor = [self._mats[i] for i in self._twotor_index]

        self._twotorrels = twotorrels
        ## A list of (2-torsion in PSL_2(Z)) matrices in Gamma_0(N) that give
        ## the orientation identification in the paths listed in twotor_index above!

        self._threetor_index = threetor_index
        ## A list of indices of the (geometric) coset representatives that
        ## form one side of an ideal triangle with an interior fixed point of
        ## a 3-torsion element of Gamma_0(N)
        self._threetor = [self._mats[i] for i in self._threetor_index]

        self._threetorrels = threetorrels
        ## A list of (3-torsion in PSL_2(Z)) matrices in Gamma_0(N) that give
        ## the interior fixed point described in threetor_index above!

        self._rels = rels
        self._rel_dict = {}
        for j, L in enumerate(rels):
            self._rel_dict[self._mats[j]] = [(d, A, self._mats[i]) for (d, A, i) in L]
        ## A list of lists of triples (d, A, i), one for each coset
        ## representative of Gamma_0(N) (ordered to correspond to the
        ## representatives of self.coset_reps) expressing the value of a
        ## modular symbol on the associated unimodular path as a sum of terms
        ##    d * (value on the i-th coset rep) | A
        ## where the index i must appear in self.gens_index, and the slash gives the
        ##  matrix action.

#        self._glue = glue_data           ## TBA... =)

    def equivalent_index(self, A):
        ky = self._P.normalize(A[t10],A[t11])
        return self._equiv_ind[ky]

    def equivalent_rep(self, A):
        """
        Returns a coset representative that is equivalent to A modulo `\Gamma_0(N)`.

        INPUT:

        - ``A`` -- a matrix in `SL_2(\ZZ)`

        OUTPUT:

        - a matrix in `SL_2(\ZZ)` congruent to ``A`` modulo `\Gamma_0(N)`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
            sage: M2Z = MatrixSpace_ZZ_2x2()
            sage: A = M2Z([5,3,38,23])
            sage: ManinRelations(60).equivalent_rep(A)
            [-7 -3]
            [26 11]
        """
        ky = self._P.normalize(A[t10],A[t11])
        return self._equiv_rep[ky]

    def gens_index(self):
        return self._gens_index

    def find_coset_rep(self, A):
        #fix this to use a dict
        i = self._P.index(A[t10],A[t11])
        m = self._P1_to_mats[i]
        return self._mats[m]

    def P1(self):
        r"""
        Returns the Sage representation of `P^1(\ZZ/N\ZZZ)`.

        OUTPUT:

        - `P^1(Z/NZ)` where N is the level of the relations.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.P1()
            The projective line over the integers modulo 11
        """
        return self._P

    def P1_to_coset_index(self,n=None):
        r"""
        Takes the n-th element of Sage's `P^1(Z/NZ)` and returns the
        index of the associated element in A.coset_reps().

        Here by associated we mean the unique coset rep whose bottom
        row corresponds to the n-th element of P^1.  If n is not
        specified the entire translation table between the Sage P^1
        and the coset reps P^1 is returned.

        INPUT:

        - ``n`` -- integer (default: None)

        OUTPUT:

        - The unique integer j satisfying that the bottom row of A.coset_reps(j) is equivalent to A.P1()[n]

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: P = A.P1()
            sage: ind = 6
            sage: a = P[ind]; a
            (1, 5)
            sage: ind2 = A.P1_to_coset_index(ind); ind2
            7
            sage: b = A.coset_reps(ind2); b
            [ 1  0]
            [-2  1]
            sage: P.index(a[0],a[1]) == P.index(b[1,0],b[1,1])
            True
        """
        if n is None:
            return self._P1_to_mats
        else:
            return self._P1_to_mats[n]

    def generator_indices(self,n=None):
        r"""
        Returns the indices of coset reps which were chosen as our
        generators.

        In particular, the associate divisors of these coset reps
        generator all divisors over Z[Gamma_0(N)], and thus a modular
        symbol is uniquely determined by its values on these divisors.

        INPUT:

        - ``n`` -- integer (default: None)

        OUTPUT:

        - The list of indices in self.coset_reps() of our generating
          set.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.generator_indices()
            [0, 2, 3]
            sage: A.generator_indices(2)
            3
            sage: A = ManinRelations(13)
            sage: A.generator_indices()
            [0, 2, 3, 4, 5]
            sage: A = ManinRelations(101)
            sage: A.generator_indices()
            [0, 2, 3, 4, 5, 6, 8, 9, 11, 13, 14, 16, 17, 19, 20, 23, 24, 26, 28]
        """
        if n is None:
            return self._gens_index
        else:
            return self._gens_index[n]

    def two_torsion_indices(self,n=None):
        r"""
        Returns indices of coset rep generators which are fixed by
        Gamma_0(N) 2-torsion.

        INPUT:

        - ``n`` -- integer (default: None)

        OUTPUT:

        - The list of indices in self.coset_reps() whose associated
          unimodular path contains a point fixed by a Gamma_0(N)
          element of order 2 (where the order is computed in
          `PSL_2(Z)`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.two_torsion_indices()
            []
            sage: A = ManinRelations(13)
            sage: A.two_torsion_indices()
            [3, 4]
            sage: A = ManinRelations(17)
            sage: A.two_torsion_indices()
            [5, 7]
        """
        if n is None:
            return self._twotor_index
        else:
            return self._twotor_index[n]

    def two_torsion_relation_matrices(self,n=None):
        r"""
        Returns the order 2 matrices corresponding to the indices in
        self.two_torsion_indices()

        INPUT:

        - ``n`` -- integer (default: None)

        OUTPUT:

        - The list of order 2 matrices which correspond to the indices
          in self.two_torsion_indices().

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(13)
            sage: A.two_torsion_relation_matrices()
            [
            [  5   2]  [  8   5]
            [-13  -5], [-13  -8]
            ]
        """
        if n is None:
            return self._twotorrels
        else:
            return self._twotorrels[n]

    def three_torsion_indices(self,n=None):
        r"""
        Returns indices of coset rep generators which are fixed by
        `\Gamma_0(N)` 3-torsion.

        INPUT:

        - ``n`` -- integer (default: None)

        OUTPUT:

        - The list of indices in self.coset_reps() whose associated
          unimodular path contains a point fixed by a `\Gamma_0(N)`
          element of order 3 in the ideal triangle directly below that
          path.  Here the order is actually computed in `PSL_2(Z)`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.three_torsion_indices()
            []
            sage: A = ManinRelations(13)
            sage: A.three_torsion_indices()
            [2, 5]
            sage: A = ManinRelations(17)
            sage: A.three_torsion_indices()
            []
            sage: A = ManinRelations(103)
            sage: A.three_torsion_indices()
            [16, 17]
        """
        if n is None:
            return self._threetor_index
        else:
            return self._threetor_index[n]

    def three_torsion_relation_matrices(self,n=None):
        r"""
        Returns the order 3 matrices corresponding to the indices in
        self.three_torsion_indices()

        INPUT:

        - ``n`` -- integer (default: None)

        OUTPUT:

        - The list of order 3 matrices which correspond to the indices
          in ``self.three_torsion_indices()``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(13)
            sage: A.three_torsion_relation_matrices()
            [
            [-4 -1]  [-10  -7]
            [13  3], [ 13   9]
            ]

        """
        if n is None:
            return self._threetorrels
        else:
            return self._threetorrels[n]

    def coset_relations(self,n=None):
        r"""
        Expresses the divisor attached to the n-th coset rep in terms
        of our chosen generators.

        Returns a list of triples `(d, A, i)` such that the divisor
        attached to the n-th coset rep equals the sum over these
        triples of:

            `d * A^(-1) * (divisor attached to i-th coset rep)`

        Here the index `i` must appear in self.generator_indices().
        This formula will allow us to recover the value of a modular
        symbol on any coset rep in terms of its values on our
        generating set.

        If n is not specified, a list containing this info for all n
        is returned.

        INPUT:

        - ``n`` -- integer (default: None)

        OUTPUT:

        - A `\ZZ[\Gamma_0(N)]`-relation expressing the divisor
          attached to the `n`-th coset rep in terms of our generating
          set.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.generator_indices()
            [0, 2, 3]
            sage: A.coset_relations(0)
            [(1, [1 0]
            [0 1], 0)]
            sage: A.coset_relations(2)
            [(1, [1 0]
            [0 1], 2)]
            sage: A.coset_relations(3)
            [(1, [1 0]
            [0 1], 3)]
            sage: A.coset_relations(4)
            [(-1, [-3 -2]
            [11  7], 2)]
            sage: B=A.coset_relations(4)[0][1]; B
            [-3 -2]
            [11  7]
            sage: B^(-1)*A.coset_reps(2)
            [ 2 -1]
            [-3  2]
            sage: A.coset_reps(4)
            [-1 -2]
            [ 2  3]
            sage: from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
            sage: M2Z = MatrixSpace_ZZ_2x2()
            sage: sig = M2Z([0,1,-1,0])
            sage: B^(-1)*A.coset_reps(2) == A.coset_reps(4)*sig
            True

        """
        if n is None:
            return self._rels
        else:
            return self._rels[n]

    def form_list_of_cusps(self):
        r"""
        Returns the intersection of a fundamental domain for
        `\Gamma_0(N)` with the real axis.

        The construction of this fundamental domain follows the
        arguments of [PS] Section 2.  The boundary of this fundamental
        domain consists entirely of unimodular paths when
        `\Gamma_0(N)` has no elements of order 3.  (See [PS] Section
        2.5 for the case when there are elements of order 3.)

        OUTPUT:

        - A sorted list of rational numbers marking the intersection
          of a fundamental domain for `\Gamma_0(N)` with the real
          axis.


        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.form_list_of_cusps()
            [-1, -2/3, -1/2, -1/3, 0]
            sage: A = ManinRelations(13)
            sage: A.form_list_of_cusps()
            [-1, -2/3, -1/2, -1/3, 0]
            sage: A = ManinRelations(101)
            sage: A.form_list_of_cusps()
            [-1, -6/7, -5/6, -4/5, -7/9, -3/4, -11/15, -8/11, -5/7, -7/10, -9/13, -2/3, -5/8, -13/21, -8/13, -3/5, -7/12, -11/19, -4/7, -1/2, -4/9, -3/7, -5/12, -7/17, -2/5, -3/8, -4/11, -1/3, -2/7, -3/11, -1/4, -2/9, -1/5, -1/6, 0]
        """
        ## Get the level
        N = self.level()

        ## Checks that the level N is a positive integer > 1
        if not ((N in ZZ) and (N > 1)):
            raise TypeError, "Error in form_list_of_cusps: level should be > 1"

        ## Some convenient shortcuts

        P = self.P1()
        sP = len(P.list())   ## Size of P^1(Z/NZ)

        ## Initialize some lists

        C = [QQ(-1),"?",QQ(0)]

        ## Initialize the list of cusps at the bottom of the fund. domain.
        ## The ? denotes that it has not yet been checked if more cusps need
        ## to be added between the surrounding cusps.

        full_domain = False     ## Says that we're not done yet!

        v = [False for r in range(sP)]
        ## This initializes a list indexed by P^1(Z/NZ) which keeps track of
        ## which right coset representatives we've found for Gamma_0(N)/SL_2(Z)
        ## thru the construction of a fundamental domain

        ## Includeds the coset repns formed by the original ideal triangle
        ## (with corners at -1, 0, infty)

        v[P.index(0,1)] = True
        v[P.index(1,-1)] = True
        v[P.index(-1,0)] = True


        ## Main Loop -- Ideal Triangle Flipping
        ## ====================================
        while (not full_domain):
            full_domain = True

            ## This loop runs through the current set of cusps
            ## and checks to see if more cusps should be added
            ## -----------------------------------------------
            for s in range(1, len(C), 2):  ## range over odd indices in the
                                           ## final list C
                if C[s] == "?":

                    ## Single out our two cusps (path from cusp2 to cusp1)
                    cusp1 = C[s-1]
                    cusp2 = C[s+1]

                    ## Makes the unimodular transform for the path from cusp2
                    ## to cusp1

                    b1 = cusp1.denominator()
                    b2 = cusp2.denominator()

                    ## This is the point where it is determined whether
                    ## or not the adjacent triangle should be added
                    ## ------------------------------------------------
                    pos = P.index(b2,b1)   ## The Sage index of the bottom
                                                 ## row of our unimodular
                                           ## transformation gam

                    ## Check if we need to flip (since this P1 element has not
                    ## yet been accounted for!)
                    if v[pos] == False:
                        v[pos] = True      ## Say this P1 element now occurs
                        v[P.index(b1,-(b1+b2))] = True ## Say that the other
                                                       ## two ideal triangle
                                                       ## edges also occur!
                        v[P.index(-(b1+b2),b2)] = True

                        ## Check to see if this triangle contains a fixed
                        ## point by an element of Gamma_0(N).  If such an
                        ## element is present, the fundamental domain can be
                        ## extended no further.

                        if (b1**2 + b2**2 + b1*b2)%N != 0:

                        ## this congruence is exactly equivalent to
                        ## gam * [0 -1; 1 -1] * gam^(-1) is in Gamma_0(N)
                        ## where gam is the matrix corresponding to the
                        ## unimodular path connecting cusp1 to cusp2

                            C[s] = "i"  ## The '?' is changed to an 'i'
                                            ## indicating that a new cusp needs to
                                        ##  be inserted here
                            full_domain = False
                        else:
                            C[s] = "x"  ## The '?' is changed to an 'x' and no
                                               ##  more checking below is needed! =)
                    else:
                        C[s] = "x"  ## The '?' is changed to an 'x' and no more
                                           ## checking below is needed! =)


            ## Now insert the missing cusps (where there is an 'i' in the
            ## final list)
            ## This will keep the fundamental domain as flat as possible!
            ## ---------------------------------------------------------------

            s=1
            while s < len(C):    ## range over odd indices in the final list C
                if C[s] == "i":
                    C[s]="?"

                    ## Single out our two cusps (path from cusp2 to cusp1)
                    cusp1 = C[s-1]
                    cusp2 = C[s+1]

                    ## Makes the unimodular transform for the path from cusp2
                    ## to cusp1
                    a1 = cusp1.numerator()
                    b1 = cusp1.denominator()
                    a2 = cusp2.numerator()
                    b2 = cusp2.denominator()

                    ## Inserts the Farey center of these two cusps!
                    a = a1 + a2
                    b = b1 + b2
                    C.insert(s+1, a/b)
                    C.insert(s+2, "?")
                    s = s+2
                s = s+2

        ## Remove the (now superfluous) extra string characters that appear
        ## in the odd list entries
        C = [QQ(C[s]) for s in range(0,len(C),2)]
        return C

    def is_unimodular_path(self, r1, r2):
        r"""
        Determines whether two (non-infinite) cusps are connected by a
        unimodular path.

        INPUT:

        - ``r1, r2`` -- rational numbers

        OUTPUT:

        - A boolean expressing whether or not a unimodular path
          connects r1 to r2.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.is_unimodular_path(0,1/3)
            True
            sage: A.is_unimodular_path(1/3,0)
            True
            sage: A.is_unimodular_path(0,2/3)
            False
            sage: A.is_unimodular_path(2/3,0)
            False
        """
        a = r1.numerator()
        b = r2.numerator()
        c = r1.denominator()
        d = r2.denominator()
        return (a*d - b*c)**2 == 1


    def unimod_to_matrices(self, r1, r2):
        r"""
        Returns the two matrices whose associated unimodular paths
        connect `r1 -> r2` and `r2 -> r1`, respectively.

        INPUT:

        - ``r1, r2`` -- rational numbers (that are assumed to be
          related by a unimodular path)

        OUTPUT:

        - a pair of `2 x 2` matrices of determinant 1

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.unimod_to_matrices(0,1/3)
            (
            [ 0  1]  [1 0]
            [-1  3], [3 1]
            )

        """
        a = r1.numerator()
        b = r2.numerator()
        c = r1.denominator()
        d = r2.denominator()
        if (a*d-b*c)==1:
            ans = M2Z([a,b,c,d]), M2Z([-b,a,-d,c])
        else:
            ans = M2Z([-a,b,-c,d]), M2Z([b,a,d,c])
        ans[0].set_immutable()
        ans[1].set_immutable()
        return ans

    def fd_boundary(self,C):
        r"""
        Finds matrices whose associated unimodular paths give the
        boundary of a fundamental domain.

        Here the fundamental domain is for `\Gamma_0(N)`.  (In the
        case when `\Gamma_0(N)` has elements of order three the shape
        cut out by these unimodular matrices is a little smaller than
        a fundamental domain.  See `\S2.5` of Pollack-Stevens.)

        INPUT:

        - a list of rational numbers coming from
          self.form_list_of_cusps()

        OUTPUT:

        - a list of `2 x 2` integer matrices of determinant 1 whose
          associated unimodular paths give the boundary of a
          fundamental domain for `Gamma_0(N)` (or nearly so in the
          case of 3-torsion).

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: C = A.form_list_of_cusps(); C
            [-1, -2/3, -1/2, -1/3, 0]
            sage: A.fd_boundary(C)
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]
            [0 1], [-1  0], [ 1  3], [ 3  2], [ 2  3], [ 3  1]
            ]
            sage: A = ManinRelations(13)
            sage: C = A.form_list_of_cusps(); C
            [-1, -2/3, -1/2, -1/3, 0]
            sage: A.fd_boundary(C)
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]
            [0 1], [-1  0], [ 1  3], [ 3  2], [ 2  3], [ 3  1]
            ]
            sage: A = ManinRelations(101)
            sage: C = A.form_list_of_cusps(); C
            [-1, -6/7, -5/6, -4/5, -7/9, -3/4, -11/15, -8/11, -5/7, -7/10, -9/13, -2/3, -5/8, -13/21, -8/13, -3/5, -7/12, -11/19, -4/7, -1/2, -4/9, -3/7, -5/12, -7/17, -2/5, -3/8, -4/11, -1/3, -2/7, -3/11, -1/4, -2/9, -1/5, -1/6, 0]
            sage: A.fd_boundary(C)
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]  [-1 -3]  [-3 -2]
            [0 1], [-1  0], [ 1  6], [ 6  5], [ 5  9], [ 9  4], [ 4 11], [11  7],
            <BLANKLINE>
            [-2 -1]  [-1 -4]  [-4 -3]  [-3 -2]  [-2 -7]  [-7 -5]  [-5 -3]  [-3 -4]
            [ 7  3], [ 3 11], [11  8], [ 8  5], [ 5 17], [17 12], [12  7], [ 7  9],
            <BLANKLINE>
            [-4 -1]  [-1 -4]  [ -4 -11]  [-11  -7]  [-7 -3]  [-3 -8]  [ -8 -13]
            [ 9  2], [ 2  7], [  7  19], [ 19  12], [12  5], [ 5 13], [ 13  21],
            <BLANKLINE>
            [-13  -5]  [-5 -2]  [-2 -9]  [-9 -7]  [-7 -5]  [-5 -8]  [ -8 -11]
            [ 21   8], [ 8  3], [ 3 13], [13 10], [10  7], [ 7 11], [ 11  15],
            <BLANKLINE>
            [-11  -3]  [-3 -7]  [-7 -4]  [-4 -5]  [-5 -6]  [-6 -1]
            [ 15   4], [ 4  9], [ 9  5], [ 5  6], [ 6  7], [ 7  1]
            ]
        """

        C.reverse() ## Reverse here to get clockwise orientation of boundary

        ## These matrices correspond to the paths from infty to 0 and -1 to infty
        mats = [M2Z([1,0,0,1]),M2Z([1,1,-1,0])]
        mats[0].set_immutable()
        mats[1].set_immutable()

        ## Now find SL_2(Z) matrices whose associated unimodular paths connect
        ## the cusps listed in C.
        ## --------------------------------------------------------
        for j in range(len(C)-1):
            a = C[j].numerator()
            b = C[j+1].numerator()
            c = C[j].denominator()
            d = C[j+1].denominator()
            new_mat = M2Z([a,b,c,d])
            new_mat.set_immutable()
            mats.append(new_mat)

        return mats
