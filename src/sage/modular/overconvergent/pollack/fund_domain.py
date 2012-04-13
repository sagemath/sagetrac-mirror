
######################################################################
##
##  Code to create the Manin Relations class, which solves the 
##  "Manin relations".  That is, a description of Div^0(P^1(Q))
##  as a Z[Gamma_0(N)]-module in terms of generators and relations
##  is found.  The method used is geometric, constructing a nice
##  fundamental domain for Gamma_0(N) and reading the relevant Manin
##  relations off of that picture.  The algorithm follows the paper 
##  of Pollack and Stevens "Overconvergent modular symbols and p-adic 
##  L-functions"
##
##  Copyright (c) 2012, Rob Pollack and Jonathan Hanke
##      <rpollack@math.bu.edu>
##      <jonhanke@gmail.com>
##
##  Released under the GNU Public License, 2012.
##
######################################################################

from sage.matrix.all import Matrix
from sage.modular.modsym.all import P1List
from sage.rings.integer_ring import ZZ
     
	
######################################
##  Define the Manin Relation Class ##
######################################

class manin_relations(SageObject):
    def __init__(self, N):
        r"""
        
        INPUT:
           ``N`` -- a positive integer

        OUTPUT:
            none

        EXAMPLES:

        """
        ## Make the identity matrix
        Id = Matrix(2,2,[1,0,0,1])

        ## Store the level
        self.__N = N

        ## Creates and stores the Sage representation of P^1(Z/NZ)
        P = P1List(N)
        self.__P = P

        ## Solve the Manin Relations!
        ## ==========================

        ## Creates a fundamental domain for Gamma_0(N), recording the cusps 
	## on the real line and the associated subset of the (geometrically 
	## constructed) right coset representatives for Gamma_0(N) / SL_2(Z)
        ## which correspond to the boundary of the fundamental domain
        cusps = self.form_list_of_cusps()
        coset_reps = self.mat_boundary_of_fd(cusps)

        ## Gives names to matrices of order 2 and 3 (in PSL_2)
        sig = Matrix(2,2,[0,1,-1,0])
        tau = Matrix(2,2,[0,-1,1,-1])

        ## Creates P^1 from the geometric coset representatives arising from 
	## fundamental domain
        p1s = [(coset_reps[j])[1] for j in range(len(coset_reps))]  
	## this gets the bottom rows of the 2x2 matrices

        ## Initializes relevant Manin data
        gens = []
        twotor = []
        twotorrels = []
        threetor = []
        threetorrels = [] 


        ## --------------------------------------------------------------------
        rels = [0 for i in range(0,len(coset_reps))]  

	## entries of rel will be lists of elements of the form (c,A,r) 
	## with c a constant, A a Gamma_0(N) matrix, and r the index of a 
	## generator.  The meaning is that the value of a modular symbol phi 
	## on the (associated divisor) of the i-th element of coset_reps 
        ## will be the sum of c * phi (r-th genetator) | A 
	## as one varies over the tuples in rel[i]

        boundary_checked = [False for i in range(0,len(coset_reps))]  
	
	## This list keeps track of which boundary pieces of the 
	## fundamental domain have been already used as we are picking 
	## our generators
 
        glue_data = [0 for i in range(0,len(coset_reps))]   ## ????

        ## The following loop will choose our generators by picking one edge 
	## out of each pair of edges that are glued to each other and picking 
	## each edge glued to itself (arising from two-torsion)
        
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

                    gens = gens + [r]    ## the index r is adding to our list
		    	   	  	 ## of indexes of generators
                    twotor = twotor + [r]  ## the index r is adding to our 
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

                    glue_data[r]=(r,gam)  

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

                        gens = gens + [r]   
			## the index r is adding to our list of indexes 
			## of generators
                        
			rels[r] = [(1,Id,r)]  

			## this relation expresses the fact that coset_reps[r]
			## is one of our basic generators
                        threetor = threetor + [r]    
			
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
                        A = Matrix(2,2,[-b,a,-d,c])
                        coset_reps = coset_reps + [A]  
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
                                
				gens = gens + [r]         
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
                               
			        glue_data[r] = (s,gam)  ## ????
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
                    
		    for j in range(s-r):  
		    ## Running between the cusps between cusp1 and cusp2
                        vA = vA + rels[r+j+2]  ## Edge relation added
                        t = (-rels[r+j+2][0][0],rels[r+j+2][0][1],rels[r+j+2][0][2]) ## This is simply the negative of the above edge relation.
                        vB = vB + [t]  ## Negative of edge relation added
                    rels = rels + [vA,vB]  ## Relations for A and B adding to relations list
####        return [coset_reps,gens,twotor,twotorrels,threetor,threetorrels,rels,glue_data]


        ## Store the data coming from solving the Manin Relations
        ## ======================================================
        
	self.__mats = coset_reps              
	
	## Coset representatives of Gamma_0(N) coming from the geometric 
	## fundamental domain algorithm

        ## Make the translation table between the Sage and Geometric 
	## descriptions of P^1
        
	v = zero_vector(len(coset_reps))
        for r in range(len(coset_reps)):  
            v[P.index(coset_reps[r][1,0], coset_reps[r][1,1])] = r
        self.__P1_to_mats = v  

        self.__gens = gens            
	## This is a list of indices of the (geometric) coset representatives 
	## whose values (on the associated degree zero divisors) determine the
	##  modular symbol.

        self.__twotor = twotor              
	## A list of indices of the (geometric) coset representatives whose 
	## paths are identified by some 2-torsion element (which switches the 
	## path orientation)

        self.__twotorrels = twotorrels      
	## A list of (2-torsion in PSL_2(Z)) matrices in Gamma_0(N) that give 
        ## the orientation identification in the paths listed in twotor above!

        self.__threetor = threetor            
	## A list of indices of the (geometric) coset representatives that 
	## form one side of an ideal triangle with an interior fixed point of 
	## a 3-torsion element of Gamma_0(N)        

        self.__threetorrels = threetorrels    
	## A list of (3-torsion in PSL_2(Z)) matrices in Gamma_0(N) that give 
        ## the interior fixed point described in threetor above!

        self.__rels = rels      
	## A list of lists of triples (d, A, i), one for each coset 
	## representative of Gamma_0(N) (ordered to correspond to the 
	## representatives of self.coset_reps) expressing the value of a 
        ## modular symbol on the associated unimodular path as a sum of terms 
        ##    d * (value on the i-th coset rep) | A
        ## where the index i must appear in self.gens, and the slash gives the
	##  matrix action.

        self.__glue = glue_data           ## TBA... =)


    def level(self):
	r"""
	Returns the level `N` of `\Gamma_0(N)` that we work with.

        INPUT:
            - none

        OUTPUT:
        
        The integer `N` of the group `\Gamma_0(N)` for which the Manin Relations are being computed.

        EXAMPLES:                        
        ::        
            sage: A = manin_relations(11)
            sage: A.level()
            11

	"""
	return self.__N

    def coset_reps(self,n=None):
	r"""
	Returns the n-th coset rep associated with our fundamental domain or all coset reps if n is not specified.

        INPUT:
            - ``n`` -- integer (default: None) 

        OUTPUT:
        
        If n is given then the n-th coset representative is returned and otherwise all coset reps are returned.

        EXAMPLES:                        
        ::        
            sage: A.coset_reps(0)
            [1 0]
            [0 1]
            sage: A.coset_reps(1)
            [ 1  1]
            [-1  0]
            sage: A.coset_reps(2)
            [ 0 -1]
            [ 1  3]
            sage: A.coset_reps()
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]  [ 1 -1]  [-1 -1]
            [0 1], [-1  0], [ 1  3], [ 3  2], [ 2  3], [ 3  1], [-1  2], [ 2  1],

            [ 1  0]  [ 0 -1]  [ 1  0]  [ 0 -1]
            [-1  1], [ 1  1], [-2  1], [ 1  2]
            ]

	"""
        if n is None:
            return deepcopy(self.__mats)
        else:
            return deepcopy(self.__mats[n])


    def P1(self):
        r"""
        Returns the Sage representation of `P^1(Z/NZ)` that we work with.
        """
        return deepcopy(self.__P)

    def P1_to_mats(self,n=None):
        r"""
        Returns the translation table (or individual element) between the Sage and geometric 
	descriptions of `P^1(Z/NZ)`.
        """
        if n is None:
            return deepcopy(self.__P1_to_mats)
        else:
            return deepcopy(self.__P1_to_mats[n])

    def generator_indices(self,n=None):
        r"""
        Returns a list of indices of the (geometric) coset representatives 
	whose values (on the associated degree zero divisors) determine the 
	modular symbol.
        """
        if n is None:
            return deepcopy(self.__gens)
        else:
            return deepcopy(self.__gens[n])

    def two_torsion_indices(self,n=None):
        r"""
        Returns a list of indices pointing to the coset generators whose 
	paths are identified by some 2-torsion element (which switches the 
	path orientation).
        """
        if n is None:
            return deepcopy(self.__twotor)
        else:
            return deepcopy(self.__twotor[n])

    def two_torsion_relation_matrices(self,n=None):
        r"""
        Returns a list of (2-torsion in `PSL_2(Z)`) matrices in `\Gamma_0(N)` 
        that give the orientation identification in the paths indexed by 
	two_torsion_indices().
        """
        if n is None:
            return deepcopy(self.__twotorrels)
        else:
            return deepcopy(self.__twotorrels[n])

    def three_torsion_indices(self,n=None):
        r"""
        Returns a list of indices pointing to (geometric) coset 
	representatives that form one side of an ideal triangle with an 
	interior fixed point of a 3-torsion element of `\Gamma_0(N)`.
        """
        if n is None:
            return deepcopy(self.__threetor)
        else:
            return deepcopy(self.__threetor[n])

    def three_torsion_relation_matrices(self,n=None):
        r"""
        A list of (3-torsion in `PSL_2(Z)`) matrices in `\Gamma_0(N)` that 
        give the interior fixed point indexed by three_torsion_indices().
        """
        if n is None:
            return deepcopy(self.__threetor)
        else:
            return deepcopy(self.__threetor[n])

    def coset_relations(self,n=None):
        r"""
        Returns a list of lists of triples `(d, A, i)`, one for each coset 
	representative of `\Gamma_0(N)` (ordered to correspond to the 
	representatives of self.mats) expressing the value of a 
        modular symbol on the associated unimodular path as a sum of terms 

            `d * (value on the i-th representative) | A`

        where the index `i` must appear in self.gens, and the slash gives 
	the matrix action.
        """
        if n is None:
            return deepcopy(self.__rels)
        else:
            return deepcopy(self.__rels[n])

    def form_list_of_cusps(self):
        r"""
        Constructs a fundamental domain for `\Gamma_0(N)` as in [PS] Section 2 
	(and Section 2.5 for the case of 3-torsion) made up of unimodular 
	paths and returns the list of rational numbers which mark where 
	this fundamental domain meets the real axis.
    
        INPUT:
            none
        
        OUTPUT:
            a list of rational numbers

        EXAMPLES:
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
            full_domain=true
    
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
                            full_domain = false
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
        C = [Rational(C[s]) for s in range(0,len(C),2)]
        return C

    def is_unimodular_path(self, r1, r2):
        r"""
        Determines whether two (non-infinite) cusps are connected by a 
	unimodular path.
    
        INPUT:
            ``r1, r2`` -- rational numbers
        
        OUTPUT:
            boolean

        EXAMPLES:
        """
        a = r1.numerator()
        b = r2.numerator()
        c = r1.denominator()
        d = r2.denominator()
        return (a*d - b*c)**2 == 1


    def unimod_to_matrices(self, r1, r2):
        r"""
	Returns the two matrices whose associated unimodular paths connect 
	`r1 -> r2` and `r2 -> r1`, respectively.
    
        INPUT:
            ``r1, r2`` -- rational numbers (that are assumed to be related by 
	    a unimodular path!)
        
        OUTPUT:
            a pair of `2 x 2` matrix of determinant 1

        EXAMPLES:
        """
        a = r1.numerator()
        b = r2.numerator()
        c = r1.denominator()
        d = r2.denominator()
        if (a*d-b*c)==1:
            return Matrix(2,2,[a,b,c,d]), Matrix(2,2,[-b,a,-d,c])
        else:
            return Matrix(2,2,[-a,b,-c,d]), Matrix(2,2,[b,a,d,c])


    def mat_boundary_of_fd(self,C):
        r"""
        Returns a list of matrices whose associated unimodular paths gives 
	the boundary of a fundamental domain for `\Gamma_0(N)`.  
        (In the case when `\Gamma_0(N)` has three torsion the shape is slightly
        smaller than a fundamental domain.  See `\S2.5` of Pollack-Stevens.)
	
        INPUT:
            a list of rational numbers coming from form_list_of_cusps
        
        OUTPUT:
            a list of `2 x 2` integer matrices of determinant 1
	      
        EXAMPLES:
        """
        ## Make the identity matrix
        Id = Matrix(2,2,[1,0,0,1])

        C = self.form_list_of_cusps()
			
        ## Now convert the finished list of cusps to matrices whose 
        ## associated unimodular paths form the fundamental domain
        ## --------------------------------------------------------

        C.reverse() ## ??????  Reverse the for reasons that I no longer 
	            ## remember but without this it crashes.
		    
	## these matrices correspond to the paths from infty to 0
	## and -1 to infty

        mats = [Id,Matrix(2,2,[1,1,-1,0])]  
        for j in range(len(C)-1):
            a = C[j].numerator()
            b = C[j+1].numerator()
            c = C[j].denominator()
            d = C[j+1].denominator()
            mats = mats + [Matrix(2,2,[a,b,c,d])]

        return mats
