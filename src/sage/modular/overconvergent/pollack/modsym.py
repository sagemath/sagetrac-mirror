from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0

@cached_function
def invert(a,b,c,d):
    """
    Takes the numbers a,b,c,d and returns the inverse to the matrix [a,b;c,d].  
    Here ad-bc is assumed to be 1
        
    INPUT:
        a,b,c,d -- integers such that ad-bc=1
        
    OUTPUT:
        a 2x2 integer matrix

    EXAMPLES:
    """
    return Matrix(2,2,[d,-b,-c,a])

@cached_function
def unimod_matrices_to_infty(r,s):
    """
    Returns a list of matrices whose associated unimodular paths connect r/s to infty.
    (This is Manin's continued fraction trick.)
        
    INPUT:
        r,s -- rational numbers
        
    OUTPUT:
        a list of SL_2(Z) matrices

    EXAMPLES:
    """
    if s!=0:
        v = []  ## Initializes the list of matrices
# the function contfrac_q in https://github.com/williamstein/psage/blob/master/psage/modform/rational/modular_symbol_map.pyx
# is very, very relevant to massively optimizing this.
        list = convergents(r/s)  ## Computes the continued fraction convergents of r/s
        ## This loop 
        for j in range(0,len(list)-1):
            a = list[j].numerator()
            c = list[j].denominator()
            b = list[j+1].numerator()
            d = list[j+1].denominator()
            v = v + [Matrix(ZZ,[[(-1)**(j+1)*a,b],[(-1)**(j+1)*c,d]])]  ## The matrix connecting two consecutive convergents is added on
        return [Matrix(ZZ,[[1,list[0].numerator()],[0,list[0].denominator()]])]+v
    else:
        return []
        
def flip(A):
    """
    Takes the matrix [a,b;c,d] and returns [-c,a;-d,b] -- i.e. reverses the orientation of the associated unimodular path        

    INPUT:
        A -- 2 x 2 matrix
        
    OUTPUT:
        2 x 2 matrix

    EXAMPLES:
    """
    return Matrix(2,2,[-A[0,1],A[0,0],-A[1,1],A[1,0]])

@cached_function
def unimod_matrices_from_infty(r,s):
    """
    Returns a list of matrices whose associated unimodular paths connect r/s to infty.
    (This is Manin's continued fraction trick.)
        
    INPUT:
        r,s -- rational numbers
        
    OUTPUT:
        a list of SL_2(Z) matrices

    EXAMPLES:
    """        
    if s != 0:
        v = []
# the function contfrac_q in https://github.com/williamstein/psage/blob/master/psage/modform/rational/modular_symbol_map.pyx
# is very, very relevant to massively optimizing this.
        list = convergents(QQ(r)/s)  
        for j in range(0,len(list)-1):
            a = list[j].numerator()
            c = list[j].denominator()
            b = list[j+1].numerator()
            d = list[j+1].denominator()
            v = v + [flip(Matrix(ZZ,[[(-1)**(j+1)*a,b],[(-1)**(j+1)*c,d]]))]
        return [flip(Matrix(ZZ,[[1,list[0].numerator()],[0,list[0].denominator()]]))]+v
    else:
        return []


def basic_hecke_matrix(a,ell):
    """
    Returns the matrix [1,a;0,ell] (if a<ell) and [ell,0;0,1] if a>=ell

    INPUT:
        a -- an integer or Infinity
        ell -- a prime
        
    OUTPUT:
        a 2 x 2 matrix of determinant ell

    EXAMPLES:
    """
    if a < ell:
        return Matrix(2,2,[1,a,0,ell])
    else:
        return Matrix(2,2,[ell,0,0,1])

@cached_function
def prep_hecke_individual(ell,M,m):
    """
    This function does some precomputations needed to compute T_ell.  In particular,
    if phi is a modular symbol and D_m is the divisor associated to our m-th chosen 
    generator, to compute (phi|T_ell)(D_m) one needs to compute phi(gam_a D_m)|gam_a where
    gam_a run thru the ell+1 matrices defining T_ell.  One then takes gam_a D_m and writes it
    as a sum of unimodular divisors.  For each such unimodular divisor, say [M] where M is a
    SL_2 matrix, we then write M=gam*M_i where gam is in Gamma_0(N) and M_i is one of our
    chosen coset representatives.  Then phi([M]) = phi([M_i]) | gam^(-1).  Thus, one has
    
        (phi | gam_a)(D_m) = sum_i sum_j phi([M_i]) | gam_{ij}^(-1) * gam_a

    as i runs over the indices of all coset representatives and j simply runs over however many
    times M_i appears in the above computation.  

    Finally, the output of this function is a list L enumerated by the coset representatives 
    in M.coset_reps() where each element of this list is a list of matrices, and the entries of L
    satisfy:

        L[i][j] = gam_{ij} * gam_a  

    INPUT:
        ell -- a prime
        M -- Manin relations of level N
        m -- index of a generator 
        
    OUTPUT:
	A list of lists (see above).

    EXAMPLES:
    """

    N = M.level()
    ans = [[] for a in range(len(M.coset_reps()))]  ## this will be the list L above enumerated by coset reps

    ##  This loop will runs thru the ell+1 (or ell) matrices defining T_ell of the form [1,a,0,ell] and carry out the computation
    ##  described above.
    ##  -------------------------------------
    for a in range(ell+1):
	if (a<ell) or (N%ell!=0):   ## if the level is not prime to ell the matrix [ell,0,0,1] is avoided.
           gama = basic_hecke_matrix(a,ell)
           t = gama*M.coset_reps(M.generator_indices(m))  ##  In the notation above this is gam_a * D_m
           v = unimod_matrices_from_infty(t[0,0],t[1,0]) + unimod_matrices_to_infty(t[0,1],t[1,1])  ##  This expresses t as a sum of unimodular divisors

           ## This loop runs over each such unimodular divisor
           ## ------------------------------------------------
           for b in range(len(v)):
               A = v[b]    ##  A is the b-th unimodular divisor
               i = M.P1().index(A[1,0],A[1,1])      ##  i is the index in SAGE P1 of the bottom row of A
               j = M.P1_to_mats(i)                  ##  j is the index of our coset rep equivalent to A
               B = M.coset_reps(j)                    ##  B is that coset rep
               C = invert(A[0,0],A[0,1],A[1,0],A[1,1])   ##  C equals A^(-1).  This is much faster than just inverting thru SAGE
               gaminv = B * C                              ##  gaminv = B*A^(-1)
               ans[j] = ans[j] + [gaminv * gama]             ##  The matrix gaminv * gama is added to our list in the j-th slot (as described above)

    return ans

@cached_function
def prep_hecke(ell,M):
    """
    Carries out prep_hecke_individual for each generator index and puts all of the answers in a long list.

    INPUT:
        ell -- a prime
        M -- Manin relations of level N
        
    OUTPUT:
	A list of lists of lists

    EXAMPLES:
    """
    ans = []
    for m in range(len(M.generator_indices())):
        ans = ans + [prep_hecke_individual(ell,M,m)]
    return ans
            

##############################
##  Define the modsym class ##
##############################

##  This class represents V-valued modular symbols on Gamma_0(N) where M satisfies:
##      1) V is a Z-module
##      2) V has a RIGHT action by some collections of matrices large enough that it contains
##         a congruence subgroup and the matrices [1,a,0,q] for all primes q and, 
##         [q,0,0,1] for q not dividing N
##
##  A modular symbol contains the info of the Manin Relations (say M) for Gamma_0(N).  
##  A symbol phi is stored by its values on our chosen generators -- i.e. those listed 
##  in M.generator_indices().  

class modsym(SageObject):
    def __init__(self,data,manin,full_data=None):
        r"""
        
        INPUT:
           data -- the list of values of the modular symbol on our specified generator set
           manin -- Manin Relations
           full_data -- the list of values of the modular symbol on all coset representatives

        OUTPUT:
            none

        EXAMPLES:

        """
        self.__manin = manin  
        self._data = data      
        if full_data!=None:
            self._full_data = full_data
        else:
            self._full_data = 0
            
    def manin(self):
        r"""
        
        INPUT:
            none

        OUTPUT:
            The manin relations defining the symbol.  

        EXAMPLES:

        """
        return deepcopy(self.__manin)

    def data(self,n=None):
        r"""
        Returns the list of values of self on our generating set or simply the n-th value if n is specified.
        
        INPUT:
            n -- integer

        OUTPUT:
            The list of values of self on our generating set or simply the n-th value if n is specified.

        EXAMPLES:

        """
        if n is None:
	        return deepcopy(self._data)
        else:
	        return deepcopy(self._data[n])

    def full_data(self,n=None):
        r"""
        Returns the value of self on all coset reps (or simply on the n-th coset rep if n is specified).

        INPUT:
            n -- integer

        OUTPUT:
            The list of values of self on all coset reps or simply the n-th value if n is specified.

        EXAMPLES:

        """
        if n is None:
	        return deepcopy(self._full_data)
        else:
	        return deepcopy(self._full_data[n])

    def __repr__(self):
        r"""
        
        INPUT:
            none

        OUTPUT:
            The representation of self.data.  

        EXAMPLES:

        """
        return repr(self.data())

    def ngens(self):
        r"""
        Returns the number of generators defining self.
        
        INPUT:
            none

        OUTPUT:
            The number of generators defining our modular symbols.  

        EXAMPLES:

        """
        return len(self.manin().generator_indices())

    def ncoset_reps(self):
        r"""
        Returns the number of coset representatives defining the full_data of phi.
        
        INPUT:
            none

        OUTPUT:
            The number of coset representatives stored in the manin relations.  (Just the size of P^1(Z/NZ))

        EXAMPLES:

        """
        return len(self.manin().coset_reps())

    def __add__(self,right):
        r"""
        Returns the sum of self and right.

        INPUT:
            right -- a modular symbol of the same level 

        OUTPUT:
            the sum of self and right

        EXAMPLES:

        """
        v = []  ## This list will store the values of the modular symbol on a generating set.
        ## This loop runs thru each of the values of the modular symbol and adds them together.
        for j in range(0,self.ngens()):
            v = v + [self.data(j)+right.data(j)]
        ## If both of the symbols already have their values computed on all coset reps then these
        ## these values are added together as well.
        if self.full_data()!=0 and right.full_data()!=0:
            w = []  ## This list will store the values of the modular symbol on all coset reps.
            for j in range(0,self.ncoset_reps()):
                w = w + [self.full_data(j) + right.full_data(j)]
        else:
            w = 0
            
        C = type(self)
        return C(v,self.manin(),w) ##  Coercing into C keeps track of which kind of modular symbol we are working with

    def normalize(self):
        r"""
        Normalizes all of the values of the symbol self.

        INPUT:
            none

        OUTPUT:
            none

        EXAMPLES:

        """
        for j in range(self.ngens()):
            self._data[j].normalize()

        if self.full_data()!=0:
            for j in range(self.ncoset_reps()):
                self._full_data[j].normalize()

    def scale(self,left):
        r"""
        Returns left * self

        INPUT:
            left -- a scalar

        OUTPUT:
            left * self

        EXAMPLES:

        """
        v = []    ##  This will be the list of values of the answer.
        for j in range(0,self.ngens()):
            v = v + [self.data(j).scale(left)]
        if self.full_data()!=0:
            w=[]
            for j in range(0,self.ncoset_reps()):
                w = w + [self.full_data(j).scale(left)]
        else:
            w = 0

        C = type(self)
        return C(v,self.manin(),w) ##  Coercing into C keeps track of which kind of modular symbol we are working with
        
    def __sub__(self,right):
        r"""
        Returns self minus right

        INPUT:
            right -- a modular symbol of the same level 

        OUTPUT:
            self minus right

        EXAMPLES:

        """
        return self + right.scale(-1)

    def __cmp__(self,right):
        return cmp((self.level,self.data()),(right.level,right.data()))

    def zero_elt(self):
        r"""
        Returns the zero element of the space where self takes values.

        INPUT:
            none

        OUTPUT:
            zero element in the space where self takes values

        EXAMPLES:

        """

        return self.data(0).zero()

    def zero(self):
        r"""
        Returns the modular symbol all of whose values are zero (thought of in the space where self takes values)

        INPUT:
            none

        OUTPUT:
            zero modular symbol

        EXAMPLES:

        """
        v = [self.zero_elt() for i in range(0,self.ngens())]
        C = type(self)
        return C(v,self.manin())

    def compute_full_data_from_gen_data(self):
        r"""
        Computes (and stores) the values of self on all coset representatives
        from its values on our generating set.  (Useful to have this precomputed
        before doing larger computations -- like hecke.)

        INPUT:
            none

        OUTPUT:
            none

        EXAMPLES:

        """
        ans = []   ## This will be the list contained all values of phi on coset reps
        ## This loop runs through all coset reps
        for m in range(self.ncoset_reps()):
            v = self.manin().coset_relations(m)  ##  This data is a list with each element of the form (c,A,j) where c is a scalar
                                                 ##  A is in Gamma_0(N) and j is an index of a generator.  The value of our modular symbol
                                                 ##  self on the m-th coset rep is given by the sum of c * self(j-th coset rep)|A as we run over
                                                 ##  all tuples in v

            t = self.zero_elt()  ## this will denote the value of self on the m-th coset rep
            ## This loops runs thru each tuple in v and adds up the appropriate values of self
            for k in range(len(v)):
                j = v[k][2]    ##  we need to examine the generator of index j (in the list of coset reps)
                A = v[k][1]    ##  This is the Gamma_0(N) matrix A as in the comment above
                c = v[k][0]    ##  This is the constant c as in the comment above
                r = self.manin().generator_indices().index(j)    ##  so the j-th coset rep corresponds to the r-th generator (in our list of gens)
                t = t + self.data(r).act_right(A).scale(c)
            ans = ans + [t]   ## the value of self on the m-th coset rep is added to our list
 
        self._full_data = ans   ## This data is now recorded in the modular symbol
    
    def eval_sl2(self,A):
        r"""
        Returns the value of self on the (unimodular) divisor corresponding to A -- i.e. on the divisor {A(0)} - {A(infty)} 

        INPUT:
            A an element of SL_2(Z) 

        OUTPUT:
            The value of self on the divisor corresponding to A -- i.e. on the divisor {A(0)} - {A(infty)} 

        EXAMPLES:

        """
        a = A[0,0]
        b = A[0,1]
        c = A[1,0]
        d = A[1,1]
        i = self.manin().P1().index(c,d)   ##  Finds the index in the SAGE P1 of the bottom row of A
        m = self.manin().P1_to_mats(i)     ##  Converts this index to the index in our list of coset reps
        B = self.manin().coset_reps(m)     ##  B is the corresponding coset rep -- so B and A are Gamma_0(N) equivalent
        C = invert(a,b,c,d)      ##  C = A^(-1)
        gaminv=B*C               ##  So A = gam B  (where gaminv = gam^(-1))
        ## Checks if the value of self on all coset reps is already precomputed
        if self.full_data()!=0: 
            return self.full_data(m).act_right(gaminv)   ##  Here we have self([A]) = self([gam B]) = self([B])|gam^(-1) 
                                                         ##     = self(m-th coset rep)|gamivn
        else:
            ##  We need to compute the value of self on the j-th coset rep from the values of self on our generators
            v = self.manin().coset_relations(m)    ##  This data is a list with each element of the form (c,A,j) where c is a scalar
                                                   ##  A is in Gamma_0(N) and j is an index of a generator.  The value of our modular symbol
                                                   ##  self on the m-th coset rep is given by the sum of c * self(j-th coset rep)|A as we run over
                                                   ##  all tuples in v

            t = self.zero_elt()    ## this will denote the value of self on the m-th coset rep acted on by gaminv
            ## This loops runs thru each tuple in v and adds up the appropriate values of self
            for k in range(len(v)):
                j = v[k][2]    ##  we need to examine the generator of index j (in the list of coset reps)
                A = v[k][1]    ##  This is the Gamma_0(N) matrix A as in the comment above
                c = v[k][0]    ##  This is the constant c as in the comment above
                r = self.manin().generator_indices().index(j)    ##  so the j-th coset rep corresponds to the r-th generator (in our list of gens)
                t = t + self.data(r).act_right(A*gaminv).scale(c)
            return t

    def eval(self,A):
        r"""
        Returns the value of self on the divisor corresponding to A -- i.e. on the divisor {A(0)} - {A(infty)} 

        INPUT:
            A = any 2x2 integral matrix  (not necessarily unimodular)

        OUTPUT:
            The value of self on the divisor corresponding to A -- i.e. on the divisor {A(0)} - {A(infty)} 

        EXAMPLES:

        """
        a = A[0,0]
        b = A[0,1]
        c = A[1,0]
        d = A[1,1]
        v1 = unimod_matrices_to_infty(b,d)   ## Creates a list of unimodular matrices whose divisors add up to {b/d} - {infty}
        v2 = unimod_matrices_to_infty(a,c)   ## Creates a list of unimodular matrices whose divisors add up to {a/c} - {infty}
        ans = self.zero_elt()   ## This will denote the value of self on A
        ## This loops compute self({b/d}-{infty}) by adding up the values of self on elements of v1
        for j in range(0,len(v1)):
            ans = ans + self.eval_sl2(v1[j])

        ## This loops subtracts away the value self({a/c}-{infty}) from ans by subtracting away the values of self on elements of v2
        ## and so in the end ans becomes self({b/d}-{a/c}) = self({A(0)} - {A(infty)}
        for j in range(0,len(v2)):
            ans = ans - self.eval_sl2(v2[j])
        return ans
 
    def act_right(self,gamma):
        r"""
        Returns self | gamma.
        
        This action is defined by (self | gamma)(D) = self(gamma D)|gamma

    	For this to work gamma must normalize Gamma_0(N) and be able to act on the values of self.
    	However, it can also be used to define Hecke operators.  Even if each individual self | gamma is not 
    	really defined on Gamma_0(N), the sum over acting by the appropriate double coset reps will be defined
    	over Gamma_0(N)

        INPUT:
            gamma = 2 x 2 matrix which acts on the values of self

        OUTPUT:
	    self | gamma
	    
        EXAMPLES:

        """
        v = []  ## This will be the data defining self | gamma
        ##  This loop rungs over all generators
        for j in range(0,self.ngens()):
            rj = self.manin().generator_indices(j)   ## rj is the index of the coset rep corresponding to j-th gen
            ## The value self(gamma*(rj-th coset rep))| gamma is added to our list
    	    v = v + [self.eval(gamma*self.manin().coset_reps(rj)).act_right(gamma)]  
        
        C = type(self)        
        ans = C(v,self.manin())
        ans.normalize()
        
        return ans
    
    def plus_part(self):
        r"""
        Returns the plus part of self -- i.e. self + self | [1,0,0,-1].
        
        Note that we haven't divided by 2.  Is this a problem?
        

        INPUT:
            none

        OUTPUT:
	    self | [1,0,0,-1] + self
	    
        EXAMPLES:
        """
        return self.act_right(Matrix(2,2,[1,0,0,-1])) + self

    def minus_part(self):
        r"""
        Returns the minus part of self -- i.e. self | [1,0,0,-1] - self
        
        Note that we haven't divided by 2.  Is this a problem?
        

        INPUT:
            none

        OUTPUT:
	    self | [1,0,0,-1] - self
	    
	    
        EXAMPLES:
        """
        return self.act_right(Matrix(2,2,[1,0,0,-1])) - self

    def normalize_full_data(self):
        r"""
        Normalizes the values of self on all coset reps (if they are already computed)
        
        INPUT:
            none

        OUTPUT:
	    none
	    
	    
        EXAMPLES:
        """
        if (self.full_data() != 0):
            for j in range(self.ncoset_reps()):
                self._full_data[j].normalize()

    def hecke_from_defn(self,ell):
        r"""
    	Computes self | T_ell directly from the definition of acting by double coset reps
	
        That is, it computes sum_a self | [1,a,0,ell] + self | [ell,0,0,1] where the last
        term occurs only if the level is prime to ell.
        
        INPUT:
            ell = prime

        OUTPUT:
    	    self | T_ell
	    
	    
        EXAMPLES:
        """
        ## If needed, precomputes the value of self on all coset reps
        if self.full_data() == 0:
            self.compute_full_data_from_gen_data()
            
        psi = self.zero()
        for a in range(0,ell): 
            psi = psi + self.act_right(Matrix(ZZ,[[1,a],[0,ell]]))
        if self.level()%ell<>0:
            psi = psi + self.act_right(Matrix(ZZ,[[ell,0],[0,1]]))
        
        psi.normalize()
        
        return psi

    def hecke(self,ell):
        r"""
        Returns self | T_ell by making use of the precomputations in prep_hecke
        
        INPUT:
            ell = prime

        OUTPUT:
            self | T_ell
            
        EXAMPLES:

        """         
        ## The values of self on all coset reps are computed and normalized if this hasn't been done yet.
        if self.full_data()==0:
            self.compute_full_data_from_gen_data()
            self.normalize_full_data()
            
        psi = self.zero()   ## psi will denote self | T_ell
        v = prep_hecke(ell,self.manin())  ## v is a long list of lists of lists with the property that the value
                                          ## of self | T_ell on the m-th generator is given by
                                          ## 
                                          ## sum_j sum_r self(j-th coset rep) | v[m][j][r]
                                          ## 
                                          ## where j runs thru all coset reps and r runs thru all entries of v[m][j]
        ## This loop computes (self | T_ell)(m-th generator)
        for m in range(self.ngens()):
            for j in range(self.ncoset_reps()):
                for r in range(len(v[m][j])):
                    psi._data[m] = psi.data(m) + self.full_data(j).act_right(v[m][j][r])
        
        psi.normalize()
        
        return psi

##  These two procedures are used to check that the symbol really does satisfy the Manin relations loop (for debugging)
    def grab_relations(self):
        v=[]
        for r in range(self.ngens()):
            for j in range(self.ncoset_reps()):
                R=self.manin().coset_relations[j]
                if (len(R)==1) and (R[0][2]==self.manin().generator_indices(r)):
                    if R[0][0]<>-1 or R[0][1]<>Id:
                        v=v+[R]
        return v

    def check_loop(self):
        list=self.grab_relations()
        t=self.zero_elt()
        for j in range(2,len(list)):
            R=list[j]
            index=R[0][2]
            rj=self.manin().generator_indices().index(index)
            t=t+self.data(rj).act_right(R[0][1]).scale(R[0][0])
        return self.data(0)-self.data(0).act_right(Matrix(2,2,[1,1,0,1]))+t

