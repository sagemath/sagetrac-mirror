##  A -- a matrix
##  t -- a list with length equal to the number of columns of A
##
##  Returns whether or not t is in the span of the rows of A modulo mod
def in_span(A,t,p,n):
	mod = p^n
	B = A.insert_row(A.nrows(),t)
	M = MatrixSpace(Integers(mod),B.nrows(),B.ncols())
	B = M(B)
	B,E,r = row_reduce_mod_pn(B,p,n)
#	B = (B.echelon_form())%(mod)
#	B = (B.echelon_form())%(mod)  ## This is a silly trick to force it to clear out multiples of mod
	
	
	return r < A.rank() + 1

##  A -- a matrix
##  t -- a list with length equal to the number of columns of A
##
##  Returns a vector expressing t as a linear combination of the rows of A modulo mod -- (in_span(A,t,mod) should be true)
def linear_combo(A,t,p,n):
	mod = p^n	
#	print "In linear_combo with",A,t,mod
	B = A.insert_row(A.nrows(),t)%(mod)
	M = MatrixSpace(Integers(mod),B.nrows(),B.ncols())
	B = M(B)
	z1,z2,r = row_reduce_mod_pn(B,p,n)
	rel = z2.row(-1)  ## This is the bottom row which expresses the relation

	assert (Matrix(rel)*B)%mod == 0, "relation not holding in linear_combo"

	c = rel[len(rel)-1]^(-1)
	rel = (c * rel)
	rel = (-rel)

	rel = Sequence(rel)
	rel.pop()


	return rel

##  N -- tame level
##  p -- prime
##  k -- weight
##  M -- accuracy we are working with (i.e. p^M)
##  d -- dimension of the ordinary space (which is computed beforehand via Hida theory)
##  sign -- 1 or -1 (indicates which subspace of modular symbols we are using
##
##  OUTPUT -- a list of OMSs of tame level Gamma_1(N), character chi, weight k which form
##  a basis of of the space of such OMSs modulo p^M.  
##
##  The method is produce a random symbol and project to the ordinary subspace.  Then find another one
##  and see if the two together span a free Z/p^M module of rank 2.  (This is done by looking at the
##  the vector of total measures.)  Then a third random symbol is formed and check if these three
##  symbols span a rank 3 free Z/p^M-module.  If not, this vector is thrown away.  And repeat.
def form_basis(N,p,k,M,chi,d,sign):
	Phis = []
	v = []

	total = 0
	while total < d:
		done = false
		while not done:
			Phi = random_ordinary_OMS(N,p,k,M,chi)
			if sign == 1:
				Phi = Phi.plus_part()
			else:
				Phi = Phi.minus_part()
			v = v + [Phi.vector_of_total_measures()]
			Phis = Phis + [Phi];
			B = Matrix(v)
			Mat = MatrixSpace(Integers(p^M),B.nrows(),B.ncols())
			B = Mat(B)
			B,E,r = row_reduce_mod_pn(B,p,M)
			print B
			print "free rank: ",r

			if r == B.nrows():
				print "Keeping it"
				total = total + 1
				done = True
			else:
				print "Failed"
				done = False
				Phis.pop()
				v.pop()

	return Phis

##  N -- tame level
##  p -- prime
##  k -- weight
##  M -- accuracy we are working with (i.e. p^M)
##  d -- dimension of the ordinary space (which is computed beforehand via Hida theory)
##  sign -- 1 or -1 (indicates which subspace of modular symbols we are using
##
##  OUTPUT -- a list of OMSs of tame level Gamma_1(N), character chi, weight k which form
##  a basis of of the space of such OMSs modulo p^M.  
##
##  The method here is an improvement over form_basis in that rather than picking a new random
##  symbol each time, the U_p-span of each symbol is considered until they no longer span a free module.
def form_basis2(N,p,k,M,chi,d,sign):
	Phis = []
	v = []
	total = 0

	new_seed_needed = True
	while (total < d):
		if new_seed_needed:
			print "NEW SEED"
			Phi = random_ordinary_OMS(N,p,k,M,chi)
			new_seed_needed = False
		else:
			Phi = Phi.hecke(p)
		if sign == 1:
			Phi = Phi.plus_part()
		else:
			Phi = Phi.minus_part()
		v = v + [Phi.vector_of_total_measures()]
		Phis = Phis + [Phi];
		B = Matrix(v)
		Mat = MatrixSpace(Integers(p^M),B.nrows(),B.ncols())
		B = Mat(B)
		B,E,r = row_reduce_mod_pn(B,p,M)
		print B
		print "free rank: ",r
		if r == B.nrows():
			print "Keeping it"
			total = total + 1
			done = True
		else:
			print "Failed"
			done = False
			Phis.pop()
			v.pop()
			new_seed_needed = True

	return Phis
	

## Phis -- a list of OMSs which generate the ordinary subspace
##
## Computes the matrix of the q-th Hecke operator acting on the span of these OMSs
def hecke_matrix(Phis,q):
	p = Phis[0].p()
	M = Phis[0].num_moments()
	A = Matrix([Phis[r].vector_of_total_measures() for r in range(len(Phis))])

	Tq = []
	for r in range(len(Phis)):
		print r
		t = Phis[r].hecke(q).vector_of_total_measures()
		Tq = Tq + [linear_combo(A,t,p,M)]
		
	Tq = Matrix(Tq)
	Mat = MatrixSpace(ZZ,Tq.nrows(),Tq.ncols())
	return Mat(Tq)
		
##  chi - Dirichlet character
##  p - prime
##  k - weight (really the weight minus 2)
##  sign - 1 or -1
##
##  Returns the dimension of the ordinary subspace of overconvergent (p-adic) modular symbols of weight k, 
##  character chi, and sign sign.  (This procedure uses Hida theory to do this computation.)
##
##  SHOULD BE REPROGRAMMED TO JUST USE HIDA THEORY TO REDUCE TO WEIGHT 2 -- ALREADY DONE IN OTHER WEIGHT 1 CODE
##  ALSO SHOULD INCLUDE P=2.
def dimension_of_ordinary_subspace(chi,p,k,sign):
	r = k % (p-1) + 2
	M = ModularSymbols(chi,r,sign,GF(p))
	
	hecke_poly = M.hecke_polynomial(p)
	R = hecke_poly.parent()
	x = R.gen()
	return hecke_poly.degree() - hecke_poly.ord(x)

##  N - integer (the tame level; assumes that p does not divide N)
##	or a mod p Dirichlet character (the Nebentypus; assumes of level
#	not divisible by p^2, i.e. of first kind)
##  p - prime
##  k - weight greater than or equal to 2
##
##  Returns the dimension of the S_k(Gamma_0(Np))^ord if N is an integer
##  or S_k(Gamma_1(Np), chi)^ord if chi = N is a Dirichlet character
##
##  SHOULD BE REPROGRAMMED TO INCLUDE P=2
def dimension_of_cuspidal_ordinary_subspace(N, p, k):
	r = ((k - 2) % (p - 1))
	if N in ZZ:
		if r == 0:
			M = ModularSymbols(Gamma0(N * p), 2, 1, GF(p)).cuspidal_subspace()
		else:
			DG = DirichletGroup(N * p, GF(p))
			chi = [GF(p)(u) ** r for u in DG.unit_gens()]	#mod p Teichmuller^r
			chi = DG(chi)
			M = ModularSymbols(chi, 2, 1, GF(p)).cuspidal_subspace()
	else:
		n = N
		if not p.divides(N.level()):
			n = N.extend(p * N.level())
		if r == 0:
			M = ModularSymbols(n, 2, 1, GF(p)).cuspidal_subspace()
		else:
			DG = DirichletGroup(n.level(), GF(p))
			chi = [GF(p)(u) ** r for u in DG.unit_gens()]
			chi = n * DG(chi)
			M = ModularSymbols(chi, 2, 1, GF(p)).cuspidal_subspace()
	
	hecke_poly = M.hecke_polynomial(p)
	R = hecke_poly.parent()
	x = R.gen()
	return hecke_poly.degree() - hecke_poly.ord(x)

##  N - integer
##  p - prime
##  k - weight greater than or equal to 2
##  sign - 1 or -1
##
##  Returns the dimension of the S_k(Gamma_0(Np))^ord
##
##  SHOULD BE REPROGRAMMED TO JUST USE HIDA THEORY TO REDUCE TO WEIGHT 2 AND INCLUDE P=2
def dimension_of_cuspidal_ordinary_subspace_old(N,p,k):
	r = ((k-2) % (p-1)) + 2
	if r == 2:
		M = ModularSymbols(N*p,r,1,GF(p)).cuspidal_subspace()
	else:
		M = ModularSymbols(N,r,1,GF(p)).cuspidal_subspace()
	
	hecke_poly = M.hecke_polynomial(p)
	R = hecke_poly.parent()
	x = R.gen()
	return hecke_poly.degree() - hecke_poly.ord(x)
		
##  B -- matrix
##  p -- prime
##  n -- integer
##
##  OUTPUT - Three bits data are returned -- (A,E,r)
##
##  This procedure attempts to row reduce B modulo p^n.  That is, it is always true that E*B = A.  However,
##  I never finished programming the case when one of the elementary factors is not 1.  I think the function
##  just ends the row reduction there.  The variable r is the rank of the matrix when all of the elementary
##  factors are 1.
def row_reduce_mod_pn(B,p,n):
	A = copy(B)
	rows = A.nrows()
	cols = A.ncols()
	M = MatrixSpace(Integers(p^n),rows,cols)
	A = M(A)
	curr_col = 0
	curr_row = 0
	E = Matrix(Integers(p^n),rows,rows,1)
	while (curr_row < rows) and (curr_col < cols):
		#print (curr_row,curr_col)
		v = [A[r][curr_col] for r in range(curr_row,rows)]
		#print v
		vals = [ZZ(v[j]).valuation(p)  for j in range(len(v))]
		#print vals
		if min(vals) == 0:
			least_ind = curr_row
			while vals[least_ind-curr_row] > 0:
				least_ind = least_ind + 1
			#print "least_ind = ",least_ind
			A.swap_rows(curr_row,least_ind)
			F = Matrix(Integers(p^n),rows,rows,1)
			F.swap_rows(curr_row,least_ind)
			E = F * E	
			for r in range(rows):
				if r != curr_row:
					#print A[curr_row,curr_col]
					c = -A[r,curr_col]/A[curr_row,curr_col]
					A.add_multiple_of_row(r,curr_row,c)
					F = Matrix(Integers(p^n),rows,rows,1)
					F[r,curr_row] = c
					E = F * E
			curr_row = curr_row + 1
		curr_col = curr_col + 1
		#print A
		#print "----"
		#print E*A
		#print "-------------"
	return A,E,curr_row
	
		