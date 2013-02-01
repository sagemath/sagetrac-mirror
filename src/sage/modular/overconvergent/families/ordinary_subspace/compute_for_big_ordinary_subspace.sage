# Works
def pwvaluation(a,p,w):
	return pvaluate(a,p)+wvaluate(a,w)

# Works
def pvaluate(a,p):
	return 0.valuation(2) if a == 0 else min([ZZ(c).valuation(ZZ(p)) for c in PolynomialRing(ZZ,'w')(a).list()])

def pwf_valuation(f,p,w):
	return 0.valuation(2) if f == 0 else min([pwvaluation(c,p,w) for c in f.list()])

# Works
def wvaluate(b,w):
	#print b.list()
	for i in range(0,len(b.list())):
		if b.list()[i] != 0:
			return i
	return -1

#def wvaluate_old(a,w):
#	return 0 if a.parent() != w.parent() else a.valuation(w)

# Somehow sage in Z/27[w]/w^3 seems to think that 1/(1+3w) = 1-w+w^2
# This function computes 1/a in Z/p^m[w]/w^n for unit a (or raises 
# ValueError if passed a non-unit

# Appears to work on units

def pwn_invert(p,m,w,n,a):
	print p,m,w,n,a
	a = w.parent()(a)
	#print("inverting {0}".format(a))
	R.<ww> = PolynomialRing(Zmod(p^m),'x').quotient_ring(x^n)
	f = a.parent().hom([ww],R)
	b = f(a)
	if(b.list()[0].valuation(p) > 0): 
		raise ValueError("{0} not invertible".format(a))
	a0 = Zmod(p^n)(1)/(Zmod(p^n)(b.list()[0]))
	b = -((b*a0)-1)
	ans = 1
	while(b != 0):
		ans += b
		b *= b
		#print b,ans
	ans *= a0
	return f_preimage_stupid(ans,ww,w)

# Given f a polynomial in variable www, return a polynomial 
# with the same coefficients except in the variable w

# Works

def f_preimage_stupid(f,www,w):
	ans = 0;
	f = www.parent()(f)
	for i in range(0,len(f.list())):
		ans += w^i*ZZ(f.list()[i])
	return ans

# Solve ax=b mod p^m,w^n if possible (raise ValueError if not)

# Works when it can reduce a to a unit, but cannot solve e.g. (w+3)x=w^2+3w

def pwn_solve(p,m,w,n,a,b):
	#print("solving ({0})x = {1} mod {2},{3}".format(a,b,p^m,w^n))
	a = w.parent()(a)
	b = w.parent()(b)
	if(pwvaluation(a,p,w) > pwvaluation(b,p,w)): 
		raise ValueError("({0})x = {1} mod {2},{3} not solvable".format(a,b,p^m,w^n))
	v = pvaluate(a,p)
	vv = wvaluate(a,w)
	#print v,vv
	S.<ww> = PolynomialRing(Zmod(p^n))
	I = S.ideal(ww^n)
	R.<www> = S.quotient_ring(I)
	#print w.parent().gens()
	f = S.hom([www],R)
	#print a,b,v,vv,b/(p^v*w^vv),a/(p^v*w^vv)
	return f_preimage_stupid(f(S(b/(p^v*w^vv)))*f(pwn_invert(p,m,ww,n,S(a/(p^v*w^vv)))),www,w)

# Given r a root mod (p,w) of the given f in (Z/p^m[w]/w^n)[x], check if r
# can be lifted to an honest root of a in this ring, and if so, return 
# that.  Otherwise, raise a ValueError

# Appears to work

def pwn_hensel(p,m,w,n,f,r):
	pwval = 2*pwvaluation(f.derivative()(r),p,w)+1
	print f.derivative()(r)
	print f(r)
	print pwval
	print pwvaluation(f(r),p,w)
	
	if(pwvaluation(f(r),p,w) < pwval):
		raise ValueError("r cannot be lifted to a root of f")

	print r,f,f.derivative()
	ans = r
	while(pwval < n+m):
		ans += pwn_solve(p,n,w,m,f.derivative()(ans),-f(ans))
		pwval += 1
		print ans
	
	return ans

# Swap columns i,j of the matrix represented by the 2D array L

# Works 

def pwn_swap_cols(L,i,j):
	for ii in range(len(L)):
	    temp = L[ii][j]
	    L[ii][j] = L[ii][i]
	    L[ii][i] = temp


# Swap rows i,j of the matrix represented by the 2D array L

# Works 

def pwn_swap_rows(L,i,j):
	if(i != j):
		A = L[i]
		L[i] = L[j]
		L[j] = A

# Change L by L[i] -= x*L[j]

# Works 

def pwn_row_op(p,m,w,n,L,i,j,x):
	#print "ROW_OP",p,m,w,n,L,i,j,x
	if(i == j): return
	for a in range(0,len(L[i])):
		L[i][a] -= x*L[j][a]
		L[i][a] = w.parent(PolynomialRing(ZZ,'w')(L[i][a])%(p^m))%(w^n)
		# (w.parent()(L[i][a])%(p^m))%(w^n)
	#print "RESULT",L

# Return the n x n 2D array representing the identity	

# Works 

def pwn_id(n):
	E = [[0 for x in range(n)] for y in range(n)]
	for i in range(n): E[i][i] = 1
	return E

# row reduce a matrix with entries in ZZ[w]/(p^m,w^n)

# Appears to work?

def pwn_row_reduce(p,m,w,n,L,s,E):
	#print "ROW_REDUCE",p,m,w,n,s
	if(s >= len(L) or s >= len(L[0])): return E

	# Put minimal $p$-valuation guy in top left
	mv = 9999
	rmv = s
	cmv = s
	for i in range(s,len(L)):
		for j in range(s,len(L[i])):
			v = pwvaluation(L[i][j],p,w)
			if(v < mv):
				mv = v
				rmv = i
				cmv = j
	#print("Swapping {0},{1} in \n{2}".format(rmv,cmv,L))
	#print("Swapping rows {0},{1}".format(s,rmv))
	pwn_swap_rows(L,s,rmv)
	pwn_swap_rows(E,s,rmv)
	pwn_swap_cols(L,s,cmv)
	#print(L)
	#print "E"
	#print E
	
	# Clean out first column
	a = L[s][s]
	if(a == 0): 
		return E
	for i in range(s+1,len(L)):
		x = pwn_solve(p,m,w,n,a,L[i][s])
		#print "row op row({0}) -= ({2})*row({1})".format(i,s,x)
		#print "E:"
		#print E
		pwn_row_op(p,m,w,n,L,i,s,x)
		pwn_row_op(p,m,w,n,E,i,s,x)
		#print "E:"
		#print E
	return pwn_row_reduce(p,m,w,n,L,s+1,E)

# (The difference between .mod() and % caused this to fail before--
# (w^3).mod(w^3) is not zero, apparently)
# Appears to work 

def is_all_zero(p,m,w,n,L):
	#print w.parent(),p,m,w,n,L
	for x in L:
		#print w.parent(PolynomialRing(ZZ,'w')(x)%(p^m))%(w^n)
		if(w.parent(PolynomialRing(ZZ,'w')(x)%(p^m))%(w^n) != 0):
			return false
	#print "ZERO"
	return true

# Probably works, subject to pwn_row_reduce

def get_rank(p,m,w,n,L):
	print "BLAH"
	print L
	E = pwn_row_reduce(p,m,w,n,L,0,pwn_id(len(L)))
	print L
	r = len(L)
	while(r > 0):
		print r,p^m,w^n,L[r-1]
		if(not is_all_zero(p,m,w,n,L[r-1])):
			print "DONE",r
			break
		r -= 1
	print "r=",r
	return [r,E]

# Multiply everything up so that the precisions match.  
# At the start, know things with precisions: 
# (p^M,p^(M-1),...,p)
# so multiply entries by 1,p,p^2,...,p^(M-1) respectively

# Works

def OMS_to_list(phi):
	p = phi.p()
	L = []
	M = phi.data[0].num_moments()
	for i in range(0,len(phi.data)):
		for j in range(0,phi.data[i].num_moments()):
			L.append(p^j*phi.data[i].moment(j))
	return L

# Return the power of p modulo which we know the jth moment in Fil^M

def fil_pow(p,M,j):
	#return j
	return ceil(((M-j)*(p-2))/(p-1))

# This function takes a family given like [(m1,m2m...),(n1,n2,...)]
# and converts it to a single list like [m1,m2,...,n1,n2,...], except
# that the precisions are all different on these, so we need to
# multiply everything up so that the precisions match.  At the start,
# know things with precisions: 
#
# (p^(ceil(M c_p)),p^(ceil((M-1)c_p)),...,p) 
#
# so multiply jth entry by p^(ceil(M c_p)-ceil((M-j)c_p)) to normalize
# all precisions

# Unsure if the normalisation to force all moments into the same precision is correct

def OMS_fam_to_list(phi):
	p = phi.p()
	L = []
	M = phi.data[0].num_moments()
	for i in range(0,len(phi.data)):
		for j in range(0,phi.data[i].num_moments()):
			#L.append(p^j*phi.data[i].moment(j))
			L.append(p^(fil_pow(p,M,0) - fil_pow(p,M,j))*phi.data[i].moment(j))
	return L

# Test if phi is in the span of the families in L up to precision M
# Returns [b,E] where b is a boolean answering the `is it in the
# span?'  question and E is the elementary matrix row-reducing [L,phi]
# (Hence if phi is in the span E[-1] will provide a linear dependence
# on [L,phi])

# Works, I think, subject to dependencies OMS_fam_to_list, fil_pow,
# and get_rank

def is_in_span(phi,L,p,N,k,M,w):
	print "checking if {0} is in span {1}".format(phi,L)
	if(len(L) == 0): 
		return [false,None]
	A = []
	for x in L:
		A.append(OMS_fam_to_list(x))
	print "lenA",len(A)
	r = get_rank(p,fil_pow(p,M,0),w,M,A)[0]
	#print(A)
	A = []
	for x in L:
		A.append(OMS_fam_to_list(x))
	A.append(OMS_fam_to_list(phi))
	#print(A)
	print "lenA",len(A)
	s = get_rank(p,fil_pow(p,M,0),w,M,A)
	print("ranks")
	print(r)
	print(s[0])
	print r == s[0],s[1]
	return [r == s[0],s[1]]

# Arguments: 
#    L,p,N,k -- a list of forms already known to be in the ordinary subapce for 
#               level N, prime p, weight k
#    M       -- number of moments to store for the distributions forming the ordinary
#               subspace
#    w       -- the variable used in the power series defining the family
#    d       -- the actual dimension of the ordinary subspace we're computing
#
# Returns: a list of two elements [[phi, Up(phi), Up^2(phi), ...,
# Up^a(phi)],E] where the first element is a cyclic subspace of the
# ordinary subspace where phi is not in the span of the given elements
# in L, and where E is the elementary matrix row-reducing the matrix
# [phi, Up(phi), Up^2(phi), ..., Up^a(phi), Up^(a+1)(phi)] (Hence the
# E[-1] gives the coefficients of the characteristic polynomial of Up
# on the returned cyclic subspace of the ordinary subspace.)

# Works subject to dependencies: is_in_span

def build_ordinary_subspace_by_cyclics(L,p,N,k,M,w,d):
	cs = []
	ans = [cs,0]
	phi=random_OMS_fam(p,N,M,M,k,w)
	for i in range(0,M):
		phi = phi.hecke(p)
	phi = phi-phi.hecke(p)
	if(is_in_span(phi,L,p,N,k,M,w)[0]): return []
	cs.append(phi)
	for i in range(0,d):
		print "i={0}".format(i)
		phi = phi.hecke(p)
		print(phi)
		s = is_in_span(phi,cs,p,N,k,M,w)
		if(s[0]): 
			print s[1]
			ans[1] = s[1]
			break
		elif(i+1 == d):
			raise Exception("SOMETHING IS WRONG: Getting bigger ordinary subspace than we should")
		cs.append(phi)
	return ans

# Arguments: 
#    p,N,k -- compute the ordinary subapce for modular forms of level N, prime p, 
#             weight k
#    M     -- number of moments to store for the distributions forming the ordinary
#             subspace
# 
# Returns: A 2-element list [L,m,f] where L is a list contianing a
#    basis of ordimary families, m is the matrix describing the action
#    of Up given as a list of blocks, and f is the characteristic
#    polynomial of Up on this space
#
# Repeatedly use build_ordinary_subspace_by_cyclics to build the full
# ordinary subspace

# Possible precision issues

def compute_ordinary_subspace(p,N,k,M,w):
	d = dimension_of_cuspidal_ordinary_subspace(N,p,k)
	
	ordinary_subspace = []
	ordinary_subspace_by_cyclics=[]
	while(sum(map(len,ordinary_subspace_by_cyclics)) < d):
		temp = build_ordinary_subspace_by_cyclics(ordinary_subspace,p,N,k,M,w,d)
		if(len(temp[0]) > 0):
			ordinary_subspace_by_cyclics.append(temp)
			ordinary_subspace += temp[0]

	print("Ordinary subspace: {0}".format(ordinary_subspace))

	S = PolynomialRing(ZZ,'w')
	Up_char_poly = PolynomialRing(S,'X')(1)
	Up_matrix = []
	for i in ordinary_subspace_by_cyclics:
		Up_char_poly *= PolynomialRing(S,'X')(i[1][-1])
		Up_block = [[0 for j in range(len(i[1]))] for j in range(len(i[1]))]
		for j in range(len(i[1])-1):
			Up_block[j][j+1] = 1
		Up_block[-1] = i[1][-1]
		Up_matrix.append(Up_block)

	print("Char poly: {0}".format(Up_char_poly))	

	return [ordinary_subspace,Up_matrix,Up_char_poly]

def pwn_modpw(f,p,w):
	L = []
	print f.list()
	for i in range(len(f.list())):
		L.append(PolynomialRing(Zmod(p),'w')((f.list()[i]))(0))
	print L
	return PolynomialRing(Zmod(p),'X')(L)

# TODO
# Maybe broken?  More likely seems it's being fed nonsense

def find_basis_of_eigenforms(p,N,k,M,w):
	O = compute_ordinary_subspace(p,N,k,M,w)
	eigenbasis = []
	U = O[1]
	f = O[2]
	Rs = pwn_modpw(f,p,w).roots(multiplicities=False)
	print "Roots of ",pwn_modpw(f,p,w)," mod p,w: ",Rs
	roots = []
	for r in Rs:
		roots.append(pwn_hensel(p,M,w,M,f,r))
	for i in range(len(roots)):
		ans = O[0]
		for j in range(len(roots)):
			if(i != j):
				ans = ans.hecke(p)-ans.scale(roots[j])
		eigenbasis.append([ans,roots[i]])
	return eigenbasis

def x011test():
	S.<w> = PolynomialRing(QQ)
	return find_basis_of_eigenforms(3,11,2-2,6,w)

def rr_test():
	R.<w> = PolynomialRing(ZZ)
	pwn_solve(3,3,w,3,3*w+1,9*w)

	R.<x> = PolynomialRing(ZZ)		
	L = [[4*x+x^2,2,3*x],[4*x+x^2+2*x^8,2-4*x,3*x],[x^7,-2,0]]
	E = pwn_row_reduce(3,3,x,3,L,0,pwn_id(len(L)))
	print(L)
	print E

