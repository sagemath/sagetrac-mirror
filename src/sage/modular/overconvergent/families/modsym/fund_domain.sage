from sage.matrix.matrix_integer_2x2 import Matrix_integer_2x2 as mi2x2
import sage.matrix.all as matrix

M2Z = matrix.MatrixSpace(ZZ,2)

Id = mi2x2(M2Z, [1, 0, 0, 1], copy = False, coerce = False)

def p1_equiv(v,w,N):
	"""determines whether the row vectors v and w are equal in P^1(Z/NZ)"""
	equiv=false
	for a in range(N):
		if gcd(a,N)==1:
			if ((v[0]*a-w[0])%N==0) and ((v[1]*a-w[1])%N==0):
				equiv=true
	return equiv
		
def bottom_row(A):
	"""returns the bottom row of the 2x2 matrix A"""
	return [A[1,0],A[1,1]]

def in_gamma0(A,N):
	"""determines if the 2x2 matrix is in Gamma_0(N)"""
	return A[1,0]%N==0

def p1_flip(v):
	"""takes the vector [a,b] and returns [-b,a] --- related to 2 term manin relation"""
	return [-v[1],v[0]]

def form_list_of_cusps(N):
	assert N>1, "Error in form_list_of_cusps: level should be > 1"
	P=P1List(N)
	sP=len(P.list())
	matrices=[mi2x2(M2Z, [1, 0, 0, 1], copy = False, coerce = False)]
	C=[-1,"?",0]
	tau=mi2x2(M2Z, [0, -1, 1, -1], copy = False, coerce = False)
	full_domain=false

	v=[0 for r in range(sP)]
	ans=P.index(0,1)
	v[ans]=1
	ans=P.index(1,-1)
	v[ans]=2
	ans=P.index(-1,0)
	v[ans]=3
	position=4

	s=0
	while (not full_domain):
		s=s+1		
		full_domain=true
	
#this loop runs through the current set of cusps and checks to see
#if more cusps should be added

		for r in range(floor((len(C)-1)/2)):
			if C[2*r+1] == "?":
				cusp1=C[2*r]
				cusp2=C[2*r+2]
				a1=cusp1.numerator()
				b1=cusp1.denominator()
				a2=cusp2.numerator()
				b2=cusp2.denominator()
				#a=a1+a2
				#b=b1+b2
				gam = mi2x2(M2Z, [a2, a1, b2, b1], copy = False, coerce = False)
				
			#this is the point where it is determined whether
			#or not the adjacent triangle should be added

				pos=P.index(b2,b1)
                                if v[pos] == 0:
                                        v[pos]=position
                                        ans=P.index(b1,-(b1+b2))
                                        v[ans]=position+1
                                        ans=P.index(-(b1+b2),b2)
                                        v[ans]=position+2
                                        position=position+3
					needed=true
				else:
					needed = false

				if needed:
					if (gam*tau*gam^(-1))[1][0] % N == 0:
						matrices=matrices+[gam]
						C[2*r+1]="t"
					else:
						C[2*r+1]="i"
						matrices=matrices+[gam]
						full_domain=false
				else:
					C[2*r+1]="x"
		r = 0
		while 2*r+1 < len(C):
			if C[2*r+1] == "i":
				C[2*r+1]="?"
				cusp1=C[2*r]
				cusp2=C[2*r+2]
				a1=cusp1.numerator()
				b1=cusp1.denominator()
				a2=cusp2.numerator()
				b2=cusp2.denominator()
				a=a1+a2
				b=b1+b2
				C.insert(2*r+2,a/b)
				C.insert(2*r+3,"?")
			r = r + 1
#	print "Finished level",s,"at",Cputime(t)

	return C



def weird_cusp_list_to_mats(C):
	v=[Rational(C[2*r]) for r in range((len(C)+1)/2)]
	v.reverse()
	mats = [Id, mi2x2(M2Z, [1,1,-1,0], copy = False, coerce = False)]
	for j in range(len(v)-1):
		a=v[j].numerator()
		b=v[j+1].numerator()
		c=v[j].denominator()
		d=v[j+1].denominator()
		mats.append(mi2x2(M2Z, [a,b,c,d], copy = False, coerce = False))
	return mats
		
def unimod(r1,r2):
	a=r1.numerator()
	b=r2.numerator()
	c=r1.denominator()
	d=r2.denominator()
	return (a*d-b*c)^2==1

def unimod_to_matrices(r1,r2):
	a=r1.numerator()
	b=r2.numerator()
	c=r1.denominator()
	d=r2.denominator()
	if (a*d-b*c)==1:
		return mi2x2(M2Z, [a,b,c,d], copy = False, coerce = False),mi2x2(M2Z, [-b,a,-d,c], copy = False, coerce = False)
	else:
		return mi2x2(M2Z, [-a,b,-c,d], copy = False, coerce = False),mi2x2(M2Z, [b,a,d,c], copy = False, coerce = False)


def solve_manin_relations(N):
	C=form_list_of_cusps(N)
	P=P1List(N)
	cusps=[Rational(C[2*r]) for r in range((len(C)+1)/2)]
	cusps.reverse()
	mats=weird_cusp_list_to_mats(C)
	sig=mi2x2(M2Z, [0,1,-1,0], copy = False, coerce = False)
	tau=mi2x2(M2Z, [0,-1,1,-1], copy = False, coerce = False)
	p1s=[bottom_row(mats[j]) for j in range(len(mats))]
	gens=[]
	twotor=[]
	twotorrels=[]
	threetor=[]
	threetorrels=[]	
#this loop is choosing our generators from the edges of the fundamental domain
	rels = [0 for i in range(0,len(mats))]
	boundary_checked=[0 for i in range(0,len(mats))]
	glue_data=[0 for i in range(0,len(mats))]
	for r in range(0,len(mats)):	
		if boundary_checked[r]==0:
			if P.index(p1s[r][0],p1s[r][1]) == P.index(-p1s[r][1],p1s[r][0]):
				#two-torsion!  mats[r] is made a generator too
				rels[r]=[(1,Id,r)]
				gens=gens+[r]
				twotor=twotor+[r]
				gam=mats[r]*sig*mats[r]^(-1)
				twotorrels=twotorrels+[gam]
		#a matrix gam appearing in twotorrels means that 1+gam kills the value on this divisor
				glue_data[r]=(r,gam)
				boundary_checked[r]=1
			else:
				c=mats[r][1,0]
				d=mats[r][1,1]
				if (c^2+d^2+c*d)%N==0:
					#three torsion!
					rels[r]=[(1,Id,r)]
					gens=gens+[r]
					threetor=threetor+[r]
					gam=mats[r]*tau*mats[r]^(-1)
					threetorrels=threetorrels+[gam]
		#a matrix gam appearing in threetorrels means that 1+gam+gam^2 kills the value on this divisor
					#GLUE DATA???
					boundary_checked[r]=1
		#need to include the three-torsion path in the other direction
					a=mats[r][0][0]
					b=mats[r][0][1]
					A=mi2x2(M2Z, [-b,a,-d,c], copy = False, coerce = False)
					mats=mats+[A]
					rels=rels+[[(-1,Id,r)]]
		#the value of phi of A is just the negative of the value of phi on mats[r]
				else:
					found=false
					s=r+1
					while (not found):
						if P.index(p1s[s][0],p1s[s][1]) == P.index(-p1s[r][1],p1s[r][0]):
	#mats[r] will now be made a generator and we need to express phi(mats[s]) in terms of phi(mats[r])
							A=mats[s]*sig 
							gam=mats[r]*A^(-1)		
							rels[r]=[(1,Id,r)]
							gens=gens+[r]
#a relation in the r-th slot of the form (1,Id,0) means that mat[r] is a basic generator with its value directly stored.
							rels[s]=[(-1,gam,r)]
#a relation in the s-th slot of the form (c,gam,r) means that phi evaluated on mats[s] equals c*phi(mats[r])|gam
							glue_data[r]=(s,gam)
							found=true
							boundary_checked[r]=1
							boundary_checked[s]=1
						s=s+1
				
#this loop expresses the relations of all of the interior paths in terms of the boundary,
	for r in range(len(cusps)-2):
		for s in range(r+2,len(cusps)):
			r1=cusps[r]
			r2=cusps[s]
			if unimod(r1,r2):
				A,B=unimod_to_matrices(r1,r2)
				mats=mats+[A,B]
				vA=[]
				vB=[]
				for j in range(s-r):
					vA=vA+rels[r+j+2]
					t=(-rels[r+j+2][0][0],rels[r+j+2][0][1],rels[r+j+2][0][2])
					vB=vB+[t]
				rels=rels+[vA,vB]
	return [mats,gens,twotor,twotorrels,threetor,threetorrels,rels,glue_data]

def form_P1_to_mat_list(mats,P):
	v=zero_vector(len(mats))
	for r in range(len(mats)):	
		v[P.index(mats[r][1,0],mats[r][1,1])]=r
	return v


class manin_relations(SageObject):
	def __init__(self,N):
		A=solve_manin_relations(N)
		self.P=P1List(N)
		self.mats=A[0]
		self.P1_to_mats=form_P1_to_mat_list(A[0],self.P)
		self.gens=A[1]
		self.twotor=A[2]
		self.twotorrels=A[3]
		self.threetor=A[4]
		self.threetorrels=A[5]
		self.rels=A[6]
		self.glue=A[7]

	############################################################################################################
	## Takes the i-th generator defining self and returns the associated matrix gamma_i -- i.e. if we are not in 	
	## two or three torsion case, we return the matrix for which 
	##
	##  phi(D_infty) | Delta = \sum phi(D_i) | (gamma_i - 1)  + two and three torsion stuff
	##
	## In the two torsion case, we return the matrix such that phi(D'_i)|(gamma_i+1) = 0 
	## In the three torsion case, we return the matrix such that phi(D''_i)|(gamma_i^2+gamma_i+1) = 0 
	############################################################################################################
	def gen_rel_mat(self,i):
		i = self.gens[i]
		if self.twotor.count(i) > 0:
			index = self.twotor.index(i)
			return self.twotorrels[index]
		elif self.threetor.count(i) > 0:
			index = self.threetor.index(i)
			return self.threetorrels[index]
		else:
			R = self.rels
			for a in range(len(R)):
				if (len(R[a]) == 1) and (R[a][0][0] == -1) and (R[a][0][2] == i):
					return R[a][0][1]