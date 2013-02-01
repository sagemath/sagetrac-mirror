from sage.structure.sage_object import SageObject

class dist(SageObject):
	def __init__(self,p,k,moments,char=None):
		"""A distribution is stored as a vector whose j-th entry is the j-th moment of the distribution.  The j-th entry is stored modulo p^(N-j) where N is the total number of moments.  (This is the accuracy that is maintained after acting by Gamma_0(p).)

Inputs: 
	p -- prime 
	k -- weight (used to determine the action of matices)
	moments	-- the list of moments given as a vector
	char -- optional Dirichlet character"""
		self.p=p
		self.weight=k
		self.moments=moments
		if char != None:
			self._char = char
		else:
		 	self._char = DirichletGroup(1,QQ).0
			
	def moment(self,n):
		"""returns the n-th moment"""
		return self.moments[n]

	def char(self):
		"""returns the character of the distribution"""
		return self._char

	def __repr__(self):
		"""Displays the moments of the distribution"""
		return repr(self.moments)

	def num_moments(self):
		"""Returns the number of moments of the distribution"""
		return len(self.moments)

	def __add__(self,right):
		"""Adds the distributions self and right"""
		assert self.weight==right.weight, "the weights are different"
		assert self.num_moments()==right.num_moments(), "the accuracies are different"
		assert self.char() == right.char(), "the characters are different"
		return dist(self.p,self.weight,self.moments+right.moments,self.char())

	#should try to figure out __lmul__
	def _lmul_(self,other):
		"""Scales self by the constant left"""
		return dist(self.p,self.weight,other*self.moments,self.char())

	def scale(self,left):
		"""Scales self by the constant left"""
		return dist(self.p,self.weight,left*self.moments,self.char())

	def __sub__(self,right):
		"""Returns self-right"""
		return self+right.scale(-1)

	def __cmp__(self,right):
		"""Returns true is both the moments and the weight agree"""
		return cmp((self.weight,self.moments,self.char()),(right.weight,right.moments,right.char()))

	def zero(self):
		"""returns a distribution with all zero moments (and the same prime and weight)"""
		return dist(self.p,self.weight,zero_vector(QQ,self.num_moments()),self.char())
	
	def is_zero(self):
		"""Returns whether all of self's moments are zero."""
		return self.moments.is_zero()

	def valuation(self):
		"""returns the highest power of p which divides all moments of the distribution"""
		p=self.p
		return min([self.moment(a).valuation(p) for a in range(self.num_moments())])

	def normalize(self):
		"""reduced modulo Fil^N -- that is the i-th moments is reduced modulo p^(N-i)"""
		p=self.p
		if self.valuation() >= 0:		
			v=vector([self.moment(i)%(p^(self.num_moments()-i)) for i in range(0,self.num_moments())])
			return dist(self.p,self.weight,v,self.char())
		else:
			return self

	def normalize_aws(self):
		"""Arizona Winter School filtration (for families) -- reduces the j-th moment modulo p^(ceil((N-j)*(p-2)/(p-1))"""
		p=self.p
		N=self.num_moments()
		assert self.valuation() >= 0, "moments not integral in normalization"
		
		v=vector([self.moment(j)%(p^(ceil((N-j)*(p-2)/(p-1)))) for j in range(0,self.num_moments())])	
		return dist(self.p,self.weight,v,self.char())

	def change_precision(self,M):
		"""only holds on to M moments"""
		assert M<=self.num_moments(),"not enough moments"

		v=[self.moment(i) for i in range(M)]
		mu=dist(self.p,self.weight,vector(v),self.char())
		return mu.normalize()

	def specialize(self):
		"""specializes to weight k -- i.e. projects to Sym^k"""
		k = self.weight
		R.<X,Y>=PolynomialRing(QQ,2)
		P = 0
		for j in range(0,k+1):
			P=P+binomial(k,j)*(-1)^j*self.moment(j)*X^j*Y^(k-j)
		return symk(k,P)	

	def act_right(self,gam):
		"""return self|gam"""
		a=gam[0,0]
		b=gam[0,1]
		c=gam[1,0]	
		d=gam[1,1]
		new_moments = (self.char())(a) * Matrix(self.moments)*form_acting_matrix_on_dist(self.p,self.num_moments(),self.weight,a,b,c,d)
			
		return dist(self.p,self.weight,new_moments[0],self.char())

	def lift_to_dist_fam(self,deg,r,w):
		v=[self.moment(j)+0*w for j in range(0,self.num_moments())]
		return dist_fam(self.p,deg,r,vector(v),self.char())
		
	def solve_diff_eqn(self):
		#assert self.moments[0][0]==0, "not total measure zero"
		#print "result accurate modulo p^",self.moment(0).valuation(self.p)
		mu=self.zero()
		for m in range(1,self.num_moments()):
			mu=mu+eta(m-1,self.p,self.weight,self.num_moments(),self.char()).scale(self.moment(m)/m)
		return mu

	def lift(self):
		"""increases the number of moments by 1"""
		w=self.moments.list()+[0]
		return dist(self.p,self.weight,vector(w),self.char())
		


@cached_function				
def form_acting_matrix_on_dist(p,M,k,a,b,c,d):
	"""forms a large M x M matrix say G such that if v is the vector of moments of a distribution mu, then v*G is the vector of moments of mu|[a,b;c,d]"""
	assert (a%p != 0) and (c%p == 0), "acting by bad matrix"

	R=PowerSeriesRing(QQ,'y',default_prec=M)
	y=R.gen()

	scale=(b+d*y)/(a+c*y)
	t=((a+c*y)^k).truncate(M)

	A = []
	for i in range(0,M):
		temp1=t.list();
		d=len(temp1)
		for j in range(d,M):
			temp1 = temp1 + [0]	
		#while len(temp1)>M:
		#	temp1.pop()
		A = A + [temp1]	
		t=(t*scale).truncate(M)
	q=p^M
	B=Matrix(QQ,A).transpose()
	for r in range(0,M):
		for c in range(0,M):
			#B[r,c]=B[r,c]%(p^(M-c))
			B[r,c]=B[r,c]%(q)

	return B

def eta(i,p,k,M,char=None):
	"""helper function in solving the difference equation -- see Lemma 4.4 of [PS]"""
	v=[0 for j in range(0,i)]+[binomial(j,i)*bernoulli(j-i) for j in range(i,M)]
	return dist(p,k,vector(v),char)

def random_dist(p,k,M,char=None):	
	"""Returns a random distribution with prime p, weight k, and M moments"""
	moments=vector([ZZ(floor(random()*p^M)) for i in [1..M]])
	
	return dist(p,k,moments,char)


