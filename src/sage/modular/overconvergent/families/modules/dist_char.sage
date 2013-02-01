from sage.structure.sage_object import SageObject

class dist_char(SageObject):
	def __init__(self,p,k,chi,moments):
		"""A distribution with character is stored as a vector whose j-th entry is the j-th moment of the distribution.  The j-th entry is stored modulo p^(N-j) where N is the total number of moments.  (This is the accuracy that is maintained after acting by Gamma_0(p).)  

Inputs: 
	p -- prime 
	k -- weight (used to determine the action of matices)
	j -- power of the mod p cyclotomic character
	moments	-- the list of moments"""
		self.p=p
		self.weight=k
		self.char=chi
		self.moments=vector(Sequence(moments))
			
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
		return dist_char(self.p,self.weight,self.char,self.moments+right.moments)

	#should try to figure out __lmul__
	def _lmul_(self,other):
		"""Scales self by the constant left"""
		return dist_char(self.p,self.weight,self.char,other*self.moments)

	def scale(self,left):
		"""Scales self by the constant left"""
		return dist_char(self.p,self.weight,self.char,left*self.moments)

	def __sub__(self,right):
		return self+right.scale(-1)

	def __cmp__(self,right):
		"""Returns true is both the moments and the weight agree"""
		return cmp((self.weight,self.moments,self.char),(right.weight,right.moments,right.char))

	def zero(self):
		"""return a distribution with all zero moments (and the same prime and weight"""
		return dist_char(self.p,self.weight,self.char,vector([0 for i in range(0,len(self.moments))]))

	def series(self,m):
		"""returns a polynomial in the variable m whose j-th coefficient is the j-th moment of the distribution"""
		ans=0
		for j in range(0,self.num_moments()):
			ans=ans+m**j*self.moments[j]
		return ans
		
	def valuation(self):
		"""returns the highest power of p which divides all moments of the distribution"""
		p=self.p
		coefs=Sequence(self.moments)
		return min([coefs[i].valuation(p) for i in range(0,self.num_moments())])

	def normalize(self):
		"""reduced modulo Fil^N"""
		p=self.p
		assert self.valuation() >= 0, "moments not integral in normalization"
		
		v=vector([self.moments[i]%(p**(self.num_moments()-i)) for i in range(0,self.num_moments())])	
		return dist_char(self.p,self.weight,self.char,v)

	def change_precision(self,M):
		"""reduces modulo Fil^M"""
		assert M<=self.num_moments(),"not enough moments"

		v=[self.moments[i] for i in range(M)]
		mu=dist_char(self.p,self.weight,self.char,v)
		return mu.normalize()

	def specialize(self):
		"""specializes to weight k -- i.e. projects to Sym^k  -- NO CHARACTER HERE!!"""
		assert 0==1, "didn't program character here yet"
		k=self.weight
		if k==0:
			return symk(0,self.moments[0])
		else:
			R.<X,Y>=PolynomialRing(QQ,2)
			P=0
			for j in range(0,k+1):
				P=P+binomial(k,j)*(-1)^j*self.moments[j]*X^j*Y^(k-j)
			return symk(k,P)	

	def act_right(self,gam):
		a=gam[0,0]
		b=gam[0,1]
		c=gam[1,0]	
		d=gam[1,1]
		new_moments=(self.moments*form_acting_matrix_on_dist(self.p,self.num_moments(),self.weight,a,b,c,d))*ZZ(self.char(a))
		return dist_char(self.p,self.weight,self.char,new_moments)

	def solve_diff_eqn(self):
		assert self.moments[0]==0, "not total measure zero"
		mu=self.zero()
		for m in range(1,self.num_moments()):
			mu=mu+eta_char(m-1,self.p,self.weight,self.char,self.num_moments()).scale(self.moments[m]/m)
		return mu	

def eta_char(i,p,k,j,M):
	v=[0 for j in range(0,i)]+[binomial(j,i)*bernoulli(j-i) for j in range(i,M)]
	return dist_char(p,k,j,v)

def random_dist_char(p,k,chi,M):	
	"""Returns a random distribution with prime p, weight k, character chi, and M moments"""
	moments=[ZZ(floor(random()*p^M)) for i in [1..M]]
	
	return dist_char(p,k,chi,moments)

