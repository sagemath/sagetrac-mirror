from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.polynomial.all import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.arith import binomial, bernoulli
from sage.modules.free_module_element import vector, zero_vector
from sage.matrix.all import Matrix
from sage.misc.prandom import random
from sage.functions.other import floor


class dist(SageObject):
	def __init__(self,p,k,moments):
		"""A distribution is stored as a vector whose j-th entry is the j-th moment of the distribution.  The j-th entry is stored modulo p^(N-j) where N is the total number of moments.  (This is the accuracy that is maintained after acting by Gamma_0(p).)

Inputs:
	p -- prime
	k -- weight (used to determine the action of matices)
	moments	-- the list of moments given as a vector"""
		self.p=p
		self.weight=k
		self.moments=moments

	def moment(self,n):
		"""returns the n-th moment"""
		return self.moments[n]

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
		return dist(self.p,self.weight,self.moments+right.moments)

	#should try to figure out __lmul__
	def _lmul_(self,other):
		"""Scales self by the constant left"""
		return dist(self.p,self.weight,other*self.moments)

	def scale(self,left):
		"""Scales self by the constant left"""
		return dist(self.p,self.weight,left*self.moments)

	def __sub__(self,right):
		"""Returns self-right"""
		return self+right.scale(-1)

	def __cmp__(self,right):
		"""Returns true is both the moments and the weight agree"""
		return cmp((self.weight,self.moments),(right.weight,right.moments))

	def zero(self):
		"""returns a distribution with all zero moments (and the same prime and weight)"""
		return dist(self.p,self.weight,zero_vector(QQ,self.num_moments()))

	def valuation(self):
		"""returns the highest power of p which divides all moments of the distribution"""
		p=self.p
		return min([self.moment(a).valuation(p) for a in range(self.num_moments())])

	def normalize(self):
		"""reduced modulo Fil^N -- that is the i-th moments is reduced modulo p^(N-i)"""
		p=self.p
		assert self.valuation() >= 0, "moments not integral in normalization"

		v=vector([self.moment(i)%(p**(self.num_moments()-i)) for i in range(0,self.num_moments())])
		return dist(self.p,self.weight,v)

	def normalize_aws(self):
		"""Arizona Winter School filtration (for families) -- reduces the j-th moment modulo p^(floor((N+1-j)*(p-2)/(p-1))"""
		p=self.p
		N=self.num_moments()
		assert self.valuation() >= 0, "moments not integral in normalization"

		v=vector([self.moment(j)%(p**(floor((N+1-j)*(p-2)/(p-1)))) for j in range(0,self.num_moments())])
		return dist(self.p,self.weight,v)

	def change_precision(self,M):
		"""only holds on to M moments"""
		assert M<=self.num_moments(),"not enough moments"

		v=[self.moment(i) for i in range(M)]
		mu=dist(self.p,self.weight,vector(v))
		return mu.normalize()

	def specialize(self):
		"""specializes to weight k -- i.e. projects to Sym^k"""
		k=self.weight
		if k==0:
			# R.<X,Y>=PolynomialRing(QQ,2)
			R = PolynomialRing(QQ,('X','Y'))
			X,Y = R.gens()
			P=0
			for j in range(0,k+1):
				P=P+binomial(k,j)*((-1)**j)*self.moment(j)*(X**j)*(Y**(k-j))
			return symk(k,P)

	def act_right(self,gam):
		"""return self|gam"""
#		print(gam)
		a=gam[0,0]
		b=gam[0,1]
		c=gam[1,0]
		d=gam[1,1]
		new_moments=Matrix(self.moments)*form_acting_matrix_on_dist(self.p,self.num_moments(),self.weight,a,b,c,d)
		return dist(self.p,self.weight,new_moments[0])

	def lift_to_dist_fam(self,deg,w):
	        v=[self.moment(j)+0*w for j in range(0,self.num_moments())]
		return dist_fam(self.p,deg,vector(v))

	def solve_diff_eqn(self):
		# assert self.moments[0][0]==0, "not total measure zero"
		# print "result accurate modulo p^",self.moment(0).valuation(self.p)
		mu=self.zero()
		for m in range(1,self.num_moments()):
			mu=mu+eta(m-1,self.p,self.weight,self.num_moments()).scale(self.moment(m)/m)
		return mu

	def lift(self):
		"""increases the number of moments by 1"""
		w=self.moments.list()+[0]
		return dist(self.p,self.weight,vector(w))



# @cached_function
def form_acting_matrix_on_dist(p,M,k,a,b,c,d):
	"""forms a large M x M matrix say G such that if v is the vector of moments of a distribution mu, then v*G is the vector of moments of mu|[a,b;c,d]"""

#	print("Checking...")
#	print(a,b,c,d)
#	print(p)

	assert (a%p != 0) and (c%p == 0), "acting by bad matrix"

	R=PowerSeriesRing(QQ,'y',default_prec=M)
	y=R.gen()

	scale=(b+d*y)/(a+c*y)
	t=((a+c*y)**k).truncate(M)

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
	q=p**M
	B=Matrix(QQ,A).transpose()
	for r in range(0,M):
		for c in range(0,M):
			#B[r,c]=B[r,c]%(p**(M-c))
			B[r,c]=B[r,c]%(q)

	return B

def eta(i,p,k,M):
	"""helper function in solving the difference equation -- see Lemma 4.4 of [PS]"""
	v=[0 for j in range(0,i)]+[binomial(j,i)*bernoulli(j-i) for j in range(i,M)]
	return dist(p,k,vector(v))

def random_dist(p,k,M):
	"""Returns a random distribution with prime p, weight k, and M moments"""
	moments=vector([ZZ(floor(random()*(p**M))) for i in range(1,M+1)])

	return dist(p,k,moments)
