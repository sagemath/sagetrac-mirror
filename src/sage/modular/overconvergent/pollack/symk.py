
import sage.rings.polynomial.multi_polynomial_element as MPol
from sage.structure.sage_object import SageObject
from sage.rings.rational_field import QQ
from sage.rings.polynomial.all import PolynomialRing
from sage.rings.arith import valuation, binomial
from sage.modules.free_module_element import vector

class symk(SageObject):
	def __init__(self,k,poly=None,base_ring=QQ):
		"""A symk object is stored as a homogeneous polynomial of degree k

Inputs:
	k -- weight
	poly -- homogeneous of weight k in X and Y
	base_ring -- ring containing the coefficients of poly"""

		#assert (poly == None) or (poly == 0) or (k==0) or (poly.is_homogeneous()), "must enter a homogeneous polynomial"
		#if (poly != None) and (poly !=0) and (k<>0):
		#	assert poly.total_degree() == k, "the weight is incorrect"
		if poly != None:
			self.poly=poly
		else:
			# R.<X,Y>=PolynomialRing(base_ring,2)
			# self.poly=0*X
			R = PolynomialRing(base_ring,('X','Y'))
			self.poly=0*R.gens()[0]
		self.weight=k
		self.base_ring=base_ring

	def __repr__(self):
		return repr(self.poly)

	def coef(self,j):
		"""returns coefficient of X^j*Y^(k-j)"""
		v=self.vars()
		X=v[0]
		Y=v[1]
		k=self.weight

		return self.poly[j,k-j]

	def __add__(self,right):
		"""returns self+right"""
		#assert self.weight == right.weight, "weights are unequal"
		return symk(self.weight,self.poly+right.poly)

	def scale(self,left):
		"""returns left*self"""
		return symk(self.weight,left*self.poly)

	def __sub__(self,right):
		"""returns self-right"""
		#assert self.weight == right.weight, "weights are unequal"
		return self+right.scale(-1)

	def __cmp__(self,right):
		return cmp((self.weight,self.poly),(right.weight,right.poly))

	def zero(self):
		return symk(self.weight)

	def vars(self):
		"""Returns the variables defining Sym^k"""
		return self.poly.parent().gens()

	def valuation(self,p):
		"""returns the exponent of the highest power of p which divides all coefficients of self"""
		#assert self.base_ring==QQ, "need to be working over Q in valuation"
		k=self.weight
		v=self.vars()
		X=v[0]
		Y=v[1]
		v=[]
		for j in range(k+1):
			v=v+[valuation(QQ(self.poly.coefficient((X**j)*(Y**(k-j)))),p)]
		return min(v)

	def normalize(self):
		return self

	def map(self,psi):
		"""psi is a map from the base_ring to Qp and this function applies psi to all polynomial coefficients and then lifts them to QQ"""
		#assert psi.domain()==self.base_ring
		k=self.weight
		# S.<X,Y>=PolynomialRing(QQ)
		S = PolynomialRing(QQ,('X','Y'))
		X,Y = S.gens()

		ans=0*X
		for j in range(k+1):
			ans=ans+psi(self.coef(j)).lift()*(X**j)*(Y**(k-j))
		temp=copy(self)
		temp.poly=ans
		temp.base_ring=psi.codomain()

		return temp


	def act_right(self,gam):
		"""returns self|gam where (P|[a,b;c,d])(X,Y) = P(dX-cY,-bX+aY)"""
		if self.weight==0:
			return self
		else:
			a=gam[0,0]
			b=gam[0,1]
			c=gam[1,0]
			d=gam[1,1]
			v=self.vars()
			X=v[0]
			Y=v[1]
			return symk(self.weight,self.poly(d*X-c*Y,-b*X+a*Y))

	def lift_to_dist(self,p,N):
		from sage.modular.overconvergent.pollack.dist import dist

		"""returns some (p-adic) distribution with N moments which specializes to self"""
		k=self.weight
		#assert N >= k+1, "need more moments"
		if k==0:
			moments=[QQ(self.poly)]
		else:
			v=self.vars()
			X=v[0]
			Y=v[1]
			moments=[]
			for j in range(k+1):
				moments=moments+[QQ(self.poly.coefficient((X**j)*(Y**(k-j)))/(binomial(k,j)*(-1)^j))]
		while len(moments)<N:
			moments=moments+[0]
		mu=dist(p,k,vector(moments))
		if mu.valuation()<0:
			print "scaling by ",p,"^",-mu.valuation()," to keep things integral"
			mu=mu.scale(p**(-mu.valuation()))
		return mu
