from sage.structure.sage_object import SageObject

class ps_fam(SageObject):
	def __init__(self,p,series):
		"""p denotes the prime, self denotes the power series in z stored as a vector, and deg denotes that we are working modulo z^deg"""
		self.p=p
		self.series=series
		self.deg=len(series)
			
	def __repr__(self):
		return repr(self.series)
	
	def gen(self):
		return self.series[0].parent().gen()
		
	def specialize(self,k,r):
		"""evaluates at ((1+p)^k-1)/p^r"""
		v=[]
		w=self.gen
		for j in range(0,len(self.series)):
			v=v+[self.series[j].substitute(w=((1+self.p)**k-1))/(self.p)^r]
		return ps_fam(self.p,v)

	def normalize(self):
		"""reduces all coefficients modulo p^[N*(p-1)/p] where N=self.deg"""
		N=self.deg
		p=self.p
		r=floor((N+1)*(p-1)/p)
		v=[self.series[a]%(p^r) for a in range(self.deg)]
		return ps_fam(p,v)

def ps_normalize(f,p,N):
	"""reduces all of the coefficients of the power series modulo p^[N*(p-1)/p]"""
	v=Sequence(f)
	r=floor((N+1)*(p-1)/(p))
	v=[v[a]%(p^r) for a in range(len(v))]
	S=f.parent()
	f=S(v)

	return f
		
def logp_fcn(p,N,z):
	"""this is the *function* on Z_p^* which sends z to log_p(z) using a power series truncated at N terms"""
	R=pAdicField(p)
	z=z/R.teichmuller(z)
	ans=0
	for m in range(1,N-1):
		ans=ans+(-1)**(m-1)*(z-1)**m/m

	return ans

def logpp(p,N,z):
	"""returns the (integral) power series for log_p(1+p*z) -- extra p here!"""
	ans=0
	for m in range(1,N-1):
		ans=ans+(-1)^(m-1)*(p*z)^m/m

	return ans

def logpp_gam(p,N,z):
	"""returns the (integral) power series log_p(1+p*z)*(1/log_p(1+p)) where the denominator is computed with some accuracy"""
	L=logpp(p,N,z)
	loggam=ZZ(logp_fcn(p,N*p^2,1+p))
	return ps_normalize(L/loggam,p,N)

@cached_function
def logpp_binom(n,p,N,z):
	"""returns the (integral) power series p^n*(log_p(1+p*z)/log_p(1+p) choose n)"""
	prod=1
	L=logpp_gam(p,N,z)
	for j in range(0,n):
		prod=prod*(L-j)
	prod=p^n*prod/factorial(n)
	
	return ps_normalize(prod.truncate(N),p,N)
	
@cached_function
def aut(p,Mw,Mm,a,c,w):
	#extra p built in here!!
	R=w.parent()
	S=PolynomialRing(R,'zz')
	SS=PolynomialRing(QQ,'yy')
	yy=SS.gen()
	zz=S.gens()[0]
	
	ans=1+0*zz
	for n in range(1,Mw): 
		LB=logpp_binom(n,p,Mm,yy)
		ta=teich(a,p,2*max(Mw,Mm))
		v=(a/ta-1)/p+c/(p*ta)*zz
		ans=ans+w^n*LB(v)

	v=Sequence(ans)
	while len(v)<Mm:
		v=v+[0*w]
	return ps_fam(p,v)


	
	