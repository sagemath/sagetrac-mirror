from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0

@cached_function
def invert(a,b,c,d):
	return Matrix(2,2,[d,-b,-c,a])

@cached_function
def unimod_matrices(r,s):
	"""connects the rational number r/s to 1/0"""
	if s<>0:
		v=[]
		list=convergents(r/s)
		for j in range(0,len(list)-1):
			a=list[j].numerator()
			c=list[j].denominator()
			b=list[j+1].numerator()
			d=list[j+1].denominator()
			v=v+[Matrix(ZZ,[[(-1)**(j+1)*a,b],[(-1)**(j+1)*c,d]])]
		return [Matrix(ZZ,[[1,list[0].numerator()],[0,list[0].denominator()]])]+v
	else:
		return []
		
def flip(A):
	return Matrix(2,2,[-A[0,1],A[0,0],-A[1,1],A[1,0]])

@cached_function
def unimod_matrices2(r,s):
	"""connects the rational number 1/0 to r/s"""
	if s<>0:
		v=[]
		list=convergents(r/s)
		for j in range(0,len(list)-1):
			a=list[j].numerator()
			c=list[j].denominator()
			b=list[j+1].numerator()
			d=list[j+1].denominator()
			v=v+[flip(Matrix(ZZ,[[(-1)**(j+1)*a,b],[(-1)**(j+1)*c,d]]))]
		return [flip(Matrix(ZZ,[[1,list[0].numerator()],[0,list[0].denominator()]]))]+v
	else:
		return []

def basic_hecke_matrix(a,ell):
	if a<=ell:
		return Matrix(2,2,[1,a,0,ell])
	else:
		return Matrix(2,2,[ell,0,0,1])

@cached_function
def prep_hecke_individual(ell,N,M,m):
	ans=[[] for a in range(len(M.mats))]
	for a in range(ell):
		gama=basic_hecke_matrix(a,ell)
		t=gama*M.mats[M.gens[m]]
		v=unimod_matrices2(t[0,0],t[1,0])+unimod_matrices(t[0,1],t[1,1])
		for b in range(len(v)):
			A=v[b]
			i=M.P.index(A[1,0],A[1,1])
			j=M.P1_to_mats[i]
			B=M.mats[j]
			C=invert(A[0,0],A[0,1],A[1,0],A[1,1])
			gaminv=B*C
			ans[j]=ans[j]+[gaminv*gama]
	if N%ell<>0:
		gama=basic_hecke_matrix(ell+1,ell)
		t=gama*M.mats[M.gens[m]]
		v=unimod_matrices2(t[0,0],t[1,0])+unimod_matrices(t[0,1],t[1,1])
		for b in range(len(v)):
			A=v[b]
			i=M.P.index(A[1,0],A[1,1])
			j=M.P1_to_mats[i]
			B=M.mats[j]
			C=invert(A[0,0],A[0,1],A[1,0],A[1,1])
			gaminv=B*C
			ans[j]=ans[j]+[gaminv*gama]

	return ans

@cached_function
def prep_hecke(ell,N,M):
	ans=[]
	for m in range(len(M.gens)):
		ans=ans+[prep_hecke_individual(ell,N,M,m)]
	return ans
			

class modsym(SageObject):
	def __init__(self,level,data,manin,full_data=None):
		self.level=level
		self.data=data
		self.manin=manin
		if full_data<>None:
			self.full_data=full_data
		else:
			self.full_data=0
			
	def __repr__(self):
		return repr(self.data)

	def ngens(self):
		return len(self.manin.gens)

	def __add__(self,right):
		assert self.level==right.level, "the levels are different"
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j]+right.data[j]]
		if self.full_data<>0 and right.full_data<>0:
			w=[]
			for j in range(0,len(self.full_data)):
				w=w+[self.full_data[j]+right.full_data[j]]
		else:
			w=0
			
		C=type(self)
		return C(self.level,v,self.manin,w).normalize()

	def normalize(self):
		for j in range(len(self.data)):
			self.data[j]=self.data[j].normalize()
		return self

	def scale(self,left):
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].scale(left)]
		if self.full_data<>0:
			w=[]
			for j in range(0,len(self.full_data)):
				w=w+[self.full_data[j].scale(left)]
		else:
			w=0

		C=type(self)
		return C(self.level,v,self.manin,w)
		
	def __sub__(self,right):
		return self+right.scale(-1)

	def __cmp__(self,right):
		return cmp((self.level,self.data),(right.level,right.data))

	def zero_elt(self):
		return self.data[0].zero()

	def zero(self):
		v=[self.zero_elt() for i in range(0,len(self.data))]
		C=type(self)
		return C(self.level,v,self.manin)

	def compute_full_data_from_gen_data(self):
		ans=[]
		for m in range(len(self.manin.mats)):
			v=self.manin.rels[m]
			t=self.data[0].zero()
			for k in range(len(v)):
				j=v[k][2]
				r=self.manin.gens.index(j)
				t=t+self.data[r].act_right(v[k][1]).scale(v[k][0])
			ans=ans+[t]
		self.full_data=ans
	
	def eval_sl2(self,A):
		i=self.manin.P.index(A[1,0],A[1,1])
		j=self.manin.P1_to_mats[i]
		B=self.manin.mats[j]
		C=invert(A[0,0],A[0,1],A[1,0],A[1,1])
		gaminv=B*C
		if self.full_data<>0:
			return self.full_data[j].act_right(gaminv)
		else:
			v=self.manin.rels[j]
			t=self.data[0].zero()
			for k in range(len(v)):
				m=v[k][2]
				r=self.manin.gens.index(m)
				t=t+self.data[r].act_right(v[k][1]*gaminv).scale(v[k][0])
			return t

	def eval(self,A):
		a=A[0,0]
		b=A[0,1]
		c=A[1,0]
		d=A[1,1]
		v1=unimod_matrices(b,d)
		v2=unimod_matrices(a,c)
		ans=self.zero_elt()
		for j in range(0,len(v1)):
			ans=ans+self.eval_sl2(v1[j])
		for j in range(0,len(v2)):
			ans=ans-self.eval_sl2(v2[j])
		return ans
				
	def act_right(self,gamma):
		v=[]
		for j in range(0,len(self.data)):
			rj=self.manin.gens[j]
			v=v+[self.eval(gamma*self.manin.mats[rj]).act_right(gamma)]

		C=type(self)		
		return C(self.level,v,self.manin).normalize()
	
	def plus_part(self):
		return self.act_right(Matrix(2,2,[1,0,0,-1]))+self

	def minus_part(self):
		return self.act_right(Matrix(2,2,[1,0,0,-1]))-self

	def normalize_full_data(self):
		if (self.full_data != 0):
			for j in range(len(self.full_data)):
				self.full_data[j]=self.full_data[j].normalize()

	def hecke2(self,ell):
		if self.full_data==0:
			self.compute_full_data_from_gen_data()
		psi=self.zero()
		for a in range(0,ell): 
			psi=psi+self.act_right(Matrix(ZZ,[[1,a],[0,ell]]))
		if self.level%ell<>0:
			psi=psi+self.act_right(Matrix(ZZ,[[ell,0],[0,1]]))
		return psi.normalize()

	def hecke(self,ell):
		if self.full_data==0:
			self.compute_full_data_from_gen_data()
			self.normalize_full_data()
		psi=self.zero()
		v=prep_hecke(ell,self.level,self.manin)
		for m in range(len(self.manin.gens)):
			for j in range(len(self.manin.mats)):
				for r in range(len(v[m][j])):
					psi.data[m]=psi.data[m]+self.full_data[j].act_right(v[m][j][r])
		return psi.normalize()

	def grab_relations(self):
		v=[]
		for r in range(len(self.manin.gens)):
			for j in range(len(self.manin.rels)):
				R=self.manin.rels[j]
				if (len(R)==1) and (R[0][2]==self.manin.gens[r]):
					if R[0][0]<>-1 or R[0][1]<>Id:
						v=v+[R]
		return v

	def check_loop(self):
		list=self.grab_relations()
		t=self.zero_elt()
		for j in range(2,len(list)):
			R=list[j]
			index=R[0][2]
			rj=self.manin.gens.index(index)
			t=t+self.data[rj].act_right(R[0][1]).scale(R[0][0])
		return self.data[0]-self.data[0].act_right(Matrix(2,2,[1,1,0,1]))+t

def zero_to_ratl(r,p):
	"""returns a Gamma_0(p) matrix which takes 0 to r"""
	c=numerator(r)
	d=denominator(r)
	assert d%p<>0, "not Gamma_0(p)-equivalent to 0"
	g,x,y=xgcd(d,-p*c)
	return Matrix(2,2,[x,c,p*y,d])
	
def eisen_gamma0p(p,M):
	assert M.act_right(Matrix(2,2,[1,0,p,1]))==M, "not a good element to use"
	manin=manin_relations(p)
	v=[]
	for j in range(0,len(manin.gens)):
		rj=manin.gens[j]
		A=manin.mats[rj]
		a=A[0,0]
		b=A[0,1]
		c=A[1,0]
		d=A[1,1]
		t=M.zero()
		if d%p<>0:
			gam1=zero_to_ratl(b/d,p)
			t=M.act_right(gam1**(-1))
		if c%p<>0:
			gam2=zero_to_ratl(a/c,p)
			t=t-M.act_right(gam2**(-1))
		v=v+[t]
	return modsym(p,v,manin)	
		

