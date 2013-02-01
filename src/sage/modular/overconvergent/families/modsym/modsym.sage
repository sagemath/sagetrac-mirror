from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
from sage.matrix.matrix_integer_2x2 import Matrix_integer_2x2 as mi2x2
import sage.matrix.all as matrix

M2Z = matrix.MatrixSpace(ZZ,2)

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
		return self.act_right(mi2x2(M2Z, [1, 0, 0, -1], copy = False, coerce = False)) + self

	def minus_part(self):
		return self.act_right(mi2x2(M2Z, [1, 0, 0,- 1], copy = False, coerce = False)) - self

	def normalize_full_data(self):
		if (self.full_data != 0):
			for j in range(len(self.full_data)):
				self.full_data[j]=self.full_data[j].normalize()

	def hecke2(self,ell):
		if self.full_data==0:
			self.compute_full_data_from_gen_data()
		psi=self.zero()
		for a in range(0,ell): 
			psi=psi+self.act_right(mi2x2(M2Z, [1, a, 0, ell], copy = False, coerce = False))
		if self.level%ell<>0:
			psi=psi+self.act_right(mi2x2(M2Z, [ell, 0, 0, 1], copy = False, coerce = False))
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

	#for testing speed
#	def hecke(self,ell):
#		print "In hecke({0})".format(ell)
#		cur_time = walltime()
#		if self.full_data==0:
#			print "Computing full data"
#			self.compute_full_data_from_gen_data()
#			self.normalize_full_data()
#			print "Time full data:", walltime() - cur_time
#			cur_time = walltime()
#		psi=self.zero()
#		v=prep_hecke(ell,self.level,self.manin)
#		print "Time prep_hecke:", walltime() - cur_time
#		cur_time = walltime() 
#		for m in range(len(self.manin.gens)):
#			for j in range(len(self.manin.mats)):
#				for r in range(len(v[m][j])):
#					psi.data[m]=psi.data[m]+self.full_data[j].act_right(v[m][j][r])
##			psi.data[m]=psi.data[m]+sum([self.full_data[j].act_right(v[m][j][r]) for j in range(len(self.manin.mats)) for r in range(len(v[m][j]))])
#		print "Time loop:", walltime() - cur_time
#		cur_time = walltime() 
#		return psi.normalize()

	def grab_relations(self):
		v=[]
		for r in range(len(self.manin.gens)):
			for j in range(len(self.manin.rels)):
				R=self.manin.rels[j]
				if (len(R)==1) and (R[0][2]==self.manin.gens[r]):
					if R[0][0]<>-1 or R[0][1]<>Id:
						v=v+[R]
		return v
	
	def aplist(self, max_ell):
		ells = prime_range(max_ell + 1)
		a_ells = []
		for ell in ells:
			b, Lambda, m = self.is_Tq_eigen(ell)
			if not b:
				raise NotImplementedError("q-expansion only implemented for eigensymbols.")
			a_ells.append(Lambda)
		return a_ells

	def check_loop(self):
		t = self.zero_elt()
		for j in range(self.ngens()):
			rj = self.manin.gens[j]
			if (self.manin.twotor.count(rj) == 0) and (self.manin.threetor.count(rj) == 0):
				t = t + self.data[j].act_right(self.gen_rel_mat(j)) - self.data[j]
			else:
				t = t - self.data[j]
		return t

	############################################################################################################
	## Takes the i-th generator defining self and returns the associated matrix gamma_i -- i.e. if we are not in 	## two or three torsion case, we return the matrix for which 
	##
	##  phi(D_infty) | Delta = \sum phi(D_i) | (gamma_i - 1)  + two and three torsion stuff
	##
	## In the two torsion case, we return the matrix such that phi(D'_i)|(gamma_i+1) = 0 
	## In the three torsion case, we return the matrix such that phi(D''_i)|(gamma_i^2+gamma_i+1) = 0 
	############################################################################################################
	def gen_rel_mat(self,i):
		return self.manin.gen_rel_mat(i)
