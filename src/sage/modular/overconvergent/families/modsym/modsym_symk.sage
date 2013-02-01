class modsym_symk(modsym):
	def ms(self):
		"""demotes to a regular modular symbol"""
		return modsym(self.level,self.data,self.manin)

	def weight(self):
		"""returns the weight of any value of self"""
		return self.data[0].weight

	def base_ring(self):
		return self.data[0].base_ring

	def valuation(self,p):
		"""returns the valuation of self at p -- i.e. the exponent of the largest power of p which divides all of the coefficients of all values of self"""
		v=self.data
		return min([v[j].valuation(p) for j in range(0,len(v))])

	def p_stabilize(self,p,alpha):
		"""p-stablizes self to form an eigensymbol for U_p with eigenvalue alpha"""
		N=self.level
		assert N%p<>0, "The level isn't prime to p"

		pp=Matrix(ZZ,[[p,0],[0,1]])
		manin=manin_relations(N*p)
		v=[]
		for j in range(0,len(manin.gens)):
			rj=manin.gens[j]
			v=v+[self.eval(manin.mats[rj])-self.eval(pp*manin.mats[rj]).act_right(pp).scale(1/alpha)]
		return modsym_symk(N*p,v,manin)


	def p_stabilize_ordinary(self,p,ap,M):
		"""returns the unique p-ordinary p-stabilization of self"""
		N=self.level
		k=self.data[0].weight
		assert N%p<>0, "The level isn't prime to p"
		assert ap%p<>0, "Not ordinary!"

		"""makes alpha the unit root of Hecke poly"""
		R=PolynomialRing(pAdicField(p,M),'y') 
		y=R.gen()
		f=y^2-ap*y+p^(k+1)
		v=f.roots()
		if Integer(v[0][0])%p<>0:
			alpha=Integer(v[0][0])
		else:
			alpha=Integer(v[1][0])

		if self.full_data == 0:
			self.compute_full_data_from_gen_data()

		return self.p_stabilize(p,alpha)

	def p_stabilize_critical(self,p,ap,M):
		N=self.level
		assert N%p<>0, "The level isn't prime to p"
		assert ap%p<>0, "Not ordinary!"

		"""makes alpha the non-unit root of Hecke poly"""
		R=PolynomialRing(Qp(p,M),'y') 
		y=R.gen()
		f=y**2-ap*y+p
		v=f.roots()
		if Integer(v[0][0])%3==0:
			alpha=Integer(v[0][0])
		else:
			alpha=Integer(v[1][0])

		if self.full_data == 0:
			self.compute_full_data_from_gen_data()

		return self.p_stabilize(p,alpha)

	def is_Tq_eigen_old(self,q,p,M):
		"""determines of self is an eigenvector for T_q modulo p^M"""
		selfq=self.hecke(q)
		r=0
		while ((self.data[r].coef(0))==0) and (r<=self.ngens()):
			r=r+1
		if r>self.ngens():
			#all coefficients of Y^k are zero which I think forces it not to be eigen
			return False
		else:	
			c=selfq.data[r].coef(0)/self.data[r].coef(0)
			val=(selfq-self.scale(c)).valuation(p)
			if val>=M:
				return True,c
			else:
				return False

	def is_Tq_eigen(self,q,p,M):
		"""determines if self is an eigenvector for T_q modulo p^M"""
		assert self.valuation(p) >= 0, "is_Tq_eigen: Not integral at p"

		selfq=self.hecke(q)
		r=0
		while ((self.data[r].coef(0))%(p)==0) and (r < self.ngens()):
			r=r+1
		if r>=self.ngens():
			#all coefficients of Y^k are zero which I think forces it not to be eigen
			return False
		else:	
			c=selfq.data[r].coef(0)/self.data[r].coef(0)
			val=(selfq-self.scale(c)).valuation(p)
			if val>=M:
				return [True, c, M]
			else:
				return [False, None, M]

	def aplist(self, max_ell, p, M):
		ells = prime_range(max_ell + 1)
		a_ells = []
		for ell in ells:
			b, Lambda, m = self.is_Tq_eigen(ell, p, M)
			if not b:
				raise NotImplementedError("q-expansion only implemented for eigensymbols.")
			a_ells.append(Lambda)
		return a_ells

	def lift_to_OMS(self,p,M):
		"""returns a (p-adic) overconvergent modular symbols with M moments which lifts self up to an Eisenstein error"""
		v=[]
		## this loop runs through each generator and lifts the value of self on that generator to D
		for j in range(1,len(self.manin.gens)):
			rj = self.manin.gens[j]
			if (self.manin.twotor.count(rj) == 0) and (self.manin.threetor.count(rj) == 0):
				v = v + [self.data[j].lift_to_dist(p,M)]
			elif (self.manin.twotor.count(rj) != 0):
				## Case of two torsion (See [PS] section 4.1)
				gam = self.manin.gen_rel_mat(j)
				mu = self.data[j].lift_to_dist(p,M)
				v = v + [(mu - mu.act_right(gam)).scale(1/2)]
			else:
				## Case of three torsion (See [PS] section 4.1)	
				gam = self.manin.gen_rel_mat(j)
				mu = self.data[j].lift_to_dist(p,M)
				v = v + [(mu.scale(2) - mu.act_right(gam) - mu.act_right(gam^2)).scale(1/3)]

		t = v[0].zero()
		## This loops adds up around the boundary of fundamental domain except the two verticle lines
		for j in range(1,len(self.manin.gens)):
			rj = self.manin.gens[j]
			if (self.manin.twotor.count(rj) == 0) and (self.manin.threetor.count(rj) == 0):
				t = t + v[j-1].act_right(self.manin.gen_rel_mat(j)) - v[j-1]
			else:
				t = t - v[j-1]


		## t now should be sum Phi(D_i) | (gamma_i - 1) - sum Phi(D'_i) - sum Phi(D''_i)
		## (Here I'm using the opposite sign convention of [PS1] regarding D'_i and D''_i)

		assert t.normalize().moment(0) == 0, "Not total measure 0 in lifting OMS" 

		mu = t.solve_diff_eqn()

		v = [mu.scale(-1)] + v
	 
		return modsym_dist(self.level,v,self.manin)	

	def lift_to_OMS_eigen(self,p,M,verbose=True):
		"""returns Hecke-eigensymbol OMS lifting phi -- phi must be a p-ordinary eigensymbol"""
		v=self.is_Tq_eigen(p,p,M)
		assert v[0], "not eigen at p!"
		ap=v[1]
		k=self.weight()
		Phi=self.lift_to_OMS(p,M)
		s=-Phi.valuation()
		if s>0:
			if verbose:
				print "Scaling by ",p,"^",s
			Phi=Phi.scale(p^(-Phi.valuation()))
		Phi=Phi.normalize()
		if verbose:
			print "Applying Hecke"
		Phi=Phi.hecke(p).scale(1/ap)
		if verbose:
			print "Killing eisenstein part"
		if (ap%(p^M))<>1:
			Phi=(Phi-Phi.hecke(p)).scale(1/(1-ap))
			e=(1-ap).valuation(p)
			if e>0:
				Phi=Phi.change_precision(M-e)
				print "change precision to",M-e
		else:
			q=2
			v=self.is_Tq_eigen(q,p,M)
			assert v[0],"not eigen at q"
			aq=v[1]
			while (q<>p) and (aq-q^(k+1)-1)%(p^M)==0:
				q=next_prime(q)
				v=self.is_Tq_eigen(q,p,M)
				assert v[0],"not eigen at q"
				aq=v[1]
			Phi=(Phi.scale(q^(k+1)+1)-Phi.hecke(q)).scale(1/(q^(k+1)+1-aq))
			e=(q^(k+1)+1-aq).valuation(p)
			if e>0:
				Phi=Phi.change_precision(M-e)
				print "change precision to",M-e
		if verbose:
			print "Iterating U_p"
		Psi=Phi.hecke(p).scale(1/ap)
		err=(Psi-Phi).valuation()
		Phi=Psi
		while err<Infinity:
			if (Phi.valuation()>=s) and (s>0):
				Phi=Phi.scale(1/p^s)
				Phi=Phi.change_precision(Phi.num_moments()-s).normalize()
				print "unscaling by p^",s
				s=Infinity
			Psi=Phi.hecke(p).scale(1/ap)
			err=(Psi-Phi).valuation()
			if verbose:
				print "error is zero modulo p^",err
			Phi=Psi
		return Phi.normalize()

	"""self is a modular symbol taking values in Sym^k(K^2), where K is a finite extension of Q, psi is a map from the K to Qp and the below function 'map' applies psi to all polynomial coefficients and then lifts them to QQ"""
	def map(self,psi):
		v=[self.data[j].map(psi) for j in range(len(self.data))]
		return modsym_symk(self.level,v,self.manin)

	def coerce_to_Qp(self,p,M):
		"""If K is the base_ring of self, this function takes all maps K-->Q_p and applies them to self return a vector of <modular symbol,map:K-->Q_p> as map varies over all such maps.  M is the accuracy"""
		K=self.base_ring()
		f=K.defining_polynomial()	
		R.<x>=PolynomialRing(pAdicField(p,M+10))
		v=R(f).roots()
		if len(v)==0:
			print "No coercion possible -- no prime over p has degree 1"
			return []
		else:
			ans=[]
			for j in range(len(v)):
				root=v[j][0]
				psi=K.hom([root],pAdicField(p,M))
				ans=ans+[[self.map(psi),psi]]
		return ans

	def lift(self,p,ap):
		"""M.Greenberg method of lifting one step at a time -- slower in these programs because of how I optimized things"""
		k=self.weight()
		ans=[]
		for a in range(self.ngens()):
			v=[]
			for j in range(k+1):
				v=v+[phi.data[a].coef(j)]
			v=dist(p,k,vector(v+[0]))
			ans=ans+[v]	
		new=modsym_dist(self.level,ans,self.manin)
		return new.hecke(p).scale(1/ap).normalize()


def form_modsym_from_elliptic_curve(E):
	"""returns the modular symbol (in the sense of Stevens) associated to E"""
	ms_plus = E.modular_symbol()
	ms_minus = E.modular_symbol(sign = -1)
	N=E.conductor()
	manin=manin_relations(N)
	v=[]
	R.<X,Y>=PolynomialRing(QQ,2)
	for j in range(0,len(manin.gens)):
		rj=manin.gens[j]
		g=manin.mats[rj]
		a=g[0,0]
		b=g[0,1]
		c=g[1,0]
		d=g[1,1]
		if c<>0:
			a1=ms_plus(a/c)+ms_minus(a/c)
		else:
			a1=ms_plus(oo)+ms_minus(oo)
		if d<>0:
			a2=ms_plus(b/d)+ms_minus(b/d)
		else:
			a2=ms_plus(oo)+ms_minus(oo)
		v=v+[symk(0,R(a1))-symk(0,R(a2))]
	return modsym_symk(N,v,manin)

def form_modsym_from_decomposition(A):
	"""A is a piece of the result from a command like ModularSymbols(---).decomposition()"""
	M=A.ambient_module()
	N=A.level()
	k=A.weight()
	manin=manin_relations(N)

	w=A.dual_eigenvector()
	K=w.parent().base_field()
	v=[]
	R.<X,Y>=PolynomialRing(K,2)
	for j in range(0,len(manin.gens)):
		rj=manin.gens[j]
		g=manin.mats[rj]
		a=g[0,0]
		b=g[0,1]
		c=g[1,0]
		d=g[1,1]
		ans=0
		if c<>0:
			r1=a/c
		else:
			r1=oo
		if d<>0:
			r2=b/d
		else:
			r2=oo
		for j in range(k-1):
			coef=w.dot_product(M.modular_symbol([j,r1,r2]).element())
			ans=ans+X^j*Y^(k-2-j)*binomial(k-2,j)*coef
		v=v+[symk(k-2,ans,K)]
	return modsym_symk(N,v,manin)



