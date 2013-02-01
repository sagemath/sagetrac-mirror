class modsym_dist(modsym):
	def ms(self):
		"""demotes to a regular modular symbol"""
		return modsym(self.level,self.data,self.manin)

	def p(self):
		"""returns the underlying prime"""
		return self.data[0].p

	def weight(self):
		"""returns the underlying weight"""		
		return self.data[0].weight

	def num_moments(self):
		"""returns the number of moments of each value of self"""
		return self.data[0].num_moments()
	def is_zero(self):
		"""Return whether all of self's moments are zero."""
		for d in self.data:
			if not d.is_zero():
				return False
		return True

	def eval(self,A):
		"""here A is a 2x2 matrix and this function returns self evaluated and the divisor attached to A = A(\infty)-A(0)"""
		ans=self.ms().eval(A)
		return ans.normalize()

	def specialize(self):
		"""returns the underlying classical symbol of weight k -- i.e. applies the canonical map D_k --> Sym^k to all values of self"""
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].specialize()]
		return modsym_symk(self.level,v,self.manin)
	
	def valuation(self):
		"""returns the exponent of the highest power of p which divides all moments of all values of self"""
		return min([self.data[j].valuation() for j in range(0,len(self.data))])

	def normalize(self):
		"""normalized every values of self -- e.g. reduces each values j-th moment modulo p^(N-j)"""
		if self.valuation()>=0:
			v=[]
			for j in range(0,len(self.data)):
				v=v+[self.data[j].normalize()]
			ans=modsym_dist(self.level,v,self.manin,self.full_data)
			ans.normalize_full_data()
			
			return ans
		else:
			return self

	def normalize_aws(self):
		"""normalized every values of self -- e.g. reduces each values j-th moment modulo p^ceil((N-j)*(p-2)/(p-1)"""
		assert self.valuation()>=0, "moments are not integral"
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].aws_normalize()]
		ans=modsym_dist(self.level,v,self.manin,self.full_data)
		ans.normalize_full_data()
		
		return ans

	def change_precision(self,M):
		"""only holds on to M moments of each value of self"""
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].change_precision(M)]
		return modsym_dist(self.level,v,self.manin)

	def lift(self,phi,ap):
		"""Greenberg trick of lifting and applying U_p --- phi is the (exact) classical symbol"""
		v=[]
		for a in range(self.ngens()):
			v=v+[self.data[a].lift()]
		Phi=modsym_dist(self.level,v,self.manin)
		k=self.weight()
		for a in range(self.ngens()):
			for j in range(k+1):
				Phi.data[a].moments[j]=(phi.data[a].coef(j))%(p^(Phi.num_moments()))
		return Phi.hecke(self.p()).scale(1/ap).normalize()



	##  This function forms a family of OMSs which specialize to the given OMS
	def lift_to_modsym_dist_fam(self,w):
		p = self.p()
		M = self.num_moments()
		k = self.weight()
		r = self.weight() % (p-1)
		chi = self.data[0].char()
		deg = M
		v = []

		## this loop runs through each generator and lifts the value of self on that generator to D
		for j in range(1,len(self.manin.gens)):
			rj = self.manin.gens[j]
			if (self.manin.twotor.count(rj) == 0) and (self.manin.threetor.count(rj) == 0):
				v = v + [self.data[j].lift_to_dist_fam(deg,r,w)]
			elif (self.manin.twotor.count(rj) != 0):
				## Case of two torsion (See [PS] section 4.1)
				gam = self.manin.gen_rel_mat(j)
				mu = self.data[j].lift_to_dist_fam(deg,r,w)
				v = v + [(mu - mu.act_right(gam)).scale(1/2)]
			else:
				## Case of three torsion (See [PS] section 4.1)	
				gam = self.manin.gen_rel_mat(j)
				mu = self.data[j].lift_to_dist_fam(deg,r,w)
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

		## We now need to make some adjustment of Phi(D_i) to make t have total measure 0

		j = 1
		rj = self.manin.gens[j]
		gam = self.manin.gen_rel_mat(j)
		a = gam[0,0]
		c = gam[1,0]

		while (j < len(self.manin.gens)-1) and ((self.manin.twotor.count(rj) != 0) or (self.manin.threetor.count(rj) != 0) or (r == 0 and c == 0) or (r != 0 and (chi(a)*(a^r))%p == 1)):
			j = j + 1
			rj = self.manin.gens[j]
			gam = self.manin.gen_rel_mat(j)
			a = gam[0,0]
			c = gam[1,0]

		assert j < len(self.manin.gens) - 1, "Everything is 2 or 3 torsion -- or other problems!  NOT YET IMPLEMENTED IN THIS CASE"

		gam = self.manin.gen_rel_mat(j)

		a = gam[0,0]
		c = gam[1,0]
		chi = self.data[0].char()
		K = aut(p,deg,M,a,c,r,chi,w)
		K0 = K[0]  ## K0 is the coefficient of z^0 in K
		K1 = K[1]  ## K1 is the coefficient of z^1 in K
		t0 = t.moment(0)
		T = PowerSeriesRing(QQ,'ww',M+1)
		err = self.zero_elt().lift_to_dist_fam(deg,r,w)

		if r != 0:
			## The following code simply computes -t0/(K0-1)
			temp = T(T(-t0)/T(K[0]-1))
			temp = temp.truncate(deg)
			R = w.parent()
			temp = R(temp.padded_list())

			print "The result will be a lifting modulo {0}^{1}".format(p, valuation(temp.substitute(w=((1+p)^k-1)/p),p))
			err.moments[0] = temp
		else:
			## The following code simply computes -t0/(K1)
			temp=T(T(-t0)/T(K1))
			temp=temp.truncate(deg)
			R = w.parent()
			temp=R(temp.padded_list())
			print "The result will be a lifting modulo {0}^{1}".format(p, valuation(temp.substitute(w=((1+p)^k-1)/p),p))
			err.moments[1] = temp

		v[j-1] = v[j-1] + err
		t = t + err.act_right(gam)-err
		print "Total measure of t: {0} (should be zero)".format(t.normalize().moment(0))

		mu = t.solve_diff_eqn()
		print "Verify solution of difference equation: {0} (should be all zeroes)".format((mu.act_right(Matrix(2,2,[1,1,0,1]))-mu-t).normalize())

		v = [mu.scale(-1)] + v
	
		Phis = modsym_dist_fam(self.level,v,self.manin)

		return Phis

#	def check_loop(self):
#		return self.ms().check_loop().normalize()

	def is_Tq_eigen(self,q):
		p = self.data[0].p
		M = self.data[0].num_moments()

		a = 0
		m = 0
		done = false
		while (self.data[a].moment(m) % (p^M) == 0) and (not done):
			if a < self.ngens() - 1:
				a = a + 1
			else:
				if m < self.num_moments() - 1:
					m = m + 1
				else:
					done = True
		if done:	#Here the symbol is zero
			return [False, None, None]
		else:
			Phiq = self.hecke(q)
			c = Phiq.data[a].moment(m)/self.data[a].moment(m)
			checker = (Phiq - self.scale(c)).is_zero()
			if checker:
				return [True, c % (p^(M-m-self.valuation())), M-m-self.valuation()]
			return [False, None, None]
	
	def qexp(self, prec = 6):
		"""
		qexp of modsym_dist
		prec - a positive integer, the q-expansion will be computed up to O(q^prec)
		
		"""
		
		a_ells = self.aplist(prec - 1)
		a_ns = aplist_to_anlist(None, a_ells, prec - 1, data = [self.level, self.weight()])
		print a_ns
		#a_ns = [a_n % (self.p() ** self.num_moments()) for a_n in a_ns]
		R = PowerSeriesRing(QQ, 'q')
		a_ns = [0] + a_ns
		return R(a_ns)
			
	def vector_of_total_measures(self):
		"""returns the vector comprising of the total measure of each distribution defining Phi"""
		v=[]
		for j in range(0,self.ngens()):
			v=v+[self.data[j].moments[0]]
		return v

	def vector_of_total_measures_mod_pn(self,n):
		"""returns the vector comprising of the total measure of each distribution defining Phi"""
		v=[]
		R=Integers(self.p()^n)
		for j in range(0,self.ngens()):
			v=v+[R(self.data[j].moments[0])]
		return v


	def project_to_ordinary_subspace(self):
		"""Iterates U_p M+1 times where M is the number of moments"""
		p = self.p()
		Phi = self
		for r in range(self.num_moments()*2+2):
			Phi = Phi.hecke(p)

		return Phi
	
def random_OMS(N,p,k,M,char=None):
	"""Returns a random OMS with tame level N, prime p, weight k, and M moments"""
	## There must be a problem here with that +1 -- should be variable depending on a c of some matrix
	## We'll need to divide by some power of p and so we add extra accuracy here.
	if k != 0:
		MM = M + valuation(k,p) + 1 + floor(log(M)/log(p)) 
	else:
		MM = M + floor(log(M)/log(p)) + 1

	manin = manin_relations(N*p)
	v = []

#this loop runs thru all of the generators (except (0)-(infty)) and randomly chooses a distribution to assign  to this generator (in the 2,3-torsion cases care is taken to satisfy the relevant relation)
	for j in range(1,len(manin.gens)):
		rj = manin.gens[j]
		if (manin.twotor.count(rj) != 0):
			gam = manin.gen_rel_mat(j)
			mu=random_dist(p,k,MM,char)
			v=v+[(mu-mu.act_right(gam))]
		elif (manin.threetor.count(rj)!=0):
			gam = manin.gen_rel_mat(j)
			mu=random_dist(p,k,MM,char)
			v=v+[(mu.scale(2)-mu.act_right(gam)-mu.act_right(gam^2))]
		else:
			v=v+[random_dist(p,k,MM,char)]

#now we compute nu_infty of Prop 5.1 of [PS1]
	t=v[0].zero()
	for j in range(1,len(manin.gens)):
		rj = manin.gens[j]
		if (manin.twotor.count(rj) == 0) and (manin.threetor.count(rj) == 0):
			t = t + v[j-1].act_right(manin.gen_rel_mat(j)) - v[j-1]
		else:
			t = t - v[j-1]
		
## If k = 0, then t has total measure zero.  However, this is not true when k != 0  
## (unlike Prop 5.1 of [PS1] this is not a lift of classical symbol).  
## So instead we simply add (const)*mu_1 to some (non-torsion) v[j] to fix this
## here since (mu_1 |_k ([a,b,c,d]-1))(trival char) = chi(a) k a^{k-1} c , 
## we take the constant to be minus the total measure of t divided by (chi(a) k a^{k-1} c)

	if k != 0:
		j = 1
		rj = manin.gens[j]
		while (j < len(manin.gens)-1) and ((manin.twotor.count(rj) != 0) or (manin.threetor.count(rj) != 0)):
			j = j + 1
			rj = manin.gens[j]
			print j
		assert j < len(manin.gens) - 1, "Everything is 2 or 3 torsion!  NOT YET IMPLEMENTED IN THIS CASE"

		gam = manin.gen_rel_mat(j)
		a = gam[0,0]
		c = gam[1,0]
		
		if char != None:
			chara = char(a)
		else:
			chara = 1
		err = -t.moments[0]/(chara*k*a^(k-1)*c)
		mu_1 = v[rj-1].zero()
		mu_1.moments[1] = err
		v[j-1] = v[j-1] + mu_1
		t = t + mu_1.act_right(gam) - mu_1
		
	mu = t.solve_diff_eqn()
	v = [mu.scale(-1)] + v
	Phi = modsym_dist(N*p,v,manin)
	
	Psi = Phi.change_precision(M)
	e = max(-Psi.valuation(),0)
	Phi = Phi.scale(p^e).change_precision(M)

	return Phi

## NOT YET REWRITTEN AS ABOVE
def random_OMS_char(N,p,k,chi,M):
	"""Returns a random OMS with tame level N, prime p, weight k, and M moments --- requires no 2 or 3-torsion"""
	manin=manin_relations(N*p)
	v=[]
	for j in range(1,len(manin.gens)):
		g=manin.gens[j]
		print g
		print manin.glue[g][1]
		gam=manin.glue[g][1]
		a=gam[0,0]
		c=gam[1,0]
		if manin.twotor.count(g)==0:
			mu=random_dist_char(p,k,chi,M).zero()
			mu.moments[0]=c
			mu.moments[1]=(chi(a)-a)*a/chi(a)  #formula to force difference equation to work 
			v=v+[mu]
		else:
			v=v+[random_dist_char(p,k,chi,M).zero()]
#			rj=manin.twotor.index(g)
#			gam=manin.twotorrels[rj]
#			mu=random_dist_char(p,k,chi,M)
#			v=v+[(mu.act_right(gam)-mu).scale(1/2)]
	t=v[0].zero()
	for j in range(2,len(manin.rels)):
		R=manin.rels[j]
		if len(R)==1:
			if R[0][0]==1:
				rj=manin.gens.index(j)
				t=t+v[rj-1]
				print v[rj-1]
			else:
				index=R[0][2]
				rj=manin.gens.index(index)
				mu=v[rj-1]
				print mu.act_right(R[0][1]).scale(R[0][0]).normalize()
				t=t+mu.act_right(R[0][1]).scale(R[0][0])
	t=t.normalize()
	mu=t.solve_diff_eqn()
	v=[mu.scale(-1)]+v

	return modsym_dist(N*p,v,manin)	

def random_ordinary_OMS(N,p,k,M,char=None):
	Phi = random_OMS(N,p,k,M,char)
#	print "first check",Phi.check_loop().normalize()
	Phi = Phi.project_to_ordinary_subspace()
#	print "second check",Phi.check_loop().normalize()
	Phi = Phi.change_precision(M)

	return Phi

## U_p is iterated r times and the resulting list of vectors of moments is returned
def hecke_span(Phi,q,r):
	p = Phi.p()
	M = Phi.num_moments()

	Psi = Phi
	list = []
	list = list + [Psi.vector_of_total_measures()]
	for j in range(r):
		Psi = Psi.hecke(q)
		list = list + [Psi.vector_of_total_measures()]
	
	A=Matrix(list)
	
	return list,A,(A.echelon_form())%(p^M)

def aplist_to_anlist(selfe, aps, max_n, data = None):
    if data is None:
        eps = selfe.character()
        k = selfe.weight()
        N = selfe.level()
    else:
        if data[0] in ZZ:
            N = data[0]
            eps = DirichletGroup(N)[0]	#trivial character
        else:
            eps = data[0]
            N = eps.level()
        k = data[1]
    anlist = [0] * max_n
    anlist[0] = 1
    ps = prime_range(max_n + 1)
    if len(aps) < len(ps):
        raise ValueError("Not enough ap's given. Must reach at least {0}.".format(max_n))
    ap_to_the_rs = {}
    for i in range(len(ps)):
        p = ps[i]
        if p > max_n:
            break
        powers = [aps[i]]
        cur_pow = p ** 2
        j = 2
        while cur_pow <= max_n:
            if ZZ(p).divides(N):
                powers.append(aps[i] * powers[-1])
            else:
                if j == 2:
                    powers.append(aps[i] ** 2 - eps(p) * (p ** (k-1)))
                else:
                    powers.append(aps[i] * powers[-1] - eps(p) * (p ** (k-1)) * powers[-2])
            cur_pow = cur_pow * p
            j += 1
        ap_to_the_rs[p] = powers
    
    for n in range(2, max_n + 1):
        fs = ZZ(n).factor()
        an = 1
        for f in fs:
            an = an * ap_to_the_rs[f[0]][f[1] - 1]
        anlist[n - 1] = an
    
    return anlist