class modsym_dist_fam(modsym):
	def ms(self):
		"""demotes to a regular modular symbol"""
		return modsym(self.level,self.data,self.manin)

	def num_moments(self):
		return self.data[0].num_moments()

	def p(self):
		return self.data[0].p

	def deg(self):
		return self.data[0].deg

	## This function returns a number between 0 and p-2 which represents which 	
	## disc in weight the family is supported on -- the integer i corresponds to the
	## disc with characters with finite order part omega^i
	def disc(self):
		return self.data[0].disc()

	def change_deg(self,new_deg):
		v=[self.data[r].change_deg(new_deg) for r in range(self.ngens())]
		return modsym_dist_fam(self.level,v,self.manin)

	def specialize(self,k):
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].specialize(k)]
		return modsym_dist_aws(self.level,v,self.manin)
	
	def valuation(self):
		#print [self.data[j].valuation() for j in range(len(self.data))]
		return min([self.data[j].valuation() for j in range(len(self.data))])

	def normalize(self):
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].normalize()]
		return modsym_dist_fam(self.level,v,self.manin)
	
	def change_precision(self,M):
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].change_precision(M)]
		return modsym_dist_fam(self.level,v,self.manin)
	
	def is_zero(self):
		"""Return true if all of selfs moments are zero."""
		for d in self.data:
			if not d.is_zero():
				return False
		return True

	## This procedure tries to find a power series c(w) such that 
	##      self | T_q = c(w) self
	## Returns a triple consisting of a boolean and if true c(w) and the precision (if false, None and None)
	def is_Tq_eigen(self,q,verbose=False):
		p = self.p()
		M = self.num_moments()
		R = self.data[0].moment(0).parent()
		T = PowerSeriesRing(QQ,'y')

		a = 0
		done = false
		while (self.data[a].moment(0) == 0) and (not done):
			if a < self.ngens() - 1:
				a = a + 1
			else:
				done = True
		if done:	
			print "All of the total measures are zero!"
			return [False, None, None]
		else:
			Phiq = self.hecke(q)
			aq = R(T(T(Phiq.data[a].moment(0))/T(self.data[a].moment(0))).padded_list())

			checker = (self.scale(aq) - Phiq).normalize().is_zero()
			if checker:
				return [True, ps_normalize(aq,p,M-self.valuation()), [M-self.valuation(), self.deg()]]

			if verbose:
				print "The difference Phi | T_q - (potential eval) * Phi is:",(self.scale(aq) - Phiq).normalize()

			return [False, None, None]
	

	def vector_of_total_measures(self):
		"""returns the vector comprising of the total measure of each distribution defining Phi"""
		v=[]
		for j in range(0,self.ngens()):
			v=v+[self.data[j].moments[0]]
		return v


#######################################################################################################################
##  This function produces a random family of OMSs.
##
##  p -- prime
##  N -- tame level
##  char -- character of conductor dividing N (or maybe Np?)
##  M -- number of moments
##  r -- integer between 0 and p-2 indicating which disc on weight space we are working
##  w -- variable which is just carried around because I had trouble with symbols being defined over different rings
#####################################################################################################################
def random_OMS_fam(p,N,M,deg,r,w,char=None,TestTorsion=False):
	if char == None:
		char = DirichletGroup(1,QQ).0
	manin = manin_relations(N*p)
	t2 = 0
	t3 = 0
	if (p == 2) and (len(manin.twotor) > 0):
		t2 = 1
	if (p == 3) and (len(manin.threetor) > 0):
		t3 = 1
	MM = M + t2 + t3 + 1 + floor(log(M)/log(p))

	v = []
	## this loop runs through each generator (different from D_infty) and picks a random value for that generator
	for j in range(1,len(manin.gens)):
		rj = manin.gens[j]
		mus = random_dist_fam(p,M,deg,r,w,char)
		if (manin.twotor.count(rj) == 0) and (manin.threetor.count(rj) == 0):
			v = v + [mus]
		elif (manin.twotor.count(rj) != 0):
			## Case of two torsion (See [PS] section 4.1)
			gam = manin.gen_rel_mat(j)
			v = v + [(mus - mus.act_right(gam)).scale(1/2)]
		else:
			## Case of three torsion (See [PS] section 4.1)	
			gam = manin.gen_rel_mat(j)
			v = v + [(mus.scale(2) - mus.act_right(gam) - mus.act_right(gam^2)).scale(1/3)]

	t = v[0].zero()
	## This loops adds up around the boundary of fundamental domain except the two verticle lines
	for j in range(1,len(manin.gens)):
		rj = manin.gens[j]
		if (manin.twotor.count(rj) == 0) and (manin.threetor.count(rj) == 0):
			t = t + v[j-1].act_right(manin.gen_rel_mat(j)) - v[j-1]
		else:
			t = t - v[j-1]

	## t now should be sum Phi(D_i) | (gamma_i - 1) - sum Phi(D'_i) - sum Phi(D''_i)
	## (Here I'm using the opposite sign convention of [PS1] regarding D'_i and D''_i)

	## We now need to make some adjustment of Phi(D_i) to make t have total measure 0

	j = 1
	rj = manin.gens[j]
	while (j < len(manin.gens)-1) and ((manin.twotor.count(rj) != 0) or (manin.threetor.count(rj) != 0)):
		j = j + 1
		rj = manin.gens[j]
	if j < len(manin.gens) - 1:
		print "Found a non-torsion generator"

	gam = manin.gen_rel_mat(j)
	a = gam[0,0]
	c = gam[1,0]
	K = aut(p,deg,M,a,c,r,char,w)
	K0 = K[0]  ## K0 is the coefficient of z^0 in K
	K1 = K[1]  ## K1 is the coefficient of z^1 in K
	t0 = t.moment(0)
	T = PowerSeriesRing(QQ,'ww')
	err = mus.zero()

	if r != 0:
		## The following code simply computes -t0/(K0-1)
		temp = T(T(-t0)/T(K[0]-1))
		temp = temp.truncate(deg)
		R = w.parent()
		temp = R(temp.padded_list())

		err.moments[0] = temp
	else:
		## The following code simply computes -t0/(K1)
		temp=T(T(-t0)/T(K1))

		temp=temp.truncate(deg)
		R = w.parent()
		temp=R(temp.padded_list())

		err.moments[1] = temp

	if manin.twotor.count(rj) > 0 or manin.threetor.count(rj) >0:
		print "All generators are two and three torsion"
		v[j-1] = v[j-1] + err.act_right(gam) - err
		t = t + err.act_right(gam) - err
	else:
		v[j-1] =v[j-1] + err
		t = t + err.act_right(gam) - err

	if TestTorsion:
		print "Total Measure", t.moments[0],t.moments[0].valuation(p)
		print "Is t zero?", t.is_zero()
		print "is two/three tor?", manin.twotor.count(rj) > 0, manin.threetor.count(rj) >0



	#v[j-1] = v[j-1] + err
	#t = t + err.act_right(gam)-err

	mus = t.solve_diff_eqn()

	v = [mus.scale(-1)] + v

	Phis = modsym_dist_fam(N*p,v,manin)

	Psis = Phis.change_precision(M)
	e = max(-Psis.valuation(),0)
	Phis = Phis.scale(p^e).change_precision(M)

	return Phis

def random_OMS_fam_new(p,N,M,deg,r,w,char=None):
	if char == None:
		char = DirichletGroup(1,QQ).0

	manin = manin_relations(N*p)

	t2 = 0
	t3 = 0
	if (p == 2) and (len(manin.twotor) > 0):
		t2 = 1
	if (p == 3) and (len(manin.threetor) > 0):
		t3 = 1
	MM = M + t2 + t3 + 1 + floor(log(M)/log(p))

	v = []
	first_time = True
	## this loop runs through each generator (different from D_infty) and picks a random value for that generator
	## Also, it adds up around the boundary of fundamental domain except the two vertical lines
	for j in range(1,len(manin.gens)):
		rj = manin.gens[j]
		mus = random_dist_fam(p,M,deg,r,w,char)
		if first_time:
			t = mus.zero()
			first_time = False
		if (manin.twotor.count(rj) == 0) and (manin.threetor.count(rj) == 0):
			v = v + [mus]
			t = t + v[j-1].act_right(manin.gen_rel_mat(j)) - v[j-1]
		elif (manin.twotor.count(rj) != 0):
			## Case of two torsion (See [PS] section 4.1)
			gam = manin.gen_rel_mat(j)
			v = v + [(mus - mus.act_right(gam)).scale(1/2)]
			t = t - v[j-1]
		else:
			## Case of three torsion (See [PS] section 4.1)	
			gam = manin.gen_rel_mat(j)
			v = v + [(mus.scale(2) - mus.act_right(gam) - mus.act_right(gam^2)).scale(1/3)]
			t = t - v[j-1]

	## t now should be sum Phi(D_i) | (gamma_i - 1) - sum Phi(D'_i) - sum Phi(D''_i)
	## (Here I'm using the opposite sign convention of [PS1] regarding D'_i and D''_i)

	## We now need to make some adjustment of Phi(D_i) to make t have total measure 0

	j = 1
	rj = manin.gens[j]
	while (j < len(manin.gens)-1) and ((manin.twotor.count(rj) != 0) or (manin.threetor.count(rj) != 0)):
		j = j + 1
		rj = manin.gens[j]
	assert j < len(manin.gens) - 1, "Everything is 2 or 3 torsion!  NOT YET IMPLEMENTED IN THIS CASE"

	gam = manin.gen_rel_mat(j)
	a = gam[0,0]
	c = gam[1,0]
	K = aut(p,deg,M,a,c,r,char,w)
	K0 = K[0]  ## K0 is the coefficient of z^0 in K
	K1 = K[1]  ## K1 is the coefficient of z^1 in K
	t0 = t.moment(0)
	T = PowerSeriesRing(QQ,'ww')
	err = mus.zero()

	if r != 0:
		## The following code simply computes -t0/(K0-1)
		temp = T(T(-t0)/T(K[0]-1))
		temp = temp.truncate(deg)
		R = w.parent()
		temp = R(temp.padded_list())

		err.moments[0] = temp
	else:
		## The following code simply computes -t0/(K1)
		temp=T(T(-t0)/T(K1))

		temp=temp.truncate(deg)
		R = w.parent()
		temp=R(temp.padded_list())

		err.moments[1] = temp

	v[j-1] = v[j-1] + err
	t = t + err.act_right(gam)-err

	mus = t.solve_diff_eqn()

	v = [mus.scale(-1)] + v

	Phis = modsym_dist_fam(N*p,v,manin)

	Psis = Phis.change_precision(M)
	e = max(-Psis.valuation(),0)
	Phis = Phis.scale(p^e).change_precision(M)

	return Phis

def aplist_to_anlist_fam(selfe, aps, max_n, data = None):
    """
    To use this function: suppose you have a family of eigensymbols Phis (i.e. a modsym_dist_fam object)
    and you want to know its formal Fourier coefficients for n=1 to max_n (inclusively) (warning: in the future,
    this will start at n=0). Here is what parameters should be passed:
        
        INPUT::
            
            selfe: set this to Phis
            aps: set this to the output of Phis.aplist(max_n)
            max_n: highest n such that you want a_n
            data: set this to the pair [N*p, weight], where N is the tame level, p is the prime, weight is the
                weight (normalized so that elliptic curves correspond to weight 2).
        
        OUPTUT::
            
            The list [a_1, a_2, ..., a_n].
        
        WARNING::
            
            In a future implementation the output will begin with a_0.
    """
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
    p_power = None
    for i in range(len(ps)):
        p = ps[i]
        if p > max_n:
            break
        powers = [aps[i]]
        cur_pow = p ** 2
        j = 2
	R = aps[0].parent()
        while cur_pow <= max_n:
            if ZZ(p).divides(N):
                powers.append(aps[i] * powers[-1])
            else:
                if p_power is None:
                    p_power = ell_power_fam(selfe.p(), p, k, 2 * selfe.num_moments(), 2 * selfe.deg(), selfe.data[0].moments[0].variables()[0])
                if j == 2:
                    powers.append(aps[i] ** 2 - R(eps(p) * p_power))
                else:
                    powers.append(aps[i] * powers[-1] - R(eps(p) * p_power) * powers[-2])
            cur_pow = cur_pow * p
            j += 1
        ap_to_the_rs[p] = powers
        p_power = None
    
    for n in range(2, max_n + 1):
        fs = ZZ(n).factor()
        an = 1
        for f in fs:
            an = an * ap_to_the_rs[f[0]][f[1] - 1]
        anlist[n - 1] = an
    
    return anlist
