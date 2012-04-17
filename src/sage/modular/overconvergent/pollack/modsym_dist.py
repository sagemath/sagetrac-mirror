
from sage.rings.integer_ring import ZZ

import modsym

class modsym_dist(modsym.modsym):
    
    def ms(self):
        r"""
	Demotes to a regular modular symbol
	"""
	return modsym(self.data(),self.manin())

    def p(self):
	r"""
	Returns the underlying prime
	"""
	return self.data()[0].p

    def weight(self):
	r"""
        Returns the underlying weight
        """		
        return self.data()[0].weight

    def num_moments(self):
        r"""
        Returns the number of moments of each value of self
        """
	return self.data()[0].num_moments()

    def eval(self,A):
        r"""
        Here `A` is a `2x2` matrix and this function returns self evaluated 
	and the divisor attached to `A = A(\infty)-A(0)`
        """
	ans = self.ms().eval(A)
	return ans.normalize()

    def specialize(self):
        r"""
        Returns the underlying classical symbol of weight `k` -- i.e. applies 
        the canonical map `D_k --> Sym^k` to all values of self
        """
	#v = []
	#for j in range(0,len(self.data())):
        #    v=v+[self.data()[j].specialize()]

	v = [j.specialize() for j in self.data()]
	return modsym_symk(self.level,v,self.manin)
	
    def valuation(self):
	r"""
        Returns the exponent of the highest power of `p` which divides all 
        moments of all values of self
        """
#	return min([self.data()[j].valuation() for j in range(len(self.data()))])
        return min([j.valuation() for j in self.data()])

    def normalize(self):
        r"""
        Normalizes every value of self -- e.g. reduces each value's `j`-th 
        moment modulo `p^(N-j)`
        """
        assert self.valuation()>=0, "moments are not integral"

        #v = []
        #for j in range(0,len(self.data())):
        #     v=v+[self.data()[j].normalize()]

        v = [j.normalize() for j in self.data()]
	ans = modsym_dist(v,self.manin(),self._full_data)
        ans.normalize_full_data()
        return ans

    def change_precision(self,M):
        r"""
        Only holds on to `M` moments of each value of self
        """
	#v=[]
	#for j in range(0,len(self.data())):
	#    v=v+[self.data()[j].change_precision(M)]

	v = [j.change_precision(M) for j in self.data()]
	return modsym_dist(v,self.manin())

    def lift(self,phi,ap):
	r"""
        Greenberg trick of lifting and applying `U_p` --- `\phi` is the 
        (exact) classical symbol
        """
#        v=[]
#	for a in range(self.ngens()):
#            v=v+[self.data()[a].lift()]

        v = [self.data()[a].lift() for a in range(self.ngens())]
	Phi = modsym_dist(v,self.manin())
	k = self.weight()
	for a in range(self.ngens()):
            for j in range(k+1):
		Phi.data[a].moments[j]=(phi.data[a].coef(j))%(p**(Phi.num_moments()))
	return Phi.hecke(self.p()).scale(1/ap).normalize()

    def is_Tq_eigen(Phi,q):
        r"""
	"""
	Phiq = Phi.hecke(q)
	c = Phiq.data[0].moment(0)/Phi.data[0].moment(0)
	M = Phi.data[0].num_moments()
	p = Phi.data[0].p
	#print Phiq-Phi.scale(c)
        return c%(p**M)
		
	
    def random_OMS(N,p,k,M):
        r"""
        Returns a random OMS with tame level `N`, prime `p`, weight `k`, and 
        `M` moments --- requires no `2` or `3`-torsion
        """
	manin = manin_relations(N*p)
	v = []
	for j in range(1,len(manin.gens)):
	    g = manin.gens[j]
	    if manin.twotor.count(g)==0:
                v = v+[random_dist(p,k,M)]
	    else:
	        rj = manin.twotor.index(g)
		gam = manin.twotorrels[rj]
	   	mu = random_dist(p,k,M)
		v = v+[(mu.act_right(gam)-mu).scale(ZZ(1)/ZZ(2))]
	t = v[0].zero()
	for j in range(2,len(manin.rels)):
	    R = manin.rels[j]
            if len(R) == 1:
		if R[0][0] == 1:
	            rj = manin.gens.index(j)
		    t = t+v[rj-1]
                else:
	            index = R[0][2]
	            rj = manin.gens.index(index)
		    mu = v[rj-1]
		    t = t+mu.act_right(R[0][1]).scale(R[0][0])
	v = [mu.scale(-1)]+v
	Phi = modsym_dist(v,manin)	
	return Phi

