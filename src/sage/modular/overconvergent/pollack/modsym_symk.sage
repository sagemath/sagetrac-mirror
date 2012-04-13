from sage.rings.arith import binomial
from sage.rings.padics.all import pAdicField
from sage.rings.integer_ring import ZZ
from sage.matrix.all import Matrix

class modsym_symk(modsym):

    def ms(self):
        r"""
        Demotes to a regular modular symbol
        """
        return modsym(self.level,self.data,self.manin)
        
    def weight(self):
        r"""
        Returns the weight of any value of self
        """
        return self.data[0].weight

    def base_ring(self):
        r"""
        Returns the base ring of self
        """
        return self.data[0].base_ring

    def valuation(self,p):
        r"""
        Returns the valuation of self at `p` -- i.e. the exponent of the 
        largest power of `p` which divides all of the coefficients of all 
        values of self
        """
        v = self.data
        return min([v[j].valuation(p) for j in range(len(v))])

    def p_stabilize(self,p,alpha):
        r"""
        `p`-stablizes self to form an eigensymbol for `U_p` with eigenvalue alpha
        """
        N = self.level
        assert N%p<>0, "The level isn't prime to p"
        pp = Matrix(ZZ,[[p,0],[0,1]])
        manin = manin_relations(N*p)
        v = []
        for j in range(len(manin.gens)):
            rj = manin.gens[j]
            v = v+[self.eval(manin.mats[rj])-self.eval(pp*manin.mats[rj]).act_right(pp).scale(1/alpha)]
        return modsym_symk(N*p,v,manin)

    def p_stabilize_ordinary(self,p,ap,M):
        r"""
        Returns the unique `p`-ordinary `p`-stabilization of self
        """
        N = self.level
        k = self.data[0].weight
        assert N%p<>0, "The level isn't prime to p"
        assert ap%p<>0, "Not ordinary!"

        # makes alpha the unit root of Hecke poly
        R = pAdicField(p,M)['y'] 
        y = R.gen()
        f = y**2-ap*y+p**(k+1)
        v = f.roots()
        if ZZ(v[0][0])%p<>0:
            alpha = ZZ(v[0][0])
        else:
            alpha = ZZ(v[1][0])
        if self.full_data == 0:
            self.compute_full_data_from_gen_data()
        return self.p_stabilize(p,alpha)

    def p_stabilize_critical(self,p,ap,M):
        r"""
        """
        N = self.level
        assert N%p<>0, "The level isn't prime to p"
        assert ap%p<>0, "Not ordinary!"

        # makes alpha the non-unit root of Hecke poly
        R = pAdicField(p,M)['y'] 
        y = R.gen()
        f = y**2-ap*y+p
        v = f.roots()
        if ZZ(v[0][0])%p == 0:
            alpha = ZZ(v[0][0])
        else:
            alpha = ZZ(v[1][0])
        if self.full_data == 0:
            self.compute_full_data_from_gen_data()
        return self.p_stabilize(p,alpha)

    def is_Tq_eigen(self,q,p,M):
        r"""
        Determines if self is an eigenvector for `T_q` modulo `p^M`
        """
        selfq = self.hecke(q)
        r = 0
        while (self.data[r].coef(0)==0) and (r<=self.ngens()):
            r = r+1
        if r > self.ngens():
        #all coefficients of Y^k are zero which I think forces it 
        #not to be eigen
            return False
        else:   
            c = selfq.data[r].coef(0)/self.data[r].coef(0)
            val = (selfq-self.scale(c)).valuation(p)
            if val >= M:
                return True,c
            else:
                return False
                
    def lift_to_OMS(self,p,M):
        r"""
        Returns a (`p`-adic) overconvergent modular symbol with `M` moments 
        which lifts self up to an Eisenstein error
        """
        v = []
        # this loop runs through each generator and lifts the value of 
        # self on that generator to D
        
        for j in range(1,len(self.manin.gens)):
            g = self.manin.gens[j]
            if (self.manin.twotor.count(g)==0) and (self.manin.threetor.count(g)==0):
                #not two or three torsion
                v = v + [self.data[j].lift_to_dist(p,M)]
            else:
                if (self.manin.twotor.count(g)<>0):
                    #case of two torsion (See [PS] section 4.1)
                    rj = self.manin.twotor.index(g)
                    gam = self.manin.twotorrels[rj]
                    mu = self.data[j].lift_to_dist(p,M)
                    v = v+ [(mu.act_right(gam)-mu).scale(ZZ(1)/ZZ(2))]
                else:
                    #case of three torsion (See [PS] section 4.1)       
                    rj = self.manin.threetor.index(g)
                    gam = self.manin.threetorrels[rj]
                    mu = self.data[j].lift_to_dist(p,M)
                    v = v+[(mu.scale(2)-mu.act_right(gam)-mu.act_right(gam**2)).scale(ZZ(1)/ZZ(3))]
        
        t = v[0].zero()
                
        # this loops adds up around the boundary of fundamental domain 
        # except the two verticle lines
        for j in range(2,len(self.manin.rels)):
            R = self.manin.rels[j]
            if len(R) == 1:
                if R[0][0] == 1:
                    rj = self.manin.gens.index(j)
                    t = t + self.data[rj].lift_to_dist(p,M)
                else:
                    if R[0][1]<>Id:
                    #rules out extra three torsion terms
                        index = R[0][2]
                        rj = self.manin.gens.index(index)
                        t = t+self.data[rj].lift_to_dist(p,M).act_right(R[0][1]).scale(R[0][0])
        mu = t.solve_diff_eqn()
        v = [mu] + v
        return modsym_dist(self.level,v,self.manin)     

    def lift_to_OMS_eigen(self,p,M,verbose=True):
        r"""
        Returns Hecke-eigensymbol OMS lifting `phi` -- `phi` must be a 
        `p`-ordinary eigensymbol
        """
        v = self.is_Tq_eigen(p,p,M)
        assert v[0], "not eigen at p!"
        ap = v[1]
        k = self.weight()
        Phi = self.lift_to_OMS(p,M)
        s = -Phi.valuation()
        if s > 0:
            if verbose:
                print "Scaling by ",p,"^",s
            Phi = Phi.scale(p**(-Phi.valuation()))
        Phi = Phi.normalize()
        if verbose:
            print "Applying Hecke"
        Phi = Phi.hecke(p).scale(1/ap)
        if verbose:
            print "Killing eisenstein part"
        if (ap%(p**M))<>1:
            Phi = (Phi-Phi.hecke(p)).scale(1/(1-ap))
            e = (1-ap).valuation(p)
            if e > 0:
                Phi = Phi.change_precision(M-e)
                print "change precision to",M-e
        else:
            q = 2
            v = self.is_Tq_eigen(q,p,M)
            assert v[0],"not eigen at q"
            aq = v[1]
            while (q<>p) and (aq-q**(k+1)-1)%(p**M)==0:
                q = next_prime(q)
                v = self.is_Tq_eigen(q,p,M)
                assert v[0],"not eigen at q"
                aq = v[1]
            Phi = (Phi.scale(q**(k+1)+1)-Phi.hecke(q)).scale(1/(q**(k+1)+1-aq))
            e = (q**(k+1)+1-aq).valuation(p)
            if e > 0:
                Phi = Phi.change_precision(M-e)
                print "change precision to",M-e
        if verbose:
            print "Iterating U_p"
        Psi = Phi.hecke(p).scale(1/ap)
        err = (Psi-Phi).valuation()
        Phi = Psi
        while err < Infinity:
            if (Phi.valuation()>=s) and (s>0):
                Phi = Phi.scale(1/p**s)
                Phi = Phi.change_precision(Phi.num_moments()-s).normalize()
                print "unscaling by p^",s
                s = Infinity
            Psi = Phi.hecke(p).scale(1/ap)
            err = (Psi-Phi).valuation()
            if verbose:
                print "error is zero modulo p^",err
            Phi = Psi
        return Phi.normalize()
        
    ### self is a modular symbol taking values in Sym^k(K^2), where 
    ### K is a finite extension of Q, psi is a map from the K to Qp and 
    ### the below function 'map' applies psi to all polynomial 
    ### coefficients and then lifts them to QQ"""
        

    def map(self,psi):
        r"""
        """
        v = [self.data[j].map(psi) for j in range(len(self.data))]
        return modsym_symk(self.level,v,self.manin)

    def coerce_to_Qp(self,p,M):
        r"""
        If `K` is the base_ring of self, this function takes all maps 
        `K-->Q_p` and applies them to self return a vector of 
        <modular symbol,map: `K-->Q_p`> as map varies over all such maps.  
        `M` is the accuracy
        """
        K = self.base_ring()
        f = K.defining_polynomial()
        R = pAdicField(p,M+10)['x']
        x = R.gen()
        v = R(f).roots()
        if len(v) == 0:
            print "No coercion possible -- no prime over p has degree 1"
            return []
        else:
            ans = []
            for j in range(len(v)):
                root = v[j][0]
                psi = K.hom([root],pAdicField(p,M))
                ans = ans+[[self.map(psi),psi]]
        return ans

    def form_modsym_from_elliptic_curve(E):
        r"""
        Returns the modular symbol (in the sense of Stevens) associated to `E`
        """
        L = E.padic_lseries(3)  #should this be p?
        N = E.conductor()
        manin = manin_relations(N)
        v = []
        R = PolynomialRing(QQ,2,names='X,Y')
        X,Y = R.gens()
        for j in range(len(manin.gens)):
            rj = manin.gens[j]
            g = manin.mats[rj]
            a,b,c,d = g.list()
            if c != 0:
                a1 = L.modular_symbol(a/c,1)+L.modular_symbol(a/c,-1)
            else:
                a1 = L.modular_symbol(oo,1)+L.modular_symbol(oo,-1)
            if d != 0:
                a2 = L.modular_symbol(b/d,1)+L.modular_symbol(b/d,-1)
            else:
                a2 = L.modular_symbol(oo,1)+L.modular_symbol(oo,-1)
            v = v + [symk(0,R(a1)) - symk(0,R(a2))]
        return modsym_symk(N,v,manin)

    def form_modsym_from_decomposition(A):
        r"""
        `A` is a piece of the result from a command like 
        ModularSymbols(---).decomposition()
        """
        M = A.ambient_module()
        N = A.level()
        k = A.weight()
        manin = manin_relations(N)
        w = A.dual_eigenvector()
        K = w.parent().base_field()
        v = []
        R = PolynomialRing(K,2,names='X,Y')
        X,Y = R.gens()
        for j in range(0,len(manin.gens)):
            rj = manin.gens[j]
            g = manin.mats[rj]
            a,b,c,d = g.list()
            ans = 0
            if c<>0:
                r1 = a/c
            else:
                r1 = oo
            if d<>0:
                r2 = b/d
            else:
                r2 = oo
            for j in range(k-1):
                coef = w.dot_product(M.modular_symbol([j,r1,r2]).element())
                ans = ans+X**j*Y**(k-2-j)*binomial(k-2,j)*coef
            v = v+[symk(k-2,ans,K)]
        return modsym_symk(N,v,manin)

    def lift(self,p,ap):
        """
        M.Greenberg method of lifting one step at a time -- slower in these 
        programs because of how I optimized things
        """
        k = self.weight()
        ans = []
        for a in range(self.ngens()):
            #v = []
            #for j in range(k+1):
            #    v = v + [phi.data[a].coef(j)]
            v = [phi.data[a].coeff(j) for j in range(k+1)]
            v = dist(p,k,vector(v+[0]))
            ans = ans + [v]     
        new = modsym_dist(self.level,ans,self.manin)
        return new.hecke(p).scale(1/ap).normalize()
