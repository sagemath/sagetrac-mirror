from sage.rings.arith import binomial
from sage.rings.padics.all import pAdicField
from sage.rings.integer_ring import ZZ
from sage.matrix.all import Matrix

class modsym_symk(modsym):

    ### self is a modular symbol taking values in Sym^k(K^2), where
    ### K is a finite extension of Q

    def ms(self):
        r"""
        Demotes self from a modsym_symk type to a modsym type.

        INPUT:
        
            - ``phi`` - modsym_symk

        OUTPUT:

        The same modular symbol phi but now simply thought of as a modsym type.

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: type(phi)
        <class '__main__.modsym_symk'>
        sage: psi=phi.ms(); psi
        [-1/5, 3/2, -1/2]
        sage: type(psi)
        <class '__main__.modsym'>

        """
        return modsym(self.data(),self.manin())
        
    def weight(self):
        r"""
        Returns the weight of any value of self.

        INPUT:

            - ``phi`` - A Sym^k-valued modular symbols

        OUTPUT:

        The weight of a value of self.

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi.weight()
        0

        """
        return self.data(0).weight()

    def base_ring(self):
        r"""
        Returns the base ring of self

        INPUT:
   
            none
        
        OUTPUT:

        The ring containing the coefficients of the values of self.

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi.base_ring()
        Rational Field

        """
        return self.data(0).base_ring()

    def valuation(self,p):
        r"""
        Returns the valuation of self at `p`.

        Here the valuation if the exponent of the largest power of `p` 
        which divides all of the coefficients of all values of self.

        INPUT:
            - ``p`` - prime

        OUTPUT:

        The valuation of self at `p`

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi.valuation(2)
        -1
        sage: phi.valuation(3)
        0
        sage: phi.valuation(5)
        -1
        sage: phi.valuation(7)
        0

        """
        return min([self.data(j).valuation(p) for j in range(self.ngens())])

    def p_stabilize(self,p,alpha):
        r"""
        Returns the `p`-stablization of self to level `N*p` on which `U_p` acts by `alpha`.

        Note that alpha is p-adic and so the resulting symbol is just an approximation to the
        true p-stabilization (depending on how well alpha is approximated).

        INPUT:
            - ``p`` -- prime not dividing the level of self.
            - ``alpha`` -- eigenvalue for `U_p` 

        OUTPUT:

        A modular symbol with the same Hecke-eigenvalues as self away from `p` and eigenvalue `alpha` at `p`.

        EXAMPLES:
 
        ::

        sage: p = 3
        sage: M = 100
        sage: R = pAdicField(p,M)['y'] 
        sage: y = R.gen()
        sage: f = y**2-ap*y+p
        sage: v = f.roots()
        sage: alpha = v[0][0]
        sage: alpha
        2 + 3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 3^8 + 2*3^9 + 2*3^11 + 2*3^13 + 3^14 + 2*3^15 + 3^18 + 3^19 + 2*3^20 + 2*3^23 + 2*3^25 + 3^28 + 3^29 + 2*3^32 + 3^33 + 3^34 + 3^36 + 3^37 + 2*3^38 + 3^39 + 3^40 + 2*3^43 + 3^44 + 3^45 + 3^46 + 2*3^48 + 3^49 + 2*3^50 + 3^51 + 3^52 + 3^53 + 3^55 + 3^56 + 3^58 + 3^60 + 3^62 + 3^63 + 2*3^65 + 3^67 + 3^70 + 3^71 + 2*3^73 + 2*3^75 + 2*3^78 + 2*3^79 + 3^82 + 2*3^84 + 3^86 + 3^87 + 3^88 + 3^89 + 2*3^91 + 2*3^92 + 3^93 + 3^94 + 3^95 + O(3^100)
        sage: alpha^2 - E.ap(3)*alpha + 3
        O(3^100)
        sage: alpha=ZZ(alpha)
        sage: phi_alpha.hecke(2)-phi_alpha.scale(E.ap(2))
        [0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: (phi_alpha.hecke(3)-phi_alpha.scale(alpha)).valuation(3)
        101
        sage: phi_alpha.hecke(5)-phi_alpha.scale(E.ap(5))
        [0, 0, 0, 0, 0, 0, 0, 0, 0]

        """
        N = self.level()
        assert N % p != 0, "The level isn't prime to p"

        pp = Matrix(ZZ,[[p,0],[0,1]])

        manin = manin_relations(N*p)

        v = [] ## this list will contain the data of the p-stabilized symbol of level N*p

        ##  This loop runs through each generator at level Np and computes the value of the
        ##  p-stabilized symbol on each generator.  Here the following formula is being used:
        ##
        ##    (p-stabilize phi with U_p-eigenvalue alpha) = phi - 1/alpha phi | [p,0;0,1]
        ##  ---------------------------------------------------------------------------------
        for j in range(len(manin.generator_indices())):
            rj = manin.generator_indices(j)
            v = v+[self.eval(manin.coset_reps(rj))-self.eval(pp*manin.coset_reps(rj)).act_right(pp).scale(1/alpha)]
        return modsym_symk(v,manin)

    def p_stabilize_ordinary(self,p,ap,M):
        r"""
        Returns the unique `p`-ordinary `p`-stabilization of self

        INPUT:
            - ``p`` -- good ordinary prime
            - ``ap`` -- Hecke eigenvalue
            - ``M`` -- precision of `Q_p`

        OUTPUT:

        The unique modular symbol of level `N*p` with the same Hecke-eigenvalues as self away from `p` 
        and unit eigenvalue at `p`.

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]        
        sage: phi_ord = phi.p_stabilize_ordinary(3,E.ap(3),10); phi_ord
        [-47611/238060, 47615/95224, 95221/95224, 47611/238060, -95217/476120, 95225/95224, -5951/11903, -95229/95224, 95237/476120]
        sage: phi_ord.hecke(2) - phi_ord.scale(E.ap(2))
        [0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: phi_ord.hecke(5) - phi_ord.scale(E.ap(5))
        [0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: p = 3
        sage: M = 100
        sage: R = pAdicField(p,M)['y'] 
        sage: y = R.gen()
        sage: f = y**2-ap*y+p
        sage: v = f.roots()
        sage: alpha = v[0][0]
        sage: alpha
        2 + 3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 3^8 + 2*3^9 + 2*3^11 + 2*3^13 + 3^14 + 2*3^15 + 3^18 + 3^19 + 2*3^20 + 2*3^23 + 2*3^25 + 3^28 + 3^29 + 2*3^32 + 3^33 + 3^34 + 3^36 + 3^37 + 2*3^38 + 3^39 + 3^40 + 2*3^43 + 3^44 + 3^45 + 3^46 + 2*3^48 + 3^49 + 2*3^50 + 3^51 + 3^52 + 3^53 + 3^55 + 3^56 + 3^58 + 3^60 + 3^62 + 3^63 + 2*3^65 + 3^67 + 3^70 + 3^71 + 2*3^73 + 2*3^75 + 2*3^78 + 2*3^79 + 3^82 + 2*3^84 + 3^86 + 3^87 + 3^88 + 3^89 + 2*3^91 + 2*3^92 + 3^93 + 3^94 + 3^95 + O(3^100)
        sage: alpha^2 - E.ap(3)*alpha + 3
        O(3^100)
        sage: alpha=ZZ(alpha)
        sage: (phi_ord.hecke(3) - phi_ord.scale(alpha)).valuation(3)
        11

        """
        N = self.level()
        k = self.data(0).weight()
        assert N%p!=0, "The level isn't prime to p"
        assert (ap%p)!=0, "Not ordinary!"

        # makes alpha the unit root of Hecke poly
        R = pAdicField(p,M)['y'] 
        y = R.gen()
        f = y**2-ap*y+p**(k+1)
        v = f.roots()
        if ZZ(v[0][0])%p<>0:
            alpha = ZZ(v[0][0])
        else:
            alpha = ZZ(v[1][0])
        if self.full_data() == 0:
            self.compute_full_data_from_gen_data()
        return self.p_stabilize(p,alpha)

    def p_stabilize_critical(self,p,ap,M):
        r"""

        INPUT:
            - ``p`` -- good ordinary prime
            - ``ap`` -- Hecke eigenvalue
            - ``M`` -- precision of `Q_p`

        OUTPUT:

        The unique modular symbol of level `N*p` with the same Hecke-eigenvalues as self away from `p` 
        and non-unit eigenvalue at `p`.

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]        
        sage: phi_ss = phi.p_stabilize_critical(3,E.ap(3),10)
        sage: phi_ord.hecke(2) - phi_ord.scale(E.ap(2))
        [0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: phi_ord.hecke(5) - phi_ord.scale(E.ap(5))
        [0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: p = 3
        sage: M = 10
        sage: R = pAdicField(p,M)['y'] 
        sage: y = R.gen()
        sage: f = y**2-ap*y+p
        sage: v = f.roots()
        sage: beta = v[1][0]
        sage: beta^2 - E.ap(3)*beta + 3
        O(3^10)
        sage: beta=ZZ(beta)
        sage: (phi_ss.hecke(3) - phi_ss.scale(beta)).valuation(3)
        9

        """
        N = self.level()
        k = self.data(0).weight()
        assert N%p<>0, "The level isn't prime to p"
        assert (ap%p)!=0, "Not ordinary!"

        # makes alpha the non-unit root of Hecke poly
        R = pAdicField(p,M)['y'] 
        y = R.gen()
        f = y**2-ap*y+p**(k+1)
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

        INPUT:
            - ``q`` -- prime of the Hecke operator
            - ``p`` -- prime we are working modulo
            - ``M`` -- degree of accuracy of approximation

        OUTPUT:

        If there is some constant `c` such that `self|T_q - c * self` has valuation greater than 
        or equal to M, then it returns (True, c).  Otherwise, it returns False.

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi_ord = phi.p_stabilize_ordinary(3,E.ap(3),10)
        sage: phi_ord.is_Tq_eigen(2,3,10)
        (True, -2)
        sage: phi_ord.is_Tq_eigen(2,3,100)
        (True, -2)
        sage: phi_ord.is_Tq_eigen(2,3,1000)
        (True, -2)
        sage: phi_ord.is_Tq_eigen(3,3,10)
        (True, -95227/47611)
        sage: phi_ord.is_Tq_eigen(3,3,100)
        False
        
        """
        selfq = self.hecke(q)
        r = 0
        while (self.data(r).coef(0)==0) and (r<=self.ngens()):
            r = r+1
        if r > self.ngens():
        #all coefficients of Y^k are zero which I think forces it 
        #not to be eigen
            return False
        else:   
            c = selfq.data(r).coef(0)/self.data(r).coef(0)
            val = (selfq-self.scale(c)).valuation(p)
            if val >= M:
                return True,c
            else:
                return False
                
    def lift_to_OMS(self,p,M):
        r"""
        Returns a (`p`-adic) overconvergent modular symbol with `M` moments which lifts self up to an Eisenstein error

        INPUT:
            - ``p`` -- 
            - ``M`` --

        OUTPUT:

        EXAMPLES:
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
        # except the two vertical lines
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

        INPUT:
            - ``p`` --
            - ``M`` --

        OUTPUT:

        EXAMPLES:

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

        Applies psi to all polynomial coefficients of self and lifts them to `QQ`

        INPUT:
            - ``psi`` -- map from `K` to `Q_p` 
        
        OUTPUT:

        EXAMPLES:

        """
        v = [self.data[j].map(psi) for j in range(len(self.data))]
        return modsym_symk(self.level,v,self.manin)

    def coerce_to_Qp(self,p,M):
        r"""
        If `K` is the base_ring of self, this function takes all maps 
        `K-->Q_p` and applies them to self return a vector of 
        <modular symbol,map: `K-->Q_p`> as map varies over all such maps.  
        `M` is the accuracy

        INPUT:
            - ``p`` --
            - ``M`` --

        OUTPUT:

        EXAMPLES:
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

    def form_modsym_from_decomposition(A):
        r"""
        `A` is a piece of the result from a command like 
        ModularSymbols(---).decomposition()
        
        INPUT:
            - ``A`` -- 

        OUTPUT:

        EXAMPLES:

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

        INPUT:
            - ``p`` --
            - ``ap`` --

        OUTPUT:

        EXAMPLES:
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


def form_modsym_from_elliptic_curve(E):
    r"""
    Returns the modular symbol (in the sense of Stevens) associated to `E`

    INPUT:
        - ``E`` --


    OUTPUT:

    EXAMPLES:
    """
    L = E.padic_lseries(3)  #initializing the L-function at 3 to access mod sym data
    N = E.conductor()
    manin = manin_relations(N)
    v = []
    R = PolynomialRing(QQ,2,names='X,Y')
    X,Y = R.gens()
    for j in range(len(manin.generator_indices())):
        rj = manin.generator_indices(j)
        g = manin.coset_reps(rj)
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
    return modsym_symk(v,manin)

