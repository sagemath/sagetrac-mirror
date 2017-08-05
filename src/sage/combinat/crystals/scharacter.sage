"""
Characters of typical representations of gl(m|n). The characters
are computed for irreducible representations of atypicality <= 1.
The characters of the Kac modules are computed in every case.

The weight lam = (lam1, ... ,lam_m; lam_{m+1}, ..., lam_{m+n}) is
dominant if lam1 >= ... >= lam_m an lam_{m+1} >= ... >= lam_{m+n}.

First an example of a typical weight. In this case, the
Kac module and the irreducible coincide, and we can also
compare the degree of the representation to the cardinality
of the crystal. Note that the "sharp" operation turns the
shape [5,2,1,1] into the weight [5,2,2,0]

sage: G0=SuperWeylCharacterRing(2,2)
sage: G0.atypicality((5,2,2,0))
0
sage: G0.character((5,2,2,0))
A1xA1(3,1,3,2) + A1xA1(3,1,4,1) + A1xA1(3,2,3,1) + A1xA1(3,0,4,2) + A1xA1(4,1,2,2) + 2*A1xA1(4,1,3,1) + A1xA1(4,1,4,0) + A1xA1(4,2,2,1) + A1xA1(4,2,3,0) + A1xA1(4,0,3,2) + A1xA1(4,0,4,1) + A1xA1(5,1,2,1) + A1xA1(5,1,3,0) + A1xA1(5,2,2,0) + A1xA1(5,0,3,1)
sage: G0.character((5,2,2,0))
A1xA1(3,1,3,2) + A1xA1(3,1,4,1) + A1xA1(3,2,3,1) + A1xA1(3,0,4,2) + A1xA1(4,1,2,2) + 2*A1xA1(4,1,3,1) + A1xA1(4,1,4,0) + A1xA1(4,2,2,1) + A1xA1(4,2,3,0) + A1xA1(4,0,3,2) + A1xA1(4,0,4,1) + A1xA1(5,1,2,1) + A1xA1(5,1,3,0) + A1xA1(5,2,2,0) + A1xA1(5,0,3,1)
sage: G0.character((5,2,2,0)).degree()
192
sage: crystals.Tableaux(['A',[1,1]],shape=[5,2,1,1]).cardinality()
192

sage: G0.atypicality((5,1,1,0))
1
sage: G0.character((5,1,1,0))
A1xA1(3,1,2,1) + A1xA1(3,0,3,1) + A1xA1(4,1,1,1) + A1xA1(4,1,2,0) + A1xA1(4,0,2,1) + A1xA1(4,0,3,0) + A1xA1(5,1,1,0) + A1xA1(5,0,2,0)
sage: G0.character((5,1,1,0)).degree()
92
sage: crystals.Tableaux(['A',[1,1]],shape=[5,1,1]).cardinality()
92
"""

class SuperWeylCharacterRing(WeylCharacterRing):
    """
    Character Ring for gl(m|n)
    """
    @staticmethod
    def __classcall__(cls, m, n):
        return super(WeylCharacterRing, cls).__classcall__(cls, m, n)
    
    def __init__(self, m, n):
        WeylCharacterRing.__init__(self, [['A',m-1],['A',n-1]])
        self.m = m
        self.n = n
        self._space = self.space()
        self.delt = {i : self._space.basis()[i-1] for i in [1..m]}
        self.epsi = {i : self._space.basis()[m+i-1] for i in [1..n]}
        self.even_positive_roots = self._space.positive_roots()
        self.odd_positive_roots = [self._space(self.delt[i]-self.epsi[j]) for i in [1..m] for j in [1..n]]
        self.rho0 = sum(self.even_positive_roots)/2
        self.rho1 = sum(self.odd_positive_roots)/2
        self.rho = self.rho0-self.rho1
        self.g0 = WeightRing(self)
        g0 = self.g0
        self.l1 = prod(g0.one()+g0(-a) for a in self.odd_positive_roots)
        self.L1 = self(self.l1)
    
    def KacCharacter(self, lam):
        return self(lam)*self.L1

    def inner_product(self, lam, mu):
        [m,n] = [self.m, self.n]
        return sum(lam[i]*mu[i] for i in [0..m-1])-sum(lam[i]*mu[i] for i in [m,m+n-1])

    def atypical_roots(self, lam):
        lam = self._space(lam)
        return [a for a in self.odd_positive_roots if self.inner_product(a,lam+self.rho)==0]

    def atypicality(self, lam):
        """
        The number of independent atypical roots.
        """
        return Matrix([x.to_vector() for x in self.atypical_roots(lam)]).rank()

    def weyl_group(self):
        return self._space.weyl_group(prefix="s")

    def weyl_action(self, expr, w):
        """
        expr: a Weight ring element
        w: a Weyl group element
        Returns w(expr)
        """
        d = expr.monomial_coefficients()
        return sum(d[k]*k.weyl_action(w) for k in d)

    def weight_to_character(self, x):
        """
        x: Weight ring element.
        returns sum (-1)^w w(x)/weyl_denominator
        as an element of self.
        """
        d = x.monomial_coefficients()
        ret = self.zero()
        for k in d:
            [epsilon,v]=self._reduce(k)
            if epsilon != 0:
                ret += d[k]*epsilon*(self(v-self.rho0))
        return ret
        
    def _reduce(self, x):
        alphacheck = self._space.simple_coroots()
        alpha = self._space.simple_roots()
        sr = self._space.weyl_group().simple_reflections()
        [epsilon, ret] = [1,x]
        done = False
        while not done:
            done = True
            for i in self._space.index_set():
                c = ret.inner_product(alphacheck[i])
                if c == 0:
                    return [0, self._space.zero()]
                elif c < 0:
                    epsilon = -epsilon
                    ret -= c*alpha[i]
                    done = False
                    break
        return [epsilon, ret]

    def character(self, lam):
        """
        This gives the character of the irreducible module with highest weight
        lam. If lam is typical, this is the Kac crystal. If lam has atypicality
        index 1, then the formula for the character is from Van der Jeugt, J.;
        Hughes, J. W. B.; King, R. C.; Thierry-Mieg, J, A character formula for
        singly atypical modules of the Lie superalgebra sl(m/n).

        In another paper, J. van der Jeugt, J.W.B. Hughes, R.C. King,
        J. Thierry-Mieg, Character formulas for irreducible modules of the Lie
        superalgebras sl(m|n), J. Math. Phys. 31 (1990) 2278–2304 conjectured a
        more complicated general formula for the character with arbitrary
        atypicality. This formula was eventually proved by Su and Zhang,
        Character and dimension formulae for general linear superalgebra,
        Adv. Math. 211 (2007), no. 1, 1–33, making use of important results of
        Serganova and Brundan. We have not yet implemented this formula.
        """
        g0 = WeightRing(self)
        lam = self._space(lam)
        if self.atypicality(lam) == 0:
            return self.KacCharacter(lam)
        elif self.atypicality(lam) == 1:
            W = self.weyl_group()
            beta = self.atypical_roots(lam)[0]
            ex = g0(lam+self.rho0)*prod(g0.one()+g0(-a) for a in self.odd_positive_roots if a is not beta)
            return self.weight_to_character(ex)
        else:
            print "not implemented for atypicality >1."
            return None
            
    def irreducible_socle(self, lam):
        """
        If lam is a dominant weight with atypicality 1, we observe
        that the kernel of the canonical map K(lam) -> L(lam)
        from the Kac module to its irreducible quotient is
        often irreducible. I this is so, this method
        attempts to find mu such that we have a short
        exact sequence L(mu) -> K(lam) -> L(lam).

        sage: G0=SuperWeylCharacterRing(2,2)
        sage: G0.atypicality([2,1,1,0])
        1
        sage: G0.irreducible_socle([2,1,1,0])
        (2, 0, 1, 1)
        sage: G0.irreducible_socle([2,0,1,1])
        (2, -2, 2, 2)
        """
        f = self.KacCharacter(lam)-self.character(lam)
        for mu in f.monomial_coefficients().keys():
            if f == self.character(mu):
                return mu
        print "irreducible socle not found."
        return None



