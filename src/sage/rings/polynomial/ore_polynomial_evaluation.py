from sage.rings.morphism import RingHomomorphism
from sage.categories.homset import Hom


class OrePolynomialEvaluation(RingHomomorphism):
    def __init__(self, P, c):
        """
        Define the evaluation map for P and c.
        
        TESTS::

            sage: A.<t> = GF(5^3)[]
            sage: A = A.fraction_field()
            sage: der = A.derivation()
            sage: K = OrePolynomialRing(A, der, 'X')
            sage: P = K.random_element()
            sage: c = A.random_element()
            sage: eval = P.eval(c)
            sage: TestSuite(eval).run()
            
        """
        self.P = P
        self.c = c
        self.coeffs = P.list()
        parent = P.parent()
        self.domain = parent.base_ring()
        RingHomomorphism.__init__(self, Hom(self.domain, self.domain))
        if parent.twisting_morphism() is None:
            self._sigma = Hom(self.domain, self.domain).identity()
        else:
            self._sigma = parent.twisting_morphism()
        if parent.twisting_derivation() is None:
            self._derivation = lambda x: parent.base_ring()(0)
        else:
            self._derivation = parent.twisting_derivation()

    def _call_(self, a):
        coeffs = self.P.list()
        result = coeffs[0]
        for i in coeffs[1:]:
            a = self.c*self._sigma(a) + self._derivation(a)
            result += i * a
        return result

    def _repr_(self):
        if self.P.parent().twisting_derivation() is None:
            der = ""
        else:
            der = "der"
        if self.c == 0:
            sig = ""
        else:
            if self.P.parent().twisting_derivation() is not None:
                der += " + "
            sig = "{}*sig".format(self.c)
        evc = der + sig
        s = "Ore evaluation endomorphism of {} by the pseudo-linear morphism ".format(self.P) + evc + " with :"
        if self.P.parent().twisting_derivation() is not None:
            s += "\n der is {}".format(self._derivation)
        if self.c != 0:
            s += "\n sig is {}".format(self._sigma)
        s += "\n Defn:   a |----> "
        n = len(self.coeffs)
        if n == 1:
            s += "{}".format(self.coeffs[0])
            return s
        for i in reversed(self.coeffs):
            n -= 1
            if n > 1:
                if i == 1:
                    s += "(" + evc + ")^{}(a) + ".format(n)
                elif i != 0:
                    s += "{}*(".format(i) + evc + ")^{}(a) + ".format(n)
            elif n == 1:
                if i == 1:
                    s += "(" + evc + ")(a) + "
                elif i != 0:
                    s += "{}*(".format(i) + evc + ")(a) + "
            elif n == 0:
                if i != 0:
                    s += "{}+ ".format(i)
        s = s[:-2]
        return s

    #def _update_slots(self):
        
    def _hash_(self):
        return hash((self.P.parent(), self.P, self.c))

    def domain(self):
        return self.domain

    def category(self):
        return Hom(self.domain, self.domain).category()
    
    def _repr_short(self):
        if self.P.parent().twisting_derivation() is None:
            der = ""
        else:
            der = "der"
        if self.c == 0:
            sig = ""
        else:
            if self.P.parent().twisting_derivation() is not None:
                der += " + "
            sig = "{}*sig".format(self.c)
        evc = der + sig
        s = "a |----> "
        n = len(self.coeffs)
        if n == 1:
            s += "{}".format(self.coeffs[0])
            return s
        for i in reversed(self.coeffs):
            n -= 1
            if n > 1:
                if i == 1:
                    s += "(" + evc + ")^{}(a) + ".format(n)
                elif i != 0:
                    s += "{}*(".format(i) + evc + ")^{}(a) + ".format(n)
            elif n == 1:
                if i == 1:
                    s += "(" + evc + ")(a) + "
                elif i != 0:
                    s += "{}*(".format(i) + evc + ")(a) + "
            elif n == 0:
                if i != 0:
                    s += "{}+ ".format(i)
        s = s[:-2]
        return s

    #def _latex_(self):

    
