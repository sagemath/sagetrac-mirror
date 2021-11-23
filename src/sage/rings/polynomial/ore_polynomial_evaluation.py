import sage
from sage.categories.morphism import IdentityMorphism, Morphism
from sage.categories.homset import Hom


class OrePolynomialEvaluation(Morphism):
    def __init__(self, P, c):
        self.P = P
        self.c = c
        parent = P.parent()
        super().__init__(Hom(parent.base_ring(), parent.base_ring()))
        if parent.twisting_morphism() is None:
            self._sigma = Hom(parent.base_ring(), parent.base_ring()).identity()
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
        s += "\n    a |----> "
        coeffs = self.P.list()
        n = len(coeffs)
        if n == 1:
            s += "{}".format(coeffs[0])
            return s
        for i in reversed(coeffs):
            n -= 1
            if i == 1:
                s += "(" + evc + ")^{}(a) + ".format(n)
            elif i != 0:
                s += "{}*(".format(i) + evc + ")^{}(a) + ".format(n)
            s = s[:-2]
        return s
