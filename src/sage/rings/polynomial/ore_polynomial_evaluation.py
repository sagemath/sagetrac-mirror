import sage
from sage.misc.prandom import randint
from sage.misc.cachefunc import cached_method
from sage.rings.infinity import Infinity
from sage.structure.category_object import normalize_names

from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import Algebra
from sage.rings.integer import Integer
from sage.structure.element import Element

from sage.categories.commutative_rings import CommutativeRings
from sage.categories.algebras import Algebras
from sage.rings.ring import _Fields


class OrePolynomialEvaluation():
    def __init__(self, P, c):
        self.P = P
        self.c = c
        self._sigma = P.parent().twisting_morphism()
        self._derivation = P.parent().twisting_derivation()

    def __call__(self, a):
        coeffs = self.P.list()
        result = coeffs[0]
        #On peut factoriser facilement cette partie en initialisant sigma à l'identité et derivation à 0. Mais je n'ai pas trouvé comment faire.
        if self._sigma == None and self._derivation == None:
            for i in coeffs[1:]:
                a = self.c*a
                result += i * a
        if self._sigma == None and self._derivation != None:
            for i in coeffs[1:]:
                a = self.c*a + self._derivation(a)
                result += i * a
        if self._sigma != None and self._derivation == None:
            for i in coeffs[1:]:
                a = self.c*self._sigma(a)
                result += i * a
        if self._sigma != None and self._derivation != None:
            for i in coeffs[1:]:
                a = self.c*self._sigma(a) + self._derivation(a)
                result += i * a
        return result

    def __repr__(self):
        s = "Ore Polynomial Evaluation P(c*sigma + der) with :\n"
        s += "P is %s \n" % self.P._repr_()
        s += "c is %s\n" % str(self.c)
        if self._sigma == None:
            s += "sigma is the identity morphism on %s \n" % self.P.parent().base_ring()._repr_()
        else:
            s += "sigma is %s \n" % self._sigma._repr_()
        if self._derivation == None:
            s += "der is 0"
        else:
            s += "der is %s" % self._derivation._repr_()
        return s
