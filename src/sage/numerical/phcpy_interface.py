import warnings

class PolynomialSystem(SageObject):
    """
    A class for systems of polynomials, primarily for numerical purposes.
    """
    def __init__(self, polys):
        """
        This is a constructor that takes a list of polynomials and
        returns an object of class PolynomialSystem.
        """
        if not isinstance(polys,list):
             raise TypeError("incorrect input")
        L=list(set([p.parent() for p in polys]))
        if len(L) != 1 or not isinstance(L[0],Ring):
            raise TypeError("polynomials don't have same parent ring")
        # better error handling for coefficient field NEEDED
        if L[0].base_ring() in set([RR,CC,QQ,ZZ]):
            self.polys = polys
            self.ring = polys[0].parent() # not strictly necessary
        elif isinstance(L[0].base_ring(),sage.symbolic.ring.SymbolicRing):
            warnings.warn("SymbolicRing expressions not checked for consistency.",RuntimeWarning)
            myvars=list(set(flatten([list(p.variables()) for p in polys])))
            self.ring = PolynomialRing(CC,len(myvars),myvars)
            print("hi")
            self.polys = [(self.ring)(p)  for p in polys]
        else:
            raise TypeError("coefficient ring")



# class Homotopy(PolynomialSystem):
    """
    A doc string describing this class
    """
#    def __init__(self, polySys, params):
