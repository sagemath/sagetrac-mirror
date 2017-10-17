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
    def evaluate(self, npoint):
        if isinstance(npoint,list):
            npoint=NumericalPoint(npoint,ring=self.ring)
        if not isinstance(npoint,NumericalPoint):
            raise TypeError("point provided must be of type NumericalPoint")
        return([f.subs(npoint.to_dict()) for f in self.polys])

class NumericalPoint():
    """
    A class for representing points numerically
    """
    def __init__(self, coords, ring=None, multiplicity=None, condition_number=None):
        """
        Construct from list of coordinates

        EXAMPLES
        R.<x,y,z> =PolynomialRing(CC,3)
        p = NumericalPoint([2,3,4],ring=R)
        p.coordinates
        p.to_dict()
        """
        self.coordinates = coords
        self.ring = ring
        self.multiplicity = multiplicity
        self.condition_number = condition_number
        # and so on as more args are added
    def to_dict(self):
        if self.ring != None:
            return(dict([(self.ring.gens()[i],self.coordinates[i]) for i in range(0,len(self.coordinates))]))
        else:
            raise AttributeError("please set a ring")


# class Homotopy(PolynomialSystem):
    """
    A doc string describing this class
    """
#    def __init__(self, polySys, params):


# STAND-


# RUNNING EXAMPLES
