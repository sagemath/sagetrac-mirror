class Polynomial_System(SageObject):
    """
    A doc string describing this class
    """
    def __init__(self, polys):
        """
        this is a constructor that takes a list of polynomials and 
        returns an object of class polynomial syste
        """
        if not isinstance(polys,list):
             raise TypeError("incorrect input")
        L=list(set([p.parent() for p in polys]))
        if len(L) != 1 or not isinstance(L[0],Ring):
            raise TypeError("polynomials don't have same parent ring")
        self.polys = polys
        self.ring = polys[0].parent() # not strictly necessary
        
class Homotopy(Polynomial_System):
    """
    A doc string describing this class
    """
    def __init__(self, polySys, params):

