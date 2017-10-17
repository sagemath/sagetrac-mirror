## TODO
# 1) Most classes and methods need examples, tests, I/O documentation

import warnings

## CLASSES

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
            self.polys = [(self.ring)(p)  for p in polys]
        else:
            raise TypeError("coefficient ring")
    def evaluate(self, npoint):
        if isinstance(npoint,list):
            npoint=NumericalPoint(npoint,ring=self.ring)
        if not isinstance(npoint,NumericalPoint):
            raise TypeError("point provided must be of type NumericalPoint")
        return([f.subs(npoint.to_dict()) for f in self.polys])
    def __str__(self):
        return("%s over %s. " %(self.polys, self.ring))
        

class NumericalPoint(SageObject):
    """
    A class for representing points numerically
    """
    def __init__(self, coords, ring=None, multiplicity=None, condition_number=None):
        """
        Construct from list of coordinates

        EXAMPLES::

            sage: R.<x,y,z> =PolynomialRing(CC,3)
            sage: p = NumericalPoint([2,3,4],ring=R)
            sage: p.coordinates
            [2, 3, 4]
            sage: p.to_dict()
            {z: 4, y: 3, x: 2}        
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
    def __str__(self):
        return("A numerical point in CC^%s." %(len(self.coordinates)))

        
class WitnessSet(SageObject):
    def __init__(self, polySys, forms, points):
        """
        INPUT: 
            *) polySys, an object of type PolynomialSystem
            *) forms, an object of type PolynomialSystem consisting of linear forms w/ the same ring as polySys
            *) points --- a list of objects of type NumericalPoint 
        """
        if not isinstance(polySys, PolynomialSystem):
            raise TypeError("first argument should be a PolynomialSystem")
        if not isinstance(forms, PolynomialSystem):
            raise TypeError("second argument should be a polynomial system")
        if not (isinstance(points, list) and len(set([p.parent() for p in points]))==1 and points[0].parent() == NumericalPoint):
            raise TypeError("third argument should be a list of NumericalPoints")
        if not polySys.ring == forms.ring:
            raise TypeError("make sure first two arguments share a common ring")
        self.system = polySys
        self.linear_forms = forms
        self.witness_points = points
        self.dimension = len(forms.polys)
    def __str__(self):
        return("A witness set for a dimension-%s component with %s points." %(self.dimension, len(self.witness_points)))

        
class NumericalIrreducibleDecomposition(SageObject):
    """
    a class which organizes the witness sets appearing in the NID of a variety
    """
    def __init__(self):
        self.components =dict()
    def append_witness_set(self, wset):
        if not isinstance(wset,WitnessSet):
            raise TypeError("must append with a witness set")
        if wset.dimension in self.components.keys():
            self.components[wset.dimension].append(wset)
        else:
            self.components[wset.dimension] = [wset]
    def __str__(self):
        return_string = ""
        for i in self.components.keys():
            return_string += "Dimension "+str(i) + ":" + "\n"
            for j in self.components[i]:
                return_string += "    Component of degree "+str(len(j.witness_points)) + "\n "
        return(return_string)

# class Homotopy(PolynomialSystem):
    """
    A doc string describing this class
    """
#    def __init__(self, polySys, params):


## STANDALONE FUNCTION



