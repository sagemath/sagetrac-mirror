from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from itertools import permutations, product

class SignedPermutation(Element):
    r"""
    A signed permutation.
    
    Converts ``l`` to a signed permutation on `\{\pm 1,\pm 2,\ldots,\pm n\}`.
    
    INPUT:
    
    - ``l`` -- a list of integers, viewed as one-line permutation notation.
    """
    
    __metaclass__ = ClasscallMetaclass
    
    def __init__(self, parent, l):
        """Constructor 
        
        INPUT:
        
        - ``l`` -- a list of integers
        """
        self.w = l
        self.n = len(l)
        
    def __call__(self, i):
        if i > 0:
            return self.w[i-1]
        return -self.w[-i-1]
        
class SignedPermutations(UniqueRepresentation, Parent):
    r"""
    Class of signed permutations of `\{\pm 1, \pm 2,\ldots, \pm n\}`
    """
    
    Element = SignedPermutation
    
    def __init__(self, n):
        self.n = n
        
    def __iter__(self):
        for p in itertools.permutations(range(1,self.n+1)):
            for s in itertools.product([1,-1],repeat=self.n):
                yield(SignedPermutation([p[i]*s[i] for i in range(self.n)]))