from sage.structure.element import Element
from sage.structure.parent import Parent
import itertools

class SignedPermutation(Element):
    r"""
    A signed permutation.
    
    Converts ``l`` to a signed permutation on `\{\pm 1,\pm 2,\ldots,\pm n\}`.
    
    INPUT:
    
    - ``l`` -- a list of integers, viewed as one-line permutation notation.
    """
    
    def __init__(self, l):
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
        
    def __str__(self):
        return repr(self)
        
    def _repr_(self):
        return repr(self.w)
        
class SignedPermutations(Parent):
    r"""
    Class of signed permutations of `\{\pm 1, \pm 2,\ldots, \pm n\}`
    """
    
    Element = SignedPermutation
    
    def __init__(self, n):
        self.n = n
    
    def _repr_(self):
        return "Signed Permutations on %s",self.n
        
    def __iter__(self):
        for p in itertools.permutations(range(1,self.n+1)):
            for s in itertools.product([1,-1],repeat=self.n):
                yield(SignedPermutation([p[i]*s[i] for i in range(self.n)]))