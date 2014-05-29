from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.combinat.permutation import Permutation
from sage.misc.lazy_attribute import lazy_attribute
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

    def abs(self):
        """Returns a Permtuation object `w` such that `w(i)=|self(i)|`"""
        return Permutation([abs(i) for i in self.w])

    @lazy_attribute
    def _descents(self):
        """Calculates where the descents are"""
        descents = []
        if self.w[0] < 0:
            descents.append(0)
        for i in range(self.n):
            a = self.w[i-1]
            b = self.w[i]
            if a > 0 and b < 0:
                if a > b:
                    descents.append(i)
            elif a < 0 and b > 0:
                if a < b:
                    descents.append(i)
        if self.w[-1] < 0:
            descents.append(self.n)
        return descents

    def descents(self):
        """Returns a list of positions of descents"""
        return self._descents

    def number_of_descents(self):
        """Returns the number of descents of `self`"""
        return len(self._descents)

    @lazy_attribute
    def _neg_entries(self):
        return [i+1 for (i,j) in enumerate(self.w) if j < 0]
    
    def negative_entries(self):
        """Returns list of `i` such that `self(i) < 0`"""
        return self._neg_entries

    def number_of_negative_entries(self):
        """Returns the number of `i` such that `self(i) < 0`"""
        return len(self._neg_entries)
        
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
