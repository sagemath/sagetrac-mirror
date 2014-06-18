from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.combinat.permutation import Permutation
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import Sets
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
        Element.__init__(self,SignedPermutations(len(l)))
        self.w = l
        self.n = len(l)
        #self._parent = SignedPermutations
    
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
        descents.extend([i+1 for i in range(self.n) if self.w[i] > self.w[i+1]])
        return descents
    
    def descents(self):
        """Returns a list of positions of descents"""
        return self._descents
    
    def number_of_descents(self):
        """Returns the number of descents of `self`"""
        return len(self._descents)
    
    def des_a(self):
        if 0 in self.descents_:
            return len(self.descents_) - 1
        else:
            return len(self.descents_)
    
    @lazy_attribute
    def _neg_entries(self):
        return [i+1 for (i,j) in enumerate(self.w) if j < 0]
    
    def negative_entries(self):
        """Returns list of `i` such that `self(i) < 0`"""
        return self._neg_entries
    
    def neg(self):
        """Returns the number of `i` such that `self(i) < 0`"""
        return len(self._neg_entries)
    
    def maj(self):
        """Returns the major index of `self`"""
        return sum(self._descents)
    
    def nsum(self):
        """Returns negative number sum of `self`"""
        return -sum(self._neg_entries)
    
    def length(self):
        return self.inv() + self.nsum()
    
    def inversions(self):
        return self._inversions
    
    @lazy_attribute
    def _inversions(self):
        return [(i,j) for i in range(1,n+1) for j in range(i+1,n+1) if\
                self.w[i] > self.w[j]]
    
    def inv(self):
        return len(self._inversions)
    
    def length(self):
        if self.length() % 2 == 0:
            return 1
        else:
            return -1
    
    def fmaj(self):
        #note that flag-major and fmaj are different!
        #this is truly fmaj from Adin, Brenti, Roichman p. 218
        return 2*self.maj() + self.neg()
    
    def fdes(self):
        return 2*self.des_a() + (self.w[0] < 0)
        
class SignedPermutations(UniqueRepresentation,Parent):
    r"""
    Class of signed permutations of `\{\pm 1, \pm 2,\ldots, \pm n\}`
    """
    Element = SignedPermutation
    def __init__(self, n):
        self.n = n
        Parent.__init__(self, category = Sets())

    def _repr_(self):
        return "Signed Permutations on %s" % self.n

    def __iter__(self):
        for p in itertools.permutations(range(1,self.n+1)):
            for s in itertools.product([1,-1],repeat=self.n):
                yield(SignedPermutation([p[i]*s[i] for i in range(self.n)]))
    
    def _element_constructor(self,l):
        return self.element_class(l)