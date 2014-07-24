from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.global_options import GlobalOptions
from sage.combinat.permutation import Permutation
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import Sets
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.matrix.all import matrix
import itertools

SignedPermutationOptions = GlobalOptions(name='signed permutations',
    doc=r"""
    Set the global options for elements of the permutation class. The
    defaults are for permutations to be displayed in list notation and
    the multiplication done from left to right (like in GAP) -- that
    is, `(\pi \psi)(i) = \psi(\pi(i))` for all `i`.

    .. NOTE::

        These options have no effect on signed permutation group elements.
    """,
    end_doc="""
    EXAMPLES::

        sage: p213 = Permutation([2,1,3])
        sage: p312 = Permutation([3,1,2])
        sage: Permutations.global_options(mult='l2r', display='list')
        sage: Permutations.global_options['display']
        'list'
        sage: p213
        [2, 1, 3]
        sage: Permutations.global_options(display='cycle')
        sage: p213
        (1,2)
        sage: Permutations.global_options(display='singleton')
        sage: p213
        (1,2)(3)
        sage: Permutations.global_options(display='list')

    ::

        sage: Permutations.global_options['mult']
        'l2r'
        sage: p213*p312
        [1, 3, 2]
        sage: Permutations.global_options(mult='r2l')
        sage: p213*p312
        [3, 2, 1]
        sage: Permutations.global_options.reset()
    """,
    display=dict(default="list",
                 description="Specifies how the permutations should be printed",
                 values=dict(list="the permutations are displayed in list notation"
                                  " (aka 1-line notation)",
                             cycle="the permutations are displayed in cycle notation"
                                  " (i. e., as products of disjoint cycles)",
                             singleton="the permutations are displayed in cycle notation"
                                       " with singleton cycles shown as well",
                             ),
                 alias=dict(word="reduced_word", reduced_expression="reduced_word"),
                 case_sensitive=False),
    latex=dict(default="list",
               description="Specifies how the permutations should be latexed",
               values=dict(list="latex as a list in one-line notation",
                           twoline="latex in two-line notation",
                           cycle="latex in cycle notation",
                           singleton="latex in cycle notation with singleton cycles shown as well",
                           reduced_word="latex as reduced words"),
               alias=dict(word="reduced_word", reduced_expression="reduced_word", oneline="list"),
               case_sensitive=False),
    latex_empty_str=dict(default="1",
                         description='The LaTeX representation of a reduced word when said word is empty',
                         checker=lambda char: isinstance(char,str)),
    generator_name=dict(default="s",
                        description="the letter used in latexing the reduced word",
                        checker=lambda char: isinstance(char,str)),
    mult=dict(default="l2r",
              description="The multiplication of permutations",
              values=dict(l2r="left to right: `(p_1 \cdot p_2)(x) = p_2(p_1(x))`",
                          r2l="right to left: `(p_1 \cdot p_2)(x) = p_1(p_2(x))`"),
              case_sensitive=False)
)

class SignedPermutation(Element):
    r"""
    A signed permutation.
    
    Converts ``l`` to a signed permutation on `\{\pm 1,\pm 2,\ldots,\pm n\}`.
    
    INPUT:
    
    - ``l`` -- a list of integers, viewed as one-line permutation notation.
    """
    __metaclass__ = ClasscallMetaclass
    @staticmethod
    def __classcall_private__(cls,l):
        """
        Using the metaclass to parse all the input.

        I'm not sure this actually *ought* to be done here.
        """
        if isinstance(l,SignedPermutation):
            return l
        if isinstance(l,str):
            #got given cycle notation, need to get one-line
            pass
        return type.__call__(cls,l)

    def __init__(self, l, natural_order = True):
        """Initializer
        
        INPUT:
        
        - ``l`` -- a list of integers
        """
        Element.__init__(self,SignedPermutations(len(l)))
        self.w = l
        self.n = len(l)
        self.natural_order = natural_order
    
    def __call__(self, i):
        #note that this is permutations acting on the left!
        if i > 0:
            return self.w[i-1]
        return -self.w[-i-1]
    
    def __str__(self):
        return repr(self)
    
    def _repr_(self):
        return repr(self.w)

    def __len__(self):
        return len(self.w)

    def __getitem__(self,key):
        return self.w[key]

    def __iter__(self):
        return iter(self.w)

    ######
    def _comp(self,a,b):
        """
        To avoid having to write everything twice for the two orders on \pm[n]
        used for signed permutations, this function compares a and b according
        to the appropriate order based on self.natural_order.
        """
        if self.natural_order:
            return a < b
        #now we're using the -1 < -2 < ... order 
        if a < 0 and b > 0:
            return True
        if a > 0 and b < 0:
            return False
        return abs(a) <= abs(b)
    
    def abs(self):
        """Returns a Permtuation object `w` such that `w(i)=|self(i)|`
        
        EXAMPLES::

            sage: w = SignedPermutation([1,-2,3])
            sage: w.abs()
            [1,2,3]
        """
        return Permutation([abs(i) for i in self.w])

    ####Multiplication####
    def __mul__(self, rp):
        if self.parent().global_options['mult'] == 'l2r':
            return self.right_action_product(rp)
        else:
            return self.left_action_product(rp)

    def __rmul__(self, lp):
        if self.parent().global_options['mult'] == 'l2r':
            return self.left_action_product(lp)
        else:
            return self.right_action_product(lp)

    def left_action_product(self,lp):
        """
        Apply lp first then self.
        """
        #make the two permutations the same length
        new_lp = lp.w[:] +[i+1 for i in range(lp.n,self.n)]
        new_p1 = self.w[:] + [i+1 for i in range(self.n,lp.n)]
        return SignedPermutation([new_p1[i-1] for i in\
            new_lp],self.natural_order)

    def right_action_product(self,rp):
        """
        Apply rp first then self.
        """
        new_rp = rp.w[:] + [i+1 for i in range(rp.n,self.n)]
        new_p1 = self.w[:] + [i+1 for i in range(self.n,rp,n)]
        return SignedPermutation([new_p1[i-1] for i in\
            new_rp],self.natural_order)

    ####Permutation Statistics####

    @lazy_attribute
    def _descents(self):
        """Calculates where the descents are"""
        descents = []
        if self.natural_order:
            if self.w[0] < 0:
                descents.append(0)
        else:
            if self.w[n] < 0:
                descents.append(n)
        descents.extend([i+1 for i in range(self.n-1) if\
            self._comp(self.w[i],self.w[i+1])])
        return descents
    
    def descent_set(self):
        r"""Returns a list of positions of descents.
        
        A signed permutation `\pi` has a descent at `i \in [0,n-1]` if `\pi(i)
        > \pi(i+1)`.
        """
        return self._descents
    
    def des(self):
        """Returns the number of descents of `self`"""
        return len(self._descents)
    
    def des_a(self):
        if self.natural_order and 0 in self.descents_:
            return len(self.descents_) - 1
        elif (not self.natural_order) and n in self.descents_:
            return len(self.descents_)-1
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
    
    def sign(self):
        if self.length() % 2 == 0:
            return 1
        else:
            return -1
    
    def fmaj(self):
        return 2*self.maj() + self.neg()
    
    def fdes(self):
        return self.des()+self.des_a()

    def reverse(self):
        """
        Returns the signed permutation obtained by reversing the one-line
        notation of `self`.
        """
        return SignedPermutation([w[i] for i in\
            range(self.n-1,-1,-1)],self.natural_order)

    def complement(self):
        """
        Returns the signed permutation which is the complement of `self`.
        """
        c = []
        for i in self.w:
            if i > 0:
                c.append(self.n+1-i)
            else:
                c.append(-(self.n+1)-i)
        return SignedPermutation(c,self.natural_order)

    def to_matrix(self):
        m = Matrix(self.n)
        for i in range(self.n):
            if self.w[i] < 0:
                m[i][abs(w[i])-1] = -1
            else:
                m[i][w[i]-1] = 1
        return m

    ####cycle representations####
    def to_cycles(self, singletons=True):
        cycles = []
        seen = dict([(i,False) for i in range(-self.n,self.n+1) if not i == 0])
        for i in range(self.n):
            if seen[self.w[i]]:
                continue
            cycleFirst = i+1
            cycle = [cycleFirst]
            seen[cycleFirst] = True
            next = self(cycleFirst)
            while next != cycleFirst:
                cycle.append(next)
                seen[next] = True
                next = self(next)
            if singletons or len(cycle) > 1:
                cycles.append(tuple(cycle))
        return cycles

    def cycle_string(self, singletons = False):
        """
        Returns a string of the permutation in cycle notation.

        If ``singletons = True``, it includes 1-cycles in the string.
        """
        cycles = self.to_cycles(singletons = singletons)
        if cycles == []:
            return "()"
        else:
            return "".join(["("+",".join([str(l) for l in x])+")" for x in cycles])
            
        
class SignedPermutations(UniqueRepresentation,Parent):
    r"""
    Class of signed permutations of `\{\pm 1, \pm 2,\ldots, \pm n\}`
    """
    Element = SignedPermutation
    global_options = SignedPermutationOptions
    def __init__(self, n, natural_order=True):
        self.n = n
        Parent.__init__(self, category = Sets())

    def _repr_(self):
        return "Signed Permutations on %s" % self.n

    def __iter__(self):
        for p in itertools.permutations(range(1,self.n+1)):
            for s in itertools.product([1,-1],repeat=self.n):
                yield(SignedPermutation([p[i]*s[i] for i in\
                    range(self.n)],self.natural_order))
    
    def _element_constructor(self,l):
        return self.element_class(l)
