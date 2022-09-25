from sage.structure.element import CommutativeRingElement

class WittVector(CommutativeRingElement):
    def __init__(self, parent, vec=None):
        self.prec = parent.precision()
        B = parent.base()
        if vec is not None:
            if len(vec) != self.prec:
                raise ValueError(f'{vec} is not the correct length. Expected length to be {self.prec}.')
            self.vec = tuple(B(x) for x in vec)
        else:
            self.vec = (B(0) for i in range(self.prec))
        CommutativeRingElement.__init__(self, parent)
    
    def _richcmp_(self, other, op):
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            return self.vec == other.vec
        elif op == op_NE:
            return self.vec != other.vec
        else:
            return NotImplemented
    
    def _repr_(self):
        return '(' + ', '.join(map(str, self.vec)) + ')'
    
    def _add_(self, other):
        P = self.parent()
        C = self.__class__
        
        # As a slight optimization, we'll check for zero ahead of time.
        # This has the benefit of allowing us to create polynomials, 
        # even if ``P._algorithm`` is 'none'.
        if other == P.zero():
            return self
        elif self == P.zero():
            return other
        
        alg = P._algorithm
        if alg == 'standard':
            s = P.sum_polynomials
            # note here this is tuple addition, i.e. concatenation
            sum_vec = tuple(s[i](*(self.vec + other.vec)) for i in range(self.prec))
            return C(P, vec=sum_vec)
        else:
            return NotImplemented
    
    def _mul_(self, other):
        P = self.parent()
        C = self.__class__
        
        # As a slight optimization, we'll check for zero or one ahead of time.
        # This has the benefit of allowing us to create polynomials,
        # even if ``P._algorithm`` is 'none'.
        if self == P.zero() or other == P.zero():
            return P.zero()
        elif other == P.one():
            return self
        elif self == P.one():
            return other
        
        alg = P._algorithm
        if alg == 'standard':
            p = P.prod_polynomials
            # note here this is tuple addition, i.e. concatenation
            prod_vec = tuple(p[i](*(self.vec + other.vec)) for i in range(self.prec))
            return C(P, vec=prod_vec)
        else:
            return NotImplemented