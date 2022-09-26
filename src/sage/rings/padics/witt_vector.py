from sage.structure.element import CommutativeRingElement
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

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
    
    def _neg_(self):
        P = self.parent()
        C = self.__class__
        # If p == 2, -1 == (1, 1, 1, ...)
        # Otherwise, -1 == (-1, 0, 0, ...)
        if P.prime == 2:
            all_ones = P(tuple(1 for _ in range(self.prec)))
            return all_ones*self
        neg_vec = tuple(-self.vec[i] for i in range(self.prec))
        return C(P, vec=neg_vec)
    
    def _div_(self, other):
        P = self.parent()
        # As a slight optimization, we'll check for one ahead of time.
        # This has the benefit of allowing us to create polynomials,
        # even if ``P._algorithm`` is 'none'.
        if other == P.one():
            return self
        elif self == P.one():
            return other._invert_()
        
        return self * other._invert_()
    
    def _invert_(self):
        if not self.vec[0].is_unit():
            raise ZeroDivisionError(f"Inverse of {self} does not exist.")
        P = self.parent()
        C = self.__class__
        
        if self == P.one():
            return self
        if self.prec == 1:
            return P((self.vec[0]**-1, ))
        
        # Strategy: Multiply ``self`` by ``(Y_0, Y_1, ...)``, set equal
        # to (1, 0, 0, ...), and solve.
        var_names = [f'Y{i}' for i in range(1, self.prec)]
        poly_ring = PolynomialRing(P.base(), var_names)
        inv_vec = list( (self.vec[0]**-1, ) + poly_ring.gens()) # We'll fill this in one-by-one
        
        # TODO: Remove the algorithm argument once other algs are implemented
        from sage.rings.padics.witt_ring_constructor import WittRing
        W = WittRing(poly_ring, prec=P.prec, algorithm=P._algorithm)
        prod_vec = (W(self.vec)*W(inv_vec)).vec
        for i in range(1, self.prec):
            poly = prod_vec[i](inv_vec[1:])
            Y_i = poly.parent().gens()[i-1]
            inv_vec[i] = -poly.constant_coefficient() / poly.monomial_coefficient(Y_i)
        
        return C(P, vec=inv_vec)