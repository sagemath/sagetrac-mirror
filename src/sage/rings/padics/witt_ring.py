from .witt_vector import WittVector
from sage.rings.ring import CommutativeRing
from sage.rings.integer_ring import ZZ
from sage.categories.commutative_rings import CommutativeRings
from sage.structure.unique_representation import UniqueRepresentation

class WittRing_base(CommutativeRing, UniqueRepresentation):
    Element = WittVector
    def __init__(self, base_ring, prec, prime, generate_op_polys=False, category=None):
        self.prec = prec
        self.prime = prime
        
        self._algorithm = 'none'
        self.sum_polys = None
        self.prod_polys = None
        
        if generate_op_polys:
            self._algorithm = 'standard'
            self.generate_sum_and_product_polynomials(base)
        
        if category is None:
            category = CommutativeRings()
        CommutativeRing.__init__(self, base_ring, category=category)
    
    def _repr_(self):
        return f"Ring of Witt Vectors of length {self.prec} over {self.base()}"
    
    def _algorithm(self):
        return self._algorithm
    
    def _coerce_map_from_(self, S):
        # Question: do we return True is S == self.base()?
        # We have the teichmuller lift, but I don't think that's
        # a "coercion" map, per se.
        return (S is ZZ)
    
    def _element_constructor_(self, x):
        if x in ZZ:
            return self.element_class(self, self._int_to_vector(x))
        elif isinstance(x, tuple) or isinstance(x, list):
            return self.element_class(self, x)
        else:
            return NotImplemented
    
    def _int_to_vector(self, k):
        p = self.prime
        
        should_negate = False
        if k < 0:
            k = -k
            should_negate = True
        
        vec_k = [k]
        for n in range(1, self.prec):
            total = k - k**(p**n) - sum(p**(n-i) * vec_k[n-i]**(p**i) for i in range(1, n))
            total //= p**n
            vec_k.append(total)
        
        if should_negate:
            if p == 2:
                return NotImplemented
            else:
                vec_k = [-x for x in vec_k]
        
        return vec_k
    
    def characteristic(self):
        if self.base(p).is_unit():
            # If p is invertible, W_n(R) is isomorphic to R^n.
            return self.base().characteristic()
        else:
            # This is a conjecture. It's known for char(R) == p.
            return p**(n-1) * self.base().characteristic()
    
    def precision(self):
        return self.prec
    
    def random_element(self):
        return self.element_class(self, tuple(self.base().random_element() for _ in range(self.prec)))
        
    def teichmuller_lift(self, x):
        if x not in self.base():
            raise Exception(f'{x} not in {self.base()}')
        else:
            return self.element_class(self, (x,) + tuple(0 for _ in range(self.prec-1)))

class WittRing_p_typical(WittRing_base):
    def __init__(self, base_ring, prec, prime, algorithm='none', category=None):
        if algorithm == 'none' or algorithm == 'finotti':
            generate_op_polys = False
        elif algorithm == 'standard':
            generate_op_polys = True
        else:
            raise ValueError(f'Algorithm {algorithm} is not a valid option.')
        
        self._algorithm = algorithm
        
        WittRing_base.__init__(self, base_ring, prec, prime, 
            generate_op_polys=generate_op_polys, category=category)

class WittRing_finite_field(WittRing_p_typical):
    def __init__(self, base_ring, prec, prime, category=None):
        self._algorithm = 'Zq_isomorphism'
        WittRing_p_typical.__init__(self, base_ring, prec, prime, 
            algorithm='none', category=category)

class WittRing_non_p_typical(WittRing_base):
    def __init__(self, base_ring, prec, prime, generate_op_polys=False, category=None):
        WittRing_base.__init__(self, base_ring, prec, prime, 
            generate_op_polys=generate_op_polys, category=category)

class WittRing_p_invertible(WittRing_non_p_typical):
    def __init__(self, base_ring, prec, prime, category=None):
        self._algorithm = 'standard_otf'
        WittRing_non_p_typical.__init__(self, base_ring, prec, prime,
            generate_op_polys=False, category=category)