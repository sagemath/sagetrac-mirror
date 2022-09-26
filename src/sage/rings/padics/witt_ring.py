from .witt_vector import WittVector
from sage.rings.ring import CommutativeRing
from sage.rings.integer_ring import ZZ
from sage.categories.commutative_rings import CommutativeRings
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class WittRing_base(CommutativeRing, UniqueRepresentation):
    Element = WittVector
    def __init__(self, base_ring, prec, prime, algorithm='none', category=None):
        self.prec = prec
        self.prime = prime
        
        self._algorithm = algorithm
        self.sum_polynomials = None
        self.prod_polynomials = None
        
        if algorithm == 'standard':
            self._generate_sum_and_product_polynomials(base_ring)
        
        if category is None:
            category = CommutativeRings()
        CommutativeRing.__init__(self, base_ring, category=category)
    
    def _repr_(self):
        return f"Ring of Witt Vectors of length {self.prec} over {self.base()}"
    
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
    
    def _generate_sum_and_product_polynomials(self, base):
        p = self.prime
        prec = self.prec
        x_var_names = ['X{}'.format(i) for i in range(prec)]
        y_var_names = ['Y{}'.format(i) for i in range(prec)]
        var_names = x_var_names + y_var_names
        
        # Okay, what's going on here? Sage, by default, relies on 
        # Singular for Multivariate Polynomial Rings, but Singular uses
        # only SIXTEEN bits (unsigned) to store its exponents. So if we 
        # want exponents larger than 2^16 - 1, we have to use the 
        # generic implementation. However, after some experimentation,
        # it seems like the generic implementation is faster?
        #
        # After trying to compute S_4 for p=5, it looks like generic is
        # faster for  very small polys, and MUCH slower for large polys.
        # So we'll default to singular unless we can't use it.
        # 
        # Remark: Since when is SIXTEEN bits sufficient for anyone???
        #
        if p**(prec-1) >= 2**16:
            implementation = 'generic'
        else:
            implementation = 'singular'
            
        # We first generate the "universal" polynomials and then project
        # to the base ring.
        R = PolynomialRing(ZZ, var_names, implementation=implementation)
        x_y_vars = R.gens()
        x_vars = x_y_vars[:prec]
        y_vars = x_y_vars[prec:]
        
        self.sum_polynomials = [0]*(self.prec)
        for n in range(prec):
            s_n = x_vars[n] + y_vars[n]
            for i in range(n):
                s_n += (x_vars[i]**(p**(n-i)) + y_vars[i]**(p**(n-i)) - self.sum_polynomials[i]**(p**(n-i))) // p**(n-i)
            self.sum_polynomials[n] = s_n
        
        self.prod_polynomials = [x_vars[0] * y_vars[0]] + [0]*(self.prec)
        for n in range(1, prec):
            x_poly = sum([p**i * x_vars[i]**(p**(n-i)) for i in range(n+1)])
            y_poly = sum([p**i * y_vars[i]**(p**(n-i)) for i in range(n+1)])
            p_poly = sum([p**i * self.prod_polynomials[i]**(p**(n-i)) for i in range(n)])
            p_n = (x_poly*y_poly - p_poly) // p**n
            self.prod_polynomials[n] = p_n
        
        # We have to use generic here, because Singular doesn't support
        # Polynomial Rings over Polynomial Rings. For example, 
        # ``PolynomialRing(GF(5)['x'], ['X', 'Y'], implementation='singular')``
        # will fail.
        S = PolynomialRing(base, x_y_vars, implementation='generic')
        for n in range(prec):
            self.sum_polynomials[n] = S(self.sum_polynomials[n])
            self.prod_polynomials[n] = S(self.prod_polynomials[n])
    
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
    def __init__(self, base_ring, prec, prime, algorithm=None, category=None):
        WittRing_base.__init__(self, base_ring, prec, prime, 
            algorithm=algorithm, category=category)

class WittRing_finite_field(WittRing_p_typical):
    def __init__(self, base_ring, prec, prime, category=None):
        WittRing_p_typical.__init__(self, base_ring, prec, prime, 
            algorithm='Zq_isomorphism', category=category)

class WittRing_non_p_typical(WittRing_base):
    def __init__(self, base_ring, prec, prime, algorithm=None, category=None):
        WittRing_base.__init__(self, base_ring, prec, prime, 
            algorithm=algorithm, category=category)

class WittRing_p_invertible(WittRing_non_p_typical):
    def __init__(self, base_ring, prec, prime, category=None):
        WittRing_non_p_typical.__init__(self, base_ring, prec, prime,
            algorithm='standard_otf', category=category)