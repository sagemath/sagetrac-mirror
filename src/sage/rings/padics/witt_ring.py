from .witt_vector import WittVector
from sage.rings.ring import CommutativeRing
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
        WittRing_p_typical.__init__(self, base_ring, prec, prime, 
            algorithm='none', category=category)

class WittRing_non_p_typical(WittRing_base):
    def __init__(self, base_ring, prec, prime, generate_op_polys=False, category=None):
        WittRing_base.__init__(self, base_ring, prec, prime, 
            generate_op_polys=generate_op_polys, category=category)

class WittRing_p_invertible(WittRing_non_p_typical):
    def __init__(self, base_ring, prec, prime, category=None):
        WittRing_non_p_typical.__init__(self, base_ring, prec, prime,
            generate_op_polys=False, category=category)