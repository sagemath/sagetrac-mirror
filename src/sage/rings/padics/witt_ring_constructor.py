import sage.rings.ring as ring

from sage.categories.commutative_rings import CommutativeRings
_CommutativeRings = CommutativeRings()
from sage.sets.primes import Primes
_Primes = Primes()

from sage.rings.padics.witt_ring import *

def WittRing(base_ring, prec=1, p=None, *args, **kwds):
    
    if not ring.is_Ring(base_ring):
        raise TypeError(f'Base ring {base_ring} must be a ring.')
    
    if base_ring not in _CommutativeRings:
        raise TypeError(f'Cannot create Ring of Witt Vectors over {base_ring}: {base_ring} is not a commutative ring.')
    
    char = base_ring.characteristic()
    prime = None
    if p is None:
        if char not in _Primes:
            raise ValueError(f'Cannot create Ring of Witt Vectors over {base_ring}: {base_ring} has non-prime characteristic and no prime was supplied.')
        else:
            prime = char
    else:
        prime = p
    
    # Look for 'algorithm' and 'category' in kwds?
    
    # Throw error if given a keyword that's not used.
    
    if prime == char: # p-typical
        if base_ring.is_field() and base_ring.is_finite():
            return WittRing_finite_field(base_ring, prec, prime, category=_CommutativeRings)
        else:
            return WittRing_p_typical(base_ring, prec, prime, category=_CommutativeRings)
    else: # non-p-typical
        if base_ring(prime).is_unit():
            return WittRing_p_invertible(base_ring, prec, prime, category=_CommutativeRings)
        else:
            return WittRing_non_p_typical(base_ring, prec, prime, category=_CommutativeRings)