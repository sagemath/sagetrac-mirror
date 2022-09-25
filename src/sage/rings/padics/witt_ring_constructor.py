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
        raise TypeError(f'Cannot create Ring of Witt Vectors over {base}: {base} is not a commutative ring.')
    
    char = base_ring.characteristic()
    prime = None
    if p is None:
        if char not in _Primes:
            raise ValueError(f'Cannot create Ring of Witt Vectors over {base}: {base} has non-prime characteristic, and no prime was supplied.')
        else:
            prime = char
    else:
        prime = p
    
    if prime == char: # p-typical
        pass
    else: # non-p-typical
        pass
    
    return WittRing_general()