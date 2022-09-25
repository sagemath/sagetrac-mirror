from .witt_vector import WittVector
from sage.rings.ring import CommutativeRing
from sage.structure.unique_representation import UniqueRepresentation

class WittRing_general(CommutativeRing, UniqueRepresentation):
    Element = WittVector
    def __init__(self):
        pass

class WittRing_p_typical(WittRing_general):
    pass

class WittRing_finite_field(WittRing_p_typical):
    pass

class WittRing_non_p_typical(WittRing_general):
    pass

class WittRing_p_invertible(WittRing_non_p_typical):
    pass