"""
Abstract Base Class for Axioms
"""

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation


class Axiom(UniqueRepresentation, SageObject):
    
    @staticmethod
    def __classcall__(cls, *args, **kwds):
        if not hasattr(cls, '_sort_key'):
            raise TypeError('you must call axioms.register(new_axiom) first')
        return super(Axiom, cls).__classcall__(cls, *args, **kwds)
    
    def __cmp__(self, other):
        return cmp(self._sort_key, other._sort_key)

    def _repr_(self):
        return self.__class__.__name__
