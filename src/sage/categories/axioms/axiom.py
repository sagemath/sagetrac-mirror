"""
Abstract Base Class for Axioms
"""

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation


class Axiom(SageObject, UniqueRepresentation):

    def __init__(self):
        from factory import axioms
        self._sort_key = axioms._add_instance(self)

    def __cmp__(self, other):
        return cmp(self._sort_key, other._sort_key)

    def _repr_(self):
        return self.__class__.__name__
