r"""
Morphisms of sets of combinatorial polyhedra
"""

from sage.categories.morphism import Morphism
from sage.categories.homset import Hom

class CombinatorialPolyhedraMorphism(Morphism):
    """
    Element class for morphisms of sets of combinatorial polyhedra.
    """
    def __init__(self, Vrep_dict=None, parent=None, domain=None, codomain=None, check=True):
        self._Vrep_dict = Vrep_dict
        if parent is None:
            parent = Hom(domain, codomain)
        super().__init__(parent)

    def _repr_type(self):
        return "Combinatorial polyhedral set"
