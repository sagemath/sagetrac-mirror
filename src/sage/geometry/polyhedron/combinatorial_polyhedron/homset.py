r"""
Homsets between sets of combinatorial polyhedra
"""

from sage.categories.homset import Homset
from .morphism import CombinatorialPolyhedraMorphism

class CombinatorialPolyhedraHomset(Homset):

    Element = CombinatorialPolyhedraMorphism

    def _element_constructor_(self, Vrep_dict, check=None, **options):

        return self.element_class(Vrep_dict, parent=self)
