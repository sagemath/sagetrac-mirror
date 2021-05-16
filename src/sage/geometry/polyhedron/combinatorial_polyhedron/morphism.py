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

    def __call__(self, p):
        """
        Apply the morphism to a combinatorial polyhedron.

        EXAMPLES::

            sage: C = polytopes.cube().combinatorial_polyhedron()
            sage: C.Vrepresentation()
            (A vertex at (1, -1, -1),
            A vertex at (1, 1, -1),
            A vertex at (1, 1, 1),
            A vertex at (1, -1, 1),
            A vertex at (-1, -1, 1),
            A vertex at (-1, -1, -1),
            A vertex at (-1, 1, -1),
            A vertex at (-1, 1, 1))
            sage: C.Hrepresentation()
            (An inequality (-1, 0, 0) x + 1 >= 0,
            An inequality (0, -1, 0) x + 1 >= 0,
            An inequality (0, 0, -1) x + 1 >= 0,
            An inequality (1, 0, 0) x + 1 >= 0,
            An inequality (0, 0, 1) x + 1 >= 0,
            An inequality (0, 1, 0) x + 1 >= 0)

        Abstraction morphism::

            sage: Vrep_dict = {vertex: index for index, vertex in enumerate(C.vertices())}
            sage: phi = C.parent().hom(Vrep_dict)
            sage: phi_C = phi(C); phi_C
            A 3-dimensional combinatorial polyhedron with 6 facets
            sage: phi_C.Vrepresentation()
            (0, 1, 2, 3, 4, 5, 6, 7)
            sage: phi_C.Hrepresentation()
            (0, 1, 2, 3, 4, 5)

        """
        # For now, degenerations are not implemented
        incidence_matrix = p.incidence_matrix()

        Vrep = p.Vrepresentation()
        def map_vrep_element(v):
            return self._Vrep_dict[v]
        Vrep = [map_vrep_element(v) for v in Vrep]

        return self.codomain()(incidence_matrix, Vrep=Vrep)
