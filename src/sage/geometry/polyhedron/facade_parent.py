r"""
Facade parents representing infinite families of polyhedra
"""

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.geometry.polyhedron.parent import Polyhedra

class Polyhedra_facade_base(UniqueRepresentation, Parent):

    def __init__(self, base_ring, ambient_dim, backend):
        self._polyhedra = Polyhedra(base_ring, ambient_dim, backend)
        Parent.__init__(self, facade=self._polyhedra)

    def _element_constructor_(self, *args, **kwds):
        r"""
        Construction of elements

        This is also used for the membership test.
        """
        polyhedron = self._polyhedra._element_constructor_(*args, **kwds)
        self._check_polyhedron(polyhedron)
        return polyhedron

    def _check_polyhedron(self, polyhedron):
        pass

    def _repr_(self):
        return "Polyhedra in " + self._polyhedra._repr_ambient_module()

    def universe(self):
        return self._polyhedra.universe()

    def ambient_dim(self):
        return self._polyhedra.ambient_dim()

class FullDimensionalPolyhedra(Polyhedra_facade_base):

    def _repr_(self):
        return "Full-dimensional " + super(FullDimensionalPolyhedra, self)._repr_()

    def _check_polyhedron(self, polyhedron):
        super(FullDimensionalPolyhedra, self)._check_polyhedron(polyhedron)
        if not polyhedron.is_full_dimensional():
            raise ValueError("{} should be full-dimensional".format(polyhedron))

class LowerDimensionalPolyhedra(Polyhedra_facade_base):

    def _repr_(self):
        return "Lower-dimensional " + super(LowerDimensionalPolyhedra, self)._repr_()

    def _check_polyhedron(self, polyhedron):
        super(LowerDimensionalPolyhedra, self)._check_polyhedron(polyhedron)
        if polyhedron.is_full_dimensional():
            raise ValueError("{} should be lower-dimensional".format(polyhedron))

class NonPointedPolyhedra(Polyhedra_facade_base):
    """
    EXAMPLES::

        sage: from sage.geometry.polyhedron.facade_parent import NonPointedPolyhedra
        sage: Polyhedron(vertices=[[0, 0]]) in NonPointedPolyhedra(QQ, 2, 'field')
        False
        sage: Polyhedron(lines=[[1, 1]]) in NonPointedPolyhedra(QQ, 2, 'field')
        True
    """
    def _repr_(self):
        return "Non-pointed " + super(NonPointedPolyhedra, self)._repr_()

    def _check_polyhedron(self, polyhedron):
        super(NonPointedPolyhedra, self)._check_polyhedron(polyhedron)
        if not polyhedron.lines():
            raise ValueError("{} should be non-pointed".format(polyhedron))
