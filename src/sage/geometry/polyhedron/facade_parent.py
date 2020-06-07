r"""
Facade parents representing infinite families of polyhedra
"""

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.geometry.polyhedron.parent import Polyhedra

class Polyhedra_facade_base(UniqueRepresentation, Parent):
    """
    Base class for facade parents representing infinite families of polyhedra.
    """

    def __init__(self, base_ring, ambient_dim, backend):
        """
        TESTS::

            sage: from sage.geometry.polyhedron.facade_parent import Polyhedra_facade_base
            sage: CP = Polyhedra_facade_base(QQ, 2, 'ppl')
            sage: TestSuite(CP).run()
        """
        self._polyhedra = polyhedra = Polyhedra(base_ring, ambient_dim, backend)
        Parent.__init__(self, facade=polyhedra)

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
    """
    Facade parent representing the family of full-dimensional polyhedra.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.facade_parent import FullDimensionalPolyhedra
        sage: Polyhedron(vertices=[[0, 0]]) in FullDimensionalPolyhedra(QQ, 2, 'field')
        False
    """

    def _repr_(self):
        return "Full-dimensional " + super(FullDimensionalPolyhedra, self)._repr_()

    def _check_polyhedron(self, polyhedron):
        super(FullDimensionalPolyhedra, self)._check_polyhedron(polyhedron)
        if not polyhedron.is_full_dimensional():
            raise ValueError("{} should be full-dimensional".format(polyhedron))

class LowerDimensionalPolyhedra(Polyhedra_facade_base):
    """
    Facade parent representing the family of lower-dimensional polyhedra.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.facade_parent import LowerDimensionalPolyhedra
        sage: Polyhedron(vertices=[[0, 0]]) in LowerDimensionalPolyhedra(QQ, 2, 'field')
        True
    """

    def _repr_(self):
        return "Lower-dimensional " + super(LowerDimensionalPolyhedra, self)._repr_()

    def _check_polyhedron(self, polyhedron):
        super(LowerDimensionalPolyhedra, self)._check_polyhedron(polyhedron)
        if polyhedron.is_full_dimensional():
            raise ValueError("{} should be lower-dimensional".format(polyhedron))

class NonPointedPolyhedra(Polyhedra_facade_base):
    """
    Facade parent representing the family of lower-dimensional polyhedra.

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
