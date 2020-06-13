r"""
Affine Vertex Algebra

AUTHORS:

- Reimundo Heluani (08-09-2019): Initial implementation.
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .universal_enveloping_vertex_algebra import \
                                                UniversalEnvelopingVertexAlgebra
from sage.sets.family import Family


class AffineVertexAlgebra(UniversalEnvelopingVertexAlgebra):
    def __init__(self, R, ct, k, names=None):
        r"""The universal affine vertex algebra over the ring `R` of the Lie
        algebra of finite Cartan Type `ct` at level `k`

        EXAMPLES::

            sage: V = AffineVertexAlgebra(QQ, 'A1', 1); V
            The universal affine vertex algebra of CartanType ['A', 1] at level 1
            sage: V = AffineVertexAlgebra(QQ, 'B3', 1); V
            The universal affine vertex algebra of CartanType ['B', 3] at level 1

        """
        from sage.algebras.lie_conformal_algebras.affine_lie_conformal_algebra\
                    import AffineLieConformalAlgebra
        if names is not None:
            prefix = ''
            bracket = ''
        else:
            prefix = 'E'
            bracket = '('
        ML = AffineLieConformalAlgebra(R, ct, names=names, prefix=prefix,
                                       bracket=bracket)
        cp = Family({ML.central_elements()[0]: k})
        super(AffineVertexAlgebra,self).__init__(R, ML,
            central_parameters = cp, names=names)

        self._level = k
        if type(ct) is str:
            from sage.combinat.root_system.cartan_type import CartanType
            ct = CartanType(ct)
        self._ct = ct
        self._c = k*self._ngens/(k+ct.dual_coxeter_number())

    def level(self):
        r"""The level of this Affine vertex algebra

        EXAMPLES:

        sage: V = AffineVertexAlgebra(QQ, 'B3', 1); V
        The universal affine vertex algebra of CartanType ['B', 3] at level 1
        sage: V.level()
        1

        """
        return self._level

    def cartan_type(self):
        r"""The Cartan Type of this Affine vertex algebra

        EXAMPLES::

            sage: V = AffineVertexAlgebra(QQ, 'B3', 1); V
            The universal affine vertex algebra of CartanType ['B', 3] at level 1
            sage: V.cartan_type()
            ['B', 3]

        """
        return self._ct

    def is_critical(self):
        r"""True if the level equals minus the dual Coxeter number of its Cartan
        Type

        EXAMPLES::

            sage: V = AffineVertexAlgebra(QQ, 'A1', -2); V
            The universal affine vertex algebra of CartanType ['A', 1] at critical level
            sage: V.is_critical()
            True

        """
        return self.level() == -self.cartan_type().dual_coxeter_number()

    def _repr_(self):
        if self.is_critical():
            return "The universal affine vertex algebra of CartanType {} at"\
                   " critical level over {}".format(self.cartan_type(),
                   self.base_ring())
        else:
            return "The universal affine vertex algebra of CartanType {} at "\
                   "level {} over {}".format(self.cartan_type(),
                                             self.level(), self.base_ring())

