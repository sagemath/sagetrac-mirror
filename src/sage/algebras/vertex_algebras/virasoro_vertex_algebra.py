"""
Virasoro Vertex Algebra

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

class VirasoroVertexAlgebra(UniversalEnvelopingVertexAlgebra):
    def __init__(self, R, c=0, names="L"):
        r"""
        The universal Virasoro vertex algebra

        INPUT:

        - ``R`` a commutative ring; the base ring of this vertex
           algebra. Undefined behaviour if
           this is not a Field of characteristic zero.

        - ``c`` a number (default: ``0``); the central charge of this
           vertex algebra.

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2); V
            The Virasoro vertex algebra of central charge 1/2

        """
        from sage.algebras.lie_conformal_algebras.\
             virasoro_lie_conformal_algebra import VirasoroLieConformalAlgebra
        ML = VirasoroLieConformalAlgebra(R)
        cp = Family({ML.gen(1):c})
        super(VirasoroVertexAlgebra,self).__init__(R, ML,
                 central_parameters=cp, names=names)
        self._c = c

    def _repr_(self):
        return "The Virasoro vertex algebra of central charge {} over {}".\
                format(self.central_charge(),self.base_ring())



