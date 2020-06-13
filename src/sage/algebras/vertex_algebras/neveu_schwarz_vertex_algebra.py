"""
Neveu-Schwarz Super Vertex Algebra

AUTHORS:

- Reimundo Heluani (06-09-2020): Initial implementation.
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

class NeveuSchwarzVertexAlgebra(UniversalEnvelopingVertexAlgebra):
    def __init__(self,R,c=0,names=('L','G')):
        """
        The Neveu-Schwarz super vertex algebra.

        INPUT:

        - ``R`` a commutative ring; the base ring. 
        - ``c`` a number (default: ``0``); the central charge.
        
        EXAMPLES::

        """
        from sage.algebras.lie_conformal_algebras.\
             neveu_schwarz_lie_conformal_algebra import \
             NeveuSchwarzLieConformalAlgebra

        ML = NeveuSchwarzLieConformalAlgebra(R,names=names)
        cp = Family({ML.gen(2):c})
        super(NeveuSchwarzVertexAlgebra,self).__init__(R, ML,
                 central_parameters=cp, names=names)
        self._c = c

    def _repr_(self):
        return "The Neveu-Schwarz super vertex algebra of central charge {}"\
               " over {}".format(self.central_charge(),self.base_ring())



