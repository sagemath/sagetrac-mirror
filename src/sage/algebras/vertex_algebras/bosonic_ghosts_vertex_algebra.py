"""
Bosonic Ghosts Vertex Algebra

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

class BosonicGhostsVertexAlgebra(UniversalEnvelopingVertexAlgebra):
    def __init__(self,R,ngens=2,names=None,index_set=None):
        r"""
        The Bosonic ghosts or `\beta-\gamma`-system vertex algebra.

        INPUT:

        - ``R``: a commutative ring.
        - ``ngens``: an even positive Integer (default: ``2``); the
          number generators of this vertex algebra.
        """
        from sage.algebras.lie_conformal_algebras.\
                    bosonic_ghosts_lie_conformal_algebra import \
                     BosonicGhostsLieConformalAlgebra

        ML = BosonicGhostsLieConformalAlgebra(R, ngens=ngens, names=names,
                                              index_set=index_set)

        cp = Family({ML.gen(-1):R.one()})
        names = ML.variable_names()[:-1]
        super(BosonicGhostsVertexAlgebra,self).__init__(R, ML,
                                            central_parameters=cp, names=names)
        self._c = self._ngens

    def _repr_(self):
        return "The Bosonic Ghosts vertex algebra with generators {}".\
                format(self.gens())





