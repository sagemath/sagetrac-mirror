"""
Free Bosons Vertex Algebra

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

class FreeBosonsVertexAlgebra(UniversalEnvelopingVertexAlgebra):
    def __init__(self,R,gram_matrix=None,ngens=None,names=None,
                 index_set=None):
        r"""
        The Free Bosons Lie conformal algebra.

        INPUT:

        - ``R``: a commutative ring.
        - ``gram_matrix``: a matrix (default: ``[1]``); a symmetric
          square matrix with coefficients in ``R``.
        - ``ngens``: a positive Integer (default ``1``); the number of
          generators of this vertex algebra.

        OUTPUT:

        The Free Bosons Lie conformal algebra with generators
         `\alpha_i`, `i=1,...,n` and `\lambda`-brackets

         .. MATH::

            [{\alpha_i}_{\lambda} \alpha_j] = \lambda M_{ij} |0\rangle,

        where `n` is the number of generators ``ngens`` and `M` is
        the ``gram_matrix``. This vertex
        algebra is `H`-graded where every generator has conformal weight
        `1`.
        """
        from sage.algebras.lie_conformal_algebras.\
                    free_bosons_lie_conformal_algebra import \
                     FreeBosonsLieConformalAlgebra

        ML = FreeBosonsLieConformalAlgebra(R,gram_matrix=gram_matrix,
                        ngens=ngens,names=names,index_set=index_set)

        cp = Family({ML.gen(-1):R.one()})

        names = ML.variable_names()[:-1]

        super(FreeBosonsVertexAlgebra,self).__init__(R, ML,
                 central_parameters=cp, names=names)
        self._c = self._ngens

    def _repr_(self):
        return "The Free Bosons vertex algebra with generators {}".\
                format(self.gens())





