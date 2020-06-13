"""
Free Fermions Vertex Algebra

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

class FreeFermionsVertexAlgebra(UniversalEnvelopingVertexAlgebra):
    def __init__(self,R,ngens=None,gram_matrix=None,names=None,
                 index_set=None):
        r"""
        The Free Fermions Super vertex algebra.

        INPUT:

        - ``R``: a commutative ring.
        - ``ngens``: a positive Integer (default ``1``); the number of
          generators of this vertex algebra.
        - ``gram_matrix``: a matrix (default: ``[1]``); a symmetric
          square matrix with coefficients in ``R``.

        OUTPUT:

        The Free Fermions super vertex algebra with generators
         `\psi_i`, `i=1,...,n` and `\lambda`-brackets

         .. MATH::

            [{\psi_i}_{\lambda} \psi_j] = \lambda M_{ij} |0\rangle,

        where `n` is the number of generators ``ngens`` and `M` is the
        ``gram_matrix``. This vertex
        algebra is `H`-graded where every generator has conformal weight
        `1/2`.

        EXAMPLES::

        """
        from sage.algebras.lie_conformal_algebras.\
                    free_fermions_lie_conformal_algebra import \
                     FreeFermionsLieConformalAlgebra

        ML = FreeFermionsLieConformalAlgebra(R,gram_matrix=gram_matrix,
                        ngens=ngens,names=names,index_set=index_set)

        cp = Family({ML.gen(-1):R.one()})

        names = ML.variable_names()[:-1]

        super(FreeFermionsVertexAlgebra,self).__init__(R, ML,
                 central_parameters=cp, names=names)
        self._c = self._ngens/2

    def _repr_(self):
        return "The Free Fermions super vertex algebra with generators {} over"\
               " {}".format(self.gens(),self.base_ring())



