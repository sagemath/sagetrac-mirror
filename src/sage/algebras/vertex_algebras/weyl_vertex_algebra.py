"""
Weyl Vertex Algebra

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

class WeylVertexAlgebra(UniversalEnvelopingVertexAlgebra):

    def __init__(self,R,gram_matrix=None,ngens=None,names=None,
                 index_set=None):
        r"""
        The Weyl vertex algebra.

        INPUT:

        - ``R``: a commutative ring.
        - ``gram_matrix``: a matrix (Default: `None`); A non-singular  
          skew-symmetric square matrix with coefficients in `R`. 
        - ``ngens``: an even positive Integer (Default `2`); The number
          of generators of this vertex algebra. 

        OUTPUT: 

        The Weyl Lie conformal algebra with generators
         `\alpha_i`, `i=1,...,ngens` and `\lambda`-brackets

        .. MATH::

            [{\alpha_i}_{\lambda} \alpha_j] = M_{ij} |0\rangle, 

        where `M` is the ``gram_matrix`` above. 

        .. NOTE::

            The returned vertex algebra is not `H`-graded. For
            a related `H`-graded vertex algebra see
            :class:`BosonicGhostsVertexAlgebra`
        """
        from sage.algebras.lie_conformal_algebras.weyl_lie_conformal_algebra import\
                                                    WeylLieConformalAlgebra
        ML = WeylLieConformalAlgebra(R, gram_matrix=gram_matrix, ngens=ngens,
                                      names=names, index_set=index_set)
        cp = Family({ML.gen(-1):R.one()})
        names = ML.variable_names()[:-1]
        super(WeylVertexAlgebra,self).__init__(R, ML,
                                            central_parameters=cp, names=names)

    def _repr_(self):
        return "The Weyl vertex algebra with generators {} over {}".format(
                self.gens(),self.base_ring())



