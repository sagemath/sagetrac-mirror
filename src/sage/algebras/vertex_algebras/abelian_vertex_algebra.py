"""
Abelian Vertex Algebra.

AUTHORS:

- Reimundo Heluani (06-15-2020): Initial implementation.
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

class AbelianVertexAlgebra(UniversalEnvelopingVertexAlgebra):
    def __init__(self, R, ngens=1, weights=None, parity=None, names=None,
                 index_set=None):
        r"""
        The Abelian vertex algebra.

        INPUT:

        - ``R`` -- a commutative ring.
        - ``ngens`` -- a positive integer (default: ``1``); the number
          of generators of this Lie conformal algebra
        - ``weights`` -- a list of numbers (default: ``1`` for each
          generator); the weights of the generators. The resulting 
          vertex algebra is `H`-graded.
        - ``parity`` -- ``None`` or a list of ``0`` or ``1`` (default:
           ``None``); The parity of the generators. If not ``None`` the
           resulting Lie Conformal algebra is a Super Lie conformal
           algebra


        OUTPUT:

        The Abelian vertex algebra with generators `a_i`,
        `i=1,...,ngens` and vanishing `\lambda`-brackets.

        EXAMPLES::

            sage: F = AbelianVertexAlgebra(QQ,2,weights=(2,3/2),parity=(0,1),names=('L','G'))
            sage: F.inject_variables()
            Defining L, G
            sage: (L*L)*L*G
            L_-2L_-2L_-2G_-3/2|0>

        .. TODO::

            implement its own class to speed up arithmetics in this
            case.
        """
        from sage.algebras.lie_conformal_algebras.abelian_lie_conformal_algebra\
             import AbelianLieConformalAlgebra

        ML = AbelianLieConformalAlgebra(R,ngens,weights,parity,names,index_set)
        super(AbelianVertexAlgebra,self).__init__(R,ML,names=names)

    def _repr_(self):
        return "The Abelian vertex algebra over {} with generators {}".\
               format(self.base_ring(), self.gens())

    from .vertex_algebra_element import UniversalEnvelopingVertexAlgebraElement
    class Element(UniversalEnvelopingVertexAlgebraElement):

        def _bracket_(self,other):
            return {}

