r"""
Poisson Vertex Algebras.

AUTHORS:

- Reimundo Heluani (08-09-2019): Initial implementation.

.. include:: ../../../poisson_vertex_algebras/poisson_vertex_algebra_desc.rst
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


from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.vertex_algebras import VertexAlgebras
from sage.categories.poisson_vertex_algebras import PoissonVertexAlgebras
from sage.categories.commutative_rings import CommutativeRings

class PoissonVertexAlgebra(UniqueRepresentation):
    """
    class
    """
    @staticmethod
    def __classcall_private__(cls, R=None, arg0=None, category=None, **kwds):
        if not R in CommutativeRings():
            raise ValueError("R must be a commutative ring, got {}".format(R))

        #This is the only exposed class so we clean keywords here
        known_keywords = ['term_order']

        for key in kwds:
            if key not in known_keywords:
                raise ValueError("PoissonVertexAlgebra(): got an unexpected " +
                                "keyword argument '%s'"%key)

        if kwds:
            raise NotImplementedError("{} is not implemented yet".format(
                                      kwds.keys()))

        category = PoissonVertexAlgebras(R).or_subcategory(category)

        #Until ModulesWithBasis(R) has a working `change_ring` method
        #We need to check the base ring of arg0 is the same as R
        if arg0 in VertexAlgebras(R).Graded().FinitelyGenerated().WithBasis():
            from .vertex_algebra_classical_limit import \
                                                    VertexAlgebraClassicalLimit
            category = category.Graded().FinitelyGenerated().WithBasis()
            if arg0.is_super():
                category = category.Super()
            return VertexAlgebraClassicalLimit(R, arg0, category=category)

        raise ValueError ("arg0 needs to be a finitely generated H-graded"\
                        "vertex algebra with basis over R, got {}".format(arg0))

