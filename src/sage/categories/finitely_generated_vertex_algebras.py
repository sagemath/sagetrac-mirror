"""
Finitely Generated Vertex Algebras

AUTHORS:

- Reimundo Heluani (2019-10-09): Initial implementation.
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

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from .graded_vertex_algebras import GradedVertexAlgebrasCategory
from sage.categories.super_modules import SuperModulesCategory
from .vertex_algebras import VertexAlgebras

class FinitelyGeneratedVertexAlgebras(CategoryWithAxiom_over_base_ring):
    """
    The subcategory of finitely generated vertex algebras.

    EXAMPLES::

        sage: VertexAlgebras(QQ).FinitelyGenerated()
        Category of finitely generated vertex algebras over Rational Field
    """
    _base_category_class_and_axiom = (VertexAlgebras, "FinitelyGeneratedAsProtoVertexAlgebra")
    class WithBasis(CategoryWithAxiom_over_base_ring):
        """
        The category of finitely generated vertex algebras with basis.

        EXAMPLES::

            sage: VertexAlgebras(QQ).FinitelyGenerated().WithBasis()
            Category of finitely generated vertex algebras with basis over Rational Field
        """

        class Super(SuperModulesCategory):
            """
            The subcategory of super finitely generated vertex algebras
            with basis.

            EXAMPLES::

                sage: VertexAlgebras(AA).FinitelyGenerated().WithBasis().Super()
                Category of super finitely generated vertex algebras with basis over Algebraic Real Field
            """

        class Graded(GradedVertexAlgebrasCategory):
            """
            The subcategory of H-graded finitely generated vertex
            algebras with basis.

            .. NOTE::

                To implement a finitely generated H-graded
                vertex algebra with basis. The basis indices
                should implement a method `subset` which
                admits at least the keyword parameter
                ``energy`` such that
                ``_indices.subset(energy=n)`` enumerates
                the indices of the basis of monomials with
                conformal weight ``n``. See for example
                :mod:`EnergyPartitionTuples<sage.algebras.vertex_algebras.energy_partition_tuples>`.

            EXAMPLES::

                sage: VertexAlgebras(QQbar).FinitelyGenerated().WithBasis().Graded()
                Category of H-graded finitely generated vertex algebras with basis over Algebraic Field
            """
    class Super(SuperModulesCategory):
        """
        The subcategory of super finitely generated vertex algebras.

        EXAMPLES::

            sage: VertexAlgebras(AA).FinitelyGenerated().Super()
            Category of super finitely generated vertex algebras over Algebraic Real Field
        """

    class Graded(GradedVertexAlgebrasCategory):
        """
        The subcategory of H-graded finitely generated vertex algebras.

        EXAMPLES::

            sage: VertexAlgebras(QQbar).FinitelyGenerated().Graded()
            Category of H-graded finitely generated vertex algebras over Algebraic Field
        """

