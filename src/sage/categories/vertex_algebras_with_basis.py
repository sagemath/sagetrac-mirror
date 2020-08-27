"""
Vertex Algebras With Basis

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

class VertexAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The subcategory of vertex algebras with a preferred basis.

    EXAMPLES::

        sage: VertexAlgebras(ZZ).WithBasis()
        Category of vertex algebras with basis over Integer Ring
    """
    class Super(SuperModulesCategory):
        """
        The subcategory of super vertex algebras with basis.

        EXAMPLES::

            sage: VertexAlgebras(AA).WithBasis().Super()
            Category of super vertex algebras with basis over Algebraic Real Field
        """
        class Graded(GradedVertexAlgebrasCategory):
            """
            The category of H-graded super vertex algebras with basis.

            EXAMPLES::

                sage: VertexAlgebras(QQbar).WithBasis().Super().Graded()
                Category of H-graded super vertex algebras with basis over Algebraic Field
            """

    class Graded(GradedVertexAlgebrasCategory):
        """
        The subcategory of H-graded vertex algebras with basis.

        EXAMPLES::

            sage: VertexAlgebras(QQbar).WithBasis().Graded()
            Category of H-graded vertex algebras with basis over Algebraic Field
        """

