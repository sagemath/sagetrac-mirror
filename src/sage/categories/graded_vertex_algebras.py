r"""
Graded Vertex Algebras

AUTHORS:

- Reimundo Heluani (2020-08-26): Initial implementation.
"""
#******************************************************************************
#       Copyright (C) 2020 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .graded_modules import GradedModulesCategory
from sage.misc.cachefunc import cached_method

class GradedVertexAlgebrasCategory(GradedModulesCategory):
    @cached_method
    def Super(self, base_ring=None):
        r"""
        Return the super-analogue category of ``self``.

        INPUT:

        - ``base_ring`` -- this is ignored

        EXAMPLES::

            sage: C = VertexAlgebras(QQbar)
            sage: C.Graded().Super() is C.Super().Graded()
            True
            sage: Cp = C.WithBasis()
            sage: Cp.Graded().Super() is Cp.Super().Graded()
            True
        """
        return self.base_category().Super(base_ring).Graded()

    def _repr_object_names(self):
        """
        The names of the objects of ``self``.

        EXAMPLES::

            sage: VertexAlgebras(QQbar).Graded()
            Category of H-graded vertex algebras over Algebraic Field

            sage: VertexAlgebras(QQbar).WithBasis().FinitelyGenerated().Graded()
            Category of H-graded finitely generated vertex algebras with basis over Algebraic Field
        """
        return "H-graded {}".format(self.base_category()._repr_object_names())

class GradedVertexAlgebras(GradedVertexAlgebrasCategory):
    """
    The category of H-graded vertex algebras.

    EXAMPLES::

        sage: VertexAlgebras(QQ).Graded()
        Category of H-graded vertex algebras over Rational Field
    """
