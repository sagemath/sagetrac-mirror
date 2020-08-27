r"""
Super Vertex Algebras

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

from .super_modules import SuperModulesCategory
from .proto_vertex_algebras import ProtoVertexAlgebras
from .graded_modules import GradedModulesCategory

class SuperVertexAlgebras(SuperModulesCategory):
    """
    The subcategory of super vertex algebras.

    EXAMPLES::

        sage: VertexAlgebras(QQbar).Super()
        Category of super vertex algebras over Algebraic Field
    """
    #Need to do all this to make Super commute with Graded.
    def extra_super_categories(self):
        """
        The extra super categories of this category.

        EXAMPLES::

            sage: VertexAlgebras(QQ).Super().super_categories()
            [Category of super Lie conformal algebras over Rational Field,
             Category of proto vertex algebras over Rational Field]
        """
        return [ProtoVertexAlgebras(self.base_ring())]

    def example(self):
        """
        An example of a super vertex algebra.

        EXAMPLES::

            sage: VertexAlgebras(QQ).Super().example()
            The Free Fermions super vertex algebra with generators (psi_-1/2|0>,) over Rational Field
        """
        from sage.algebras.vertex_algebras.free_fermions_vertex_algebra import \
                                                  FreeFermionsVertexAlgebra
        return FreeFermionsVertexAlgebra(self.base_ring())

    class Graded(GradedModulesCategory):
        """
        The category of H-graded super vertex algebras.

        EXAMPLES::

            sage: VertexAlgebras(QQ).Super().Graded()
            Category of H-graded super vertex algebras over Rational Field
        """
        def _repr_object_names(self):
            """
            The names of the objects of this category.

            EXAMPLES::

                sage: VertexAlgebras(QQbar).Super().Graded()
                Category of H-graded super vertex algebras over Algebraic Field
            """
            return "H-graded {}".format(self.base_category()._repr_object_names())
