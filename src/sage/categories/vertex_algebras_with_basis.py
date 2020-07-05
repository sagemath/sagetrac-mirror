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
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.super_modules import SuperModulesCategory

class VertexAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The subcategory of vertex algebras with a preferred basis.

    EXAMPLES::

        sage: VertexAlgebras(ZZ).WithBasis()
        Category of vertex algebras with basis over Integer Ring
    """
    class ElementMethods:
        """
        Base class of elements of a vertex algebra with basis.
        """
        def is_monomial(self):
            """
            Whether this element is a monomial.

            EXAMPLES::

                sage: V = vertex_algebras.Virasoro(QQ,1); V.inject_variables()
                Defining L
                sage: L.T(2)*L
                2*L_-4L_-2|0>
                sage: (L.T(2)*L).is_monomial()
                True
                sage: L*L.T(2)
                2*L_-4L_-2|0> + 4*L_-6|0>
                sage: (L*L.T(2)).is_monomial()
                False
            """
            return (len(self.monomial_coefficients()) == 1 or\
                    self.is_zero())

        def index(self):
            """
            The basis index parametrizing this monomial element.

            INPUT:

            - ``self`` -- a monomial in this vertex algebra.

            EXAMPLES::

                sage: V = vertex_algebras.Affine(QQ, 'A1', 1, names = ('e','h', 'f'));
                sage: V.inject_variables()
                Defining e, h, f
                sage: v = e.T(3)*(e.T(2)*(e*(h.T()*(h*f)))); v
                12*e_-4e_-3e_-1h_-2h_-1f_-1|0>
                sage: v.index()
                ([4, 3, 1], [2, 1], [1])

                sage: v = f*e; v.index()
                Traceback (most recent call last):
                ...
                ValueError: index can only be computed for monomials, got e_-1f_-1|0> - h_-2|0>
            """
            if self.is_zero():
                return None
            if not self.is_monomial():
                raise ValueError ("index can only be computed for monomials,"\
                                  " got {}".format(self))
            return next(iter(self.monomial_coefficients()))

    class Super(SuperModulesCategory):
        """
        The subcategory of super vertex algebras with basis.

        EXAMPLES::

            sage: VertexAlgebras(AA).WithBasis().Super()
            Category of super vertex algebras with basis over Algebraic Real Field
        """
        def extra_super_categories(self):
            """
            The extra super categories of this category

            EXAMPLES::

                sage: VertexAlgebras(QQ).WithBasis().Super().super_categories()
                [Category of super Lie conformal algebras with basis over Rational Field,
                Category of vertex algebras with basis over Rational Field,
                Category of super vertex algebras over Rational Field]
            """
            return [self.base_category(),]

    class Graded(GradedModulesCategory):
        """
        The subcategory of H-graded vertex algebras with basis.

        EXAMPLES::

            sage: VertexAlgebras(QQbar).WithBasis().Graded()
            Category of H-graded vertex algebras with basis over Algebraic Field
        """
        def _repr_object_names(self):
            """
            The names of the objects of this category.

            EXAMPLES::

                sage: VertexAlgebras(QQbar).WithBasis().Graded()
                Category of H-graded vertex algebras with basis over Algebraic Field
            """
            return "H-graded {}".format(self.base_category()._repr_object_names())

        class ParentMethods:
            """
            Base class of Parents of an H-graded vertex algebra with
            basis.
            """
            def degree_on_basis(self, m):
                """
                The conformal weight of the basis element indexed by
                ``m``.

                EXAMPLES::

                    sage: V = vertex_algebras.Virasoro(QQ,1)
                    sage: B = V._indices
                    sage: V.degree_on_basis(B([[3,2,1]]))
                    9                    
                """
                assert m in self._indices
                return m.energy()

        class ElementMethods:
            """
            Base class for elements of an H-graded vertex algebra with
            basis.
            """
            def homogeneous_terms(self):
                """
                Return a tuple with the homogeneous terms of this
                element with respect to conformal weight.

                EXAMPLES::

                    sage: V = vertex_algebras.Virasoro(QQ,1/2)
                    sage: V.inject_variables()
                    Defining L
                    sage: (L+L.T(3) + L*L.T() + L.T(7)/factorial(7)).homogeneous_terms()
                    (L_-2|0>, L_-3L_-2|0> + 7*L_-5|0>, L_-9|0>)

                TESTS::

                    sage: V = vertex_algebras.Virasoro(QQ,1); V(0).homogeneous_terms()
                    (0,)
                """
                if self.is_zero():
                    return tuple([self])
                S = {}
                p = self.parent()
                for m in self.terms():
                    w = m.weight()
                    S[w] = S.get(w,p.zero()) + m
                return tuple(S.values())

        class Super(SuperModulesCategory):
            """
            The subcategory of super H-graded vertex algebras with
            basis.

            EXAMPLES::

                sage: C = VertexAlgebras(QQbar).WithBasis()
                sage: C.Graded().Super()
                Category of super H-graded vertex algebras with basis over Algebraic Field
                sage: C.Graded().Super() is C.Super().Graded()
                True
            """
            def extra_super_categories(self):
                """
                The extra super categories of this category.

                EXAMPLES::

                    sage: VertexAlgebras(QQ).WithBasis().Graded().Super().super_categories()
                    [Category of super vertex algebras with basis over Rational Field,
                    Category of super H-graded Lie conformal algebras with basis over Rational Field,
                    Category of H-graded vertex algebras with basis over Rational Field,
                    Category of super H-graded vertex algebras over Rational Field]
                """
                return [self.base_category(),]
