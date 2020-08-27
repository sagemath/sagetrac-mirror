r"""
Graded Proto Vertex Algebras

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
from sage.misc.abstract_method import abstract_method
from sage.rings.all import ZZ

class GradedProtoVertexAlgebras(GradedModulesCategory):
    """
    The category of H-graded proto vertex algebras.

    This is an abstract base category for H-graded vertex and super
    vertex algebras.
    """
    class ElementMethods:
        """
        Base class for elements of an H-graded vertex algebra.
        """
        @abstract_method
        def degree(self):
            """
            The conformal weight of this element.

            EXAMPLES::

                sage: V = vertex_algebras.Virasoro(QQ,1/2); L = V.0
                sage: L.degree()
                2
                sage: W = vertex_algebras.Affine(QQ, 'A1', 2); E = W.0
                sage: E.degree()
                1
                sage: L.T().degree()
                3
                sage: (L + L.T()).degree()
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous
            """

        def weight(self):
            """
            The conformal weight of this element.

            EXAMPLES::

                sage: V = vertex_algebras.NeveuSchwarz(QQ,1);
                sage: V.gens()
                (L_-2|0>, G_-3/2|0>)
                sage: [g.weight() for g in V.gens()]
                [2, 3/2]
            """
            return self.degree()

        @abstract_method
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

        def nmodeproduct(self, other, n):
            r"""
            The shifted product of these two elements.

            INPUT:

            - ``other`` -- an element of this vertex algebra
            - ``n`` -- a rational.

            OUTPUT:

            The shifted `n`-product of ``self`` with ``other``, which is
            defined as follows.

            For an element `a` of degree `p`, that is `a \in V_p`,
            then `a_n b` is defined as `a_{(n+p-1)}b`.

            .. NOTE::

                For this method to be defined, the element ``self``
                needs to be homogeneous of rational conformal weight
                say ``w``. In addition, ``n+w`` needs to be an
                integer number. If these conditions are not met
                this method raises ``ValueError``.

            EXAMPLES::

                sage: V = vertex_algebras.Virasoro(QQ, 1/2); V.register_lift()
                sage: V.inject_variables()
                Defining L
                sage: L.nmodeproduct(L.T(),0)
                3*L_-3|0>

                sage: (L + V.vacuum()).nmodeproduct(L,0)
                Traceback (most recent call last):
                ...
                ValueError: couldn't compute weight of |0> + L_-2|0>, is it not homogeneous?
            """
            try:
                weight = self.weight()
            except ValueError:
                raise ValueError("couldn't compute weight of {}, "\
                                "is it not homogeneous?".format(self))
            if n+weight not in ZZ:
                raise ValueError("{} does not have a mode {}".format(self,
                                                                         n))
            return self._nproduct_(other, n+weight-1)
