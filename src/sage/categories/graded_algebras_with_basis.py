r"""
Graded algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import GradedAlgebras, GradedModulesWithBasis, AlgebrasWithBasis
from sage.misc.cachefunc import cached_method

class GradedAlgebrasWithBasis(Category_over_base_ring):
    """
    The category of graded algebras with a distinguished basis

    EXAMPLES::

        sage: GradedAlgebrasWithBasis(ZZ)
        Category of graded algebras with basis over Integer Ring
        sage: GradedAlgebrasWithBasis(ZZ).super_categories()
        [Category of graded modules with basis over Integer Ring, Category of graded algebras over Integer Ring, Category of algebras with basis over Integer Ring]

    TESTS::

        sage: TestSuite(GradedAlgebrasWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GradedAlgebrasWithBasis(QQ).super_categories()
            [Category of graded modules with basis over Rational Field, Category of graded algebras over Rational Field, Category of algebras with basis over Rational Field]
        """
        R = self.base_ring()
        return [GradedModulesWithBasis(R),GradedAlgebras(R), AlgebrasWithBasis(R)]

    def example(self, generators = ('a','b','c'), degrees = 1):
        """
        Returns an example of a graded algebra with basis as per
        :meth:`Category.example <sage.categories.category.Category.example>`.

        EXAMPLES::

            sage: GradedAlgebrasWithBasis(QQ).example()
            An example of a graded algebra with basis: the polynomial algebra on generators ('a', 'b', 'c') of degrees (1, 1, 1) over Rational Field

        Another set of generators, and their degrees, can be specified
        as optional arguments::

            sage: GradedAlgebrasWithBasis(QQ).example(generators=("x", "z"), degrees=(3, 10))
            An example of a graded algebra with basis: the polynomial algebra on generators ('x', 'z') of degrees (3, 10) over Rational Field
        """
        from sage.categories.examples.graded_algebras_with_basis import Example
        return Example(self.base_ring(), generators, degrees)

    class ParentMethods:

        # TODO: which syntax do we prefer?
        # A.basis(degree = 3)
        # A.basis().subset(degree=3)

        # This is related to the following design question:
        # If F = (f_i)_{i\in I} is a family, should ``F.subset(degree = 3)``
        # be the elements of F of degree 3 or those whose index is of degree 3?

        def basis(self, degree=None):
            """
            Returns the basis for (an homogeneous component of) this graded algebra

            INPUT:

            - `degree` -- non negative integer or ``None``, optional (default: ``None``)

            If `degree` is None, returns a basis of the algebra.
            Otherwise, returns the basis of the homogeneous component of degree
            `degree`.

            EXAMPLES::

                sage: A = GradedAlgebrasWithBasis(QQ).example(generators=("x", "y"), degrees=(2, 3))
                sage: A.basis(4)
                Lazy family (Term map from Integer vectors weighted by [2, 3] to An example of a graded algebra with basis: the polynomial algebra on generators ('x', 'y') of degrees (2, 3) over Rational Field(i))_{i in Integer vectors of 4 weighted by [2, 3]}
                sage: A.basis(4) # todo: not implemented (output)
                Family (x^{2},)
                sage: A.basis(6) # todo: not implemented (output)
                Family (y^{2}, x^{3})
                sage: A.basis(-10) # todo: not implemented (output)
                Family ()

            Without arguments, the full basis is returned::

                sage: A.basis()
                Lazy family (Term map from Integer vectors weighted by [2, 3] to An example of a graded algebra with basis: the polynomial algebra on generators ('x', 'y') of degrees (2, 3) over Rational Field(i))_{i in Integer vectors weighted by [2, 3]}
            """
            from sage.sets.family import Family
            if degree is None:
                return Family(self._basis_keys, self.monomial)
            else:
                return Family(self._basis_keys.subset(size=degree), self.monomial)

        @cached_method
        def homogeneous_component(self, degree):
            """
            Returns the degree `degree` homogeneous component of this graded algebra.

            EXAMPLES::

                sage: A = GradedAlgebrasWithBasis(QQ).example(generators=("x", "y"), degrees=(2, 3))
                sage: A.homogeneous_component(1)
                Free module generated by Integer vectors of 1 weighted by [2, 3] over Rational Field
                sage: A.homogeneous_component(6) # todo: not implemented (output)
                Free module generated by (y^{2}, x^{3}) over Rational Field

            As a shorthand, graded algebras may allow to access their
            homogeneous components using square brackets; but this is
            not systematically enforced::

                sage: A[9]
                Free module generated by Integer vectors of 9 weighted by [2, 3] over Rational Field

            TODO:

             - How do we want homogeneous components to be printed out?

             - If `A` and `B` are two graded algebras over the same index
               set and ground ring, should ``A.homogeneous_component(3)``
               and `B.homogeneous_component(3)`` coincide? Or should they
               instead be respectively aware that they are subsets of `A`
               and `B`?
            """
            from sage.combinat.free_module import CombinatorialFreeModule
            M = CombinatorialFreeModule(self.base_ring(),
                            self.basis().keys().subset(size=degree),
                            element_class=self.Element)
            #M._name = "Free module generated by %s"%(tuple(self.monomial(a) for a in basis),)
            return M

        # TODO: degree_on_basis should be an abstract method; write doc and doctests
        # from sage.misc.abstract_method import abstract_method
        #@abstract_method
        #def degree_on_basis(self):
        #    pass

    class ElementMethods:
        def is_homogeneous(self):
            """
            Return whether this element is homogeneous.

            EXAMPLES::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: (3*x).is_homogeneous()
                True
                sage: (x^3 - y^2).is_homogeneous()
                True
                sage: ((x + y)^2).is_homogeneous()
                False
            """
            degree_on_basis = self.parent().degree_on_basis
            degree = None
            for m in self.support():
                if degree is None:
                    degree = degree_on_basis(m)
                else:
                    if degree != degree_on_basis(m):
                        return False
            return True

        def homogeneous_degree(self):
            """
            The degree of this element.

            .. note::

               This raises an error if the element is not homogeneous.
               To obtain the maximum of the degrees of the homogeneous
               summands, use :meth:`maximal_degree`

            .. seealso: :meth:`maximal_degree`

            EXAMPLES::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: x.homogeneous_degree()
                2
                sage: (x^3 + 4*y^2).homogeneous_degree()
                6
                sage: ((1 + x)^3).homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: Element is not homogeneous.

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S.zero().degree()
                Traceback (most recent call last):
                ...
                ValueError: The zero element does not have a well-defined degree.
            """
            if self.is_zero():
                raise ValueError("The zero element does not have a well-defined degree.")
            try:
                assert self.is_homogeneous()
                return self.parent().degree_on_basis(self.leading_support())
            except AssertionError:
                raise ValueError("Element is not homogeneous.")

        # default choice for degree; will be overridden as necessary
        degree = homogeneous_degree

        def maximal_degree(self):
            """
            The maximum of the degrees of the homogeneous summands.

            .. seealso: :meth:`homogeneous_degree`

            EXAMPLES::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: x.maximal_degree()
                2
                sage: (x^3 + 4*y^2).maximal_degree()
                6
                sage: ((1 + x)^3).maximal_degree()
                6

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S.zero().degree()
                Traceback (most recent call last):
                ...
                ValueError: The zero element does not have a well-defined degree.
            """
            if self.is_zero():
                raise ValueError("The zero element does not have a well-defined degree.")
            else:
                degree_on_basis = self.parent().degree_on_basis
                return max(degree_on_basis(m) for m in self.support())
