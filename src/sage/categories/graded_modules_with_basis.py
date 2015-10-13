r"""
Graded modules with basis
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring

class GradedModulesWithBasis(GradedModulesCategory):
    """
    The category of graded modules with a distinguished basis.

    EXAMPLES::

        sage: C = GradedModulesWithBasis(ZZ); C
        Category of graded modules with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of graded modules over Integer Ring,
         Category of modules with basis over Integer Ring]
        sage: C is ModulesWithBasis(ZZ).Graded()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    class ParentMethods:

        # TODO: which syntax do we prefer?
        # A.basis(degree = 3)
        # A.basis().subset(degree=3)

        # This is related to the following design question:
        # If F = (f_i)_{i\in I} is a family, should ``F.subset(degree = 3)``
        # be the elements of F of degree 3 or those whose index is of degree 3?

        def basis(self, d=None):
            """
            Return the basis for (an homogeneous component of) ``self``.

            INPUT:

            - `d` -- non negative integer or ``None``, optional (default: ``None``)

            If `d` is None, returns a basis of the module.
            Otherwise, returns the basis of the homogeneous component of degree `d`.

            EXAMPLES::

                sage: A = GradedModulesWithBasis(ZZ).example()
                sage: A.basis(4)
                Lazy family (Term map from Partitions to An example of a graded module with basis: the free module on partitions over Integer Ring(i))_{i in Partitions of the integer 4}

            Without arguments, the full basis is returned::

                sage: A.basis()
                Lazy family (Term map from Partitions to An example of a graded module with basis: the free module on partitions over Integer Ring(i))_{i in Partitions}
                sage: A.basis()
                Lazy family (Term map from Partitions to An example of a graded module with basis: the free module on partitions over Integer Ring(i))_{i in Partitions}
            """
            if d is None:
                from sage.sets.family import Family
                return Family(self._indices, self.monomial)
            else:
                return self.graded_component_basis(d)

        def graded_component_basis(self, d):
            """
            Return a basis for the ``d``-th graded component of ``self``.

            EXAMPLES::

                sage: A = GradedModulesWithBasis(ZZ).example()
                sage: A.graded_component_basis(4)
                Lazy family (Term map from Partitions to An example of a graded module with basis: the free module on partitions over Integer Ring(i))_{i in Partitions of the integer 4}
            """
            from sage.sets.family import Family
            return Family(self._indices.subset(size=d), self.monomial)

        def graded_component(self, d):
            """
            Return the ``d``-th graded component of ``self``.

            EXAMPLES::

                sage: A = GradedModulesWithBasis(ZZ).example()
                sage: A.graded_component(4)
                Degree 4 component of An example of a graded module with
                 basis: the free module on partitions over Integer Ring
            """
            from sage.categories.modules_with_basis import ModulesWithBasis
            category = ModulesWithBasis(self.category().base_ring())
            M = self.submodule(self.graded_component_basis(d),
                               category=category,
                               already_echelonized=True)
            M.rename("Degree {} component of {}".format(d, self))
            return M

    class ElementMethods:
        def is_homogeneous(self):
            """
            Return whether this element is homogeneous.

            EXAMPLES::

                sage: A = GradedModulesWithBasis(ZZ).example()
                sage: x=A(Partition((3,2,1)))
                sage: y=A(Partition((4,4,1)))
                sage: z=A(Partition((2,2,2)))
                sage: (3*x).is_homogeneous()
                True
                sage: (x - y).is_homogeneous()
                False
                sage: (x+2*z).is_homogeneous()
                True
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

        def degree(self):
            """
            The degree of this element in the graded module.

            .. note::

                This raises an error if the element is not homogeneous.
                Another implementation option would be to return the
                maximum of the degrees of the homogeneous summands.

            EXAMPLES::

                sage: A = GradedModulesWithBasis(ZZ).example()
                sage: x = A(Partition((3,2,1)))
                sage: y = A(Partition((4,4,1)))
                sage: z = A(Partition((2,2,2)))
                sage: x.degree()
                6
                sage: (x + 2*z).degree()
                6
                sage: (y - x).degree()
                Traceback (most recent call last):
                ...
                ValueError: Element is not homogeneous.
            """
            if not self.support():
                raise ValueError("The zero element does not have a well-defined degree.")
            if self.is_homogeneous():
                return self.parent().degree_on_basis(self.leading_support())
            else:
                raise ValueError("Element is not homogeneous.")

        def homogeneous_component(self, n):
            """
            Return the homogeneous component of degree ``n`` of this
            element.

            EXAMPLES::

                sage: A = GradedModulesWithBasis(ZZ).example()
                sage: x = A.an_element(); x
                2*P[] + 2*P[1] + 3*P[2]
                sage: x.homogeneous_component(-1)
                0
                sage: x.homogeneous_component(0)
                2*P[]
                sage: x.homogeneous_component(1)
                2*P[1]
                sage: x.homogeneous_component(2)
                3*P[2]
                sage: x.homogeneous_component(3)
                0

            TESTS:

            Check that this really return ``A.zero()`` and not a plain ``0``::

                sage: x.homogeneous_component(3).parent() is A
                True
            """
            degree_on_basis = self.parent().degree_on_basis
            return self.parent().sum_of_terms((i, c)
                                              for (i, c) in self
                                              if degree_on_basis(i) == n)

        def truncate(self, n):
            """
            Return the sum of the homogeneous components of degree ``< n`` of this element

            EXAMPLES::

                sage: A = GradedModulesWithBasis(ZZ).example()
                sage: x = A.an_element(); x
                2*P[] + 2*P[1] + 3*P[2]
                sage: x.truncate(0)
                0
                sage: x.truncate(1)
                2*P[]
                sage: x.truncate(2)
                2*P[] + 2*P[1]
                sage: x.truncate(3)
                2*P[] + 2*P[1] + 3*P[2]

            TESTS:

            Check that this really return ``A.zero()`` and not a plain ``0``::

                sage: x.truncate(0).parent() is A
                True
            """
            degree_on_basis = self.parent().degree_on_basis
            return self.parent().sum_of_terms((i, c) for (i, c) in self
                                              if degree_on_basis(i) < n)

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):
        class ParentMethods:
            def graded_component_basis(self, d):
                """
                Return a basis for the ``d``-th graded component of ``self``.

                EXAMPLES::

                    sage: cat = GradedModulesWithBasis(ZZ).FiniteDimensional()
                    sage: C = CombinatorialFreeModule(ZZ, ['a', 'b'], category=cat)
                    sage: C.degree_on_basis = lambda x: 1 if x == 'a' else 2
                    sage: C.graded_component_basis(1)
                    Finite family {'a': B['a']}
                    sage: C.graded_component_basis(2)
                    Finite family {'b': B['b']}
                """
                from sage.sets.family import Family
                try:
                    S = self._indices.subset(size=d)
                except (AttributeError, ValueError, TypeError):
                    S = [i for i in self._indices if self.degree_on_basis(i) == d]
                return Family(S, self.monomial)

