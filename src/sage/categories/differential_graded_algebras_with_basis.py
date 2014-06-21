"""
Differential graded algebras with basis
"""
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.category_types import ChainComplexes
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute

class DifferentialGradedAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The category of differential graded algebras with a distinguished basis.

    EXAMPLES::

        sage: from sage.categories.differential_graded_algebras_with_basis import DifferentialGradedAlgebrasWithBasis
        sage: DifferentialGradedAlgebrasWithBasis(QQ)
        Category of differential graded algebras with basis over Rational Field
        sage: DifferentialGradedAlgebrasWithBasis(QQ).super_categories()
        [Category of graded algebras with basis over Rational Field,
         Category of chain complexes over Rational Field]

    TESTS::

        sage: TestSuite(DifferentialGradedAlgebrasWithBasis(QQ)).run()
    """
    _base_category_class_and_axiom = (GradedAlgebrasWithBasis, "Differential")

    def extra_super_categories(self):
        r"""
        Return the :class:`ChainComplexes` category.

        This method specifies that a differential graded algebra with
        a basis is a chain complex.

        .. SEEALSO::

            The :ref:`axioms-deduction-rules` section in the
            documentation of axioms

        EXAMPLES::

            sage: C = Algebras(QQ).WithBasis().Graded().Differential()
            sage: C.extra_super_categories()
            (Category of chain complexes over Rational Field,)
        """
        return (ChainComplexes(self.base_ring()),)

    def example(self):
        """
        An example of differential graded algebra with basis.

        EXAMPLES::

            sage: from sage.categories.differential_graded_algebras_with_basis import DifferentialGradedAlgebrasWithBasis
            sage: DifferentialGradedAlgebrasWithBasis(QQ).example()
            Free commutative differential graded algebra over Rational Field on generators in degrees 1, 1, 1, 1, 1
        """
        from sage.algebras.differential_graded_algebras.commutative_dga import CommutativeDGA
        return CommutativeDGA(degrees=(1,1,1,1,1), differential=(None, None, None, (((1,1,0,0,0), 1),), (((0,1,1,0,0), 1),)))

    class ParentMethods:
        @abstract_method(optional=True)
        def differential_on_basis(self, t):
            """
            The differential of the algebra on the basis (optional).

            INPUT:

            - ``t`` -- the indices of an element of the basis of ``self``

            Returns the differential of the corresponding basis element.
            If implemented, the differential of the algebra is defined
            from it by linearity.

            EXAMPLES::

                sage: from sage.algebras.differential_graded_algebras.commutative_dga import CommutativeDGA
                sage: D = CommutativeDGA(degrees=(1,2), differential=(None, (((1,1), 1),)))
                sage: D.gens()
                (a_1, b_2)
                sage: D.differential_on_basis((0,1))
                a_1 * b_2
            """
            pass

        @lazy_attribute
        def differential(self):
            """
            The differential of this algebra.

            If :meth:`.differential_basis` is available, this
            constructs the differential morphism from ``self``
            to ``self`` by extending it by linearity. Otherwise,
            :meth:`self.differential_by_coercion` is used, if
            available.

            EXAMPLES::

                sage: from sage.algebras.differential_graded_algebras.commutative_dga import CommutativeDGA
                sage: D = CommutativeDGA(degrees=(1,2), differential=(None, (((1,1), 1),)))
                sage: D.gens()
                (a_1, b_2)
                sage: D.differential(D.gen(1))
                a_1 * b_2
            """
            if self.differential_on_basis is not NotImplemented:
                return self._module_morphism(self.differential_on_basis, codomain = self)
            elif hasattr(self, "differential_by_coercion"):
                return self.differential_by_coercion

    class ElementMethods:
        pass

