r"""
Lie Algebras With Basis

AUTHORS:

- Travis Scrimshaw (07-15-2013): Initial implementation
"""
#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
#from sage.misc.cachefunc import cached_method
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.lie_algebras import LieAlgebras

class LieAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    Category of Lie algebras with a basis.
    """
    _base_category_class_and_axiom = [LieAlgebras, "WithBasis"]

    def example(self, gens=None):
        """
        Return an example of a Lie algebra as per
        :meth:`Category.example <sage.categories.category.Category.example>`.

        EXAMPLES::

            sage: LieAlgebras(QQ).WithBasis().example()
            An example of a Lie algebra: the abelian Lie algebra on the
             generators indexed by Partitions over Rational Field

        Another set of generators can be specified as an optional argument::

            sage: LieAlgebras(QQ).WithBasis().example(Compositions())
            An example of a Lie algebra: the abelian Lie algebra on the
             generators indexed by Compositions of non-negative integers
             over Rational Field
        """
        if gens is None:
            from sage.combinat.partition import Partitions
            gens = Partitions()
        from sage.categories.examples.lie_algebras_with_basis import Example
        return Example(self.base_ring(), gens)

    class ParentMethods:
        def _basis_cmp(self, x, y):
            """
            Compare two basis element indices. The default is to call ``cmp``.

            TESTS::

                sage: L = LieAlgebras(QQ).WithBasis().example()
                sage: L._basis_cmp(Partition([3,1]), Partition([2,1,1]))
                1
            """
            return cmp(x, y)

        @abstract_method(optional=True)
        def bracket_on_basis(self, x, y):
            """
            Return the bracket of basis elements indexed by ``x`` and ``y``
            where ``x < y``. If this is not implemented, then the method
            ``_bracket_()`` for the elements must be overwritten.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).WithBasis().example()
                sage: L.bracket_on_basis(Partition([3,1]), Partition([2,2,1,1]))
                0
            """

        def module(self):
            """
            Return an `R`-module which is isomorphic to the
            underlying `R`-module of ``self``.

            See
            :meth:`sage.categories.lie_algebras.LieAlgebras.module` for
            an explanation.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).WithBasis().example()
                sage: L.module()
                Free module generated by Partitions over Rational Field
            """
            from sage.combinat.free_module import CombinatorialFreeModule
            try:
                # Try to see if it has an indexing set
                return CombinatorialFreeModule(self.base_ring(), self.basis().keys())
            except AttributeError:
                # Otherwise just index by the basis of ``self`` as a fallback
                return CombinatorialFreeModule(self.base_ring(), self.basis())

        def from_vector(self, v):
            """
            Return the element of ``self`` corresponding to the
            vector ``v`` in ``self.module()``.

            Implement this if you implement :meth:`module`; see the
            documentation of
            :meth:`sage.categories.lie_algebras.LieAlgebras.module`
            for how this is to be done.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: u = L.from_vector(vector(QQ, (1, 0, 0))); u
                (1, 0, 0)
                sage: parent(u) is L
                True
            """
            B = self.basis()
            return self.sum(v[i] * B[i] for i in v.support())

    class ElementMethods:
        def _bracket_(self, y):
            """
            Return the Lie bracket ``[self, y]``, where ``y`` is an
            element of the same Lie algebra as ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).WithBasis().example()
                sage: G = L.lie_algebra_generators()
                sage: x = G[Partition([4,3,3,1])]
                sage: y = G[Partition([6,1])]
                sage: x.bracket(y)
                0
            """
            P = self.parent()
            def term(ml,mr):
                comp = P._basis_cmp(ml,mr)
                if comp == 0:
                    return P.zero()
                if comp < 0:
                    return P.bracket_on_basis(ml, mr)
                return -P.bracket_on_basis(mr, ml)
            return P.sum(cl*cr * term(ml,mr) for ml,cl in self for mr,cr in y)

        def to_vector(self):
            """
            Return the vector in ``g.module()`` corresponding to the
            element ``self`` of ``g`` (where ``g`` is the parent of
            ``self``).

            Implement this if you implement ``g.module()``.
            See :meth:`sage.categories.lie_algebras.LieAlgebras.module`
            for how this is to be done.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.an_element().to_vector()
                (0, 0, 0)

            .. TODO::

                Doctest this implementation on an example not overshadowed.
            """
            M = self.parent().module()
            B = M.basis()
            return M.sum(self[i] * B[i] for i in self.support())

