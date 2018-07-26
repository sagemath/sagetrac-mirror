"""
This patch implement some tools for testing primitivity.
"""
from sage.combinat.cha.patches.monkey_patching import MonkeyPatch

from sage.categories.graded_hopf_algebras_with_basis import \
    GradedHopfAlgebrasWithBasis


class _(MonkeyPatch, GradedHopfAlgebrasWithBasis.ParentMethods):

    def get_primitive_elements(self, max_elem=100):
        """
        Return the list of the primitive elements of the basis between the
        first ``max_elem`` explored on the basis.

        TESTS::

            sage: F = FQSym(QQ).F()
            sage: F.get_primitive_elements()
        """
        it = self.basis().keys().__iter__()
        return filter(
            lambda elem: self(elem).is_primitive(),
            map(lambda _: it.next(), range(max_elem))
        )

    def get_group_type_elements(self, max_elem=100):
        """
        Return the list of the group type elements of the basis between the
        first ``max_elem`` explored on the basis.

        TESTS::

            sage: F = FQSym(QQ).F()
            sage: F.get_group_type_elements()
            [[]]
        """
        it = self.basis().keys().__iter__()
        return filter(
            lambda elem: self(elem).is_group_type(),
            map(lambda _: it.next(), range(max_elem))
        )


class _(MonkeyPatch, GradedHopfAlgebrasWithBasis.ElementMethods):

    def is_primitive(self):
        """
        Indicate if the element *self* is primitive.

        TESTS::

            sage: F = FQSym(QQ).F()
            sage: F[1].is_primitive()
            True
            sage: F[1].coproduct()
            F[] # F[1] + F[1] # F[]
            sage: F[1,2].is_primitive()
            False
            sage: F[1,2].coproduct()
            F[] # F[1, 2] + F[1] # F[1] + F[1, 2] # F[]
        """
        from sage.categories.tensor import tensor
        return self.coproduct() == tensor((self, self.parent().one())) + \
            tensor((self.parent().one(), self))

    def is_group_type(self):
        """
        Indicate if the element *self* is group type.

            ...
        """
        from sage.categories.tensor import tensor
        return self.coproduct() == tensor((self, self))
