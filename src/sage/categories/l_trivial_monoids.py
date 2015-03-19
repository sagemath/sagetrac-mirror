from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.h_trivial_monoids import HTrivialMonoids

class LTrivialMonoids(Category):
    """
    The category of `L`-trivial monoids

    Let `M` be a monoid. The *`L`-preorder* is defined by `x\leq_L y`
    if `x \in My`.  The *`L`-classes* are the the equivalence classes
    for the associated equivalence relation.  A monoid is *`L`-trivial*
    if all its `L`-classes are trivial, that is of cardinality `1`, or
    equivalently if the `L`-preoder is in fact an order.

    EXAMPLES::

        sage: C = LTrivialMonoids(); C
        Category of l trivial monoids
        sage: C.super_categories()
        [Category of h trivial monoids]

    .. seealso:: :class:`RTrivialMonoids`, :class:`HTrivialMonoids`, :class:`JTrivialMonoids`
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES:

        An L-trivial monoid is also H-trivial::

            sage: LTrivialMonoids().super_categories()
            [Category of h trivial monoids]
        """
        return [HTrivialMonoids()]

    class Finite(CategoryWithAxiom):

        class ParentMethods:

            def index_of_regular_j_class(self, idempotent):
                """
                Returns the index that should be used for an idempotent in the transversal

                In this implementation, each idempotent e is indexed
                by the subset of the indices `i` of the generators
                `s_i` such that `es_i=e` (that is `s_1` acts by `1` on
                the corresponding simple module).

                .. seealso:: :meth:`FiniteSemigroups.ParentMethods.j_transversal_of_idempotents`

                .. todo::

                    This is mostly a duplicate of
                    :meth:`RTrivialMonoids.Finite.ParentMethods.j_transversal_of_idempotents`

                    Instead this should be generalized to
                    DASemigroups.Finite, by testing if idempotent *
                    s[i] is in the same J-class. And recycled to build
                    the corresponding simple module.

                EXAMPLES::

                    TODO!
                """
                s = self.semigroup_generators()
                return tuple(i for i in s.keys() if s[i] * idempotent == idempotent)
