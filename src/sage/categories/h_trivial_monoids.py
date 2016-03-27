from sage.misc.cachefunc import cached_method
from sage.categories.monoids import Monoids
from sage.categories.category import Category
from sage.categories.category_with_axiom import CategoryWithAxiom

class HTrivialMonoids(Category):
    """
    The category of `H`-trivial monoids

    Let `M` be a monoid. The *`H`-preorder* is defined by `x\leq_H y`
    if `x \in My` and `x \in yM`. The *`H`-classes* are the the
    equivalence classes for the associated equivalence relation.  A
    monoid is *`H`-trivial* if all its `H`-classes are trivial, that
    is of cardinality `1`, or equivalently if the `H`-preoder is in
    fact an order.

    EXAMPLES::

        sage: from sage.categories.h_trivial_monoids import *
        sage: C = HTrivialMonoids(); C
        Category of h trivial monoids
        sage: C.super_categories()
        [Category of monoids]
        sage: C.example()
        NotImplemented


    .. seealso:: :class:`LTrivialMonoids`, :class:`RTrivialMonoids`, :class:`JTrivialMonoids`
    """

    @cached_method
    def super_categories(self):
        """

        """
        return [Monoids()]

    class Finite(CategoryWithAxiom):

        class ParentMethods:
            pass

        class ElementMethods:

            def pow_omega(self):
                """
                The omega power of ``self``.
                """
                res_old = self
                res_new = res_old*res_old
                while res_old != res_new:
                    res_old = res_new
                    res_new = res_old*res_old
                return res_new

            pow_infinity = pow_omega # for backward compatibility
