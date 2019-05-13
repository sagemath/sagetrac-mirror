r"""
Filtered Modules

A *filtered module* over a ring `R` with a totally ordered
indexing set `I` (typically `I = \NN`) is an `R`-module `M` equipped
with a family `(F_i)_{i \in I}` of `R`-submodules satisfying
`F_i \subseteq F_j` for all `i,j \in I` having `i \leq j`, and
`M = \bigcup_{i \in I} F_i`. This family is called a *filtration*
of the given module `M`.

.. TODO::

    Implement a notion for decreasing filtrations: where `F_j \subseteq F_i`
    when `i \leq j`.

.. TODO::

    Implement filtrations for all concrete categories.

.. TODO::

    Implement `\operatorname{gr}` as a functor.
"""
# ****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory


class FilteredModulesCategory(RegressiveCovariantConstructionCategory, Category_over_base_ring):
    def __init__(self, base_category):
        """
        EXAMPLES::

            sage: C = Algebras(QQ).Filtered()
            sage: C
            Category of filtered algebras over Rational Field
            sage: C.base_category()
            Category of algebras over Rational Field
            sage: sorted(C.super_categories(), key=str)
            [Category of algebras over Rational Field,
             Category of filtered modules over Rational Field]

            sage: AlgebrasWithBasis(QQ).Filtered().base_ring()
            Rational Field
            sage: HopfAlgebrasWithBasis(QQ).Filtered().base_ring()
            Rational Field
        """
        super(FilteredModulesCategory, self).__init__(base_category, base_category.base_ring())

    _functor_category = "Filtered"

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: AlgebrasWithBasis(QQ).Filtered()  # indirect doctest
            Category of filtered algebras with basis over Rational Field
        """
        return "filtered {}".format(self.base_category()._repr_object_names())

class FilteredModules(FilteredModulesCategory):
    r"""
    The category of filtered modules over a given ring `R`.

    A *filtered module* over a ring `R` with a totally ordered
    indexing set `I` (typically `I = \NN`) is an `R`-module `M` equipped
    with a family `(F_i)_{i \in I}` of `R`-submodules satisfying
    `F_i \subseteq F_j` for all `i,j \in I` having `i \leq j`, and
    `M = \bigcup_{i \in I} F_i`. This family is called a *filtration*
    of the given module `M`.

    EXAMPLES::

        sage: Modules(ZZ).Filtered()
        Category of filtered modules over Integer Ring
        sage: Modules(ZZ).Filtered().super_categories()
        [Category of modules over Integer Ring]

    TESTS::

        sage: TestSuite(Modules(ZZ).Filtered()).run()

    REFERENCES:

    - :wikipedia:`Filtration_(mathematics)`
    """
    def extra_super_categories(self):
        r"""
        Add :class:`VectorSpaces` to the super categories of ``self`` if
        the base ring is a field.

        EXAMPLES::

            sage: Modules(QQ).Filtered().extra_super_categories()
            [Category of vector spaces over Rational Field]
            sage: Modules(ZZ).Filtered().extra_super_categories()
            []

        This makes sure that ``Modules(QQ).Filtered()`` returns an
        instance of :class:`FilteredModules` and not a join category of
        an instance of this class and of ``VectorSpaces(QQ)``::

            sage: type(Modules(QQ).Filtered())
            <class 'sage.categories.filtered_modules.FilteredModules_with_category'>

        .. TODO::

            Get rid of this workaround once there is a more systematic
            approach for the alias ``Modules(QQ)`` -> ``VectorSpaces(QQ)``.
            Probably the latter should be a category with axiom, and
            covariant constructions should play well with axioms.
        """
        from sage.categories.modules import Modules
        from sage.categories.fields import Fields
        base_ring = self.base_ring()
        if base_ring in Fields:
            return [Modules(base_ring)]
        else:
            return []

    class SubcategoryMethods:

        @cached_method
        def Connected(self):
            r"""
            Return the full subcategory of the connected objects of ``self``.

            A filtered `R`-module `M` with filtration
            `(F_0, F_1, F_2, \ldots)` (indexed by `\NN`)
            is said to be *connected* if `F_0` is isomorphic
            to `R`.

            EXAMPLES::

                sage: Modules(ZZ).Filtered().Connected()
                Category of filtered connected modules over Integer Ring
                sage: Coalgebras(QQ).Filtered().Connected()
                Join of Category of filtered connected modules over Rational Field
                    and Category of coalgebras over Rational Field
                sage: AlgebrasWithBasis(QQ).Filtered().Connected()
                Category of filtered connected algebras with basis over Rational Field

            TESTS::

                sage: TestSuite(Modules(ZZ).Filtered().Connected()).run()
                sage: Coalgebras(QQ).Filtered().Connected.__module__
                'sage.categories.filtered_modules'
            """
            return self._with_axiom("Connected")

    class Connected(CategoryWithAxiom_over_base_ring):
        pass

