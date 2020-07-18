r"""
Normed Vector Spaces
"""
#*****************************************************************************
#  Copyright (C) 2014, 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#                2020 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory

class NormedSpacesCategory(RegressiveCovariantConstructionCategory, Category_over_base_ring):

    _functor_category = "Normed"

    def __init__(self, base_category):
        """
        EXAMPLES::

            sage: VectorSpaces(QQ).Filtered().Normed().super_categories()


            sage: FreeModules(ZZ).Normed().super_categories()

            sage: VectorSpaces(QQ).Filtered().Normed().super_categories()
            [Category of , Category of metric spaces]

        """
        super(NormedSpacesCategory, self).__init__(base_category, base_category.base_ring())

    @classmethod
    def default_super_categories(cls, category):
        """
        Return the default super categories of ``category.Normed()``.

        Mathematical meaning: if `A` is a normed space in the
        category `C`, then `A` is also a metric space.

        INPUT:

        - ``cls`` -- a subclass of ``NormedSpacesCategory``
        - ``category`` -- a category `C`

        OUTPUT:

        A (join) category

        In practice, this returns ``category.Metric()``, joined
        together with the result of the method
        :meth:`RegressiveCovariantConstructionCategory.default_super_categories()
        <sage.categories.covariant_functorial_construction.RegressiveCovariantConstructionCategory.default_super_categories>`
        (that is the join of ``category`` and ``cat.Normed()`` for
        each ``cat`` in the super categories of ``category``).

        EXAMPLES:

        Consider ``category=Groups()``. Then, a group `G` with a metric
        is simultaneously a topological group by itself, and a
        metric space::

        This resulted from the following call::

            sage: sage.categories.normed_vector_spaces.NormedSpacesCategory.default_super_categories(VectorSpaces(QQ).Filtered())
            Join of Category of topological groups and Category of metric spaces
        """
        return Category.join([category.Metric(),
                              super(NormedSpacesCategory, cls).default_super_categories(category)])

    # We currently don't have a use for this, but we probably will
    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: VectorSpaces(QQ).Normed()  # indirect doctest
            Join of Category of topological groups and Category of metric spaces
        """
        return "normed {}".format(self.base_category()._repr_object_names())

class NormedVectorSpaces(NormedSpacesCategory):
    r"""
    The category of normed vector spaces.

    A norm on a vector space `V` is a function `\|\cdot\| : V \to \RR` such that:

    - `\|x\| \geq 0`,
    - `\|x\| = 0` if and only if `x = 0`,
    - `\|x+y\| \leq \|x\| + \|y\|`,
    - `\|\lambda x\| = |\lambda| \cdot \|x\|`.

    """

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Sets().Metric()  # indirect doctest
            Category of metric spaces
        """
        return "normed vector spaces"

    class ElementMethods:

        def abs_squared(self):
            """
            Return the square of the norm of ``self``.

            Classes implementing a norm should implement one of the methods
            ``abs_squared`` and ``abs``.  The other method is supplied
            automatically.

            """
            return self.norm() ** 2

        def abs(self):
            """
            Return the norm of ``self``.

            Classes implementing a norm should implement one of the methods
            ``abs_squared`` and ``abs``.  The other method is supplied
            automatically.
            """
            return sqrt(self.abs_squared())
