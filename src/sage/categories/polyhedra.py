r"""
Category of sets of polyhedra
"""
#*****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_types import Category_over_base_ring

class PolyhedralSets(Category_over_base_ring):
    r"""
    The category of sets of polyhedra over a ring.

    A set of polyhedra is a parent.

    A polyhedron is an element.

    EXAMPLES:

    We create the category of sets of polyhedra over `\QQ`::

        sage: PolyhedralSets(QQ)
        Category of polyhedral sets over Rational Field

    TESTS::

        sage: TestSuite(PolyhedralSets(RDF)).run()

        sage: P = Polyhedron()
        sage: P.parent().category().element_class
        <class 'sage.categories.category.JoinCategory.element_class'>
        sage: P.parent().category().element_class.mro()
        [<class 'sage.categories.category.JoinCategory.element_class'>,
         <class 'sage.categories.polyhedra.PolyhedralSets.element_class'>,
         <class 'sage.categories.magmas.Magmas.Commutative.element_class'>,
         <class 'sage.categories.magmas.Magmas.element_class'>,
         <class 'sage.categories.additive_monoids.AdditiveMonoids.element_class'>,
         <class 'sage.categories.additive_magmas.AdditiveMagmas.AdditiveUnital.element_class'>,
         <class 'sage.categories.additive_semigroups.AdditiveSemigroups.element_class'>,
         <class 'sage.categories.additive_magmas.AdditiveMagmas.element_class'>,
         <class 'sage.categories.finite_enumerated_sets.FiniteEnumeratedSets.element_class'>,
         <class 'sage.categories.enumerated_sets.EnumeratedSets.element_class'>,
         <class 'sage.categories.combinatorial_polyhedral_sets.CombinatorialPolyhedralSets.Finite.element_class'>,
         <class 'sage.categories.finite_sets.FiniteSets.element_class'>,
         <class 'sage.categories.combinatorial_polyhedral_sets.CombinatorialPolyhedralSets.element_class'>,
         <class 'sage.categories.sets_cat.Sets.element_class'>,
         <class 'sage.categories.sets_with_partial_maps.SetsWithPartialMaps.element_class'>,
         <class 'sage.categories.objects.Objects.element_class'>,
         <class 'object'>]
        sage: isinstance(P, P.parent().category().element_class)
        True
    """

    def __init__(self, R):
        """
        TESTS::

            sage: PolyhedralSets(AA)
            Category of polyhedral sets over Algebraic Real Field
        """
        Category_over_base_ring.__init__(self, R)

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: PolyhedralSets(QQ).super_categories()
            [Category of commutative magmas,
             Category of additive monoids,
             Category of combinatorial polyhedral sets]
        """
        from sage.categories.magmas import Magmas
        from sage.categories.additive_monoids import AdditiveMonoids
        from sage.categories.combinatorial_polyhedral_sets import CombinatorialPolyhedralSets
        return [Magmas().Commutative(), AdditiveMonoids(), CombinatorialPolyhedralSets()]




