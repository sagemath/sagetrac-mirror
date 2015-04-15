r"""
Finite semigroups
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Florent Hivert <florent.hivert at univ-rouen.fr>
#                2008-2014 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.misc import attrcall
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

class FiniteSemigroups(CategoryWithAxiom):
    r"""
    The category of finite (multiplicative) semigroups.

    A finite semigroup is a :class:`finite set <FiniteSets>` endowed
    with an associative binary operation `*`.

    .. WARNING::

        Finite semigroups in Sage used to be automatically endowed
        with an :class:`enumerated set <EnumeratedSets>` structure;
        the default enumeration is then obtained by iteratively
        multiplying the semigroup generators. This forced any finite
        semigroup to either implement an enumeration, or provide
        semigroup generators; this was often inconvenient.

        Instead, finite semigroups that provide a distinguished finite
        set of generators with :meth:`semigroup_generators` should now
        explicitly declare themselves in the category of
        :class:`finitely generated semigroups
        <Semigroups.FinitelyGeneratedSemigroup>`::

            sage: Semigroups().FinitelyGenerated()
            Category of finitely generated semigroups

        This is a backward incompatible change.

    EXAMPLES::

        sage: C = FiniteSemigroups(); C
        Category of finite semigroups
        sage: C.super_categories()
        [Category of semigroups, Category of finite sets]
        sage: sorted(C.axioms())
        ['Associative', 'Finite']
        sage: C.example()
        An example of a finite semigroup: the left regular band generated by ('a', 'b', 'c', 'd')

    TESTS::

        sage: TestSuite(C).run()
    """

    class ParentMethods:
        def idempotents(self):
            r"""
            Returns the idempotents of the semigroup

            EXAMPLES::

                sage: S = FiniteSemigroups().example(alphabet=('x','y'))
                sage: sorted(S.idempotents())
                ['x', 'xy', 'y', 'yx']
            """
            return [x for x in self if x.is_idempotent()]

        @cached_method
        def j_classes(self):
            r"""
            Returns the $J$-classes of the semigroup.

            Two elements $u$ and $v$ of a monoid are in the same $J$-class
            if $u$ divides $v$ and $v$ divides $u$.

            OUTPUT:

             All the $J$-classes of self, as a list of lists.

            EXAMPLES::

                sage: S = FiniteSemigroups().example(alphabet=('a','b', 'c'))
                sage: sorted(map(sorted, S.j_classes()))
                [['a'], ['ab', 'ba'], ['abc', 'acb', 'bac', 'bca', 'cab', 'cba'], ['ac', 'ca'], ['b'], ['bc', 'cb'], ['c']]
            """
            return self.cayley_graph(side="twosided", simple=True).strongly_connected_components()

        @cached_method
        def j_classes_of_idempotents(self):
            r"""
            Returns all the idempotents of self, grouped by J-class.

            OUTPUT:

             a list of lists.

            EXAMPLES::

                sage: S = FiniteSemigroups().example(alphabet=('a','b', 'c'))
                sage: sorted(map(sorted, S.j_classes_of_idempotents()))
                [['a'], ['ab', 'ba'], ['abc', 'acb', 'bac', 'bca', 'cab', 'cba'], ['ac', 'ca'], ['b'], ['bc', 'cb'], ['c']]
            """
            return [l for l in map(lambda cl: list(filter(attrcall('is_idempotent'), cl)),
                                   self.j_classes()) if len(l) > 0]

        @cached_method
        def j_transversal_of_idempotents(self):
            r"""
            Returns a list of one idempotent per regular J-class

            EXAMPLES::

                sage: S = FiniteSemigroups().example(alphabet=('a','b', 'c'))
                sage: sorted(S.j_transversal_of_idempotents())
                ['a', 'ab', 'ac', 'acb', 'b', 'c', 'cb']
            """
            def first_idempotent(l):
                for x in l:
                    if x.is_idempotent():
                        return x
                return None
            return [x for x in map(first_idempotent, self.j_classes()) if not x is None]

        # TODO: compute eJe, where J is the J-class of e
        # TODO: construct the action of self on it, as a permutation group
