# -*- coding: utf-8 -*-
"""
Category of species of structures
---------------------------------


Definition (**Species** _[BLL])
-------------------------------

A *species of structures* is a rule `F` which produces

    i) for each finite set `U`, a finite set `F[U]`,
   ii) for each bijection `\sigma: U \to V`, a function `F[\sigma]: F[U] \to
       F[V]`.

The functions `F[\sigma]` should further satisfy the following functorial
properties:

    a) for all bijections `\sigma : U \to V` and `tau :V \to W`,

    MATH::

        F[\tau \circ \sigma] = F[\tau] \circ F[\sigma]\,,

    b) for the identity map `Id_U : U \to U`,

    MATH::

        F[Id_U] = Id_F[U].

An element `s \in F[U]` is called an *`F`-structure on `U`* (or even a
*structure of species `F` on `U`*).
The function `F [\sigma]` is called the *transport of `F`-structures along
`\sigma`*.

References:
-----------

.. [BLL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle, and Pierre Leroux

"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.partition import Partition
from sage.misc.misc_c import prod
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.combinat.sf.sf import SymmetricFunctions
from sage.misc.cachefunc import cached_method
from sage.rings.arith import factorial
from sage.sets.set import Set
from sage.categories.category import Category
from sage.categories.sets_cat import Sets
from sage.misc.abstract_method import abstract_method
from sage.categories.combinatorial_structures import CombinatorialStructures
from sage.rings.all import Integers


class Species(Category):
    """
    Category of species of structures (see _[BLL])
    """

    def super_categories(self):
        return [CombinatorialStructures()]

    class ParentMethods:

        @abstract_method(optional=False)
        def transport(self, sigma):
            """
            @param sigma: a permutation or a bijection of set `\sigma : U \to V`
            @return: `F[\sigma]`, the *transport* of `F`-structures along
            `\sigma`.

            `F[\sigma]` has to be a (sage) morphism from `F[U]` to `F[V]` such
            that can be compose automatically.
            """


        def structures(self, U):
            """
            @param U: a finite set or a non-negative integer
            @return: `F[U]`, the set of `F`-structures on `U` (or on `[U]` if
                     `U` is an integer).
            """
            if U in Sets():
                return self.Structures(U)
            else:
                assert(U in Integers()), "*U* must be a finite set or an integer"
                return self.Structures(set(range(1, U+1)))

        def Fix(self, sigma):
            """
            TESTS::

                sage: from sage.categories.examples.species_permutations import \
                    Permutations
                sage: P = Permutations()
                sage: from sage.combinat.permutation import Permutations as Perm
                sage: for sig in Perm(3):
                    P.Fix(sig)
                {[(1, 3), (2,)], [(1,), (2,), (3,)], [(1, 2), (3,)], [(1, 2, 3)], [(1,), (2, 3)]}
                {[(1, 2, 3)], [(1,), (2, 3)]}
                {[(1, 2), (3,)], [(1, 2, 3)]}
                {[(1, 2, 3)]}
                {[(1, 2, 3)]}
                {[(1, 3), (2,)], [(1, 2, 3)]}
            """
            from sage.combinat.permutation import Permutations
            assert(sigma in Permutations())
            l = []
            trans = self.transport(sigma)
            for obj in self.graded_component(len(sigma)):
                if trans(obj) == obj:
                    l.append(obj)
            return Set(l)

        @cached_method
        def fix(self, sigma):
            """
            TESTS::

                sage: from sage.categories.examples.species_permutations import \
                    Permutations
                sage: P = Permutations()
                sage: from sage.combinat.permutation import Permutations as Perm
                sage: for sig in Perm(3):
                    P.fix(sig)
                5
                2
                2
                1
                1
                2
            """
            return Integer(len(self.Fix(sigma)))

        def cycle_index_series(self, n):
            """
            The *cycle index series* of a species of structures `F` is the
            formal series in symmetric functions::

            MATH::

                Z_F = \sum_{n \geqslant 0} \frac{1}{n!} \left(\sum_{\sigma \in S_n}
                \mathtt{fix}\ F[\sigma] p_1^{\sigma_1} p_2^{\sigma_2} p_3^{\sigma_3}
                \cdots \right)\,,

            with `S_n` the group of permutations of `[n]` and
            `\mathtt{fix}\ F[\sigma]` the number of `F`-structures on `[n]`
            fixed by `F[\sigma]`.

            """
            # TODO :: use LazyPowerSeries?
            return self.graded_component(n).cycle_index_series()

        #def _test_transport_composition(self):
            # TODO implement somewhere how to test the composition
        #    pass

        #def _test_transport_identity(self):
            # TODO   ---------------------------   the identity
        #    pass


    class ElementMethods:

        def underlying_set(self):
            pass

        def transport(self, sigma):
            """
            TESTS::

                sage: Permutation([1,3,2]).transport([3,1,2])
                [3, 2, 1]

            """
            return self.parent().transport(sigma)(self)

    class Structures(Category):

        def super_categories(self):
            return [CombinatorialStructures.GradedComponents()]

        class ParentMethods:

            def _repr_(self):
                return "Structures of " + repr(self.ambient()) + \
                       " over the set `" + repr(self.underlying_set()) + "`"

            @abstract_method(optional=False)
            def underlying_set(self):
                """

                """

            def cycle_index_series(self):
                """

                """
                from sage.combinat.permutation import Permutations
                p = SymmetricFunctions(QQ).p()
                n = self.grading()
                return Integer(1)/factorial(n) * sum(self.fix(sig) *
                        prod(p.monomial(Partition([i+1])) ** Integer(cy)
                             for i, cy in enumerate(sig.cycle_type().to_exp())
                        ) for sig in Permutations(n))

            def Fix(self, sigma):
                return self.ambient().Fix(sigma)

            def fix(self, sigma):
                return self.ambient().fix(sigma)

            def grading(self):
                return len(self.underlying_set())