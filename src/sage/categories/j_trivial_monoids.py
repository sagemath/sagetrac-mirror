r"""
Finite J-Trivial Monoids
"""
#*****************************************************************************
#  Copyright (C) 2009-2010 Florent Hivert <florent.hivert at univ-rouen.fr>
#                2009-2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.l_trivial_monoids import LTrivialMonoids
from sage.categories.r_trivial_monoids import RTrivialMonoids
from sage.misc.cachefunc import cached_in_parent_method
from sage.sets.family import Family


class JTrivialMonoids(Category):
    r"""
    The category of `J`-trivial monoids

    Let `M` be a monoid. The `J`-relation on `M` is given by
    `a \sim b iff MaM = MbM`. A monoid is `J`-trivial if all its `J`-classes
    are of cardinality one.
    """

    @cached_method
    def super_categories(self):
        return [LTrivialMonoids(), RTrivialMonoids()]


    class Finite(CategoryWithAxiom):
        """
        The category of finite `J`-trivial monoids
        """
        class ParentMethods:

            @cached_method
            def semigroup_generators(self):
                """
                Returns the canonical minimal set of generators. It
                consists of the irreducible elements, that is elements
                which are not of the form `x*y` with `x` and `y` in
                ``self`` distinct from `x*y`.
                """
                res = []
                G = self.cayley_graph(side="twosided")
                for x in G:
                    if x == self.one().value:
                        continue
                    incoming = set(G.incoming_edges(x))
                    if all(l == x for u,v,l in incoming):
                        res.append(self(x))
                return Family(res)


        class ElementMethods:

            @cached_in_parent_method
            def symbol(self, side = "left"):
                """
                INPUT:

                 - ``self`` -- a monoid element `x`
                 - ``side`` -- "left", "right"

                Returns the unique minimal idempotent `e` (in J-order)
                such that `e x = x` (resp. `xe = x`).
                """
                monoid = self.parent()

                if side == "left":
                    fix = [ s for s in monoid.semigroup_generators() if (s * self == self) ]
                else:
                    fix = [ s for s in monoid.semigroup_generators() if (self * s == self) ]
                return (monoid.prod(fix)).pow_omega()
