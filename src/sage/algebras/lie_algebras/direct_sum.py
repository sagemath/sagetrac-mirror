"""
Direct Sum of Lie Algebras

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element_wrapper import ElementWrapper
from sage.categories.category import Category
from sage.algebras.lie_algebras.lie_algebra import LieAlgebra

class DirectSumOfLieAlgebras(LieAlgebra):
    r"""
    The direct sum of lie algebras.

    Given Lie algebras `\mathfrak{g}_1, \mathfrak{g}_2, \ldots \mathfrak{g}_m`,
    the direct sum `\mathfrak{d} = \bigoplus_i \mathfrak{g}_i` is given a Lie
    algebra structure by

    .. MATH::

        \bigl[(x_1, x_2, \ldots, x_n), (y_1, y_2, \ldots, y_n)]
        = ([x_1, y_1], [x_2, y_2], \ldots, [x_n, y_n]).
    """
    @staticmethod
    def __classcall_private__(cls, algs, base_ring=None):
        """
        Normalize arguments and flatten direct sums.

        EXAMPLES::
        """
        if base_ring is None:
            base_ring = g[0].base_ring()
        final = []
        for g in algs:
            # FIXME: Should we instead try to create a common base ring?
            if g.base_ring() != base_ring:
                raise ValueError("inconsistant base rings")
            if isinstance(g, DirectSumOfLieAlgebras):
                final += list(g._g)
            else:
                final.append(g)
        category = Category.meet([g.category() for g in algs])
        return super(DirectSumOfLieAlgebras, cls).__classcall__(cls, tuple(final), base_ring, category)

    def __init__(self, algs, R, category):
        """
        Initialize ``self``.

        EXAMPLES::
        """
        self._g = algs
        LieAlgebra.__init__(self, R, category=category)

    class Element(ElementWrapper):
        def _add_(self, y):
            """
            Add ``self`` and ``y``.

            EXAMPLES::
            """
            return self.__class__(self.parent(), tuple(x + y.value[i]
                                                       for i,x in enumerate(self.value)))

        def _sub_(self, rhs):
            """
            Subtract ``self`` and ``rhs``.

            EXAMPLES::
            """
            return self.__class__(self.parent(), tuple(x - y.value[i]
                                                       for i,x in enumerate(self.value)))

        def __neg__(self):
            """
            Return the negation of ``self``.

            EXAMPLES::
            """
            return self.__class__(self.parent(), tuple(-x for x in self.value))

        def _bracket_(self, y):
            """
            Return the Lie bracket ``[self, y]``.

            EXAMPLES::
            """
            return self.__class__(self.parent(), tuple(x._bracket_(y.value[i])
                                                       for i,x in enumerate(self.value)))

        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of a scalar on ``self``.

            EXAMPLES::
            """
            # With the current design, the coercion model does not have
            # enough information to detect apriori that this method only
            # accepts scalars; so it tries on some elements(), and we need
            # to make sure to report an error.
            if hasattr( scalar, 'parent' ) and scalar.parent() != self.base_ring():
                # Temporary needed by coercion (see Polynomial/FractionField tests).
                if self.base_ring().has_coerce_map_from(scalar.parent()):
                    scalar = self.base_ring()( scalar )
                else:
                    return None
            if self_on_left:
                return self.__class__(self.parent(), tuple(x * scalar for x in self.value))
            return self.__class__(self.parent(), tuple(scalar * x for x in self.value))

        def __getitem__(self, i):
            """
            Return the ``i``-th component of ``self``.

            EXAMPLES::
            """
            return self.value.__getitem__(i)

