r"""
Ordered Sets
"""

#*****************************************************************************
#  Copyright (C) 2020 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.lattice_posets import LatticePosets

class OrderedSets(Category):

    r"""
    The category of (totally) ordered sets.

    EXAMPLES::

        sage: from sage.categories.ordered_sets import OrderedSets
        sage: OrderedSets()
        Category of ordered sets
        sage: OrderedSets().super_categories()
        [Category of lattice posets]

    """

    @cached_method
    def super_categories(self):
        r"""
        Return a list of the (immediate) super categories of
        ``self``, as per :meth:`Category.super_categories`.

        EXAMPLES::

            sage: from sage.categories.ordered_sets import OrderedSets
            sage: OrderedSets().super_categories()
            [Category of lattice posets]
        """
        return [LatticePosets()]

    class ParentMethods:

        def meet(self, x, y):
            """
            Returns the meet of ``x`` and ``y`` in this lattice

            Because ``self`` is totally ordered, this is just the maximum
            of ``x`` and ``y``.

            INPUT:

             - ``x``, ``y`` -- elements of ``self``

            """
            if self.le(x, y):
                return y
            else:
                return x

        def join(self, x, y):
            """
            Returns the join of `x` and `y` in this lattice

            Because ``self`` is totally ordered, this is just the minimum
            of ``x`` and ``y``.

            INPUT:

             - ``x``, ``y`` -- elements of ``self``

            """
            if self.le(x, y):
                return x
            else:
                return y
