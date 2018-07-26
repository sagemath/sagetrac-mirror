# -*- coding: utf-8 -*-
r"""
Dendriform algebras

AUTHORS:

 - Rémi Maurice & Jean-Baptiste Priez (first version)
"""
#*****************************************************************************
#  Copyright (C) 2013      Rémi Maurice <maurice@univ-mlv.fr>
#                          Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.category_types import Category_over_base_ring

from dendriform_algebras import DendriformAlgebras
from dendriform_coalgebras import DendriformCoalgebras


class BidendriformBialgebras(Category_over_base_ring):
    """
    The category of bidendriform bialgebras over a given base ring.

    A bidendriform bialgebra over a ring `R` is an bialgebra which is a
    dendriform algebra and a dendriform coalgebra such that respect these
    following conditions:

    .. MATH::

        TODO::

    TODO: should `R` be a commutative ring?

    EXAMPLES::

        sage: DendriformAlgebras(ZZ)
        Category of dendriform algebras over Integer Ring
        sage: Algebras(ZZ).super_categories()
        [Category of non unital algebras over Integer Ring, Category of rings]

    TESTS::

        sage: TestSuite(Algebras(ZZ)).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: Algebras(ZZ).super_categories()
            [Category of non unital algebras over Integer Ring,
             Category of rings]
        """
        R = self.base_ring()
        return [DendriformAlgebras(R), DendriformCoalgebras(R)]

    class WithBasis(Category_over_base_ring):

        def __repr__(self):
            return "Category of bidendriform bialgebras with basis over %s" % \
                (self.base_ring())

        def super_categories(self):
            """
            EXAMPLES::

                sage: Algebras(ZZ).super_categories()
                [Category of non unital algebras over Integer Ring,
                 Category of rings]
            """
            R = self.base_ring()
            return [DendriformAlgebras.WithBasis(R),
                    DendriformCoalgebras.WithBasis(R)]
