# -*- coding: utf-8 -*-
"""
Implement the cartesian product of combinatorial classes of structures.

Let `F` and `G` be combinatorial structures classes.

We denote `H = F \cdot G` (or just `H = FG`) the cartesian product of `F` and `G`.

The generating series of `FG` is given the *Hadamard product* of the series:

MATH::

    (fg)(x) := f(x) \square g(x) = \sum_{n \geqslant 0} f_n g_n x^n

AUTHOR:

- Jean-Baptiste Priez (2014)
"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from itertools import product, imap
from sage.combinat.structures import Structures
from __init__ import _Operations
from sage.combinat.structures.operations.product import CauchyProduct


class CartesianProduct(_Operations):

    def _repr_(self):
        """
        TESTS::

        """
        return "Cartesian product of structures : '" + ", ".join(map(repr, self._structures)) + "'"

    def generating_series(self):
        """
        Return the generating serie of ``self``. This is the Hadamard product of
        the series ``F`` and ``G``.

        MATH::

            (fg)(x) := f(x) \square g(x) = \sum_{n \geqslant 0} f_n g_n x^n

        """
        return # TODO


    class GradedComponent(Structures.GradedComponent):

        def __iter__(self):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: C = Compositions()
                sage: BcC = B.cartesian_product(C)
                sage: BcC.graded_component(3).list()
                [[[., [., [., .]]], [1, 1, 1]],
                 [[., [., [., .]]], [1, 2]],
                 [[., [., [., .]]], [2, 1]],
                 [[., [., [., .]]], [3]],
                 [[., [[., .], .]], [1, 1, 1]],
                 [[., [[., .], .]], [1, 2]],
                 [[., [[., .], .]], [2, 1]],
                 [[., [[., .], .]], [3]],
                 [[[., .], [., .]], [1, 1, 1]],
                 [[[., .], [., .]], [1, 2]],
                 [[[., .], [., .]], [2, 1]],
                 [[[., .], [., .]], [3]],
                 [[[., [., .]], .], [1, 1, 1]],
                 [[[., [., .]], .], [1, 2]],
                 [[[., [., .]], .], [2, 1]],
                 [[[., [., .]], .], [3]],
                 [[[[., .], .], .], [1, 1, 1]],
                 [[[[., .], .], .], [1, 2]],
                 [[[[., .], .], .], [2, 1]],
                 [[[[., .], .], .], [3]]]

            """
            k = self.grading()
            return imap(
                self._element_constructor_,
                product(*imap(
                    lambda F: F.graded_component(k),
                    self.ambient()._structures
            )  ))

    class Element(CauchyProduct.Element):
        pass
