# -*- coding: utf-8 -*-
r"""
Formal Power Series

Reference
---------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux
  Cambridge University Press, 1998
}
"""
#*****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.category import Category
from sage.misc.cachefunc import cached_method


class FormalPowerSeries(Category):
    """
    A *formal power series*

    MATH::

        F(t) = \sum_{n \geqslant 0} f_n t^n\,.

    The idea of this category is to formalize common methods use by power series.
    (We define only few thing to use exponential or ordinary power series associated to cycles index series _[BBL],
    .. see :mod:`sage.categories.cycle_index_series`)
    """

    def super_categories(self):
        return []

    class ParentMethods:

        @cached_method
        def coefficient(self, n):
            """
            The coefficient of degree `n`

            MATH::

            :param n: a non-negative integer

            (or may be just an integer cf Laurent series...)
            """