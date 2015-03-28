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
from itertools import imap
from sage.categories.category import Category
from sage.categories.objects import Objects
from sage.misc.cachefunc import cached_method
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer


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
        return [Objects()]

    def zero(self):
        from sage.combinat.species2.formal_power_series.some_characteristic_fps import ZeroFPS
        return ZeroFPS(category=self)

    def one(self):
        from sage.combinat.species2.formal_power_series.some_characteristic_fps import OneFPS
        return OneFPS(category=self)

    def singletons(self):
        from sage.combinat.species2.formal_power_series.some_characteristic_fps import SingletonsFPS
        return SingletonsFPS(category=self)

    def sets(self):
        from sage.combinat.species2.formal_power_series.some_characteristic_fps import SetsFPS
        return SetsFPS(category=self)

    def recursive_formal_power_series(self, name="f"):
        from sage.combinat.species2.formal_power_series.operations.recursive_fps import RecursiveFormalPowerSeries
        return RecursiveFormalPowerSeries(name=name, category=self)

    class ParentMethods:

        @cached_method
        def coefficient(self, n):
            """
            The coefficient of degree `n`

            MATH::

            :param n: a non-negative integer

            (or may be just an integer cf Laurent series...)
            """

        def coefficients(self, n):
            """
            Iterator of coefficients between `0` and `n`
            """
            return imap(self.coefficient, range(n+1))

        def add(f, g):
            """
            Sum of formal power series

            MATH::

                (f + g)(t) = f(t) + g(t) = \sum_{n \geqslant 0} (f_n + g_n)t^n\,.

            """
            from sage.combinat.species2.formal_power_series.operations.add import Add
            return Add((f,1), (g,1), category=f.category())

        __add__ = _add_ = add

        def __pow__(self, power):
            if power == 0:
                return self.category().one()
            elif power == 1:
                return self
            else:
                return self * self ** (power-1)

        def restricted(self, min=Integer(0), max=Infinity):
            """
            Restriction of formal power series
            """
            from sage.combinat.species2.formal_power_series.operations.restriction import Restriction
            return Restriction(self, min=min, max=max)

        def hadamard_product(f, g):
            """
            Hadamard product of formal power series

            MATH::

                f(t) \times g(t) = \sum_{n \geqslant 0} f_n g_n t^n\,.

            """
            from sage.combinat.species2.formal_power_series.operations.hadamard_product import HadamardProduct
            return HadamardProduct((f,1), (g,1), category=f.category())

        def derivative(f):
            """
            Derivative of formal power series

            MATH::

                f'(t) = \sum_{n \geqslant 1} n * f_n t^{n-1}\,.

            """
            from sage.combinat.species2.formal_power_series.operations.derivative import Derivative
            return Derivative(f)

        def pointing(self):
            """
            Pointing of formal power series

            MATH::

                f^\bullet (t) = t \frac{d}{dt}f(t)\,.

            """
            return self.category().singletons() * self.derivative()


class ExponentialPowerSeries(FormalPowerSeries):

    def super_categories(self):
        return [FormalPowerSeries()]

    class ParentMethods:

        def product(f, g):
            """
            Product of exponential formal power series

            MATH::

                (f \cdot g) = f(t) \cdot g(t) = \sum_{n \geqslant 0}
                \left(\sum_{i + j = n} \binom{n}{i} \cdot f_i \cdot g_j \right) t^n\,.

            with `f, g` are both exponential generating series.
            """
            from sage.combinat.species2.formal_power_series.operations.product import ExponentialProd
            return ExponentialProd((f,1), (g,1))

        __mul__ = _mul_ = product

        def substitution(f, g):
            """
            The plethystic substition of exponential formal power series

            MATH::

                (f\circ g)(t) = f(g(t))


            """
            from sage.combinat.species2.formal_power_series.operations.substitution import Substitution
            return Substitution(f, g)

        composition = substitution

        def functorial_composite(f, g):
            """
            Functorial composite of formal power series

            MATH::

                (f \Box g)(t) = \sum_{n \seq 0} f_{g_n} t^n\,.

            """
            from sage.combinat.species2.formal_power_series.operations.functorial_composite import FunctorialComposite
            return FunctorialComposite(f, g)

class OrdinaryPowerSeries(FormalPowerSeries):

    def super_categories(self):
        return [FormalPowerSeries()]

    class ParentMethods:

        def product(f, g):
            """
            Product of ordinary power series

            MATH::

                (f \cdot g) = f(t) \cdot g(t) = \sum_{n \geqslant 0}
                \left(\sum_{i + j = n} f_i \cdot g_j \right) t^n\,.

            with `f, g` are both ordinary generating series.
            """
            from sage.combinat.species2.formal_power_series.operations.product import OrdinaryProd
            return OrdinaryProd((f,1), (g,1))

        __mul__ = _mul_ = product
