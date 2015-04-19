# -*- coding: utf-8 -*-
"""
Substitution of formal power series

Reference
---------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux
  Cambridge University Press, 1998

"""
# *****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ******************************************************************************
from sage.categories.formal_power_series import ExponentialPowerSeries
from sage.combinat.species2.formal_power_series import FPS
from sage.combinat.species2.formal_power_series.operations.add import Add
from sage.combinat.species2.formal_power_series.operations.product import ExponentialProd
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.arith import factorial
from sage.rings.integer import Integer


class Substitution(FPS):

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args):
        f, g = args
        # neutral element
        if f == ExponentialPowerSeries().singletons():
            return g
        if g == ExponentialPowerSeries().singletons():
            return f
        # associative: f(g(h)) |--> f(g)(h)
        if isinstance(g, Substitution):
            g, h = g._f_, g._g_
            return Substitution(Substitution(f, g), h)
        # distributive 1
        if isinstance(f, ExponentialProd):
            return ExponentialProd(*[(Substitution(h, g), nh) for (h, nh) in f._dic_fs_.iteritems()])
        # distributive 2
        if isinstance(f, Add):
            return Add(*[(Substitution(h, g), nh) for (h, nh) in f._dic_fs_.iteritems()])
        # # otherwise ##
        return super(Substitution, cls).__classcall__(cls, f, g)

    def __init__(self, f, g):
        assert(g._valuation_() > 0), "The formal power series should satisfy: `[t^0]g(t) = 0`."
        FPS.__init__(self, category=ExponentialPowerSeries())
        self._f_ = f
        self._g_ = g

        self._gen = self._gen_()
        self._gen.send(None)

    def _gen_(self):
        f = self._f_
        g = self._g_

        self._index_ = [(f.coefficient(0), ExponentialPowerSeries().one())] * (self._valuation_()+1)
        self._last_index_ = self._valuation_()

        while True:
            n = (yield)
            while self._last_index_ < n:
                self._last_index_ += 1

                # f_n
                coef = f.coefficient(self._last_index_) / factorial(self._last_index_)
                # t^n --> g(t)^n
                self._index_.append((coef, g**self._last_index_))
                # [t^n]f(g(t)) |---> [t^n](f_0 + f_1 g(t) + f_2 g(t)^2 + f_3 g(t)^3 + ...)

    @cached_method
    def coefficient(self, n):
        """
        MATH::

            \begin{align*}
                f(g(t)) &= \sum_{n \geqslant 0} f_n g(t)^n\\
                        &= \sum_{n \geqslant 0} \left(
                            \sum_{k | n} f_k [t^n]g(t)^k
                        \right) t^n
            \end{align*}

        """
        n = Integer(n)
        if n == 0:
            return self._f_.coefficient(0)

        self._gen.send(n)
        return sum(self._index_[k][0] * self._index_[k][1].coefficient(n) for k in range(n+1))

    def _valuation_(self):
        # FIXME: naive implementation... is it ok?
        return self._f_._valuation_() * self._g_._valuation_()