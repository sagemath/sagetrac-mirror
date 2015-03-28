# -*- coding: utf-8 -*-
"""
Restriction of formal power series

"""
# *****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ******************************************************************************
from sage.categories.formal_power_series import ExponentialPowerSeries
from sage.combinat.species2.formal_power_series import FPS
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer


class Restriction(FPS):

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        f = args[0]

        mi, ma = opts["min"], opts["max"]

        if mi == ma:
            if mi == 0 and f.coefficient(mi) == 1:
                return f.category().one()
            if mi == 1 and f.coefficient(mi) == 1:
                return f.category().singletons()
        elif mi == 0 and ma == Infinity:
            return f

        return super(Restriction, cls).__classcall__(cls, f, **opts)

    def __init__(self, f, min=0, max=Infinity):
        FPS.__init__(self, category=f.category())
        self._f_ = f
        self._min_ = min
        self._max_ = max

    def coefficient(self, n):
        """
        MATH::

            [t^n]f_{[min, max]}(t) = \begin{dcases*}
                [t^n]f(t) & if `min \leqslant n \leqslant max`,\\
                 0        & otherwise.
            \end{dcases*}

        """
        if self._min_ <= n <= self._max_:
            return self._f_.coefficient(n)
        return Integer(0)

    def _valuation_(self):
        return max(self._f_._valuation_(), self._min_)  # FIXME It seems to be reasonable for human use...
